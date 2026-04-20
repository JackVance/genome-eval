"""Call mtDNA haplogroup for a subject using HaploGrep3 + PhyloTree.

Pipeline:
  1. Extract MT variants from the subject's chip parquet (imputation panels
     don't cover MT, so chip is authoritative).
  2. Generate a VCF with rCRS as reference; for each position where subject's
     allele differs from rCRS, emit ALT + genotype.
  3. Run HaploGrep3 classify with --chip flag (handles chip sparsity and
     alignment rules for partial data).
  4. Parse the top-hit haplogroup + quality metrics.
  5. Record a finding with the haplogroup, defining variants observed,
     quality score, and the full list of hits for transparency.

Input file conventions:
  standardized-genomes/<subject>.parquet     (chip parquet with MT)
  reference/haplogroups/mtdna/rCRS.fasta      (NC_012920.1)
  reference/haplogroups/mtdna/haplogrep3/     (HaploGrep3 JAR + tree)

Output files:
  standardized-genomes/haplogroups/<subject>.mtdna.vcf
  standardized-genomes/haplogroups/<subject>.mtdna.haplogrep.txt

Ledger entry: topic=`mtdna_haplogroup`, evidence_class=`mendelian_trait`.

Usage:
    python scripts/run_mtdna_haplogroup.py alice
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
sys.path.insert(0, str(SCRIPT_DIR))

import ledger_io
from tier_rules import compute_tier, TIER_RULE_VERSION

MT_REF_FASTA = PROJECT_ROOT / "reference" / "haplogroups" / "mtdna" / "rCRS.fasta"
HAPLOGREP_DIR = PROJECT_ROOT / "reference" / "haplogroups" / "mtdna" / "haplogrep3"
HAPLOGREP_JAR = HAPLOGREP_DIR / "haplogrep3.jar"
JAVA_EXE = PROJECT_ROOT / "reference" / "imputation" / "jdk" / "jdk-21.0.6+7" / "bin" / "java.exe"
OUT_DIR = PROJECT_ROOT / "standardized-genomes" / "haplogroups"
TREE_ID = "phylotree-rcrs@17.2"  # latest PhyloTree (rCRS version)


def load_rcrs_sequence() -> str:
    if not MT_REF_FASTA.exists() or MT_REF_FASTA.stat().st_size < 15_000:
        _download_rcrs()
    lines = MT_REF_FASTA.read_text().splitlines()
    return "".join(l for l in lines if not l.startswith(">"))


def _download_rcrs() -> None:
    """One-time fetch of the rCRS mtDNA reference from NCBI.

    NC_012920.1 is the revised Cambridge Reference Sequence — the canonical
    mtDNA reference PhyloTree is built against. 16,569 bp; the download is
    a tiny text file (~16 KB gzipped over HTTPS). Only called if the local
    copy is missing — subsequent runs read from disk.
    """
    import urllib.request
    url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        "?db=nuccore&id=NC_012920.1&rettype=fasta&retmode=text"
    )
    MT_REF_FASTA.parent.mkdir(parents=True, exist_ok=True)
    print(f"Fetching rCRS from NCBI → {MT_REF_FASTA.name} …")
    data = urllib.request.urlopen(url, timeout=60).read().decode()
    if len(data) < 15_000:
        raise RuntimeError(
            f"rCRS download looks too small ({len(data)} chars). Expected ~16.8 KB."
        )
    MT_REF_FASTA.write_text(data)
    print(f"  wrote {len(data)} chars ({MT_REF_FASTA.stat().st_size} bytes)")


def extract_mt_variants(subject_id: str) -> pd.DataFrame:
    """Load MT-chromosome rows from the chip parquet.

    Returns a dataframe with columns: rsid, pos, a1, a2 (for MT these are
    haploid so a1 == a2 always). The chip parquet is authoritative for MT;
    imputation panels don't cover it.
    """
    chip_path = PROJECT_ROOT / "standardized-genomes" / f"{subject_id}.parquet"
    df = pd.read_parquet(chip_path)
    df["chrom"] = df["chrom"].astype(str)
    mt = df[df["chrom"].isin(["MT", "chrM", "M"])].copy()
    mt = mt[mt["a1"].notna() & (mt["a1"] != "-")]
    # Collapse indel probes (23andMe D/I encoding isn't directly usable for mtDNA);
    # PhyloTree defining variants are SNPs and short indels represented canonically.
    mt = mt[mt["a1"].isin(["A", "C", "G", "T"])]
    mt = mt.sort_values("pos").reset_index(drop=True)
    return mt


def make_mt_vcf(mt_df: pd.DataFrame, rcrs: str, sample_id: str, out_path: Path) -> int:
    """Emit a minimal VCF of MT variants relative to rCRS.

    Returns the number of variant rows written. Only positions where the
    subject's allele differs from rCRS are emitted; unvaried positions would
    be treated as reference by HaploGrep3 regardless. The --chip flag on
    HaploGrep3 handles sparse input correctly.

    VCF header uses `chrMT` as the contig name — HaploGrep3 accepts this
    plus several aliases via its alignment rules.
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)
    n_written = 0
    with out_path.open("w", encoding="utf-8") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.write("##source=run_mtdna_haplogroup.py\n")
        fh.write("##contig=<ID=chrMT,length=16569>\n")
        fh.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        fh.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_id}\n")

        for row in mt_df.itertuples(index=False):
            pos = int(row.pos)
            if pos < 1 or pos > len(rcrs):
                continue
            ref_base = rcrs[pos - 1].upper()
            alt_base = row.a1.upper()
            if alt_base == ref_base:
                continue
            if alt_base not in ("A", "C", "G", "T"):
                continue
            # mtDNA is haploid → GT encoded as the single allele index "1".
            fh.write(
                f"chrMT\t{pos}\t{row.rsid}\t{ref_base}\t{alt_base}\t"
                f".\tPASS\t.\tGT\t1\n"
            )
            n_written += 1
    return n_written


def run_haplogrep(vcf_path: Path, out_path: Path) -> None:
    cmd = [
        str(JAVA_EXE),
        "-jar",
        str(HAPLOGREP_JAR),
        "classify",
        "--in",
        str(vcf_path),
        "--out",
        str(out_path),
        "--tree",
        TREE_ID,
        "--chip",
        "--extend-report",
    ]
    print(f"  $ {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print("HaploGrep3 stdout:", result.stdout[-1500:])
        print("HaploGrep3 stderr:", result.stderr[-1500:])
        raise RuntimeError(f"HaploGrep3 exited with code {result.returncode}")


def _haplogroup_context(haplogroup: str) -> dict:
    """Map a mtDNA haplogroup to a short descriptor + one-paragraph context.

    Keyed by the initial trunk letter + subclade prefix when common. The
    dataset is intentionally tiny — just the 10-15 macrohaplogroups a reader
    of this project is likely to see. For finer resolution the notes include
    PhyloTree defining variants so the user can investigate independently.
    """
    h = (haplogroup or "").upper()
    default = {"short": "haplogroup details in notes",
               "context": "See PhyloTree (https://www.phylotree.org) for lineage context."}
    table = [
        ("H", "pan-European common lineage",
         "~40-50% of modern Europeans; expanded from SW Europe after the Last Glacial Maximum ~15-20 kya."),
        ("U5", "oldest extant European mtDNA lineage",
         "U5 is the oldest mtDNA lineage still found in Europe today, diverged ~35-50 kya. "
         "U5a is more common in northern/northeastern Europe, U5b in western Europe. "
         "Enriched in pre-Neolithic hunter-gatherer remains."),
        ("U4", "ancient European hunter-gatherer lineage",
         "Common in Mesolithic remains across northern Europe; persisted in Baltic, Finno-Ugric populations."),
        ("U", "European or Near-Eastern ancient lineage",
         "U-branch mtDNA is ~45 kya; subclades (U1-U8) trace specific ancient migrations."),
        ("K", "post-LGM European / Near Eastern lineage",
         "K diverged from U8 ~35 kya; enriched in Ashkenazi Jewish populations (~30-40%)."),
        ("T", "Neolithic farmer-associated lineage",
         "T and T2 expanded in Europe with Neolithic agriculture from the Near East ~8-9 kya."),
        ("J", "Near Eastern / European lineage",
         "J and subclades expanded from the Near East into Europe in the Neolithic; J1c and J2b are European-enriched."),
        ("V", "Iberian / Scandinavian lineage",
         "V peaks in the Basque country, Sami, and northern Scandinavia; derived from pre-V (HV0)."),
        ("W", "European / Near Eastern low-frequency lineage", "~1-3% of Europeans."),
        ("X", "rare trans-continental lineage",
         "Found in Europe, Near East, and some Native American populations (X2a); distribution suggests complex deep history."),
        ("I", "northern European / North African lineage", "Low-frequency across Europe."),
        ("N", "macrohaplogroup out-of-Africa root",
         "Parent of most non-African mtDNA branches. Rarely seen directly in Europeans as a terminal call."),
        ("L", "African-origin lineage",
         "L0-L6 are the African macrohaplogroups. A European subject calling L-terminal would be unusual."),
        ("A", "East Asian / Native American lineage",
         "A peaks in Native American and East Asian populations."),
        ("B", "East Asian / Austronesian lineage", "B peaks in East Asia, Polynesia, and the Americas."),
        ("M", "Asian macrohaplogroup",
         "Parent of many Asian and Oceanian branches; rare in Europeans."),
    ]
    for prefix, short, context in table:
        if h.startswith(prefix):
            return {"short": short, "context": context}
    return default


def _strip_quoted(cell: str) -> str:
    """HaploGrep3 wraps every TSV cell in double-quotes; strip them."""
    c = cell.strip()
    if len(c) >= 2 and c[0] == '"' and c[-1] == '"':
        return c[1:-1]
    return c


def parse_haplogrep_output(path: Path) -> dict:
    """HaploGrep3 TSV output columns (extend-report mode):

    "SampleID" "Haplogroup" "Rank" "Quality" "Range" "Not_Found_Polys"
    "Found_Polys" "Remaining_Polys" "AAC_In_Remainings" "Input_Sample"

    Every cell is double-quoted in the output, so we strip quotes when
    reading — otherwise the keys lookup with raw names would silently miss.
    """
    rows = []
    with path.open(encoding="utf-8") as fh:
        header = [_strip_quoted(c) for c in fh.readline().rstrip("\n").split("\t")]
        for line in fh:
            parts = [_strip_quoted(c) for c in line.rstrip("\n").split("\t")]
            rows.append(dict(zip(header, parts)))
    if not rows:
        raise RuntimeError(f"HaploGrep3 output empty: {path}")
    return {"header": header, "rows": rows}


def main(subject_id: str) -> None:
    print(f"=== mtDNA haplogroup: {subject_id} ===")
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    rcrs = load_rcrs_sequence()
    print(f"rCRS: {len(rcrs)} bp")

    mt_df = extract_mt_variants(subject_id)
    print(f"MT chip variants: {len(mt_df)}")

    vcf_path = OUT_DIR / f"{subject_id}.mtdna.vcf"
    n = make_mt_vcf(mt_df, rcrs, subject_id, vcf_path)
    print(f"Wrote {vcf_path.name} with {n} non-reference variants")

    hg_path = OUT_DIR / f"{subject_id}.mtdna.haplogrep.txt"
    run_haplogrep(vcf_path, hg_path)
    print(f"HaploGrep3 output: {hg_path.name}")

    parsed = parse_haplogrep_output(hg_path)
    top = parsed["rows"][0]
    haplogroup = top.get("Haplogroup", "unknown")
    quality = top.get("Quality", "unknown")
    range_covered = top.get("Range", "unknown")
    print(f"\nTop haplogroup: {haplogroup}")
    print(f"Quality:        {quality}")
    print(f"Range:          {range_covered}")

    # Extract defining polymorphisms from the extend-report columns. Column
    # names vary by HaploGrep3 version; fall back to best alternatives.
    found_polys = top.get("Found_Polys") or ""
    not_found_polys = top.get("Not_Found_Polys") or top.get("Remaining_Polys") or ""
    input_sample = top.get("Input_Sample") or ""
    aac = top.get("AAC_In_Remainings") or ""

    # Also pull any lower-ranked alternatives for transparency. HaploGrep3
    # only emits the top rank unless we pass --hits, but the report may
    # contain multiple rows if the user did.
    all_hits = [(r.get("Haplogroup"), r.get("Quality")) for r in parsed["rows"]]

    # Look up haplogroup background to contextualize the call. Keep to a
    # short summary — deep ancestry context belongs in the notes, not a
    # reference table.
    hg_context = _haplogroup_context(haplogroup)

    claim = (
        f"mtDNA haplogroup: {haplogroup} ({hg_context['short']}, quality {quality})"
    )
    notes = (
        f"mtDNA haplogroup call via HaploGrep3 3.2.2 against PhyloTree "
        f"{TREE_ID}. Subject chip MT variants: {len(mt_df)} positions on "
        f"23andMe v5; {n} non-reference relative to rCRS. Quality score: "
        f"{quality} (1.00 = all defining variants of the haplogroup observed; "
        f"lower values indicate missing expected polymorphisms due to chip "
        f"sparsity). Range covered: {range_covered}. Found defining "
        f"polymorphisms ({len(found_polys.split()) if found_polys else 0}): "
        f"{found_polys or 'n/a'}. Missing expected polymorphisms "
        f"({len(not_found_polys.split()) if not_found_polys else 0}): "
        f"{not_found_polys or 'n/a'}. Background: {hg_context['context']}"
    )

    # Supersede any prior active mtDNA finding for this subject — re-running
    # should produce one canonical active record, not a pile of duplicates.
    prior = [
        r for r in ledger_io.load_active_findings(subject_id=subject_id)
        if r.get("topic") == "mtdna_haplogroup"
    ]
    supersedes_id = prior[0]["finding_id"] if prior else None

    rec = {
        "subject_id": subject_id,
        "topic": "mtdna_haplogroup",
        "supersedes": supersedes_id,
        "claim": claim,
        "variants": [{"chrom": "MT", "n_variants": n, "tool": "HaploGrep3",
                      "tree": TREE_ID, "rsid": "mtdna_bundle"}],
        "effect": {
            "type": "haplogroup",
            "lineage": "mtDNA",
            "haplogroup": haplogroup,
            "quality": quality,
            "range": range_covered,
            "found_polys": found_polys,
            "missing_polys": not_found_polys,
            "alternative_hits": all_hits[:5],
            "context_short": hg_context["short"],
        },
        "cohort_ancestry": "global",
        "subject_ancestry_match": "match",
        "source_ids": ["haplogrep3:phylotree-rcrs-17.2"],
        "notes": notes,
        "evidence_class": "mendelian_trait",
        "replication_count": 0,
        "inference_confidence": "high",
        "clinvar_significance": None,
        "clinvar_review_stars": None,
        "pvalue": None,
        "n_cases": None,
        "n_controls": None,
        "odds_ratio": None,
        "investigation_id": None,
    }
    tier, downgrade = compute_tier(rec)
    rec["tier_computed"] = tier
    rec["tier_computed_at"] = datetime.now(timezone.utc).isoformat()
    rec["tier_rule_version"] = TIER_RULE_VERSION
    rec["ancestry_downgrade"] = downgrade

    ledger_io.append_source(
        source_id="haplogrep3:phylotree-rcrs-17.2",
        kind="tool",
        url="https://github.com/genepi/haplogrep3",
        accessed_at=datetime.now(timezone.utc).isoformat(),
        citation="Weissensteiner H et al. HaploGrep 3 - Phylogenetic analysis of mitochondrial DNA. NAR (2024).",
        ancestry_cohort="global",
    )
    fid = ledger_io.append_finding(**rec)
    print(f"\nAppended finding {fid} (tier {tier})")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("subject_id")
    args = parser.parse_args()
    main(args.subject_id)
