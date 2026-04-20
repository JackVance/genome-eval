"""Call Y-chromosome haplogroup for a male subject using yhaplo (ISOGG tree).

Pipeline:
  1. Extract Y-chromosome variants from the subject's chip parquet (imputation
     panels don't cover Y, so chip is authoritative).
  2. Generate a yhaplo sample-major text input (.genos.txt) — simpler than
     indexed VCF for single-subject single-chromosome data.
  3. Run yhaplo; it walks the ISOGG tree and emits short-form + long-form
     haplogroup calls.
  4. Parse output, record finding with the classification + defining SNPs
     observed + alternative hits.

Input:
  standardized-genomes/<subject>.parquet        (chip parquet with Y)

Output:
  standardized-genomes/haplogroups/<subject>.ychr.genos.txt
  standardized-genomes/haplogroups/yhaplo_out/<subject>/*.txt   (yhaplo dir)

Usage:
    python scripts/run_y_haplogroup.py alice
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

OUT_DIR = PROJECT_ROOT / "standardized-genomes" / "haplogroups"


def extract_y_variants(subject_id: str) -> pd.DataFrame:
    """Return Y-chromosome rows from the chip parquet, genotype as a single
    letter per site. Chip Y SNPs are effectively haploid (male); 23andMe
    reports the single allele in both a1 and a2 columns.
    """
    chip_path = PROJECT_ROOT / "standardized-genomes" / f"{subject_id}.parquet"
    df = pd.read_parquet(chip_path)
    df["chrom"] = df["chrom"].astype(str)
    y = df[df["chrom"].isin(["Y", "chrY", "24"])].copy()
    y = y[y["a1"].notna() & (y["a1"] != "-")]
    y = y[y["a1"].isin(["A", "C", "G", "T"])]
    y = y.sort_values("pos").reset_index(drop=True)
    return y


def write_genos_txt(y_df: pd.DataFrame, sample_id: str, out_path: Path) -> int:
    """Write yhaplo sample-major input (.genos.txt format).

    File format (discovered from the yhaplo repo fixtures):
        Line 0: "ID" + tab-separated physical positions (GRCh37).
        Line 1+: sample_id + tab-separated genotype letters (A/C/G/T or ".").

    The format docs in --help say "whitespace-separated" but the reference
    example is tab-delimited; whitespace-split fails to parse positions
    correctly on multi-character input.
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)
    positions = y_df["pos"].astype(int).tolist()
    genotypes = y_df["a1"].tolist()
    with out_path.open("w", encoding="utf-8") as fh:
        fh.write("ID\t" + "\t".join(str(p) for p in positions) + "\n")
        fh.write(sample_id + "\t" + "\t".join(genotypes) + "\n")
    return len(positions)


def run_yhaplo(genos_path: Path, yhaplo_out_dir: Path) -> None:
    yhaplo_out_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        "yhaplo",
        "--input",
        str(genos_path),
        "--out_dir",
        str(yhaplo_out_dir),
        "--all_aux_output",    # emits paths, derived_snps, ancestral_snps, counts — everything we need
    ]
    print(f"  $ {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print("yhaplo stdout:", result.stdout[-1500:])
        print("yhaplo stderr:", result.stderr[-1500:])
        raise RuntimeError(f"yhaplo exited with code {result.returncode}")
    # Dump stdout for the log
    if result.stdout:
        print(result.stdout[-800:])


def parse_yhaplo_output(yhaplo_out_dir: Path, sample_id: str) -> dict:
    """yhaplo writes several files in the output dir:
      - haplogroups.<tag>.txt        (per-sample short-form + YCC-long-form)
      - paths.<tag>.txt              (the tree path traversed)
      - derived_snps.<tag>.txt       (defining SNPs called derived)
      - ancestral_snps.<tag>.txt     (defining SNPs called ancestral)

    Returns a dict with the parsed fields.
    """
    # Find the haplogroups file (there should be exactly one)
    hg_files = list(yhaplo_out_dir.glob("haplogroups*.txt"))
    if not hg_files:
        raise RuntimeError(f"No haplogroups output in {yhaplo_out_dir}")
    hg_file = hg_files[0]
    lines = hg_file.read_text().strip().splitlines()
    # Format: "<sample> <hg_short_form> <hg_intermediate> <hg_long_form>" —
    # yhaplo's haplogroups.txt uses padded-space alignment, not tab or single
    # space, so use unlimited whitespace-splitting.
    parts = lines[0].split()
    if len(parts) < 4 or parts[0] != sample_id:
        raise RuntimeError(f"Unexpected yhaplo output line: {lines[0]!r}")
    short = parts[1]
    intermediate = parts[2]
    long_form = parts[3]

    # Path traversal (which nodes in the ISOGG tree the call visited)
    path_lines = []
    path_files = list(yhaplo_out_dir.glob("paths*.txt"))
    if path_files:
        pt = path_files[0].read_text().strip().splitlines()
        # Path files have one line per sample, first col is sample_id
        for line in pt:
            if line.startswith(sample_id + " ") or line.startswith(sample_id + "\t"):
                path_lines.append(line.split(None, 1)[1] if len(line.split(None, 1)) > 1 else "")
                break

    # Derived / ancestral SNPs
    derived = []
    anc = []
    for f in yhaplo_out_dir.glob("derived_snps*.txt"):
        for line in f.read_text().strip().splitlines():
            p = line.split()
            if p and p[0] == sample_id:
                derived = p[1:]
                break
    for f in yhaplo_out_dir.glob("ancestral_snps*.txt"):
        for line in f.read_text().strip().splitlines():
            p = line.split()
            if p and p[0] == sample_id:
                anc = p[1:]
                break

    return {
        "short": short,
        "intermediate": intermediate,
        "long_form": long_form,
        "path": path_lines[0] if path_lines else "",
        "n_derived_snps": len(derived),
        "n_ancestral_snps": len(anc),
        "derived_snps_sample": derived[:30],
    }


def _y_haplogroup_context(short: str, long_form: str = "") -> dict:
    """Brief one-paragraph context per major Y haplogroup for the notes block.

    Keyed by the YCC long-form prefix (R1b, R1a, I1, I2, J2, E1b, G2a, etc.).
    Short-forms like "R-M269" don't start with the macro-haplogroup letters,
    so we match against the long form when present and fall back to the
    short form otherwise. Intentionally compact — deeper ancestry context
    is a rabbit hole that belongs in an external ISOGG / YFull lookup.
    """
    candidates = [c for c in (long_form, short) if c]
    h_candidates = [c.upper() for c in candidates]
    table = [
        ("R1B", "Western European patrilineage",
         "R1b-M269 is the most common Y haplogroup in Western Europe (>60% in Britain, Ireland, France, Spain). "
         "Subclades U106 (Germanic/Anglo-Saxon), P312 (Celtic/Iberian), L21 (Insular Celtic), and U152 (Alpine) "
         "trace specific Bronze-Age expansions."),
        ("R1A", "Eastern European / South Asian patrilineage",
         "R1a expanded from Eastern Europe / Central Asia ~4-5 kya. Z283 (Slavic), Z282 (Baltic), "
         "Z93 (South Asian / Indo-Iranian) are the major subclades."),
        ("I1", "Northern European patrilineage",
         "I1-M253 is found at highest frequency in Scandinavia (~35-40%), spread across northern Europe with "
         "Germanic migrations. Relatively young (~5 kya) starburst expansion."),
        ("I2", "European patrilineage, ancient",
         "I2-M438 is one of the oldest Y lineages still common in Europe. I2a (Balkan), I2c2 (British Isles) "
         "trace pre-Bronze-Age European populations."),
        ("J2", "Near Eastern / Mediterranean patrilineage",
         "J2-M172 expanded from Anatolia / Levant during the Neolithic farming dispersal; common in "
         "southern Europe, the Near East, and parts of South Asia."),
        ("J1", "Arabian / Semitic patrilineage",
         "J1-M267 peaks in the Arabian peninsula and among Semitic-speaking populations."),
        ("E1B", "North African / European patrilineage",
         "E1b1b-M215 originated in Africa; subclade M81 is Berber/North African, V13 spread into Europe."),
        ("G2A", "Neolithic farmer patrilineage",
         "G2a-P15 is associated with early Neolithic European farmers (Cardial/LBK cultures); enriched in "
         "ancient European remains from ~7-8 kya."),
        ("N", "Finno-Ugric / Siberian patrilineage",
         "N-M231 is dominant among Finno-Ugric populations and in eastern Eurasia."),
        ("Q", "Native American / Central Asian patrilineage",
         "Q-M242 is the primary Y lineage of Native Americans; Q1b is broadly Eurasian."),
        ("T", "Near Eastern / African patrilineage",
         "T-M184 is low-frequency in Europe, enriched in parts of the Near East and East Africa."),
        ("A", "oldest African patrilineage",
         "A is the deepest-rooted Y lineage, primarily African. A European subject calling A-terminal would "
         "suggest deep African ancestry or a recent Sub-Saharan paternal line."),
        ("B", "Sub-Saharan African patrilineage", "B-M181 is primarily African."),
        ("C", "Asian / Oceanian patrilineage",
         "C-M130 is found in East Asia, Central Asia, and Oceania; C2 is associated with Mongolic/Turkic expansions."),
    ]
    for prefix, short_d, context in table:
        for h in h_candidates:
            if h.startswith(prefix):
                return {"short": short_d, "context": context}
    return {"short": "see ISOGG tree for lineage context",
            "context": "See https://isogg.org/tree/ for the full Y-chromosome phylogenetic tree."}


def main(subject_id: str) -> None:
    print(f"=== Y haplogroup: {subject_id} ===")
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    y_df = extract_y_variants(subject_id)
    print(f"Y chip variants: {len(y_df)}")

    genos_path = OUT_DIR / f"{subject_id}.ychr.genos.txt"
    n = write_genos_txt(y_df, subject_id, genos_path)
    print(f"Wrote {genos_path.name} with {n} Y positions")

    yhaplo_out = OUT_DIR / "yhaplo_out" / subject_id
    run_yhaplo(genos_path, yhaplo_out)

    parsed = parse_yhaplo_output(yhaplo_out, subject_id)
    short = parsed["short"]
    print(f"\nShort-form haplogroup: {short}")
    print(f"YCC long-form:          {parsed['long_form']}")
    print(f"Intermediate:           {parsed['intermediate']}")
    print(f"Derived SNPs called:    {parsed['n_derived_snps']}")
    print(f"Ancestral SNPs called:  {parsed['n_ancestral_snps']}")

    context = _y_haplogroup_context(short, long_form=parsed["long_form"])
    claim = f"Y-chromosome haplogroup: {short} ({context['short']})"
    notes = (
        f"Y haplogroup call via yhaplo 2.1.17 against ISOGG tree. Subject "
        f"chip Y variants: {n}. Derived (variant-allele) SNPs: "
        f"{parsed['n_derived_snps']}; ancestral (reference-allele) SNPs: "
        f"{parsed['n_ancestral_snps']}. Short-form: {short}. YCC long-form: "
        f"{parsed['long_form']}. Background: {context['context']} "
        f"Representative derived SNPs observed: "
        f"{', '.join(parsed['derived_snps_sample'][:15]) or 'n/a'}."
    )

    # Supersede any prior active Y finding
    prior = [r for r in ledger_io.load_active_findings(subject_id=subject_id)
             if r.get("topic") == "y_haplogroup"]
    supersedes_id = prior[0]["finding_id"] if prior else None

    rec = {
        "subject_id": subject_id,
        "topic": "y_haplogroup",
        "supersedes": supersedes_id,
        "claim": claim,
        "variants": [{"chrom": "Y", "n_variants": n, "tool": "yhaplo",
                      "tree": "ISOGG", "rsid": "ychr_bundle"}],
        "effect": {
            "type": "haplogroup",
            "lineage": "Y",
            "haplogroup_short": short,
            "haplogroup_long": parsed["long_form"],
            "haplogroup_intermediate": parsed["intermediate"],
            "n_derived_snps": parsed["n_derived_snps"],
            "n_ancestral_snps": parsed["n_ancestral_snps"],
            "context_short": context["short"],
        },
        "cohort_ancestry": "global",
        "subject_ancestry_match": "match",
        "source_ids": ["yhaplo:ISOGG"],
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
        source_id="yhaplo:ISOGG",
        kind="tool",
        url="https://github.com/23andMe/yhaplo",
        accessed_at=datetime.now(timezone.utc).isoformat(),
        citation="Poznik GD. Identifying Y-chromosome haplogroups in arbitrarily large samples of sequenced or genotyped men. bioRxiv (2016).",
        ancestry_cohort="global",
    )
    fid = ledger_io.append_finding(**rec)
    print(f"\nAppended finding {fid} (tier {tier})")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("subject_id")
    args = parser.parse_args()
    main(args.subject_id)
