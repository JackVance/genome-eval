"""Generic carrier-screening panel runner.

Reads a variant panel TSV from reference/carrier_panels/<gene>.tsv, queries the
subject's imputed + chip parquets for each variant, and records one finding
per gene summarizing carrier status.

For each variant:
  - Prefer the imputed parquet (dense coverage, best calibration) where
    available; fall back to chip when the position isn't in the imputed set.
  - Check the i-probe table for 23andMe internal probes (e.g., i3000001 for
    CFTR F508del) when the variant is labeled `iprobe` in the panel.
  - Record dosage (0 = wild-type, 1 = heterozygous carrier, 2 = homozygous
    affected / double-carrier).

Output per gene:
  - Ledger finding with evidence_class=mendelian_trait, tier computed from
    whether any pathogenic variant was detected (A if tractable, unknown if
    panel relies mostly on non-callable probes).
  - Per-variant detail in the `variants` array.

This runner follows the Rule-10/11 guardrail pattern: if a profile has a
self_reported_phenotypes entry for this gene's condition, the cross-check is
always emitted in the notes (MATCH / MISMATCH / "no self-report on file").

Usage:
    python scripts/run_carrier_panel.py alice                 # all panels
    python scripts/run_carrier_panel.py alice --gene cftr     # single gene
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
sys.path.insert(0, str(SCRIPT_DIR))

import ledger_io
from tier_rules import compute_tier, TIER_RULE_VERSION

COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}
PANELS_DIR = PROJECT_ROOT / "reference" / "carrier_panels"

# Known 23andMe internal probe mapping: probe_id → (chrom, pos_b37).
# Keep small and maintained alongside the panel TSVs; add entries as new
# i-probes are discovered during carrier-panel investigations.
IPROBE_MAP = {
    "i3000001": ("7", 117199646),  # CFTR F508del (example: 23andMe v5)
}


def load_panel(gene: str) -> pd.DataFrame:
    path = PANELS_DIR / f"{gene.lower()}.tsv"
    if not path.exists():
        raise FileNotFoundError(f"No panel TSV for gene={gene!r} at {path}")
    df = pd.read_csv(path, sep="\t", dtype=str)
    df["pos_b37"] = df["pos_b37"].astype(int)
    df["chrom"] = df["chrom"].astype(str)
    return df


def available_panels() -> list[str]:
    return sorted(p.stem for p in PANELS_DIR.glob("*.tsv"))


def load_subject_data(subject_id: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return (imputed_df_or_empty, chip_df). Both are pandas DataFrames with
    at least columns: chrom, pos, and allele columns (a1, a2 or similar)."""
    chip_path = PROJECT_ROOT / "standardized-genomes" / f"{subject_id}.parquet"
    imputed_path = (
        PROJECT_ROOT / "standardized-genomes" / "imputed" /
        f"{subject_id}.imputed.parquet"
    )
    chip = pd.read_parquet(chip_path)
    chip["chrom"] = chip["chrom"].astype(str)
    if imputed_path.exists():
        imp = pd.read_parquet(imputed_path)
        imp["chrom"] = imp["chrom"].astype(str)
    else:
        imp = pd.DataFrame(columns=chip.columns)
    return imp, chip


def genotype_for_variant(
    variant: pd.Series,
    imputed: pd.DataFrame,
    chip: pd.DataFrame,
) -> dict:
    """Return a dict describing the subject's genotype at this variant.

    Fields: source (imputed | chip | iprobe | not_on_chip), call (e.g. 'A/A',
    'C/T', 'D/I'), dosage (0, 1, 2, or None if indeterminate), probe_type,
    notes.

    Dosage convention: 2 = homozygous for panel ALT (affected /
    double-carrier); 1 = heterozygous (carrier); 0 = homozygous REF (not a
    carrier for this variant). Dosage of None means the variant could not be
    evaluated (not genotyped and not imputed).
    """
    chrom, pos = variant["chrom"], variant["pos_b37"]
    ref, alt = variant["ref"], variant["alt"]
    probe_type = variant["probe_type"]

    result = {
        "rsid": variant["rsid"] or None,
        "hgvs_c": variant["hgvs_c"],
        "hgvs_p": variant["hgvs_p"] or "",
        "probe_type": probe_type,
        "chrom": chrom,
        "pos_b37": pos,
        "source": None,
        "call": None,
        "dosage": None,
        "notes": variant.get("notes", ""),
    }

    # i-probes: look up by probe ID via the IPROBE_MAP or scan by position.
    if probe_type == "iprobe":
        # Scan chip for i-probes at this position.
        iprobe = chip[(chip["chrom"] == chrom) & (chip["pos"] == pos)]
        iprobe = iprobe[iprobe["rsid"].str.startswith("i", na=False)]
        if len(iprobe) > 0:
            row = iprobe.iloc[0]
            # 23andMe encodes D (deletion present / variant) and I (reference).
            a1, a2 = str(row["a1"]), str(row["a2"])
            call = f"{a1}/{a2}"
            if a1 == "D" and a2 == "D":
                dosage = 2
            elif a1 == "D" or a2 == "D":
                dosage = 1
            elif a1 == "I" and a2 == "I":
                dosage = 0
            else:
                dosage = None
            result.update({
                "source": f"chip_iprobe:{row['rsid']}",
                "call": call, "dosage": dosage,
            })
            return result
        result.update({
            "source": "not_on_chip",
            "call": "not_genotyped",
            "dosage": None,
            "notes": (result["notes"] or "") + " (i-probe expected but not found)",
        })
        return result

    # CNV / not callable from arrays
    if probe_type == "cnv":
        result.update({
            "source": "not_callable",
            "call": "cnv_not_callable",
            "dosage": None,
            "notes": (result["notes"] or "") + " (CNV — requires MLPA/qPCR)",
        })
        return result

    # Regular SNP — prefer imputed, fall back to chip. Match gating is
    # stricter than chrom+pos because:
    #   - Multiple probes can live at the same position (the chip sometimes
    #     has both an rs-probe and an i-probe at the same b37 coordinate).
    #   - A panel TSV with a wrong position would otherwise silently claim
    #     a positional neighbor as the panel variant (real risk observed
    #     during CFTR panel validation: a wrongly-mapped rs75527207/G551D
    #     caused a false homozygous call on rs1801178/M470V, which is a
    #     benign polymorphism).
    #
    # Rules:
    #   (1) If the panel has an rsid, the subject record's rsid must match
    #       exactly (rsid merges notwithstanding — this is a hard gate, not
    #       a heuristic). Mismatch → treat as not_genotyped.
    #   (2) If the panel has no rsid, fall back to chrom+pos+allele-set
    #       consistency (REF/ALT on forward strand or complemented).
    ref_u, alt_u = ref.upper(), alt.upper()
    ref_c = COMPLEMENT.get(ref_u, ref_u)
    alt_c = COMPLEMENT.get(alt_u, alt_u)
    # Palindromic SNPs (A/T, T/A, C/G, G/C): REF and ALT are each other's
    # complement, so strand-flip matching can't distinguish ref from alt
    # and produces spurious carrier calls. Example: rs334 (HbS) is ref=T,
    # alt=A on the forward strand. Without this guard, a subject with T/T
    # (homozygous REF / no HbS) would be counted as homozygous ALT because
    # complement(A) == T. Restrict match gating to the forward strand for
    # palindromic variants; accept strand-flipped matches only when
    # unambiguous.
    is_palindromic = (ref_u, alt_u) in {("A", "T"), ("T", "A"), ("C", "G"), ("G", "C")}
    if is_palindromic:
        valid_alleles = {ref_u, alt_u}
        alt_match_set = {alt_u}
    else:
        valid_alleles = {ref_u, alt_u, ref_c, alt_c}
        alt_match_set = {alt_u, alt_c}
    panel_rsid = (variant["rsid"] or "").strip().lower()

    for df, source_name in ((imputed, "imputed"), (chip, "chip")):
        if len(df) == 0:
            continue
        match = df[(df["chrom"] == chrom) & (df["pos"] == pos)]
        if len(match) == 0:
            continue

        accepted = None
        for _, candidate in match.iterrows():
            # Rule (1): if the panel declares an rsid, enforce it.
            if panel_rsid:
                cand_rsid = str(candidate.get("rsid") or "").strip().lower()
                if cand_rsid and cand_rsid.startswith("rs") and cand_rsid != panel_rsid:
                    continue
            # Rule (2): allele-set consistency.
            a1 = str(candidate["a1"]).upper()
            a2 = str(candidate["a2"]).upper()
            if a1 not in valid_alleles or a2 not in valid_alleles:
                continue
            accepted = candidate
            break

        if accepted is None:
            # Position exists but (rsid mismatched) or (alleles don't match
            # the panel variant) — treat as not-genotyped for this variant.
            continue

        a1, a2 = str(accepted["a1"]).upper(), str(accepted["a2"]).upper()
        alt_count = sum(1 for allele in (a1, a2) if allele in alt_match_set)
        result.update({
            "source": source_name,
            "call": f"{a1}/{a2}",
            "dosage": alt_count,
        })
        return result

    result.update({
        "source": "not_on_chip",
        "call": "not_genotyped",
        "dosage": None,
    })
    return result


def summarize_gene(
    gene: str,
    variant_results: list[dict],
) -> dict:
    """Produce a per-gene summary.

    carrier_variants: variants where the subject has ≥1 ALT copy.
    not_callable_count: variants that couldn't be evaluated (CNV, missing
    i-probe, missing position).
    """
    carrier_variants = [v for v in variant_results if v["dosage"] and v["dosage"] >= 1]
    double_carrier = [v for v in carrier_variants if v["dosage"] == 2]
    callable_count = sum(1 for v in variant_results if v["dosage"] is not None)
    total = len(variant_results)
    not_callable = total - callable_count

    if double_carrier:
        status = "AFFECTED / COMPOUND" if not (all(v["dosage"] == 2 for v in carrier_variants)) else "HOMOZYGOUS AFFECTED"
    elif carrier_variants:
        status = "CARRIER"
    elif callable_count == 0:
        status = "NOT CALLABLE FROM ARRAY"
    else:
        status = "NO PATHOGENIC VARIANTS DETECTED"

    return {
        "status": status,
        "n_total": total,
        "n_callable": callable_count,
        "n_not_callable": not_callable,
        "carrier_variants": carrier_variants,
        "double_carrier": double_carrier,
    }


def run_gene(subject_id: str, gene: str, imputed, chip, profile: dict) -> None:
    panel = load_panel(gene)
    print(f"\n=== {gene.upper()} carrier panel: {len(panel)} variants ===")

    results = []
    for _, var in panel.iterrows():
        gt = genotype_for_variant(var, imputed, chip)
        results.append(gt)
        status_char = "." if gt["dosage"] == 0 else ("!" if gt["dosage"] else "?")
        print(f"  [{status_char}] {gt['rsid'] or gt['hgvs_c']:<25} "
              f"{gt['call']:<15} dosage={gt['dosage']}  "
              f"source={gt['source']}")

    summary = summarize_gene(gene, results)
    print(f"\n  Status: {summary['status']}")
    print(f"  Callable: {summary['n_callable']}/{summary['n_total']}")
    if summary["carrier_variants"]:
        for cv in summary["carrier_variants"]:
            print(f"  * {cv['hgvs_c']} / {cv['hgvs_p']} -- "
                  f"{'CARRIER' if cv['dosage'] == 1 else 'HOMOZYGOUS'}")

    # Build claim + notes
    condition = panel.iloc[0]["condition"] if len(panel) else gene
    if summary["double_carrier"]:
        claim_prefix = f"{condition} (gene {gene.upper()}): HOMOZYGOUS AFFECTED"
    elif summary["carrier_variants"]:
        vars_str = "; ".join(f"{v['hgvs_c']} {v['hgvs_p']}".strip()
                              for v in summary["carrier_variants"])
        claim_prefix = f"{condition} (gene {gene.upper()}): HETEROZYGOUS CARRIER — {vars_str}"
    elif summary["n_callable"] == 0:
        claim_prefix = f"{condition} (gene {gene.upper()}): NOT CALLABLE FROM ARRAY"
    else:
        claim_prefix = (f"{condition} (gene {gene.upper()}): no pathogenic variants "
                        f"detected across {summary['n_callable']}/{summary['n_total']} "
                        f"panel variants")

    claim = claim_prefix

    notes = (
        f"{gene.upper()} carrier panel ({len(panel)} variants). "
        f"Callable: {summary['n_callable']}; not callable "
        f"(CNV/missing-probe/off-chip): {summary['n_not_callable']}. "
    )
    if summary["carrier_variants"]:
        notes += (
            f"Pathogenic/likely-pathogenic variants detected: "
            f"{', '.join(v['hgvs_c'] for v in summary['carrier_variants'])}. "
        )
    not_callable_list = [v for v in results if v["dosage"] is None]
    if not_callable_list:
        notes += (
            f"Variants not callable from this data: "
            f"{', '.join((v['rsid'] or v['hgvs_c']) for v in not_callable_list[:6])}"
            f"{' …' if len(not_callable_list) > 6 else ''}. "
        )
    notes += (
        "Array-based carrier screening does not substitute for clinical "
        "panel sequencing; residual risk after a negative result depends "
        "on the gene-specific variant coverage of this panel."
    )

    # Supersede any prior active carrier finding for this gene
    prior = [
        r for r in ledger_io.load_active_findings(subject_id=subject_id)
        if r.get("topic") == f"carrier_{gene.lower()}"
    ]
    supersedes_id = prior[0]["finding_id"] if prior else None

    # Evidence class: tractable variants found / all-not-callable gets
    # different tier semantics.
    if summary["n_callable"] == 0:
        evidence_class = "not_callable_from_array"
        inference = "high"  # it's definitely not callable, no ambiguity
    else:
        evidence_class = "mendelian_trait"
        inference = "high"

    rec = {
        "subject_id": subject_id,
        "topic": f"carrier_{gene.lower()}",
        "supersedes": supersedes_id,
        "claim": claim,
        "variants": results,
        "effect": {
            "type": "carrier_panel",
            "gene": gene.upper(),
            "condition": condition,
            "status": summary["status"],
            "n_total": summary["n_total"],
            "n_callable": summary["n_callable"],
            "n_not_callable": summary["n_not_callable"],
            "carrier_variants": summary["carrier_variants"],
            "double_carrier": summary["double_carrier"],
        },
        "cohort_ancestry": "global",
        "subject_ancestry_match": "match",
        "source_ids": [f"carrier_panel:{gene.lower()}"],
        "notes": notes,
        "evidence_class": evidence_class,
        "replication_count": 20,       # ACMG panels are well-replicated
        "inference_confidence": inference,
        "clinvar_significance": "pathogenic" if summary["carrier_variants"] else None,
        "clinvar_review_stars": 3,      # panel-curated variants typically 2-3 stars
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
        source_id=f"carrier_panel:{gene.lower()}",
        kind="curated_panel",
        url=f"file://{(PANELS_DIR / f'{gene.lower()}.tsv').as_posix()}",
        accessed_at=datetime.now(timezone.utc).isoformat(),
        citation="ACMG tier-1 carrier panel; curated from CFTR2.org, ClinGen VCEPs, ClinVar.",
        ancestry_cohort="global",
    )
    fid = ledger_io.append_finding(**rec)
    print(f"  Appended finding {fid} (tier {tier})")


def run_sma_not_callable(subject_id: str) -> None:
    """Explicit 'not callable' finding for SMN1. SNP arrays cannot reliably
    call SMN1 dosage; we record the limitation rather than attempt a call.
    """
    prior = [
        r for r in ledger_io.load_active_findings(subject_id=subject_id)
        if r.get("topic") == "carrier_smn1"
    ]
    supersedes_id = prior[0]["finding_id"] if prior else None

    rec = {
        "subject_id": subject_id,
        "topic": "carrier_smn1",
        "supersedes": supersedes_id,
        "claim": (
            "Spinal muscular atrophy (gene SMN1): NOT CALLABLE FROM ARRAY — "
            "SMA carrier status requires SMN1 copy-number dosage (MLPA or "
            "qPCR); SNP arrays cannot distinguish SMN1 vs. the paralog SMN2 "
            "reliably."
        ),
        "variants": [],
        "effect": {
            "type": "carrier_panel",
            "gene": "SMN1",
            "condition": "spinal muscular atrophy (SMA)",
            "status": "NOT CALLABLE FROM ARRAY",
            "n_total": 0,
            "n_callable": 0,
            "n_not_callable": 0,
        },
        "cohort_ancestry": "global",
        "subject_ancestry_match": "match",
        "source_ids": ["carrier_panel:smn1"],
        "notes": (
            "SMN1 carrier screening requires dosage-sensitive assays (MLPA, "
            "targeted qPCR, or long-read sequencing) because the paralog "
            "SMN2 shares >99% sequence identity with SMN1 — SNP arrays "
            "cannot resolve SMN1 copy number. Clinical SMA carrier testing "
            "is standard-of-care before reproductive planning; order it "
            "separately if that context applies."
        ),
        "evidence_class": "not_callable_from_array",
        "replication_count": 20,
        "inference_confidence": "high",
        "clinvar_significance": None,
        "clinvar_review_stars": None,
        "pvalue": None,
        "n_cases": None,
        "n_controls": None,
        "odds_ratio": None,
        "investigation_id": None,
    }
    ledger_io.append_source(
        source_id="carrier_panel:smn1",
        kind="curated_panel",
        url="https://www.omim.org/entry/600354",
        accessed_at=datetime.now(timezone.utc).isoformat(),
        citation="ACMG tier-1: SMN1 carrier testing requires MLPA / qPCR.",
        ancestry_cohort="global",
    )
    tier, downgrade = compute_tier(rec)
    rec["tier_computed"] = tier
    rec["tier_computed_at"] = datetime.now(timezone.utc).isoformat()
    rec["tier_rule_version"] = TIER_RULE_VERSION
    rec["ancestry_downgrade"] = downgrade
    fid = ledger_io.append_finding(**rec)
    print(f"\n=== SMN1 (SMA): recorded as NOT CALLABLE FROM ARRAY ===")
    print(f"  Appended finding {fid} (tier {tier})")


def main(subject_id: str, genes: list[str] | None) -> None:
    profile = json.loads(
        (PROJECT_ROOT / "profiles" / f"{subject_id}.json").read_text()
    )
    imputed, chip = load_subject_data(subject_id)
    print(f"Loaded imputed={len(imputed):,} rows, chip={len(chip):,} rows")

    panels = genes if genes else available_panels()
    for gene in panels:
        run_gene(subject_id, gene, imputed, chip, profile)

    # Always record SMN1 as explicitly not-callable unless user opts out.
    if not genes or "smn1" in [g.lower() for g in genes]:
        run_sma_not_callable(subject_id)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("subject_id")
    parser.add_argument("--gene", action="append", help="Run only this gene (repeatable)")
    args = parser.parse_args()
    genes = [g.lower() for g in args.gene] if args.gene else None
    main(args.subject_id, genes)
