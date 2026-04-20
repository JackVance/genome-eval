"""Derived tier classification for findings.

Tiers are computed from finding metrics, never assigned directly.
Bump TIER_RULE_VERSION if the rules change; old tiers can then be recomputed
via scripts/reclassify.py or topic-specific recomputation.

v2 (2026-04-17):
  - Added evidence_class-driven tiering for well-established trait findings
    that lack formal GWAS N/p metrics (Mendelian traits, well-replicated
    common variants, weakly predictive variants).
  - Added inference_confidence handling for gene-presence / multi-SNP
    composite / other derived findings where the finding's confidence is
    itself the headline evidence axis.
  - Accept both older and newer field names (p_value / pvalue, study_n /
    n_cases+n_controls) so earlier findings remain classifiable.

Tier captures evidence quality for the general claim (how well-established
the genotype-phenotype association is). It does NOT capture penetrance,
genotype-call reliability, or per-subject confidence — those are surfaced
separately per Rule 10.
"""
from __future__ import annotations

from typing import Any

TIER_RULE_VERSION = "v2"


def _get_p(finding: dict[str, Any]):
    for k in ("p_value", "pvalue"):
        v = finding.get(k)
        if v is not None:
            return v
    return None


def _get_n(finding: dict[str, Any]):
    """Return a best-effort study N — legacy study_n field, or n_cases+n_controls."""
    n = finding.get("study_n")
    if n is not None:
        return n
    nc = finding.get("n_cases")
    nctl = finding.get("n_controls")
    if nc is not None and nctl is not None:
        return nc + nctl
    if nc is not None:
        return nc
    return None


def compute_tier(finding: dict[str, Any]) -> tuple[str, bool]:
    """Return (tier, ancestry_downgrade).

    Tier is one of: A, B, C, D, E, unknown.
    ancestry_downgrade is True when we dropped a tier because the subject's
    declared ancestry doesn't match the source cohort ancestry.
    """
    evidence_class = (finding.get("evidence_class") or "").lower()
    clinvar_stars = finding.get("clinvar_review_stars")
    clinvar_sig = (finding.get("clinvar_significance") or "").lower()
    cpic_level = (finding.get("cpic_level") or "").upper()
    inference_conf = (finding.get("inference_confidence") or "").lower()

    n = _get_n(finding)
    p = _get_p(finding)
    rep = finding.get("replication_count")
    ancestry_match = (finding.get("subject_ancestry_match") or "unknown").lower()

    # --- Tier A: clinical-grade evidence ---
    if cpic_level in ("A", "B"):
        return _apply_ancestry("A", ancestry_match)
    if (
        clinvar_sig in ("pathogenic", "likely pathogenic", "pathogenic/likely pathogenic")
        and isinstance(clinvar_stars, (int, float))
        and clinvar_stars >= 2
    ):
        return _apply_ancestry("A", ancestry_match)
    # Near-Mendelian / classical single-gene traits with strong replication.
    if evidence_class in ("mendelian_trait", "near_mendelian_trait"):
        if isinstance(rep, (int, float)) and rep >= 5:
            return _apply_ancestry("A", ancestry_match)
        # Even without replication_count, Mendelian classification implies A.
        return _apply_ancestry("A", ancestry_match)

    # --- Tier B: well-replicated common-variant trait, OR high-confidence inference,
    #     OR GWAS-significant replicated ancestry-matched study ---
    if evidence_class == "well_replicated_common_variant":
        if isinstance(rep, (int, float)) and rep >= 10:
            return _apply_ancestry("B", ancestry_match)
        if isinstance(rep, (int, float)) and rep >= 3:
            return _apply_ancestry("C", ancestry_match)
        # Insufficient replication count even if labeled this way.
        return _apply_ancestry("D", ancestry_match)

    if evidence_class in ("multi_snp_composite", "gene_presence_inference"):
        if inference_conf == "high":
            return _apply_ancestry("B", ancestry_match)
        if inference_conf == "moderate":
            return _apply_ancestry("C", ancestry_match)
        if inference_conf == "low":
            return _apply_ancestry("D", ancestry_match)
        # inference_class set but no confidence → unknown is honest.
        return "unknown", False

    if evidence_class == "weakly_predictive_variant":
        return _apply_ancestry("D", ancestry_match)

    if evidence_class in ("suspected_miscall", "array_limitation", "not_callable_from_array"):
        # Findings that explicitly acknowledge the chip can't give a reliable answer.
        # Tier "unknown" is appropriate — evidence quality is not the issue; calling reliability is.
        return "unknown", False

    # GWAS-based rules (need N and p).
    if n is None or p is None:
        # Nothing left to key on except replication count for a rough fallback.
        if isinstance(rep, (int, float)) and rep >= 10:
            return _apply_ancestry("C", ancestry_match)
        if isinstance(rep, (int, float)) and rep >= 3:
            return _apply_ancestry("D", ancestry_match)
        return "unknown", False

    # Tier B: replicated, large, GWAS-significant, ancestry-matched.
    if (
        isinstance(rep, (int, float)) and rep >= 2
        and n >= 10_000
        and p < 5e-8
        and ancestry_match == "match"
    ):
        return "B", False

    # Tier C: single large GWAS OR multiple small consistent studies.
    if n >= 5_000 and p < 5e-8:
        return _apply_ancestry("C", ancestry_match)
    if (
        isinstance(rep, (int, float)) and rep >= 2
        and 500 <= n <= 5_000
    ):
        return _apply_ancestry("C", ancestry_match)

    # Tier D: nominal significance or candidate-gene.
    if p < 0.05:
        return _apply_ancestry("D", ancestry_match)

    # Tier E: weak.
    if evidence_class in ("snpedia", "press", "community"):
        return "E", False
    if n < 500:
        return "E", False

    return "unknown", False


def _apply_ancestry(tier: str, ancestry_match: str) -> tuple[str, bool]:
    if ancestry_match == "mismatch":
        order = ["A", "B", "C", "D", "E"]
        if tier in order:
            idx = order.index(tier)
            downgraded = order[min(idx + 1, len(order) - 1)]
            return downgraded, True
    return tier, False


if __name__ == "__main__":
    # Smoke tests covering old + new paths.
    examples = [
        ("CPIC A drug", {"cpic_level": "A", "subject_ancestry_match": "match"}),
        ("ClinVar pathogenic 3-star", {"clinvar_significance": "Pathogenic",
                                        "clinvar_review_stars": 3,
                                        "subject_ancestry_match": "match"}),
        ("Mendelian trait, well-replicated", {"evidence_class": "mendelian_trait",
                                               "replication_count": 20,
                                               "subject_ancestry_match": "match"}),
        ("Well-replicated common variant (eye color)", {
            "evidence_class": "well_replicated_common_variant",
            "replication_count": 15,
            "subject_ancestry_match": "match",
        }),
        ("Weakly predictive (CYP1A2 caffeine)", {
            "evidence_class": "weakly_predictive_variant",
            "subject_ancestry_match": "match",
        }),
        ("Gene-presence inference, high confidence (Rh+)", {
            "evidence_class": "gene_presence_inference",
            "inference_confidence": "high",
            "subject_ancestry_match": "match",
        }),
        ("Multi-SNP composite, high confidence (eye color)", {
            "evidence_class": "multi_snp_composite",
            "inference_confidence": "high",
            "subject_ancestry_match": "match",
        }),
        ("Not callable from array (CYP2D6)", {
            "evidence_class": "not_callable_from_array",
            "subject_ancestry_match": "match",
        }),
        ("GWAS replicated large matched", {
            "n_cases": 30_000, "n_controls": 30_000, "pvalue": 1e-20,
            "replication_count": 5, "subject_ancestry_match": "match",
        }),
    ]
    for label, ex in examples:
        tier, dg = compute_tier(ex)
        print(f"  {label:60s} -> tier={tier}  downgrade={dg}")
