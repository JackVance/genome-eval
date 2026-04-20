"""One-off backfill: record 'not genotyped' findings for Tier 2 pathogenic
loci that the earlier investigation silently skipped because the variant
wasn't present on the chip.

Skill invariant: 'not on chip' ≠ 'wild-type'. Absence must be recorded as a
distinct finding so future queries against the ledger see the gap.

Idempotent: checks for an active finding on each topic with on_chip=False
before appending.
"""
from __future__ import annotations

import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

import ledger_io
from tier_rules import compute_tier, TIER_RULE_VERSION
from datetime import datetime, timezone

SUBJECT = "alice"

# Topic → (record template). Only include loci where ClinVar backs the
# evidence strength of the locus.
NOT_ON_CHIP_RECORDS = [
    {
        "topic": "prothrombin_g20210a",
        "rsid": "rs1799963",
        "gene": "F2",
        "chrom": "11",
        "pos": 46761055,
        "ref": "G",
        "alt": "A",
        "claim": "Prothrombin G20210A (rs1799963): NOT GENOTYPED — variant absent from this 23andMe v5 chip",
        "clinvar_significance": "Pathogenic",
        "clinvar_review_stars": 3,
        "cohort_ancestry": "European-dominant",
        "subject_ancestry_match": "match",
        "source_ids": ["clinvar:prothrombin-g20210a"],
        "notes": (
            "Locus is ClinVar Pathogenic (3-star) but is NOT on the chip. "
            "Absence of the variant is therefore unknown; the subject's status cannot be inferred. "
            "Carrier frequency ~2% in European populations. "
            "If clinically relevant (unexplained VTE, family history), order targeted testing."
        ),
    },
]


def active_topics_with_variant(findings, topic, rsid):
    """Return list of active (non-superseded) findings with this topic and rsid."""
    superseded = {f["supersedes"] for f in findings if f.get("supersedes")}
    out = []
    for f in findings:
        if f["finding_id"] in superseded:
            continue
        if f.get("topic") != topic:
            continue
        for v in f.get("variants", []):
            if v.get("rsid") == rsid:
                out.append(f)
                break
    return out


def main():
    findings = ledger_io.load_findings()
    added = []
    for rec in NOT_ON_CHIP_RECORDS:
        existing = active_topics_with_variant(findings, rec["topic"], rec["rsid"])
        if existing:
            print(f"skip {rec['topic']} / {rec['rsid']} — already has active finding(s)")
            continue

        variants = [{
            "rsid": rec["rsid"],
            "gene": rec["gene"],
            "chrom": rec["chrom"],
            "pos": rec["pos"],
            "ref": rec["ref"],
            "alt": rec["alt"],
            "genotype": None,
            "on_chip": False,
        }]

        kwargs = {
            "subject_id": SUBJECT,
            "topic": rec["topic"],
            "claim": rec["claim"],
            "variants": variants,
            "effect": None,
            "clinvar_significance": rec["clinvar_significance"],
            "clinvar_review_stars": rec["clinvar_review_stars"],
            "cohort_ancestry": rec["cohort_ancestry"],
            "subject_ancestry_match": rec["subject_ancestry_match"],
            "source_ids": rec["source_ids"],
            "notes": rec["notes"],
            "supersedes": None,
        }
        tier, downgrade = compute_tier(kwargs)
        kwargs["tier_computed"] = tier
        kwargs["tier_computed_at"] = datetime.now(timezone.utc).isoformat()
        kwargs["tier_rule_version"] = TIER_RULE_VERSION
        kwargs["ancestry_downgrade"] = downgrade

        fid = ledger_io.append_finding(**kwargs)
        added.append({"topic": rec["topic"], "rsid": rec["rsid"],
                      "tier": tier, "finding_id": fid})

    if not added:
        print("no new findings added")
    else:
        import json
        print(json.dumps(added, indent=2))


if __name__ == "__main__":
    main()
