"""Append superseding records for PRS findings whose numeric output is
suspect. Does NOT rerun the pipeline — just annotates the finding with an
explicit caveat block so future readers know the number is not load-bearing.

Current cases:
  - PGS000889 (LDL, Graham 2021 P+T): z = -7.5 sits far outside the
    empirical 1000G EUR range (46-55 raw score, SD 1.4). The weight file is
    marked `weight_type=NR` (not reported) — coefficients are on a
    non-standard scale. Calibration + scoring pipeline likely have a scale
    mismatch. Needs a different LDL PGS (e.g., the pure-EUR PGS000891 variant)
    before the number can be trusted.
  - PGS003971 (SBP, Shetty 2023 PRS-CS): only ~2% of the 1.1M panel variants
    matched the imputed parquet. The Shetty weights may use a coordinate
    system (hg38?) that differs from our GRCh37 reference despite the file
    name suggesting hg19. Needs investigation.
"""
from __future__ import annotations

import sys
from datetime import datetime, timezone
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

import ledger_io

SUSPECT = {
    "prs_ldl_pgs000889": {
        "reason": (
            "z-score of -7.5 is outside the empirical 1000G EUR distribution "
            "range (46.1 - 55.5 raw score, SD = 1.41). PGS000889 uses "
            "weight_type=NR (not reported) P+T coefficients on a "
            "non-standard scale; there's a scoring/calibration mismatch "
            "between scripts/calibrate_prs_empirical.py and "
            "scripts/run_prs.py that produces systematic divergence for "
            "this particular scoring scale."
        ),
        "followup": (
            "Try the pure-EUR Graham 2021 variant (PGS000891 candidate) or "
            "a different LDL PGS (e.g., Klarin 2018, Ripatti 2010 derivatives "
            "if available). Also audit the calibrate vs run_prs variant-"
            "matching logic for this specific weight-type case."
        ),
    },
    "prs_sbp_pgs003971": {
        "reason": (
            "Only ~2% of PGS003971's 1.1M variants contributed (22,379 out "
            "of 1,115,520), despite file metadata claiming hg19/GRCh37 "
            "harmonization. The z-score of +2.0 is therefore computed from "
            "a tiny subset of the intended panel and is not representative "
            "of the score's true design."
        ),
        "followup": (
            "Inspect PGS003971 weights against the imputed-parquet coordinates "
            "directly; chromosome prefix conventions (chr1 vs 1) or build "
            "mismatch are the most likely causes. Re-download the score via "
            "prs_download.py after confirming metadata. As fallback, try "
            "Evangelou 2018-derived PGP000283 scores once a specific ID is "
            "confirmed."
        ),
    },
}


def main() -> None:
    active = ledger_io.load_active_findings(subject_id="alice")
    topic_to_row = {r["topic"]: r for r in active}

    for topic, info in SUSPECT.items():
        old = topic_to_row.get(topic)
        if old is None:
            print(f"No active row for topic={topic!r}; skipping.")
            continue

        new = dict(old)
        new.pop("finding_id", None)
        new.pop("timestamp", None)
        new["supersedes"] = old["finding_id"]
        caveat = (
            f" \n\nSUSPECT NUMERIC OUTPUT — do not use as-is. "
            f"Reason: {info['reason']} "
            f"Follow-up: {info['followup']}"
        )
        new["notes"] = (new.get("notes") or "").rstrip() + caveat
        # Don't edit the claim — original numbers stay on the claim line so
        # a reader sees exactly what was computed; the notes explain why
        # it's not load-bearing.
        # Mark inference_confidence as low to reflect the suspicion.
        new["inference_confidence"] = "low"
        # Recompute tier with the updated inference_confidence.
        from tier_rules import compute_tier, TIER_RULE_VERSION
        tier, downgrade = compute_tier(new)
        new["tier_computed"] = tier
        new["tier_computed_at"] = datetime.now(timezone.utc).isoformat()
        new["tier_rule_version"] = TIER_RULE_VERSION
        new["ancestry_downgrade"] = downgrade
        fid = ledger_io.append_finding(**new)
        print(f"Flagged {topic} ({old['finding_id'][:8]} -> {fid[:8]}) tier={tier}")


if __name__ == "__main__":
    main()
