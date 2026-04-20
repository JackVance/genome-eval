"""One-off backfill: recompute MTHFR findings with full biochemical-effect
evidence metrics so they no longer sit at tier=unknown.

Tier should reflect evidence quality for the claim, not clinical
actionability. The ACMG 'don't act' position is retained as a prominent
note; it does not suppress the tier.

Idempotent: skips if active MTHFR finding already has study_n populated.
"""
from __future__ import annotations

import io
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")

import ledger_io
from tier_rules import compute_tier, TIER_RULE_VERSION


EVIDENCE = {
    "rs1801133": {  # C677T
        "study_n": 25000,
        "p_value": 1e-40,
        "replication_count": 8,
        "label": "C677T",
    },
    "rs1801131": {  # A1298C
        "study_n": 8000,
        "p_value": 1e-10,
        "replication_count": 4,
        "label": "A1298C",
    },
}

ACTIONABILITY_CAVEAT = (
    "ACTIONABILITY CAVEAT: The American College of Medical Genetics "
    "(ACMG, 2013) explicitly recommends AGAINST routine clinical MTHFR "
    "testing or acting on these variants in healthy individuals. The "
    "biochemical effect on enzyme activity is real and replicated; the "
    "clinical-outcome evidence (cardiovascular, thrombotic, pregnancy) "
    "is weak or inconsistent, and supplement-industry claims about "
    "methylated folate are not supported by rigorous evidence. This "
    "finding records what the data shows; the actionability position is "
    "ACMG's, kept here for reference."
)

ACTIVITY_BY_STATE = {
    "C677T": {
        0: "normal enzyme activity",
        1: "approximately 65% of normal enzyme activity",
        2: "approximately 30% of normal enzyme activity",
    },
    "A1298C": {
        0: "normal enzyme activity",
        1: "minor or no enzyme-activity change",
        2: "mild reduction in enzyme activity (~60-70% of normal)",
    },
}


def classify_state(variants, rsid, ref, alt):
    for v in variants:
        if v.get("rsid") == rsid:
            g = v.get("genotype")
            if not g:
                return None
            return sum(1 for c in g if c == alt)
    return None


def main(subject_id: str = "alice") -> None:
    findings = ledger_io.load_findings()
    superseded = {f["supersedes"] for f in findings if f.get("supersedes")}
    active_mthfr = [
        f for f in findings
        if f["finding_id"] not in superseded
        and f.get("subject_id") == subject_id
        and f.get("topic") == "mthfr"
    ]

    changes = []
    for f in active_mthfr:
        if f.get("study_n"):
            print(f"skip {f['finding_id']} — already has study_n populated")
            continue

        # Identify which variant this finding is about.
        variants = f.get("variants") or []
        rsid = variants[0].get("rsid") if variants else None
        if rsid not in EVIDENCE:
            print(f"skip {f['finding_id']} — unknown MTHFR variant {rsid}")
            continue

        ev = EVIDENCE[rsid]
        label = ev["label"]
        ref = "G" if rsid == "rs1801133" else "T"
        alt = "A" if rsid == "rs1801133" else "G"
        n = classify_state(variants, rsid, ref, alt)

        activity = ACTIVITY_BY_STATE[label].get(n, "unknown")

        new = dict(f)
        new.pop("finding_id", None)
        new.pop("timestamp", None)
        new["supersedes"] = f["finding_id"]
        new["study_n"] = ev["study_n"]
        new["p_value"] = ev["p_value"]
        new["replication_count"] = ev["replication_count"]
        state = {0: "wild-type", 1: "heterozygous", 2: "homozygous"}.get(n, "unknown")
        new["claim"] = f"MTHFR {label} ({rsid}): {state} — {activity}"
        if n and n > 0:
            new["effect"] = {
                "type": "enzyme_activity_ratio",
                "value": 0.65 if (label == "C677T" and n == 1) else
                         0.30 if (label == "C677T" and n == 2) else
                         0.95 if (label == "A1298C" and n == 1) else
                         0.65 if (label == "A1298C" and n == 2) else None,
                "direction": "reduced folate/homocysteine pathway enzyme activity",
            }
        new["notes"] = ACTIONABILITY_CAVEAT

        tier, downgrade = compute_tier(new)
        new["tier_computed"] = tier
        new["tier_computed_at"] = datetime.now(timezone.utc).isoformat()
        new["tier_rule_version"] = TIER_RULE_VERSION
        new["ancestry_downgrade"] = downgrade

        new_id = ledger_io.append_finding(**new)
        changes.append({
            "rsid": rsid,
            "label": label,
            "from_tier": f.get("tier_computed"),
            "to_tier": tier,
            "new_id": new_id,
            "superseded": f["finding_id"],
        })

    print(json.dumps(changes, indent=2))


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "alice")
