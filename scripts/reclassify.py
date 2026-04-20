"""Re-derive subject_ancestry_match + tier across the ledger and append
superseding records only where the tier changes.

Honest-history approach: tiers are derived, but we preserve the history by
writing new rows with `supersedes: <old_id>` rather than editing.
"""
from __future__ import annotations

import json
import sys
import uuid
from datetime import datetime, timezone
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
sys.path.insert(0, str(SCRIPT_DIR))

import ledger_io
from tier_rules import compute_tier, TIER_RULE_VERSION


def match_ancestry(declared: dict | None, cohort: str | None) -> str:
    """Return match / partial / mismatch / unknown."""
    if not declared or not cohort:
        return "unknown"
    top = (declared.get("top_level") or "").lower()
    c = cohort.lower()

    # Explicit top-level match.
    if top and (top in c or c in top):
        return "match"
    # "European-dominant" / "european" / derivatives.
    if top == "european" and (
        "european" in c or c == "eur" or c in ("nfe", "nwe")
    ):
        return "match"
    # Broadly inclusive cohorts — treat as match (imperfect but not a
    # mismatch).
    if c in ("multi-ethnic", "multiethnic", "global", "n/a", ""):
        return "match"
    # Otherwise assume mismatch.
    return "mismatch"


def main(subject_id: str = "alice") -> None:
    profile_path = PROJECT_ROOT / "profiles" / f"{subject_id}.json"
    profile = json.loads(profile_path.read_text(encoding="utf-8"))
    declared = profile.get("declared_ancestry")

    # Use the canonical active-filter from ledger_io — chain_root + tombstone
    # aware, so we never try to refresh metadata on a dead chain.
    active = ledger_io.load_active_findings(subject_id=subject_id)

    changes = []
    for f in active:
        new_match = match_ancestry(declared, f.get("cohort_ancestry"))
        candidate = dict(f)
        candidate["subject_ancestry_match"] = new_match
        new_tier, downgrade = compute_tier(candidate)

        tier_changed = new_tier != f.get("tier_computed")
        match_changed = new_match != f.get("subject_ancestry_match")
        downgrade_changed = downgrade != f.get("ancestry_downgrade", False)
        version_changed = TIER_RULE_VERSION != f.get("tier_rule_version")

        if not (tier_changed or match_changed or downgrade_changed or version_changed):
            continue

        new_rec = dict(f)
        new_rec["finding_id"] = str(uuid.uuid4())
        new_rec["timestamp"] = datetime.now(timezone.utc).isoformat()
        new_rec["supersedes"] = f["finding_id"]
        new_rec["subject_ancestry_match"] = new_match
        new_rec["tier_computed"] = new_tier
        new_rec["tier_computed_at"] = new_rec["timestamp"]
        new_rec["tier_rule_version"] = TIER_RULE_VERSION
        new_rec["ancestry_downgrade"] = downgrade

        # Compose a concise explanation of what moved. A metadata-only refresh
        # (tier unchanged) deserves a different line than a tier movement, so
        # a reader can tell at a glance why the superseder exists.
        parts = []
        if tier_changed:
            parts.append(
                f"tier {f.get('tier_computed')} → {new_tier}"
            )
        if match_changed:
            parts.append(
                f"ancestry_match {f.get('subject_ancestry_match')!r} → {new_match!r}"
            )
        if downgrade_changed:
            parts.append(
                f"ancestry_downgrade {f.get('ancestry_downgrade', False)} → {downgrade}"
            )
        if version_changed and not (tier_changed or match_changed or downgrade_changed):
            parts.append(
                f"tier_rule_version {f.get('tier_rule_version')!r} → {TIER_RULE_VERSION!r}"
            )
        add = (
            f"Metadata refresh: {'; '.join(parts)} "
            f"(declared_ancestry={declared.get('top_level')!r}, "
            f"cohort={f.get('cohort_ancestry')!r})."
        )
        note = (new_rec.get("notes") or "").strip()
        new_rec["notes"] = (note + " " + add).strip() if note else add

        ledger_io.append_finding(**{k: v for k, v in new_rec.items()
                                    if k not in ("finding_id", "timestamp")})
        changes.append({
            "topic": f["topic"],
            "from_tier": f.get("tier_computed"),
            "to_tier": new_tier,
            "ancestry_match_was": f.get("subject_ancestry_match"),
            "ancestry_match_now": new_match,
            "downgrade_was": f.get("ancestry_downgrade", False),
            "downgrade_now": downgrade,
            "superseded_id": f["finding_id"],
            "kind": "tier_movement" if tier_changed else "metadata_refresh",
        })

    return changes


if __name__ == "__main__":
    subj = sys.argv[1] if len(sys.argv) > 1 else "alice"
    changes = main(subj)
    print(json.dumps(changes, indent=2))
