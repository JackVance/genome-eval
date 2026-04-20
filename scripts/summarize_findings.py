"""Pretty-print the current active findings (post-supersede) for a subject.

Usage:
    python scripts/summarize_findings.py [subject_id]
"""
from __future__ import annotations

import io
import sys
from collections import Counter
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")

import ledger_io


def main(subject_id: str = "alice") -> None:
    findings = ledger_io.load_findings()
    active = ledger_io.load_active_findings(subject_id=subject_id, findings=findings)
    total_subject_rows = sum(1 for f in findings if f.get("subject_id") == subject_id)

    print(f"=== {len(active)} active findings for {subject_id} "
          f"({total_subject_rows} ledger rows for this subject, "
          f"{total_subject_rows - len(active)} superseded / tombstoned) ===\n")

    tier_counts = Counter(f.get("tier_computed", "?") for f in active)
    print("Tier distribution:", dict(tier_counts))
    print()

    tier_order = {"A": 0, "B": 1, "C": 2, "D": 3, "E": 4, "unknown": 5, "?": 6}
    active.sort(key=lambda f: (tier_order.get(f.get("tier_computed", "?"), 9),
                                f.get("topic", "")))

    for f in active:
        variants = ", ".join(
            f"{v.get('rsid')}="
            + (v.get("genotype") or ("not-on-chip" if not v.get("on_chip") else "no-call"))
            for v in f.get("variants", [])
        )
        tier = f.get("tier_computed", "?")
        downgrade = " (ancestry-downgrade)" if f.get("ancestry_downgrade") else ""
        print(f"[{tier}]{downgrade} {f['topic']}")
        print(f"       {f['claim']}")
        print(f"       variants: {variants}")
        print()


if __name__ == "__main__":
    subj = sys.argv[1] if len(sys.argv) > 1 else "alice"
    main(subj)
