"""One-time migration: backfill `supersede_chain_root` + `is_tombstone` on every
existing findings.jsonl row.

The append-only invariant of the ledger is about *content* — you don't revise a
finding by editing an existing row, you append a superseder. Adding derived
metadata fields (chain_root computed from the existing supersedes pointers,
is_tombstone inferred from the claim prefix) is not content revision, so an
in-place rewrite is acceptable. The script makes a timestamped backup first.

After this migration, `ledger_io.load_active_findings()` is the canonical way to
read the ledger; callers no longer need to build their own superseded_ids set.

Usage:
    python scripts/migrate_ledger_schema.py [--dry-run]
"""
from __future__ import annotations

import json
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
sys.path.insert(0, str(SCRIPT_DIR))

import ledger_io

LEDGER_PATH = PROJECT_ROOT / "ledger" / "findings.jsonl"


def main(dry_run: bool = False) -> int:
    rows = ledger_io.load_findings()
    if not rows:
        print("No findings to migrate.")
        return 0

    rows_by_id = {r["finding_id"]: r for r in rows}

    migrated = []
    changed_chain_root = 0
    changed_is_tombstone = 0
    for r in rows:
        new = dict(r)
        if "supersede_chain_root" not in new or not new.get("supersede_chain_root"):
            new["supersede_chain_root"] = ledger_io.compute_chain_root(new, rows_by_id)
            changed_chain_root += 1
        if "is_tombstone" not in new:
            new["is_tombstone"] = ledger_io._infer_tombstone(new)
            changed_is_tombstone += 1
        migrated.append(new)

    print(f"Rows total:              {len(rows)}")
    print(f"chain_root backfilled:   {changed_chain_root}")
    print(f"is_tombstone backfilled: {changed_is_tombstone}")
    tombstones = sum(1 for r in migrated if r.get("is_tombstone"))
    print(f"Tombstones (flagged):    {tombstones}")

    # Sanity cross-check: after migration, load_active_findings should return
    # the same count as the legacy (naive-filter + startswith-SUPERSEDED)
    # pre-migration path.
    legacy_superseded = {r["supersedes"] for r in rows if r.get("supersedes")}
    legacy_active = [
        r for r in rows
        if r["finding_id"] not in legacy_superseded
        and not (r.get("claim") or "").startswith("SUPERSEDED")
    ]
    # Count by subject for both paths.
    from collections import Counter
    legacy_by_subj = Counter(r.get("subject_id", "?") for r in legacy_active)
    new_active = ledger_io.load_active_findings(findings=migrated)
    new_by_subj = Counter(r.get("subject_id", "?") for r in new_active)
    print(f"Legacy-filter active (by subject): {dict(legacy_by_subj)}")
    print(f"New-filter active    (by subject): {dict(new_by_subj)}")
    if legacy_by_subj != new_by_subj:
        print("WARNING: new-filter active count differs from legacy-filter.")
        # Show the differences.
        legacy_ids = {r["finding_id"] for r in legacy_active}
        new_ids = {r["finding_id"] for r in new_active}
        only_in_legacy = legacy_ids - new_ids
        only_in_new = new_ids - legacy_ids
        if only_in_legacy:
            print(f"  In legacy but not new ({len(only_in_legacy)}):")
            for fid in list(only_in_legacy)[:10]:
                r = rows_by_id[fid]
                print(f"    {fid[:8]} {r.get('topic')} | {(r.get('claim') or '')[:80]}")
        if only_in_new:
            print(f"  In new but not legacy ({len(only_in_new)}):")
            for fid in list(only_in_new)[:10]:
                r = rows_by_id[fid]
                print(f"    {fid[:8]} {r.get('topic')} | {(r.get('claim') or '')[:80]}")

    if dry_run:
        print("Dry run — no file written.")
        return 0

    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    backup = LEDGER_PATH.with_name(f"findings.jsonl.bak.{stamp}")
    shutil.copy2(LEDGER_PATH, backup)
    print(f"Backup written: {backup.name}")

    with LEDGER_PATH.open("w", encoding="utf-8") as fh:
        for r in migrated:
            fh.write(json.dumps(r) + "\n")
    print(f"Rewrote {LEDGER_PATH.name} with {len(migrated)} rows.")

    # Final verification: re-load and count.
    final = ledger_io.load_active_findings()
    print(f"Post-migration active count: {len(final)}")
    return 0


if __name__ == "__main__":
    dry = "--dry-run" in sys.argv
    sys.exit(main(dry_run=dry))
