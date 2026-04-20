"""Append-only helpers for the ledger JSONL files.

Never edits existing rows. To revise a finding, append a new row with
`supersedes` set to the superseded finding_id.

Schema fields (findings.jsonl):
  - finding_id: stable UUID for this row.
  - timestamp: ISO-8601 write time.
  - supersedes: UUID of the parent row this one replaces, or None.
  - supersede_chain_root: UUID of the earliest row in this chain. Equals
      finding_id for a brand-new root. Inherited from the parent otherwise.
      Computed at write time so readers never have to walk the chain.
  - is_tombstone: bool. True iff the row exists only to announce that an
      older row has been superseded and carries no new content. The
      append-only model has always permitted content-carrying superseders
      (which is what you almost always want); is_tombstone is reserved for
      the rare case where the actual canonical lives elsewhere and this row
      is pure bookkeeping. A tombstone chain whose head is a tombstone is
      considered inactive — the chain is dead.
  - (plus all the finding fields: subject_id, topic, claim, variants, effect,
    notes, tiering metadata, etc.)
"""
from __future__ import annotations

import json
import uuid
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Iterator

PROJECT_ROOT = Path(__file__).resolve().parent.parent
LEDGER_DIR = PROJECT_ROOT / "ledger"


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _append(path: Path, record: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "a", encoding="utf-8") as fh:
        fh.write(json.dumps(record) + "\n")


def _infer_tombstone(record: dict[str, Any]) -> bool:
    """Best-effort inference of tombstone status for legacy rows.

    True when the row is explicitly flagged, or its claim begins with
    "SUPERSEDED" — the pattern used by older one-off supersede scripts.
    Returning the heuristic here (rather than at read sites) means every
    consumer sees the same definition.
    """
    if record.get("is_tombstone") is True:
        return True
    claim = record.get("claim") or ""
    return claim.lstrip().startswith("SUPERSEDED")


def compute_chain_root(
    record: dict[str, Any],
    rows_by_id: dict[str, dict[str, Any]],
) -> str:
    """Walk `supersedes` pointers back to the earliest row in the chain.

    Returns the finding_id of the oldest ancestor. If `record` has no parent
    (or the parent is missing from rows_by_id), the chain root is the
    record's own finding_id. A cycle guard prevents infinite loops if the
    ledger ever contains a malformed chain.
    """
    seen: set[str] = set()
    cur = record
    while True:
        fid = cur.get("finding_id")
        if fid is None:
            # Record doesn't have an id yet (being assembled). Walking is
            # still possible via `supersedes`.
            pass
        elif fid in seen:
            # Cycle detected — bail out at the current cur.
            break
        else:
            seen.add(fid)

        parent_id = cur.get("supersedes")
        if not parent_id:
            return cur["finding_id"]
        parent = rows_by_id.get(parent_id)
        if parent is None:
            # Orphan chain — parent was never loaded. Treat cur as the root.
            return cur["finding_id"]
        cur = parent


def append_finding(**kwargs) -> str:
    """Append a finding row.

    Auto-computes `supersede_chain_root` and `is_tombstone` when absent.
    Callers may pass them explicitly to skip the derivation (e.g., during
    bulk migration where the ledger is held in memory).
    """
    rec: dict[str, Any] = {
        "finding_id": str(uuid.uuid4()),
        "timestamp": _now(),
        "supersedes": None,
        **kwargs,
    }
    rec.setdefault("is_tombstone", _infer_tombstone(rec))
    if "supersede_chain_root" not in rec:
        parent_id = rec.get("supersedes")
        if not parent_id:
            rec["supersede_chain_root"] = rec["finding_id"]
        else:
            # Lookup parent in the current ledger to inherit its chain_root.
            existing = load_findings()
            rows_by_id = {r["finding_id"]: r for r in existing}
            parent = rows_by_id.get(parent_id)
            if parent and parent.get("supersede_chain_root"):
                rec["supersede_chain_root"] = parent["supersede_chain_root"]
            else:
                # Fall back to walking the chain. For a legacy ledger where
                # chain_root hasn't been backfilled, this still yields the
                # correct root.
                rec["supersede_chain_root"] = compute_chain_root(rec, rows_by_id)
    _append(LEDGER_DIR / "findings.jsonl", rec)
    return rec["finding_id"]


def append_source(**kwargs) -> str:
    source_id = kwargs.get("source_id") or str(uuid.uuid4())
    rec = {
        "source_id": source_id,
        "accessed_at": _now(),
        **kwargs,
    }
    rec["source_id"] = source_id
    _append(LEDGER_DIR / "sources.jsonl", rec)
    return source_id


def append_investigation(**kwargs) -> str:
    rec = {
        "investigation_id": str(uuid.uuid4()),
        "timestamp": _now(),
        **kwargs,
    }
    _append(LEDGER_DIR / "investigations.jsonl", rec)
    return rec["investigation_id"]


def read_jsonl(path: Path) -> Iterator[dict[str, Any]]:
    if not path.exists():
        return
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            yield json.loads(line)


def load_findings() -> list[dict[str, Any]]:
    return list(read_jsonl(LEDGER_DIR / "findings.jsonl"))


def load_sources() -> list[dict[str, Any]]:
    return list(read_jsonl(LEDGER_DIR / "sources.jsonl"))


def load_investigations() -> list[dict[str, Any]]:
    return list(read_jsonl(LEDGER_DIR / "investigations.jsonl"))


def _chain_root_for(
    record: dict[str, Any],
    rows_by_id: dict[str, dict[str, Any]],
) -> str:
    """Return the stored chain_root if present, else compute on the fly.

    Keeps this logic local so both `load_active_findings` and the migration
    tool can share a definition.
    """
    stored = record.get("supersede_chain_root")
    if stored:
        return stored
    return compute_chain_root(record, rows_by_id)


def load_active_findings(
    subject_id: str | None = None,
    findings: Iterable[dict[str, Any]] | None = None,
) -> list[dict[str, Any]]:
    """Return the current active findings per the chain_root + tombstone model.

    Grouping: one group per `supersede_chain_root`. Within each group, the
    "head" is the latest row by timestamp. A group is active iff its head is
    not a tombstone — otherwise the chain is dead (no live canonical). The
    result is one row per active chain.

    When `subject_id` is set, the filter is applied before grouping so chains
    belonging to other subjects are never considered.

    This is the one correct way to read "active findings". Consumers that
    previously rolled their own `superseded_ids = {...}` filter should switch
    to this helper — it handles tombstones and supersede-of-supersede chains
    in one place.
    """
    rows = list(findings) if findings is not None else load_findings()
    if subject_id is not None:
        rows = [r for r in rows if r.get("subject_id") == subject_id]
    rows_by_id = {r["finding_id"]: r for r in rows}

    by_chain: dict[str, list[dict[str, Any]]] = {}
    for r in rows:
        root = _chain_root_for(r, rows_by_id)
        by_chain.setdefault(root, []).append(r)

    active: list[dict[str, Any]] = []
    for chain in by_chain.values():
        # Latest by timestamp (fall back to insertion order if ties / missing).
        chain_sorted = sorted(chain, key=lambda r: r.get("timestamp") or "")
        head = chain_sorted[-1]
        if _infer_tombstone(head):
            continue
        active.append(head)
    return active
