"""Normalize a raw genotype file to standardized parquet + profile JSON.

Usage:
    python scripts/normalize.py <subject_id> <raw_file_path> [--display-name NAME]

Detects provider from header, dispatches to the correct parser, writes:
  - standardized-genomes/<id>.parquet
  - profiles/<id>.json

Analysis code must never read from raw-source-genomes/ — only from these outputs.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
sys.path.insert(0, str(SCRIPT_DIR))

import parse_23andme
import parse_ancestry
import parse_myheritage
import parse_ftdna

PARSERS = {
    "23andMe": (parse_23andme.parse, parse_23andme.extract_header_meta),
    "AncestryDNA": (parse_ancestry.parse, parse_ancestry.extract_header_meta),
    "MyHeritage": (parse_myheritage.parse, parse_myheritage.extract_header_meta),
    "FamilyTreeDNA": (parse_ftdna.parse, parse_ftdna.extract_header_meta),
}


def detect_provider(path: Path) -> tuple[str, float]:
    """Return (provider, confidence 0..1). Raises on failure."""
    with open(path, "rt", encoding="utf-8", errors="replace") as fh:
        head = "".join(fh.readline() for _ in range(50))
    lower = head.lower()
    if "23andme" in lower:
        return "23andMe", 1.0
    if "ancestrydna" in lower or "ancestry.com" in lower:
        return "AncestryDNA", 1.0
    if "myheritage" in lower:
        return "MyHeritage", 1.0
    if "familytreedna" in lower or "ftdna" in lower:
        return "FamilyTreeDNA", 1.0

    # Structural fallback.
    data_lines = [
        ln for ln in head.splitlines()
        if ln and not ln.startswith("#") and not ln.lower().startswith("rsid")
    ]
    if data_lines:
        sample = data_lines[0]
        tab_cols = len(sample.split("\t"))
        comma_cols = len(sample.split(","))
        if tab_cols == 4:
            return "23andMe", 0.5
        if tab_cols == 5:
            return "AncestryDNA", 0.5
        if comma_cols == 4:
            return "MyHeritage", 0.4

    raise ValueError(
        f"Could not detect provider. First 500 chars:\n{head[:500]}"
    )


def detect_chip_version_23andme(n_snps: int) -> str:
    if n_snps > 900_000:
        return "v3 (~960k)"
    if 700_000 <= n_snps <= 900_000:
        return "v4 or v5 (~600-700k)"
    if 500_000 <= n_snps < 700_000:
        return "v5 (~640k) / v4 (~600k)"
    return f"unknown (N={n_snps})"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def infer_sex(df: pd.DataFrame) -> str:
    y_calls = df[(df["chrom"] == "Y") & df["a1"].notna()]
    y_total = df[df["chrom"] == "Y"]
    if len(y_total) == 0:
        return "unknown (no Y markers on chip)"
    call_rate = len(y_calls) / len(y_total) if len(y_total) else 0
    if call_rate > 0.5:
        return "male (Y call rate high)"
    if call_rate < 0.05:
        return "female (Y effectively no-call)"
    return f"ambiguous (Y call rate {call_rate:.1%})"


def parse_stats(df: pd.DataFrame) -> dict:
    total = len(df)
    no_call = int(((df["a1"].isna()) & (df["a2"].isna())).sum())
    chroms = sorted(df["chrom"].unique().tolist())
    has_mt = "MT" in chroms
    has_x = "X" in chroms
    has_y = "Y" in chroms
    return {
        "total_snps": total,
        "no_call_count": no_call,
        "no_call_rate": round(no_call / total, 6) if total else 0,
        "chromosomes_present": chroms,
        "has_mitochondrial": has_mt,
        "has_x": has_x,
        "has_y": has_y,
        "inferred_sex": infer_sex(df),
    }


def ingest(subject_id: str, raw_path: Path, display_name: str | None = None) -> dict:
    provider, confidence = detect_provider(raw_path)
    parser_fn, header_fn = PARSERS[provider]

    df = parser_fn(raw_path)

    out_parquet = PROJECT_ROOT / "standardized-genomes" / f"{subject_id}.parquet"
    out_parquet.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(out_parquet, index=False)

    stats = parse_stats(df)
    header_meta = header_fn(raw_path)

    declared_build = header_meta.get("declared_build", "GRCh37")

    chip_version = None
    if provider == "23andMe":
        chip_version = detect_chip_version_23andme(stats["total_snps"])

    profile = {
        "subject_id": subject_id,
        "display_name": display_name or subject_id,
        "provider": provider,
        "provider_detection_confidence": confidence,
        "chip_version": chip_version,
        "build": declared_build,
        "raw_file": {
            "path": str(raw_path.relative_to(PROJECT_ROOT)).replace("\\", "/"),
            "sha256": sha256_file(raw_path),
            "size_bytes": raw_path.stat().st_size,
            "provider_header": header_meta,
        },
        "standardized_file": {
            "path": str(out_parquet.relative_to(PROJECT_ROOT)).replace("\\", "/"),
        },
        "parse_stats": stats,
        "declared_ancestry": None,
        "current_medications": None,
        "family_history": None,
        "sharing_sensitivity": "unset",
        "ingested_at": datetime.now(timezone.utc).isoformat(),
    }

    profile_path = PROJECT_ROOT / "profiles" / f"{subject_id}.json"
    profile_path.parent.mkdir(parents=True, exist_ok=True)
    profile_path.write_text(json.dumps(profile, indent=2), encoding="utf-8")

    return profile


def ensure_ledger():
    ledger_dir = PROJECT_ROOT / "ledger"
    ledger_dir.mkdir(parents=True, exist_ok=True)
    for name in ("findings.jsonl", "sources.jsonl", "investigations.jsonl"):
        p = ledger_dir / name
        if not p.exists():
            p.write_text("", encoding="utf-8")


def main():
    ap = argparse.ArgumentParser(description="Normalize a raw genotype file.")
    ap.add_argument("subject_id", help="Short filename-safe ID, e.g. 'alice'")
    ap.add_argument("raw_path", help="Path to raw provider file")
    ap.add_argument("--display-name", default=None)
    args = ap.parse_args()

    raw = Path(args.raw_path).resolve()
    if not raw.exists():
        raise SystemExit(f"Raw file not found: {raw}")

    ensure_ledger()
    profile = ingest(args.subject_id, raw, args.display_name)
    print(json.dumps(profile, indent=2))


if __name__ == "__main__":
    main()
