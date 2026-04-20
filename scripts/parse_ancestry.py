"""Parser for AncestryDNA raw data TSV files.

Expected format: `#`-prefixed header comments, optional `rsid...` header row,
then 5 tab-separated columns:
    rsid \t chromosome \t position \t allele1 \t allele2

Alleles are in separate columns. No-call is `0` (not `--`).
Chromosome encoding uses 23/24/25/26 for X/Y/PAR/MT — remapped on parse.
Build: GRCh37.
"""
from __future__ import annotations

import gzip
from pathlib import Path

import pandas as pd

_CHROM_MAP = {"23": "X", "24": "Y", "25": "X", "26": "MT"}


def parse(path: Path) -> pd.DataFrame:
    opener = gzip.open if str(path).endswith(".gz") else open
    rows = []
    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            line = line.rstrip("\n").rstrip("\r")
            if not line:
                continue
            if line.lower().startswith("rsid"):
                continue
            parts = line.split("\t")
            if len(parts) != 5:
                continue
            rsid, chrom, pos, a1, a2 = parts
            try:
                pos_i = int(pos)
            except ValueError:
                continue
            chrom = _CHROM_MAP.get(chrom, chrom)
            a1 = None if a1 in ("0", "", "-") else a1
            a2 = None if a2 in ("0", "", "-") else a2
            rows.append((rsid, chrom, pos_i, a1, a2))

    df = pd.DataFrame(rows, columns=["rsid", "chrom", "pos", "a1", "a2"])
    df["chrom"] = df["chrom"].astype(str)
    return df


def extract_header_meta(path: Path) -> dict:
    return {}
