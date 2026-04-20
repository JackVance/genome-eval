"""Parser for MyHeritage raw data CSV files.

Expected format: quoted CSV with columns
    "RSID","CHROMOSOME","POSITION","RESULT"

Genotype in one 2-char concatenated field. Newer files (2022+) may use GRCh38 —
caller should check the header. If GRCh38, either liftover to 37 or track the
build separately in the profile metadata.
"""
from __future__ import annotations

import gzip
from pathlib import Path

import pandas as pd


def parse(path: Path) -> pd.DataFrame:
    opener = gzip.open if str(path).endswith(".gz") else open
    rows = []
    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
        header_seen = False
        for line in fh:
            if not line or line.startswith("#"):
                continue
            stripped = line.rstrip("\n").rstrip("\r")
            if not stripped:
                continue
            cells = [c.strip().strip('"') for c in stripped.split(",")]
            if not header_seen and cells and cells[0].lower() in ("rsid", "rs"):
                header_seen = True
                continue
            if len(cells) != 4:
                continue
            rsid, chrom, pos, geno = cells
            try:
                pos_i = int(pos)
            except ValueError:
                continue

            if not geno or geno in ("--", "00"):
                a1, a2 = None, None
            elif len(geno) == 1:
                c = geno if geno != "-" else None
                a1, a2 = c, c
            else:
                c1 = geno[0] if geno[0] != "-" else None
                c2 = geno[1] if geno[1] != "-" else None
                a1, a2 = c1, c2

            rows.append((rsid, chrom, pos_i, a1, a2))

    df = pd.DataFrame(rows, columns=["rsid", "chrom", "pos", "a1", "a2"])
    df["chrom"] = df["chrom"].astype(str)
    return df


def extract_header_meta(path: Path) -> dict:
    opener = gzip.open if str(path).endswith(".gz") else open
    meta = {}
    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
        head = "".join(fh.readline() for _ in range(30))
    if "GRCh38" in head or "hg38" in head.lower() or "build 38" in head.lower():
        meta["declared_build"] = "GRCh38"
    return meta
