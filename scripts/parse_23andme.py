"""Parser for 23andMe raw data TSV files.

Expected format: `#`-prefixed header comments, then 4 columns:
    rsid \t chromosome \t position \t genotype

Genotype is a 2-char concatenated string (AG, CC), `--` for no-call,
single char for X/Y/MT hemizygous positions in males.
"""
from __future__ import annotations

import gzip
from pathlib import Path

import pandas as pd


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
            parts = line.split("\t")
            if len(parts) != 4:
                continue
            rsid, chrom, pos, geno = parts
            try:
                pos_i = int(pos)
            except ValueError:
                continue

            if not geno or geno == "--":
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
        for _ in range(50):
            line = fh.readline()
            if not line:
                break
            if not line.startswith("#"):
                break
            body = line.lstrip("#").strip()
            for key in ("file_id", "signature", "timestamp"):
                prefix = key + ":"
                if body.startswith(prefix):
                    meta[key] = body[len(prefix):].strip()
    return meta
