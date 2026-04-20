"""Extract EUR_AF per variant from canonical 1000G Phase 3 release VCFs.

Streams each chromosome's VCF (~500 MB to 2 GB) and pulls out:
    chrom, pos, ref, alt, rsid, eur_af

from the INFO field. Writes a compact `reference/population_cache/1kg_eur_afs.parquet`
with all autosomes merged. That file becomes the local AF lookup source
for all future PRS runs — no Ensembl API needed.

Parallelizes across chromosomes with multiprocessing.

After running, the raw `reference/imputation/1kg_ebi_release/` VCFs can be
deleted to reclaim disk.

Usage:
    python scripts/extract_eur_afs.py                 # all autosomes
    python scripts/extract_eur_afs.py --chr 22        # single
"""
from __future__ import annotations

import argparse
import gzip
import multiprocessing as mp
import sys
from pathlib import Path

import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
SRC_DIR = PROJECT_ROOT / 'reference' / 'imputation' / '1kg_ebi_release'
OUT_DIR = PROJECT_ROOT / 'reference' / 'population_cache'

FILENAME_TPL = 'ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz'


def parse_info_for_eur_af(info: str) -> float | None:
    """Pull EUR_AF out of the INFO field. Returns None if missing or multi-allelic."""
    for kv in info.split(';'):
        if kv.startswith('EUR_AF='):
            val = kv.split('=', 1)[1]
            if ',' in val:
                # Multi-allelic — would need per-alt parsing; skip.
                return None
            try:
                return float(val)
            except ValueError:
                return None
    return None


def process_chrom(chrom: str) -> str:
    src = SRC_DIR / FILENAME_TPL.format(chrom=chrom)
    if not src.exists():
        return f'chr{chrom}: missing {src.name}; skipped'
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    dst = OUT_DIR / f'1kg_eur_afs.chr{chrom}.parquet'
    if dst.exists():
        return f'chr{chrom}: already extracted ({dst.stat().st_size:,} bytes); skipped'

    rows = []
    n_in = 0
    n_multi = 0
    n_no_eur = 0
    with gzip.open(src, 'rt', encoding='utf-8') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            n_in += 1
            try:
                p_chrom, p_pos, p_rsid, p_ref, p_alt, _q, _f, info, *_ = line.split('\t', 8)
            except ValueError:
                continue
            if ',' in p_alt:
                n_multi += 1
                continue
            if len(p_ref) != 1 or len(p_alt) != 1:
                continue  # indel
            eur_af = parse_info_for_eur_af(info)
            if eur_af is None:
                n_no_eur += 1
                continue
            rows.append((p_chrom, int(p_pos), p_ref, p_alt,
                         p_rsid if p_rsid != '.' else None, eur_af))

    df = pd.DataFrame(rows, columns=['chrom', 'pos', 'ref', 'alt', 'rsid', 'eur_af'])
    df.to_parquet(dst, index=False)
    return (f'chr{chrom}: in={n_in:,} kept={len(df):,} '
            f'(multi={n_multi:,}, no_eur={n_no_eur:,}) -> {dst.name}')


def merge_all():
    """Combine per-chromosome parquets into one sorted lookup file."""
    parts = sorted(OUT_DIR.glob('1kg_eur_afs.chr*.parquet'))
    if not parts:
        print('No per-chromosome parquets to merge.')
        return
    frames = [pd.read_parquet(p) for p in parts]
    merged = pd.concat(frames, ignore_index=True)
    merged['chrom'] = merged['chrom'].astype(str)
    merged = merged.sort_values(['chrom', 'pos']).reset_index(drop=True)
    out = OUT_DIR / '1kg_eur_afs.parquet'
    merged.to_parquet(out, index=False)
    print(f'Merged: {len(merged):,} variants -> {out} ({out.stat().st_size:,} bytes)')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--chr', default='all')
    ap.add_argument('--workers', type=int, default=6)
    ap.add_argument('--skip-merge', action='store_true')
    args = ap.parse_args()

    chroms = [str(i) for i in range(1, 23)] if args.chr == 'all' else args.chr.split(',')

    if len(chroms) == 1:
        print(process_chrom(chroms[0]))
    else:
        with mp.Pool(args.workers) as pool:
            for result in pool.imap_unordered(process_chrom, chroms):
                print(result, flush=True)

    if not args.skip_merge:
        merge_all()


if __name__ == '__main__':
    main()
