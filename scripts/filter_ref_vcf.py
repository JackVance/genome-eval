"""Filter 1000G reference VCFs to biallelic SNPs with unique non-empty IDs.

Required because conform-gt fails on multiallelic indels with duplicate
or empty ID fields. Writes a `<chr>.filtered.vcf.gz` alongside the
original; downstream scripts prefer the filtered file when present.

Parallelizes across chromosomes (multiprocessing) and writes output with
low gzip compression level (fast) since the filtered file is an
intermediate artifact, not long-term storage.

Usage:
    python scripts/filter_ref_vcf.py            # all chromosomes, parallel
    python scripts/filter_ref_vcf.py --chr 2    # just one
"""
from __future__ import annotations

import argparse
import gzip
import multiprocessing as mp
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
REF_DIR = PROJECT_ROOT / 'reference' / 'imputation' / '1kg_ref_b37'


def filter_vcf(src: Path, dst: Path):
    n_in = 0
    n_out = 0
    n_multi = 0
    n_indel = 0
    n_dup_id = 0
    seen_ids = set()
    # Low compression (level=1) for speed; the file is intermediate.
    # Batch writes for further speed.
    with gzip.open(src, 'rt', encoding='utf-8') as fin, gzip.open(dst, 'wt', compresslevel=1, encoding='utf-8') as fout:
        out_buffer = []
        for line in fin:
            if line.startswith('#'):
                out_buffer.append(line)
                if len(out_buffer) >= 1000:
                    fout.write(''.join(out_buffer))
                    out_buffer.clear()
                continue
            n_in += 1
            # Parse only the first 5 columns (faster than full split)
            try:
                p1, p2, p3, p4, p5, rest = line.split('\t', 5)
            except ValueError:
                continue
            rsid, ref, alt = p3, p4, p5
            if ',' in alt:
                n_multi += 1
                continue
            if len(ref) != 1 or len(alt) != 1 or ref not in 'ACGT' or alt not in 'ACGT':
                n_indel += 1
                continue
            if rsid == '.' or ';' in rsid:
                rsid = f'{p1}:{p2}'
            if rsid in seen_ids:
                n_dup_id += 1
                continue
            seen_ids.add(rsid)
            out_buffer.append(f'{p1}\t{p2}\t{rsid}\t{ref}\t{alt}\t{rest}')
            n_out += 1
            if len(out_buffer) >= 5000:
                fout.write(''.join(out_buffer))
                out_buffer.clear()
        if out_buffer:
            fout.write(''.join(out_buffer))
    return {'in': n_in, 'out': n_out, 'multi': n_multi, 'indel': n_indel, 'dup_id': n_dup_id}


def _worker(chrom: str) -> str:
    src_candidates = list(REF_DIR.glob(f'chr{chrom}.*.vcf.gz'))
    src_candidates = [p for p in src_candidates if '.filtered.' not in p.name]
    if not src_candidates:
        return f'chr{chrom}: no source VCF; skipped'
    src = src_candidates[0]
    dst = REF_DIR / f'chr{chrom}.1kg.phase3.v5a.filtered.vcf.gz'
    if dst.exists() and dst.stat().st_size > 0:
        return f'chr{chrom}: filtered already exists ({dst.stat().st_size:,} bytes); skipped'
    try:
        stats = filter_vcf(src, dst)
    except Exception as e:
        return f'chr{chrom}: FAILED {e!r}'
    return (f'chr{chrom}: in={stats["in"]:,} out={stats["out"]:,} '
            f'(multi={stats["multi"]:,}, indel={stats["indel"]:,}, dup_id={stats["dup_id"]:,})')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--chr', default='all')
    ap.add_argument('--workers', type=int, default=6,
                    help='Parallel worker processes (default 6; CPU-bound + I/O-bound mix)')
    args = ap.parse_args()

    if args.chr == 'all':
        chroms = [str(i) for i in range(1, 23)]
    else:
        chroms = args.chr.split(',')

    if len(chroms) == 1:
        print(_worker(chroms[0]), flush=True)
    else:
        with mp.Pool(args.workers) as pool:
            # Process in submission order but collect results as they finish
            for result in pool.imap_unordered(_worker, chroms):
                print(result, flush=True)


if __name__ == '__main__':
    main()
