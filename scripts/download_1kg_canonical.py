"""Download canonical 1000G Phase 3 release VCFs from EBI.

Unlike the bochet imputation-panel VCFs we already have (INFO field is `.`),
these canonical release VCFs include per-population AF in the INFO field
(EUR_AF, AFR_AF, AMR_AF, EAS_AF, SAS_AF). We parse EUR_AF out of them once
and cache locally, eliminating the Ensembl API bottleneck for every future
PRS run.

Downloads to `reference/imputation/1kg_ebi_release/`. ~15 GB across 22 autosomes.
Supports retry + resume. Safe to delete the raw VCFs after running
`scripts/extract_eur_afs.py` (which builds the compact AF lookup parquet).

Usage:
    python scripts/download_1kg_canonical.py              # all autosomes
    python scripts/download_1kg_canonical.py --chr 22     # single
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import requests

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
DEST_DIR = PROJECT_ROOT / 'reference' / 'imputation' / '1kg_ebi_release'

BASE_URL = 'http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502'
FILENAME_TPL = 'ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz'


def stream_download(url: str, target: Path, max_retries: int = 6) -> Path:
    target.parent.mkdir(parents=True, exist_ok=True)
    try:
        head = requests.head(url, timeout=30, allow_redirects=True)
        total_expected = int(head.headers.get('Content-Length', 0))
    except Exception:
        total_expected = 0

    if target.exists() and total_expected and target.stat().st_size == total_expected:
        print(f'  Already present: {target.name} ({target.stat().st_size:,} bytes)')
        return target

    attempt = 0
    while attempt < max_retries:
        attempt += 1
        existing = target.stat().st_size if target.exists() else 0
        headers = {'Range': f'bytes={existing}-'} if existing > 0 else {}
        action = 'Resuming' if existing > 0 else 'Downloading'
        print(f'  {action} {target.name} at {existing:,} bytes (attempt {attempt})')
        try:
            resp = requests.get(url, headers=headers, stream=True, timeout=(30, 180))
            if resp.status_code not in (200, 206):
                resp.raise_for_status()
            mode = 'ab' if existing > 0 and resp.status_code == 206 else 'wb'
            if mode == 'wb':
                existing = 0
            written = existing
            total = total_expected or (existing + int(resp.headers.get('Content-Length', 0)))
            last_print = existing
            with open(target, mode) as fh:
                for chunk in resp.iter_content(chunk_size=1 << 20):
                    if chunk:
                        fh.write(chunk)
                        written += len(chunk)
                        if total > 0 and written - last_print >= (50 << 20):
                            print(f'    {written / (1 << 20):,.0f} / {total / (1 << 20):,.0f} MB', flush=True)
                            last_print = written
            if total > 0 and target.stat().st_size < total:
                print(f'  File short ({target.stat().st_size:,}<{total:,}); retrying.')
                continue
            print(f'  Done: {target.name} ({target.stat().st_size:,} bytes)')
            return target
        except (requests.exceptions.ChunkedEncodingError,
                requests.exceptions.ConnectionError,
                requests.exceptions.ReadTimeout,
                ConnectionResetError) as e:
            print(f'  Connection error on attempt {attempt}: {e}')
            time.sleep(5 * attempt)
            continue
    raise RuntimeError(f'Failed to download {url} after {max_retries} attempts')


def download_chrom(chrom: str):
    filename = FILENAME_TPL.format(chrom=chrom)
    url = f'{BASE_URL}/{filename}'
    target = DEST_DIR / filename
    print(f'[chr{chrom}] {url}')
    stream_download(url, target)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--chr', default='all', help='Comma-separated chromosomes or "all"')
    args = ap.parse_args()
    chroms = [str(i) for i in range(1, 23)] if args.chr == 'all' else args.chr.split(',')
    for c in chroms:
        download_chrom(c)


if __name__ == '__main__':
    main()
