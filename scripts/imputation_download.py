"""Download Beagle, reference panel, and genetic maps for imputation.

Files:
  - Beagle JAR (~5 MB)
  - Genetic maps for all autosomes (small zip, ~80 MB uncompressed)
  - 1000 Genomes Phase 3 bref3 reference panels (one per chromosome, ~200-500 MB each, ~10 GB total)

Use --chr to download only specific chromosomes (e.g., --chr 22 for validation).

Usage:
    python scripts/imputation_download.py --beagle          # Beagle JAR only
    python scripts/imputation_download.py --maps            # Genetic maps only
    python scripts/imputation_download.py --chr 22          # One chromosome of reference
    python scripts/imputation_download.py --all             # Everything
"""
from __future__ import annotations

import argparse
import sys
import zipfile
from pathlib import Path

import requests

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
IMPUTATION_DIR = PROJECT_ROOT / 'reference' / 'imputation'
BEAGLE_DIR = IMPUTATION_DIR / 'beagle'
REF_DIR = IMPUTATION_DIR / '1kg_ref_b37'
MAPS_DIR = IMPUTATION_DIR / 'genetic_maps'

# Known-working URLs as of 2026-04. Update if dead links surface.
BEAGLE_URL = 'https://faculty.washington.edu/browning/beagle/beagle.22Jul22.46e.jar'
BEAGLE_FILENAME = 'beagle.22Jul22.46e.jar'
CONFORM_GT_URL = 'https://faculty.washington.edu/browning/conform-gt/conform-gt.24May16.cee.jar'
CONFORM_GT_FILENAME = 'conform-gt.24May16.cee.jar'
MAPS_ZIP_URL = 'https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh37.map.zip'
REF_BASE_URL = 'https://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/b37.bref3'
REF_VCF_BASE_URL = 'https://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/b37.vcf'


def stream_download(url: str, target: Path, description: str = '', max_retries: int = 5) -> Path:
    """Download with retry-on-connection-reset and HTTP Range resume."""
    import time
    target.parent.mkdir(parents=True, exist_ok=True)
    # Check if already fully downloaded by issuing a HEAD
    head = None
    try:
        head = requests.head(url, timeout=30, allow_redirects=True)
    except Exception:
        pass
    total_expected = int(head.headers.get('Content-Length', 0)) if head is not None else 0

    if target.exists() and total_expected and target.stat().st_size == total_expected:
        print(f'  Already present: {target.name} ({target.stat().st_size:,} bytes)')
        return target

    attempt = 0
    while attempt < max_retries:
        attempt += 1
        existing = target.stat().st_size if target.exists() else 0
        headers = {}
        if existing > 0:
            headers['Range'] = f'bytes={existing}-'
            print(f'  Resuming {description or target.name} at {existing:,} bytes (attempt {attempt})')
        else:
            print(f'  Downloading {description or target.name} (attempt {attempt})')
        try:
            resp = requests.get(url, headers=headers, stream=True, timeout=(30, 120))
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
                        if total > 0 and written - last_print >= (20 << 20):
                            print(f'    {written / (1 << 20):,.0f} / {total / (1 << 20):,.0f} MB')
                            last_print = written
            # Verify size
            if total > 0 and target.stat().st_size < total:
                print(f'  File short: {target.stat().st_size:,} < {total:,}. Retrying.')
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


def download_beagle():
    print('[beagle]')
    stream_download(BEAGLE_URL, BEAGLE_DIR / BEAGLE_FILENAME, 'Beagle 5.4 JAR')
    print('[conform-gt]')
    return stream_download(CONFORM_GT_URL, BEAGLE_DIR / CONFORM_GT_FILENAME, 'conform-gt (allele aligner)')


def download_maps():
    print('[maps]')
    zip_path = MAPS_DIR / 'plink.GRCh37.map.zip'
    stream_download(MAPS_ZIP_URL, zip_path, 'Genetic map archive')
    # Extract
    print(f'  Extracting to {MAPS_DIR}...')
    with zipfile.ZipFile(zip_path, 'r') as zf:
        zf.extractall(MAPS_DIR)
    maps = sorted(MAPS_DIR.glob('plink.chr*.GRCh37.map'))
    print(f'  Extracted {len(maps)} per-chromosome map files.')


def download_ref(chrom: str, include_vcf: bool = True):
    """Download one chromosome's reference panel.

    Both bref3 (for Beagle's `ref=` argument) and VCF (for conform-gt's
    allele-alignment step) are needed for the full pipeline.
    """
    # bref3 format for Beagle
    filename = f'chr{chrom}.1kg.phase3.v5a.b37.bref3'
    url = f'{REF_BASE_URL}/{filename}'
    print(f'[ref chr{chrom} bref3]')
    stream_download(url, REF_DIR / filename, f'1000G Phase 3 EUR reference bref3, chr{chrom}')

    if include_vcf:
        vcf_filename = f'chr{chrom}.1kg.phase3.v5a.vcf.gz'
        vcf_url = f'{REF_VCF_BASE_URL}/{vcf_filename}'
        print(f'[ref chr{chrom} vcf]')
        stream_download(url=vcf_url, target=REF_DIR / vcf_filename,
                        description=f'1000G Phase 3 EUR reference VCF, chr{chrom}')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--beagle', action='store_true')
    ap.add_argument('--maps', action='store_true')
    ap.add_argument('--chr', help='Comma-separated chromosome numbers (e.g., 22 or 1,2,22) or "all"')
    ap.add_argument('--all', action='store_true', help='Beagle + maps + all autosomes')
    args = ap.parse_args()

    if args.all:
        download_beagle()
        download_maps()
        for c in [str(i) for i in range(1, 23)]:
            download_ref(c)
        return

    if args.beagle:
        download_beagle()
    if args.maps:
        download_maps()
    if args.chr:
        chroms = [str(i) for i in range(1, 23)] if args.chr == 'all' else args.chr.split(',')
        for c in chroms:
            download_ref(c)

    if not any([args.beagle, args.maps, args.chr, args.all]):
        ap.print_help()
        sys.exit(1)


if __name__ == '__main__':
    main()
