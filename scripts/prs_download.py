"""Download PGS Catalog scoring files (harmonized to a given build).

Usage:
    python scripts/prs_download.py PGS000297
    python scripts/prs_download.py PGS000297 --build GRCh37

Pulls harmonized file from the EBI FTP mirror. Saves to
`reference/prs_weights/<pgs_id>_hmPOS_<build>.txt.gz`.

Files are generally small (a few MB) for clumped scores; can be 100+ MB
for LDpred2-auto / PRS-CS scores with millions of variants.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import requests

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
WEIGHTS_DIR = PROJECT_ROOT / 'reference' / 'prs_weights'


def download_pgs(pgs_id: str, build: str = 'GRCh37') -> Path:
    """Download a PGS Catalog harmonized scoring file."""
    WEIGHTS_DIR.mkdir(parents=True, exist_ok=True)
    filename = f'{pgs_id}_hmPOS_{build}.txt.gz'
    target = WEIGHTS_DIR / filename
    if target.exists() and target.stat().st_size > 0:
        print(f'Already present: {target} ({target.stat().st_size:,} bytes)')
        return target

    url = (
        f'https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/'
        f'{pgs_id}/ScoringFiles/Harmonized/{filename}'
    )
    print(f'Downloading: {url}')
    resp = requests.get(url, stream=True, timeout=120)
    resp.raise_for_status()
    total = 0
    with open(target, 'wb') as fh:
        for chunk in resp.iter_content(chunk_size=65536):
            if chunk:
                fh.write(chunk)
                total += len(chunk)
    print(f'Saved: {target} ({total:,} bytes)')
    return target


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('pgs_id', help='PGS Catalog ID, e.g., PGS000297')
    ap.add_argument('--build', default='GRCh37', choices=['GRCh37', 'GRCh38'])
    args = ap.parse_args()
    try:
        download_pgs(args.pgs_id, args.build)
    except Exception as e:
        print(f'FAILED: {e}', file=sys.stderr)
        sys.exit(1)
