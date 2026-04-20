"""Fetch 1000 Genomes European (EUR) allele frequencies for PRS variants
that overlap with a subject's chip. Uses Ensembl REST API single-variant
endpoint (batch endpoint doesn't return population data even with pops=1).

Caches results to reference/population_cache/1kg_eur_afs.json keyed by rsID.

Usage:
    python scripts/prs_fetch_afs.py PGS000297 alice
"""
from __future__ import annotations

import argparse
import json
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd
import requests

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
sys.path.insert(0, str(SCRIPT_DIR))

from prs_pipeline import load_weights

CACHE_DIR = PROJECT_ROOT / 'reference' / 'population_cache'
CACHE_FILE = CACHE_DIR / '1kg_eur_afs.json'
WEIGHTS_DIR = PROJECT_ROOT / 'reference' / 'prs_weights'


def load_cache() -> dict:
    if CACHE_FILE.exists():
        return json.loads(CACHE_FILE.read_text(encoding='utf-8'))
    return {}


def save_cache(cache: dict) -> None:
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    CACHE_FILE.write_text(json.dumps(cache, indent=2), encoding='utf-8')


def fetch_eur_af_single(rsid: str, timeout: int = 15) -> dict:
    """Fetch EUR allele frequencies for a single rsID from Ensembl."""
    url = f'https://rest.ensembl.org/variation/human/{rsid}'
    headers = {'Accept': 'application/json'}
    resp = requests.get(url, headers=headers, params={'pops': 1}, timeout=timeout)
    if resp.status_code == 404:
        return {}
    resp.raise_for_status()
    data = resp.json()
    out = {}
    for pop in data.get('populations', []):
        if pop.get('population') == '1000GENOMES:phase_3:EUR':
            allele = pop.get('allele')
            freq = pop.get('frequency')
            if allele is not None and freq is not None:
                out[allele] = float(freq)
    return out


def _fetch_one_with_retry(rsid: str) -> tuple[str, dict]:
    """Worker function for the thread pool — returns (rsid, af_dict)."""
    try:
        return rsid, fetch_eur_af_single(rsid)
    except Exception:
        time.sleep(1)
        try:
            return rsid, fetch_eur_af_single(rsid)
        except Exception:
            return rsid, {}


def main(pgs_id: str, subject_id: str, workers: int = 10, save_every: int = 100):
    weights_path = WEIGHTS_DIR / f'{pgs_id}_hmPOS_GRCh37.txt.gz'
    weights, _meta = load_weights(weights_path)

    # Prefer imputed parquet if present (it has ~50x more overlap with PRS panels)
    imputed = PROJECT_ROOT / 'standardized-genomes' / 'imputed' / f'{subject_id}.imputed.parquet'
    chip = PROJECT_ROOT / 'standardized-genomes' / f'{subject_id}.parquet'
    subject_parquet = imputed if imputed.exists() else chip
    print(f'Using {"IMPUTED" if subject_parquet is imputed else "CHIP"} parquet: {subject_parquet}')
    subject = pd.read_parquet(subject_parquet)
    subject['chrom'] = subject['chrom'].astype(str)

    overlap = weights.merge(subject, left_on=['chr', 'pos'], right_on=['chrom', 'pos'], how='inner')
    rsids = [r for r in overlap['rsID'].dropna().tolist() if isinstance(r, str) and r.startswith('rs')]
    rsids = list(dict.fromkeys(rsids))
    print(f'{pgs_id} x {subject_id}: {len(rsids)} overlap rsIDs to look up')

    cache = load_cache()
    to_fetch = [r for r in rsids if r not in cache]
    print(f'Cached: {len(rsids) - len(to_fetch)}; fetching: {len(to_fetch)} with {workers} parallel workers')

    cache_lock = threading.Lock()
    completed = 0
    start = time.time()

    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(_fetch_one_with_retry, rsid): rsid for rsid in to_fetch}
        for fut in as_completed(futures):
            rsid, af = fut.result()
            with cache_lock:
                cache[rsid] = af
                completed += 1
                if completed % save_every == 0:
                    save_cache(cache)
                    elapsed = time.time() - start
                    rate = completed / elapsed
                    eta = (len(to_fetch) - completed) / rate if rate > 0 else 0
                    print(f'  [{completed}/{len(to_fetch)}] {rate:.1f} req/s, ETA {eta/60:.1f} min; '
                          f'{sum(1 for v in cache.values() if v)} non-empty in cache',
                          flush=True)

    save_cache(cache)
    populated = sum(1 for r in rsids if cache.get(r))
    elapsed = time.time() - start
    print()
    print(f'Done in {elapsed:.0f}s. {populated} / {len(rsids)} overlap rsIDs have EUR AF data.')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('pgs_id')
    ap.add_argument('subject_id')
    ap.add_argument('--workers', type=int, default=10)
    args = ap.parse_args()
    main(args.pgs_id, args.subject_id, args.workers)
