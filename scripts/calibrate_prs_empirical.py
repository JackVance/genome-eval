"""Compute empirical PRS reference distributions on 1000G Phase 3 EUR samples.

For each PGS, applies the score to each of the 503 EUR individuals in 1000G
(CEU + FIN + GBR + IBS + TSI) and saves the empirical mean and SD of the
resulting distribution. These are then used by run_prs.py to produce
properly-calibrated z-scores, replacing the theoretical independence-formula
approximation that fails for LDpred2-style scores.

Output: reference/population_cache/prs_empirical/<PGS_id>.json, with:
  - mean, SD, min, max, n_samples
  - per-sample scores (for future re-analysis)
  - reference description

Parallelizes across chromosomes. Expects the canonical EBI 1000G release VCFs
at reference/imputation/1kg_ebi_release/ (already downloaded) and the sample
panel metadata file at that directory's metadata/ subfolder.

Usage:
    python scripts/calibrate_prs_empirical.py PGS002804 PGS002231 PGS002135
    python scripts/calibrate_prs_empirical.py all          # every PGS in prs_traits.py
"""
from __future__ import annotations

import argparse
import gzip
import json
import multiprocessing as mp
import statistics
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
sys.path.insert(0, str(SCRIPT_DIR))

from prs_pipeline import load_weights, COMPLEMENT

VCF_DIR = PROJECT_ROOT / 'reference' / 'imputation' / '1kg_ebi_release'
PANEL_FILE = VCF_DIR / 'metadata' / 'integrated_call_samples_v3.20130502.ALL.panel'
OUT_DIR = PROJECT_ROOT / 'reference' / 'population_cache' / 'prs_empirical'
WEIGHTS_DIR = PROJECT_ROOT / 'reference' / 'prs_weights'

REFERENCE_DESCRIPTION = (
    '1000 Genomes Phase 3 EUR samples (N=503; sub-populations: CEU=99, FIN=99, '
    'GBR=91, IBS=107, TSI=107). Standard publicly-available European reference; '
    'weighted toward Northern/Western European sub-populations, limited '
    'representation of Southern/Eastern/Ashkenazi. Empirical mean/SD from this '
    'sample is an approximation of the true EUR-population distribution, subject '
    'to sampling error (~3% on SD) and sub-population weighting bias.'
)


def load_eur_sample_ids() -> list[str]:
    eur = []
    with open(PANEL_FILE, 'r', encoding='utf-8') as f:
        for i, line in enumerate(f):
            if i == 0:
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 3 and parts[2] == 'EUR':
                eur.append(parts[0])
    return eur


def compute_chromosome(args):
    """Worker: parse one chromosome VCF, return per-sample partial scores per PGS."""
    chrom, pgs_weights_list, eur_sample_ids = args

    # Build per-PGS position lookup for this chromosome
    combined_positions = set()
    per_pgs: dict[str, dict] = {}
    for pgs_id, weights in pgs_weights_list:
        pw = {}
        for r in weights[weights['chr'].astype(str) == chrom].itertuples():
            key = (str(r.chr), int(r.pos))
            pw[key] = (r.effect_allele, r.effect_weight, r.other_allele)
            combined_positions.add(key)
        per_pgs[pgs_id] = pw

    scores = {pgs_id: {s: 0.0 for s in eur_sample_ids} for pgs_id, _ in pgs_weights_list}
    eur_set = set(eur_sample_ids)

    vcf_path = VCF_DIR / f'ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz'
    if not vcf_path.exists():
        return f'chr{chrom}: VCF missing', {pgs_id: {} for pgs_id, _ in pgs_weights_list}

    sample_col_idx = None
    n_matched = 0
    n_skipped_allele_mismatch = 0

    with gzip.open(vcf_path, 'rt', encoding='utf-8') as fh:
        for line in fh:
            if line.startswith('#CHROM'):
                headers = line.rstrip('\n').split('\t')
                all_samples = headers[9:]
                sample_col_idx = [(i, sid) for i, sid in enumerate(all_samples) if sid in eur_set]
                continue
            if line.startswith('#'):
                continue

            # Fast position check before full parse
            try:
                i1 = line.index('\t')
                i2 = line.index('\t', i1 + 1)
            except ValueError:
                continue
            chrom_v = line[:i1]
            pos_str = line[i1 + 1:i2]
            try:
                pos = int(pos_str)
            except ValueError:
                continue
            key = (chrom_v, pos)
            if key not in combined_positions:
                continue

            parts = line.rstrip('\n').split('\t')
            if len(parts) < 10:
                continue
            _chrom_v, _pos_s, _rsid, ref, alt, _qual, _filter, _info, fmt, *samples = parts
            if ',' in alt:
                # multiallelic — skip (PGS positions are typically biallelic)
                continue
            if len(ref) != 1 or len(alt) != 1:
                continue  # indel

            # Find GT field index in FORMAT
            fmt_fields = fmt.split(':')
            try:
                gt_idx = fmt_fields.index('GT')
            except ValueError:
                continue

            for pgs_id, _ in pgs_weights_list:
                pw = per_pgs[pgs_id]
                if key not in pw:
                    continue
                ea, beta, oa = pw[key]

                if ea == alt:
                    ea_char = '1'
                elif ea == ref:
                    ea_char = '0'
                elif COMPLEMENT.get(ea) == alt:
                    ea_char = '1'
                elif COMPLEMENT.get(ea) == ref:
                    ea_char = '0'
                else:
                    n_skipped_allele_mismatch += 1
                    continue

                # Per-sample genotype parse
                for col_idx, sid in sample_col_idx:
                    try:
                        sample_col = samples[col_idx]
                    except IndexError:
                        continue
                    if ':' in sample_col:
                        gt_str = sample_col.split(':')[gt_idx]
                    else:
                        gt_str = sample_col
                    # Phased ("0|1") or unphased ("0/1")
                    if '|' in gt_str:
                        a, b = gt_str.split('|', 1)
                    elif '/' in gt_str:
                        a, b = gt_str.split('/', 1)
                    else:
                        continue
                    if a == '.' or b == '.':
                        continue
                    count = (1 if a == ea_char else 0) + (1 if b == ea_char else 0)
                    scores[pgs_id][sid] += count * beta

                n_matched += 1

    status = f'chr{chrom}: matched {n_matched} site-PGS pairs (skipped {n_skipped_allele_mismatch} allele-mismatch)'
    return status, scores


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('pgs_ids', nargs='+', help='PGS IDs to calibrate, or "all"')
    ap.add_argument('--workers', type=int, default=6)
    ap.add_argument(
        '--subject-observed',
        help='Optional subject ID. Restricts calibration to variants present '
             'in <subject>.imputed.parquet (or chip parquet). This produces '
             'an apples-to-apples empirical distribution for any score run '
             'against a subject whose coverage is sub-100%. Output file is '
             'suffixed with ".<subject>" to avoid overwriting the full '
             'calibration.',
    )
    args = ap.parse_args()

    # Resolve "all"
    pgs_ids = args.pgs_ids
    if pgs_ids == ['all']:
        from prs_traits import TRAITS
        pgs_ids = sorted({cfg['pgs_id'] for cfg in TRAITS.values()})

    # Load EUR sample list
    if not PANEL_FILE.exists():
        sys.exit(f'Panel file missing: {PANEL_FILE}')
    eur = load_eur_sample_ids()
    print(f'EUR reference samples: {len(eur)}')

    # Optional subject-observed-set filter for per-subject comparable z-scores.
    observed_set = None
    subject_suffix = ''
    if args.subject_observed:
        import pandas as pd
        subj_imputed = PROJECT_ROOT / 'standardized-genomes' / 'imputed' / f'{args.subject_observed}.imputed.parquet'
        subj_chip = PROJECT_ROOT / 'standardized-genomes' / f'{args.subject_observed}.parquet'
        src = subj_imputed if subj_imputed.exists() else subj_chip
        print(f'Subject-observed filter: reading {src.name}')
        sdf = pd.read_parquet(src, columns=['chrom', 'pos'])
        sdf['chrom'] = sdf['chrom'].astype(str)
        observed_set = set(zip(sdf['chrom'], sdf['pos'].astype(int)))
        subject_suffix = f'.{args.subject_observed}'
        print(f'Subject positions available: {len(observed_set):,}')

    # Load weights for each PGS
    weights_list = []
    for pgs in pgs_ids:
        wpath = WEIGHTS_DIR / f'{pgs}_hmPOS_GRCh37.txt.gz'
        if not wpath.exists():
            print(f'  SKIP {pgs}: weights file missing at {wpath}')
            continue
        w, _ = load_weights(wpath)
        w = w.copy()
        w['chr'] = w['chr'].astype(str)
        if observed_set is not None:
            before = len(w)
            keys = list(zip(w['chr'], w['pos'].astype(int)))
            w = w[[k in observed_set for k in keys]].reset_index(drop=True)
            print(f'  {pgs}: {len(w):,} variants ({before - len(w):,} filtered out as not in subject observation)')
        else:
            print(f'  {pgs}: {len(w):,} variants loaded')
        weights_list.append((pgs, w))

    if not weights_list:
        sys.exit('No PGS weights to process.')

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # Parallelize across chromosomes; every worker sees all PGS at once
    chroms = [str(i) for i in range(1, 23)]
    args_list = [(c, weights_list, eur) for c in chroms]

    final = {pgs_id: {s: 0.0 for s in eur} for pgs_id, _ in weights_list}

    with mp.Pool(args.workers) as pool:
        for status, chrom_scores in pool.imap_unordered(compute_chromosome, args_list):
            print(' ', status, flush=True)
            for pgs_id, sample_scores in chrom_scores.items():
                if pgs_id not in final:
                    continue
                for sid, s in sample_scores.items():
                    final[pgs_id][sid] += s

    # Save empirical distributions
    print()
    for pgs_id, sample_scores in final.items():
        values = list(sample_scores.values())
        empirical = {
            'pgs_id': pgs_id,
            'reference': REFERENCE_DESCRIPTION,
            'n_samples': len(values),
            'mean': statistics.mean(values),
            'sd': statistics.stdev(values) if len(values) > 1 else 0.0,
            'min': min(values),
            'max': max(values),
            'per_sample': sample_scores,
        }
        out_path = OUT_DIR / f'{pgs_id}{subject_suffix}.json'
        out_path.write_text(json.dumps(empirical, indent=2), encoding='utf-8')
        print(f'  {pgs_id}: empirical mean={empirical["mean"]:.4f}  SD={empirical["sd"]:.4f}  '
              f'(range {empirical["min"]:.3f} to {empirical["max"]:.3f})  -> {out_path.name}')


if __name__ == '__main__':
    main()
