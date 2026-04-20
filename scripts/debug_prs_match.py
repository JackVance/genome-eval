"""Diagnose why PRS coverage is low and what's causing allele mismatches."""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
sys.path.insert(0, str(SCRIPT_DIR))

from prs_pipeline import load_weights

weights, meta = load_weights(PROJECT_ROOT / 'reference' / 'prs_weights' / 'PGS000297_hmPOS_GRCh37.txt.gz')
subject = pd.read_parquet(PROJECT_ROOT / 'standardized-genomes' / 'alice.parquet')
subject['chrom'] = subject['chrom'].astype(str)

print('Weight file head:')
print(weights.head())
print()
print(f'Weight chr values: {sorted(weights["chr"].unique())[:15]}')
print(f'Subject chrom values: {sorted(subject["chrom"].unique())[:15]}')
print()

# Try match by rsID first
rsid_match = weights.merge(subject, left_on='rsID', right_on='rsid', how='left')
print(f'rsID-based match: {rsid_match["rsid"].notna().sum()} / {len(weights)}')

# Now (chr, pos) match
subject['chrom_str'] = subject['chrom'].astype(str)
pos_match = weights.merge(
    subject,
    left_on=['chr', 'pos'],
    right_on=['chrom_str', 'pos'],
    how='left',
)
print(f'(chr, pos)-based match: {pos_match["a1"].notna().sum()} / {len(weights)}')

# Look at examples of allele mismatches for the chr,pos-matched ones
matched = pos_match[pos_match['a1'].notna()].copy()
print()
print(f'{len(matched)} SNPs matched by position. Examining allele match:')
COMPLEMENT = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
direct = 0
flipped = 0
palindromic = 0
mismatch = 0
examples = []
for r in matched.itertuples():
    ea = r.effect_allele
    oa = r.other_allele
    is_palindromic = {ea, oa} in [{'A', 'T'}, {'C', 'G'}]
    if ea in (r.a1, r.a2):
        direct += 1
    elif is_palindromic:
        palindromic += 1
    elif COMPLEMENT.get(ea) in (r.a1, r.a2):
        flipped += 1
    else:
        mismatch += 1
        if len(examples) < 10:
            examples.append(f'  chr{r.chr}:{r.pos} ea={ea} oa={oa} subject={r.a1}{r.a2}')

print(f'  direct match: {direct}')
print(f'  strand flip:  {flipped}')
print(f'  palindromic:  {palindromic}')
print(f'  mismatch:     {mismatch}')
print()
print('Example mismatches:')
for e in examples:
    print(e)
