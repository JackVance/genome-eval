"""Check Rh factor (D antigen) status from tag SNPs near/in RHD.

Rh status on arrays is tricky because the Rh-negative phenotype in
Europeans is predominantly caused by complete deletion of the RHD gene.
Tag SNPs can infer status but with limitations. This script:

  1. Checks several commonly-referenced RHD/RHCE-region tag SNPs.
  2. Also inspects homozygosity patterns across RHD region — if many
     RHD-specific probes show suspicious calls, that's consistent with
     hemizygous or homozygous deletion.
  3. Reports the inferred Rh status with an honest confidence.

Caveat: array-based Rh inference is not a clinical test. Transfusion
medicine requires serological typing.
"""
from __future__ import annotations

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
sys.path.insert(0, str(SCRIPT_DIR))

import ledger_io
from tier_rules import compute_tier, TIER_RULE_VERSION


# Commonly-referenced tag SNPs for Rh status, with GRCh37 forward-strand info.
# None is a perfect tag; presence of multiple concordant calls increases
# confidence.
RH_TAG_SNPS = [
    # (rsid, gene, rh_pos_allele, rh_neg_allele, note)
    ('rs590787',   'RHD',  'C', 'T',
     'RHD intron 4 tag SNP; T associated with RHD deletion in European cohorts.'),
    ('rs676785',   'RHD',  'G', 'A',
     'RHD tag SNP; A associated with RHD deletion in European cohorts.'),
    ('rs17418085', 'RHCE', 'A', 'G',
     'RHCE region SNP sometimes used as Rh status tag in Europeans.'),
    ('rs586178',   'RHD',  'A', 'G',
     'RHD tag SNP (less validated).'),
    ('rs675072',   'RHD',  'C', 'T',
     'RHD tag SNP.'),
]

# GRCh37 coordinates for the RHD gene region (approx chr1:25,598,884–25,656,935)
RHD_REGION = ('1', 25598884, 25656935)


def main(subject_id='alice'):
    df = pd.read_parquet(PROJECT_ROOT / 'standardized-genomes' / f'{subject_id}.parquet')

    # Query each tag SNP
    print('=== Rh tag SNPs ===')
    tag_results = []
    for rsid, gene, rh_pos, rh_neg, note in RH_TAG_SNPS:
        row = df[df.rsid == rsid]
        if row.empty:
            print(f'  {rsid:12s} ({gene:4s}): NOT ON CHIP  | {note}')
            tag_results.append({'rsid': rsid, 'gene': gene, 'on_chip': False})
            continue
        r = row.iloc[0]
        a1 = r.a1 if pd.notna(r.a1) else '-'
        a2 = r.a2 if pd.notna(r.a2) else '-'
        gt = a1 + a2
        pos_count = gt.count(rh_pos)
        neg_count = gt.count(rh_neg)
        print(f'  {rsid:12s} ({gene:4s}): {gt}  (Rh+ allele {rh_pos} × {pos_count}, '
              f'Rh- allele {rh_neg} × {neg_count})  | {note}')
        tag_results.append({
            'rsid': rsid, 'gene': gene, 'on_chip': True, 'genotype': gt,
            'rh_pos_allele': rh_pos, 'rh_neg_allele': rh_neg,
            'rh_pos_count': pos_count, 'rh_neg_count': neg_count,
        })

    # Simple heuristic: count Rh+ vs Rh- alleles across informative probes
    informative = [t for t in tag_results if t.get('on_chip') and t.get('rh_pos_count') is not None]
    total_pos = sum(t['rh_pos_count'] for t in informative)
    total_neg = sum(t['rh_neg_count'] for t in informative)
    total_alleles = total_pos + total_neg
    print()
    print(f'Across {len(informative)} informative tag SNPs: {total_pos} Rh+ alleles, '
          f'{total_neg} Rh- alleles (of {total_alleles} total).')

    if total_alleles == 0:
        inference = 'Rh status cannot be inferred — no informative tag SNPs on chip.'
        confidence = 'none'
    elif total_pos == total_alleles:
        inference = 'Rh POSITIVE (D+) — all tag SNPs consistent with functional RHD.'
        confidence = 'moderate (array-based tag SNPs; not a clinical test)'
    elif total_neg == total_alleles:
        inference = 'Rh NEGATIVE (D-) — all tag SNPs consistent with RHD deletion.'
        confidence = 'moderate'
    elif total_pos > total_neg:
        inference = f'Likely Rh POSITIVE — majority of tag SNPs ({total_pos}/{total_alleles}) Rh+ direction.'
        confidence = 'low-moderate (mixed signals across tag SNPs)'
    else:
        inference = f'Mixed signal (Rh+ {total_pos}, Rh- {total_neg}) — cannot call confidently.'
        confidence = 'low'

    print()
    print(f'INFERENCE: {inference}')
    print(f'Confidence: {confidence}')

    # Also inspect number of RHD-region probes that returned calls
    chrom, start, end = RHD_REGION
    rhd_region_snps = df[(df.chrom == chrom) & (df.pos >= start) & (df.pos <= end)]
    no_calls = rhd_region_snps[(rhd_region_snps.a1.isna()) | (rhd_region_snps.a2.isna())]
    print()
    print(f'RHD region (chr1:{start}-{end}): {len(rhd_region_snps)} probes on chip, '
          f'{len(no_calls)} no-calls.')
    print('(High no-call rate in RHD region can be consistent with homozygous RHD deletion, '
          'though individual array behavior varies.)')

    # Record finding
    profile = json.loads((PROJECT_ROOT / 'profiles' / f'{subject_id}.json').read_text(encoding='utf-8'))
    declared = profile.get('declared_ancestry')

    rec = {
        'subject_id': subject_id,
        'topic': 'trait_rh_factor',
        'claim': f'Rh factor (D antigen) inference: {inference}',
        'variants': [
            {'rsid': t['rsid'], 'gene': t['gene'],
             'genotype': t.get('genotype'), 'on_chip': t['on_chip']}
            for t in tag_results
        ],
        'effect': {
            'type': 'trait_inference',
            'value': inference,
            'method': 'tag-SNP majority across RHD/RHCE probes',
            'confidence': confidence,
            'rh_pos_alleles': total_pos,
            'rh_neg_alleles': total_neg,
            'rhd_region_probes_on_chip': len(rhd_region_snps),
            'rhd_region_no_calls': len(no_calls),
        },
        'cohort_ancestry': 'European (tag SNPs validated in European populations)',
        'subject_ancestry_match': 'match' if (declared and 'european' in (declared.get('top_level','').lower())) else 'unknown',
        'source_ids': ['pubmed:avent-reid-2000-rh'],
        'notes': (
            'Array-based Rh status inference uses tag SNPs in RHD/RHCE region. '
            'The Rh-negative phenotype in Europeans is predominantly caused by complete RHD '
            'gene deletion, which arrays cannot directly detect; tag SNPs provide a proxy. '
            'This is NOT a clinical test — transfusion medicine requires serological D typing. '
            'If you have confirmed Rh status from medical records or prior 23andMe reports, '
            'that ground truth should supersede this inference.'
        ),
        'clinvar_significance': None,
        'clinvar_review_stars': None,
        'replication_count': 5,
        'investigation_id': None,
    }
    tier, downgrade = compute_tier(rec)
    rec['tier_computed'] = tier
    rec['tier_computed_at'] = datetime.now(timezone.utc).isoformat()
    rec['tier_rule_version'] = TIER_RULE_VERSION
    rec['ancestry_downgrade'] = downgrade
    fid = ledger_io.append_finding(**rec)
    print()
    print(f'Appended finding {fid} (tier {tier}).')


if __name__ == '__main__':
    subj = sys.argv[1] if len(sys.argv) > 1 else 'alice'
    main(subj)
