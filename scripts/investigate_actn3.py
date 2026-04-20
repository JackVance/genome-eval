"""ACTN3 rs1815739 (R577X) — sprint vs. endurance predisposition.

Single well-studied SNP in alpha-actinin-3. XX homozygotes make no functional
ACTN3 protein (fast-twitch muscle fibers); RR homozygotes make the full
protein. XX is overrepresented in elite endurance athletes; RR in sprint/power.
Effect size is small (~2-3% of elite-athlete performance variance), undetectable
at recreational level. Population frequency: ~18% XX in Europeans.

Reports genotype and writes a trait finding to the ledger.
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

RSID = 'rs1815739'
CHROM = '11'
POS = 66560624  # GRCh37
# REF=C (R577, Arg577) produces functional alpha-actinin-3 — sprint/power direction
# ALT=T (X577, stop codon)  — no functional protein — endurance direction
REF_ALLELE = 'C'
ALT_ALLELE = 'T'


def interpret(gt: str) -> tuple[str, str]:
    """Return (plain-language result, technical explanation)."""
    c = sorted(gt)
    if c == ['C', 'C']:
        return (
            'RR homozygote (sprint/power direction)',
            'Both copies produce functional alpha-actinin-3. Overrepresented in elite sprinters '
            'and power athletes (~30% of Europeans); at the recreational level the effect is '
            'too small to matter.'
        )
    if c == ['C', 'T']:
        return (
            'RX heterozygote (carrier — no clear direction)',
            'One functional copy, one stop-codon copy. Intermediate. Most common genotype in '
            'Europeans (~50%).'
        )
    if c == ['T', 'T']:
        return (
            'XX homozygote (endurance direction)',
            'Both copies have the premature stop codon — no functional alpha-actinin-3 protein. '
            'Overrepresented in elite endurance athletes (~18% of Europeans). The effect is real '
            'but small (~2-3% of elite-athlete performance variance); for recreational training, '
            'environmental factors (training, sleep, diet, age) dominate.'
        )
    return ('Uninterpretable genotype', f'Observed: {gt}')


def main(subject_id='alice'):
    # Prefer imputed parquet if present
    imputed = PROJECT_ROOT / 'standardized-genomes' / 'imputed' / f'{subject_id}.imputed.parquet'
    chip = PROJECT_ROOT / 'standardized-genomes' / f'{subject_id}.parquet'
    parquet_path = imputed if imputed.exists() else chip
    print(f'Using {"IMPUTED" if parquet_path is imputed else "CHIP"} parquet: {parquet_path.name}')

    df = pd.read_parquet(parquet_path)
    df['chrom'] = df['chrom'].astype(str)

    # Match by position (most reliable)
    row = df[(df['chrom'] == CHROM) & (df['pos'] == POS)]
    if row.empty:
        # Fallback to rsID
        row = df[df['rsid'] == RSID]

    if row.empty:
        print(f'ACTN3 {RSID} not found in {parquet_path.name}.')
        return

    r = row.iloc[0]
    a1 = r['a1'] if 'a1' in r.index else None
    a2 = r['a2'] if 'a2' in r.index else None
    if not a1 or not a2 or (isinstance(a1, float) and pd.isna(a1)):
        print(f'{RSID}: genotype missing (a1={a1}, a2={a2}).')
        return
    gt = f'{a1}{a2}'
    imputed_flag = bool(r.get('imputed', False))

    label, explanation = interpret(gt)

    print()
    print(f'=== ACTN3 rs1815739 ===')
    print(f'  Position: chr{CHROM}:{POS} (GRCh37)')
    print(f'  Genotype: {gt} ({"imputed" if imputed_flag else "directly measured"})')
    print(f'  Result:   {label}')
    print(f'  Meaning:  {explanation}')

    # Log finding
    profile = json.loads((PROJECT_ROOT / 'profiles' / f'{subject_id}.json').read_text(encoding='utf-8'))
    declared = profile.get('declared_ancestry')

    rec = {
        'subject_id': subject_id,
        'topic': 'trait_actn3_sprint_endurance',
        'claim': f'Athletic predisposition (ACTN3 rs1815739): {label} — genotype {gt}',
        'variants': [{
            'rsid': RSID, 'gene': 'ACTN3', 'chrom': CHROM, 'pos': POS,
            'ref': REF_ALLELE, 'alt': ALT_ALLELE,
            'genotype': gt, 'on_chip': not imputed_flag, 'imputed': imputed_flag,
        }],
        'effect': {
            'type': 'trait',
            'value': label,
            'explanation': explanation,
            'effect_size_note': (
                '~2-3% of elite-athlete performance variance. Undetectable at recreational level. '
                'Training, sleep, diet, and age dominate actual performance. Worth knowing but not '
                'actionable for training decisions.'
            ),
        },
        'cohort_ancestry': 'multi-ethnic',
        'subject_ancestry_match': 'match' if (declared and 'european' in (declared.get('top_level','').lower())) else 'unknown',
        'source_ids': ['snpedia:trait-curated', 'pmid:14523377'],  # Yang 2003 original ACTN3 paper
        'notes': (
            'Classical single-SNP athletic-trait association. Well-replicated genotype-frequency '
            'differences between elite endurance vs sprint/power athletes, but effect magnitude at '
            'the individual recreational level is modest. Reported for completeness; not a training '
            'prescription.'
        ),
        'clinvar_significance': None,
        'clinvar_review_stars': None,
        'pvalue': None,
        'n_cases': None,
        'n_controls': None,
        'odds_ratio': None,
        'replication_count': 20,
        'evidence_class': 'well_replicated_common_variant',
        'investigation_id': None,
    }
    tier, downgrade = compute_tier(rec)
    rec['tier_computed'] = tier
    rec['tier_computed_at'] = datetime.now(timezone.utc).isoformat()
    rec['tier_rule_version'] = TIER_RULE_VERSION
    rec['ancestry_downgrade'] = downgrade
    fid = ledger_io.append_finding(**rec)
    print(f'  Tier: {tier}')
    print(f'  Appended finding {fid}')


if __name__ == '__main__':
    subj = sys.argv[1] if len(sys.argv) > 1 else 'alice'
    main(subj)
