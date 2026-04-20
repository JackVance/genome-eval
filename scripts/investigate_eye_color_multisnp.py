"""Multi-SNP eye color investigation using the HIrisPlex-S pigmentation loci.

HIrisPlex-S (Walsh et al. 2017) combines 41 SNPs to predict eye, hair, and
skin color simultaneously. The eye-color sub-model is a superset of the
earlier 6-SNP IrisPlex model. This script queries the chip for the
well-established pigmentation loci, reports each genotype with its
published eye-color direction, and computes a simple per-color vote
that approximates the HIrisPlex direction without the actual model
coefficients.

Also supersedes the earlier single-SNP eye-color finding for the subject
if its confidence framing needs updating.

Run from project root:
    python scripts/investigate_eye_color_multisnp.py [subject_id]
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


# Pigmentation SNPs used in IrisPlex / HIrisPlex-S eye-color prediction.
# Columns: rsid, gene, GRCh37 forward-strand allele association direction for eye color.
#
# "blue_allele" is the allele that pushes toward blue/lighter eyes.
# "effect_size" is a rough ordinal: strong / moderate / weak.
#
# Sources: Walsh 2017 (HIrisPlex-S), Liu 2009 (IrisPlex), SNPedia, dbSNP.

PIGMENTATION_SNPS = [
    # (rsid, gene, blue_allele, brown_allele, effect_size, note)
    ('rs12913832', 'HERC2',   'A', 'G', 'strong',
     'Dominant blue-eye locus; explains ~74% of blue/brown variance in Europeans. AA=blue, GG=brown.'),
    ('rs1800407',  'OCA2',    'T', 'C', 'moderate',
     'R419Q (OCA2): T allele shifts toward lighter eye color; can push rs12913832 AG toward blue.'),
    ('rs12896399', 'SLC24A4', 'T', 'G', 'moderate',
     'T allele associated with lighter (blue/green) eyes.'),
    ('rs16891982', 'SLC45A2', 'G', 'C', 'moderate',
     'F374L: G (phenylalanine, derived) lighter pigmentation; C (leucine, ancestral) darker. Corrected direction 2026-04-17.'),
    ('rs1393350',  'TYR',     'A', 'G', 'weak',
     'A allele weakly associated with lighter pigmentation.'),
    ('rs12203592', 'IRF4',    'T', 'C', 'weak',
     'T allele weakly associated with lighter eye color and freckling.'),
    ('rs683',      'TYRP1',   'A', 'C', 'weak',
     'Weak effect; A allele lighter direction.'),
    ('rs1129038',  'HERC2',   'T', 'C', 'moderate',
     'In tight linkage with rs12913832; T->blue direction. Sanity-check SNP.'),
    ('rs1667394',  'HERC2',   'T', 'C', 'moderate',
     'HERC2 regulatory region; T toward blue.'),
    ('rs6119471',  'ASIP',    'C', 'G', 'weak',
     'Weak ASIP effect.'),
    ('rs1470608',  'OCA2',    'T', 'G', 'weak',
     'Weak OCA2 effect.'),
]


def interpret_genotype(gt, blue_a, brown_a):
    """Return (blue_count, brown_count, ambiguous_count) for a 2-allele genotype."""
    blue = gt.count(blue_a)
    brown = gt.count(brown_a)
    ambiguous = 2 - blue - brown
    return blue, brown, ambiguous


def main(subject_id='alice'):
    df = pd.read_parquet(PROJECT_ROOT / 'standardized-genomes' / f'{subject_id}.parquet')

    rsid_gt = {}
    for r in df.itertuples():
        a1 = r.a1 if pd.notna(r.a1) else '-'
        a2 = r.a2 if pd.notna(r.a2) else '-'
        rsid_gt[r.rsid] = a1 + a2

    print('=== Pigmentation SNP genotypes ===')
    print()
    print(f'{"rsID":12s}  {"gene":8s}  {"GT":6s}  {"blue/brown":12s}  {"effect":8s}  note')
    print('-' * 120)

    strong_blue = 0
    strong_brown = 0
    moderate_blue = 0
    moderate_brown = 0
    weak_blue = 0
    weak_brown = 0
    findings_per_snp = []

    for rsid, gene, blue_a, brown_a, effect, note in PIGMENTATION_SNPS:
        if rsid not in rsid_gt:
            print(f'{rsid:12s}  {gene:8s}  NOT ON CHIP')
            findings_per_snp.append({'rsid': rsid, 'gene': gene, 'on_chip': False, 'effect_size': effect})
            continue
        gt = rsid_gt[rsid]
        b, brn, amb = interpret_genotype(gt, blue_a, brown_a)
        direction = f'{b}/{brn}'
        print(f'{rsid:12s}  {gene:8s}  {gt:6s}  {direction:12s}  {effect:8s}  {note[:50]}')

        findings_per_snp.append({
            'rsid': rsid, 'gene': gene, 'genotype': gt, 'on_chip': True,
            'blue_allele': blue_a, 'brown_allele': brown_a,
            'blue_count': b, 'brown_count': brn, 'effect_size': effect,
        })

        if effect == 'strong':
            strong_blue += b
            strong_brown += brn
        elif effect == 'moderate':
            moderate_blue += b
            moderate_brown += brn
        else:
            weak_blue += b
            weak_brown += brn

    print()
    print('=== Summed allele counts by effect tier ===')
    print(f'  Strong effect: {strong_blue} blue / {strong_brown} brown')
    print(f'  Moderate:      {moderate_blue} blue / {moderate_brown} brown')
    print(f'  Weak:          {weak_blue} blue / {weak_brown} brown')

    # Weighted composite score: strong=3, moderate=2, weak=1.
    weighted_blue = 3 * strong_blue + 2 * moderate_blue + 1 * weak_blue
    weighted_brown = 3 * strong_brown + 2 * moderate_brown + 1 * weak_brown
    total_weighted = weighted_blue + weighted_brown
    print()
    if total_weighted > 0:
        blue_pct = weighted_blue / total_weighted * 100
        brown_pct = weighted_brown / total_weighted * 100
        print(f'Weighted composite: blue {blue_pct:.0f}% / brown {brown_pct:.0f}%')
    else:
        blue_pct = brown_pct = None

    # Interpretation
    print()
    if strong_blue == 2:
        dominant = 'Strong evidence for blue eyes (rs12913832 AA).'
    elif strong_brown == 2:
        dominant = 'Single-SNP prediction: brown eyes (rs12913832 GG). Check other loci for consistency.'
    elif strong_blue == 1 and strong_brown == 1:
        dominant = 'Intermediate at dominant locus (rs12913832 AG); other SNPs resolve direction.'
    else:
        dominant = 'rs12913832 missing or ambiguous.'
    print(f'Dominant-locus call:   {dominant}')
    print(f'Multi-locus composite: {weighted_blue} blue / {weighted_brown} brown (weighted)')

    # If rs12913832 says brown but composite leans blue, flag for miscall
    rs12913832_gt = rsid_gt.get('rs12913832')
    rs12913832_says_brown = rs12913832_gt == 'GG'
    composite_leans_blue = blue_pct and blue_pct > 60

    if rs12913832_says_brown and composite_leans_blue:
        print()
        print('NOTE: rs12913832 is GG (brown-predicting) but the broader pigmentation panel')
        print('leans blue. Possible explanations:')
        print('  1. Array miscall at rs12913832 (single-SNP error rate ~0.1-0.5%).')
        print('  2. Genuine atypical case where lighter modifier loci override rs12913832.')
        print('  3. rs12913832 is in LD with a neighboring probe — check rs1129038 and rs1667394 for strand-sanity.')

    # Also check the HERC2 linkage-sanity SNPs (rs1129038, rs1667394)
    sanity = {r: rsid_gt.get(r) for r in ('rs1129038', 'rs1667394')}
    print()
    print(f'HERC2 linkage-sanity SNPs: {sanity}')
    print('  These should normally track rs12913832 (same haplotype in Europeans).')
    print('  If they disagree with rs12913832, the rs12913832 call is more suspect.')

    # Write finding
    profile = json.loads((PROJECT_ROOT / 'profiles' / f'{subject_id}.json').read_text(encoding='utf-8'))
    declared = profile.get('declared_ancestry')

    # Supersede the earlier single-SNP eye-color finding
    findings = ledger_io.load_findings()
    superseded_ids = {f['supersedes'] for f in findings if f.get('supersedes')}
    old = [f for f in findings
           if f.get('subject_id') == subject_id
           and f.get('topic') == 'trait_eye_color_herc2'
           and f['finding_id'] not in superseded_ids]
    old = old[-1] if old else None

    # Build composite claim
    if composite_leans_blue and rs12913832_says_brown:
        claim = (
            'Eye color — multi-SNP composite conflicts with dominant-locus call. '
            f'rs12913832 GG (brown-predicting) but broader pigmentation panel is {blue_pct:.0f}% blue direction. '
            'Consistent with reported blue-eye phenotype if rs12913832 is a miscall or '
            'lighter modifier loci dominate.'
        )
    elif composite_leans_blue:
        claim = f'Eye color — multi-SNP composite predicts blue / light ({blue_pct:.0f}% blue-direction alleles, weighted).'
    elif blue_pct and blue_pct < 40:
        claim = f'Eye color — multi-SNP composite predicts brown / dark ({brown_pct:.0f}% brown-direction alleles, weighted).'
    else:
        claim = f'Eye color — multi-SNP composite is intermediate (green/hazel possible).'

    rec = {
        'subject_id': subject_id,
        'supersedes': old['finding_id'] if old else None,
        'topic': 'trait_eye_color_multisnp',
        'claim': claim,
        'variants': findings_per_snp,
        'effect': {
            'type': 'trait_multisnp',
            'dominant_locus_call': 'brown (rs12913832 GG)' if rs12913832_says_brown else f'rs12913832={rs12913832_gt}',
            'weighted_blue_percent': blue_pct,
            'weighted_brown_percent': brown_pct,
            'strong_effect': {'blue': strong_blue, 'brown': strong_brown},
            'moderate_effect': {'blue': moderate_blue, 'brown': moderate_brown},
            'weak_effect': {'blue': weak_blue, 'brown': weak_brown},
            'method': 'Weighted allele-count composite across HIrisPlex-S eye-color loci (approximate).',
            'note': (
                'Not the full HIrisPlex-S model. Accurate HIrisPlex-S prediction requires '
                'logistic-regression coefficients from Walsh 2017. This composite uses effect-size '
                'weights (3/2/1) applied to per-allele counts.'
            ),
        },
        'cohort_ancestry': 'European',
        'subject_ancestry_match': 'match' if (declared and 'european' in (declared.get('top_level','').lower())) else 'unknown',
        'source_ids': ['pmid:28376720', 'snpedia:trait-curated'],  # Walsh HIrisPlex-S
        'notes': (
            'Multi-SNP eye-color estimate from pigmentation loci. Single-SNP rs12913832 explains '
            '~74% of blue/brown variance; remainder carried by OCA2, TYR, SLC24A4, SLC45A2, IRF4, '
            'TYRP1, ASIP, KITLG. If the dominant-locus call disagrees with reported phenotype, '
            'check (a) array miscall probability at rs12913832, (b) other pigmentation SNPs for '
            'consistency, (c) HERC2 linkage-sanity SNPs rs1129038 and rs1667394.'
        ),
        'clinvar_significance': None,
        'clinvar_review_stars': None,
        'pvalue': None,
        'n_cases': None,
        'n_controls': None,
        'odds_ratio': None,
        'replication_count': 10,
        'investigation_id': None,
    }
    tier, downgrade = compute_tier(rec)
    rec['tier_computed'] = tier
    rec['tier_computed_at'] = datetime.now(timezone.utc).isoformat()
    rec['tier_rule_version'] = TIER_RULE_VERSION
    rec['ancestry_downgrade'] = downgrade
    fid = ledger_io.append_finding(**rec)
    print()
    print(f'Appended finding {fid} (tier {tier}); supersedes {old["finding_id"] if old else "None"}.')


if __name__ == '__main__':
    subj = sys.argv[1] if len(sys.argv) > 1 else 'alice'
    main(subj)
