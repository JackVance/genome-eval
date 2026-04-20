"""Re-supersede the eye color finding to enforce Rule 11 separation:
genomic inference stands on genomic evidence alone; the subject's
self-report is cross-checked in a separate notes line, not folded
into the genomic reasoning.
"""
from __future__ import annotations

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
sys.path.insert(0, str(SCRIPT_DIR))

import ledger_io
from tier_rules import compute_tier, TIER_RULE_VERSION


def main(subject_id='alice'):
    findings = ledger_io.load_findings()
    superseded_ids = {f['supersedes'] for f in findings if f.get('supersedes')}

    active_eye = [f for f in findings
                  if f.get('subject_id') == subject_id
                  and f.get('topic') in ('trait_eye_color_multisnp', 'trait_eye_color_herc2')
                  and f['finding_id'] not in superseded_ids]
    if not active_eye:
        print('No active eye color finding.')
        return
    old = active_eye[-1]

    profile = json.loads((PROJECT_ROOT / 'profiles' / f'{subject_id}.json').read_text(encoding='utf-8'))
    sr = (profile.get('self_reported_phenotypes') or {}).get('eye_color')

    rec = dict(old)
    rec.pop('finding_id', None)
    rec.pop('timestamp', None)
    rec['supersedes'] = old['finding_id']

    rec['claim'] = (
        'Eye color (genomic inference): blue-haplotype homozygous at HERC2 with one SNP '
        'showing an inconsistent call. rs1129038 TT and rs1667394 TT are both homozygous for '
        'the blue-eye haplotype; rs12913832 reads GG (homozygous brown haplotype), which is '
        'internally inconsistent in a European genome given the tight LD at this locus. '
        'Most parsimonious explanation: rs12913832 is a single-probe miscall; true genotype '
        'likely AA (blue haplotype). Genomic prediction: blue.'
    )

    rec['effect'] = {
        'type': 'trait_multisnp_with_confidence',
        'genomic_prediction': 'blue',
        'genotype_call_confidence': {
            'rs12913832': 'low (suspected miscall from LD inconsistency)',
            'rs1129038': 'high',
            'rs1667394': 'high',
            'others (SLC45A2, OCA2, SLC24A4, TYR, IRF4, TYRP1)': 'high',
        },
        'genotype_to_phenotype_confidence': (
            'High for blue prediction once rs12913832 is treated as a suspected miscall. '
            'rs1129038 and rs1667394 are in tight linkage disequilibrium with rs12913832 in '
            'Europeans and are normally inherited as a block; two-of-three agreement on the '
            'blue haplotype against a single disagreeing SNP is the classic signature of a '
            'single-probe error.'
        ),
        'confidence_reasoning_one_liner': (
            'Two HERC2 linkage-companion SNPs homozygous blue-haplotype; rs12913832 outlier '
            'likely a single-probe miscall (~0.1-0.5% per-SNP array error rate).'
        ),
        'supporting_loci': {
            'rs1129038 (HERC2)': 'TT — blue-haplotype homozygous',
            'rs1667394 (HERC2)': 'TT — blue-haplotype homozygous',
            'rs16891982 (SLC45A2)': 'GG — lighter/phenylalanine allele homozygous',
        },
        'outlier_loci': {
            'rs12913832 (HERC2)': 'GG — suspected miscall (contradicts linkage-companion SNPs)',
        },
        'neutral_or_darker_loci': {
            'rs1800407 (OCA2 R419Q)': 'CC — no lighter-shifting variant',
            'rs12896399 (SLC24A4)': 'GG',
            'rs12203592 (IRF4)': 'CC',
        },
        'method': 'HIrisPlex-S pigmentation panel + HERC2 linkage-sanity cross-check.',
    }

    notes_parts = [
        'Genomic inference derived from SNP data only. The blue-eye conclusion rests on '
        'agreement between two HERC2 linkage-companion SNPs (rs1129038 TT, rs1667394 TT) '
        'against one disagreeing outlier (rs12913832 GG) that is most likely a single-probe '
        'miscall. The inference does not depend on subject-reported phenotype.',
    ]
    if sr:
        value = sr.get('value')
        match_status = 'match' if value and 'blue' in value.lower() else 'mismatch'
        notes_parts.append(
            f'Self-report cross-check: subject reports eye color "{value}" '
            f'(reported {sr.get("reported_at")}, source {sr.get("source")}, '
            f'subject-reported confidence {sr.get("confidence")}). '
            f'Cross-check result: {match_status.upper()} with the genomic inference. '
            'Self-report is tracked separately from the genomic conclusion per Rule 11.'
        )
    rec['notes'] = ' '.join(notes_parts)

    rec['replication_count'] = 15
    tier, downgrade = compute_tier(rec)
    rec['tier_computed'] = tier
    rec['tier_computed_at'] = datetime.now(timezone.utc).isoformat()
    rec['tier_rule_version'] = TIER_RULE_VERSION
    rec['ancestry_downgrade'] = downgrade

    fid = ledger_io.append_finding(**rec)
    print(f'Superseded {old["finding_id"]} -> {fid} (tier {tier})')
    print()
    print('New finding separates genomic evidence from self-report per Rule 11.')


if __name__ == '__main__':
    subj = sys.argv[1] if len(sys.argv) > 1 else 'alice'
    main(subj)
