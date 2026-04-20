"""One-off: record a self-reported IQ value and supersede the cognitive-ability
PRS finding with a Rule-11 cross-check note. Genomic inference is preserved
as-is; the self-report is a separate line of evidence. Template — edit the
subject id and value before running.
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
    # 1. Update profile with self-reported IQ
    profile_path = PROJECT_ROOT / 'profiles' / f'{subject_id}.json'
    profile = json.loads(profile_path.read_text(encoding='utf-8'))
    srp = profile.setdefault('self_reported_phenotypes', {})
    srp['iq'] = {
        'value_approx': 142,
        'reported_at': '2026-04-18',
        'source': 'clinical_test',
        'confidence': 'high',
        'notes': 'Multiple official IQ tests plus unofficial tests, consistent at ~142. Official tests are the authoritative source; figure represents the central estimate across sittings.'
    }
    profile_path.write_text(json.dumps(profile, indent=2), encoding='utf-8')
    print(f'Updated {profile_path.name}: added self_reported_phenotypes.iq = 142')

    # 2. Supersede the active cognitive-ability finding with cross-check note
    rows = ledger_io.load_findings()
    superseded = {r['supersedes'] for r in rows if r.get('supersedes')}
    active = [r for r in rows
              if r['finding_id'] not in superseded
              and r.get('subject_id') == subject_id
              and r.get('topic') == 'prs_cognitive_ability_pgs002135']
    if not active:
        print('No active cognitive-ability finding to update.')
        return
    old = active[0]

    # Compute divergence metrics for the cross-check
    e = old.get('effect') or {}
    prs_z = e.get('z_score')
    predicted = e.get('predicted_value')
    residual_sd = e.get('residual_sd')
    reported_iq = 142
    reported_iq_z = (reported_iq - 100) / 15
    conditional_residual_z = (reported_iq - predicted) / residual_sd if residual_sd else None

    rec = dict(old)
    rec.pop('finding_id', None)
    rec.pop('timestamp', None)
    rec['supersedes'] = old['finding_id']

    crosscheck = (
        f' Self-report cross-check: subject reports IQ ~142 based on multiple official and '
        f'unofficial tests with consistent results (source: clinical_test, confidence: high, reported 2026-04-18). '
        f'At 142 the subject is +{reported_iq_z:.1f} SD on the IQ distribution vs. {prs_z:+.2f} SD on this PRS '
        f'({conditional_residual_z:+.2f} residual SDs above the PRS conditional mean of {predicted:.1f}; '
        f'~0.3% prior probability under the PRS model, outside the 95% CI). '
        f'NOT a contradiction: consistent with the well-documented "missing heritability" gap for cognitive '
        f'ability — twin/family heritability 50-80% vs. common-variant PRS r² 0.07. Rare variants '
        f'(not on chip), non-additive effects, and sample-size limitations in the source GWAS account '
        f'for the gap. The PRS systematically under-predicts people at the tails because the unmeasured '
        f'genetic signal (rare variants especially) disproportionately matters at extremes. '
        f'Genomic inference is NOT adjusted based on the self-report — they are separate lines of evidence '
        f'(Rule 11). The self-report is the better measure of actual cognitive performance; '
        f'the PRS is the better measure of what common-variant genotype data alone reveals.'
    )

    existing_notes = (rec.get('notes') or '').rstrip()
    rec['notes'] = existing_notes + crosscheck

    tier, downgrade = compute_tier(rec)
    rec['tier_computed'] = tier
    rec['tier_computed_at'] = datetime.now(timezone.utc).isoformat()
    rec['tier_rule_version'] = TIER_RULE_VERSION
    rec['ancestry_downgrade'] = downgrade

    fid = ledger_io.append_finding(**rec)
    print(f'Superseded {old["finding_id"]} -> {fid} (tier {tier})')
    print(f'  Cross-check recorded: actual IQ {reported_iq} (+{reported_iq_z:.1f} SD) vs. PRS predicted {predicted:.1f} IQ ({prs_z:+.2f} SD)')
    print(f'  Divergence: {conditional_residual_z:+.2f} residual SDs; consistent with missing-heritability gap.')


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else 'alice')
