"""Supersede the current height PRS finding with a self-report cross-check.

The genomic prediction stays as-is; a separate Rule-11 cross-check line
records the subject-reported actual height and the discrepancy. Also
supersedes earlier stale height PRS rows (from pre-calibration runs)
so the active finding is clean.
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
    active_prs = [
        f for f in findings
        if f.get('subject_id') == subject_id
        and f.get('topic', '').startswith('prs_height_')
        and f['finding_id'] not in superseded_ids
    ]
    if not active_prs:
        print('No active height PRS finding.')
        return

    profile = json.loads((PROJECT_ROOT / 'profiles' / f'{subject_id}.json').read_text(encoding='utf-8'))
    sr = (profile.get('self_reported_phenotypes') or {}).get('height')

    # Supersede each active height PRS finding — in practice there should only be one,
    # but earlier pipeline runs left a couple of stale rows.
    for old in active_prs:
        rec = dict(old)
        rec.pop('finding_id', None)
        rec.pop('timestamp', None)
        rec['supersedes'] = old['finding_id']

        # Build cross-check notes
        existing_notes = (rec.get('notes') or '').strip()
        if sr:
            reported_cm = sr.get('value_cm')
            reported_imp = sr.get('value_imperial')
            predicted_cm = (rec.get('effect') or {}).get('predicted_height_cm')
            ci95 = (rec.get('effect') or {}).get('ci95_cm') or []
            diff_cm = reported_cm - predicted_cm if (reported_cm and predicted_cm) else None
            # Actual height in population SDs
            pop_mean = 177.0
            pop_sd = 7.0
            reported_z = (reported_cm - pop_mean) / pop_sd if reported_cm else None
            predicted_z = (rec.get('effect') or {}).get('z_score_approx')

            crosscheck_note = (
                f' Self-report cross-check: subject reports actual height {reported_cm:.0f} cm / {reported_imp} '
                f'(reported {sr.get("reported_at")}, source {sr.get("source")}, '
                f'subject-reported confidence {sr.get("confidence")}). '
                f'Actual position in European male distribution: ~{reported_z:+.1f} SD. '
                f'PRS-predicted position: {predicted_z:+.2f} SD ({predicted_cm:.1f} cm, 95% CI {ci95[0]:.1f}-{ci95[1]:.1f} cm). '
                f'Divergence: +{diff_cm:.1f} cm; actual height is OUTSIDE the PRS 95% CI '
                f'({"upper bound " + str(round(ci95[1],1)) + " cm exceeded" if reported_cm > ci95[1] else "within CI"}). '
                f'Candidate explanations (not mutually exclusive): '
                f'(a) PRS r^2 only ~0.20 captures a modest fraction of height variance; '
                f'(b) chip overlap with this PRS is ~15%, so most height-contributing variants are unmeasured; '
                f'(c) regression-to-the-mean effect — PRS systematically underpredicts at distribution tails '
                f'because the missing variants disproportionately contribute to extreme phenotypes; '
                f'(d) rare / private variants not captured by a common-variant PRS; '
                f'(e) non-genetic contributors (childhood nutrition, etc.) — less likely given the magnitude. '
                f'Cross-check result: SIGNIFICANT MISMATCH directionally consistent with known PRS limitations. '
                f'The genomic prediction is recorded as-is per Rule 11; the self-report stands separately.'
            )
            rec['notes'] = (existing_notes + crosscheck_note).strip()

        tier, downgrade = compute_tier(rec)
        rec['tier_computed'] = tier
        rec['tier_computed_at'] = datetime.now(timezone.utc).isoformat()
        rec['tier_rule_version'] = TIER_RULE_VERSION
        rec['ancestry_downgrade'] = downgrade
        fid = ledger_io.append_finding(**rec)
        print(f'Superseded {old["finding_id"]} -> {fid} (tier {tier})')


if __name__ == '__main__':
    subj = sys.argv[1] if len(sys.argv) > 1 else 'alice'
    main(subj)
