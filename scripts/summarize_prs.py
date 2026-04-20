"""Print the cross-trait PRS table for a subject per SKILL.md Rule 10.1.

Reads active PRS findings from `ledger/findings.jsonl` and emits a markdown
table with one row per trait. Use this for session reviews, report generation,
or any time multiple PRS are being surfaced together.

Usage:
    python scripts/summarize_prs.py alice
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
sys.path.insert(0, str(SCRIPT_DIR))

import ledger_io
from run_prs import fmt_value, fmt_anchor, _get_unit_preference


def _suspect_with_reason(row: dict) -> tuple[bool, str]:
    """Return (is_suspect, reason). "Suspect" means the number shouldn't be
    trusted because of a pipeline / calibration issue — NOT because the
    underlying PGS has low r² for individual prediction. Low r² is an
    inherent limitation of the PRS, visible from the r² column itself, and
    doesn't need a separate SUSPECT tag.

    Suspect criteria:
    - `effect.suspect` is truthy (explicit flag set by superseder or runner).
    - |z| > 5 — outside the plausible empirical range for a common-variant
      PRS; almost always a scoring/calibration artifact.
    - Panel coverage < 10% — coordinate / strand mismatch between PGS weights
      and the imputed parquet; z-score derived from a tiny subset of the
      intended panel.
    - Notes explicitly include "SUSPECT NUMERIC OUTPUT" — the marker written
      by `supersede_suspect_prs.py` when a score is flagged post-hoc.
    """
    e = row.get('effect') or {}
    v = (row.get('variants') or [{}])[0]
    if e.get('suspect'):
        return True, str(e.get('suspect_reason', 'flagged suspect by runner'))
    z = e.get('z_score')
    if z is not None and abs(z) > 5:
        return True, (
            f'|z| = {abs(z):.2f} > 5 — outside the plausible empirical range '
            f'for a common-variant PRS; likely scoring/calibration artifact.'
        )
    cov = v.get('coverage_fraction')
    if cov is not None and cov < 0.10:
        return True, (
            f'Panel coverage only {cov:.1%} — most of the PGS weights did '
            f'not map to the imputed parquet (coordinate/strand mismatch). '
            f'The z-score is computed from a tiny subset of the intended panel.'
        )
    notes = row.get('notes') or ''
    if 'SUSPECT NUMERIC OUTPUT' in notes:
        return True, (
            'Flagged in notes as SUSPECT NUMERIC OUTPUT — see finding notes '
            'for the specific caveat.'
        )
    return False, ''


def _is_suspect(row: dict) -> bool:
    ok, _ = _suspect_with_reason(row)
    return ok


def main(subject_id='alice'):
    active = ledger_io.load_active_findings(subject_id=subject_id)
    prs = [r for r in active if (r.get('effect') or {}).get('type') == 'prs']
    if not prs:
        print(f'No active PRS findings for {subject_id}.')
        return

    prs.sort(key=lambda r: (r.get('effect') or {}).get('r_squared_from_paper') or 0, reverse=True)
    pref = _get_unit_preference(subject_id)

    print()
    n_suspect = sum(1 for r in prs if _is_suspect(r))
    suspect_note = f"  [{n_suspect} flagged as SUSPECT]" if n_suspect else ""
    print(f'PRS summary for subject `{subject_id}` ({len(prs)} active, display={pref}){suspect_note}')
    print()
    print('| Status | Trait (PGS id) | PRS z | %ile | Phenotype anchor | Predicted | 95% CI | r² | Naive (if r²=1) |')
    print('|---|---|---|---|---|---|---|---|---|')
    suspect_reasons = []
    for r in prs:
        e = r.get('effect') or {}
        v = (r.get('variants') or [{}])[0]
        trait = e.get('trait') or v.get('trait') or r.get('topic', '?')
        pgs = e.get('pgs_id') or v.get('pgs_id') or ''
        z = e.get('z_score')
        pctile = e.get('percentile')
        predicted = e.get('predicted_value')
        ci = e.get('ci95') or [None, None]
        r2 = e.get('r_squared_from_paper')
        sex = e.get('sex_assumed', '')

        units = ''
        anchor_str = 'n/a'
        naive_str = 'n/a'
        try:
            from prs_traits import TRAITS
            cfg = None
            for _, c in TRAITS.items():
                if c['pgs_id'] == pgs:
                    cfg = c
                    break
            if cfg:
                units = cfg.get('units', '')
                anchor = cfg.get('anchors', {}).get(sex) or cfg.get('anchors', {}).get('any')
                if anchor:
                    anchor_str = fmt_anchor(anchor['mean'], anchor['sd'], units, pref)
                    if z is not None:
                        naive = anchor['mean'] + anchor['sd'] * z
                        naive_str = fmt_value(naive, units, pref)
        except Exception:
            pass

        predicted_str = fmt_value(predicted, units, pref)
        ci_str = f'{fmt_value(ci[0], units, pref)} – {fmt_value(ci[1], units, pref)}' if ci and ci[0] is not None else 'n/a'
        z_str = f'{z:+.2f}' if z is not None else 'n/a'
        pctile_str = f'{pctile:.1f}' if pctile is not None else 'n/a'
        r2_str = f'{r2:.2f}' if r2 is not None else 'n/a'

        # Status column: visible marker for suspect results so they don't
        # slip past a casual reader of the table. Reasons are collected
        # below the table for traceability.
        suspect, reason = _suspect_with_reason(r)
        status = '**SUSPECT**' if suspect else 'ok'
        if suspect:
            suspect_reasons.append((trait, pgs, reason))

        print(f'| {status} | {trait} ({pgs}) | {z_str} | {pctile_str} | {anchor_str} | {predicted_str} | {ci_str} | {r2_str} | {naive_str} |')
    print()

    if suspect_reasons:
        print('### Suspect-result reasons')
        print()
        for trait, pgs, reason in suspect_reasons:
            print(f'- **{trait} ({pgs}):** {reason}')
        print()


if __name__ == '__main__':
    subj = sys.argv[1] if len(sys.argv) > 1 else 'alice'
    main(subj)
