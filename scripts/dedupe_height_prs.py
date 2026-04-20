"""Chain the 5 active height PRS rows into a single canonical active one.

Keep the most recent well-calibrated finding; mark all others superseded
by pointing each one back via a chain of supersede records.
"""
import json
import sys
import uuid
from datetime import datetime, timezone
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
sys.path.insert(0, str(SCRIPT_DIR))

import ledger_io


def main():
    findings = ledger_io.load_findings()
    superseded = {f['supersedes'] for f in findings if f.get('supersedes')}
    active_h = [
        f for f in findings
        if f['finding_id'] not in superseded
        and f.get('topic', '').startswith('prs_height_')
    ]
    if len(active_h) <= 1:
        print('Already clean.')
        return

    # Sort by timestamp so we chain in order
    active_h.sort(key=lambda f: f.get('timestamp', ''))
    canonical = active_h[-1]
    print(f'Canonical: {canonical["finding_id"]} ({canonical["claim"][:60]}...)')
    print(f'Superseding {len(active_h) - 1} earlier rows into the canonical chain.')

    # Each earlier row gets marked superseded by appending a tombstone row
    # pointing to itself as the predecessor of the next one in chronological order.
    for i in range(len(active_h) - 1):
        old = active_h[i]
        next_row = active_h[i + 1]
        tombstone = {
            'subject_id': old['subject_id'],
            'topic': old['topic'],
            'supersedes': old['finding_id'],
            'claim': f'SUPERSEDED (pipeline iteration cleanup): rolled into {next_row["finding_id"]} (later calibration run)',
            'variants': [],
            'effect': None,
            'cohort_ancestry': old.get('cohort_ancestry'),
            'subject_ancestry_match': old.get('subject_ancestry_match'),
            'source_ids': old.get('source_ids', []),
            'notes': (
                f'Tombstone row. The original finding was from an intermediate pipeline run '
                f'(calibration or coverage state). Superseded by {next_row["finding_id"]} which '
                f'uses the final calibrated pipeline. No tiering claim; refer to the active chain head.'
            ),
            'clinvar_significance': None,
            'clinvar_review_stars': None,
            'pvalue': None,
            'n_cases': None,
            'n_controls': None,
            'odds_ratio': None,
            'replication_count': None,
            'tier_computed': 'unknown',
            'tier_computed_at': datetime.now(timezone.utc).isoformat(),
            'tier_rule_version': 'v2',
            'ancestry_downgrade': False,
            'evidence_class': 'pipeline_intermediate',
        }
        fid = ledger_io.append_finding(**tombstone)
        print(f'  {old["finding_id"]} -> tombstone {fid}')


if __name__ == '__main__':
    main()
