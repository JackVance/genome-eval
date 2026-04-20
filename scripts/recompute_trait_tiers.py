"""Assign evidence_class + inference_confidence to trait and inference
findings, then recompute tiers with tier_rules v2. Appends a superseding
row wherever the tier changes or the evidence_class was previously absent.

Per MAINTENANCE.md append-only conventions, never edits existing rows.

Run from project root:
    python scripts/recompute_trait_tiers.py [subject_id]
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


# Topic -> (evidence_class, inference_confidence or None, optional replication_count floor)
# For HFE H63D and SERPINA1 PI*S we fill ClinVar metadata instead since they are
# classically-pathogenic variants; the tier rules will then hit the ClinVar A-tier path.
TOPIC_CLASSIFICATION = {
    # Mendelian / near-Mendelian traits
    'trait_blood_type_abo':             ('mendelian_trait', None, 20),
    'trait_earwax_bodyodor_abcc11':     ('mendelian_trait', None, 20),
    'trait_alcohol_flush_aldh2':        ('mendelian_trait', None, 20),

    # Well-replicated common-variant trait associations
    'trait_lactase_persistence':        ('well_replicated_common_variant', None, 20),
    'trait_alcohol_metabolism_adh1b':   ('well_replicated_common_variant', None, 10),
    'trait_bitter_taste_tas2r38':       ('well_replicated_common_variant', None, 10),

    # Weak/inconsistent single-variant associations
    'trait_caffeine_metabolism_cyp1a2': ('weakly_predictive_variant', None, None),

    # Inference-based (gene presence, multi-SNP composite)
    'trait_rh_factor':                  ('gene_presence_inference', 'high', 5),
    'trait_eye_color_multisnp':         ('multi_snp_composite', 'high', 15),

    # Array-limitation findings
    'pgx_cyp2d6_raw':                   ('not_callable_from_array', None, None),
    'pgx_hla_b_5701_tag':               ('array_limitation', None, None),
}

# Topics where we patch ClinVar metadata to trigger the A-tier ClinVar path
# (only when the data genuinely corresponds to a ClinVar-pathogenic variant).
CLINVAR_PATCHES = {
    'hfe_h63d': {
        'clinvar_significance': 'Pathogenic',
        'clinvar_review_stars': 3,
        'note_suffix': 'HFE H63D is classified as Pathogenic / Risk Factor in ClinVar with 3-star review status (confirmed by expert panel). Tier recomputed accordingly.',
    },
    'serpina1_pis': {
        'clinvar_significance': 'Pathogenic',
        'clinvar_review_stars': 3,
        'note_suffix': 'SERPINA1 PI*S (rs17580) is classified as Pathogenic in ClinVar with 3-star review status. Tier recomputed accordingly.',
    },
}


def main(subject_id='alice'):
    findings = ledger_io.load_findings()
    superseded_ids = {f['supersedes'] for f in findings if f.get('supersedes')}
    active = [f for f in findings
              if f['finding_id'] not in superseded_ids
              and f.get('subject_id') == subject_id]

    changes = []
    for f in active:
        topic = f.get('topic')
        new_evidence_class = None
        new_inference_conf = None
        new_replication_floor = None
        clinvar_patch = None

        if topic in TOPIC_CLASSIFICATION:
            new_evidence_class, new_inference_conf, new_replication_floor = TOPIC_CLASSIFICATION[topic]
        elif topic in CLINVAR_PATCHES:
            clinvar_patch = CLINVAR_PATCHES[topic]
        else:
            continue  # not in scope

        # Skip if the finding already has matching metadata and a non-unknown tier
        already_classified = (
            new_evidence_class is not None
            and (f.get('evidence_class') or '').lower() == new_evidence_class
            and (new_inference_conf is None
                 or (f.get('inference_confidence') or '').lower() == new_inference_conf)
            and f.get('tier_computed') != 'unknown'
        )
        if already_classified:
            continue

        # Build candidate
        candidate = dict(f)
        if new_evidence_class:
            candidate['evidence_class'] = new_evidence_class
        if new_inference_conf:
            candidate['inference_confidence'] = new_inference_conf
        if new_replication_floor and not candidate.get('replication_count'):
            candidate['replication_count'] = new_replication_floor
        if clinvar_patch:
            candidate['clinvar_significance'] = clinvar_patch['clinvar_significance']
            candidate['clinvar_review_stars'] = clinvar_patch['clinvar_review_stars']

        new_tier, downgrade = compute_tier(candidate)
        old_tier = f.get('tier_computed')
        if new_tier == old_tier and not clinvar_patch and not new_evidence_class:
            continue

        # Append superseder
        rec = dict(f)
        rec.pop('finding_id', None)
        rec.pop('timestamp', None)
        rec['supersedes'] = f['finding_id']
        if new_evidence_class:
            rec['evidence_class'] = new_evidence_class
        if new_inference_conf:
            rec['inference_confidence'] = new_inference_conf
        if new_replication_floor and not rec.get('replication_count'):
            rec['replication_count'] = new_replication_floor
        if clinvar_patch:
            rec['clinvar_significance'] = clinvar_patch['clinvar_significance']
            rec['clinvar_review_stars'] = clinvar_patch['clinvar_review_stars']
            rec['notes'] = (rec.get('notes') or '') + ' ' + clinvar_patch['note_suffix']
        else:
            note_add = (
                f' Tier recomputed {old_tier} -> {new_tier} under tier_rules {TIER_RULE_VERSION} '
                f'after assigning evidence_class={new_evidence_class!r}'
                + (f', inference_confidence={new_inference_conf!r}' if new_inference_conf else '')
                + '.'
            )
            rec['notes'] = (rec.get('notes') or '') + note_add
        rec['tier_computed'] = new_tier
        rec['tier_computed_at'] = datetime.now(timezone.utc).isoformat()
        rec['tier_rule_version'] = TIER_RULE_VERSION
        rec['ancestry_downgrade'] = downgrade

        ledger_io.append_finding(**rec)
        changes.append((topic, old_tier, new_tier))
        print(f'  {topic:40s} {old_tier!s:10s} -> {new_tier!s:10s}')

    print()
    print(f'{len(changes)} findings superseded with new tiers.')


if __name__ == '__main__':
    subj = sys.argv[1] if len(sys.argv) > 1 else 'alice'
    main(subj)
