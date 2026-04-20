"""Generic PRS runner — dispatches to the right trait-specific config.

Trait configs live in `scripts/prs_traits.py` as a dict keyed by trait name.
Each config contains: PGS ID, trait label, r^2 source, population anchors
(if any), and narrative notes.

Usage:
    python scripts/run_prs.py <subject> <trait>
    # e.g. python scripts/run_prs.py alice height_yengo2022
    # e.g. python scripts/run_prs.py alice educational_attainment
"""
from __future__ import annotations

import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
sys.path.insert(0, str(SCRIPT_DIR))

import ledger_io
from tier_rules import compute_tier, TIER_RULE_VERSION
from prs_pipeline import load_weights, compute_score, zscore_and_percentile, COMPLEMENT
from prs_traits import TRAITS

WEIGHTS_DIR = PROJECT_ROOT / 'reference' / 'prs_weights'


def _get_unit_preference(subject_id: str) -> str:
    """Return 'us_customary' (default) or 'si' based on profile override or env var.

    Lookup order: GENOME_EVAL_UNITS env var > profile.display_preferences.units > default.
    """
    import os
    env = os.environ.get('GENOME_EVAL_UNITS')
    if env in ('us_customary', 'us', 'si'):
        return 'us_customary' if env.startswith('us') else 'si'
    profile_path = PROJECT_ROOT / 'profiles' / f'{subject_id}.json'
    if profile_path.exists():
        try:
            prof = json.loads(profile_path.read_text(encoding='utf-8'))
            pref = (prof.get('display_preferences') or {}).get('units')
            if pref in ('us_customary', 'si'):
                return pref
        except Exception:
            pass
    return 'us_customary'


def fmt_value(value, units: str, preference: str) -> str:
    """Format a phenotype value in the user's preferred unit system.

    Storage in the ledger stays in the canonical unit (typically SI) regardless.
    """
    if value is None:
        return 'n/a'
    if preference == 'us_customary':
        if units == 'cm':
            total_inches = value / 2.54
            feet = int(total_inches // 12)
            inches = total_inches - feet * 12
            return f'{feet}\'{inches:.1f}"'
        if units == 'kg':
            return f'{value / 0.453592:.1f} lb'
        if units in ('C', '°C'):
            return f'{value * 9 / 5 + 32:.1f} °F'
    return f'{value:.2f} {units}' if isinstance(value, float) else f'{value} {units}'


def fmt_anchor(mean, sd, units: str, preference: str) -> str:
    """Format the 'mean ± SD per SD' anchor."""
    if preference == 'us_customary' and units == 'cm':
        m_in = mean / 2.54
        m_ft = int(m_in // 12)
        m_rem = m_in - m_ft * 12
        sd_in = sd / 2.54
        return f'{m_ft}\'{m_rem:.1f}" ± {sd_in:.1f}" per SD'
    return f'{mean:g} {units} ± {sd:g} {units} per SD'


def theoretical_from_afs_json(details, af_lookup: dict):
    """Legacy JSON-cache path: keyed by rsid -> {allele: freq}."""
    contributing = []
    mean = 0.0
    var = 0.0
    for d in details:
        rsid = d.get('rsid')
        ea = d.get('effect_allele_panel')
        af = None
        entry = af_lookup.get(rsid, {}) if rsid else {}
        if entry:
            af = entry.get(ea)
            if af is None:
                af = entry.get(COMPLEMENT.get(ea))
        if af is None:
            continue
        d['eur_af'] = af
        contributing.append(d)
        beta = d['effect_weight']
        mean += 2 * af * beta
        var += 2 * af * (1 - af) * beta * beta
    return mean, math.sqrt(var) if var > 0 else 0.0, contributing


def theoretical_from_afs_parquet(details, af_df: pd.DataFrame):
    """Local-parquet path: keyed by (chrom, pos) with ref/alt/eur_af columns.

    Matches (chrom, pos) and validates that the effect_allele_panel matches the
    parquet's ref or alt (direct strand) or its complement (flipped strand). The
    AF returned is frequency of the effect allele (not always the ALT).
    """
    # Build lookup keyed by (chrom_str, pos)
    af_df = af_df.copy()
    af_df['chrom'] = af_df['chrom'].astype(str)
    key_idx = {(r.chrom, r.pos): (r.ref, r.alt, r.eur_af) for r in af_df.itertuples(index=False)}

    contributing = []
    mean = 0.0
    var = 0.0
    for d in details:
        key = (str(d.get('chr')), d.get('pos'))
        entry = key_idx.get(key)
        if entry is None:
            continue
        ref, alt, alt_af = entry
        ea = d.get('effect_allele_panel')
        # Determine effect-allele frequency (ALT AF if ea matches alt, else 1-alt_af if ea matches ref)
        if ea == alt:
            af = alt_af
        elif ea == ref:
            af = 1.0 - alt_af
        elif COMPLEMENT.get(ea) == alt:
            af = alt_af
        elif COMPLEMENT.get(ea) == ref:
            af = 1.0 - alt_af
        else:
            continue
        d['eur_af'] = af
        contributing.append(d)
        beta = d['effect_weight']
        mean += 2 * af * beta
        var += 2 * af * (1 - af) * beta * beta
    return mean, math.sqrt(var) if var > 0 else 0.0, contributing


def theoretical_uniform(betas, assumed_af: float = 0.35):
    s = pd.Series(betas)
    mean = float((2 * assumed_af * s).sum())
    var = float((2 * assumed_af * (1 - assumed_af) * (s ** 2)).sum())
    return mean, math.sqrt(var)


def predicted_value(z: float, trait_cfg: dict, sex: str):
    """Return (predicted, ci95_low, ci95_high, residual_sd) in the trait's units."""
    anchors = trait_cfg.get('anchors', {})
    anchor = anchors.get(sex) or anchors.get('any')
    if not anchor:
        return None, None, None, None
    r = math.sqrt(trait_cfg['r_squared'])
    predicted = anchor['mean'] + r * anchor['sd'] * z
    residual_sd = anchor['sd'] * math.sqrt(1 - trait_cfg['r_squared'])
    return predicted, predicted - 1.96 * residual_sd, predicted + 1.96 * residual_sd, residual_sd


def build_self_report_crosscheck(
    trait_cfg: dict,
    profile: dict,
    predicted,
    ci_low,
    ci_high,
    residual_sd,
    anchor,
) -> tuple[str, dict]:
    """Return (notes_block, structured_effect_fields) for the self-report cross-check.

    Guardrail per Rules 10/11: if the trait config declares a `self_report_key`,
    this function ALWAYS returns a non-empty block — MATCH, MISMATCH, or an
    explicit "no self-report on file" statement. You cannot silently skip it.

    If the trait config does not declare a self_report_key, returns ("", {}) —
    the caller decides whether to still note the omission. For PRS traits
    this should not happen; self_report_key should be set on every trait in
    prs_traits.py. If it is missing, a loud warning is emitted so the gap is
    visible rather than silently swallowed.
    """
    if 'self_report_key' not in trait_cfg:
        # Missing entirely → warn; config author forgot to declare.
        print(
            f'[WARN] Trait {trait_cfg.get("topic")} has no self_report_key in '
            f'prs_traits.py. Rule-10/11 cross-check will be absent from the '
            f'ledger row. Add a self_report_key (or set it to None explicitly) '
            f'to silence this warning.'
        )
        return '', {}
    sr_key = trait_cfg.get('self_report_key')
    if sr_key is None:
        # Explicitly opted out — genuinely inapplicable (e.g., binary disease
        # with no self-reported status). Silent skip.
        return '', {}

    sr_all = profile.get('self_reported_phenotypes') or {}
    sr = sr_all.get(sr_key)
    if not sr:
        block = (
            f' Self-report cross-check: no self-reported value on file for '
            f'`{sr_key}` in profiles/{profile.get("subject_id", "?")}.json; '
            f'cross-check not possible. To enable, add a '
            f'`self_reported_phenotypes.{sr_key}` entry per Rule 11.'
        )
        structured = {
            'self_report_crosscheck': {
                'status': 'no_self_report',
                'self_report_key': sr_key,
            },
        }
        return block, structured

    value_field = trait_cfg.get('self_report_value_field', 'value')
    sr_value = sr.get(value_field)
    if sr_value is None:
        # Fall back to a couple of common alternate keys
        for alt in ('value', 'value_approx', 'value_cm', 'value_years'):
            if sr.get(alt) is not None:
                sr_value = sr[alt]
                value_field = alt
                break
    if sr_value is None or predicted is None or ci_low is None or ci_high is None:
        block = (
            f' Self-report cross-check: self-report present for `{sr_key}` but '
            f'could not be numerically compared (value_field=`{value_field}`, '
            f'sr_value={sr_value!r}, predicted={predicted!r}). Cross-check not '
            f'possible — investigate the profile entry or the trait config.'
        )
        structured = {
            'self_report_crosscheck': {
                'status': 'uncomparable',
                'self_report_key': sr_key,
                'self_report_value_field': value_field,
                'self_report_value': sr_value,
            },
        }
        return block, structured

    conv = trait_cfg.get('self_report_unit_conversion', 1.0)
    sr_value_ledger_units = sr_value * conv

    within_ci = ci_low <= sr_value_ledger_units <= ci_high
    gap = sr_value_ledger_units - predicted
    gap_residual_sd = gap / residual_sd if residual_sd else float('nan')
    pheno_sd = (anchor or {}).get('sd')
    gap_population_sd = gap / pheno_sd if pheno_sd else float('nan')

    source = sr.get('source', 'unspecified source')
    confidence = sr.get('confidence', 'unspecified confidence')
    reported_at = sr.get('reported_at', 'unspecified date')
    display_value = sr.get('value_imperial') or sr.get('value') or f'{sr_value} ({value_field})'

    if within_ci:
        status = 'MATCH'
        interp = (
            f'Self-report sits within the PRS 95% CI. Note the CI width reflects '
            f'r^2={trait_cfg["r_squared"]:.2f} ({(1 - trait_cfg["r_squared"]):.0%} of '
            f'variance unexplained) — "match" here is consistent with the PRS but '
            f'loosely constraining.'
        )
    else:
        status = 'MISMATCH'
        interp = (
            f'Self-report sits outside the PRS 95% CI. Gap = {gap:+.2f} '
            f'{trait_cfg.get("units", "")} ({gap_residual_sd:+.2f} residual SDs '
            f'relative to the PRS point estimate; {gap_population_sd:+.2f} '
            f'population SDs from the mean). At this r^2, common-variant PRS '
            f'systematically underestimates tails (regression-to-the-mean shrinks '
            f'extreme true values toward the population mean in the *prediction* — '
            f'not in the person). This is a documented tail limitation, not a '
            f'PRS error and not a self-report error. Per Rule 11 genomic inference '
            f'is preserved as-is; the mismatch itself is a finding about the '
            f'PRS-at-tails.'
        )

    block = (
        f' Self-report cross-check: subject reports {display_value} for `{sr_key}` '
        f'(source {source}, confidence {confidence}, reported {reported_at}). '
        f'PRS-predicted: {predicted:.2f} {trait_cfg.get("units", "")} (95% CI '
        f'{ci_low:.2f} - {ci_high:.2f}). Cross-check result: {status}. {interp}'
    )
    structured = {
        'self_report_crosscheck': {
            'status': status.lower(),
            'self_report_key': sr_key,
            'self_report_value_field': value_field,
            'self_report_value': sr_value_ledger_units,
            'self_report_value_raw': sr_value,
            'self_report_source': source,
            'self_report_confidence': confidence,
            'self_report_reported_at': reported_at,
            'gap': gap,
            'gap_residual_sd': gap_residual_sd,
            'gap_population_sd': gap_population_sd,
            'within_95ci': within_ci,
        },
    }
    return block, structured


def run(subject_id: str, trait_key: str):
    if trait_key not in TRAITS:
        print(f'Unknown trait: {trait_key}')
        print(f'Available: {", ".join(TRAITS.keys())}')
        sys.exit(1)
    cfg = TRAITS[trait_key]
    pgs_id = cfg['pgs_id']

    weights_path = WEIGHTS_DIR / f'{pgs_id}_hmPOS_GRCh37.txt.gz'
    if not weights_path.exists():
        print(f'Weight file missing: {weights_path}')
        print(f'Run: python scripts/prs_download.py {pgs_id}')
        sys.exit(1)

    # Load subject data (prefer imputed)
    imputed = PROJECT_ROOT / 'standardized-genomes' / 'imputed' / f'{subject_id}.imputed.parquet'
    chip = PROJECT_ROOT / 'standardized-genomes' / f'{subject_id}.parquet'
    subject_parquet = imputed if imputed.exists() else chip
    print(f'Subject: {subject_parquet.name}')
    subject_df = pd.read_parquet(subject_parquet)
    subject_df['chrom'] = subject_df['chrom'].astype(str)

    profile = json.loads((PROJECT_ROOT / 'profiles' / f'{subject_id}.json').read_text(encoding='utf-8'))
    sex_stats = (profile.get('parse_stats') or {}).get('inferred_sex', 'unknown').lower()
    sex = 'male' if ('male' in sex_stats and 'female' not in sex_stats) else (
        'female' if 'female' in sex_stats else 'unknown'
    )
    declared = profile.get('declared_ancestry') or {}
    ancestry_match = 'match' if 'european' in (declared.get('top_level', '').lower()) else 'unknown'

    # Compute score
    weights, meta = load_weights(weights_path)
    print(f'PGS: {pgs_id} ({cfg["label"]}) — {len(weights):,} variants')
    result = compute_score(weights, subject_df)
    print(f'  panel={result["n_panel"]:,}  on_chip={result["n_on_chip"]:,}  '
          f'contributing={result["n_contributing"]:,}  coverage={result["coverage_fraction"]:.1%}')
    print(f'  strand_flips={result["n_strand_flip"]}  palindromic_skipped={result["n_palindromic_skipped"]}  '
          f'allele_mismatch={result["n_allele_mismatch_after_flip"]}')

    # Calibration: prefer empirical reference distribution (applies the PRS to
    # 1000G EUR samples) over the theoretical independence-formula approximation.
    # Lookup order:
    #   1) <PGS>.<subject>.json  — empirical restricted to variants the subject
    #      actually has. Produces apples-to-apples z-scores when coverage < 100%
    #      by eliminating the "subject is missing variants that 1000G samples
    #      have" systematic bias (observed for PGS000889 LDL: full-calibration
    #      z = -7.5 despite 92% coverage because the 8% missing variants
    #      accumulate systematically in the +LDL direction).
    #   2) <PGS>.json  — full-panel 1000G EUR empirical distribution.
    # The empirical file is the product of scripts/calibrate_prs_empirical.py.
    pop_dir = PROJECT_ROOT / 'reference' / 'population_cache'
    subject_empirical = pop_dir / 'prs_empirical' / f'{pgs_id}.{subject_id}.json'
    full_empirical = pop_dir / 'prs_empirical' / f'{pgs_id}.json'
    if subject_empirical.exists():
        empirical_path = subject_empirical
        print(f'Using subject-observed empirical: {empirical_path.name}')
    else:
        empirical_path = full_empirical
    if empirical_path.exists():
        emp = json.loads(empirical_path.read_text(encoding='utf-8'))
        theo_mean = emp['mean']
        theo_sd = emp['sd']
        # Use the score computed over only those SNPs that contributed (no AF filter here —
        # empirical calibration makes the per-SNP AF lookup irrelevant for mean/SD).
        score_for_z = result['score']
        calibration = f'empirical_1kg_eur_n{emp["n_samples"]}'
        n_with_af = n_contributing = result['n_contributing']
        print(f'Calibration: empirical 1000G EUR (N={emp["n_samples"]}); mean={theo_mean:.4f}, SD={theo_sd:.4f}')
        z, percentile = zscore_and_percentile(score_for_z, theo_mean, theo_sd)
        details_with_af = result['contributing_details']
        predicted, ci_low, ci_high, residual_sd = predicted_value(z, cfg, sex)
        # Skip the theoretical-AF path below
        empirical_used = True
    else:
        empirical_used = False

    af_parquet_path = pop_dir / '1kg_eur_afs.parquet'
    af_json_path = pop_dir / '1kg_eur_afs.json'

    if empirical_used:
        pass  # z / predicted already set above
    else:
        if af_parquet_path.exists():
            print(f'Loading AF parquet: {af_parquet_path.name} ({af_parquet_path.stat().st_size/1e6:.0f} MB)')
            need_df = pd.DataFrame({
                'chrom': [str(d.get('chr')) for d in result['contributing_details']],
                'pos': [d.get('pos') for d in result['contributing_details']],
            }).drop_duplicates()
            chroms_needed = set(need_df['chrom'])
            af_df = pd.read_parquet(af_parquet_path, filters=[('chrom', 'in', list(chroms_needed))])
            af_df['chrom'] = af_df['chrom'].astype(str)
            af_df = af_df.merge(need_df, on=['chrom', 'pos'], how='inner')
            print(f'  matched {len(af_df):,} of {len(need_df):,} contributing positions')
            theo_mean, theo_sd, details_with_af = theoretical_from_afs_parquet(
                result['contributing_details'], af_df
            )
        else:
            af_cache = json.loads(af_json_path.read_text(encoding='utf-8')) if af_json_path.exists() else {}
            theo_mean, theo_sd, details_with_af = theoretical_from_afs_json(
                result['contributing_details'], af_cache
            )

        n_with_af = len(details_with_af)
        n_contributing = result['n_contributing']

        if n_with_af >= 0.5 * n_contributing:
            score_for_z = sum(d['dosage'] * d['effect_weight'] for d in details_with_af)
            calibration = 'per_snp_eur_afs_from_1000G_phase3'
            print(f'Calibration: per-SNP EUR AFs on {n_with_af}/{n_contributing} contributing')
        else:
            theo_mean, theo_sd = theoretical_uniform(result['contributing_betas'])
            score_for_z = result['score']
            calibration = 'theoretical_EUR_AF_0.35_uniform'
            print(f'Calibration: uniform p=0.35 (only {n_with_af}/{n_contributing} have AF)')

        z, percentile = zscore_and_percentile(score_for_z, theo_mean, theo_sd)
        predicted, ci_low, ci_high, residual_sd = predicted_value(z, cfg, sex)

    # All four numbers visible together per the reporting convention:
    #   PRS z-score (SDs on PRS distribution) | percentile | phenotype-SD unit | raw value + CI | r²
    anchor = cfg.get('anchors', {}).get(sex) or cfg.get('anchors', {}).get('any') or {}
    anchor_mean = anchor.get('mean')
    anchor_sd = anchor.get('sd')
    units = cfg.get('units', '')

    # Single-trait console output is a key-value list per SKILL.md Rule 10.1.
    # Display respects the user's unit preference (default US-customary).
    pref = _get_unit_preference(subject_id)
    naive_z = (anchor_mean + anchor_sd * z) if anchor_mean is not None else None

    print()
    print(f'{cfg["label"]} ({pgs_id})')
    print(f'  PRS z-score:        {z:+.3f} SD ({percentile:.1f}th percentile on PRS distribution)')
    if predicted is not None:
        print(f'  Phenotype anchor:   {fmt_anchor(anchor_mean, anchor_sd, units, pref)}')
        print(f'  Predicted:          {fmt_value(predicted, units, pref)}  '
              f'(95% CI {fmt_value(ci_low, units, pref)} – {fmt_value(ci_high, units, pref)})')
        print(f'  Residual SD:        {fmt_value(residual_sd, units, pref)}')
    print(f'  r²:                 {cfg["r_squared"]:.2f}  (PRS explains ~{cfg["r_squared"]:.0%} of trait variance)')
    if naive_z is not None and predicted is not None:
        print(f'  Naive (if r²=1):    {fmt_value(naive_z, units, pref)}  '
              f'(regression-to-mean attenuates to {fmt_value(predicted, units, pref)})')

    # Self-report cross-check preview (the same block gets written to the ledger
    # notes further down via the guardrail). Shown here so the runtime user
    # sees the match/mismatch status, not just the raw PRS numbers.
    _preview_anchor = cfg.get('anchors', {}).get(sex) or cfg.get('anchors', {}).get('any') or {}
    _preview_block, _preview_structured = build_self_report_crosscheck(
        trait_cfg=cfg,
        profile=profile,
        predicted=predicted,
        ci_low=ci_low,
        ci_high=ci_high,
        residual_sd=residual_sd,
        anchor=_preview_anchor,
    )
    if _preview_block:
        status = (
            _preview_structured.get('self_report_crosscheck', {}).get('status', '?').upper()
        )
        print(f'  Self-report cross-check: {status}')
    print()

    # Build finding
    claim = cfg['claim_template'].format(
        predicted=f'{predicted:.2f} {cfg.get("units", "")}' if predicted is not None else 'n/a',
        ci_low=f'{ci_low:.2f}' if ci_low is not None else 'n/a',
        ci_high=f'{ci_high:.2f}' if ci_high is not None else 'n/a',
        percentile=f'{percentile:.0f}',
        z=f'{z:+.2f}',
        r2=f'{cfg["r_squared"]:.2f}',
        label=cfg['label'],
        pheno_sd=f'{anchor_sd:g}' if anchor_sd is not None else 'n/a',
        pheno_mean=f'{anchor_mean:g}' if anchor_mean is not None else 'n/a',
    )

    # Rule-10/11 guardrail: always build the self-report cross-check block
    # before finalizing the ledger row. See build_self_report_crosscheck.
    crosscheck_block, crosscheck_effect = build_self_report_crosscheck(
        trait_cfg=cfg,
        profile=profile,
        predicted=predicted,
        ci_low=ci_low,
        ci_high=ci_high,
        residual_sd=residual_sd,
        anchor=anchor,
    )

    base_notes = cfg['notes_template'].format(
        pgs_id=pgs_id,
        label=cfg['label'],
        r2=f'{cfg["r_squared"]:.2f}',
        panel=f'{result["n_panel"]:,}',
        contributing=f'{result["n_contributing"]:,}',
        coverage=f'{result["coverage_fraction"]:.1%}',
        unexplained=f'{(1 - cfg["r_squared"]):.0%}',
    )
    # Trim trailing whitespace then append the cross-check block with a single
    # leading space (build_self_report_crosscheck already starts with a space).
    full_notes = base_notes.rstrip() + crosscheck_block

    # Hard guardrail: if a self_report_key is declared, the notes MUST contain
    # the cross-check phrase. Prevents a template refactor from silently
    # dropping it.
    if cfg.get('self_report_key') and 'Self-report cross-check' not in full_notes:
        raise RuntimeError(
            f'Rule-10/11 guardrail failed: trait {cfg.get("topic")} declares '
            f'self_report_key={cfg["self_report_key"]!r} but the assembled '
            f'notes do not contain "Self-report cross-check". Refusing to '
            f'write the finding. Check build_self_report_crosscheck().'
        )

    effect = {
        'type': 'prs',
        'pgs_id': pgs_id,
        'trait': cfg['label'],
        'z_score': z,
        'percentile': percentile,
        'predicted_value': predicted,
        'ci95': [ci_low, ci_high] if ci_low is not None else None,
        'residual_sd': residual_sd,
        'r_squared_from_paper': cfg['r_squared'],
        'calibration_method': calibration,
        'sex_assumed': sex,
        'n_contributing': n_contributing,
        'n_with_af': n_with_af,
    }
    effect.update(crosscheck_effect)

    rec = {
        'subject_id': subject_id,
        'topic': cfg['topic'],
        'claim': claim,
        'variants': [{
            'pgs_id': pgs_id,
            'trait': cfg['label'],
            'panel_size': result['n_panel'],
            'on_chip': result['n_on_chip'],
            'contributing': result['n_contributing'],
            'coverage_fraction': round(result['coverage_fraction'], 4),
            'strand_flips': result['n_strand_flip'],
            'palindromic_skipped': result['n_palindromic_skipped'],
        }],
        'effect': effect,
        'cohort_ancestry': 'European',
        'subject_ancestry_match': ancestry_match,
        'source_ids': [f'pgs_catalog:{pgs_id}'],
        'notes': full_notes,
        'evidence_class': 'well_replicated_common_variant',
        'replication_count': cfg.get('replication_count', 5),
        'inference_confidence': cfg.get('inference_confidence', 'moderate'),
        'clinvar_significance': None,
        'clinvar_review_stars': None,
        'pvalue': None,
        'n_cases': None,
        'n_controls': None,
        'odds_ratio': None,
        'investigation_id': None,
    }
    tier, downgrade = compute_tier(rec)
    rec['tier_computed'] = tier
    rec['tier_computed_at'] = datetime.now(timezone.utc).isoformat()
    rec['tier_rule_version'] = TIER_RULE_VERSION
    rec['ancestry_downgrade'] = downgrade

    ledger_io.append_source(
        source_id=f'pgs_catalog:{pgs_id}',
        kind='pgs_catalog',
        url=f'https://www.pgscatalog.org/score/{pgs_id}/',
        accessed_at=datetime.now(timezone.utc).isoformat(),
        citation=cfg.get('citation', ''),
        ancestry_cohort='European',
    )

    fid = ledger_io.append_finding(**rec)
    print()
    print(f'Appended finding {fid} (tier {tier})')


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print('Usage: python scripts/run_prs.py <subject> <trait>')
        print(f'Available traits: {", ".join(TRAITS.keys())}')
        sys.exit(1)
    run(sys.argv[1], sys.argv[2])
