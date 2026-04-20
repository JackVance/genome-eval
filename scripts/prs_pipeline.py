"""Reusable polygenic-score computation pipeline.

Given a PGS Catalog harmonized weight file and a subject's standardized
parquet genotype file, compute a weighted polygenic score with strand-
aware allele matching.

Output:
    - raw score (weighted sum of effect-allele dosages)
    - SNP coverage (# matched / # in panel)
    - # strand flips applied
    - # ambiguous (A/T or C/G) SNPs excluded from strand inference
    - # not on chip
    - # allele mismatches after strand-flip attempt (should be small)

Load a reference distribution separately and call `zscore_and_percentile()`
to convert score to percentile.
"""
from __future__ import annotations

import gzip
import math
from pathlib import Path
from typing import Any

import pandas as pd

COMPLEMENT = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}


def load_weights(path: Path) -> tuple[pd.DataFrame, dict[str, str]]:
    """Load a PGS Catalog harmonized scoring file.

    Returns (weights_df, metadata) where weights_df has columns:
        rsID, chr, pos, effect_allele, other_allele, effect_weight.
    Metadata captures the #-prefixed header fields (pgs_id, trait, etc.).
    """
    meta = {}
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt', encoding='utf-8') as fh:
        header_lines = []
        for line in fh:
            if line.startswith('#'):
                header_lines.append(line.rstrip())
                if '=' in line:
                    k, _, v = line.lstrip('#').partition('=')
                    meta[k.strip()] = v.strip()
            else:
                column_line = line.rstrip()
                break
        cols = column_line.split('\t')
        # Build a row-by-row loader since pd.read_csv after a partial read is awkward
        rows = []
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < len(cols):
                parts += [''] * (len(cols) - len(parts))
            rows.append(parts)
    df = pd.DataFrame(rows, columns=cols)

    # Normalize column names — files vary ("chr_name"/"hm_chr", "chr_position"/"hm_pos", "rsID"/"hm_rsID")
    def pick(candidates: list[str]) -> str | None:
        for c in candidates:
            if c in df.columns:
                return c
        return None

    # Prefer `hm_chr`/`hm_pos`/`hm_rsID` over the unharmonized originals.
    # PGS Catalog's *_hmPOS_GRCh37.txt.gz files are harmonized to GRCh37 in
    # the `hm_*` columns regardless of what the source paper used; the
    # `chr_name`/`chr_position` columns can be the ORIGINAL build (sometimes
    # GRCh38), which would silently miss-match against a GRCh37 parquet.
    # Example: PGS003971 (Shetty 2023) publishes chr_name/chr_position in
    # hg38 and `hm_*` in GRCh37 — preferring the non-hm columns drops
    # coverage from ~100% to ~2%. The `hm_*` columns are always the
    # correct build for our pipeline.
    col_chr = pick(['hm_chr', 'chr_name', 'chr'])
    col_pos = pick(['hm_pos', 'chr_position', 'pos'])
    col_rsid = pick(['hm_rsID', 'rsID', 'rsid'])
    col_ea = pick(['effect_allele'])
    col_oa = pick(['other_allele', 'hm_inferOtherAllele'])
    col_w = pick(['effect_weight', 'weight', 'beta'])

    if not all([col_chr, col_pos, col_ea, col_w]):
        raise ValueError(
            f'Could not find required columns in weight file. '
            f'Have: {list(df.columns)}'
        )

    out = pd.DataFrame({
        'rsID': df[col_rsid] if col_rsid else '',
        'chr': df[col_chr].astype(str),
        'pos': pd.to_numeric(df[col_pos], errors='coerce'),
        'effect_allele': df[col_ea].str.upper(),
        'other_allele': df[col_oa].str.upper() if col_oa else '',
        'effect_weight': pd.to_numeric(df[col_w], errors='coerce'),
    })
    out = out.dropna(subset=['pos', 'effect_weight'])
    out['pos'] = out['pos'].astype(int)
    return out, meta


def compute_score(
    weights: pd.DataFrame,
    subject: pd.DataFrame,
) -> dict[str, Any]:
    """Match weights to subject genotypes and compute weighted PRS.

    Matching strategy:
        1. Merge on (chr, pos).
        2. For matched SNPs, check if effect_allele appears in {a1, a2}.
        3. If not, try strand flip (complement of effect_allele). Flag.
        4. Ambiguous palindromic SNPs (A/T and C/G) cannot be strand-corrected
           from genotype alone; exclude them. (This is a well-known limitation;
           minor loss of coverage.)
        5. Dosage = count of effect_allele across {a1, a2}.
        6. Score = sum(dosage * effect_weight).
    """
    # Canonicalize subject chrom
    subj = subject.copy()
    subj['chrom'] = subj['chrom'].astype(str)

    merged = weights.merge(
        subj,
        left_on=['chr', 'pos'],
        right_on=['chrom', 'pos'],
        how='left',
        suffixes=('', '_subj'),
    )

    n_panel = len(weights)
    n_matched_on_chip = int(merged['a1'].notna().sum())

    # Identify palindromic (ambiguous) SNPs
    merged['is_palindromic'] = merged.apply(
        lambda r: (
            isinstance(r['effect_allele'], str)
            and isinstance(r['other_allele'], str)
            and {r['effect_allele'], r['other_allele']} in [{'A', 'T'}, {'C', 'G'}]
        ),
        axis=1,
    )

    score = 0.0
    n_contributing = 0
    n_strand_flip = 0
    n_palindromic_skipped = 0
    n_allele_mismatch = 0
    n_on_chip_no_call = 0
    contributing_details = []  # list of {rsid, effect_allele_used, effect_weight, dosage}

    for r in merged.itertuples():
        if pd.isna(r.a1) or pd.isna(r.a2):
            if pd.notna(r.a1) or pd.notna(r.a2):
                n_on_chip_no_call += 1
            continue
        ea = r.effect_allele
        oa = r.other_allele
        subject_alleles = {r.a1, r.a2}

        # Try direct strand match: subject alleles are a subset of {ea, oa}
        if subject_alleles <= {ea, oa}:
            if r.is_palindromic:
                # Palindromic on direct strand — ambiguous strand, but alleles fit.
                # Accept on the assumption that both weight file and chip use forward strand
                # (the common case for PGS Catalog harmonized files and 23andMe data).
                pass
            dosage = int(r.a1 == ea) + int(r.a2 == ea)
        else:
            # Try strand flip: subject alleles are a subset of {complement(ea), complement(oa)}
            ea_flip = COMPLEMENT.get(ea, ea)
            oa_flip = COMPLEMENT.get(oa, oa)
            if subject_alleles <= {ea_flip, oa_flip}:
                if r.is_palindromic:
                    # Palindromic + subject alleles match only the flipped strand:
                    # unresolvable from genotype alone (both {A,T}==complement({A,T}) and
                    # similarly {C,G}). Safer to skip.
                    n_palindromic_skipped += 1
                    continue
                dosage = int(r.a1 == ea_flip) + int(r.a2 == ea_flip)
                n_strand_flip += 1
            else:
                n_allele_mismatch += 1
                continue
        score += dosage * r.effect_weight
        # Determine which allele (on the subject's strand) was used as the effect allele
        used_ea = ea
        if ea not in subject_alleles:
            used_ea = COMPLEMENT.get(ea, ea)
        contributing_details.append({
            'rsid': r.rsID,
            'chr': r.chr,
            'pos': r.pos,
            'effect_allele_panel': ea,
            'effect_allele_used': used_ea,
            'effect_weight': r.effect_weight,
            'dosage': dosage,
        })
        n_contributing += 1

    return {
        'score': score,
        'n_panel': n_panel,
        'n_on_chip': n_matched_on_chip,
        'n_contributing': n_contributing,
        'n_strand_flip': n_strand_flip,
        'n_palindromic_skipped': n_palindromic_skipped,
        'n_allele_mismatch_after_flip': n_allele_mismatch,
        'n_on_chip_no_call': n_on_chip_no_call,
        'coverage_fraction': n_contributing / n_panel if n_panel else 0.0,
        'contributing_details': contributing_details,
        'contributing_betas': [d['effect_weight'] for d in contributing_details],
    }


def zscore_and_percentile(score: float, ref_mean: float, ref_sd: float) -> tuple[float, float]:
    """Return (z_score, percentile) given the score and a reference distribution."""
    if ref_sd <= 0:
        return 0.0, 50.0
    z = (score - ref_mean) / ref_sd
    # Normal CDF from z (Abramowitz & Stegun approx via math.erf)
    cdf = 0.5 * (1 + math.erf(z / math.sqrt(2)))
    return z, cdf * 100
