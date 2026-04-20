"""Parse imputed VCF files (Beagle output) back into a parquet with dosages.

Beagle output columns of interest:
  - GT: genotype (most-likely call)
  - DS: expected dosage of ALT (float in [0, 2])
  - AP1, AP2: allele probabilities per haplotype (float)
  - GP: genotype probabilities (three floats)

Produces `standardized-genomes/imputed/<subject>.imputed.parquet` with schema:
  rsid, chrom, pos, a1, a2, ref, alt, dosage, gp_hom_ref, gp_het, gp_hom_alt, imputed (bool)

The imputed flag is True for variants not in the original chip; False for directly-measured.

Usage:
    python scripts/vcf_to_parquet.py alice --chr 22
    python scripts/vcf_to_parquet.py alice --all
"""
from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path

import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent


def parse_format_field(format_str: str, sample_str: str) -> dict:
    keys = format_str.split(':')
    vals = sample_str.split(':')
    return dict(zip(keys, vals))


def parse_imputed_vcf(vcf_path: Path, chrom: str) -> list[dict]:
    rows = []
    opener = gzip.open if vcf_path.suffix == '.gz' else open
    with opener(vcf_path, 'rt', encoding='utf-8') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 10:
                continue
            chrom_v, pos, rsid, ref, alt, qual, filt, info, fmt, sample = parts[:10]
            # Skip multiallelic (shouldn't happen for Beagle EUR panel, but be safe)
            if ',' in alt:
                continue
            fd = parse_format_field(fmt, sample)
            gt_str = fd.get('GT', './.')
            # Parse GT — Beagle output is typically phased "0|1"
            sep = '|' if '|' in gt_str else '/'
            gt_parts = gt_str.split(sep)
            if len(gt_parts) != 2 or '.' in gt_parts:
                a1 = a2 = None
            else:
                a1 = ref if gt_parts[0] == '0' else alt
                a2 = ref if gt_parts[1] == '0' else alt
            # Dosage
            try:
                dosage = float(fd['DS']) if 'DS' in fd else None
            except Exception:
                dosage = None
            # Genotype probabilities
            gp_ref = gp_het = gp_alt = None
            if 'GP' in fd:
                try:
                    gp_vals = [float(x) for x in fd['GP'].split(',')]
                    if len(gp_vals) == 3:
                        gp_ref, gp_het, gp_alt = gp_vals
                except Exception:
                    pass

            rows.append({
                'rsid': rsid if rsid != '.' else None,
                'chrom': chrom,
                'pos': int(pos),
                'ref': ref,
                'alt': alt,
                'a1': a1,
                'a2': a2,
                'dosage': dosage,
                'gp_hom_ref': gp_ref,
                'gp_het': gp_het,
                'gp_hom_alt': gp_alt,
            })
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('subject_id')
    ap.add_argument('--chr', default='all', help='Chromosome or comma-list or "all"')
    ap.add_argument('--outdir', default=None)
    args = ap.parse_args()

    # Load original chip parquet to flag "imputed" vs "directly measured"
    original = pd.read_parquet(PROJECT_ROOT / 'standardized-genomes' / f'{args.subject_id}.parquet')
    original['chrom'] = original['chrom'].astype(str)
    original_keys = set(zip(original['chrom'].tolist(), original['pos'].tolist()))

    beagle_dir = PROJECT_ROOT / 'standardized-genomes' / 'imputed' / args.subject_id / 'beagle'
    if not beagle_dir.exists():
        sys.exit(f'No imputed directory at {beagle_dir}. Run scripts/run_imputation.py first.')

    if args.chr == 'all':
        chroms = [str(i) for i in range(1, 23)]
    else:
        chroms = args.chr.split(',')

    all_rows = []
    for chrom in chroms:
        vcf = beagle_dir / f'chr{chrom}.vcf.gz'
        if not vcf.exists():
            print(f'chr{chrom}: no imputed VCF at {vcf}, skipping')
            continue
        print(f'chr{chrom}: parsing {vcf}...')
        rows = parse_imputed_vcf(vcf, chrom)
        print(f'  {len(rows):,} variants')
        all_rows.extend(rows)

    if not all_rows:
        sys.exit('No imputed data parsed.')

    df = pd.DataFrame(all_rows)
    df['imputed'] = ~df.apply(lambda r: (r['chrom'], r['pos']) in original_keys, axis=1)

    outdir = Path(args.outdir) if args.outdir else PROJECT_ROOT / 'standardized-genomes' / 'imputed'
    outdir.mkdir(parents=True, exist_ok=True)
    out_path = outdir / f'{args.subject_id}.imputed.parquet'
    df.to_parquet(out_path, index=False)
    print()
    print(f'Wrote {len(df):,} rows to {out_path}')
    print(f'  directly measured: {(~df["imputed"]).sum():,}')
    print(f'  imputed:           {df["imputed"].sum():,}')
    chroms_processed = df['chrom'].unique().tolist()
    chip_on_same_chroms = original[original['chrom'].isin(chroms_processed)]
    if len(chip_on_same_chroms) > 0:
        print(f'  expansion on processed chromosomes ({",".join(sorted(chroms_processed, key=int))}): '
              f'{len(df) / len(chip_on_same_chroms):.1f}x vs chip directly measured')


if __name__ == '__main__':
    main()
