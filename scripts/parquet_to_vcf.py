"""Convert a subject's standardized parquet to VCF format for Beagle imputation.

Per chromosome output (Beagle imputes one chromosome at a time). VCF spec:
https://samtools.github.io/hts-specs/VCFv4.2.pdf

Strategy:
  - Looks up canonical REF / ALT for each (chrom, pos) from the 1000 Genomes
    Phase 3 VCF reference panel. This gives proper biallelic output even at
    sites where the subject is homozygous — we know which of the two possible
    alleles they carry because the reference tells us what both options are.
  - Writes genotype 0/0 (homozygous REF), 0/1 (heterozygous), or 1/1
    (homozygous ALT). If subject has a strand-flipped allele pair, we complement
    and re-check. If after flipping the alleles still don't match reference
    REF/ALT, the site is skipped (truly off-reference, rare).
  - Indels (23andMe D/I encoding) and multi-allelic sites are skipped.
  - No-calls encoded as `./.` dropped entirely (Beagle doesn't use them).

Usage:
    python scripts/parquet_to_vcf.py alice --chr 22
    python scripts/parquet_to_vcf.py alice --all
"""
from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path

import pandas as pd

COMPLEMENT = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent


def load_ref_panel_alleles(ref_vcf: Path, chrom: str) -> dict[int, tuple[str, str]]:
    """Build a {pos: (REF, ALT)} mapping from the reference panel VCF for this chromosome.

    Skips multiallelic sites (those are rare and not useful for imputation of chip data).
    """
    mapping: dict[int, tuple[str, str]] = {}
    with gzip.open(ref_vcf, 'rt', encoding='utf-8') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            # We only need the first 5 columns: CHROM POS ID REF ALT
            parts = line.split('\t', 5)
            if len(parts) < 5:
                continue
            chrom_v, pos_s, _rsid, ref, alt = parts[0], parts[1], parts[2], parts[3], parts[4]
            if chrom_v != chrom:
                continue
            if ',' in alt:  # multiallelic
                continue
            if len(ref) != 1 or len(alt) != 1:
                continue
            if ref not in 'ACGT' or alt not in 'ACGT':
                continue
            try:
                pos = int(pos_s)
            except ValueError:
                continue
            # Keep first occurrence at each position (there shouldn't be duplicates in 1kg panel)
            if pos not in mapping:
                mapping[pos] = (ref, alt)
    return mapping


def write_vcf_for_chromosome(
    df: pd.DataFrame,
    subject_id: str,
    chrom: str,
    out_path: Path,
    ref_alleles: dict[int, tuple[str, str]] | None = None,
):
    """Write VCF for a single chromosome.

    df must already be filtered to that chromosome, sorted by pos.
    If ref_alleles is provided, uses reference REF/ALT (proper biallelic output).
    If not, falls back to alphabetical ordering (may produce REF=ALT='.' at homozygous sites).
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)
    open_fn = gzip.open if out_path.suffix == '.gz' else open
    with open_fn(out_path, 'wt', encoding='utf-8', newline='\n') as fh:
        # Header
        fh.write('##fileformat=VCFv4.2\n')
        fh.write('##source=genome-eval parquet_to_vcf.py\n')
        fh.write('##reference=GRCh37\n')
        fh.write(f'##contig=<ID={chrom}>\n')
        fh.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        fh.write(f'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{subject_id}\n')

        n_written = 0
        n_skipped_indel = 0
        n_skipped_nocall = 0
        n_skipped_noncanonical = 0
        n_not_in_ref = 0
        n_strand_flipped = 0
        n_allele_mismatch = 0

        for r in df.itertuples():
            a1, a2 = r.a1, r.a2
            if pd.isna(a1) or pd.isna(a2) or a1 == '-' or a2 == '-':
                n_skipped_nocall += 1
                continue
            if a1 in ('D', 'I') or a2 in ('D', 'I'):
                n_skipped_indel += 1
                continue
            if a1 not in 'ACGT' or a2 not in 'ACGT':
                n_skipped_noncanonical += 1
                continue

            pos = int(r.pos)
            rsid = r.rsid if r.rsid and str(r.rsid) != 'nan' else '.'

            if ref_alleles is not None:
                ref_alt = ref_alleles.get(pos)
                if ref_alt is None:
                    n_not_in_ref += 1
                    continue
                ref, alt = ref_alt
                # Check whether subject alleles match {ref, alt} directly or via strand flip
                subj = {a1, a2}
                if subj <= {ref, alt}:
                    pass  # direct match
                elif subj <= {COMPLEMENT[ref], COMPLEMENT[alt]}:
                    # Strand flip: complement the subject's alleles
                    a1 = COMPLEMENT[a1]
                    a2 = COMPLEMENT[a2]
                    n_strand_flipped += 1
                else:
                    n_allele_mismatch += 1
                    continue
                # Encode GT
                if a1 == ref and a2 == ref:
                    gt = '0/0'
                elif a1 == alt and a2 == alt:
                    gt = '1/1'
                else:
                    gt = '0/1'
            else:
                # Fallback (no reference): alphabetical; homozygous -> ALT='.'
                if a1 == a2:
                    ref, alt = a1, '.'
                    gt = '0/0'
                else:
                    if a1 < a2:
                        ref, alt = a1, a2
                    else:
                        ref, alt = a2, a1
                    gt = '0/1'

            fh.write(f'{chrom}\t{pos}\t{rsid}\t{ref}\t{alt}\t.\tPASS\t.\tGT\t{gt}\n')
            n_written += 1

    return {
        'written': n_written,
        'skipped_indel': n_skipped_indel,
        'skipped_nocall': n_skipped_nocall,
        'skipped_noncanonical': n_skipped_noncanonical,
        'not_in_ref': n_not_in_ref,
        'strand_flipped': n_strand_flipped,
        'allele_mismatch': n_allele_mismatch,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('subject_id')
    ap.add_argument('--chr', help='Single chromosome (e.g., 22), or "all"', default='all')
    ap.add_argument('--outdir', default=None,
                    help='Output directory (default: standardized-genomes/imputed/<subject>/vcf/)')
    ap.add_argument('--no-gzip', action='store_true')
    args = ap.parse_args()

    parquet_path = PROJECT_ROOT / 'standardized-genomes' / f'{args.subject_id}.parquet'
    df = pd.read_parquet(parquet_path)
    df['chrom'] = df['chrom'].astype(str)

    outdir = Path(args.outdir) if args.outdir else (
        PROJECT_ROOT / 'standardized-genomes' / 'imputed' / args.subject_id / 'vcf'
    )

    if args.chr == 'all':
        chroms = [str(i) for i in range(1, 23)]  # autosomes only; Beagle handles X separately
    else:
        chroms = [args.chr]

    ref_dir = PROJECT_ROOT / 'reference' / 'imputation' / '1kg_ref_b37'
    for chrom in chroms:
        sub = df[df['chrom'] == chrom].sort_values('pos')
        if sub.empty:
            print(f'chr{chrom}: no SNPs in subject data, skipping')
            continue
        # Deduplicate by position: some chips have multiple probes at the same site
        # (e.g. 23andMe i-probe + rs-probe). Keep the first — Beagle only wants one.
        n_before = len(sub)
        sub = sub.drop_duplicates(subset='pos', keep='first')
        if len(sub) < n_before:
            print(f'chr{chrom}: dropped {n_before - len(sub)} duplicate-position probes.')

        ref_vcf_candidates = list(ref_dir.glob(f'chr{chrom}.*.vcf.gz'))
        ref_alleles = None
        if ref_vcf_candidates:
            print(f'chr{chrom}: loading reference alleles from {ref_vcf_candidates[0].name}...')
            ref_alleles = load_ref_panel_alleles(ref_vcf_candidates[0], chrom)
            print(f'  {len(ref_alleles):,} biallelic reference positions loaded.')
        else:
            print(f'chr{chrom}: no reference VCF found; using fallback alphabetical allele ordering (may produce REF=ALT=. at homozygous sites).')

        suffix = '.vcf' if args.no_gzip else '.vcf.gz'
        out_path = outdir / f'{args.subject_id}.chr{chrom}{suffix}'
        stats = write_vcf_for_chromosome(sub, args.subject_id, chrom, out_path, ref_alleles)
        print(
            f'chr{chrom:>2s}: {stats["written"]:>6d} SNPs written | '
            f'skipped: {stats["skipped_nocall"]} no-call, '
            f'{stats["skipped_indel"]} indel, {stats["skipped_noncanonical"]} non-canonical, '
            f'{stats["not_in_ref"]} not-in-ref, {stats["allele_mismatch"]} allele-mismatch | '
            f'strand-flipped: {stats["strand_flipped"]} | -> {out_path}'
        )


if __name__ == '__main__':
    main()
