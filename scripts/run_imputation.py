"""Orchestrate Beagle 5.4 imputation for a subject, one chromosome at a time.

Prerequisites (see imputation_download.py):
  - Beagle JAR in reference/imputation/beagle/
  - 1000G Phase 3 bref3 panels in reference/imputation/1kg_ref_b37/
  - Genetic maps in reference/imputation/genetic_maps/
  - Java on PATH

Workflow per chromosome:
  1. Convert subject parquet -> chr-specific VCF (via parquet_to_vcf.py)
  2. Run Beagle: java -jar beagle.jar gt=input ref=... map=... out=...
  3. Beagle emits <out>.vcf.gz with imputed genotypes + dosages (DS) + allele probabilities (AP)

Output: standardized-genomes/imputed/<subject>/beagle/chr<N>.vcf.gz

Usage:
    python scripts/run_imputation.py alice --chr 22              # one chromosome
    python scripts/run_imputation.py alice --all --mem 8g        # all autosomes
"""
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
IMPUTATION_DIR = PROJECT_ROOT / 'reference' / 'imputation'
BEAGLE_DIR = IMPUTATION_DIR / 'beagle'
REF_DIR = IMPUTATION_DIR / '1kg_ref_b37'
MAPS_DIR = IMPUTATION_DIR / 'genetic_maps'


def find_beagle_jar() -> Path:
    jars = sorted(BEAGLE_DIR.glob('beagle*.jar'))
    if not jars:
        sys.exit(f'No Beagle JAR found in {BEAGLE_DIR}. Run: python scripts/imputation_download.py --beagle')
    return jars[-1]


def find_java() -> str:
    """Return the java executable path. Falls back to 'java' on PATH."""
    # Prefer the project-local portable JDK if present
    local = list((IMPUTATION_DIR / 'jdk').glob('jdk-*/bin/java.exe'))
    if local:
        return str(local[0])
    java_env = os.environ.get('JAVA_HOME')
    if java_env:
        jpath = Path(java_env) / 'bin' / 'java.exe'
        if jpath.exists():
            return str(jpath)
    which = shutil.which('java')
    if which:
        return which
    # Common Windows install paths
    candidates = list(Path('C:/Program Files').glob('*/jdk-*/bin/java.exe'))
    candidates += list(Path('C:/Program Files/Eclipse Adoptium').glob('*/bin/java.exe'))
    candidates += list(Path('C:/Program Files/Microsoft').glob('jdk-*/bin/java.exe'))
    candidates += list(Path(os.path.expanduser('~/AppData/Local/Programs')).glob('**/bin/java.exe'))
    if candidates:
        return str(candidates[0])
    sys.exit('No java found. Run: python scripts/install_portable_jdk.py')


def find_ref_and_map(chrom: str) -> tuple[Path, Path, Path]:
    """Return (bref3_panel, vcf_panel, genetic_map) for a chromosome."""
    ref_candidates = list(REF_DIR.glob(f'chr{chrom}.*.bref3'))
    if not ref_candidates:
        sys.exit(f'No reference panel for chr{chrom} in {REF_DIR}. Run: python scripts/imputation_download.py --chr {chrom}')
    # Prefer filtered reference VCF (biallelic-SNP-only) over the raw file,
    # since conform-gt chokes on multiallelic indels with non-unique IDs.
    filtered_vcf = list(REF_DIR.glob(f'chr{chrom}.*.filtered.vcf.gz'))
    if filtered_vcf:
        vcf_candidates = filtered_vcf
    else:
        vcf_candidates = [p for p in REF_DIR.glob(f'chr{chrom}.*.vcf.gz') if '.filtered.' not in p.name]
    if not vcf_candidates:
        sys.exit(f'No VCF reference for chr{chrom} in {REF_DIR}. Needed for conform-gt allele alignment. Run: python scripts/imputation_download.py --chr {chrom} and then python scripts/filter_ref_vcf.py --chr {chrom}')
    map_candidates = list(MAPS_DIR.glob(f'plink.chr{chrom}.*.map'))
    if not map_candidates:
        sys.exit(f'No genetic map for chr{chrom} in {MAPS_DIR}. Run: python scripts/imputation_download.py --maps')
    return ref_candidates[0], vcf_candidates[0], map_candidates[0]


def find_conform_gt_jar() -> Path:
    jars = list(BEAGLE_DIR.glob('conform-gt*.jar'))
    if not jars:
        sys.exit(f'conform-gt JAR not in {BEAGLE_DIR}. Run: python scripts/imputation_download.py --beagle')
    return jars[0]


def run_conform_gt(java: str, conform_jar: Path, input_vcf: Path, ref_vcf: Path, chrom: str, out_prefix: Path) -> Path:
    """Run conform-gt to align input VCF's REF/ALT to the reference panel.

    Writes <out_prefix>.vcf.gz (conform-gt appends .vcf.gz automatically). Returns that path.
    """
    cmd = [
        java, '-Xmx4g', '-jar', str(conform_jar),
        f'gt={input_vcf}',
        f'ref={ref_vcf}',
        f'chrom={chrom}',
        f'out={out_prefix}',
    ]
    print('  Running conform-gt (allele alignment):')
    print('   ', ' '.join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
    if result.returncode != 0:
        print('  conform-gt stdout:', result.stdout[-800:])
        print('  conform-gt stderr:', result.stderr[-800:])
        raise RuntimeError(f'conform-gt failed with code {result.returncode}')
    # Log tail
    print('  conform-gt log (tail):')
    for line in (result.stdout or '').strip().splitlines()[-10:]:
        print(f'    {line}')
    expected = Path(str(out_prefix) + '.vcf.gz')
    if not expected.exists():
        raise RuntimeError(f'conform-gt did not produce {expected}')
    return expected


def run_beagle(
    java: str,
    jar: Path,
    input_vcf: Path,
    ref_panel: Path,
    genetic_map: Path,
    out_prefix: Path,
    mem: str = '8g',
    threads: int = 4,
) -> subprocess.CompletedProcess:
    cmd = [
        java, f'-Xmx{mem}', '-jar', str(jar),
        f'gt={input_vcf}',
        f'ref={ref_panel}',
        f'map={genetic_map}',
        f'out={out_prefix}',
        f'nthreads={threads}',
        'impute=true',
        'ap=true',       # output allele probabilities
        'gp=true',       # output genotype probabilities
    ]
    print('  Running Beagle:')
    print('   ', ' '.join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=7200)
    return result


def impute_chromosome(subject_id: str, chrom: str, mem: str, threads: int) -> bool:
    java = find_java()
    jar = find_beagle_jar()
    conform_jar = find_conform_gt_jar()
    ref_panel_bref3, ref_panel_vcf, genetic_map = find_ref_and_map(chrom)

    # Ensure raw input VCF exists (run parquet_to_vcf.py if needed)
    input_vcf_dir = PROJECT_ROOT / 'standardized-genomes' / 'imputed' / subject_id / 'vcf'
    input_vcf = input_vcf_dir / f'{subject_id}.chr{chrom}.vcf.gz'
    if not input_vcf.exists():
        print(f'  Creating input VCF for chr{chrom}...')
        subprocess.run(
            [sys.executable, str(SCRIPT_DIR / 'parquet_to_vcf.py'), subject_id, '--chr', chrom],
            check=True,
        )

    # Align alleles to reference panel with conform-gt
    conform_out_prefix = input_vcf_dir / f'{subject_id}.chr{chrom}.aligned'
    aligned_vcf = input_vcf_dir / f'{subject_id}.chr{chrom}.aligned.vcf.gz'
    if not aligned_vcf.exists():
        run_conform_gt(java, conform_jar, input_vcf, ref_panel_vcf, chrom, conform_out_prefix)
        if not aligned_vcf.exists():
            raise RuntimeError(f'conform-gt did not produce {aligned_vcf}')
    else:
        print(f'  Already aligned: {aligned_vcf}')

    out_dir = PROJECT_ROOT / 'standardized-genomes' / 'imputed' / subject_id / 'beagle'
    out_dir.mkdir(parents=True, exist_ok=True)
    out_prefix = out_dir / f'chr{chrom}'

    expected_out = out_prefix.with_suffix('.vcf.gz')
    if expected_out.exists() and expected_out.stat().st_size > 0:
        print(f'  chr{chrom}: imputed VCF already exists ({expected_out.stat().st_size:,} bytes). Skipping.')
        return True

    beagle_out_prefix = out_dir / f'chr{chrom}'
    result = run_beagle(java, jar, aligned_vcf, ref_panel_bref3, genetic_map, beagle_out_prefix, mem=mem, threads=threads)
    if result.returncode != 0:
        print(f'  BEAGLE FAILED for chr{chrom} (exit {result.returncode})')
        print('  stdout:', result.stdout[-500:])
        print('  stderr:', result.stderr[-500:])
        return False
    # Log the tail of Beagle's output for context
    print('  Beagle log (tail):')
    for line in result.stdout.strip().splitlines()[-15:]:
        print(f'    {line}')
    print(f'  Done: {expected_out}')
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('subject_id')
    ap.add_argument('--chr', default='22', help='Chromosome (single) or comma-separated list or "all"')
    ap.add_argument('--all', action='store_true', help='Impute all 22 autosomes')
    ap.add_argument('--mem', default='8g', help='Java heap size (e.g., 8g, 16g)')
    ap.add_argument('--threads', type=int, default=4)
    args = ap.parse_args()

    if args.all:
        chroms = [str(i) for i in range(1, 23)]
    elif args.chr == 'all':
        chroms = [str(i) for i in range(1, 23)]
    else:
        chroms = args.chr.split(',')

    for chrom in chroms:
        print(f'\n=== chr{chrom} ===')
        ok = impute_chromosome(args.subject_id, chrom, args.mem, args.threads)
        if not ok:
            print(f'  Stopping after chr{chrom} failure.')
            break


if __name__ == '__main__':
    main()
