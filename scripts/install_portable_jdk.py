"""Download and extract a portable JDK ZIP for Beagle, no admin rights needed.

Places the JDK under reference/imputation/jdk/ so it's self-contained with the
imputation tooling. Exposes the java.exe path for other scripts.
"""
from __future__ import annotations

import sys
import zipfile
from pathlib import Path

import requests

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
JDK_DIR = PROJECT_ROOT / 'reference' / 'imputation' / 'jdk'

# Temurin JDK 21 Windows x64 ZIP (portable, no installer)
ADOPTIUM_URL = (
    'https://github.com/adoptium/temurin21-binaries/releases/download/'
    'jdk-21.0.6%2B7/OpenJDK21U-jdk_x64_windows_hotspot_21.0.6_7.zip'
)


def main():
    JDK_DIR.mkdir(parents=True, exist_ok=True)
    zip_path = JDK_DIR / 'temurin21.zip'
    if not zip_path.exists() or zip_path.stat().st_size < 100_000_000:
        print(f'Downloading Temurin JDK 21 from {ADOPTIUM_URL}...')
        resp = requests.get(ADOPTIUM_URL, stream=True, timeout=120)
        resp.raise_for_status()
        total = int(resp.headers.get('Content-Length', 0))
        written = 0
        with open(zip_path, 'wb') as fh:
            for chunk in resp.iter_content(chunk_size=1 << 20):
                if chunk:
                    fh.write(chunk)
                    written += len(chunk)
        print(f'  Downloaded {written:,} bytes')
    else:
        print(f'Already downloaded: {zip_path} ({zip_path.stat().st_size:,} bytes)')

    # Check if already extracted
    extracted = list(JDK_DIR.glob('jdk-*/bin/java.exe'))
    if extracted:
        java_exe = extracted[0]
        print(f'Already extracted: {java_exe}')
    else:
        print(f'Extracting to {JDK_DIR}...')
        with zipfile.ZipFile(zip_path, 'r') as zf:
            zf.extractall(JDK_DIR)
        extracted = list(JDK_DIR.glob('jdk-*/bin/java.exe'))
        if not extracted:
            sys.exit('Extraction did not produce a jdk-*/bin/java.exe')
        java_exe = extracted[0]
        print(f'Extracted to: {java_exe.parent.parent}')

    print()
    print(f'java.exe: {java_exe}')
    print(f'JAVA_HOME: {java_exe.parent.parent}')
    print()
    # Version check
    import subprocess
    result = subprocess.run([str(java_exe), '-version'], capture_output=True, text=True)
    print(result.stderr.strip())


if __name__ == '__main__':
    main()
