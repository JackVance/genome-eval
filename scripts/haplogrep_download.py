"""Download + unpack HaploGrep3 (mtDNA haplogroup classifier).

One-time install. HaploGrep3 ships as a zip of JARs + launcher scripts.
We extract to reference/haplogroups/mtdna/haplogrep3/ so future subjects reuse.
"""
from __future__ import annotations

import io
import sys
import urllib.request
import zipfile
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
DEST = PROJECT_ROOT / "reference" / "haplogroups" / "mtdna" / "haplogrep3"

RELEASE_URL = (
    "https://github.com/genepi/haplogrep3/releases/download/v3.2.2/"
    "haplogrep3-3.2.2-windows.zip"
)


def main() -> int:
    DEST.mkdir(parents=True, exist_ok=True)
    sentinel = DEST / "haplogrep3.jar"
    if sentinel.exists() and sentinel.stat().st_size > 100_000:
        print(f"HaploGrep3 already installed at {DEST}")
        return 0

    print(f"Downloading {RELEASE_URL}")
    with urllib.request.urlopen(RELEASE_URL) as resp:
        data = resp.read()
    print(f"  downloaded {len(data)/1e6:.1f} MB")

    with zipfile.ZipFile(io.BytesIO(data)) as zf:
        zf.extractall(DEST)
    print(f"Extracted {len(data)/1e6:.1f} MB to {DEST}")

    # Normalize: some zip builds put everything under a top-level dir.
    top_level_dirs = [p for p in DEST.iterdir() if p.is_dir()]
    if len(top_level_dirs) == 1 and not sentinel.exists():
        inner = top_level_dirs[0]
        for child in inner.iterdir():
            child.rename(DEST / child.name)
        inner.rmdir()
        print(f"  flattened {inner.name}")

    for p in DEST.iterdir():
        print(f"  {p.name} ({p.stat().st_size/1e6:.2f} MB)" if p.is_file() else f"  {p.name}/")
    if not sentinel.exists():
        print(f"WARNING: expected haplogrep3.jar not found after extraction.")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
