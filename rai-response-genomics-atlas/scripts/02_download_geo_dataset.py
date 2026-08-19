#!/usr/bin/env python3
"""02 — download a GEO dataset's SOFT family file and series matrix.

Usage:  python3 scripts/02_download_geo_dataset.py <ACCESSION>

Tries GEOparse first; falls back to direct GEO FTP URLs if needed.
Outputs land in data/raw/<ACCESSION>/.
"""
from __future__ import annotations
import gzip
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw"


def geo_stub(acc: str) -> str:
    # GEO directory naming: GSE151179 → GSE151nnn
    return acc[:-3] + "nnn"


def download_with_geoparse(acc: str, out_dir: Path) -> bool:
    try:
        import GEOparse
        out_dir.mkdir(parents=True, exist_ok=True)
        gse = GEOparse.get_GEO(geo=acc, destdir=str(out_dir), include_data=True, silent=False)
        # write a small metadata summary
        n_samples = len(gse.gsms) if hasattr(gse, "gsms") else 0
        summary = out_dir / f"{acc}_summary.txt"
        summary.write_text(
            f"acc: {acc}\nn_samples: {n_samples}\nplatforms: {list(gse.gpls.keys()) if hasattr(gse,'gpls') else []}\n"
        )
        return True
    except Exception as e:
        print(f"  GEOparse failed: {e}", file=sys.stderr)
        return False


def download_fallback(acc: str, out_dir: Path) -> bool:
    import urllib.request
    out_dir.mkdir(parents=True, exist_ok=True)
    base = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{geo_stub(acc)}/{acc}"
    candidates = [
        f"{base}/matrix/{acc}_series_matrix.txt.gz",
        f"{base}/soft/{acc}_family.soft.gz",
    ]
    ok_any = False
    for url in candidates:
        fname = url.rsplit("/", 1)[1]
        out = out_dir / fname
        try:
            print(f"  fetching {url} → {out}")
            urllib.request.urlretrieve(url, out)
            ok_any = True
        except Exception as e:
            print(f"  failed: {e}", file=sys.stderr)
    return ok_any


def main() -> None:
    if len(sys.argv) < 2:
        print("usage: 02_download_geo_dataset.py <ACCESSION>", file=sys.stderr); sys.exit(1)
    acc = sys.argv[1].strip()
    out_dir = RAW / acc
    print(f"# download {acc} → {out_dir}")
    ok = download_with_geoparse(acc, out_dir)
    if not ok:
        print("  → falling back to direct GEO FTP")
        ok = download_fallback(acc, out_dir)
    if not ok:
        print(f"ERROR: could not fetch {acc}", file=sys.stderr); sys.exit(2)
    listing = sorted(p.name for p in out_dir.iterdir())
    print(f"# files in {out_dir}:")
    for f in listing:
        size = (out_dir / f).stat().st_size
        print(f"  {size:>12d}  {f}")


if __name__ == "__main__":
    main()
