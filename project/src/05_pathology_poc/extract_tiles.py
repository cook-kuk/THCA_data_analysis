#!/usr/bin/env python3
"""Spot-centered H&E tile extraction for pathology AI POC.

Output: data/processed/GSE250521/tiles/{sample_id}/{spot_id}_size{tile_size}.png
        + tile_metadata.tsv (spot_id, sample_id, stage, x, y, RAI_8, DM1_like, TDS_like)

Foundation model embedding (UNI/CONCH/Virchow2/Prov-GigaPath) is OPTIONAL.
This script extracts tiles only; embedding + Ridge regression skeleton lives in
this directory's README.md.
"""
from __future__ import annotations
import argparse, logging, sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from PIL import Image

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger("tiles")
Image.MAX_IMAGE_PIXELS = None


def extract_for_sample(adata: ad.AnnData, tile_size: int, tiles_dir: Path,
                       max_per_sample: int | None) -> pd.DataFrame:
    sid = adata.obs["sample_id"].iloc[0]
    spatial = adata.uns.get("spatial", {}).get(sid, {})
    sf = spatial.get("scalefactors", {})
    hires_sf = sf.get("tissue_hires_scalef", 1.0)
    img_path = spatial.get("images", {}).get("hires_path")
    if not img_path or not Path(img_path).exists():
        log.warning("[%s] no image; skipping", sid); return pd.DataFrame()
    img = np.asarray(Image.open(img_path))
    H, W = img.shape[:2]

    in_tissue = adata.obs["in_tissue"].fillna(0).astype(int) == 1
    sub = adata[in_tissue].copy()
    if "DM1_like_score" not in sub.obs.columns:
        log.warning("[%s] no DM1_like_score; skipping", sid); return pd.DataFrame()

    if max_per_sample and sub.n_obs > max_per_sample:
        rng = np.random.default_rng(0)
        idx = rng.choice(sub.n_obs, size=max_per_sample, replace=False)
        sub = sub[idx].copy()

    out_dir = tiles_dir / sid
    out_dir.mkdir(parents=True, exist_ok=True)
    half = tile_size // 2
    rows = []
    skipped = 0
    for i in range(sub.n_obs):
        spot_id = sub.obs.index[i]
        # in hires image space
        x = int(sub.obsm["spatial"][i, 0] * hires_sf)
        y = int(sub.obsm["spatial"][i, 1] * hires_sf)
        x0, x1 = x - half, x + half
        y0, y1 = y - half, y + half
        if x0 < 0 or y0 < 0 or x1 > W or y1 > H:
            skipped += 1; continue
        tile = img[y0:y1, x0:x1]
        if tile.shape[0] != tile_size or tile.shape[1] != tile_size:
            skipped += 1; continue
        out_png = out_dir / f"{spot_id}_size{tile_size}.png"
        if not out_png.exists():
            Image.fromarray(tile).save(out_png)
        rows.append({"spot_id": spot_id, "sample_id": sid,
                     "stage": sub.obs["stage"].iloc[i],
                     "tile_path": str(out_png),
                     "tile_size": tile_size, "x_hires": x, "y_hires": y,
                     "RAI_8_score": sub.obs["RAI_8_score"].iloc[i],
                     "DM1_like_score": sub.obs["DM1_like_score"].iloc[i],
                     "TDS_like_score": sub.obs["TDS_like_score"].iloc[i]})
    log.info("[%s] %d tiles (skipped border=%d)", sid, len(rows), skipped)
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--processed-dir", default="project/data/processed/GSE250521")
    ap.add_argument("--tiles-dir", default="project/data/processed/GSE250521/tiles")
    ap.add_argument("--meta-out", default="project/results/03_pathology_poc/tile_metadata.tsv.gz")
    ap.add_argument("--tile-size", type=int, default=224, choices=[224, 448, 672])
    ap.add_argument("--max-per-sample", type=int, default=400,
                    help="cap per sample for POC (set 0 for all)")
    args = ap.parse_args()

    cap = args.max_per_sample if args.max_per_sample > 0 else None
    tiles_dir = Path(args.tiles_dir)
    all_rows = []
    for h in sorted(Path(args.processed_dir).rglob("*.scored.h5ad")):
        adata = ad.read_h5ad(h)
        df = extract_for_sample(adata, args.tile_size, tiles_dir, cap)
        if not df.empty:
            all_rows.append(df)
    if not all_rows:
        log.error("no tiles extracted"); sys.exit(1)
    meta = pd.concat(all_rows, ignore_index=True)
    Path(args.meta_out).parent.mkdir(parents=True, exist_ok=True)
    meta.to_csv(args.meta_out, sep="\t", index=False, compression="gzip")
    log.info("metadata: %d tiles → %s", len(meta), args.meta_out)
    print(meta["stage"].value_counts())
    print("samples:", meta["sample_id"].nunique())


if __name__ == "__main__":
    sys.exit(main())
