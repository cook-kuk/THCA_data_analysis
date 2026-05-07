#!/usr/bin/env python3
"""Extract H&E tiles from 12 external Visium slides (GSE230424 + GSE248205).
Reuses existing GSE250521 pattern to get full 28-slide tile dataset for AI POC."""
from pathlib import Path
import argparse
import anndata as ad
import numpy as np
import pandas as pd
from PIL import Image

Image.MAX_IMAGE_PIXELS = None


def extract_one(adata, tile_size, out_dir, max_per_sample, score_df):
    sid = adata.obs["sample_id"].iloc[0]
    spatial = adata.uns.get("spatial", {}).get(sid, {})
    sf = spatial.get("scalefactors", {})
    hires_sf = sf.get("tissue_hires_scalef", 1.0) if sf else 1.0
    # Try image
    img_path = (spatial.get("images", {}) or {}).get("hires_path") \
               or (spatial.get("images", {}) or {}).get("HE_path")
    if not img_path or not Path(img_path).exists():
        print(f"[skip] {sid}: no image"); return pd.DataFrame()
    img = np.asarray(Image.open(img_path))
    if img.ndim == 2: img = np.stack([img]*3, axis=-1)
    H, W = img.shape[:2]
    # use full-res coords scaled by hires_sf only if hires image; for HE.jpg without scalefactors, use raw pxl
    # GSE230424 has HE.jpg full-res; coords already in fullres → sf=1.0
    # GSE248205 has hires.png → coords need *= hires_sf
    half = tile_size // 2
    sub = adata[adata.obs["in_tissue"].fillna(0).astype(int) == 1].copy()
    if max_per_sample and sub.n_obs > max_per_sample:
        rng = np.random.default_rng(0)
        idx = rng.choice(sub.n_obs, size=max_per_sample, replace=False)
        sub = sub[idx].copy()

    # Score lookup
    sc_lookup = score_df.set_index("spot_id") if score_df is not None and "spot_id" in score_df.columns else None

    rows = []
    skipped = 0
    out_dir.mkdir(parents=True, exist_ok=True)
    for i in range(sub.n_obs):
        spot_id = sub.obs.index[i]
        x = int(sub.obsm["spatial"][i, 0] * hires_sf)
        y = int(sub.obsm["spatial"][i, 1] * hires_sf)
        x0, x1 = x - half, x + half; y0, y1 = y - half, y + half
        if x0 < 0 or y0 < 0 or x1 > W or y1 > H:
            skipped += 1; continue
        tile = img[y0:y1, x0:x1]
        if tile.shape[0] != tile_size or tile.shape[1] != tile_size:
            skipped += 1; continue
        out_png = out_dir / f"{spot_id}_size{tile_size}.png"
        if not out_png.exists():
            Image.fromarray(tile).save(out_png)
        row = {"spot_id": spot_id, "sample_id": sid, "tile_path": str(out_png),
               "tile_size": tile_size,
               "stage": sub.obs.get("stage", sub.obs.get("condition_inferred")).iloc[i],
               "dataset": sub.obs["dataset"].iloc[0] if "dataset" in sub.obs.columns else "",
               "x_hires": x, "y_hires": y}
        if sc_lookup is not None and spot_id in sc_lookup.index:
            for col in ("DM1_like_score_resid","RAI_8_score_resid","THYROID_NONOVERLAP_score_resid","Epithelial_score_resid"):
                if col in sc_lookup.columns:
                    row[col] = sc_lookup.loc[spot_id, col]
        rows.append(row)
    print(f"  [{sid}] {len(rows)} tiles (skipped border={skipped})")
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--registry", default="project_external_st/results/registry/dataset_sample_registry.tsv")
    ap.add_argument("--scores-dir", default="project_external_st/results/scores")
    ap.add_argument("--tiles-dir", default="project_external_st/data/processed/tiles_external")
    ap.add_argument("--meta-out", default="project_external_st/results/extra/g1_external_tile_metadata.tsv.gz")
    ap.add_argument("--tile-size", type=int, default=224)
    ap.add_argument("--max-per-sample", type=int, default=200)
    args = ap.parse_args()

    reg = pd.read_csv(args.registry, sep="\t")
    reg = reg[reg["status"] == "ok"]
    rows = []
    for _, r in reg.iterrows():
        sid = r["sample_id"]; ds = r["dataset"]
        a = ad.read_h5ad(r["h5ad"])
        spot_tsv = Path(args.scores_dir) / f"{ds}_{sid}_spot_scores.tsv.gz"
        score_df = pd.read_csv(spot_tsv, sep="\t") if spot_tsv.exists() else None
        out_dir = Path(args.tiles_dir) / sid
        df = extract_one(a, args.tile_size, out_dir, args.max_per_sample, score_df)
        if not df.empty:
            rows.append(df)
    if not rows:
        print("no tiles"); return
    meta = pd.concat(rows, ignore_index=True)
    Path(args.meta_out).parent.mkdir(parents=True, exist_ok=True)
    meta.to_csv(args.meta_out, sep="\t", index=False, compression="gzip")
    print(f"\nMetadata: {len(meta)} tiles → {args.meta_out}")
    print(meta["dataset"].value_counts())


if __name__ == "__main__":
    main()
