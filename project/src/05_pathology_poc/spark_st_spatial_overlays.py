#!/usr/bin/env python3
"""
Per-slide spatial overlay figures — DM1_like_score + SPARK B5 cluster + IGLC2 expression
overlaid on H&E hires thumbnail. 16 slides × 4-panel = 64 figures.

Inputs:
  project/data/processed/GSE250521/{sid}/{sid}.scored.h5ad     (Visium with obsm['spatial'])
  project/data/processed/GSE250521/{sid}/spatial/tissue_hires_image.png
  project/results/03_pathology_poc/spark_analytical_per_tile.tsv.gz

Output:
  project/results/03_pathology_poc/spatial_overlays/{sid}.png      (4-panel paper figure)
  project/results/03_pathology_poc/spatial_overlays/index.html
"""
from __future__ import annotations
import sys
from pathlib import Path

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image

ROOT = Path(__file__).resolve().parent.parent.parent.parent
GS = ROOT / "project/data/processed/GSE250521"
SPARK_TILE = ROOT / "project/results/03_pathology_poc/spark_analytical_per_tile.tsv.gz"
OUT_DIR = ROOT / "project/results/03_pathology_poc/spatial_overlays"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def get_scalefactor(spatial_uns):
    if not spatial_uns:
        return 1.0, None
    libs = list(spatial_uns.keys())
    sf = spatial_uns[libs[0]].get("scalefactors", {})
    return float(sf.get("tissue_hires_scalef", 1.0)), libs[0]


def hires_path(sdir):
    cands = [sdir / "spatial/tissue_hires_image.png",
             sdir / f"{sdir.name}_tissue_hires_image.png"]
    for c in cands:
        if c.exists(): return c
    cands = list(sdir.rglob("tissue_hires_image.png"))
    return cands[0] if cands else None


def render_slide(sid, spark):
    sdir = GS / sid
    h5 = sdir / f"{sid}.scored.h5ad"
    if not h5.exists():
        return None
    a = ad.read_h5ad(h5)
    sf, lib = get_scalefactor(a.uns.get("spatial", {}))
    bg = hires_path(sdir)
    if bg is None:
        return None
    bg_img = np.array(Image.open(bg).convert("RGB"))

    # spot coords in hires pixel space
    coords = a.obsm["spatial"] * sf
    spark_s = spark[spark.sample_id == sid].set_index("spot_id")
    obs = a.obs.reset_index().rename(columns={"index": "spot_id"})
    obs = obs.merge(spark_s[["B5_thyrocyte_cluster_med", "B4_tumor_stroma_ratio"]],
                    left_on="spot_id", right_index=True, how="left")
    if "spot_id" not in obs.columns:
        obs.index.name = "spot_id"
        obs = obs.reset_index()

    # IGLC2 expression
    iglc2_idx = np.where(a.var.index == "IGLC2")[0]
    if len(iglc2_idx):
        x = a.X
        if hasattr(x, "toarray"):
            x = x.toarray()
        iglc2 = np.asarray(x[:, iglc2_idx[0]]).flatten()
    else:
        iglc2 = np.zeros(a.n_obs)

    fig, axes = plt.subplots(1, 4, figsize=(20, 5.5), constrained_layout=True)
    panels = [
        ("H&E", None, None, None),
        ("DM1_like_score", a.obs["DM1_like_score"].values, "RdBu_r", "DM1"),
        ("SPARK B5 thyrocyte_cluster_med", obs["B5_thyrocyte_cluster_med"].values, "viridis", "B5"),
        ("IGLC2 (plasma cell)", iglc2, "magma", "IGLC2"),
    ]
    for ax, (title, vals, cmap, _) in zip(axes, panels):
        ax.imshow(bg_img)
        ax.set_title(f"{sid} · {title}", fontsize=10)
        ax.axis("off")
        if vals is not None:
            keep = ~np.isnan(vals) if vals.dtype.kind == "f" else np.ones_like(vals, dtype=bool)
            v = vals[keep]
            xy = coords[keep]
            if v.size and np.nanstd(v) > 0:
                lo, hi = np.nanpercentile(v, [5, 95])
                sc = ax.scatter(xy[:, 0], xy[:, 1], c=v, cmap=cmap,
                                s=4, alpha=0.75, vmin=lo, vmax=hi, edgecolor="none")
                fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.02)
    out = OUT_DIR / f"{sid}.png"
    fig.savefig(out, dpi=140, bbox_inches="tight")
    plt.close(fig)
    return out


def main():
    spark = pd.read_csv(SPARK_TILE, sep="\t")
    sids = sorted([d.name for d in GS.iterdir() if d.is_dir() and (d / f"{d.name}.scored.h5ad").exists()])
    print(f"slides: {len(sids)}")
    out_paths = []
    for sid in sids:
        p = render_slide(sid, spark)
        print(f"  {sid}: {'ok ' + str(p.relative_to(ROOT)) if p else 'skip'}")
        if p:
            out_paths.append(p)

    # Build a small index.html gallery
    html = ['<!doctype html><html><head><meta charset="utf-8"><title>SPARK · DM spatial overlays</title>',
            '<style>body{background:#0a0d18;color:#e1e8f0;font-family:Inter,system-ui;margin:24px}',
            'h1{margin:0 0 8px} .sub{color:#94a3b8;margin-bottom:24px}',
            'figure{margin:0 0 18px;border:1px solid #1f2a44;border-radius:10px;padding:8px;background:#0e1426}',
            'figcaption{font-size:12px;color:#94a3b8;padding:6px 4px}',
            'img{width:100%;display:block;border-radius:6px}</style></head><body>',
            '<h1>SPARK · DM spatial overlays (16 slides × 4 panels)</h1>',
            '<div class="sub">H&E + DM1_like_score + SPARK B5 thyrocyte_cluster_med + IGLC2 plasma marker · per-spot Visium overlay</div>']
    for p in out_paths:
        rel = p.name
        html.append(f'<figure><img src="{rel}" /><figcaption>{p.stem}</figcaption></figure>')
    html.append("</body></html>")
    (OUT_DIR / "index.html").write_text("\n".join(html))
    print(f"\nwrote gallery {OUT_DIR / 'index.html'}")
    print(f"  ({len(out_paths)} PNG figures)")


if __name__ == "__main__":
    main()
