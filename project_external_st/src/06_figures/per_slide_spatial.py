#!/usr/bin/env python3
"""Per-slide spatial heatmaps for external datasets — RAI_8, DM1_like, THYROID_NONOVERLAP, Epithelial.
Reuses scored spot tables (already computed) — no new scoring."""
from __future__ import annotations
import argparse, sys
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image

Image.MAX_IMAGE_PIXELS = None

SCORES = ["RAI_8_score_resid", "DM1_like_score_resid", "THYROID_NONOVERLAP_score_resid", "Epithelial_score_resid"]
TITLES = {"RAI_8_score_resid":"RAI_8 (depth-resid)",
          "DM1_like_score_resid":"DM1_like (depth-resid)",
          "THYROID_NONOVERLAP_score_resid":"THYROID_NONOVERLAP (depth-resid)",
          "Epithelial_score_resid":"Epithelial mask"}
CMAPS = {"RAI_8_score_resid":"RdBu_r","DM1_like_score_resid":"RdBu_r",
         "THYROID_NONOVERLAP_score_resid":"RdBu_r","Epithelial_score_resid":"viridis"}


def plot_one(ax, df, score, *, img=None, sf=1.0, title=None):
    x = df["pxl_col_in_fullres"].to_numpy() * sf
    y = df["pxl_row_in_fullres"].to_numpy() * sf
    v = df[score].to_numpy()
    vmin, vmax = np.nanpercentile(v, 2), np.nanpercentile(v, 98)
    if vmin == vmax: vmin, vmax = vmin - 1, vmax + 1
    if img is not None:
        ax.imshow(img, alpha=0.55)
    sc = ax.scatter(x, y, c=v, s=4, cmap=CMAPS[score], vmin=vmin, vmax=vmax,
                    linewidths=0, edgecolors="none")
    ax.set_aspect("equal"); ax.invert_yaxis()
    ax.set_xticks([]); ax.set_yticks([])
    if title: ax.set_title(title, fontsize=9)
    plt.colorbar(sc, ax=ax, fraction=0.04, pad=0.02)


def load_image_and_sf(meta_row):
    """Load hires image + scalefactor from registry. GSE230424 has only HE.jpg, no scalefactors."""
    img_path = meta_row.get("image", "")
    img = None
    sf = 1.0
    if img_path and Path(img_path).exists():
        try: img = np.asarray(Image.open(img_path))
        except Exception: img = None
    # Scalefactor: read from .h5ad uns if present (GSE248205); GSE230424 has none → sf=1.0 (raw pixel coords)
    h5ad_path = meta_row.get("h5ad")
    if h5ad_path and Path(h5ad_path).exists():
        try:
            import anndata as ad
            a = ad.read_h5ad(h5ad_path, backed="r")
            sample_id = a.obs["sample_id"].iloc[0]
            sf = a.uns.get("spatial", {}).get(sample_id, {}).get("scalefactors", {}).get("tissue_hires_scalef", 1.0)
            a.file.close()
        except Exception: pass
    return img, sf


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--registry", default="project_external_st/results/registry/dataset_sample_registry.tsv")
    ap.add_argument("--scores-dir", default="project_external_st/results/scores")
    ap.add_argument("--out-dir", default="project_external_st/results/figures/per_slide")
    args = ap.parse_args()
    out = Path(args.out_dir); out.mkdir(parents=True, exist_ok=True)

    reg = pd.read_csv(args.registry, sep="\t")
    for _, r in reg[reg["status"] == "ok"].iterrows():
        sid = r["sample_id"]; ds = r["dataset"]
        spot_tsv = Path(args.scores_dir) / f"{ds}_{sid}_spot_scores.tsv.gz"
        if not spot_tsv.exists():
            print(f"[skip] no spot file for {sid}"); continue
        df = pd.read_csv(spot_tsv, sep="\t")
        img, sf = load_image_and_sf(r)
        fig, axes = plt.subplots(1, 4, figsize=(18, 5))
        for ax, sc in zip(axes, SCORES):
            plot_one(ax, df, sc, img=img, sf=sf,
                     title=f"{sid}\n[{r['condition_inferred']}] {TITLES[sc]}")
        fig.suptitle(f"{ds} · {sid} · {r['condition_inferred']} · {df.shape[0]} spots",
                     fontsize=11, y=1.02)
        fig.tight_layout()
        out_png = out / f"{ds}_{sid}_4panel.png"
        fig.savefig(out_png, dpi=130, bbox_inches="tight")
        plt.close(fig)
        print(f"→ {out_png}")


if __name__ == "__main__":
    sys.exit(main())
