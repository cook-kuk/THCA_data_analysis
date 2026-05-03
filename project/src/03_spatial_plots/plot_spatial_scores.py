#!/usr/bin/env python3
"""Per-sample spatial heatmaps for RAI_8, DM1_like, TDS_like + 2x4 advisor figure."""
from __future__ import annotations
import argparse, logging, sys
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger("plot")
Image.MAX_IMAGE_PIXELS = None

STAGE_ORDER = ["PT", "PTC", "LPTC", "ATC"]


def load_image(adata: ad.AnnData) -> tuple[np.ndarray | None, float]:
    sid = adata.obs["sample_id"].iloc[0]
    spatial = adata.uns.get("spatial", {}).get(sid, {})
    sf = spatial.get("scalefactors", {})
    hires_sf = sf.get("tissue_hires_scalef", 1.0)
    img_path = spatial.get("images", {}).get("hires_path")
    if img_path and Path(img_path).exists():
        return np.asarray(Image.open(img_path)), hires_sf
    return None, hires_sf


def plot_one(adata: ad.AnnData, score: str, ax, *, vmin=None, vmax=None, title=None, cmap="RdBu_r"):
    img, sf = load_image(adata)
    x = adata.obsm["spatial"][:, 0] * sf
    y = adata.obsm["spatial"][:, 1] * sf
    vals = adata.obs[score].to_numpy()
    if vmin is None: vmin = np.nanpercentile(vals, 2)
    if vmax is None: vmax = np.nanpercentile(vals, 98)
    if vmin == vmax: vmin, vmax = vmin - 1, vmax + 1
    if img is not None:
        ax.imshow(img, alpha=0.55)
    sc = ax.scatter(x, y, c=vals, s=4, cmap=cmap, vmin=vmin, vmax=vmax,
                    linewidths=0, edgecolors="none")
    ax.set_aspect("equal"); ax.invert_yaxis()
    ax.set_xticks([]); ax.set_yticks([])
    if title: ax.set_title(title, fontsize=9)
    plt.colorbar(sc, ax=ax, fraction=0.04, pad=0.02)


def per_sample_plots(scored_dir: Path, out_dir: Path) -> list[Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = []
    for h in sorted(scored_dir.rglob("*.scored.h5ad")):
        adata = ad.read_h5ad(h)
        sid = adata.obs["sample_id"].iloc[0]
        stage = adata.obs["stage"].iloc[0]
        fig, axes = plt.subplots(1, 3, figsize=(13, 4.5))
        for ax, sc in zip(axes, ["RAI_8_score", "DM1_like_score", "TDS_like_score"]):
            plot_one(adata, sc, ax, title=f"{sid} [{stage}]\n{sc}")
        fig.tight_layout()
        out = out_dir / f"{sid}_spatial_RAI_DM1_TDS.png"
        fig.savefig(out, dpi=150, bbox_inches="tight")
        plt.close(fig)
        log.info("[%s] → %s", sid, out)
        paths.append(out)
    return paths


def advisor_figure(scored_dir: Path, out_path: Path) -> Path | None:
    """2 rows (RAI_8, DM1_like) × 4 cols (PT, PTC, LPTC, ATC) using one rep per stage."""
    by_stage: dict[str, ad.AnnData] = {}
    for h in sorted(scored_dir.rglob("*.scored.h5ad")):
        a = ad.read_h5ad(h)
        st = a.obs["stage"].iloc[0]
        if st not in STAGE_ORDER: continue
        if st not in by_stage:  # first sample per stage = rep
            by_stage[st] = a
    if not by_stage:
        return None
    fig, axes = plt.subplots(2, 4, figsize=(18, 10))
    for j, st in enumerate(STAGE_ORDER):
        if st not in by_stage:
            for r in range(2):
                axes[r, j].set_axis_off()
                axes[r, j].set_title(f"{st}: no sample", fontsize=10)
            continue
        a = by_stage[st]
        sid = a.obs["sample_id"].iloc[0]
        plot_one(a, "RAI_8_score", axes[0, j], title=f"{st} | {sid}\nRAI_8_score")
        plot_one(a, "DM1_like_score", axes[1, j], title=f"DM1_like_score (= -RAI_8)")
    fig.suptitle("GSE250521 spatial validation: RAI_8 (row1) vs DM1-like (row2) across PT→PTC→LPTC→ATC",
                 fontsize=12, y=0.995)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    fig.savefig(out_path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    log.info("advisor 2x4 → %s", out_path)
    return out_path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--processed-dir", default="project/data/processed/GSE250521")
    ap.add_argument("--per-sample-dir", default="project/results/01_spatial_score/per_sample")
    ap.add_argument("--advisor-out", default="project/results/figures_for_advisor/fig_2x4_RAI_DM1_per_stage.png")
    args = ap.parse_args()

    per_sample_plots(Path(args.processed_dir), Path(args.per_sample_dir))
    advisor_figure(Path(args.processed_dir), Path(args.advisor_out))


if __name__ == "__main__":
    sys.exit(main())
