#!/usr/bin/env python3
"""Diagnostic plots: predicted vs observed for primary target across folds."""
from __future__ import annotations
import argparse, logging
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger("plot")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--preds", required=True)
    ap.add_argument("--metrics", required=True)
    ap.add_argument("--target", default="DM1_like_score_resid")
    ap.add_argument("--model", default="ridge")
    ap.add_argument("--out", default="project/results/03_pathology_poc/pred_vs_obs_resnet50.png")
    args = ap.parse_args()

    preds = pd.read_csv(args.preds, sep="\t")
    metrics = pd.read_csv(args.metrics, sep="\t")
    df = preds[(preds.target == args.target) & (preds.model == args.model)].dropna()
    if df.empty:
        log.error("no predictions"); raise SystemExit(2)

    sizes = sorted(df.tile_size.unique())
    fig, axes = plt.subplots(1, len(sizes), figsize=(5 * len(sizes), 5), squeeze=False)
    for ax, ts in zip(axes[0], sizes):
        sub = df[df.tile_size == ts]
        sp = spearmanr(sub.y_obs, sub.y_pred).correlation
        ax.scatter(sub.y_obs, sub.y_pred, s=3, alpha=0.3, c="steelblue")
        lim = [sub.y_obs.min(), sub.y_obs.max()]
        ax.plot(lim, lim, "k--", lw=0.7)
        ax.set_xlabel(f"observed {args.target}")
        ax.set_ylabel("predicted")
        ax.set_title(f"tile_size={ts}\npooled Spearman={sp:.3f}  n={len(sub)}")
    fig.suptitle(f"{args.model} LOSO  ({args.target})", y=1.02)
    fig.tight_layout()
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    log.info("wrote %s", out)

    # per-slide spearman bar
    fig2, ax2 = plt.subplots(1, len(sizes), figsize=(6 * len(sizes), 4), squeeze=False)
    for ax, ts in zip(fig2.axes, sizes):
        sub = df[df.tile_size == ts]
        per_slide = sub.groupby("sample_id").apply(
            lambda g: spearmanr(g.y_obs, g.y_pred).correlation if len(g) >= 5 else np.nan
        ).sort_values()
        ax.barh(per_slide.index, per_slide.values, color="steelblue")
        ax.axvline(0, color="k", lw=0.5)
        ax.set_xlabel("per-slide LOSO Spearman r")
        ax.set_title(f"tile_size={ts}")
    fig2.tight_layout()
    out2 = out.with_name(out.stem + "_per_slide.png")
    fig2.savefig(out2, dpi=150, bbox_inches="tight")
    log.info("wrote %s", out2)


if __name__ == "__main__":
    main()
