#!/usr/bin/env python3
"""8-condition DM1_like representative mosaic — one slide per condition.
Conditions: PT, PTC, LPTC, ATC (GSE250521) + CONTROL, HT, GD (GSE248205) + PTC_HT (GSE230424)."""
from __future__ import annotations
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image

Image.MAX_IMAGE_PIXELS = None

CONDITIONS = [("PT","GSE250521"),("PTC","GSE250521"),("LPTC","GSE250521"),("ATC","GSE250521"),
              ("CONTROL","GSE248205"),("HT","GSE248205"),("GD","GSE248205"),("PTC_HT","GSE230424")]


def load_data(condition: str, dataset: str):
    if dataset == "GSE250521":
        df = pd.read_csv("project/results/01_spatial_score/all_spots_scored.tsv.gz", sep="\t")
        sub = df[df["stage"] == condition]
        sid = sub["sample_id"].iloc[0]
        sub = sub[sub["sample_id"] == sid].copy()
        sub["DM1_like_resid"] = sub["DM1_like_score"]  # GSE250521 only has raw, use as proxy for visualization
        sub["score"] = sub["DM1_like_resid"]
        return sub, sid, "DM1_like (within-sample z; depth-resid not in source data)"
    else:
        # external: pick first slide of this condition
        ext = pd.read_csv("project_external_st/results/scores/all_external_spots_scored.tsv.gz", sep="\t")
        sub = ext[(ext["condition_inferred"] == condition) & (ext["dataset"] == dataset)]
        sid = sub["sample_id"].iloc[0]
        sub = sub[sub["sample_id"] == sid].copy()
        sub["score"] = sub["DM1_like_score_resid"]
        return sub, sid, "DM1_like (depth-resid)"


def main():
    fig, axes = plt.subplots(2, 4, figsize=(18, 10))
    for ax, (cond, ds) in zip(axes.flat, CONDITIONS):
        df, sid, note = load_data(cond, ds)
        x = df["pxl_col_in_fullres"].to_numpy()
        y = df["pxl_row_in_fullres"].to_numpy()
        v = df["score"].to_numpy()
        vmin, vmax = np.nanpercentile(v, 2), np.nanpercentile(v, 98)
        if vmin == vmax: vmin, vmax = vmin - 1, vmax + 1
        sc = ax.scatter(x, y, c=v, s=4, cmap="RdBu_r", vmin=vmin, vmax=vmax,
                        linewidths=0, edgecolors="none")
        ax.set_aspect("equal"); ax.invert_yaxis()
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(f"{cond} · {ds}\n{sid}\n{note}", fontsize=9)
        plt.colorbar(sc, ax=ax, fraction=0.04, pad=0.02)
    fig.suptitle("DM1_like spatial distribution — 8 conditions (one representative slide each)\n"
                 "row 1: cancer progression (PT to ATC); row 2: non-cancer + Hashimoto-overlap controls",
                 fontsize=12, y=1.0)
    fig.tight_layout()
    out = Path("project_external_st/results/figures/condition_dm1_mosaic_8panel.png")
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=140, bbox_inches="tight")
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    print(f"→ {out}")


if __name__ == "__main__":
    main()
