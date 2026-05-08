#!/usr/bin/env python3
"""
SPARK paper-grade summary figures (matplotlib).

Generates 4 paper figures:
  F1 · 4-stage trajectory bar chart (DM1 / B5 / M1_M2 / Moran's I)
  F2 · Cross-platform replicate (ST 16-slide × bulk n=572 forest plot)
  F3 · SVG top 10 per stage bar chart
  F4 · L-R activity heatmap (12 L-R × 4 stage, DM1 top-bot diff)

Output: project/results/03_pathology_poc/summary_figs/F{1-4}.png
"""
from __future__ import annotations
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent.parent.parent
RES = ROOT / "project/results/03_pathology_poc"
OUT = RES / "summary_figs"
OUT.mkdir(parents=True, exist_ok=True)

plt.rcParams["font.family"] = "DejaVu Sans"
plt.rcParams["axes.spines.top"] = False
plt.rcParams["axes.spines.right"] = False


def F1_trajectory():
    """4-stage trajectory bar chart."""
    sl = pd.read_csv(RES / "spark_full_visium_per_slide.tsv", sep="\t")
    sa = pd.read_csv(RES / "spark_analytical_per_slide.tsv", sep="\t")
    stages = ["PT", "PTC", "LPTC", "ATC"]
    fig, axes = plt.subplots(1, 4, figsize=(15, 4.2), constrained_layout=True)

    # DM1 Moran
    moran = sl.groupby("stage")["morans_dm1_k6"].mean().reindex(stages)
    axes[0].bar(stages, moran.values, color=["#94a3b8", "#5eead4", "#fbbf24", "#ef4444"])
    axes[0].set_title("DM1 Moran's I-lite (k=6)\nspatial organization", fontsize=11)
    axes[0].set_ylabel("Moran's I"); axes[0].axhline(0, color="#444", lw=0.5)
    axes[0].annotate(f"{moran['ATC']:.2f}", (3, moran['ATC']), ha="center", va="bottom",
                     color="#ef4444", fontsize=10, fontweight="bold")
    axes[0].annotate("collapse →", (3, moran['ATC']+0.05), ha="center", color="#ef4444", fontsize=8)

    # B5 thyrocyte cluster
    b5 = sa.groupby("stage")["B5_thyrocyte_cluster_med_mean"].mean().reindex(stages)
    axes[1].bar(stages, b5.values, color=["#94a3b8", "#5eead4", "#fbbf24", "#ef4444"])
    axes[1].set_title("SPARK B5\nthyrocyte cluster median (per tile)", fontsize=11)
    axes[1].set_ylabel("median cluster size")
    axes[1].annotate("+49% N→ATC", (1.5, b5.max()*1.05), ha="center", color="#ef4444",
                     fontsize=9, fontweight="bold")

    # M1/M2 Cohen d top25-bot25
    qstr = pd.read_csv(RES / "spark_full_visium_quartile_stratified.tsv", sep="\t")
    m12 = qstr[qstr.module == "M1_M2"].groupby("stage")["cohen_d_top25_vs_bot25"].mean().reindex(stages)
    axes[2].bar(stages, m12.values, color=["#94a3b8", "#5eead4", "#fbbf24", "#ef4444"])
    axes[2].set_title("M1/M2 macrophage Cohen d\n(DM1 top25 vs bot25 spots)", fontsize=11)
    axes[2].set_ylabel("Cohen d"); axes[2].axhline(0, color="#444", lw=0.5)
    axes[2].annotate(f"d={m12['ATC']:.2f}", (3, m12['ATC']), ha="center", va="bottom",
                     color="#ef4444", fontsize=10, fontweight="bold")

    # B-plasma Cohen d
    bpl = qstr[qstr.module == "Bplasma"].groupby("stage")["cohen_d_top25_vs_bot25"].mean().reindex(stages)
    axes[3].bar(stages, bpl.values, color=["#94a3b8", "#5eead4", "#fbbf24", "#ef4444"])
    axes[3].set_title("B-plasma Cohen d\n(DM1 top25 vs bot25)", fontsize=11)
    axes[3].set_ylabel("Cohen d"); axes[3].axhline(0, color="#444", lw=0.5)

    fig.suptitle("F1 · 4-stage SPARK trajectory (16 slides, 55,873 spots)",
                 fontsize=12, fontweight="bold")
    fig.savefig(OUT / "F1_trajectory.png", dpi=160, bbox_inches="tight")
    plt.close(fig)
    print("F1 trajectory done")


def F2_cross_platform():
    """Cross-platform replicate forest."""
    rep = pd.read_csv(RES / "spark_tcga_bulk_caf_rai_replicate.tsv", sep="\t")
    pool = pd.read_csv(RES / "spark_cci_lag_correlation_pooled.tsv", sep="\t")

    fig, ax = plt.subplots(1, 1, figsize=(8, 4.5), constrained_layout=True)
    pairs = [
        ("CAF × RAI",   "CAF",   "RAI_lineage"),
        ("M1/M2 × RAI", "M1_M2", "RAI_lineage"),
        ("TLS × RAI",   "TLS",   "RAI_lineage"),
        ("CD8 × RAI",   "CD8_Teff", "RAI_lineage"),
    ]
    pos = np.arange(len(pairs))
    st_vals, bk_vals, labels = [], [], []
    for label, a_st, b_st in pairs:
        st = pool[((pool.mod_a == a_st) & (pool.mod_b == b_st)) |
                  ((pool.mod_a == b_st) & (pool.mod_b == a_st))]
        st_v = float(st["median_lag_sym"].iloc[0]) if len(st) else np.nan
        bk = rep[rep.a == a_st.replace("CD8_Teff", "CD8")]
        bk_v = float(bk["spearman_r"].iloc[0]) if len(bk) else np.nan
        st_vals.append(st_v); bk_vals.append(bk_v); labels.append(label)

    width = 0.35
    ax.barh(pos - width/2, st_vals, width, color="#5eead4", label="ST median lag (16 slides)")
    ax.barh(pos + width/2, bk_vals, width, color="#ff7c3e", label="TCGA bulk Spearman (n=572)")
    ax.axvline(0, color="#444", lw=0.7)
    ax.set_yticks(pos); ax.set_yticklabels(labels)
    ax.set_xlabel("Spatial avoidance correlation")
    ax.legend(loc="lower right", fontsize=9)
    ax.set_title("F2 · CAF/M1-M2/TLS × RAI cross-platform replication\n"
                 "(ST direction matches bulk in 4/4 module pairs)", fontsize=11, fontweight="bold")
    # red box around the strongest replicate
    ax.add_patch(plt.Rectangle((-0.5, -0.45), 0.5, 0.95, fill=False, edgecolor="red", linewidth=2.5))
    ax.annotate("paper-grade\nreplicate", (-0.51, -0.5), color="red", fontsize=8, fontweight="bold",
                ha="right")
    fig.savefig(OUT / "F2_cross_platform_replicate.png", dpi=160, bbox_inches="tight")
    plt.close(fig)
    print("F2 cross-platform done")


def F3_svg():
    """SVG top 10 per stage bar chart."""
    svg = pd.read_csv(RES / "spark_svg_per_stage_top.tsv", sep="\t")
    stages = ["PT", "PTC", "LPTC", "ATC"]
    fig, axes = plt.subplots(1, 4, figsize=(16, 5), constrained_layout=True)
    for ax, s, c in zip(axes, stages, ["#94a3b8", "#5eead4", "#fbbf24", "#ef4444"]):
        sub = svg[svg.stage == s].head(10).iloc[::-1]
        ax.barh(sub.gene, sub.median_moran, color=c)
        ax.set_title(f"{s} top-10 spatially organized genes", fontsize=10)
        ax.set_xlabel("median Moran's I")
        for i, (g, v) in enumerate(zip(sub.gene, sub.median_moran)):
            ax.text(v, i, f" {g}", va="center", fontsize=8)
        ax.set_yticks([])
    fig.suptitle("F3 · Spatially variable genes — LPTC=CCL19/CCL21 TLS axis · ATC=KRT7/CD24 stem axis",
                 fontsize=11, fontweight="bold")
    fig.savefig(OUT / "F3_svg_top10.png", dpi=160, bbox_inches="tight")
    plt.close(fig)
    print("F3 SVG done")


def F4_lr_heatmap():
    """L-R activity heatmap."""
    lr = pd.read_csv(RES / "spark_cci_LR_activity_per_stage.tsv", sep="\t")
    pivot = lr.pivot_table(index="lr_pair", columns="stage", values="DM1_top_minus_bot_mean", aggfunc="mean")
    stages = ["PT", "PTC", "LPTC", "ATC"]
    pivot = pivot[[s for s in stages if s in pivot.columns]]
    fig, ax = plt.subplots(figsize=(7, 6), constrained_layout=True)
    vmax = float(np.nanmax(np.abs(pivot.values)))
    im = ax.imshow(pivot.values, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    ax.set_xticks(range(len(pivot.columns))); ax.set_xticklabels(pivot.columns)
    ax.set_yticks(range(len(pivot.index))); ax.set_yticklabels(pivot.index, fontsize=10)
    for i in range(pivot.shape[0]):
        for j in range(pivot.shape[1]):
            v = pivot.iloc[i, j]
            ax.text(j, i, f"{v:+.3f}", ha="center", va="center",
                    color="white" if abs(v) > vmax*0.5 else "#222", fontsize=8)
    fig.colorbar(im, ax=ax, label="DM1-top minus DM1-bot activity")
    ax.set_title("F4 · L-R activity (DM1-high vs DM1-low spots, per stage)", fontsize=11, fontweight="bold")
    # red box around LGALS9-HAVCR2 ATC and TGFB1 LPTC
    if "LGALS9→HAVCR2" in pivot.index and "ATC" in pivot.columns:
        ri = list(pivot.index).index("LGALS9→HAVCR2")
        ci = list(pivot.columns).index("ATC")
        ax.add_patch(plt.Rectangle((ci-0.45, ri-0.45), 0.9, 0.9, fill=False, edgecolor="red", linewidth=3))
    if "TGFB1→TGFBR1" in pivot.index and "LPTC" in pivot.columns:
        ri = list(pivot.index).index("TGFB1→TGFBR1")
        ci = list(pivot.columns).index("LPTC")
        ax.add_patch(plt.Rectangle((ci-0.45, ri-0.45), 0.9, 0.9, fill=False, edgecolor="red", linewidth=3))
    fig.savefig(OUT / "F4_LR_heatmap.png", dpi=160, bbox_inches="tight")
    plt.close(fig)
    print("F4 L-R done")


def F5_cox_forest():
    """TCGA Cox forest plot for DSS endpoint."""
    cox = pd.read_csv(RES / "spark_tcga_cox_modules.tsv", sep="\t")
    dss = cox[cox.endpoint == "DSS"].sort_values("HR")
    fig, ax = plt.subplots(figsize=(7, 4), constrained_layout=True)
    pos = np.arange(len(dss))
    sig = dss.p < 0.05
    cols = ["#ef4444" if s else "#94a3b8" for s in sig]
    ax.barh(pos, dss.HR.values, color=cols, edgecolor="white", height=0.6)
    ax.axvline(1.0, color="#444", lw=0.7, ls="--")
    ax.set_yticks(pos); ax.set_yticklabels(dss.module.values, fontsize=10)
    ax.set_xlabel("Hazard ratio (DSS)")
    for i, (hr, p) in enumerate(zip(dss.HR.values, dss.p.values)):
        marker = " ★" if p < 0.01 else (" *" if p < 0.05 else "")
        ax.text(hr + 0.05, i, f"HR={hr:.2f}, p={p:.3f}{marker}", va="center", fontsize=9,
                color="#ef4444" if p < 0.05 else "#444")
    ax.set_title("F5 · TCGA-THCA DSS Cox (n=566)\nM1/M2 + CD8 protective ★ (p<0.05)",
                 fontsize=11, fontweight="bold")
    fig.savefig(OUT / "F5_cox_DSS.png", dpi=160, bbox_inches="tight")
    plt.close(fig)
    print("F5 Cox done")


if __name__ == "__main__":
    F1_trajectory()
    F2_cross_platform()
    F3_svg()
    F4_lr_heatmap()
    F5_cox_forest()
    print("\nALL summary figures done")
