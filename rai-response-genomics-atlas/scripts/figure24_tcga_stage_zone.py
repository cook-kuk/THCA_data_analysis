#!/usr/bin/env python3
"""Figure 24 — TCGA-THCA AJCC stage × R17 zone composition.

Per-stage stacked bar of zone fraction. Demonstrates that dark-matter and
BRAF-like zones increase progressively with stage.
"""
from __future__ import annotations
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent))
from _nature_style import setup_rc, panel_letter, panel_title, IVORY, INK, MUTED, BLUE, RED

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

setup_rc()

ROOT = Path(__file__).resolve().parent.parent
TCGA_MERGED = Path("/home/seungho/personal/THCA_data_analysis/project/results/r17_tcga_panel_d4p2_reconciliation/r17_per_sample_merged.tsv")
CLIN = Path("/home/seungho/personal/THCA_data_analysis/project/results/dark_matter_phase1/tcga_dm_master_with_pfi.tsv")
OUT_PNG = ROOT / "results" / "figures" / "figure24_tcga_stage_zone.png"
OUT_PDF = ROOT / "results" / "figures" / "figure24_tcga_stage_zone.pdf"

ZONE_ORDER = ["WT-like", "RAS-like", "BRAF-like", "dark-matter"]
ZONE_COLOR = {"WT-like": "#cdc6c0", "RAS-like": RED, "BRAF-like": BLUE, "dark-matter": "#5a4470"}


def stage_simplify(s):
    if pd.isna(s): return None
    s = str(s).strip().lower()
    if "iv" in s: return "IV"
    if "iii" in s: return "III"
    if "ii" in s: return "II"
    if "i" in s and "i" not in s.replace("stage i", "", 1): return "I"
    if "stage i" in s: return "I"
    return None


def main():
    df = pd.read_csv(TCGA_MERGED, sep="\t").dropna(subset=["RAI_8"])
    clin = pd.read_csv(CLIN, sep="\t")

    df["pat"] = df["sampleId"].str.slice(0, 12)
    df = df.merge(clin[["tcga_short", "stage", "ajcc_pathologic_tumor_stage"]],
                  left_on="pat", right_on="tcga_short", how="left")

    # Build zone
    panel = df["panel_DM"].fillna("").astype(str)
    d4    = df["d4p2_DM"].fillna("").astype(str)
    zone = pd.Series("WT-like", index=df.index)
    zone[(panel == "DM1") & (d4 == "DM2")] = "BRAF-like"
    zone[(panel == "DM2") & (d4 == "DM1")] = "RAS-like"
    zone[(panel == "DM1") & (d4 == "DM1")] = "dark-matter"
    zone[(panel == "DM2") & (d4 == "DM2")] = "WT-like"
    df["zone"] = zone

    df["stage_simple"] = df["ajcc_pathologic_tumor_stage"].apply(stage_simplify)
    df = df.dropna(subset=["stage_simple"])
    stages = ["I", "II", "III", "IV"]
    stages = [s for s in stages if (df["stage_simple"] == s).sum() > 0]

    # Compute per-stage zone fractions
    counts = df.groupby(["stage_simple", "zone"]).size().unstack(fill_value=0)
    counts = counts.reindex(index=stages, columns=ZONE_ORDER, fill_value=0)
    fracs = counts.div(counts.sum(axis=1), axis=0)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.6), facecolor=IVORY,
                              gridspec_kw=dict(wspace=0.32, left=0.07, right=0.97, top=0.83, bottom=0.18))
    fig.suptitle("TCGA-THCA AJCC stage × R17 zone composition",
                 fontsize=13, fontweight="bold", color=INK, y=0.97)

    # Left: stacked bar
    ax = axes[0]
    bottom = np.zeros(len(stages))
    x = np.arange(len(stages))
    for c in ZONE_ORDER:
        vals = fracs[c].values
        ax.bar(x, vals, bottom=bottom, color=ZONE_COLOR[c], label=c, width=0.6,
               edgecolor="white", lw=0.7)
        for xi, v, b in zip(x, vals, bottom):
            if v > 0.06:
                ax.text(xi, b + v/2, f"{v*100:.0f}%", ha="center", va="center",
                        fontsize=10, color="white", fontweight="bold")
        bottom += vals
    ax.set_xticks(x); ax.set_xticklabels([f"Stage {s}\n(n={int(counts.iloc[i].sum())})" for i, s in enumerate(stages)], fontsize=10)
    ax.set_ylabel("zone fraction")
    ax.set_ylim(0, 1.02)
    ax.legend(fontsize=8.8, loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=4, frameon=False)
    panel_title(ax, "Stage × zone composition"); panel_letter(ax, "a")

    # Right: panel z by stage boxplot
    ax = axes[1]
    data = [df.loc[df["stage_simple"] == s, "RAI_8"].dropna().values for s in stages]
    bp = ax.boxplot(data, positions=range(len(stages)), widths=0.55,
                     patch_artist=True, medianprops=dict(color=INK, lw=1.3),
                     boxprops=dict(lw=0.7), whiskerprops=dict(lw=0.7), capprops=dict(lw=0.7),
                     flierprops=dict(marker="o", markersize=2, alpha=0.4))
    cols = ["#7a8b9c", "#5e7891", "#456478", "#2d4659"]
    for patch, col in zip(bp["boxes"], cols):
        patch.set_facecolor(col); patch.set_alpha(0.55); patch.set_edgecolor(col)
    rng = np.random.default_rng(11)
    for i, vals in enumerate(data):
        jit = rng.uniform(-0.10, 0.10, len(vals))
        ax.scatter(np.full(len(vals), i)+jit, vals, color=cols[i], s=10,
                   alpha=0.55, edgecolor="none", zorder=2)
    ax.set_xticks(range(len(stages))); ax.set_xticklabels([f"Stage {s}" for s in stages], fontsize=10)
    ax.set_ylabel("panel z (RAI_8)")
    ax.axhline(0, color=MUTED, ls="--", lw=0.7)
    from scipy.stats import kruskal
    valid = [v for v in data if len(v) >= 2]
    h, p = kruskal(*valid) if len(valid) >= 2 else (float("nan"), float("nan"))
    ax.text(0.02, 0.04, f"Kruskal–Wallis H = {h:.1f}, p = {p:.3g}",
            transform=ax.transAxes, fontsize=10, color=INK, fontweight="bold")
    panel_title(ax, "Panel z by AJCC stage"); panel_letter(ax, "b")

    fig.text(0.5, 0.04,
             "Higher AJCC stage trends toward elevated BRAF-like + dark-matter zone fraction and lower panel z.",
             ha="center", fontsize=9.5, color=INK, style="italic")
    fig.savefig(OUT_PNG, dpi=220, facecolor=IVORY)
    fig.savefig(OUT_PDF, facecolor=IVORY)
    plt.close(fig)
    print(f"wrote {OUT_PNG}\nwrote {OUT_PDF}")


if __name__ == "__main__":
    main()
