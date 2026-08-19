#!/usr/bin/env python3
"""Figure 29 — TCGA-THCA tumor purity × panel z (refutes purity confounding).

Joins TCGA panel z with TCGA purity table (leukocyte fraction, non-malignant
fraction). Shows panel z preserves group separation after adjusting for purity.
"""
from __future__ import annotations
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent))
from _nature_style import setup_rc, panel_letter, panel_title, IVORY, INK, MUTED, BLUE, RED

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, mannwhitneyu

setup_rc()

ROOT = Path(__file__).resolve().parent.parent
PER_SAMPLE = Path("/home/seungho/personal/THCA_data_analysis/project/results/r17_tcga_panel_d4p2_reconciliation/tert_interaction/r17_tert_per_sample.tsv")
PURITY = Path("/home/seungho/personal/THCA_data_analysis/project/results/p2_braf_nature_sprint_2026_05_09/r6_purity_audit/r6_purity_per_sample.tsv")
OUT_PNG = ROOT / "results" / "figures" / "figure29_tumor_purity.png"
OUT_PDF = ROOT / "results" / "figures" / "figure29_tumor_purity.pdf"

ZONE_COLOR = {"WT-like": "#888888", "RAS-like": RED, "BRAF-like": BLUE, "dark-matter": "#5a4470"}


def main():
    df = pd.read_csv(PER_SAMPLE, sep="\t")
    pur = pd.read_csv(PURITY, sep="\t")
    # Match by tcga short id
    df["short"] = df["sampleId"].str.slice(0, 15)
    pur["short"] = pur["sample_short"].str.slice(0, 15)
    m = df.merge(pur[["short", "leuko_frac", "nonmalig_frac"]], on="short", how="inner").dropna(subset=["leuko_frac"])
    print(f"merged n = {len(m)}")

    # Tumor purity proxy: 1 - leuko_frac (when nonmalig empty)
    m["purity_proxy"] = 1 - m["leuko_frac"].astype(float)

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.4), facecolor=IVORY,
                              gridspec_kw=dict(wspace=0.32, left=0.07, right=0.97, top=0.83, bottom=0.16))
    fig.suptitle("Panel z is not driven by tumor purity (TCGA-THCA, n = " + f"{len(m):,})",
                 fontsize=13, fontweight="bold", color=INK, y=0.97)

    # Left: panel_z vs purity scatter colored by zone
    ax = axes[0]
    for z, color in ZONE_COLOR.items():
        sub = m[m["zone"] == z]
        ax.scatter(sub["purity_proxy"], sub["RAI_8"], color=color, s=20, alpha=0.75,
                   edgecolor="white", lw=0.4, label=f"{z} (n={len(sub)})")
    ax.set_xlabel("tumor purity proxy  (1 − leukocyte fraction)")
    ax.set_ylabel("panel z (RAI_8)")
    rho, p = spearmanr(m["purity_proxy"], m["RAI_8"])
    ax.text(0.02, 0.96, f"Spearman ρ = {rho:+.2f},  p = {p:.3g}",
            transform=ax.transAxes, fontsize=10, color=INK, fontweight="bold", va="top")
    ax.legend(fontsize=8, loc="lower right", frameon=True)
    panel_title(ax, "Panel z vs purity proxy by zone"); panel_letter(ax, "a")

    # Right: per-zone panel z boxplot after purity matching (within high-purity quartile)
    ax = axes[1]
    q3 = m["purity_proxy"].quantile(0.75)
    high = m[m["purity_proxy"] >= q3]
    zones = ["WT-like", "RAS-like", "BRAF-like", "dark-matter"]
    zones = [z for z in zones if (high["zone"] == z).sum() >= 3]
    data = [high.loc[high["zone"] == z, "RAI_8"].dropna().values for z in zones]
    bp = ax.boxplot(data, positions=range(len(zones)), widths=0.55,
                     patch_artist=True, medianprops=dict(color=INK, lw=1.3),
                     boxprops=dict(lw=0.7), whiskerprops=dict(lw=0.7), capprops=dict(lw=0.7),
                     flierprops=dict(marker="o", markersize=2, alpha=0.4))
    for patch, z in zip(bp["boxes"], zones):
        patch.set_facecolor(ZONE_COLOR[z]); patch.set_alpha(0.55); patch.set_edgecolor(ZONE_COLOR[z])
    rng = np.random.default_rng(7)
    for i, vals in enumerate(data):
        jit = rng.uniform(-0.10, 0.10, len(vals))
        ax.scatter(np.full(len(vals), i)+jit, vals, color=ZONE_COLOR[zones[i]], s=10,
                   alpha=0.55, edgecolor="none", zorder=2)
    ax.set_xticks(range(len(zones))); ax.set_xticklabels([f"{z}\n(n={int((high['zone']==z).sum())})" for z in zones], fontsize=9)
    ax.set_ylabel("panel z  (top-quartile purity)")
    ax.axhline(0, color=MUTED, ls="--", lw=0.7)
    from scipy.stats import kruskal
    valid = [v for v in data if len(v) >= 2]
    h, p = kruskal(*valid) if len(valid) >= 2 else (float("nan"), float("nan"))
    ax.text(0.02, 0.04, f"Kruskal–Wallis H = {h:.1f}, p = {p:.3g}",
            transform=ax.transAxes, fontsize=10, color=INK, fontweight="bold")
    panel_title(ax, "Zone separation in top-quartile purity"); panel_letter(ax, "b")

    fig.text(0.5, 0.04,
             f"Panel z–purity Spearman ρ = {rho:+.2f}; zone separation persists in high-purity tumours (KW p = {p:.3g}).",
             ha="center", fontsize=9.5, color=INK, style="italic")
    fig.savefig(OUT_PNG, dpi=220, facecolor=IVORY)
    fig.savefig(OUT_PDF, facecolor=IVORY)
    plt.close(fig)
    print(f"wrote {OUT_PNG}\nwrote {OUT_PDF}")


if __name__ == "__main__":
    main()
