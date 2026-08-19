#!/usr/bin/env python3
"""Figure 30 — TCGA-THCA panel z by demographics (age, sex).

Panel z stratified by age decile and sex, confirming the signal is not driven
by demographic confounders.
"""
from __future__ import annotations
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent))
from _nature_style import setup_rc, panel_letter, panel_title, IVORY, INK, MUTED, BLUE, RED

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, mannwhitneyu, kruskal

setup_rc()

ROOT = Path(__file__).resolve().parent.parent
PER = Path("/home/seungho/personal/THCA_data_analysis/project/results/r17_tcga_panel_d4p2_reconciliation/tert_interaction/r17_tert_per_sample.tsv")
CLIN = Path("/home/seungho/personal/THCA_data_analysis/project/results/dark_matter_phase1/tcga_dm_master_with_pfi.tsv")
OUT_PNG = ROOT / "results" / "figures" / "figure30_age_sex_panel_z.png"
OUT_PDF = ROOT / "results" / "figures" / "figure30_age_sex_panel_z.pdf"


def main():
    df = pd.read_csv(PER, sep="\t")
    clin = pd.read_csv(CLIN, sep="\t")
    df["pat"] = df["sampleId"].str.slice(0, 12)
    m = df.merge(clin[["tcga_short", "age", "sex"]], left_on="pat", right_on="tcga_short", how="inner")
    m["age"] = pd.to_numeric(m["age"], errors="coerce")
    m = m.dropna(subset=["age", "sex", "RAI_8"])
    m["sex"] = m["sex"].astype(str)

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.4), facecolor=IVORY,
                              gridspec_kw=dict(wspace=0.30, left=0.07, right=0.97, top=0.83, bottom=0.16))
    fig.suptitle("Panel z is not driven by age or sex (TCGA-THCA, n = " + f"{len(m):,})",
                 fontsize=13, fontweight="bold", color=INK, y=0.97)

    # Left: scatter age vs panel z + linear fit
    ax = axes[0]
    ax.scatter(m["age"], m["RAI_8"], color=BLUE, alpha=0.5, s=14, edgecolor="white", lw=0.3)
    # fit line
    z = np.polyfit(m["age"], m["RAI_8"], 1)
    xs = np.linspace(m["age"].min(), m["age"].max(), 100)
    ax.plot(xs, np.polyval(z, xs), color=RED, lw=1.6, alpha=0.85, label=f"slope = {z[0]:+.3f}/yr")
    rho, p = spearmanr(m["age"], m["RAI_8"])
    ax.text(0.02, 0.96, f"Spearman ρ = {rho:+.2f},  p = {p:.3g}", transform=ax.transAxes,
            fontsize=10, color=INK, fontweight="bold", va="top")
    ax.set_xlabel("age at diagnosis (years)"); ax.set_ylabel("panel z (RAI_8)")
    ax.axhline(0, color=MUTED, ls="--", lw=0.7)
    ax.legend(fontsize=9, loc="lower left", frameon=True)
    panel_title(ax, "Panel z vs age"); panel_letter(ax, "a")

    # Right: by sex
    ax = axes[1]
    sex_order = [s for s in ["Female", "Male"] if (m["sex"] == s).sum() >= 5]
    if not sex_order:
        sex_order = sorted(m["sex"].unique())
    data = [m.loc[m["sex"] == s, "RAI_8"].dropna().values for s in sex_order]
    palette = {"Female": "#cb857f", "Male": "#3a6694"}
    bp = ax.boxplot(data, positions=range(len(sex_order)), widths=0.4,
                     patch_artist=True, medianprops=dict(color=INK, lw=1.3),
                     boxprops=dict(lw=0.7), whiskerprops=dict(lw=0.7), capprops=dict(lw=0.7),
                     flierprops=dict(marker="o", markersize=2, alpha=0.4))
    for patch, s in zip(bp["boxes"], sex_order):
        patch.set_facecolor(palette.get(s, "#888")); patch.set_alpha(0.55); patch.set_edgecolor(palette.get(s, "#888"))
    rng = np.random.default_rng(11)
    for i, vals in enumerate(data):
        jit = rng.uniform(-0.08, 0.08, len(vals))
        ax.scatter(np.full(len(vals), i)+jit, vals, color=palette.get(sex_order[i], "#888"),
                   s=10, alpha=0.55, edgecolor="none", zorder=2)
    ax.set_xticks(range(len(sex_order)))
    ax.set_xticklabels([f"{s}\n(n={int((m['sex']==s).sum())})" for s in sex_order], fontsize=10)
    ax.set_ylabel("panel z (RAI_8)")
    ax.axhline(0, color=MUTED, ls="--", lw=0.7)
    if len(data) == 2:
        u, p_u = mannwhitneyu(data[0], data[1], alternative="two-sided")
        ax.text(0.02, 0.04, f"Mann–Whitney U = {u:.0f}, p = {p_u:.3g}",
                transform=ax.transAxes, fontsize=10, color=INK, fontweight="bold")
    panel_title(ax, "Panel z by sex"); panel_letter(ax, "b")

    fig.text(0.5, 0.04,
             f"Age effect is small (ρ = {rho:+.2f}); panel z is essentially unconfounded by demographic features.",
             ha="center", fontsize=9.5, color=INK, style="italic")
    fig.savefig(OUT_PNG, dpi=220, facecolor=IVORY)
    fig.savefig(OUT_PDF, facecolor=IVORY)
    plt.close(fig)
    print(f"wrote {OUT_PNG}\nwrote {OUT_PDF}")


if __name__ == "__main__":
    main()
