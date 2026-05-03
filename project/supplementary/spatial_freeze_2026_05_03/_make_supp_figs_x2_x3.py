#!/usr/bin/env python3
"""Repackage existing meta tables into standalone Supp Fig X2 (cross-validation scatter)
and Supp Fig X3 (autoimmune negative control). Pure plotting, no new analysis."""
from __future__ import annotations
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, pearsonr

OUT = Path("project/supplementary/spatial_freeze_2026_05_03")
SUMMARY = Path("project_external_st/results/meta/sample_level_score_summary.tsv")

CONDITION_ORDER = ["CONTROL", "HT", "GD", "PTC_HT"]
PALETTE = {"CONTROL": "#4C72B0", "HT": "#DD8452", "GD": "#55A467", "PTC_HT": "#C44E52"}

s = pd.read_csv(SUMMARY, sep="\t")

# --- Supp Fig X2: DM1_like vs THYROID_NONOVERLAP cross-validation, all 12 external samples ---
fig, ax = plt.subplots(figsize=(8, 6.5))
for ds, sub in s.groupby("dataset"):
    ax.scatter(sub["mean_THYROID_NONOVERLAP_score_resid_epi25"],
               sub["mean_DM1_like_score_resid_epi25"],
               c=[PALETTE[c] for c in sub["condition"]],
               s=200, edgecolor="black", linewidths=1.2,
               marker="o" if ds == "GSE230424" else "s", label=ds)
x = s["mean_THYROID_NONOVERLAP_score_resid_epi25"].to_numpy()
y = s["mean_DM1_like_score_resid_epi25"].to_numpy()
ok = ~(np.isnan(x) | np.isnan(y))
rho, p = spearmanr(x[ok], y[ok])
r, p2 = pearsonr(x[ok], y[ok])
xs = np.sort(x[ok])
ax.plot(xs, np.poly1d(np.polyfit(x[ok], y[ok], 1))(xs),
        color="black", lw=1.2, ls="--", alpha=0.6, label="OLS fit")
for cond, color in PALETTE.items():
    if cond in s["condition"].values:
        ax.scatter([], [], c=color, s=140, label=cond, edgecolor="black")
ax.set_xlabel("THYROID_NONOVERLAP score (depth-residualized, epi top 25% sample mean)\n"
              "[SLC26A4, IYD, DUOX1, DUOX2, TFF3, HHEX, GLIS3, DIO2 — zero overlap with RAI_8]",
              fontsize=10)
ax.set_ylabel("DM1_like score (= -RAI_8; depth-residualized, epi top 25% sample mean)", fontsize=10)
ax.set_title(f"Supp Fig X2 — Cross-validation against independent thyroid-lineage axis\n"
             f"All 12 external samples · Spearman ρ = {rho:.2f} (p = {p:.1e}) · Pearson r = {r:.2f} (p = {p2:.1e})",
             fontsize=11)
ax.axhline(0, color="grey", lw=0.4, ls=":")
ax.axvline(0, color="grey", lw=0.4, ls=":")
ax.legend(loc="upper right", frameon=False, fontsize=9)
fig.tight_layout()
fig.savefig(OUT / "SuppFig_X2_DM1_thyroid_nonoverlap_crossvalidation.png", dpi=200, bbox_inches="tight")
fig.savefig(OUT / "SuppFig_X2_DM1_thyroid_nonoverlap_crossvalidation.pdf", bbox_inches="tight")
plt.close(fig)

# --- Supp Fig X3: GSE248205 autoimmune negative control (CONTROL/HT/GD DM1 distribution) ---
gse248 = s[s["dataset"] == "GSE248205"].copy()
fig, ax = plt.subplots(figsize=(7.5, 6))
order248 = ["CONTROL", "HT", "GD"]
positions = np.arange(len(order248))
for i, cond in enumerate(order248):
    sub = gse248[gse248["condition"] == cond]
    ax.scatter([i] * len(sub), sub["mean_DM1_like_score_resid_epi25"],
               c=PALETTE[cond], s=200, edgecolor="black", linewidths=1.2,
               zorder=3, label=cond if i == 0 else None)
    if len(sub) > 0:
        m = sub["mean_DM1_like_score_resid_epi25"].mean()
        sem = sub["mean_DM1_like_score_resid_epi25"].std(ddof=1) / np.sqrt(len(sub)) if len(sub) > 1 else 0
        ax.errorbar(i, m, yerr=sem, color="black", capsize=10, lw=2, zorder=4)
        ax.scatter(i, m, marker="_", s=400, color="black", lw=3, zorder=5)
        ax.text(i, m + (sem if sem > 0 else 0.02) + 0.03, f"n={len(sub)}\nmean={m:.2f}",
                ha="center", fontsize=10)
ax.axhline(0, color="grey", lw=0.5, ls=":")
ax.set_xticks(positions); ax.set_xticklabels(order248, fontsize=12)
ax.set_xlabel("")
ax.set_ylabel("DM1_like score (depth-residualized, epi top 25% sample mean)", fontsize=10)
ax.set_title("Supp Fig X3 — GSE248205 autoimmune thyroid (no cancer) control\n"
             "DM1_like is direction-consistent but underpowered by inflammation alone\n"
             "(HT n=3, GD n=3, CONTROL n=2; Mann-Whitney p > 0.10 for all pairs)",
             fontsize=11)
fig.tight_layout()
fig.savefig(OUT / "SuppFig_X3_GSE248205_autoimmune_negative_control.png", dpi=200, bbox_inches="tight")
fig.savefig(OUT / "SuppFig_X3_GSE248205_autoimmune_negative_control.pdf", bbox_inches="tight")
plt.close(fig)

print("X2 →", OUT / "SuppFig_X2_DM1_thyroid_nonoverlap_crossvalidation.png")
print("X3 →", OUT / "SuppFig_X3_GSE248205_autoimmune_negative_control.png")
