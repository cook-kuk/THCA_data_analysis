"""Regenerate hero forest + inflation-gap figs after appending Wave1."""
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

WAVE7 = Path("/home/seungho/personal/THCA_data_analysis/project/results/p_neo_bayesian_2026_05_09/wave7")
df = pd.read_csv(WAVE7 / "wave7_results.tsv", sep="\t")

# Hero forest
no_ov = df[df["subset"] == "ITSNdb_no_overlap"].copy().sort_values("AUROC", ascending=True).reset_index(drop=True)
plt.figure(figsize=(10, max(6, len(no_ov) * 0.32)))
y_pos = np.arange(len(no_ov))
aucs = no_ov["AUROC"].values
los = np.where(no_ov["AUROC_lo95"].isna(), aucs, no_ov["AUROC_lo95"].values)
his = np.where(no_ov["AUROC_hi95"].isna(), aucs, no_ov["AUROC_hi95"].values)
errs_lo = aucs - los
errs_hi = his - aucs
colors = []
for m in no_ov["model"]:
    if "W7-" in m: colors.append("#d62728")
    elif "MHCflurry" in m: colors.append("#2ca02c")
    elif "Wave1" in m: colors.append("#9467bd")
    elif "Structure_LR" in m: colors.append("#ff7f0e")
    else: colors.append("#7f7f7f")
plt.errorbar(aucs, y_pos, xerr=[errs_lo, errs_hi], fmt='none', ecolor='gray', capsize=3, alpha=0.6)
for i, c in enumerate(colors):
    plt.plot(aucs[i], y_pos[i], 'o', color=c, markersize=9, mec='black')
plt.axvline(0.5, color="gray", linestyle="--", alpha=0.5, label="chance (0.5)")
plt.axvline(0.668, color="#2ca02c", linestyle="--", alpha=0.5, label="MHCflurry alone (0.668)")
plt.yticks(y_pos, no_ov["model"], fontsize=8)
plt.xlabel("AUROC (ITSNdb no_overlap; bars = 95% bootstrap CI)")
plt.title(f"Wave 7 Hero Forest — {len(no_ov)} methods on ITSNdb no_overlap (n=106)")
plt.legend(loc="lower right", fontsize=8)
plt.xlim(0.25, 0.95)
plt.tight_layout()
plt.savefig(WAVE7 / "fig_wave7_hero_forest.png", dpi=150)
plt.savefig(WAVE7 / "fig_wave7_hero_forest.pdf")
plt.close()
print("hero forest regenerated")

# Inflation gap with Wave1 included
piv = df.pivot_table(index="model", columns="subset", values="AUROC", aggfunc="first")
piv["gap"] = piv["ITSNdb_in_master"] - piv["ITSNdb_no_overlap"]
piv2 = piv.dropna(subset=["gap"]).sort_values("gap", ascending=True)
plt.figure(figsize=(9, max(4, len(piv2) * 0.32)))
colors = []
for m in piv2.index:
    if "W7-" in m: colors.append("#d62728")
    elif "MHCflurry" in m: colors.append("#2ca02c")
    elif "Wave1" in m: colors.append("#9467bd")
    elif "Structure_LR" in m: colors.append("#ff7f0e")
    else: colors.append("#7f7f7f")
plt.barh(np.arange(len(piv2)), piv2["gap"], color=colors, edgecolor='black', alpha=0.8)
plt.yticks(np.arange(len(piv2)), piv2.index, fontsize=8)
plt.xlabel("Δ AUROC (in_master − no_overlap)  →  larger = more inflated/leaky")
plt.title("Wave 7 — Inflation gap across 17+ methods (truly held-out vs leakage-contaminated)")
plt.axvline(0, color="black", alpha=0.5)
plt.tight_layout()
plt.savefig(WAVE7 / "fig_wave7_inflation_gap.png", dpi=150)
plt.savefig(WAVE7 / "fig_wave7_inflation_gap.pdf")
plt.close()
print("inflation gap regenerated")
