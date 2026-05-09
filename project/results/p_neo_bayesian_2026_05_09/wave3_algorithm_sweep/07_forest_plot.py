#!/usr/bin/env python3
"""Forest plot: leakage-stratified AUROC across off-the-shelf + ours.

X-axis: AUROC [0, 1]
Y-axis: algorithms (off-the-shelf + ours)
Two columns: in_master=True (light) vs in_master=False (dark, headline)
Annotation: Δ inflation per row.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

ROOT = Path(__file__).resolve().parent
res = pd.read_csv(ROOT / "algorithm_sweep_results.tsv", sep="\t")
print(res.head())

COMB = "ext_itsndb_main+ext_itsndb_val"
sub = res[res["testset"] == COMB].copy()

# Get rows for in_master True/False per algorithm
piv = sub.pivot_table(
    index="algorithm",
    columns="in_master",
    values=["AUROC", "AUROC_lo95", "AUROC_hi95", "n"],
    aggfunc="first",
)
piv.columns = ["_".join(map(str, c)) for c in piv.columns]
print(piv)

# Display order: off-the-shelf first, then ours
DISPLAY_ORDER = [
    "MHCflurry_presentation",
    "MHCflurry_affinity_neg",
    "BigMHC_EL",
    "BigMHC_IM",
    "PRIME_score",
    "PRIME_rank_neg",
    "Ours_ESM2_Bayesian",
    "Ours_BiophysStructLR",
    "Ours_Wave1_plus_Struct",
    "Ours_VQC",
    "Ours_LR8d",
]
algos = [a for a in DISPLAY_ORDER if a in piv.index]
piv = piv.loc[algos]

# Compute delta inflation
piv["delta"] = piv["AUROC_True"] - piv["AUROC_False"]

# Plot
fig, ax = plt.subplots(figsize=(11, 6.5))

y = np.arange(len(algos))
height = 0.35

# in_master=True (light, training-overlapping = inflated)
auroc_T = piv["AUROC_True"].values
lo_T = piv["AUROC_lo95_True"].values
hi_T = piv["AUROC_hi95_True"].values
err_T = np.array([auroc_T - lo_T, hi_T - auroc_T])

# in_master=False (dark, no-overlap = headline truth)
auroc_F = piv["AUROC_False"].values
lo_F = piv["AUROC_lo95_False"].values
hi_F = piv["AUROC_hi95_False"].values
err_F = np.array([auroc_F - lo_F, hi_F - auroc_F])

ax.errorbar(auroc_T, y - height/2, xerr=err_T, fmt="s",
            color="#88B0D6", ecolor="#88B0D6", capsize=3, markersize=8,
            label="in_master=True (training-overlapping, n≈208)")
ax.errorbar(auroc_F, y + height/2, xerr=err_F, fmt="o",
            color="#1f4e79", ecolor="#1f4e79", capsize=3, markersize=9,
            label="in_master=False (no_overlap, n≈103)")

# Chance line
ax.axvline(0.5, color="grey", linestyle="--", linewidth=1, alpha=0.7, label="Chance")

# Annotation: delta
for i, alg in enumerate(algos):
    d = piv.loc[alg, "delta"]
    if not np.isnan(d):
        ax.text(1.01, i, f"Δ={d:+.2f}", va="center", fontsize=9,
                color="firebrick" if d > 0.15 else ("forestgreen" if d < -0.05 else "black"))

ax.set_yticks(y)
ax.set_yticklabels(algos)
ax.invert_yaxis()
ax.set_xlim(0.20, 1.05)
ax.set_xlabel("AUROC (95% bootstrap CI, 1000 resamples)")
ax.set_title("Wave 3 — Leakage Inflation in Neoantigen Immunogenicity Predictors\n(ITSNdb combined; Δ = in_master AUROC − no_overlap AUROC)")
ax.grid(axis="x", alpha=0.3)
ax.legend(loc="lower right", framealpha=0.9, fontsize=9)

plt.tight_layout()
fig.savefig(ROOT / "fig_forest_leakage_inflation.png", dpi=200, bbox_inches="tight")
fig.savefig(ROOT / "fig_forest_leakage_inflation.pdf", bbox_inches="tight")
print(f"saved fig_forest_leakage_inflation.png/.pdf")

# Also save the pivot for reporting
piv.to_csv(ROOT / "forest_pivot.tsv", sep="\t")
print(piv[["AUROC_True", "AUROC_False", "delta", "n_True", "n_False"]])
