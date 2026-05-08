#!/usr/bin/env python
"""Wave 5B — assemble final results tables, ablation table, figure, MC uncertainty
sanity, and per-allele LOSO comparison vs Wave 1.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.metrics import roc_auc_score
from scipy.stats import pearsonr, spearmanr

ROOT = Path("/home/seungho/personal/THCA_data_analysis/project/results/p_neo_bayesian_2026_05_09")
W5B = ROOT / "wave5b"

# -------------------------------------------------------------------- 1. Per-tag results table

TAGS = ["full", "A_only", "AB", "AC", "AD"]
all_rows = []
for tag in TAGS:
    f = W5B / f"results_{tag}.tsv"
    if not f.exists():
        print(f"missing {f}"); continue
    r = pd.read_csv(f, sep="\t")
    r["tag"] = tag
    all_rows.append(r)
res = pd.concat(all_rows, ignore_index=True)

# Wave 5B canonical results table (Head A only) for the brief
wave5b_results = res[(res["head"] == "A")].copy()
wave5b_results = wave5b_results[["tag", "head", "testset", "n", "n_pos", "AUROC", "CI_lo95", "CI_hi95"]]
wave5b_results.to_csv(W5B / "wave5b_results.tsv", sep="\t", index=False)
print("=== wave5b_results.tsv ===")
print(wave5b_results.to_string(index=False))

# -------------------------------------------------------------------- 2. Ablation table

ablation_rows = []
for tag in TAGS:
    sub = res[(res["tag"] == tag) & (res["head"] == "A") & (res["testset"] == "ITSNdb_no_overlap")]
    if len(sub) == 1:
        ablation_rows.append(dict(combination=tag, no_overlap_AUROC=float(sub["AUROC"].iloc[0]),
                                  CI_lo95=float(sub["CI_lo95"].iloc[0]),
                                  CI_hi95=float(sub["CI_hi95"].iloc[0])))
abl = pd.DataFrame(ablation_rows)
abl.to_csv(W5B / "ablation_results.tsv", sep="\t", index=False)
print("\n=== ablation_results.tsv ===")
print(abl.to_string(index=False))

# -------------------------------------------------------------------- 3. Wave-1 reference

WAVE1_NO_OVERLAP = 0.4106514084507042  # from wave1/itsndb_bayesian_results.tsv
WAVE1_COMBINED = 0.7840628984275395
WAVE1_IN_MASTER = 0.5638045540796963   # ITSNdb_main only — closer comp uses combined

# Use the wave1 in-domain prediction file to compute wave1 in_master AUROC for fair comparison
w1_pred = pd.read_csv(ROOT / "wave1" / "predictions_itsndb.tsv", sep="\t")
w1_in_master = w1_pred[w1_pred["in_master"] == True]
auc_w1_im = roc_auc_score(w1_in_master["label"].astype(int), w1_in_master["pred_mean"]) if len(w1_in_master) > 5 else float("nan")
w1_combined = roc_auc_score(w1_pred["label"].astype(int), w1_pred["pred_mean"])
w1_no_ov = w1_pred[w1_pred["in_master"] == False]
auc_w1_no = roc_auc_score(w1_no_ov["label"].astype(int), w1_no_ov["pred_mean"])
print(f"\n[wave1] combined={w1_combined:.4f} no_overlap={auc_w1_no:.4f} in_master={auc_w1_im:.4f}")

# -------------------------------------------------------------------- 4. Bar chart figure

fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Panel A: ablation bars (no_overlap AUROC)
ax = axes[0]
order = ["W1\n(single-task)", "A_only", "AB", "AC", "AD", "Full ABCD"]
vals = [auc_w1_no] + [abl[abl["combination"] == t]["no_overlap_AUROC"].iloc[0]
                      for t in ["A_only", "AB", "AC", "AD", "full"]]
los  = [np.nan]   + [abl[abl["combination"] == t]["CI_lo95"].iloc[0]
                      for t in ["A_only", "AB", "AC", "AD", "full"]]
his  = [np.nan]   + [abl[abl["combination"] == t]["CI_hi95"].iloc[0]
                      for t in ["A_only", "AB", "AC", "AD", "full"]]
err_lo = [v - lo for v, lo in zip(vals, los)]
err_hi = [hi - v for v, hi in zip(vals, his)]
errs = [[0 if np.isnan(v) else v for v in err_lo], [0 if np.isnan(v) else v for v in err_hi]]
colors = ["#777"] + ["#4292c6"] + ["#9ecae1"]*3 + ["#08306b"]
bars = ax.bar(order, vals, yerr=errs, capsize=4, color=colors)
ax.axhline(0.5, ls="--", color="grey", alpha=0.5)
ax.set_ylim(0.30, 0.70)
ax.set_ylabel("ITSNdb_no_overlap AUROC (Head A)")
ax.set_title("Wave 5B multi-task ablation (no_overlap)")
for b, v in zip(bars, vals):
    if not np.isnan(v):
        ax.text(b.get_x() + b.get_width()/2, v + 0.01, f"{v:.3f}",
                ha="center", fontsize=9)

# Panel B: full result on 3 strata
ax = axes[1]
strata = ["combined", "no_overlap", "in_master"]
w1 = [w1_combined, auc_w1_no, auc_w1_im]
w5b_full = wave5b_results[wave5b_results["tag"] == "full"].set_index("testset")
w5b = [w5b_full.loc[f"ITSNdb_{s}", "AUROC"] if f"ITSNdb_{s}" in w5b_full.index else np.nan
       for s in strata]
x = np.arange(3)
ax.bar(x - 0.18, w1, 0.36, label="Wave 1 (single-task)", color="#777")
ax.bar(x + 0.18, w5b, 0.36, label="Wave 5B (4-head)", color="#08306b")
ax.set_xticks(x); ax.set_xticklabels(strata)
ax.axhline(0.5, ls="--", color="grey", alpha=0.5)
ax.set_ylabel("AUROC")
ax.set_ylim(0.30, 1.00)
ax.set_title("Wave 1 vs Wave 5B by ITSNdb stratum")
ax.legend()
for xx, v1, v5 in zip(x, w1, w5b):
    ax.text(xx - 0.18, v1 + 0.01, f"{v1:.3f}", ha="center", fontsize=8)
    ax.text(xx + 0.18, v5 + 0.01, f"{v5:.3f}", ha="center", fontsize=8)

plt.tight_layout()
plt.savefig(W5B / "fig_wave5b_multitask.png", dpi=160, bbox_inches="tight")
plt.savefig(W5B / "fig_wave5b_multitask.pdf", bbox_inches="tight")
print(f"\n[fig] wrote fig_wave5b_multitask.png/pdf")

# -------------------------------------------------------------------- 5. Bayesian uncertainty sanity

pred_full = pd.read_csv(W5B / "predictions_full.tsv", sep="\t")
pred_a_only = pd.read_csv(W5B / "predictions_A_only.tsv", sep="\t")

# Mean MC-Dropout std as a proxy for predictive uncertainty
unc_full = pred_full["pred_A_std"].astype(float)
unc_aonly = pred_a_only["pred_A_std"].astype(float)
print(f"\n[unc] mean MC-std Head A — full: {unc_full.mean():.4f}  A_only: {unc_aonly.mean():.4f}")
print(f"[unc] median  MC-std Head A — full: {unc_full.median():.4f}  A_only: {unc_aonly.median():.4f}")

# Relationship of uncertainty to error: rows where prediction is wrong should have higher std
# We define "wrong" as |pred - label| > 0.5
err_full = (pred_full["pred_A_mean"] - pred_full["label_A"]).abs()
rho_full, _ = spearmanr(err_full, unc_full)
err_a = (pred_a_only["pred_A_mean"] - pred_a_only["label_A"]).abs()
rho_a, _ = spearmanr(err_a, unc_aonly)
print(f"[unc] spearman(|err|, MC-std) — full: {rho_full:.3f}  A_only: {rho_a:.3f}")

unc_summary = pd.DataFrame([
    dict(model="A_only",   mean_std=float(unc_aonly.mean()), median_std=float(unc_aonly.median()),
         spearman_err_std=float(rho_a)),
    dict(model="full4head", mean_std=float(unc_full.mean()),  median_std=float(unc_full.median()),
         spearman_err_std=float(rho_full)),
])
unc_summary.to_csv(W5B / "uncertainty_compare.tsv", sep="\t", index=False)

# -------------------------------------------------------------------- 6. Per-head MSE on no_overlap (sanity)

for tag in ["full", "AB", "AC", "AD"]:
    pf = W5B / f"predictions_{tag}.tsv"
    if not pf.exists(): continue
    p = pd.read_csv(pf, sep="\t")
    m = (p["in_master"] == False).values
    print(f"\n[sanity {tag}] MSE on no_overlap (z-scaled):")
    for h, gtcol in [("B", "label_B_z"), ("C", "label_C_z"), ("D", "label_D_z")]:
        col = f"pred_{h}_mean"
        if col in p.columns:
            mse = float(((p.loc[m, col] - p.loc[m, gtcol].astype(float)) ** 2).mean())
            print(f"   head {h}  MSE={mse:.4f}")

print("\n[done] summarize complete")
