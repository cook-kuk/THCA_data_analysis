"""Build all Wave 6 figures.

Generated:
  fig_track6a_uplift.png/pdf            - Bar chart MHCflurry-only / +PWM / +ESM2 / triple combo on no_overlap
  fig_track6b_representation_hierarchy  - DRP Table 1 forest
  fig_track6c_shortcut_panel            - 6-panel shortcut tests
  fig_track6c_counterfactual_examples   - example sequence edits
"""
from __future__ import annotations
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path("/home/seungho/personal/THCA_data_analysis/project/results/p_neo_bayesian_2026_05_09")
WAVE6 = ROOT / "wave6"

MHCFLURRY_BASELINE = 0.6683

# ---- Track 6A uplift ----
results = pd.read_csv(WAVE6 / "track6b_representation_comparison.tsv", sep="\t")
score_cal = pd.read_csv(WAVE6 / "track6a_score_calibration.tsv", sep="\t")
print("results:\n", results)

# Build uplift comparison: MHCflurry-score (0.668) vs (mhcflurry_only, +pwm, +esm2, +pwm+esm2)
key_configs = ["mhcflurry_only","mhcflurry+pwm","mhcflurry+esm2","mhcflurry+pwm+esm2"]
plot_rows = []
plot_rows.append({"name":"MHCflurry score\n(baseline)","auc":MHCFLURRY_BASELINE,"lo":0.536,"hi":0.787,"color":"#444444"})
for cfg in key_configs:
    sub = results[results["config"]==cfg]
    if len(sub):
        r = sub.iloc[0]
        plot_rows.append({"name":cfg.replace("_"," ").replace("+","\n+"),
                          "auc":r["no_overlap_auc_mean"],
                          "lo":r["no_overlap_auc_lo95"],
                          "hi":r["no_overlap_auc_hi95"],
                          "color":"#1f77b4"})

fig, ax = plt.subplots(figsize=(8,5))
xs = np.arange(len(plot_rows))
aucs = [r["auc"] for r in plot_rows]
los  = [r["auc"]-r["lo"] for r in plot_rows]
his  = [r["hi"]-r["auc"] for r in plot_rows]
colors = [r["color"] for r in plot_rows]
ax.bar(xs, aucs, color=colors, edgecolor="k", linewidth=0.8)
ax.errorbar(xs, aucs, yerr=[los, his], fmt="none", color="k", capsize=4, lw=1)
ax.axhline(0.5, color="gray", linestyle=":", lw=1)
ax.axhline(MHCFLURRY_BASELINE, color="#444444", linestyle="--", lw=1, alpha=0.6)
ax.set_xticks(xs); ax.set_xticklabels([r["name"] for r in plot_rows], fontsize=8)
ax.set_ylabel("ITSNdb_no_overlap AUROC")
ax.set_ylim(0.40, 0.85)
ax.set_title("Track 6A: MHCflurry-features-as-backbone vs published score baseline")
for i,a in enumerate(aucs):
    ax.text(i, a+0.015, f"{a:.3f}", ha="center", fontsize=8)
plt.tight_layout()
plt.savefig(WAVE6/"fig_track6a_uplift.png", dpi=150)
plt.savefig(WAVE6/"fig_track6a_uplift.pdf")
plt.close()
print("saved fig_track6a_uplift")

# ---- Track 6B representation hierarchy ----
ranking = results.sort_values("no_overlap_auc_mean", ascending=True).reset_index(drop=True)
fig, ax = plt.subplots(figsize=(9,6))
y = np.arange(len(ranking))
ax.errorbar(ranking["no_overlap_auc_mean"], y,
            xerr=[ranking["no_overlap_auc_mean"]-ranking["no_overlap_auc_lo95"],
                  ranking["no_overlap_auc_hi95"]-ranking["no_overlap_auc_mean"]],
            fmt="o", color="#1f77b4", capsize=4, ms=8)
ax.axvline(0.5, color="gray", linestyle=":", lw=1)
ax.axvline(MHCFLURRY_BASELINE, color="red", linestyle="--", lw=1.2, label=f"MHCflurry score (0.668)")
ax.set_yticks(y); ax.set_yticklabels(ranking["config"])
ax.set_xlabel("ITSNdb_no_overlap AUROC (mean over seeds)")
ax.set_xlim(0.30, 0.85)
ax.set_title("Track 6B: Representation hierarchy under identical Bayesian-MLP head")
ax.legend(loc="lower right")
for i, r in ranking.iterrows():
    ax.text(r["no_overlap_auc_mean"]+0.005, i, f"  {r['no_overlap_auc_mean']:.3f}", va="center", fontsize=7)
plt.tight_layout()
plt.savefig(WAVE6/"fig_track6b_representation_hierarchy.png", dpi=150)
plt.savefig(WAVE6/"fig_track6b_representation_hierarchy.pdf")
plt.close()
print("saved fig_track6b_hierarchy")

# ---- Track 6C shortcut panel (6 sub-plots) ----
abl   = pd.read_csv(WAVE6 / "track6c_feature_ablation.tsv", sep="\t")
hla_d = pd.read_csv(WAVE6 / "track6c_per_hla.tsv", sep="\t")
len_d = pd.read_csv(WAVE6 / "track6c_per_length.tsv", sep="\t")
pos_d = pd.read_csv(WAVE6 / "track6c_position_ablation.tsv", sep="\t")
sup_d = pd.read_csv(WAVE6 / "track6c_per_supertype.tsv", sep="\t")
src_d = pd.read_csv(WAVE6 / "track6c_per_source.tsv", sep="\t")

fig, axes = plt.subplots(2, 3, figsize=(15, 9))

# (a) feature ablation
ax = axes[0,0]
xs = np.arange(len(abl))
ax.bar(xs, abl["no_overlap_auc"], color=["#1f77b4"]+["#aaa"]*(len(abl)-1), edgecolor="k")
ax.axhline(0.5, color="gray", linestyle=":")
ax.set_xticks(xs); ax.set_xticklabels(abl["test"], rotation=20, ha="right", fontsize=7)
ax.set_ylabel("AUROC"); ax.set_ylim(0.3, 0.85)
ax.set_title("(a) Feature ablation (sanity)")
for i, a in enumerate(abl["no_overlap_auc"]):
    ax.text(i, a+0.01, f"{a:.3f}", ha="center", fontsize=7)

# (b) per HLA
ax = axes[0,1]
hla_top = hla_d.head(15)
ys = np.arange(len(hla_top))[::-1]
ax.errorbar(hla_top["auc"], ys,
            xerr=[hla_top["auc"]-hla_top["lo"], hla_top["hi"]-hla_top["auc"]],
            fmt="o", color="#2ca02c", capsize=3, ms=5)
ax.axvline(0.5, color="gray", linestyle=":"); ax.axvline(MHCFLURRY_BASELINE, color="red", linestyle="--", lw=1)
ax.set_yticks(ys); ax.set_yticklabels(hla_top["hla"], fontsize=7)
ax.set_xlabel("AUROC"); ax.set_xlim(0.0, 1.05)
ax.set_title(f"(b) Per-HLA (n={len(hla_d)} alleles ≥5 samples, top 15)")

# (c) per length
ax = axes[0,2]
xs = np.arange(len(len_d))
ax.bar(xs, len_d["auc"], color="#ff7f0e", edgecolor="k")
ax.errorbar(xs, len_d["auc"], yerr=[len_d["auc"]-len_d["lo"], len_d["hi"]-len_d["auc"]],
            fmt="none", color="k", capsize=3)
ax.axhline(0.5, color="gray", linestyle=":")
ax.axhline(MHCFLURRY_BASELINE, color="red", linestyle="--", lw=1)
ax.set_xticks(xs); ax.set_xticklabels([f"{int(L)}-mer\n(n={int(n)})" for L,n in zip(len_d["length"], len_d["n"])], fontsize=7)
ax.set_ylabel("AUROC"); ax.set_ylim(0.0, 1.05)
ax.set_title("(c) Per-length stratification")

# (d) position ablation (9-mer surrogate)
ax = axes[1,0]
pos_d2 = pos_d[pos_d["position"]!="baseline_full"]
xs = np.arange(len(pos_d2))
deltas = pos_d2["delta"].values
colors = ["#d62728" if d < 0 else "#2ca02c" for d in deltas]
ax.bar(xs, deltas, color=colors, edgecolor="k")
ax.axhline(0, color="black", lw=0.8)
ax.set_xticks(xs); ax.set_xticklabels(pos_d2["position"])
ax.set_ylabel("Δ AUROC vs baseline\n(zero this position)")
ax.set_title("(d) Position ablation (1-hot 9-mer surrogate)")
for i, d in enumerate(deltas):
    ax.text(i, d-(0.01 if d<0 else -0.005), f"{d:+.3f}", ha="center", fontsize=7,
            va="top" if d<0 else "bottom")

# (e) supertype
ax = axes[1,1]
ys = np.arange(len(sup_d))[::-1]
ax.errorbar(sup_d["auc"], ys, xerr=[sup_d["auc"]-sup_d["lo"], sup_d["hi"]-sup_d["auc"]],
            fmt="o", color="#9467bd", capsize=3)
ax.axvline(0.5, color="gray", linestyle=":")
ax.axvline(MHCFLURRY_BASELINE, color="red", linestyle="--", lw=1)
ax.set_yticks(ys); ax.set_yticklabels(sup_d["supertype"], fontsize=7)
ax.set_xlabel("AUROC"); ax.set_xlim(0,1.05)
ax.set_title("(e) Per HLA-supertype")

# (f) source
ax = axes[1,2]
ys = np.arange(len(src_d))[::-1]
ax.errorbar(src_d["auc"], ys, xerr=[src_d["auc"]-src_d["lo"], src_d["hi"]-src_d["auc"]],
            fmt="o", color="#8c564b", capsize=3)
ax.axvline(0.5, color="gray", linestyle=":")
ax.axvline(MHCFLURRY_BASELINE, color="red", linestyle="--", lw=1)
ax.set_yticks(ys); ax.set_yticklabels(src_d["source"], fontsize=7)
ax.set_xlabel("AUROC"); ax.set_xlim(0,1.05)
ax.set_title("(f) Per source")

plt.tight_layout()
plt.savefig(WAVE6/"fig_track6c_shortcut_panel.png", dpi=150)
plt.savefig(WAVE6/"fig_track6c_shortcut_panel.pdf")
plt.close()
print("saved fig_track6c_shortcut_panel")

# ---- Counterfactual examples ----
edits = pd.read_csv(WAVE6 / "track6c_counterfactual_edits.tsv", sep="\t")
top_edits = edits.sort_values("delta", ascending=False).head(8)
fig, ax = plt.subplots(figsize=(11,5))
ax.axis("off")
text_lines = ["Top counterfactual edits (single-AA flip raising prediction)\n"]
text_lines.append(f"  {'orig':12s}  {'pos':>3s}  {'aa':>3s}  {'base→new':>15s}  {'Δ':>7s}  {'flip':>5s}  {'hla':<14s}")
for _, r in top_edits.iterrows():
    if r["best_pos"] > 0:
        new_aa = r["best_aa"]
        pos_idx = int(r["best_pos_idx"])
        new_pep = r["orig_peptide"][:pos_idx] + new_aa + r["orig_peptide"][pos_idx+1:]
        text_lines.append(f"  {r['orig_peptide']:12s} → {new_pep:12s}  P{int(r['best_pos'])}  {new_aa:>3s}  "
                          f"{r['base_score']:.3f}→{r['new_score']:.3f}  {r['delta']:+.3f}  {str(r['flipped']):>5s}  {r['hla']}")
ax.text(0.01, 0.98, "\n".join(text_lines), family="monospace", fontsize=10, va="top")
plt.tight_layout()
plt.savefig(WAVE6/"fig_track6c_counterfactual_examples.png", dpi=150)
plt.savefig(WAVE6/"fig_track6c_counterfactual_examples.pdf")
plt.close()
print("saved fig_track6c_counterfactual_examples")

print("\n[done] all figures generated")
