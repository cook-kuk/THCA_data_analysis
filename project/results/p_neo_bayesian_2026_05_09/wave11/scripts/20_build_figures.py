#!/usr/bin/env python3
"""Wave 11 — build the 4 mega-figures from wave11_mega_table_expanded.tsv.

Fig 1: heatmap (algorithm × test set, color = external AUROC)
Fig 2: inflation profile per algo (Δ = in_training − external)
Fig 3: per-test-set top-3 robustness ranking
Fig 4: N-effective vs best-AUROC scatter
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path("/home/seungho/personal/THCA_data_analysis/project/results/p_neo_bayesian_2026_05_09")
WAVE11 = ROOT / "wave11"
FIGS = WAVE11 / "figs"
FIGS.mkdir(parents=True, exist_ok=True)

mega = pd.read_csv(WAVE11 / "wave11_mega_table_expanded.tsv", sep="\t")
print(f"loaded mega: {mega.shape}", flush=True)

# ============================================================
# Fig 1 — multi-test-set heatmap (external AUROC)
# ============================================================
piv = mega.pivot_table(index="algorithm", columns="test_bundle", values="external_AUROC", aggfunc="first")
# fall back to combined AUROC if external is NaN
piv_full = mega.pivot_table(index="algorithm", columns="test_bundle", values="AUROC", aggfunc="first")
piv = piv.combine_first(piv_full)

# sort algos by mean external AUROC (ignoring NaN)
order = piv.mean(axis=1).sort_values(ascending=False).index
piv = piv.loc[order]

fig, ax = plt.subplots(figsize=(max(8, piv.shape[1]*1.0), max(4, piv.shape[0]*0.45)))
im = ax.imshow(piv.values, cmap="RdYlGn", vmin=0.4, vmax=0.85, aspect="auto")
ax.set_xticks(range(piv.shape[1])); ax.set_xticklabels(piv.columns, rotation=45, ha="right")
ax.set_yticks(range(piv.shape[0])); ax.set_yticklabels(piv.index)
for i in range(piv.shape[0]):
    for j in range(piv.shape[1]):
        v = piv.values[i, j]
        if not np.isnan(v):
            ax.text(j, i, f"{v:.3f}", ha="center", va="center",
                    color="black" if 0.55 <= v <= 0.75 else "white", fontsize=8)
plt.colorbar(im, ax=ax, label="External AUROC")
ax.set_title("Wave 11 — algorithm × test-set external AUROC heatmap")
plt.tight_layout()
plt.savefig(FIGS / "fig_wave11_megaheatmap.png", dpi=150, bbox_inches="tight")
plt.savefig(FIGS / "fig_wave11_megaheatmap.pdf", bbox_inches="tight")
plt.close()
print("[fig1] heatmap done", flush=True)

# ============================================================
# Fig 2 — inflation profile (Δ = in_training - external) per algo per bundle
# ============================================================
piv2 = mega.pivot_table(index="algorithm", columns="test_bundle",
                         values="inflation_delta", aggfunc="first")
# Reindex with shared order; missing rows are filled with NaN
piv2 = piv2.reindex(piv.index)

fig, ax = plt.subplots(figsize=(max(8, piv2.shape[1]*1.0), max(4, piv2.shape[0]*0.45)))
im = ax.imshow(piv2.values, cmap="RdBu_r", vmin=-0.3, vmax=0.3, aspect="auto")
ax.set_xticks(range(piv2.shape[1])); ax.set_xticklabels(piv2.columns, rotation=45, ha="right")
ax.set_yticks(range(piv2.shape[0])); ax.set_yticklabels(piv2.index)
for i in range(piv2.shape[0]):
    for j in range(piv2.shape[1]):
        v = piv2.values[i, j]
        if not np.isnan(v):
            ax.text(j, i, f"{v:+.2f}", ha="center", va="center",
                    color="black" if abs(v) < 0.15 else "white", fontsize=8)
plt.colorbar(im, ax=ax, label="Δ AUROC (in_training − external)")
ax.set_title("Wave 11 — leakage inflation profile per algorithm × test set")
plt.tight_layout()
plt.savefig(FIGS / "fig_wave11_inflation_per_set.png", dpi=150, bbox_inches="tight")
plt.savefig(FIGS / "fig_wave11_inflation_per_set.pdf", bbox_inches="tight")
plt.close()
print("[fig2] inflation done", flush=True)

# ============================================================
# Fig 3 — per-test-set top-3 ranking
# ============================================================
test_sets = sorted([c for c in piv.columns])
n_sets = len(test_sets)
ncols = min(4, n_sets)
nrows = int(np.ceil(n_sets / ncols))
fig, axes = plt.subplots(nrows, ncols, figsize=(ncols*4.0, nrows*2.5))
axes = np.atleast_1d(axes).flatten()

for ax, ts in zip(axes, test_sets):
    sub = mega[mega["test_bundle"] == ts].copy()
    sub = sub.dropna(subset=["external_AUROC"]) if not sub["external_AUROC"].isna().all() else sub.dropna(subset=["AUROC"])
    score_col = "external_AUROC" if not sub["external_AUROC"].isna().all() else "AUROC"
    sub = sub.sort_values(score_col, ascending=True).tail(5)
    if len(sub) == 0:
        ax.axis("off"); ax.set_title(ts, fontsize=10); continue
    colors = ["#3182bd"] * len(sub)
    if len(sub) >= 1: colors[-1] = "#e34a33"
    ax.barh(sub["algorithm"], sub[score_col], color=colors)
    for i, v in enumerate(sub[score_col]):
        ax.text(v + 0.005, i, f"{v:.3f}", va="center", fontsize=8)
    ax.set_xlim(0.4, 0.95)
    ax.axvline(0.5, color="gray", lw=0.5, ls="--")
    ax.set_title(ts, fontsize=10)
    ax.tick_params(axis="y", labelsize=8)
for ax in axes[len(test_sets):]:
    ax.axis("off")
plt.suptitle("Wave 11 — per-test-set top-5 algorithm ranking (external AUROC)", y=1.005)
plt.tight_layout()
plt.savefig(FIGS / "fig_wave11_robustness_ranking.png", dpi=150, bbox_inches="tight")
plt.savefig(FIGS / "fig_wave11_robustness_ranking.pdf", bbox_inches="tight")
plt.close()
print("[fig3] ranking done", flush=True)

# ============================================================
# Fig 4 — N-effective vs best-algo AUROC
# ============================================================
ts_summary = mega.groupby("test_bundle").agg(
    n_samples=("n", "max"),
    n_pos=("n_pos", "max"),
    best_external_auroc=("external_AUROC", "max"),
    best_overall_auroc=("AUROC", "max"),
    n_algos=("algorithm", "count"),
).reset_index()
ts_summary["n_eff"] = 4 * ts_summary["n_pos"] * (ts_summary["n_samples"] - ts_summary["n_pos"]) / ts_summary["n_samples"].clip(lower=1)
ts_summary["best_auroc"] = ts_summary["best_external_auroc"].fillna(ts_summary["best_overall_auroc"])

fig, ax = plt.subplots(figsize=(8, 5.5))
sub = ts_summary.dropna(subset=["best_auroc"])
ax.scatter(sub["n_eff"], sub["best_auroc"], s=80, c="#1f77b4", edgecolor="black")
for _, r in sub.iterrows():
    ax.annotate(r["test_bundle"], (r["n_eff"], r["best_auroc"]),
                xytext=(5, 5), textcoords="offset points", fontsize=9)
ax.set_xscale("log")
ax.set_xlabel("Effective n (4·pos·neg/N)")
ax.set_ylabel("Best-algorithm AUROC")
ax.axhline(0.5, color="gray", ls="--", lw=0.5)
ax.set_title("Wave 11 — small-n test sets are unreliable")
plt.tight_layout()
plt.savefig(FIGS / "fig_wave11_neff_vs_auroc.png", dpi=150, bbox_inches="tight")
plt.savefig(FIGS / "fig_wave11_neff_vs_auroc.pdf", bbox_inches="tight")
plt.close()
print("[fig4] neff scatter done", flush=True)

print("\n[done] 4 figures written to wave11/figs/")
ts_summary.to_csv(WAVE11 / "wave11_test_set_summary.tsv", sep="\t", index=False)
