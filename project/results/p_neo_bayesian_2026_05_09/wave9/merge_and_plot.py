#!/usr/bin/env python
"""Merge wave3 (6 algorithms) + wave9 (3 new) AUROC summaries into a single
9-algorithm comparison table and forest plot.
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams["font.family"] = "DejaVu Sans"

ROOT = Path("/home/seungho/personal/THCA_data_analysis/project/results/p_neo_bayesian_2026_05_09")
OUT = ROOT / "wave9"

# ---- prior 6 algorithms (wave3) ----
prior = {
    "MHCflurry":   ROOT / "wave3_mhcflurry/auroc_summary.tsv",
    "BigMHC":      ROOT / "wave3_bigmhc/auroc_summary.tsv",
    "DeepImmuno":  ROOT / "wave3_deepimmuno/auroc_summary.tsv",
    "PRIME":       ROOT / "wave3_prime/auroc_summary.tsv",
    "TransPHLA":   ROOT / "wave3_transphla/auroc_summary.tsv",
    "NetMHCpan":   ROOT / "wave3_netmhcpan/auroc_summary.tsv",
}

# ---- wave9 new ----
new = {
    "MHCnuggets":     OUT / "auroc_mhcnuggets.tsv",
    "T-SCAPE":        OUT / "auroc_tscape.tsv",
    "NetMHCstabpan":  OUT / "auroc_netmhcstabpan.tsv",
}


def load_one(name: str, path: Path) -> pd.DataFrame:
    if not path.exists():
        return None
    d = pd.read_csv(path, sep="\t")
    d["algorithm"] = name
    # Standardize testset column names
    keep = ["algorithm", "testset", "n", "n_pos", "AUROC", "AUROC_lo95", "AUROC_hi95"]
    d = d[[c for c in keep if c in d.columns]]
    return d


rows = []
for name, p in prior.items():
    df = load_one(name, p)
    if df is None:
        print(f"[merge] missing {name}: {p}")
        continue
    rows.append(df)
for name, p in new.items():
    df = load_one(name, p)
    if df is None:
        print(f"[merge] missing {name}: {p}")
        continue
    rows.append(df)

combined = pd.concat(rows, ignore_index=True)
combined.to_csv(OUT / "wave9_combined_results.tsv", sep="\t", index=False)
print(f"[merge] wrote wave9_combined_results.tsv ({len(combined)} rows)")
print(combined.to_string(index=False))


# Pivot for forest: rows = algorithm, cols by testset
forest = combined.pivot_table(index="algorithm", columns="testset",
                              values=["AUROC", "AUROC_lo95", "AUROC_hi95"])
forest.to_csv(OUT / "wave9_forest_pivot.tsv", sep="\t")
print("\n[merge] forest pivot:")
print(forest)


# ---- Forest plot ----
# Order: by no_overlap AUROC descending (most generalizing on top)
no_ov = combined[combined["testset"] == "ITSNdb_no_overlap"].copy()
no_ov = no_ov.sort_values("AUROC", ascending=True).reset_index(drop=True)

algs_order = no_ov["algorithm"].tolist()
testsets = ["ITSNdb_in_master", "ITSNdb_combined", "ITSNdb_no_overlap"]
ts_color = {"ITSNdb_in_master": "#d62728",   # leakage-suspect (red)
            "ITSNdb_combined":  "#7f7f7f",   # combined (grey)
            "ITSNdb_no_overlap": "#1f77b4"}  # truly external (blue)
ts_label = {"ITSNdb_in_master": "in_master (leakage-suspect)",
            "ITSNdb_combined":  "combined",
            "ITSNdb_no_overlap": "no_overlap (clean)"}

fig, ax = plt.subplots(figsize=(9, 5.5))
n_alg = len(algs_order)
y_off = {"ITSNdb_in_master": -0.25, "ITSNdb_combined": 0.0, "ITSNdb_no_overlap": +0.25}

for ts in testsets:
    sub = combined[combined["testset"] == ts].set_index("algorithm")
    if sub.empty:
        continue
    xs, los, his, ys = [], [], [], []
    for i, a in enumerate(algs_order):
        if a not in sub.index:
            continue
        row = sub.loc[a]
        xs.append(row["AUROC"])
        los.append(row["AUROC"] - row["AUROC_lo95"])
        his.append(row["AUROC_hi95"] - row["AUROC"])
        ys.append(i + y_off[ts])
    ax.errorbar(xs, ys, xerr=[los, his], fmt="o", color=ts_color[ts],
                ecolor=ts_color[ts], capsize=2.5, ms=5, lw=1.2, label=ts_label[ts])

ax.axvline(0.5, color="black", lw=0.7, ls="--", alpha=0.5)
ax.set_yticks(range(n_alg))
ax.set_yticklabels(algs_order)
ax.set_xlabel("AUROC (1000-bootstrap 95% CI)")
ax.set_xlim(0.30, 1.0)
ax.set_title("Wave 9 — 9-algorithm leakage-stratified ITSNdb forest\n(prior 6 + 3 new: MHCnuggets / T-SCAPE / NetMHCstabpan)")
ax.legend(loc="lower right", frameon=True, fontsize=9)
# Highlight the no_overlap baseline horizontal lines
for i in range(n_alg):
    ax.axhline(i, color="lightgrey", lw=0.4, alpha=0.4, zorder=0)
ax.grid(axis="x", alpha=0.3)
plt.tight_layout()
plt.savefig(OUT / "fig_wave9_9algorithm_forest.png", dpi=150)
plt.savefig(OUT / "fig_wave9_9algorithm_forest.pdf")
print(f"[merge] saved fig_wave9_9algorithm_forest.{{png,pdf}}")


# Compact 9-algorithm summary
summary = combined.pivot_table(index="algorithm", columns="testset", values="AUROC")
# Order rows by no_overlap AUROC descending
summary = summary.reindex(no_ov.sort_values("AUROC", ascending=False)["algorithm"].tolist())
summary.to_csv(OUT / "wave9_summary_table.tsv", sep="\t")
print("\n[merge] summary AUROC table (rows ordered by no_overlap descending):")
print(summary.round(3))
