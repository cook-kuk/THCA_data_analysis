#!/usr/bin/env python3
"""Re-run only the figure generation (assumes wave10_mega_results.tsv exists)."""
import sys
import json
import numpy as np
import pandas as pd
from pathlib import Path

ROOT = Path("/home/seungho/personal/THCA_data_analysis/project/results/p_neo_bayesian_2026_05_09")
sys.path.insert(0, str(ROOT / "wave10"))
import build_wave10_mega as bw  # noqa: E402

OUT = bw.OUT
long_df = pd.read_csv(OUT / "wave10_mega_results.tsv", sep="\t")
print(f"loaded long_df: {long_df.shape}")

# Reload predictions for figs that need them
predictions, _ = bw.collect_all_predictions()

# Skip fig1 (already done) — but redo it anyway since it's quick
bw.fig1_megaheatmap(long_df, OUT / "fig_wave10_megaheatmap.png")
print("fig1 done")
bw.fig2_smallmultiples(long_df, OUT / "fig_wave10_smallmultiples.png")
print("fig2 done")
inflation_df = bw.fig3_inflation_diagonal(long_df, OUT / "fig_wave10_inflation_diagonal.png")
print("fig3 done")
bw.fig4_per_allele(long_df, OUT / "fig_wave10_per_allele.png")
print("fig4 done")
bw.fig5_calibration(long_df, OUT / "fig_wave10_calibration.png")
print("fig5 done")
sig_long = bw.fig6_significance(predictions, OUT / "fig_wave10_significance_matrix.png")
print(f"fig6 done — {len(sig_long)} pairs")
bw.fig7_selective_topK(predictions, OUT / "fig_wave10_selective_topK.png")
print("fig7 done")

# Recompute summary
summary = {}
no = long_df[long_df["test_set"] == "ITSNdb_no_overlap"].sort_values("AUROC", ascending=False)
summary["top3_no_overlap"] = no.head(3)[["method", "AUROC", "AUROC_lo95", "AUROC_hi95"]].to_dict("records")

CORE_TESTSETS = ["ITSNdb_no_overlap", "ITSNdb_in_master", "ITSNdb_combined",
                 "ITSNdb_main", "ITSNdb_Val"]
pivot_au = long_df.pivot_table(index="method", columns="test_set", values="AUROC", aggfunc="mean")
avail = pivot_au[[c for c in CORE_TESTSETS if c in pivot_au.columns]].notna().sum(axis=1)
keep_methods = avail[avail >= 3].index
avg_au = pivot_au.loc[keep_methods, [c for c in CORE_TESTSETS if c in pivot_au.columns]].mean(axis=1).sort_values(ascending=False)
summary["top3_robust_avg_core"] = [{"method": m, "avg_AUROC": float(avg_au[m])} for m in avg_au.head(3).index]

cal = long_df[(long_df["test_set"] == "ITSNdb_no_overlap") & long_df["ECE"].notna()].sort_values("ECE")
summary["top3_calibrated_no_overlap"] = cal.head(3)[["method", "ECE", "Brier"]].to_dict("records")

if inflation_df is not None and not inflation_df.empty:
    zero_inf = inflation_df[inflation_df["delta"].abs() < 0.05].sort_values("delta", key=abs)
    summary["zero_inflation"] = zero_inf[["method", "no_overlap", "in_master", "delta"]].to_dict("records")
    massive = inflation_df[inflation_df["delta"] > 0.4].sort_values("delta", ascending=False)
    summary["massive_inflation"] = massive[["method", "no_overlap", "in_master", "delta"]].to_dict("records")
else:
    summary["zero_inflation"] = []
    summary["massive_inflation"] = []

summary["mega_dimensions"] = {
    "n_methods": int(long_df["method"].nunique()),
    "n_test_sets": int(long_df["test_set"].nunique()),
    "n_metrics": 7,
    "n_long_rows": int(len(long_df)),
}

with open(OUT / "wave10_summary.json", "w") as f:
    json.dump(summary, f, indent=2, default=str)
print(json.dumps(summary, indent=2, default=str))
