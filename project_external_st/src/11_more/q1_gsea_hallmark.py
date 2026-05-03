#!/usr/bin/env python3
"""Q1 — GSEA Hallmark enrichment on DM1-high vs DM1-low DE."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

OUT = Path("project_external_st/results/extra")

de = pd.read_csv(OUT / "t_DM1_DE_genes.tsv", sep="\t")
print(f"DE genes: {len(de)}")

# decoupler hallmark
import decoupler as dc
print("Loading MSigDB Hallmark...")
hallmark = dc.op.hallmark()
print(f"  Hallmark: {hallmark['source'].nunique()} pathways, {len(hallmark)} gene-pathway memberships")

# Build ranked stat as t-statistic, then ULM enrichment
de_ranked = de.dropna(subset=["t","gene"]).copy()
de_ranked = de_ranked.sort_values("t", ascending=False)

# Use ULM with t-statistic as input
mat = pd.DataFrame({"DM1_high_vs_low": de_ranked.set_index("gene")["t"]}).T
# decoupler expects samples x genes, here we have 1 sample (the contrast)
print(f"  ranking matrix: {mat.shape}")

result = dc.mt.ulm(data=mat, net=hallmark)
if isinstance(result, tuple):
    estimates = result[0]
    pvals = result[1] if len(result) > 1 else None
elif isinstance(result, dict):
    estimates = result.get("ulm_estimate", list(result.values())[0])
    pvals = result.get("ulm_pvals", None)
else:
    estimates = result; pvals = None

if hasattr(estimates, "to_df"):
    estimates = estimates.to_df()
if pvals is not None and hasattr(pvals, "to_df"):
    pvals = pvals.to_df()

# decoupler v2 sometimes returns score matrix (samples x pathways)
print(f"  estimates shape: {estimates.shape}")
res = pd.DataFrame({"score": estimates.iloc[0],
                    "abs_score": estimates.iloc[0].abs()}).sort_values("score", ascending=False)
if pvals is not None:
    res["p"] = pvals.iloc[0]
print("\n=== Top 15 enriched in DM1-high (positive score) ===")
print(res.head(15).to_string())
print("\n=== Top 15 depleted in DM1-high (negative score) ===")
print(res.tail(15).to_string())
res.to_csv(OUT / "q1_gsea_hallmark.tsv", sep="\t")

# Figure
fig, ax = plt.subplots(figsize=(9, 10))
n = 20
top = pd.concat([res.head(n//2), res.tail(n//2)])
colors = ["#962E2E" if s > 0 else "#3C6B4F" for s in top["score"]]
ax.barh(top.index[::-1], top["score"][::-1], color=colors[::-1], edgecolor="black")
ax.axvline(0, color="grey", lw=0.5, ls=":")
ax.set_xlabel("ULM score (DM1-high vs DM1-low t-statistic)")
ax.set_title("Q1 — MSigDB Hallmark enrichment in DM1-high (TCGA-THCA n=261)\n"
             "Red = enriched in DM1-high; Green = enriched in DM1-low",
             fontsize=11)
fig.tight_layout()
fig.savefig(OUT / "q1_gsea_hallmark.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"→ {OUT / 'q1_gsea_hallmark.png'}")
