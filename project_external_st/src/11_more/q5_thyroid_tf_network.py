#!/usr/bin/env python3
"""Q5 — Thyroid TF network coordination check.
Test if FOXE1, NKX2-1, PAX8, HHEX (4 thyroid lineage TFs) collapse coordinately
or independently in DM1-high. Coordinated collapse = stronger mechanism story."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

OUT = Path("project_external_st/results/extra")

# Load TF activity matrix from M sprint output
diff = pd.read_csv(OUT / "m_tf_activity_diff.tsv", sep="\t")

# Re-load TF activity per sample (need to re-run minimal)
# Actually we have m_master_regulator output; the per-sample TF activity isn't saved explicitly
# Let me check if there's an intermediate TF activity matrix saved
import glob
print("Available extra files:", [Path(p).name for p in sorted(glob.glob(str(OUT / "*")))][:30])

# If not saved, derive from CCLE-like analysis on TCGA expression
# Quick approach: just use the score TSV + expression load again for these 4 TFs
import gzip
TFs = ["FOXE1","NKX2-1","NKX2_1","TITF1","PAX8","HHEX","TG","TPO","TSHR","SLC5A5","DIO1","TRPS1","HOXB3","STAT3","FOSL1","JUNB","DNMT1","DNMT3B"]
ALL_TF = set(TFs)

print("Streaming pancan for thyroid TF/target genes...")
rows_dict = {}
with gzip.open("project/data/raw/TCGA_pancan/pancan_geneExp.gz","rt") as f:
    header = next(f).strip().split("\t")
    samples = header[1:]
    for line in f:
        parts = line.rstrip("\n").split("\t")
        sym = parts[0].strip()
        if sym in ALL_TF:
            rows_dict[sym] = [float(x) if x not in ("","NA","NaN") else np.nan for x in parts[1:]]
expr = pd.DataFrame(rows_dict, index=samples).T
if "NKX2_1" in expr.index and "NKX2-1" not in expr.index:
    expr = expr.rename(index={"NKX2_1":"NKX2-1"})
elif "TITF1" in expr.index and "NKX2-1" not in expr.index:
    expr = expr.rename(index={"TITF1":"NKX2-1"})
expr = expr[~expr.index.duplicated(keep="first")]
print(f"  expr {expr.shape[0]} TF genes × {expr.shape[1]} samples")

# THCA samples
score = pd.read_csv(OUT / "s_tcga_thca_scored.tsv", sep="\t")
thca_patients = set(score["patient"])
sample_to_pt = {s: "-".join(s.split("-")[:3]) for s in expr.columns}
thca_samples = [s for s,p in sample_to_pt.items() if p in thca_patients]
expr = expr[thca_samples]
print(f"  THCA: {expr.shape}")

# Pairwise correlation matrix among thyroid TFs
thyroid_tfs = ["FOXE1","NKX2-1","PAX8","HHEX","TG","TPO","TSHR","SLC5A5","DIO1","TRPS1"]
de_tfs = ["STAT3","FOSL1","JUNB","DNMT1","DNMT3B"]
keep = [t for t in (thyroid_tfs + de_tfs) if t in expr.index]

sub = expr.loc[keep]
cor = sub.T.corr()  # samples-side correlation; need genes-side
cor_genes = sub.T.corr()

print(f"\n=== Pairwise correlations among thyroid lineage genes (TCGA-THCA n={expr.shape[1]}) ===")
print(cor_genes.round(2).to_string())

# DM1 stratified
patient_to_dm1 = score.drop_duplicates("patient").set_index("patient")["DM1_like"]
sample_dm1 = pd.Series({s: float(patient_to_dm1[p]) for s,p in sample_to_pt.items()
                        if s in expr.columns and p in patient_to_dm1.index})
hi_thr = sample_dm1.quantile(0.75); lo_thr = sample_dm1.quantile(0.25)
hi_samples = sample_dm1[sample_dm1 >= hi_thr].index.intersection(expr.columns)
lo_samples = sample_dm1[sample_dm1 <= lo_thr].index.intersection(expr.columns)

print(f"\n  DM1-high (top 25%): n={len(hi_samples)}")
print(f"  DM1-low (bot 25%): n={len(lo_samples)}")

# Mean expression per group
mean_hi = expr[hi_samples].mean(axis=1)
mean_lo = expr[lo_samples].mean(axis=1)
delta = mean_hi - mean_lo

print(f"\n=== Mean Δ (DM1-high − DM1-low) per gene ===")
result = pd.DataFrame({"mean_DM1high": mean_hi, "mean_DM1low": mean_lo, "delta": delta})
result = result.reindex(keep).sort_values("delta")
print(result.to_string())

result.to_csv(OUT / "q5_thyroid_tf_network.tsv", sep="\t")

# Figure: heatmap
fig, axes = plt.subplots(1, 2, figsize=(15, 6))
import seaborn as sns

# panel A: correlation heatmap (thyroid TFs)
ax = axes[0]
order_a = thyroid_tfs
order_a = [t for t in order_a if t in cor_genes.index]
sns.heatmap(cor_genes.loc[order_a, order_a], cmap="RdBu_r", center=0, vmin=-1, vmax=1,
            annot=True, fmt=".2f", ax=ax, cbar_kws={"label":"Pearson r"})
ax.set_title("Q5.A — Thyroid lineage TF/target co-expression\n"
             "(TCGA-THCA n=" + str(expr.shape[1]) + " samples)\n"
             "Coordinated collapse = high pairwise positive correlation",
             fontsize=11)

# panel B: bar plot of Δ DM1-high − DM1-low
ax = axes[1]
order_b = result.index.tolist()
colors = ["#3C6B4F" if d < 0 else "#962E2E" for d in result["delta"]]
ax.barh(order_b, result["delta"], color=colors, edgecolor="black")
ax.axvline(0, color="grey", lw=0.5, ls=":")
ax.set_xlabel("mean expression Δ (DM1-high − DM1-low)")
ax.set_title("Q5.B — Coordinated TF / target shift\n"
             "Green = down in DM1-high (lineage collapse), Red = up", fontsize=11)
fig.tight_layout()
fig.savefig(OUT / "q5_thyroid_tf_network.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"→ {OUT / 'q5_thyroid_tf_network.png'}")
