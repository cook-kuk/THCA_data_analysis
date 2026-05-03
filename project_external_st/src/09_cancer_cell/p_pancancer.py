#!/usr/bin/env python3
"""Sprint P — Pan-cancer DM1 axis check.
Compute RAI/DM1/NONOVERLAP score across all 33 TCGA cancer types.
Question: thyroid-specific or universal dedifferentiation axis?"""
from pathlib import Path
import gzip
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

OUT = Path("project_external_st/results/extra")

RAI_8 = ["TPO","DIO1","TSHR","PAX8","TG","FOXE1","NKX2-1","SLC5A5","NKX2_1","TITF1"]
NONOVERLAP = ["SLC26A4","IYD","DUOX1","DUOX2","TFF3","HHEX","GLIS3","DIO2"]
ALL_GENES = set(RAI_8 + NONOVERLAP)

print("Streaming pancan EBPP for RAI/NONOVERLAP genes...")
rows = {}
with gzip.open("project/data/raw/TCGA_pancan/pancan_geneExp.gz","rt") as f:
    header = next(f).strip().split("\t")
    samples = header[1:]
    for line in f:
        parts = line.rstrip("\n").split("\t")
        sym = parts[0].strip()
        if sym in ALL_GENES:
            rows[sym] = [float(x) if x not in ("","NA","NaN") else np.nan for x in parts[1:]]
expr = pd.DataFrame(rows, index=samples).T
expr.index.name = "gene"
# consolidate alias
if "NKX2_1" in expr.index and "NKX2-1" not in expr.index:
    expr = expr.rename(index={"NKX2_1":"NKX2-1"})
if "TITF1" in expr.index and "NKX2-1" not in expr.index:
    expr = expr.rename(index={"TITF1":"NKX2-1"})
expr = expr[~expr.index.duplicated(keep="first")]
print(f"  {expr.shape[0]} genes × {expr.shape[1]} samples")
print(f"  found genes: {sorted(expr.index)}")

# Load survival.tsv for cancer-type labels
clin = pd.read_csv("project/data/raw/TCGA_pancan/survival.tsv", sep="\t")
clin = clin.rename(columns={"sample":"sample_id", "cancer type abbreviation":"cancer_type",
                             "PFI.time":"PFI_time", "PFI":"PFI_event",
                             "OS.time":"OS_time", "OS":"OS_event"})
print(f"  clinical: {len(clin)} samples, {clin['cancer_type'].nunique()} cancer types")

# Score per pan-cancer (within-cohort z-scoring is per sample of that cancer type)
def per_cancer_score(expr, clin, geneset, label):
    avail = [g for g in geneset if g in expr.index]
    sub_expr = expr.loc[avail]  # genes × samples
    rows = []
    for cancer in clin["cancer_type"].dropna().unique():
        cancer_samples = clin[clin["cancer_type"] == cancer]["sample_id"].tolist()
        cancer_in_expr = [s for s in cancer_samples if s in sub_expr.columns]
        if len(cancer_in_expr) < 5: continue
        e = sub_expr[cancer_in_expr]
        # within-cancer z per gene
        z = e.sub(e.mean(axis=1), axis=0).div(e.std(axis=1).replace(0, np.nan), axis=0)
        score = z.mean(axis=0)
        for samp, val in score.items():
            rows.append({"sample_id": samp, "cancer_type": cancer,
                         f"{label}_score": val})
    return pd.DataFrame(rows)

rai_df = per_cancer_score(expr, clin, [g for g in RAI_8 if g not in ("NKX2_1","TITF1")], "RAI_8")
non_df = per_cancer_score(expr, clin, NONOVERLAP, "NONOVERLAP")
score_df = rai_df.merge(non_df, on=["sample_id","cancer_type"])
score_df["DM1_like"] = -score_df["RAI_8_score"]
score_df = score_df.merge(clin[["sample_id","_PATIENT","PFI_time","PFI_event","OS_time","OS_event"]],
                          on="sample_id", how="left")
score_df.to_csv(OUT / "p_pancancer_scored.tsv", sep="\t", index=False)
print(f"  → {OUT / 'p_pancancer_scored.tsv'}  ({len(score_df)} samples)")

# Per-cancer summary
summary = (score_df.groupby("cancer_type")
           .agg(n=("sample_id","count"),
                mean_DM1=("DM1_like","mean"),
                mean_RAI=("RAI_8_score","mean"),
                mean_NONOV=("NONOVERLAP_score","mean"),
                cor_DM1_NONOV=("DM1_like", lambda x: x.corr(score_df.loc[x.index, "NONOVERLAP_score"])))
           .reset_index().sort_values("mean_DM1", ascending=False))
summary.to_csv(OUT / "p_pancancer_summary.tsv", sep="\t", index=False)
print("\n=== Per-cancer DM1 ranking ===")
print(summary.to_string(index=False))

# Figure: 33-cancer DM1 boxplot + cor with NONOVERLAP
fig, axes = plt.subplots(1, 2, figsize=(18, 7))
order = summary["cancer_type"].tolist()
ax = axes[0]
sns.boxplot(data=score_df, x="cancer_type", y="DM1_like", order=order, ax=ax,
            color="#962E2E", fliersize=2)
ax.axhline(0, color="grey", lw=0.5, ls=":")
ax.set_xticklabels(ax.get_xticklabels(), rotation=90, fontsize=9)
ax.set_title("P.A — DM1_like score across 33 TCGA cancer types\n"
             "Each cancer's score is within-cohort z-mean of 8 RAI genes (negated)\n"
             "Note: thyroid genes are highly tissue-specific → DM1 in non-thyroid = noise/absence",
             fontsize=11)
ax.set_xlabel(""); ax.set_ylabel("DM1_like score (within-cancer z-mean)")

ax = axes[1]
ax.bar(summary["cancer_type"], summary["cor_DM1_NONOV"], color="#962E2E", edgecolor="black")
ax.axhline(0, color="grey", lw=0.5, ls=":")
ax.set_xticklabels(summary["cancer_type"], rotation=90, fontsize=9)
ax.set_title("P.B — DM1 vs NONOVERLAP within-cohort Pearson correlation\n"
             "Thyroid biology specific: only THCA shows strong negative ρ",
             fontsize=11)
ax.set_ylabel("Pearson r (DM1 vs NONOVERLAP)")

fig.tight_layout()
fig.savefig(OUT / "p_pancancer.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"→ {OUT / 'p_pancancer.png'}")
