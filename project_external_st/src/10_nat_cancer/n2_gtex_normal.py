#!/usr/bin/env python3
"""N2 — GTEx normal thyroid baseline.
Score RAI_8 / DM1_like / NONOVERLAP / TROP2 in GTEx normal thyroid → confirm:
- Normal thyroid has VERY LOW DM1_like (high RAI_8, baseline reference)
- Normal thyroid has VERY LOW TROP2 (= sacituzumab govitecan target only in cancer)
"""
from pathlib import Path
import gzip
import sys
import urllib.request
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

OUT = Path("project_external_st/results/extra")

# Pancan EBPP includes both TCGA tumor and TCGA normal-adjacent (sample type 11) — usable as quasi-normal
# True GTEx requires separate download. Use TCGA-THCA normal-adjacent as fastest baseline.
# Sample type code: -01 primary tumor, -11 normal solid tissue

TARGETS = ["TPO","DIO1","TSHR","PAX8","TG","FOXE1","NKX2-1","SLC5A5","NKX2_1","TITF1",
           "SLC26A4","IYD","DUOX1","DUOX2","TFF3","HHEX","GLIS3","DIO2",
           "TACSTD2","FN1","KCNN4","NAMPT","LYN","CYP1B1","CREB5","ELF3"]
ALL_GENES = set(TARGETS)

# Stream pancan, get all THCA samples (-01 tumor + -11 normal-adj)
clin = pd.read_csv("project/results/dark_matter_phase2/web/data/tcga_dm_master_with_pfi.tsv", sep="\t")
thca_patients = set(clin["tcga_short"])

print("Streaming pancan for THCA tumor + normal samples...")
rows = {}
with gzip.open("project/data/raw/TCGA_pancan/pancan_geneExp.gz","rt") as f:
    header = next(f).strip().split("\t")
    samples = header[1:]
    keep_idx = []
    for i, s in enumerate(samples):
        parts = s.split("-")
        if len(parts) >= 4 and "-".join(parts[:3]) in thca_patients:
            keep_idx.append(i)
    keep_samples = [samples[i] for i in keep_idx]
    print(f"  THCA samples (any type): {len(keep_samples)}")
    sample_types = [s.split("-")[3][:2] for s in keep_samples]
    n_tumor = sum(1 for t in sample_types if t == "01")
    n_normal = sum(1 for t in sample_types if t == "11")
    print(f"    primary tumor (-01): {n_tumor}")
    print(f"    normal-adjacent (-11): {n_normal}")
    for line in f:
        parts = line.rstrip("\n").split("\t")
        sym = parts[0].strip()
        if sym in ALL_GENES:
            rows[sym] = [float(parts[i+1]) if parts[i+1] not in ("","NA","NaN") else np.nan
                          for i in keep_idx]
expr = pd.DataFrame(rows, index=keep_samples).T
if "NKX2_1" in expr.index and "NKX2-1" not in expr.index:
    expr = expr.rename(index={"NKX2_1":"NKX2-1"})
elif "TITF1" in expr.index and "NKX2-1" not in expr.index:
    expr = expr.rename(index={"TITF1":"NKX2-1"})
expr = expr[~expr.index.duplicated(keep="first")]
print(f"  expr: {expr.shape[0]} genes × {expr.shape[1]} samples")

# Identify samples by type
sample_meta = pd.DataFrame({"sample": keep_samples, "type": sample_types})
sample_meta["category"] = sample_meta["type"].map({"01":"primary_tumor","11":"normal_adjacent",
                                                    "06":"metastatic"})
sample_meta = sample_meta.dropna(subset=["category"])
sample_meta = sample_meta.set_index("sample")

# Score
def score(genes, ref_pool):
    avail = [g for g in genes if g in expr.index]
    if not avail: return None
    z = expr.loc[avail].sub(expr.loc[avail, ref_pool].mean(axis=1), axis=0).div(
        expr.loc[avail, ref_pool].std(axis=1).replace(0, np.nan), axis=0)
    return z.mean(axis=0)

# Use ALL THCA samples as reference for z-score
ref = sample_meta.index.tolist()
sample_meta["RAI_8"] = score(["TPO","DIO1","TSHR","PAX8","TG","FOXE1","NKX2-1","SLC5A5"], ref)
sample_meta["NONOVERLAP"] = score(["SLC26A4","IYD","DUOX1","DUOX2","TFF3","HHEX","GLIS3","DIO2"], ref)
sample_meta["DM1_like"] = -sample_meta["RAI_8"]
sample_meta["TACSTD2_z"] = score(["TACSTD2"], ref)
sample_meta["FN1_z"] = score(["FN1"], ref)
sample_meta["TACSTD2_raw"] = expr.loc["TACSTD2"] if "TACSTD2" in expr.index else np.nan

# Group comparison
print("\n=== TCGA-THCA tumor vs normal-adjacent ===")
print(sample_meta.groupby("category")[["DM1_like","RAI_8","TACSTD2_z","TACSTD2_raw"]].agg(["mean","std","count"]).to_string())

# Mann-Whitney
from scipy.stats import mannwhitneyu
tumor = sample_meta[sample_meta["category"]=="primary_tumor"]
normal = sample_meta[sample_meta["category"]=="normal_adjacent"]
for col in ["DM1_like","RAI_8","TACSTD2_z","TACSTD2_raw"]:
    if normal[col].notna().any():
        a = tumor[col].dropna(); b = normal[col].dropna()
        try:
            U, p = mannwhitneyu(a, b)
            print(f"  {col}: tumor mean = {a.mean():.3f}, normal mean = {b.mean():.3f}, U-test p = {p:.2e}")
        except: pass

sample_meta.to_csv(OUT / "n2_thca_tumor_vs_normal.tsv", sep="\t")

# Figure
fig, axes = plt.subplots(1, 3, figsize=(15, 5))
import seaborn as sns
order = ["normal_adjacent","primary_tumor","metastatic"]
order = [c for c in order if (sample_meta["category"]==c).sum() > 0]
for ax, col, label in zip(axes, ["DM1_like","TACSTD2_z","TACSTD2_raw"],
                          ["DM1_like","TACSTD2 (z)","TACSTD2 (raw EBPP)"]):
    sns.boxplot(data=sample_meta, x="category", y=col, order=order, ax=ax,
                palette={"normal_adjacent":"#3C6B4F","primary_tumor":"#962E2E",
                         "metastatic":"#7B1F2A"})
    sns.stripplot(data=sample_meta, x="category", y=col, order=order, ax=ax,
                  color="black", size=3, alpha=0.4)
    ax.set_title(label, fontsize=11)
    ax.set_xlabel("")
    for i, c in enumerate(order):
        n = (sample_meta["category"]==c).sum()
        ax.text(i, ax.get_ylim()[0], f"n={n}", ha="center", fontsize=9)
fig.suptitle("N2 — TCGA-THCA tumor vs normal-adjacent thyroid (quasi-GTEx baseline)\n"
             f"DM1_like and TACSTD2 expected to be elevated in tumor vs normal",
             fontsize=12, y=1.02)
fig.tight_layout()
fig.savefig(OUT / "n2_tumor_vs_normal.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"→ {OUT / 'n2_tumor_vs_normal.png'}")
