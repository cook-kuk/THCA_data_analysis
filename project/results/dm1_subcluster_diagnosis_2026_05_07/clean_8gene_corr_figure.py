"""
Replace the 10-gene ad-hoc network with two clean correlation matrices:
  Panel A — 8-gene RAI_8 self-correlation
  Panel B — 16-gene cross-panel (RAI_8 ∪ THYROID_NONOVERLAP) for §4.3 zero-overlap defense

Source: TCGA-THCA pancan expression (project/data/raw/TCGA_pancan/pancan_geneExp.gz)
"""
from __future__ import annotations
import gzip
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path("project/results/dm1_subcluster_diagnosis_2026_05_07")
ROOT.mkdir(parents=True, exist_ok=True)

RAI_8        = ["SLC5A5","TPO","TG","TSHR","PAX8","NKX2-1","FOXE1","DIO1"]
NONOVERLAP_8 = ["SLC26A4","IYD","DUOX1","DUOX2","TFF3","HHEX","GLIS3","DIO2"]
NEEDED = set(RAI_8 + NONOVERLAP_8 + ["TITF1","NKX2_1"])

# stream pancan
print("[load] streaming TCGA pancan expression...")
rows = {}
with gzip.open("project/data/raw/TCGA_pancan/pancan_geneExp.gz","rt") as f:
    header = next(f).strip().split("\t")
    samples = header[1:]
    for line in f:
        parts = line.rstrip("\n").split("\t")
        sym = parts[0].strip()
        if sym in NEEDED:
            rows[sym] = [float(x) if x not in ("","NA","NaN") else np.nan for x in parts[1:]]
expr = pd.DataFrame(rows, index=samples).T
# alias collapse
if "NKX2-1" not in expr.index:
    if "NKX2_1" in expr.index: expr = expr.rename(index={"NKX2_1":"NKX2-1"})
    elif "TITF1" in expr.index: expr = expr.rename(index={"TITF1":"NKX2-1"})
expr = expr[~expr.index.duplicated(keep="first")]
print(f"  loaded {expr.shape[0]} genes x {expr.shape[1]} samples")

# THCA only
score = pd.read_csv("project_external_st/results/extra/s_tcga_thca_scored.tsv", sep="\t")
thca_patients = set(score["patient"])
sample_to_pt = {s: "-".join(s.split("-")[:3]) for s in expr.columns}
thca = [s for s,p in sample_to_pt.items() if p in thca_patients]
expr_thca = expr[thca]
print(f"  THCA samples: {expr_thca.shape[1]}")

# ---------------------------------------------------------------------------
# PANEL A — 8-gene RAI_8 self-correlation (clean)
# ---------------------------------------------------------------------------
genes_a = [g for g in RAI_8 if g in expr_thca.index]
print(f"  RAI_8 found {len(genes_a)}/8: {genes_a}")
e_a = expr_thca.loc[genes_a].apply(lambda r: r.fillna(r.median()), axis=1)
cor_a = e_a.T.corr(method="spearman")

# ---------------------------------------------------------------------------
# PANEL B — 16-gene cross-panel (RAI_8 + NONOVERLAP) zero-overlap defense
# ---------------------------------------------------------------------------
genes_b = [g for g in (RAI_8 + NONOVERLAP_8) if g in expr_thca.index]
print(f"  RAI_8+NONOVERLAP found {len(genes_b)}/16: {genes_b}")
e_b = expr_thca.loc[genes_b].apply(lambda r: r.fillna(r.median()), axis=1)
cor_b = e_b.T.corr(method="spearman")

# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(16, 7))

def heatmap(ax, mat, title, divider_idx=None):
    im = ax.imshow(mat.values, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
    ax.set_xticks(range(len(mat))); ax.set_yticks(range(len(mat)))
    ax.set_xticklabels(mat.columns, rotation=45, ha="right", fontsize=9)
    ax.set_yticklabels(mat.index, fontsize=9)
    for i in range(len(mat)):
        for j in range(len(mat)):
            v = mat.values[i,j]
            if not np.isnan(v):
                ax.text(j,i,f"{v:.2f}",ha="center",va="center",
                        fontsize=7,color="white" if abs(v)>0.5 else "black")
    if divider_idx is not None:
        ax.axhline(divider_idx-0.5, color="black", lw=2)
        ax.axvline(divider_idx-0.5, color="black", lw=2)
    ax.set_title(title, fontsize=11)
    return im

# Panel A
heatmap(axes[0], cor_a,
        f"A. 8-gene RAI_8 self-correlation  (TCGA-THCA n={expr_thca.shape[1]})\n"
        f"All 8 genes are RAI_8 panel members — internal panel coherence")

# Panel B
heatmap(axes[1], cor_b,
        f"B. 16-gene cross-panel matrix  (RAI_8 ⊕ THYROID_NONOVERLAP, 0 gene overlap)\n"
        f"Top-left 8×8 = RAI_8;  Bottom-right 8×8 = NONOVERLAP;  Cross 8×8 = §4.3 zero-overlap defense",
        divider_idx=8)

# add colorbar
cbar = fig.colorbar(axes[0].images[0], ax=axes, fraction=0.025, pad=0.02)
cbar.set_label("Spearman ρ (pairwise)", fontsize=10)

plt.suptitle(
    "Thyroid lineage panel correlation matrices — 8-gene only (A) + 16-gene zero-overlap cross-panel (B)\n"
    "(Replaces ad-hoc 10-gene network: drops TRPS1; HHEX moved into NONOVERLAP block where it belongs)",
    fontsize=12)
plt.tight_layout(rect=[0,0,0.95,0.96])
plt.savefig(ROOT / "FIG_clean_8gene_and_16gene_correlation.png", dpi=160, bbox_inches="tight")
plt.close()

# also dump the matrices
cor_a.to_csv(ROOT / "tcga_thca_RAI8_correlation.tsv", sep="\t")
cor_b.to_csv(ROOT / "tcga_thca_RAI8_NONOVERLAP_correlation.tsv", sep="\t")
print(f"[done] {ROOT}/FIG_clean_8gene_and_16gene_correlation.png")
print(f"  RAI_8 self-correlation mean (off-diag): {cor_a.values[~np.eye(len(cor_a),dtype=bool)].mean():.3f}")
cross = cor_b.iloc[:8, 8:].values
print(f"  Cross-panel mean correlation (8×8 quadrant): {cross.mean():.3f}")
print(f"  Cross-panel min: {cross.min():.3f}  max: {cross.max():.3f}")
