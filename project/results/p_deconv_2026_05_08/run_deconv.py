"""
Bulk RNA-seq cell-type deconvolution via NNLS pseudobulk regression.

Reference: Lu 2023 single-cell (GSE193581_hvg_adata.h5ad) — 8 author cell types.
Bulk: TCGA-THCA log2(TPM+1) — 572 samples.
DM labels: dark_matter_phase2/p2d_per_sample_classification.tsv (137 DM1+DM2).

Output:
- per_sample_celltype_fractions.tsv (n_samples x 8 cell types)
- dm1_vs_dm2_celltype_d.tsv (per cell type Cohen's d, MW p)
- residualized_8gene_dm1_dm2.tsv (S4 expansion: 8-gene panel d before/after residualization on cell-type fractions)
"""
from __future__ import annotations
from pathlib import Path
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import anndata as ad
from scipy.optimize import nnls
from scipy import stats

OUT_DIR = Path("/home/seungho/personal/THCA_data_analysis/project/results/p_deconv_2026_05_08")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# 8-gene panel (Paper 1)
PANEL_8G = ["SLC5A5","TPO","TG","TSHR","PAX8","NKX2-1","FOXE1","DIO1"]

print("=" * 70)
print("STEP 1 / Reference build — Lu 2023 sc → pseudobulk (8 cell types)")
print("=" * 70)
ref = ad.read_h5ad("/data/thca/repo_results/v17_lu2023/GSE193581_hvg_adata.h5ad")
ct_col = "author_celltype"
print(f"sc cells: {ref.shape[0]}, HVG: {ref.shape[1]}, {ct_col}: {ref.obs[ct_col].nunique()} types")

# Mean expression per cell type (8 × 2000)
celltypes = sorted(ref.obs[ct_col].astype(str).unique())
print(f"Cell types: {celltypes}")
import scipy.sparse as sp
X = ref.X
if sp.issparse(X):
    X = X.toarray()
pseudobulk = pd.DataFrame(index=ref.var_names, columns=celltypes, dtype=float)
for ct in celltypes:
    mask = (ref.obs[ct_col].astype(str) == ct).values
    pseudobulk[ct] = X[mask].mean(axis=0)
print(f"pseudobulk: {pseudobulk.shape} (genes x celltypes)")

print()
print("=" * 70)
print("STEP 2 / Load TCGA-THCA bulk + intersect genes")
print("=" * 70)
tcga = pd.read_csv("/data/thca/data_processed/bulk_rnaseq/TCGA-THCA_rnaseq_expression_log2.tsv",
                    sep="\t", index_col=0)
print(f"TCGA bulk: {tcga.shape} (genes x samples)")

inter = sorted(set(pseudobulk.index) & set(tcga.index))
print(f"Intersect genes: {len(inter)}")
ref_X = pseudobulk.loc[inter].values   # (G, K=8)
bulk_X = tcga.loc[inter].values         # (G, S=572)
print(f"ref matrix (G, K): {ref_X.shape}, bulk matrix (G, S): {bulk_X.shape}")

# Both reference and TCGA bulk are log-scale (Lu .X is scaled, TCGA is log2(TPM+1))
# Per-gene z-score across samples + cell types together for comparability
# For NNLS deconvolution, simpler: normalize each column (cell-type and sample) to unit scale
# Standard approach: just use as-is (both log-scale); NNLS is robust to scale shifts as long as gene-relative profiles match.

print()
print("=" * 70)
print("STEP 3 / NNLS deconvolution per sample")
print("=" * 70)
sample_ids = list(tcga.columns)
fractions = np.zeros((len(sample_ids), len(celltypes)))
for i, sid in enumerate(sample_ids):
    y = bulk_X[:, i]
    w, residual = nnls(ref_X, y)
    if w.sum() > 0:
        w = w / w.sum()  # normalize to fractions
    fractions[i, :] = w
    if (i+1) % 100 == 0:
        print(f"  done {i+1}/{len(sample_ids)}")

frac_df = pd.DataFrame(fractions, index=sample_ids, columns=celltypes)
frac_df.index.name = "tcga_sample_id"
frac_df.to_csv(OUT_DIR / "per_sample_celltype_fractions.tsv", sep="\t")
print(f"Saved: per_sample_celltype_fractions.tsv  shape={frac_df.shape}")
print(f"\nMean fractions:")
print(frac_df.mean().round(4))

print()
print("=" * 70)
print("STEP 4 / DM1 vs DM2 cell-type fraction analysis")
print("=" * 70)
dm = pd.read_csv("/data/thca/repo_results/dark_matter_phase2/p2d_per_sample_classification.tsv", sep="\t")
print(f"DM table: {dm.shape}")
print(dm['v17_dark_cluster'].value_counts(dropna=False).to_string())

# Map TCGA-XX-XXXX (4-letter) form
def map_id(sid):
    # frac index is e.g., TCGA-DJ-A2Q6-01A. dm tcga_short is e.g., TCGA-4C-A93U
    return "-".join(sid.split("-")[:3])
frac_df["tcga_short"] = [map_id(s) for s in frac_df.index]
merged = dm.merge(frac_df, on="tcga_short", how="inner")
merged_dm = merged[merged['v17_dark_cluster'].isin(['DM1','DM2'])].copy()
print(f"merged DM1+DM2 samples: {len(merged_dm)} (DM1 {(merged_dm['v17_dark_cluster']=='DM1').sum()}, DM2 {(merged_dm['v17_dark_cluster']=='DM2').sum()})")

results = []
for ct in celltypes:
    dm1 = merged_dm.loc[merged_dm['v17_dark_cluster']=='DM1', ct].values
    dm2 = merged_dm.loc[merged_dm['v17_dark_cluster']=='DM2', ct].values
    n1, n2 = len(dm1), len(dm2)
    pooled_sd = np.sqrt(((n1-1)*dm1.var(ddof=1) + (n2-1)*dm2.var(ddof=1)) / (n1+n2-2))
    cohen_d = (dm1.mean() - dm2.mean()) / pooled_sd if pooled_sd > 0 else 0
    u_stat, p = stats.mannwhitneyu(dm1, dm2, alternative='two-sided')
    results.append({
        "cell_type": ct,
        "DM1_mean": dm1.mean(),
        "DM2_mean": dm2.mean(),
        "delta_DM1_minus_DM2": dm1.mean() - dm2.mean(),
        "cohen_d": cohen_d,
        "MW_p": p,
        "n_DM1": n1,
        "n_DM2": n2,
    })
celltype_d = pd.DataFrame(results).sort_values("cohen_d", key=abs, ascending=False)
celltype_d.to_csv(OUT_DIR / "dm1_vs_dm2_celltype_d.tsv", sep="\t", index=False)
print()
print("DM1 vs DM2 cell-type fraction analysis:")
print(celltype_d.round(4).to_string(index=False))

print()
print("=" * 70)
print("STEP 5 / 8-gene panel score residualization on cell-type fractions")
print("=" * 70)
# Compute 8-gene score (mean of z-scored panel) per TCGA sample
inter_panel = [g for g in PANEL_8G if g in tcga.index]
print(f"Panel genes in TCGA: {len(inter_panel)}/{len(PANEL_8G)}: {inter_panel}")
panel_mat = tcga.loc[inter_panel].T  # samples x panel
panel_z = (panel_mat - panel_mat.mean()) / panel_mat.std()
score_8g = panel_z.mean(axis=1)
score_8g.name = "panel_8g_z"
score_df = score_8g.to_frame()
score_df["tcga_short"] = [map_id(s) for s in score_df.index]
merged2 = merged.merge(score_df, on="tcga_short", how="inner")
m_dm = merged2[merged2['v17_dark_cluster'].isin(['DM1','DM2'])].copy()
print(f"DM1+DM2 with score: {len(m_dm)}")

# Raw d
raw_d_dm1 = m_dm.loc[m_dm['v17_dark_cluster']=='DM1', 'panel_8g_z'].values
raw_d_dm2 = m_dm.loc[m_dm['v17_dark_cluster']=='DM2', 'panel_8g_z'].values
n1, n2 = len(raw_d_dm1), len(raw_d_dm2)
pooled = np.sqrt(((n1-1)*raw_d_dm1.var(ddof=1) + (n2-1)*raw_d_dm2.var(ddof=1)) / (n1+n2-2))
raw_d = (raw_d_dm1.mean() - raw_d_dm2.mean()) / pooled
print(f"Raw 8-gene panel DM1 vs DM2 Cohen's d: {raw_d:.3f} (n_DM1={n1}, n_DM2={n2})")

# Residualize on each cell type / combinations
from sklearn.linear_model import LinearRegression
def residualize(y, X):
    lr = LinearRegression().fit(X, y)
    return y - lr.predict(X)

residual_results = [{"covariate": "raw (no residualization)", "cohen_d": raw_d, "delta_d": 0.0}]

# All celltypes individually
for ct in celltypes:
    X = m_dm[[ct]].values
    y = m_dm['panel_8g_z'].values
    y_resid = residualize(y, X)
    d1 = y_resid[m_dm['v17_dark_cluster'].values=='DM1']
    d2 = y_resid[m_dm['v17_dark_cluster'].values=='DM2']
    p = np.sqrt(((len(d1)-1)*d1.var(ddof=1) + (len(d2)-1)*d2.var(ddof=1)) / (len(d1)+len(d2)-2))
    cd = (d1.mean()-d2.mean())/p
    residual_results.append({
        "covariate": f"residualize on {ct}",
        "cohen_d": cd,
        "delta_d": cd - raw_d,
    })

# All cell types together
X = m_dm[celltypes].values
y = m_dm['panel_8g_z'].values
y_resid = residualize(y, X)
d1 = y_resid[m_dm['v17_dark_cluster'].values=='DM1']
d2 = y_resid[m_dm['v17_dark_cluster'].values=='DM2']
p = np.sqrt(((len(d1)-1)*d1.var(ddof=1) + (len(d2)-1)*d2.var(ddof=1)) / (len(d1)+len(d2)-2))
cd_all = (d1.mean()-d2.mean())/p
residual_results.append({
    "covariate": "residualize on ALL 8 cell-type fractions",
    "cohen_d": cd_all,
    "delta_d": cd_all - raw_d,
})

# Immune-only (T/Myeloid/B/NK)
immune_cts = [c for c in celltypes if c in ['T cell','Myeloid cell','B cell','NK cell']]
X = m_dm[immune_cts].values
y_resid = residualize(y, X)
d1 = y_resid[m_dm['v17_dark_cluster'].values=='DM1']
d2 = y_resid[m_dm['v17_dark_cluster'].values=='DM2']
p = np.sqrt(((len(d1)-1)*d1.var(ddof=1) + (len(d2)-1)*d2.var(ddof=1)) / (len(d1)+len(d2)-2))
cd_imm = (d1.mean()-d2.mean())/p
residual_results.append({
    "covariate": f"residualize on immune ({', '.join(immune_cts)})",
    "cohen_d": cd_imm,
    "delta_d": cd_imm - raw_d,
})

# Stromal
stromal_cts = [c for c in celltypes if c in ['Fibroblast','Endothelial cell']]
if stromal_cts:
    X = m_dm[stromal_cts].values
    y_resid = residualize(y, X)
    d1 = y_resid[m_dm['v17_dark_cluster'].values=='DM1']
    d2 = y_resid[m_dm['v17_dark_cluster'].values=='DM2']
    p = np.sqrt(((len(d1)-1)*d1.var(ddof=1) + (len(d2)-1)*d2.var(ddof=1)) / (len(d1)+len(d2)-2))
    cd_str = (d1.mean()-d2.mean())/p
    residual_results.append({
        "covariate": f"residualize on stromal ({', '.join(stromal_cts)})",
        "cohen_d": cd_str,
        "delta_d": cd_str - raw_d,
    })

resid_df = pd.DataFrame(residual_results)
resid_df.to_csv(OUT_DIR / "residualized_8gene_dm1_dm2.tsv", sep="\t", index=False)
print()
print("8-gene panel residualization:")
print(resid_df.round(3).to_string(index=False))

print()
print("=" * 70)
print("DONE")
print("=" * 70)
print(f"Outputs in: {OUT_DIR}")
print("Files: per_sample_celltype_fractions.tsv, dm1_vs_dm2_celltype_d.tsv, residualized_8gene_dm1_dm2.tsv")
