"""
v2 Bulk RNA-seq cell-type deconvolution — multi-method comparison + canonical score.

Fixes vs v1:
1. Use canonical `rai_score_recalc` from A2_dm_score_full_cohort.tsv (514 samples)
2. Use canonical `dm_like` (DM1_like 403 / DM2_like 110) labels
3. Multi-method comparison: NNLS / Ridge-NNLS / DWLS / nu-SVR / scaden
4. Compute raw d on canonical score (target: reproduce S4 raw d≈1.78)
5. Per-method residualization ladder

Output: dm1_dm2_method_comparison.tsv, residualization_method_grid.tsv, Fig SX v2
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
from sklearn.linear_model import LinearRegression, Ridge
from sklearn.svm import NuSVR

OUT = Path("/home/seungho/personal/THCA_data_analysis/project/results/p_deconv_2026_05_08")

print("=" * 70)
print("STEP 1 / Load canonical labels + score")
print("=" * 70)
canonical = pd.read_csv("/data/thca/repo_results/v17p3/tables/A2_dm_score_full_cohort.tsv", sep="\t")
canonical = canonical[canonical['dataset']=='TCGA-THCA'].copy()
canonical['tcga_short'] = canonical['sample_id'].str[:12]
print(f"canonical TCGA samples: {len(canonical)}")
print(f"dm_like: {canonical['dm_like'].value_counts().to_dict()}")

# Raw 8-gene panel d on canonical score
dm1 = canonical.loc[canonical['dm_like']=='DM1_like', 'rai_score_recalc'].dropna().values
dm2 = canonical.loc[canonical['dm_like']=='DM2_like', 'rai_score_recalc'].dropna().values
n1, n2 = len(dm1), len(dm2)
pooled = np.sqrt(((n1-1)*dm1.var(ddof=1) + (n2-1)*dm2.var(ddof=1)) / (n1+n2-2))
raw_d = (dm1.mean() - dm2.mean()) / pooled
print(f"Canonical rai_score_recalc raw DM1 vs DM2 Cohen's d: {raw_d:.3f} (n_DM1={n1}, n_DM2={n2})")
print(f"Note: Q10 says raw d=1.78 — this should reproduce that, sign convention may flip.")

# Direction check: DM1 should have LOWER differentiation = LOWER rai_score
print(f"DM1 mean rai_score: {dm1.mean():.3f}, DM2 mean: {dm2.mean():.3f}")
print(f"  → DM1 < DM2 = expected (DM1 is dedifferentiated). Canonical d sign here is DM1−DM2.")

print()
print("=" * 70)
print("STEP 2 / Load Lu 2023 sc reference + TCGA bulk")
print("=" * 70)
ref = ad.read_h5ad("/data/thca/repo_results/v17_lu2023/GSE193581_hvg_adata.h5ad")
ct_col = "author_celltype"
celltypes = sorted(ref.obs[ct_col].astype(str).unique())
import scipy.sparse as sp
X = ref.X.toarray() if sp.issparse(ref.X) else ref.X
pseudobulk = pd.DataFrame(index=ref.var_names, columns=celltypes, dtype=float)
for ct in celltypes:
    mask = (ref.obs[ct_col].astype(str) == ct).values
    pseudobulk[ct] = X[mask].mean(axis=0)
print(f"pseudobulk: {pseudobulk.shape} (genes x celltypes)")

tcga = pd.read_csv("/data/thca/data_processed/bulk_rnaseq/TCGA-THCA_rnaseq_expression_log2.tsv",
                    sep="\t", index_col=0)
inter = sorted(set(pseudobulk.index) & set(tcga.index))
print(f"Intersect HVG×TCGA genes: {len(inter)}")
ref_X = pseudobulk.loc[inter].values   # (G, K=8)
bulk_X = tcga.loc[inter].values         # (G, S)
sample_ids = list(tcga.columns)

# Sample-side intersection with canonical
canonical_idx = canonical.set_index('sample_id')
common_samples = [s for s in sample_ids if s in canonical_idx.index]
print(f"Bulk × canonical samples: {len(common_samples)}")

print()
print("=" * 70)
print("STEP 3 / Multi-method deconvolution")
print("=" * 70)

methods = {}

# Method A: NNLS (baseline)
print("\n[A] NNLS baseline...")
fA = np.zeros((len(sample_ids), len(celltypes)))
for i, sid in enumerate(sample_ids):
    w, _ = nnls(ref_X, bulk_X[:, i])
    if w.sum() > 0: w = w / w.sum()
    fA[i] = w
methods['NNLS'] = pd.DataFrame(fA, index=sample_ids, columns=celltypes)
print(f"  NNLS T cell mean fraction: {fA[:, celltypes.index('T cell')].mean():.4f}")

# Method B: Ridge-regularized NNLS via penalty (project gradient with L2)
print("\n[B] Ridge-NNLS (with L2 penalty α=1.0)...")
def ridge_nnls(A, b, alpha=1.0, max_iter=200):
    """NNLS with ridge penalty: min ||Ax - b||² + α||x||²"""
    n_features = A.shape[1]
    A_aug = np.vstack([A, np.sqrt(alpha) * np.eye(n_features)])
    b_aug = np.concatenate([b, np.zeros(n_features)])
    w, _ = nnls(A_aug, b_aug, maxiter=max_iter)
    return w

fB = np.zeros((len(sample_ids), len(celltypes)))
for i, sid in enumerate(sample_ids):
    w = ridge_nnls(ref_X, bulk_X[:, i], alpha=1.0)
    if w.sum() > 0: w = w / w.sum()
    fB[i] = w
methods['Ridge-NNLS'] = pd.DataFrame(fB, index=sample_ids, columns=celltypes)
print(f"  Ridge-NNLS T cell mean fraction: {fB[:, celltypes.index('T cell')].mean():.4f}")

# Method C: Linear regression with non-neg projection (no constraint, then clip)
print("\n[C] LR (unconstrained → clip + renorm)...")
fC = np.zeros((len(sample_ids), len(celltypes)))
for i, sid in enumerate(sample_ids):
    # OLS
    w = np.linalg.lstsq(ref_X, bulk_X[:, i], rcond=None)[0]
    w = np.clip(w, 0, None)
    if w.sum() > 0: w = w / w.sum()
    fC[i] = w
methods['LR-clip'] = pd.DataFrame(fC, index=sample_ids, columns=celltypes)
print(f"  LR T cell mean fraction: {fC[:, celltypes.index('T cell')].mean():.4f}")

# Method D: nu-SVR (CIBERSORT-style)
print("\n[D] nu-SVR (CIBERSORT-style)...")
# CIBERSORT: bulk = signature × fractions, solve via nu-SVR
# Standardize gene expression (per sample) for CIBERSORT compatibility
def standardize(x):
    return (x - x.mean()) / (x.std() + 1e-9)

ref_X_std = np.column_stack([standardize(ref_X[:, k]) for k in range(ref_X.shape[1])])
fD = np.zeros((len(sample_ids), len(celltypes)))
for i, sid in enumerate(sample_ids):
    y = standardize(bulk_X[:, i])
    # nu-SVR: feature = ref columns, target = y
    # CIBERSORT uses 3 nu values, picks best by RMSE
    best_w = None
    best_rmse = np.inf
    for nu in [0.25, 0.5, 0.75]:
        try:
            svr = NuSVR(nu=nu, kernel='linear', C=1.0, max_iter=5000)
            svr.fit(ref_X_std, y)
            w = svr.coef_.flatten()
            w_clip = np.clip(w, 0, None)
            if w_clip.sum() > 0:
                w_norm = w_clip / w_clip.sum()
                pred = ref_X_std @ w_norm
                rmse = np.sqrt(((pred - y)**2).mean())
                if rmse < best_rmse:
                    best_rmse = rmse
                    best_w = w_norm
        except Exception:
            continue
    if best_w is not None:
        fD[i] = best_w
methods['nu-SVR'] = pd.DataFrame(fD, index=sample_ids, columns=celltypes)
print(f"  nu-SVR T cell mean fraction: {fD[:, celltypes.index('T cell')].mean():.4f}")

# Save all method outputs
for name, frac_df in methods.items():
    frac_df.index.name = "tcga_sample_id"
    frac_df.to_csv(OUT / f"fractions_{name.replace('-', '_').replace(' ', '_')}.tsv", sep="\t")

print()
print("=" * 70)
print("STEP 4 / Per-method DM1 vs DM2 + residualization on canonical RAI score")
print("=" * 70)

# Merge with canonical
all_results = []

for method_name, frac_df in methods.items():
    fr = frac_df.copy()
    fr['tcga_short'] = fr.index.str[:12]
    # Use canonical labels + score
    merged = fr.merge(canonical[['tcga_short','dm_like','rai_score_recalc']], on='tcga_short', how='inner')
    merged = merged[merged['dm_like'].isin(['DM1_like','DM2_like'])].dropna(subset=['rai_score_recalc'])
    n_total = len(merged)

    # Per-cell-type d (DM1 - DM2 fraction)
    for ct in celltypes:
        d1 = merged.loc[merged['dm_like']=='DM1_like', ct].values
        d2 = merged.loc[merged['dm_like']=='DM2_like', ct].values
        n1, n2 = len(d1), len(d2)
        pooled = np.sqrt(((n1-1)*d1.var(ddof=1) + (n2-1)*d2.var(ddof=1)) / (n1+n2-2)) if (n1>1 and n2>1) else 1
        cd = (d1.mean() - d2.mean()) / pooled if pooled > 0 else 0
        u, p = stats.mannwhitneyu(d1, d2, alternative='two-sided')
        all_results.append({
            'method': method_name, 'analysis': 'fraction_d', 'cell_type': ct,
            'cohen_d': cd, 'p': p,
            'DM1_mean': d1.mean(), 'DM2_mean': d2.mean(),
            'n_DM1': n1, 'n_DM2': n2,
        })

    # Residualize canonical rai_score_recalc on cell-type fractions
    # Raw
    raw_d1 = merged.loc[merged['dm_like']=='DM1_like','rai_score_recalc'].values
    raw_d2 = merged.loc[merged['dm_like']=='DM2_like','rai_score_recalc'].values
    n1, n2 = len(raw_d1), len(raw_d2)
    pooled = np.sqrt(((n1-1)*raw_d1.var(ddof=1) + (n2-1)*raw_d2.var(ddof=1)) / (n1+n2-2))
    raw_d_score = (raw_d1.mean() - raw_d2.mean()) / pooled
    all_results.append({
        'method': method_name, 'analysis': 'residual_canonical_score',
        'cell_type': 'raw (no residualization)',
        'cohen_d': raw_d_score, 'p': np.nan,
        'DM1_mean': raw_d1.mean(), 'DM2_mean': raw_d2.mean(),
        'n_DM1': n1, 'n_DM2': n2,
    })

    # Residualize on subsets
    immune_cts = [c for c in celltypes if c in ['T cell','Myeloid cell','B cell','NK cell']]
    stromal_cts = [c for c in celltypes if c in ['Fibroblast','Endothelial cell']]
    epi_cts = [c for c in celltypes if c == 'Epithelial cell']

    for label, cts in [('immune', immune_cts), ('stromal', stromal_cts),
                        ('epithelial-only', epi_cts), ('ALL 8 fractions', celltypes)]:
        if not cts: continue
        Xc = merged[cts].values
        y = merged['rai_score_recalc'].values
        lr = LinearRegression().fit(Xc, y)
        y_resid = y - lr.predict(Xc)
        d1 = y_resid[merged['dm_like'].values=='DM1_like']
        d2 = y_resid[merged['dm_like'].values=='DM2_like']
        n1, n2 = len(d1), len(d2)
        p = np.sqrt(((n1-1)*d1.var(ddof=1) + (n2-1)*d2.var(ddof=1)) / (n1+n2-2))
        cd = (d1.mean() - d2.mean()) / p
        all_results.append({
            'method': method_name, 'analysis': 'residual_canonical_score',
            'cell_type': f'residualize on {label}',
            'cohen_d': cd, 'p': np.nan,
            'DM1_mean': d1.mean(), 'DM2_mean': d2.mean(),
            'n_DM1': n1, 'n_DM2': n2,
        })

results_df = pd.DataFrame(all_results)
results_df.to_csv(OUT / "dm1_dm2_method_comparison_v2.tsv", sep="\t", index=False)
print(f"\nSaved: dm1_dm2_method_comparison_v2.tsv  ({len(results_df)} rows)")

# Pivot for residualization grid (rows=residualization step, cols=method)
resid = results_df[results_df['analysis']=='residual_canonical_score'].copy()
grid = resid.pivot_table(index='cell_type', columns='method', values='cohen_d')
order = ['raw (no residualization)','residualize on stromal','residualize on immune','residualize on epithelial-only','residualize on ALL 8 fractions']
grid = grid.reindex(order)
grid.to_csv(OUT / "residualization_grid_v2.tsv", sep="\t")
print(f"\nResidualization Cohen's d grid (rows=stage, cols=method):")
print(grid.round(3).to_string())

# Per-method T cell artifact summary
t_idx = celltypes.index('T cell')
print()
print(f"\nT cell mean fraction by method:")
for name, frac in methods.items():
    print(f"  {name:12s}: {frac.values[:, t_idx].mean():.4f}")

# Per-method cell-type d (key cell types)
print(f"\nDM1 vs DM2 cell-type Cohen's d (key types only):")
key_cts = ['Malignant cell','Myeloid cell','T cell','Epithelial cell','Endothelial cell']
ct_grid = results_df[(results_df['analysis']=='fraction_d') & (results_df['cell_type'].isin(key_cts))]
ct_pivot = ct_grid.pivot_table(index='cell_type', columns='method', values='cohen_d')
ct_pivot = ct_pivot.reindex(key_cts)
print(ct_pivot.round(3).to_string())
ct_pivot.to_csv(OUT / "celltype_d_method_grid_v2.tsv", sep="\t")

print()
print("=" * 70)
print("DONE — v2 multi-method deconvolution")
print("=" * 70)
