"""
v5 D — Pu 2021 full-transcriptome reference robustness (FINAL: NNLS + LinearSVR)

Replaces v5_d.py (nu-SVR libsvm too slow at 10K-33K rows) and v5_d_fast.py
(top-10K nu-SVR also slow). Final version uses LinearSVR (liblinear) for the
3rd method — fast enough to converge in seconds for 10K rows × 8 features.

Strategy:
  1. Pu 2021 sparse → Lu 2023 label transfer → 5K-cell pseudobulk full ~33K genes
  2. TCGA bulk re-deconv with:
     a. NNLS over Pu full transcriptome (~21K genes after intersect)
     b. LinearSVR over Pu top-10K variance genes (epsilon-insensitive linear)
  3. Compare DM d to Lu HVG nu-SVR (v2 primary) for 3-way concordance.
"""
from __future__ import annotations
from pathlib import Path
import warnings, gc, time
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp
from scipy.optimize import nnls
from scipy import stats
from sklearn.svm import LinearSVR

OUT = Path("/home/seungho/personal/THCA_data_analysis/project/results/p_deconv_2026_05_08")

def standardize(x):
    s = x.std()
    return (x - x.mean()) / (s + 1e-9)

def cohen_d(a, b):
    a, b = np.asarray(a), np.asarray(b)
    n1, n2 = len(a), len(b)
    if n1<2 or n2<2 or (a.var()+b.var())==0: return 0
    pooled = np.sqrt(((n1-1)*a.var(ddof=1)+(n2-1)*b.var(ddof=1))/(n1+n2-2))
    return (a.mean()-b.mean())/pooled if pooled>0 else 0

def nnls_deconv(ref_X, bulk_X, sample_ids, celltypes):
    fr = np.zeros((len(sample_ids), len(celltypes)))
    for i in range(bulk_X.shape[1]):
        w, _ = nnls(ref_X, bulk_X[:, i])
        if w.sum() > 0: w = w / w.sum()
        fr[i] = w
    return pd.DataFrame(fr, index=sample_ids, columns=celltypes)

def linear_svr_deconv(ref_X, bulk_X, sample_ids, celltypes):
    """LinearSVR (liblinear) - fast linear ε-SVR for CIBERSORT-style robustness."""
    ref_X_std = np.column_stack([standardize(ref_X[:, k]) for k in range(ref_X.shape[1])])
    fr = np.zeros((len(sample_ids), len(celltypes)))
    for i in range(len(sample_ids)):
        y = standardize(bulk_X[:, i])
        best_w, best_rmse = None, np.inf
        for eps in [0.001, 0.01, 0.1]:
            try:
                svr = LinearSVR(epsilon=eps, C=1.0, max_iter=2000, dual='auto')
                svr.fit(ref_X_std, y)
                w = np.clip(svr.coef_.flatten(), 0, None)
                if w.sum()>0:
                    w = w/w.sum()
                    pred = ref_X_std @ w
                    rmse = np.sqrt(((pred-y)**2).mean())
                    if rmse<best_rmse: best_rmse, best_w = rmse, w
            except Exception:
                continue
        if best_w is not None: fr[i] = best_w
    return pd.DataFrame(fr, index=sample_ids, columns=celltypes)


print("=" * 70)
print("D final — Pu 2021 → label transfer → full-transcript NNLS + top-10K LinearSVR")
print("=" * 70)

t0 = time.time()
pu_backed = ad.read_h5ad("/data/thca/scrna/raw/scrna_raw.h5ad", backed='r')
print(f"Pu 2021: {pu_backed.shape}")

np.random.seed(42)
n_target = 5000
pat_ids = np.array(pu_backed.obs.patient_id)
unique_pat = np.unique(pat_ids)
per_pat = max(1, n_target // len(unique_pat))
sub_ix = []
for pat in unique_pat:
    pat_idx = np.where(pat_ids == pat)[0]
    take = np.random.choice(pat_idx, min(len(pat_idx), per_pat * 2), replace=False)
    sub_ix.extend(take)
sub_ix = np.array(sorted(np.random.choice(sub_ix, min(n_target, len(sub_ix)), replace=False)))
print(f"Subsampled {len(sub_ix)} cells")

pu = pu_backed[sub_ix, :].to_memory()
del pu_backed; gc.collect()

if 'gene_symbol' in pu.var.columns:
    new_names = pu.var['gene_symbol'].astype(str).fillna('').values
    ensembl_ids = list(pu.var.index)
    final_names = [sym if sym and sym != 'nan' else ens for sym, ens in zip(new_names, ensembl_ids)]
    pu.var.index = pd.Index(final_names)
    pu.var_names_make_unique()
print(f"Pu in memory: {pu.shape}, t={time.time()-t0:.1f}s")

# === Lu HVG label transfer ===
lu = ad.read_h5ad("/data/thca/repo_results/v17_lu2023/GSE193581_hvg_adata.h5ad")
ct_col = "author_celltype"
celltypes_lu = sorted(lu.obs[ct_col].astype(str).unique())
Lu_X = lu.X.toarray() if sp.issparse(lu.X) else lu.X
Lu_pseudo = pd.DataFrame(index=lu.var_names, columns=celltypes_lu, dtype=float)
for ct in celltypes_lu:
    mask = (lu.obs[ct_col].astype(str)==ct).values
    Lu_pseudo[ct] = Lu_X[mask].mean(axis=0)

inter_pu_lu = sorted(set(lu.var_names) & set(pu.var_names))
print(f"Pu × Lu intersect: {len(inter_pu_lu)}")
lu_pseudo_inter = Lu_pseudo.loc[inter_pu_lu].values

pu_var_idx = {g: i for i, g in enumerate(pu.var_names)}
pu_idx_inter = [pu_var_idx[g] for g in inter_pu_lu]

if sp.issparse(pu.X):
    pu_inter = pu.X[:, pu_idx_inter].toarray()
else:
    pu_inter = pu.X[:, pu_idx_inter]
pu_lib = pu_inter.sum(axis=1, keepdims=True) + 1
pu_logn = np.log1p(pu_inter / pu_lib * 10000)
del pu_inter
pu_z = (pu_logn - pu_logn.mean(axis=0, keepdims=True)) / (pu_logn.std(axis=0, keepdims=True) + 1e-9)
lu_z = (lu_pseudo_inter - lu_pseudo_inter.mean(axis=0, keepdims=True)) / (lu_pseudo_inter.std(axis=0, keepdims=True) + 1e-9)
del pu_logn, lu_pseudo_inter; gc.collect()

pu_norm = pu_z / (np.linalg.norm(pu_z, axis=1, keepdims=True) + 1e-9)
lu_norm = lu_z.T / (np.linalg.norm(lu_z.T, axis=1, keepdims=True) + 1e-9)
sim = pu_norm @ lu_norm.T
del pu_z, lu_z, pu_norm, lu_norm; gc.collect()

pu_assigned = np.array([celltypes_lu[i] for i in sim.argmax(axis=1)])
print(f"label transfer counts: {dict(pd.Series(pu_assigned).value_counts())}")
print(f"label transfer t={time.time()-t0:.1f}s")

# === Pu full pseudobulk ===
print(f"\nBuilding Pu full pseudobulk...")
pu_X_csr = pu.X.tocsr().astype(np.float32) if sp.issparse(pu.X) else sp.csr_matrix(pu.X)
lib_full = np.array(pu_X_csr.sum(axis=1)).flatten() + 1
pseudo_data = {ct: np.zeros(pu.n_vars, dtype=np.float64) for ct in celltypes_lu}
counts = {ct: 0 for ct in celltypes_lu}
for i in range(pu.n_obs):
    row = pu_X_csr[i].toarray().flatten().astype(np.float64)
    row_norm = np.log1p(row / lib_full[i] * 10000)
    ct = pu_assigned[i]
    pseudo_data[ct] += row_norm
    counts[ct] += 1
Pu_pseudobulk_full = pd.DataFrame(index=pu.var_names, dtype=float)
for ct in celltypes_lu:
    if counts[ct] > 0:
        Pu_pseudobulk_full[ct] = pseudo_data[ct] / counts[ct]
    else:
        Pu_pseudobulk_full[ct] = 0
print(f"Pu pseudobulk: {Pu_pseudobulk_full.shape}, t={time.time()-t0:.1f}s")

# === TCGA bulk ===
tcga = pd.read_csv("/data/thca/data_processed/bulk_rnaseq/TCGA-THCA_rnaseq_expression_log2.tsv",
                   sep="\t", index_col=0)
inter_pu_tcga = sorted(set(Pu_pseudobulk_full.index) & set(tcga.index))
print(f"\nPu × TCGA gene intersect: {len(inter_pu_tcga)}")

ref_full = Pu_pseudobulk_full.loc[inter_pu_tcga].values
nonzero_cols = (ref_full.sum(axis=0) > 0)
celltypes_kept = [c for c, k in zip(celltypes_lu, nonzero_cols) if k]
ref_full_v = ref_full[:, nonzero_cols]
bulk_full = tcga.loc[inter_pu_tcga].values
sample_ids_tcga = list(tcga.columns)
print(f"non-empty cell types ({len(celltypes_kept)}): {celltypes_kept}")
print(f"ref shape: {ref_full_v.shape}, bulk shape: {bulk_full.shape}")

# === D1: NNLS full ===
print(f"\n[D1] NNLS Pu full ({len(inter_pu_tcga)} genes)...")
t1 = time.time()
tcga_full_nnls = nnls_deconv(ref_full_v, bulk_full, sample_ids_tcga, celltypes_kept)
print(f"NNLS done in {time.time()-t1:.1f}s")
tcga_full_nnls.to_csv(OUT / "fractions_TCGA_full_Pu_NNLS.tsv", sep="\t")

# === D2: LinearSVR top-10K ===
print(f"\n[D2] LinearSVR Pu top-10K variance genes...")
t2 = time.time()
gene_var = ref_full_v.var(axis=1)
top10k_idx = np.argsort(-gene_var)[:10000]
ref_10k = ref_full_v[top10k_idx, :]
bulk_10k = bulk_full[top10k_idx, :]
print(f"top-10K subset: ref {ref_10k.shape}, bulk {bulk_10k.shape}")
tcga_10k_lsvr = linear_svr_deconv(ref_10k, bulk_10k, sample_ids_tcga, celltypes_kept)
print(f"LinearSVR top-10K done in {time.time()-t2:.1f}s")
tcga_10k_lsvr.to_csv(OUT / "fractions_TCGA_top10k_Pu_LinearSVR.tsv", sep="\t")
print("Mean fractions (LinearSVR top-10K):")
print(tcga_10k_lsvr.mean().round(4))

# === DM d for both ===
canonical = pd.read_csv("/data/thca/repo_results/v17p3/tables/A2_dm_score_full_cohort.tsv", sep="\t")
canonical = canonical[canonical['dataset']=='TCGA-THCA'].copy()

def dm_d(frac_df, label):
    frac_df = frac_df.copy()
    frac_df['sample_id'] = frac_df.index
    m = frac_df.merge(canonical[['sample_id','dm_like']], on='sample_id')
    m = m[m['dm_like'].isin(['DM1_like','DM2_like'])]
    n1 = (m['dm_like']=='DM1_like').sum()
    n2 = (m['dm_like']=='DM2_like').sum()
    print(f"\n{label} DM1 vs DM2 (n_DM1={n1}, n_DM2={n2}):")
    rows = []
    for ct in celltypes_kept:
        a = m.loc[m['dm_like']=='DM1_like', ct].values
        b = m.loc[m['dm_like']=='DM2_like', ct].values
        d = cohen_d(a, b)
        rows.append({'cell_type':ct,'cohen_d':d,'DM1_mean':a.mean(),'DM2_mean':b.mean()})
    df = pd.DataFrame(rows).sort_values('cohen_d', key=abs, ascending=False)
    print(df.round(4).to_string(index=False))
    return df

D_nnls = dm_d(tcga_full_nnls, "D1 NNLS Pu full")
D_nnls.to_csv(OUT / "v5D_dm1_dm2_full_pu_NNLS.tsv", sep="\t", index=False)

D_lsvr = dm_d(tcga_10k_lsvr, "D2 LinearSVR Pu top-10K")
D_lsvr.to_csv(OUT / "v5D_dm1_dm2_top10k_pu_LinearSVR.tsv", sep="\t", index=False)

# === 3-way concordance ===
v2 = pd.read_csv(OUT / "celltype_d_method_grid_v2.tsv", sep="\t")
v2_d = dict(zip(v2['cell_type'], v2['nu_SVR']))
nnls_d = dict(zip(D_nnls['cell_type'], D_nnls['cohen_d']))
lsvr_d = dict(zip(D_lsvr['cell_type'], D_lsvr['cohen_d']))

rows = []
for ct in celltypes_lu:
    rows.append({
        'cell_type': ct,
        'd_LuHVG_nuSVR_v2': v2_d.get(ct, 0),
        'd_PuFull_NNLS': nnls_d.get(ct, 0),
        'd_PuTop10K_LinearSVR': lsvr_d.get(ct, 0),
    })
concord = pd.DataFrame(rows)
concord['sign_match_3way'] = (
    (np.sign(concord['d_LuHVG_nuSVR_v2'])
     == np.sign(concord['d_PuFull_NNLS']))
    & (np.sign(concord['d_LuHVG_nuSVR_v2'])
       == np.sign(concord['d_PuTop10K_LinearSVR']))
)
concord['informative'] = (
    (concord['d_LuHVG_nuSVR_v2'].abs() > 0.05)
    & (concord['d_PuFull_NNLS'].abs() > 0.05)
    & (concord['d_PuTop10K_LinearSVR'].abs() > 0.05)
)
concord.to_csv(OUT / "v5D_three_way_concordance.tsv", sep="\t", index=False)
print("\n=== 3-way concordance: Lu HVG nu-SVR v2 / Pu full NNLS / Pu top-10K LinearSVR ===")
print(concord.round(3).to_string(index=False))

print(f"\nDONE — v5 D final (total t={time.time()-t0:.1f}s)")
