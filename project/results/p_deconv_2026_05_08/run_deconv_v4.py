"""
v4 Extensions:
  E3. Lee/GSE213647 binary DM1/DM2 split + direct effect-size comparison with TCGA
  E4. Lee residualization mirror (TCGA S4 ladder applied to Lee)
  E5. GSE76039 (Landa 2016) advanced disease (PDTC+ATC, n≈37) deconvolution
  E6. scaden 5th method retry (gene-namespace fix)
"""
from __future__ import annotations
from pathlib import Path
import warnings
warnings.filterwarnings("ignore")
import os, subprocess

import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp
from scipy import stats
from sklearn.svm import NuSVR
from sklearn.linear_model import LinearRegression

OUT = Path("/home/seungho/personal/THCA_data_analysis/project/results/p_deconv_2026_05_08")

# ---- Reference (Lu 2023, same as v2/v3) ----
ref = ad.read_h5ad("/data/thca/repo_results/v17_lu2023/GSE193581_hvg_adata.h5ad")
ct_col = "author_celltype"
celltypes = sorted(ref.obs[ct_col].astype(str).unique())
X = ref.X.toarray() if sp.issparse(ref.X) else ref.X
pseudobulk = pd.DataFrame(index=ref.var_names, columns=celltypes, dtype=float)
for ct in celltypes:
    mask = (ref.obs[ct_col].astype(str) == ct).values
    pseudobulk[ct] = X[mask].mean(axis=0)

def standardize(x): return (x - x.mean()) / (x.std() + 1e-9)

def nu_svr_deconv(ref_X, bulk_X, sample_ids, celltypes):
    ref_X_std = np.column_stack([standardize(ref_X[:, k]) for k in range(ref_X.shape[1])])
    fr = np.zeros((len(sample_ids), len(celltypes)))
    for i in range(len(sample_ids)):
        y = standardize(bulk_X[:, i])
        best_w, best_rmse = None, np.inf
        for nu in [0.25, 0.5, 0.75]:
            try:
                svr = NuSVR(nu=nu, kernel='linear', C=1.0, max_iter=5000)
                svr.fit(ref_X_std, y)
                w = np.clip(svr.coef_.flatten(), 0, None)
                if w.sum() > 0:
                    w = w / w.sum()
                    pred = ref_X_std @ w
                    rmse = np.sqrt(((pred - y)**2).mean())
                    if rmse < best_rmse:
                        best_rmse, best_w = rmse, w
            except Exception:
                continue
        if best_w is not None: fr[i] = best_w
    return pd.DataFrame(fr, index=sample_ids, columns=celltypes)

def cohen_d(a, b):
    a, b = np.asarray(a), np.asarray(b)
    n1, n2 = len(a), len(b)
    if n1 < 2 or n2 < 2 or (a.var()+b.var())==0: return 0
    pooled = np.sqrt(((n1-1)*a.var(ddof=1) + (n2-1)*b.var(ddof=1)) / (n1+n2-2))
    return (a.mean() - b.mean()) / pooled if pooled > 0 else 0

# Load Ensembl-symbol mapping
gmap = pd.read_csv("/data/thca/repo_results/v17p3/tables/F1_gene_recovery_mapping.tsv", sep="\t")
ens_to_sym = dict(zip(gmap.ensembl, gmap.symbol))


print("=" * 70)
print("E3 + E4. Lee/GSE213647 binary DM1/DM2 split + residualization mirror")
print("=" * 70)
lee_frac = pd.read_csv(OUT / "fractions_LEE_nu_SVR.tsv", sep="\t", index_col=0)
lee_panel = pd.read_csv("/data/thca/repo_results/v17_korean/GSE213647_panel_score.tsv", sep="\t")
m_lee = lee_frac.copy()
m_lee['gsm'] = m_lee.index
m_lee = m_lee.merge(lee_panel[['gsm','panel_z','tissue_type','histology']], on='gsm')

# Tumor only
tumor_lee = m_lee[m_lee['tissue_type']!='Normal'].copy()
print(f"Lee tumor-only: {len(tumor_lee)}")

# Binary DM1/DM2 split — match TCGA proportion (DM1 ~78%, DM2 ~22%)
# So use the 22nd percentile of panel_z as threshold (LOW panel_z = DM1-like)
# Actually since panel_z is centered, simpler: top 22% as DM2 (high panel_z = differentiated = DM2)
threshold = np.quantile(tumor_lee['panel_z'], 0.78)  # 22% above this = DM2_like
tumor_lee['lee_dm_like'] = np.where(tumor_lee['panel_z'] >= threshold, 'DM2_like', 'DM1_like')
n_dm1 = (tumor_lee['lee_dm_like']=='DM1_like').sum()
n_dm2 = (tumor_lee['lee_dm_like']=='DM2_like').sum()
print(f"Lee tumor binary split (78/22 proportion-matched to TCGA):")
print(f"  threshold panel_z = {threshold:.3f}")
print(f"  DM1_like {n_dm1}, DM2_like {n_dm2}")

# Per cell-type DM1 vs DM2 d in Lee
lee_d = []
for ct in celltypes:
    d1 = tumor_lee.loc[tumor_lee['lee_dm_like']=='DM1_like', ct].values
    d2 = tumor_lee.loc[tumor_lee['lee_dm_like']=='DM2_like', ct].values
    cd = cohen_d(d1, d2)
    if len(d1)>1 and len(d2)>1:
        _, p = stats.mannwhitneyu(d1, d2, alternative='two-sided')
    else:
        p = 1
    lee_d.append({'cell_type':ct,'cohen_d':cd,'p':p,
                  'DM1_mean':d1.mean(),'DM2_mean':d2.mean(),'n_DM1':len(d1),'n_DM2':len(d2)})
lee_d_df = pd.DataFrame(lee_d).sort_values('cohen_d', key=abs, ascending=False)
lee_d_df.to_csv(OUT / "lee_celltype_dm1_dm2_binary.tsv", sep="\t", index=False)
print("\nLee DM1 vs DM2 cell-type Cohen's d (binary, 78/22 split):")
print(lee_d_df.round(4).to_string(index=False))

# E4. Lee residualization ladder (mirror of TCGA S4)
y = tumor_lee['panel_z'].values
raw_d_lee = cohen_d(y[tumor_lee['lee_dm_like'].values=='DM1_like'],
                    y[tumor_lee['lee_dm_like'].values=='DM2_like'])
print(f"\nLee panel_z raw Cohen's d (DM1−DM2 binary): {raw_d_lee:.3f}")

resid_lee = [{"covariate": "raw (no residualization)", "cohen_d": raw_d_lee, "delta_d": 0}]
immune = [c for c in celltypes if c in ['T cell','Myeloid cell','B cell','NK cell']]
stromal = [c for c in celltypes if c in ['Fibroblast','Endothelial cell']]
epi = [c for c in celltypes if c=='Epithelial cell']
for label, cts in [('stromal', stromal), ('immune', immune),
                    ('epithelial-only', epi), ('ALL 8 fractions', celltypes)]:
    if not cts: continue
    Xc = tumor_lee[cts].values
    lr = LinearRegression().fit(Xc, y)
    y_resid = y - lr.predict(Xc)
    cd_r = cohen_d(y_resid[tumor_lee['lee_dm_like'].values=='DM1_like'],
                   y_resid[tumor_lee['lee_dm_like'].values=='DM2_like'])
    resid_lee.append({"covariate": f"residualize on {label}", "cohen_d": cd_r,
                      "delta_d": cd_r - raw_d_lee})
resid_lee_df = pd.DataFrame(resid_lee)
resid_lee_df.to_csv(OUT / "lee_residualization_ladder.tsv", sep="\t", index=False)
print(f"\nLee residualization ladder (panel_z DM1 vs DM2):")
print(resid_lee_df.round(3).to_string(index=False))

# Compute retention
final_d = abs(resid_lee_df.iloc[-1]['cohen_d'])
raw_d = abs(resid_lee_df.iloc[0]['cohen_d'])
retention = final_d / raw_d * 100
print(f"\nLee retention after ALL-8 residualization: {retention:.1f}% (TCGA nu-SVR: 47%, S4 canonical: 56%)")


print()
print("=" * 70)
print("E5. GSE76039 (Landa 2016) advanced disease (PDTC+ATC, n=37) deconvolution")
print("=" * 70)
g76 = pd.read_csv("/data/thca/data_processed/microarray/GSE76039_microarray_expression_log2.tsv",
                   sep="\t", index_col=0)
print(f"GSE76039 microarray shape: {g76.shape}")
inter76 = sorted(set(pseudobulk.index) & set(g76.index))
print(f"GSE76039 × Lu HVG intersect: {len(inter76)}")

ref76 = pseudobulk.loc[inter76].values
bulk76 = g76.loc[inter76].values
sids76 = list(g76.columns)
print(f"Running nu-SVR on GSE76039 n={len(sids76)}...")
g76_frac = nu_svr_deconv(ref76, bulk76, sids76, celltypes)
print("Mean GSE76039 fractions (PDTC+ATC):")
print(g76_frac.mean().round(4))
g76_frac.to_csv(OUT / "fractions_GSE76039_nu_SVR.tsv", sep="\t")

# Compare GSE76039 (advanced) vs Lee tumor (PTC dominant, panel_z mid) vs Lee normal
print("\nMean cell-type fractions across cohorts:")
mean_compare = pd.DataFrame({
    'GSE76039 (advanced PDTC+ATC, n=37)': g76_frac.mean(),
    'Lee DM1_like (n=' + str(n_dm1) + ')': tumor_lee.loc[tumor_lee['lee_dm_like']=='DM1_like', celltypes].mean(),
    'Lee DM2_like (n=' + str(n_dm2) + ')': tumor_lee.loc[tumor_lee['lee_dm_like']=='DM2_like', celltypes].mean(),
    'Lee Normal (n=263)': m_lee.loc[m_lee['tissue_type']=='Normal', celltypes].mean(),
})
print(mean_compare.round(4).to_string())
mean_compare.to_csv(OUT / "cohort_celltype_mean_compare.tsv", sep="\t")

# GSE76039 vs Lee Normal Cohen's d (advanced disease vs normal — biggest contrast)
g76_vs_norm = []
norm_subset = m_lee.loc[m_lee['tissue_type']=='Normal', celltypes]
for ct in celltypes:
    cd = cohen_d(g76_frac[ct].values, norm_subset[ct].values)
    g76_vs_norm.append({'cell_type':ct,'cohen_d':cd,'GSE76039_mean':g76_frac[ct].mean(),
                        'Lee_normal_mean':norm_subset[ct].mean(),
                        'n_advanced':len(g76_frac),'n_normal':len(norm_subset)})
g76_norm_df = pd.DataFrame(g76_vs_norm).sort_values('cohen_d', key=abs, ascending=False)
g76_norm_df.to_csv(OUT / "gse76039_vs_lee_normal_celltype.tsv", sep="\t", index=False)
print("\nGSE76039 (advanced) vs Lee normal cell-type Cohen's d:")
print(g76_norm_df.round(4).to_string(index=False))


print()
print("=" * 70)
print("E6. scaden 5th method — gene-namespace fix retry")
print("=" * 70)
SCADEN = Path("/tmp/scaden_v4")
SCADEN.mkdir(exist_ok=True, parents=True)

# Subsample Lu 2023 + use HVG-restricted gene namespace = exactly what TCGA bulk has
import numpy as np
np.random.seed(42)
from collections import defaultdict
ref_sub = ad.read_h5ad("/data/thca/repo_results/v17_lu2023/GSE193581_hvg_adata.h5ad")
idxs = defaultdict(list)
for i, ct in enumerate(ref_sub.obs[ct_col]):
    idxs[ct].append(i)
sub = []
for ct, ids in idxs.items():
    sub.extend(np.random.choice(ids, min(len(ids), 1500), replace=False))
sub = sorted(sub)
X_sub = ref_sub.X.toarray() if sp.issparse(ref_sub.X) else ref_sub.X
X_sub = X_sub[sub]
X_sub = X_sub - X_sub.min()
# Use only genes that intersect with TCGA bulk
tcga = pd.read_csv("/data/thca/data_processed/bulk_rnaseq/TCGA-THCA_rnaseq_expression_log2.tsv",
                    sep="\t", index_col=0)
inter_tg = sorted(set(ref_sub.var_names) & set(tcga.index))
print(f"Lu HVG ∩ TCGA gene_symbol: {len(inter_tg)}")
gene_idx = [list(ref_sub.var_names).index(g) for g in inter_tg]
X_sub_g = X_sub[:, gene_idx]
print(f"sc subsampled+gene-restricted: {X_sub_g.shape}")

cells_df = pd.DataFrame(X_sub_g, columns=inter_tg)
ct_df = pd.DataFrame({'Celltype': ref_sub.obs[ct_col].values[sub]})
cells_df.to_csv(SCADEN / "lu2023_counts.txt", sep="\t", index=False)
ct_df.to_csv(SCADEN / "lu2023_celltypes.txt", sep="\t", index=False)
print(f"counts.txt: {cells_df.shape}, celltypes.txt: {ct_df.shape}")

# TCGA bulk in same gene space, samples × genes, non-negative
tcga_T = tcga.loc[inter_tg].T
tcga_T = tcga_T - tcga_T.values.min()
tcga_T.to_csv(SCADEN / "tcga_bulk.txt", sep="\t", index=True)
print(f"TCGA bulk for scaden: {tcga_T.shape}")

# scaden simulate
print("\nscaden simulate (300 samples, 100 cells each)...")
(SCADEN / "sim").mkdir(exist_ok=True)
res = subprocess.run(["scaden", "simulate", "-d", str(SCADEN), "-n", "300", "-c", "100",
                      "--pattern", "*_counts.txt", "--prefix", "lu2023",
                      "-o", str(SCADEN / "sim")],
                     capture_output=True, text=True, timeout=600)
sim_h5ad = SCADEN / "sim" / "lu2023.h5ad"
print(f"sim file exists: {sim_h5ad.exists()}, size: {sim_h5ad.stat().st_size if sim_h5ad.exists() else 0}")
if not sim_h5ad.exists():
    print("STDERR:", res.stderr[-1000:])
    print("scaden simulate FAILED")
else:
    # scaden process
    print("\nscaden process...")
    res = subprocess.run(["scaden", "process", str(sim_h5ad), str(SCADEN/"tcga_bulk.txt"),
                          "--processed_path", str(SCADEN/"sim"/"processed.h5ad")],
                         capture_output=True, text=True, timeout=300)
    proc_path = SCADEN/"sim"/"processed.h5ad"
    print(f"processed exists: {proc_path.exists()}")
    if not proc_path.exists():
        print("STDERR:", res.stderr[-1500:])
        print("scaden process FAILED — defer to future")
    else:
        # scaden train (very few epochs for speed)
        print("\nscaden train (5000 steps, fast)...")
        (SCADEN/"models").mkdir(exist_ok=True)
        res = subprocess.run(["scaden", "train", str(proc_path), "--steps", "5000",
                              "--model_dir", str(SCADEN/"models")],
                             capture_output=True, text=True, timeout=1800)
        print(f"train return code: {res.returncode}")
        # Predict
        print("\nscaden predict on TCGA...")
        res = subprocess.run(["scaden", "predict", str(SCADEN/"tcga_bulk.txt"),
                              "--model_dir", str(SCADEN/"models"),
                              "--outname", str(OUT / "scaden_predictions.tsv")],
                             capture_output=True, text=True, timeout=300)
        if (OUT / "scaden_predictions.tsv").exists():
            scaden_pred = pd.read_csv(OUT / "scaden_predictions.tsv", sep="\t", index_col=0)
            print(f"scaden predictions: {scaden_pred.shape}")
            print(f"Mean scaden fractions:")
            print(scaden_pred.mean().round(4))
        else:
            print("scaden predict FAILED:", res.stderr[-800:])

print()
print("=" * 70)
print("DONE — v4")
print("=" * 70)
