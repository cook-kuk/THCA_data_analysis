"""
v3 Extensions (fixed):
  E1. Lee/GSE213647 (n=632) nu-SVR — Ensembl→symbol mapping fix
  E2. TCGA within-DM1 fusion+ vs fusion− — aggregate cbio_sv_thca per-sample
"""
from __future__ import annotations
from pathlib import Path
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp
from scipy import stats
from sklearn.svm import NuSVR

OUT = Path("/home/seungho/personal/THCA_data_analysis/project/results/p_deconv_2026_05_08")

# ---- Load Lu 2023 reference (same as v2)
ref = ad.read_h5ad("/data/thca/repo_results/v17_lu2023/GSE193581_hvg_adata.h5ad")
ct_col = "author_celltype"
celltypes = sorted(ref.obs[ct_col].astype(str).unique())
X = ref.X.toarray() if sp.issparse(ref.X) else ref.X
pseudobulk = pd.DataFrame(index=ref.var_names, columns=celltypes, dtype=float)
for ct in celltypes:
    mask = (ref.obs[ct_col].astype(str) == ct).values
    pseudobulk[ct] = X[mask].mean(axis=0)

# Load Ensembl ↔ symbol mapping
gmap = pd.read_csv("/data/thca/repo_results/v17p3/tables/F1_gene_recovery_mapping.tsv", sep="\t")
print(f"Gene map: {gmap.shape} (Ensembl ↔ symbol)")
ens_to_sym = dict(zip(gmap.ensembl, gmap.symbol))

def standardize(x):
    return (x - x.mean()) / (x.std() + 1e-9)

def nu_svr_deconv(ref_X, bulk_X, sample_ids, celltypes):
    ref_X_std = np.column_stack([standardize(ref_X[:, k]) for k in range(ref_X.shape[1])])
    fr = np.zeros((len(sample_ids), len(celltypes)))
    for i, sid in enumerate(sample_ids):
        y = standardize(bulk_X[:, i])
        best_w, best_rmse = None, np.inf
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
            fr[i] = best_w
    return pd.DataFrame(fr, index=sample_ids, columns=celltypes)


print("=" * 70)
print("E1. Lee/GSE213647 cross-cohort deconvolution (Ensembl→symbol fix)")
print("=" * 70)
lee = pd.read_csv("/data/thca/data_processed/bulk_rnaseq/GSE213647_rnaseq_expression_log2.tsv",
                   sep="\t", index_col=0)
print(f"Lee bulk shape: {lee.shape} (Ensembl gene_id × samples)")

# Convert Lee Ensembl → symbol
lee_sym = lee.copy()
lee_sym['symbol'] = lee_sym.index.map(ens_to_sym)
lee_sym = lee_sym.dropna(subset=['symbol']).drop_duplicates('symbol').set_index('symbol')
print(f"Lee mapped to symbol: {lee_sym.shape}")

inter = sorted(set(pseudobulk.index) & set(lee_sym.index))
print(f"Lee-symbol × Lu HVG intersect: {len(inter)}")

ref_X = pseudobulk.loc[inter].values
bulk_X = lee_sym.loc[inter].values
sample_ids = list(lee_sym.columns)

print(f"ref_X shape: {ref_X.shape}, bulk_X shape: {bulk_X.shape}")
print(f"bulk_X stats: min={bulk_X.min():.2f}, max={bulk_X.max():.2f}, mean={bulk_X.mean():.2f}")

print("\nRunning nu-SVR on Lee...")
lee_frac = nu_svr_deconv(ref_X, bulk_X, sample_ids, celltypes)
print("Mean Lee fractions:")
print(lee_frac.mean().round(4))
lee_frac.to_csv(OUT / "fractions_LEE_nu_SVR.tsv", sep="\t")

# Lee panel_z join
lee_panel = pd.read_csv("/data/thca/repo_results/v17_korean/GSE213647_panel_score.tsv", sep="\t")
merged_lee = lee_frac.copy()
merged_lee['gsm'] = merged_lee.index
merged_lee = merged_lee.merge(lee_panel[['gsm','panel_z','tissue_type','histology']], on='gsm', how='inner')
print(f"\nmerged Lee: {len(merged_lee)} samples; tissue_type: {merged_lee.tissue_type.value_counts().to_dict()}")

tumor_lee = merged_lee[merged_lee['tissue_type']!='Normal'].copy()
print(f"Lee tumor-only (PTC+ATC+PDTC): {len(tumor_lee)}")

# Per-cell-type Spearman with panel_z (panel_z higher = MORE differentiated)
spearman_results = []
for ct in celltypes:
    if tumor_lee[ct].sum() == 0:
        spearman_results.append({'cohort':'Lee tumor-only','cell_type':ct,'spearman_r':0,'p':1,'n':len(tumor_lee)})
        continue
    r, p = stats.spearmanr(tumor_lee[ct], tumor_lee['panel_z'])
    spearman_results.append({'cohort':'Lee tumor-only','cell_type':ct,'spearman_r':r,'p':p,'n':len(tumor_lee)})

# Tumor + Normal
for ct in celltypes:
    if merged_lee[ct].sum() == 0:
        spearman_results.append({'cohort':'Lee tumor+normal','cell_type':ct,'spearman_r':0,'p':1,'n':len(merged_lee)})
        continue
    r, p = stats.spearmanr(merged_lee[ct], merged_lee['panel_z'])
    spearman_results.append({'cohort':'Lee tumor+normal','cell_type':ct,'spearman_r':r,'p':p,'n':len(merged_lee)})

lee_corr = pd.DataFrame(spearman_results)
lee_corr.to_csv(OUT / "lee_celltype_panelz_spearman.tsv", sep="\t", index=False)
print("\nLee per-cell-type Spearman r vs panel_z (tumor-only):")
print(lee_corr[lee_corr.cohort=='Lee tumor-only'].sort_values('spearman_r').round(4).to_string(index=False))

# Lee tumor vs normal
lee_tn = []
for ct in celltypes:
    t = merged_lee.loc[merged_lee['tissue_type']!='Normal', ct].values
    n = merged_lee.loc[merged_lee['tissue_type']=='Normal', ct].values
    if len(t)<2 or len(n)<2 or (t.var()+n.var())==0:
        lee_tn.append({'cell_type':ct,'tumor_mean':t.mean() if len(t) else 0,
                       'normal_mean':n.mean() if len(n) else 0,'cohen_d':0,'p':1,
                       'n_tumor':len(t),'n_normal':len(n)})
        continue
    pooled = np.sqrt(((len(t)-1)*t.var(ddof=1)+(len(n)-1)*n.var(ddof=1))/(len(t)+len(n)-2))
    cd = (t.mean()-n.mean())/pooled if pooled>0 else 0
    _, p = stats.mannwhitneyu(t, n, alternative='two-sided')
    lee_tn.append({'cell_type':ct,'tumor_mean':t.mean(),'normal_mean':n.mean(),
                   'cohen_d':cd,'p':p,'n_tumor':len(t),'n_normal':len(n)})
lee_tn_df = pd.DataFrame(lee_tn).sort_values('cohen_d', key=abs, ascending=False)
lee_tn_df.to_csv(OUT / "lee_tumor_vs_normal_celltype.tsv", sep="\t", index=False)
print("\nLee tumor vs normal cell-type Cohen's d:")
print(lee_tn_df.round(4).to_string(index=False))

print()
print("=" * 70)
print("E2. TCGA within-DM1 fusion+ vs fusion− cell-type fraction (cbio_sv aggregation)")
print("=" * 70)
# Aggregate cbio_sv per-sample fusion flag
cbio = pd.read_csv("/data/thca/repo_results/audit_2026_04_30/round3/cbio_sv_thca.tsv", sep="\t")
# Filter to actual fusion events (fusion connection or kinase fusion)
# Most rows are fusions; we just count any SV record per sample as fusion+
cbio_fusion_samples = set(cbio['sampleId'].unique())
print(f"cbio_sv total events: {len(cbio)}, unique samples with at least 1 event: {len(cbio_fusion_samples)}")
# But we want kinase fusion specifically — eventInfo or annotation
ki_keywords = ['RET','NTRK','ALK','BRAF','PAX8','PPARG']
kinase_fusion_mask = cbio['eventInfo'].fillna('').str.contains('|'.join(ki_keywords), case=False, regex=True)
cbio_kinase = cbio[kinase_fusion_mask]
kinase_fusion_samples = set(cbio_kinase['sampleId'].unique())
print(f"Kinase-fusion (RET/NTRK/ALK/BRAF/PAX8/PPARG) samples: {len(kinase_fusion_samples)}")

# Load existing nu-SVR TCGA fractions
tcga_frac = pd.read_csv(OUT / "fractions_nu_SVR.tsv", sep="\t", index_col=0)
canonical = pd.read_csv("/data/thca/repo_results/v17p3/tables/A2_dm_score_full_cohort.tsv", sep="\t")
canonical = canonical[canonical['dataset']=='TCGA-THCA'].copy()

# Match sampleId format. cbio sampleId: TCGA-DJ-A2Q6-01 (no letter)
# tcga_frac index:                       TCGA-DJ-A2Q6-01A (with sample-type letter)
# Strategy: normalize to 3-segment patient ID (TCGA-XX-XXXX) since each patient typically has 1 tumor sample
def patient_id(s):
    return "-".join(s.split("-")[:3])

tcga_frac['patient'] = [patient_id(s) for s in tcga_frac.index]
canonical['patient'] = [patient_id(s) for s in canonical['sample_id']]
fusion_kin_pt = set(patient_id(s) for s in kinase_fusion_samples)
fusion_any_pt = set(patient_id(s) for s in cbio_fusion_samples)
print(f"  patient-id alignment: cbio kinase {len(fusion_kin_pt)}, frac patients in fusion+ kin: {tcga_frac['patient'].isin(fusion_kin_pt).sum()}")

m1 = tcga_frac.merge(canonical[['patient','dm_like']], on='patient', how='inner')
m1['fusion_pos_kinase'] = m1['patient'].isin(fusion_kin_pt)
m1['fusion_pos_any'] = m1['patient'].isin(fusion_any_pt)
print(f"merged TCGA × DM × fusion: {len(m1)}")
print(f"  kinase fusion+: {m1['fusion_pos_kinase'].sum()}, kinase fusion−: {(~m1['fusion_pos_kinase']).sum()}")
print(f"  any fusion+:    {m1['fusion_pos_any'].sum()}, any fusion−:    {(~m1['fusion_pos_any']).sum()}")

dm1_only = m1[m1['dm_like']=='DM1_like'].copy()
print(f"\nDM1 samples: {len(dm1_only)}")
print(f"  DM1 kinase-fusion+: {dm1_only['fusion_pos_kinase'].sum()} | kinase-fusion−: {(~dm1_only['fusion_pos_kinase']).sum()}")
print(f"  DM1 any-fusion+:    {dm1_only['fusion_pos_any'].sum()} | any-fusion−:    {(~dm1_only['fusion_pos_any']).sum()}")

within_dm1 = []
for ct in celltypes:
    fp = dm1_only.loc[dm1_only['fusion_pos_kinase'], ct].values
    fn = dm1_only.loc[~dm1_only['fusion_pos_kinase'], ct].values
    if len(fp)<2 or len(fn)<2 or (fp.var()+fn.var())==0:
        within_dm1.append({'cell_type':ct,'fusion_pos_mean':fp.mean() if len(fp) else 0,
                           'fusion_neg_mean':fn.mean() if len(fn) else 0,
                           'cohen_d':0,'p':1,'n_fp':len(fp),'n_fn':len(fn)})
        continue
    pooled = np.sqrt(((len(fp)-1)*fp.var(ddof=1)+(len(fn)-1)*fn.var(ddof=1))/(len(fp)+len(fn)-2))
    cd = (fp.mean()-fn.mean())/pooled if pooled>0 else 0
    _, p = stats.mannwhitneyu(fp, fn, alternative='two-sided')
    within_dm1.append({'cell_type':ct,'fusion_pos_mean':fp.mean(),'fusion_neg_mean':fn.mean(),
                       'cohen_d':cd,'p':p,'n_fp':len(fp),'n_fn':len(fn)})

wd_df = pd.DataFrame(within_dm1).sort_values('cohen_d', key=abs, ascending=False)
wd_df.to_csv(OUT / "within_dm1_fusion_celltype.tsv", sep="\t", index=False)
print("\nWithin DM1 (n=403) fusion+ vs fusion− cell-type Cohen's d (kinase fusion):")
print(wd_df.round(4).to_string(index=False))

print()
print("=" * 70)
print("DONE — v3 cross-cohort + within-DM1 (fixed)")
print("=" * 70)
