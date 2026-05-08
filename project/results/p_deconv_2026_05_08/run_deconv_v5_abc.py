"""v5 ABC only — driver class + methylation × cell-type + sub-A/B teaser."""
from __future__ import annotations
from pathlib import Path
import warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
from scipy import stats

OUT = Path("/home/seungho/personal/THCA_data_analysis/project/results/p_deconv_2026_05_08")

def cohen_d(a, b):
    a, b = np.asarray(a), np.asarray(b)
    n1, n2 = len(a), len(b)
    if n1<2 or n2<2 or (a.var()+b.var())==0: return 0
    pooled = np.sqrt(((n1-1)*a.var(ddof=1)+(n2-1)*b.var(ddof=1))/(n1+n2-2))
    return (a.mean()-b.mean())/pooled if pooled>0 else 0

def patient_id(s): return "-".join(s.split("-")[:3])

celltypes = ['B cell','Endothelial cell','Epithelial cell','Fibroblast','Malignant cell','Myeloid cell','NK cell','T cell']

# Load common
tcga_frac = pd.read_csv(OUT / "fractions_nu_SVR.tsv", sep="\t", index_col=0)
fc = pd.read_csv("/data/thca/repo_results/v17/tables/fusion_calls_per_sample.tsv", sep="\t")
canonical = pd.read_csv("/data/thca/repo_results/v17p3/tables/A2_dm_score_full_cohort.tsv", sep="\t")
canonical = canonical[canonical['dataset']=='TCGA-THCA'].copy()

tcga_frac['sample_id'] = tcga_frac.index

# A. Per-driver-class
print("=" * 70)
print("A. Per-driver-class TCGA cell-type composition")
print("=" * 70)
m = tcga_frac.merge(fc[['sample_id','v3_anchor_6class']], on='sample_id', how='inner')
m = m.merge(canonical[['sample_id','dm_like']], on='sample_id', how='inner')
print(f"merged: {len(m)} samples; v3_anchor_6class: {m['v3_anchor_6class'].value_counts().to_dict()}")

A_results = []
for ct in celltypes:
    row = {'cell_type': ct}
    for c in ['BRAF_V600E','RAS_mutant','other']:
        v = m.loc[m['v3_anchor_6class']==c, ct].values
        row[f'{c}_mean'] = v.mean()
        row[f'{c}_n'] = len(v)
    braf = m.loc[m['v3_anchor_6class']=='BRAF_V600E', ct].values
    ras = m.loc[m['v3_anchor_6class']=='RAS_mutant', ct].values
    other = m.loc[m['v3_anchor_6class']=='other', ct].values
    row['d_other_vs_BRAF'] = cohen_d(other, braf)
    row['d_other_vs_RAS'] = cohen_d(other, ras)
    row['d_RAS_vs_BRAF'] = cohen_d(ras, braf)
    A_results.append(row)
A_df = pd.DataFrame(A_results)
A_df.to_csv(OUT / "v5A_per_driver_class.tsv", sep="\t", index=False)
print("\nPer-driver-class cell-type fraction:")
print(A_df[['cell_type','BRAF_V600E_mean','RAS_mutant_mean','other_mean',
            'd_other_vs_BRAF','d_other_vs_RAS']].round(3).to_string(index=False))

# B. Methylation × cell-type
print()
print("=" * 70)
print("B. Methylation × cell-type fraction Spearman correlation")
print("=" * 70)
meth = pd.read_csv("/data/thca/repo_results/audit_2026_04_30/round5/r5_2_sample_methylation_8gene.tsv", sep="\t")
meth_genes = ['DIO1','FOXE1','NKX2-1','PAX8','SLC5A5','TG','TPO','TSHR']
print(f"Methylation: {meth.shape}")

tcga_frac['patient'] = [patient_id(s) for s in tcga_frac.index]
meth['patient'] = meth['sample_short']
m_meth = tcga_frac.merge(meth, on='patient', how='inner')
print(f"Bulk × meth merged: {len(m_meth)}")

B_results = []
for gene in meth_genes + ['mean_8g_beta']:
    for ct in celltypes:
        if m_meth[ct].sum() == 0:
            B_results.append({'gene':gene,'cell_type':ct,'spearman_r':0,'p':1,'n':len(m_meth)})
            continue
        r, p = stats.spearmanr(m_meth[gene], m_meth[ct])
        B_results.append({'gene':gene,'cell_type':ct,'spearman_r':r,'p':p,'n':len(m_meth)})

B_df = pd.DataFrame(B_results)
B_df.to_csv(OUT / "v5B_methylation_celltype_corr.tsv", sep="\t", index=False)
B_pivot = B_df.pivot_table(index='gene', columns='cell_type', values='spearman_r')
B_pivot.to_csv(OUT / "v5B_methylation_celltype_pivot.tsv", sep="\t")
print("\nMethylation × cell-type Spearman r:")
print(B_pivot.round(3).to_string())

# C. DM1 sub-A vs sub-B
print()
print("=" * 70)
print("C. DM1 sub-A vs sub-B cell-type fraction (Paper 2 teaser)")
print("=" * 70)
sub = pd.read_csv("/data/thca/repo_results/d6p7_dm1_subcluster/dm1_subcluster_labels.tsv", sep="\t")
sub.columns = ['sample_id','sub_cluster']
print(f"sub labels: {sub.shape}; {sub.sub_cluster.value_counts().to_dict()}")

m_sub = tcga_frac.merge(sub, on='sample_id', how='inner')
print(f"merged sub × frac: {len(m_sub)}")

C_results = []
for ct in celltypes:
    a = m_sub.loc[m_sub['sub_cluster']=='sub_A', ct].values
    b = m_sub.loc[m_sub['sub_cluster']=='sub_B', ct].values
    cd = cohen_d(a, b)
    if len(a)>1 and len(b)>1:
        _, p = stats.mannwhitneyu(a, b, alternative='two-sided')
    else: p=1
    C_results.append({'cell_type':ct,'sub_A_mean':a.mean(),'sub_B_mean':b.mean(),
                      'cohen_d_A_vs_B':cd,'p':p,'n_A':len(a),'n_B':len(b)})
C_df = pd.DataFrame(C_results).sort_values('cohen_d_A_vs_B', key=abs, ascending=False)
C_df.to_csv(OUT / "v5C_dm1_subA_subB_celltype.tsv", sep="\t", index=False)
print("\nDM1 sub-A vs sub-B cell-type Cohen's d:")
print(C_df.round(4).to_string(index=False))

print()
print("=" * 70)
print("DONE — v5 ABC")
print("=" * 70)
