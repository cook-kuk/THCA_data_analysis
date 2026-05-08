"""v6: per-fusion-partner deconv + sub-A/B × methylation interaction."""
from pathlib import Path
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from scipy import stats

OUT = Path("/home/seungho/personal/THCA_data_analysis/project/results/p_deconv_2026_05_08")

def cohen_d(a, b):
    a, b = np.asarray(a), np.asarray(b)
    n1, n2 = len(a), len(b)
    if n1<2 or n2<2 or (a.var()+b.var())==0: return 0
    pooled = np.sqrt(((n1-1)*a.var(ddof=1)+(n2-1)*b.var(ddof=1))/(n1+n2-2))
    return (a.mean()-b.mean())/pooled if pooled>0 else 0

# Load nu-SVR fractions + canonical + cbio_sv (events)
frac = pd.read_csv(OUT/"fractions_nu_SVR.tsv", sep="\t", index_col=0)
canonical = pd.read_csv("/data/thca/repo_results/v17p3/tables/A2_dm_score_full_cohort.tsv", sep="\t")
canonical = canonical[canonical['dataset']=='TCGA-THCA'].copy()
cbio = pd.read_csv("/data/thca/repo_results/audit_2026_04_30/round3/cbio_sv_thca.tsv", sep="\t")

def patient_id(s): return "-".join(s.split("-")[:3])
frac['patient'] = [patient_id(s) for s in frac.index]
canonical['patient'] = [patient_id(s) for s in canonical['sample_id']]

# Per-sample fusion partner classification
# Classes: RET / NTRK / ALK / BRAF / PAX8_PPARG / no_fusion
def classify_fusion(s):
    s = str(s) if pd.notna(s) else ''
    if 'RET' in s.upper(): return 'RET'
    if 'NTRK' in s.upper(): return 'NTRK'
    if 'ALK' in s.upper(): return 'ALK'
    if 'BRAF' in s.upper() and 'Fusion' in s: return 'BRAF'
    if 'PAX8' in s.upper() and 'PPARG' in s.upper(): return 'PAX8-PPARG'
    return None

cbio['fusion_class'] = cbio['eventInfo'].apply(classify_fusion)
cbio_classified = cbio.dropna(subset=['fusion_class']).copy()
cbio_classified['patient'] = [patient_id(s) for s in cbio_classified['sampleId']]

# Per-patient: take first fusion class (most common)
patient_fusion = cbio_classified.groupby('patient')['fusion_class'].first().reset_index()
print(f"Patients with classified fusion: {len(patient_fusion)}")
print(patient_fusion['fusion_class'].value_counts().to_dict())

# Merge frac × dm_like × fusion class
m = frac.merge(canonical[['patient','dm_like']], on='patient', how='inner')
m = m.merge(patient_fusion, on='patient', how='left')
m['fusion_class'] = m['fusion_class'].fillna('no_fusion')
print(f"\nMerged: {len(m)}, fusion classes: {m['fusion_class'].value_counts().to_dict()}")

# Within DM1 only — per-fusion-partner cell-type fraction
celltypes = ['B cell','Endothelial cell','Epithelial cell','Fibroblast','Malignant cell','Myeloid cell','NK cell','T cell']
dm1 = m[m['dm_like']=='DM1_like'].copy()
print(f"\nDM1: {len(dm1)}, by fusion class: {dm1['fusion_class'].value_counts().to_dict()}")

# Per-fusion mean fraction
fusion_means = dm1.groupby('fusion_class')[celltypes].mean()
fusion_n = dm1['fusion_class'].value_counts()
print(f"\nDM1 per-fusion-partner mean cell-type fractions:")
print(fusion_means.round(4))
fusion_means.to_csv(OUT/"v6A_per_fusion_partner_means.tsv", sep="\t")

# Pairwise Cohen's d: each fusion partner vs no_fusion (within DM1)
d_results = []
for fc in ['RET','NTRK','ALK','BRAF','PAX8-PPARG']:
    if fc not in dm1['fusion_class'].values: continue
    for ct in celltypes:
        a = dm1.loc[dm1['fusion_class']==fc, ct].values
        b = dm1.loc[dm1['fusion_class']=='no_fusion', ct].values
        cd = cohen_d(a, b)
        if len(a)>1 and len(b)>1:
            _, p = stats.mannwhitneyu(a, b, alternative='two-sided')
        else: p = 1
        d_results.append({'fusion_class':fc,'cell_type':ct,'cohen_d_vs_no':cd,'p':p,'n':len(a)})
d_df = pd.DataFrame(d_results)
d_df.to_csv(OUT/"v6A_per_fusion_vs_no_fusion_d.tsv", sep="\t", index=False)
print(f"\nPer-fusion-partner vs no-fusion within DM1 Cohen's d:")
pivot = d_df.pivot_table(index='fusion_class', columns='cell_type', values='cohen_d_vs_no')
print(pivot.round(3))
pivot.to_csv(OUT/"v6A_pivot.tsv", sep="\t")

# B. Sub-A/B × methylation interaction
print("\n=== B. Sub-A vs sub-B × methylation interaction ===")
sub = pd.read_csv("/data/thca/repo_results/d6p7_dm1_subcluster/dm1_subcluster_labels.tsv", sep="\t")
sub.columns = ['sample_id','sub_cluster']
meth = pd.read_csv("/data/thca/repo_results/audit_2026_04_30/round5/r5_2_sample_methylation_8gene.tsv", sep="\t")
sub['patient'] = sub['sample_id'].apply(patient_id)
meth['patient'] = meth['sample_short']

sub_meth = sub.merge(meth, on='patient', how='inner')
print(f"Sub × methylation: {len(sub_meth)}")
print(sub_meth['sub_cluster'].value_counts().to_dict())

meth_genes = ['DIO1','FOXE1','NKX2-1','PAX8','SLC5A5','TG','TPO','TSHR','mean_8g_beta']
B_results = []
for g in meth_genes:
    a = sub_meth.loc[sub_meth['sub_cluster']=='sub_A', g].values
    b = sub_meth.loc[sub_meth['sub_cluster']=='sub_B', g].values
    cd = cohen_d(a, b)
    if len(a)>1 and len(b)>1: _, p = stats.mannwhitneyu(a, b, alternative='two-sided')
    else: p = 1
    B_results.append({'gene':g,'sub_A_mean':a.mean(),'sub_B_mean':b.mean(),
                      'cohen_d_A_vs_B':cd,'p':p,'n_A':len(a),'n_B':len(b)})
B_df = pd.DataFrame(B_results)
B_df.to_csv(OUT/"v6B_subA_subB_methylation.tsv", sep="\t", index=False)
print(f"\nSub-A vs sub-B per-gene methylation Cohen's d:")
print(B_df.round(4).to_string(index=False))

print("\n=== DONE v6 ===")
