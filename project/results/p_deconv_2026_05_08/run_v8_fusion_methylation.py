"""v8: per-fusion-partner × 8-gene methylation within DM1."""
from pathlib import Path
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from scipy import stats

OUT = Path("/home/seungho/personal/THCA_data_analysis/project/results/p_deconv_2026_05_08")

def cohen_d(a,b):
    a,b = np.asarray(a),np.asarray(b)
    n1,n2=len(a),len(b)
    if n1<2 or n2<2: return 0
    pooled = np.sqrt(((n1-1)*a.var(ddof=1)+(n2-1)*b.var(ddof=1))/(n1+n2-2))
    return (a.mean()-b.mean())/pooled if pooled>0 else 0

def patient_id(s): return "-".join(s.split("-")[:3])

# Load methylation + canonical + cbio
meth = pd.read_csv("/data/thca/repo_results/audit_2026_04_30/round5/r5_2_sample_methylation_8gene.tsv", sep="\t")
canonical = pd.read_csv("/data/thca/repo_results/v17p3/tables/A2_dm_score_full_cohort.tsv", sep="\t")
canonical = canonical[canonical['dataset']=='TCGA-THCA']
cbio = pd.read_csv("/data/thca/repo_results/audit_2026_04_30/round3/cbio_sv_thca.tsv", sep="\t")

meth['patient'] = meth['sample_short']
canonical['patient'] = canonical['sample_id'].apply(patient_id)
cbio['patient'] = cbio['sampleId'].apply(patient_id)

def classify(s):
    s = str(s) if pd.notna(s) else ''
    if 'RET' in s.upper(): return 'RET_fusion'
    if 'NTRK' in s.upper(): return 'NTRK_fusion'
    if 'ALK' in s.upper(): return 'ALK_fusion'
    if 'BRAF' in s.upper() and 'Fusion' in s: return 'BRAF_fusion'
    return None
cbio['fusion_class'] = cbio['eventInfo'].apply(classify)
fusion = cbio.dropna(subset=['fusion_class']).groupby('patient')['fusion_class'].first().reset_index()

m = canonical[['patient','dm_like']].merge(meth, on='patient', how='inner')
m = m.merge(fusion, on='patient', how='left')
m['fusion_class'] = m['fusion_class'].fillna('no_fusion')

dm1 = m[m['dm_like']=='DM1_like'].copy()
print(f"DM1 × meth: n={len(dm1)}")
print(f"  fusion classes: {dm1['fusion_class'].value_counts().to_dict()}")

genes = ['DIO1','FOXE1','NKX2-1','PAX8','SLC5A5','TG','TPO','TSHR','mean_8g_beta']
fc_order = ['RET_fusion','NTRK_fusion','BRAF_fusion','ALK_fusion']

rows = []
for fc in fc_order:
    n_fc = (dm1['fusion_class']==fc).sum()
    if n_fc < 3: continue
    for g in genes:
        a = dm1.loc[dm1['fusion_class']==fc, g].values
        b = dm1.loc[dm1['fusion_class']=='no_fusion', g].values
        cd = cohen_d(a, b)
        if len(a)>1 and len(b)>1:
            _, p = stats.mannwhitneyu(a, b, alternative='two-sided')
        else: p=1
        rows.append({'fusion_class':fc,'gene':g,'fc_mean':a.mean(),'no_mean':b.mean(),
                    'cohen_d_vs_no':cd,'p':p,'n_fc':len(a),'n_no':len(b)})
df = pd.DataFrame(rows)
df.to_csv(OUT/"v8_fusion_methylation_d.tsv", sep="\t", index=False)

print("\nFusion-partner vs no-fusion methylation Cohen's d (within DM1):")
pivot = df.pivot_table(index='fusion_class', columns='gene', values='cohen_d_vs_no')
print(pivot.round(3))
pivot.to_csv(OUT/"v8_fusion_methylation_pivot.tsv", sep="\t")

# Also driver-class — BRAF V600E vs RAS vs other on methylation (whole cohort)
print("\n=== Driver class × methylation (full cohort) ===")
clin = pd.read_csv("/data/thca/repo_results/dark_matter_phase1/tcga_dm_master_with_pfi.tsv", sep="\t")
clin['patient'] = clin['tcga_short']
m2 = clin[['patient','has_braf_v600e','has_ras_mut','driver_anchor_v17']].merge(meth, on='patient', how='inner')
print(f"n={len(m2)}")

def driver(r):
    if r['has_braf_v600e']: return 'BRAF_V600E'
    if r['has_ras_mut']: return 'RAS_mut'
    return 'other'
m2['driver'] = m2.apply(driver, axis=1)
print(f"  driver: {m2['driver'].value_counts().to_dict()}")

rows2 = []
for g in genes:
    for fc in ['BRAF_V600E','RAS_mut']:
        a = m2.loc[m2['driver']==fc, g].values
        b = m2.loc[m2['driver']=='other', g].values
        cd = cohen_d(a, b)
        rows2.append({'driver':fc,'gene':g,'mean':a.mean(),'other_mean':b.mean(),
                      'cohen_d_vs_other':cd,'n':len(a)})
df2 = pd.DataFrame(rows2)
pivot2 = df2.pivot_table(index='driver', columns='gene', values='cohen_d_vs_other')
print(f"\nDriver vs other (driver-neg) methylation Cohen's d:")
print(pivot2.round(3))
df2.to_csv(OUT/"v8_driver_methylation_d.tsv", sep="\t", index=False)
pivot2.to_csv(OUT/"v8_driver_methylation_pivot.tsv", sep="\t")

print("\n=== DONE v8 ===")
