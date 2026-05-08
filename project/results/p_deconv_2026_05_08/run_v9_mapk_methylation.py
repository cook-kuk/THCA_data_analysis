"""v9: BRAF V600E vs RET-fusion methylation — is the BRAF→hyper-methylation MAPK-driven or BRAF-specific?"""
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

meth = pd.read_csv("/data/thca/repo_results/audit_2026_04_30/round5/r5_2_sample_methylation_8gene.tsv", sep="\t")
clin = pd.read_csv("/data/thca/repo_results/dark_matter_phase1/tcga_dm_master_with_pfi.tsv", sep="\t")
cbio = pd.read_csv("/data/thca/repo_results/audit_2026_04_30/round3/cbio_sv_thca.tsv", sep="\t")

meth['patient'] = meth['sample_short']
clin['patient'] = clin['tcga_short']
cbio['patient'] = cbio['sampleId'].apply(patient_id)

def fc(s):
    s = str(s) if pd.notna(s) else ''
    if 'RET' in s.upper(): return 'RET_fusion'
    if 'NTRK' in s.upper(): return 'NTRK_fusion'
    if 'ALK' in s.upper(): return 'ALK_fusion'
    return None
cbio['fusion_class'] = cbio['eventInfo'].apply(fc)
fusion = cbio.dropna(subset=['fusion_class']).groupby('patient')['fusion_class'].first().reset_index()

m = clin[['patient','has_braf_v600e','has_ras_mut']].merge(meth, on='patient', how='inner')
m = m.merge(fusion, on='patient', how='left')

# Build composite class:
def build_class(r):
    if r['has_braf_v600e']: return 'BRAF_V600E'
    if r['has_ras_mut']: return 'RAS_mut'
    if r['fusion_class']=='RET_fusion': return 'RET_fusion'
    if r['fusion_class']=='NTRK_fusion': return 'NTRK_fusion'
    if r['fusion_class']=='ALK_fusion': return 'ALK_fusion'
    return 'driver_neg'
m['driver'] = m.apply(build_class, axis=1)
counts = m['driver'].value_counts().to_dict()
print(f"Class composition (n={len(m)}): {counts}")

genes = ['DIO1','FOXE1','NKX2-1','PAX8','SLC5A5','TG','TPO','TSHR','mean_8g_beta']

# Pairwise contrasts:
contrasts = [
    ('BRAF_V600E', 'driver_neg'),
    ('RAS_mut',    'driver_neg'),
    ('RET_fusion', 'driver_neg'),
    ('NTRK_fusion','driver_neg'),
    ('BRAF_V600E', 'RAS_mut'),
    ('BRAF_V600E', 'RET_fusion'),
    ('RET_fusion', 'RAS_mut'),
]
rows = []
for ca, cb in contrasts:
    if (m['driver']==ca).sum()<3 or (m['driver']==cb).sum()<3: continue
    for g in genes:
        a = m.loc[m['driver']==ca, g].values
        b = m.loc[m['driver']==cb, g].values
        cd = cohen_d(a, b)
        _, p = stats.mannwhitneyu(a, b, alternative='two-sided')
        rows.append({'contrast':f"{ca}_vs_{cb}",'gene':g,'a_mean':a.mean(),'b_mean':b.mean(),
                     'cohen_d':cd,'p':p,'n_a':len(a),'n_b':len(b)})
df = pd.DataFrame(rows)
df.to_csv(OUT/"v9_mapk_methylation_pairwise.tsv", sep="\t", index=False)
pivot = df.pivot_table(index='contrast', columns='gene', values='cohen_d')
print(f"\nCohen's d (each contrast):")
print(pivot.round(3))
pivot.to_csv(OUT/"v9_mapk_methylation_pivot.tsv", sep="\t")

# Test: is BRAF V600E vs RET_fusion small (MAPK-driven) or large (BRAF-specific)?
print("\n--- Key contrast: BRAF V600E vs RET_fusion on mean_8g_β ---")
b_braf = m.loc[m['driver']=='BRAF_V600E','mean_8g_beta'].values
b_ret  = m.loc[m['driver']=='RET_fusion','mean_8g_beta'].values
b_ras  = m.loc[m['driver']=='RAS_mut','mean_8g_beta'].values
print(f"  BRAF V600E (n={len(b_braf)}): mean β = {b_braf.mean():.4f}")
print(f"  RET fusion  (n={len(b_ret)}): mean β = {b_ret.mean():.4f}")
print(f"  RAS mut     (n={len(b_ras)}): mean β = {b_ras.mean():.4f}")
print(f"  d(BRAF − RET) = {cohen_d(b_braf, b_ret):.3f}")
print(f"  d(BRAF − RAS) = {cohen_d(b_braf, b_ras):.3f}")
print(f"  d(RET − RAS)  = {cohen_d(b_ret, b_ras):.3f}")

print("\nInterpretation:")
print("  - If d(BRAF−RET) is small AND both d(BRAF−RAS), d(RET−RAS) are large:")
print("    → MAPK-pathway-driven (both MAPK-activating drivers hyper-methylate)")
print("  - If d(BRAF−RET) is large:")
print("    → BRAF-specific mechanism (not just MAPK output)")

print("\n=== DONE v9 ===")
