"""v7: within-DM1 OS by fusion class + TERT additive (NC paper-strengthening)."""
from pathlib import Path
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test, multivariate_logrank_test

OUT = Path("/home/seungho/personal/THCA_data_analysis/project/results/p_deconv_2026_05_08")

# Load TCGA-THCA clinical (OS)
clin = pd.read_csv("/data/thca/repo_results/dark_matter_phase1/tcga_dm_master_with_pfi.tsv", sep="\t")
print(f"Clinical: {len(clin)} samples, columns: {list(clin.columns)[:10]}")
print(clin.head(2))

# Load fractions + fusion class + canonical
frac = pd.read_csv(OUT/"fractions_nu_SVR.tsv", sep="\t", index_col=0)
canonical = pd.read_csv("/data/thca/repo_results/v17p3/tables/A2_dm_score_full_cohort.tsv", sep="\t")
canonical = canonical[canonical['dataset']=='TCGA-THCA']
cbio = pd.read_csv("/data/thca/repo_results/audit_2026_04_30/round3/cbio_sv_thca.tsv", sep="\t")
# tert_pos defined after patient col is assigned (below)
tert = None

def patient_id(s): return "-".join(s.split("-")[:3])

clin['patient'] = clin['tcga_short'].astype(str)
tert = clin[['patient','tert_pos']].copy()
tert['TERT_pos'] = tert['tert_pos'].astype(int)
canonical['patient'] = canonical['sample_id'].apply(patient_id)
cbio['patient'] = cbio['sampleId'].apply(patient_id)
print(f"TERT n={len(tert)}, pos={int(tert['TERT_pos'].sum())}")

# Fusion class
def classify(s):
    s = str(s) if pd.notna(s) else ''
    if 'RET' in s.upper(): return 'RET_fusion'
    if 'NTRK' in s.upper(): return 'NTRK_fusion'
    if 'ALK' in s.upper(): return 'ALK_fusion'
    if 'BRAF' in s.upper() and 'Fusion' in s: return 'BRAF_fusion'
    return None
cbio['fusion_class'] = cbio['eventInfo'].apply(classify)
fusion = cbio.dropna(subset=['fusion_class']).groupby('patient')['fusion_class'].first().reset_index()

# Merge
m = canonical[['patient','dm_like']].merge(clin[['patient','os_days','os_event']].rename(columns={'os_days':'OS_time','os_event':'OS_event'}), on='patient', how='inner')
m = m.merge(fusion, on='patient', how='left')
m['fusion_any'] = m['fusion_class'].notna().astype(int)
m['fusion_class'] = m['fusion_class'].fillna('no_fusion')
m = m.merge(tert[['patient','TERT_pos']].drop_duplicates(subset='patient'), on='patient', how='left')
m['TERT_pos'] = m['TERT_pos'].fillna(0).astype(int)

m = m.dropna(subset=['OS_time','OS_event']).reset_index(drop=True)
print(f"\nFinal merge: n={len(m)}, DM1={int((m['dm_like']=='DM1_like').sum())}, DM2={int((m['dm_like']=='DM2_like').sum())}")
print(f"  fusion_any={int(m['fusion_any'].sum())}, TERT_pos={int(m['TERT_pos'].sum())}")
print(f"  events={int(m['OS_event'].sum())}")

# Within DM1 — fusion class vs no_fusion OS
dm1 = m[m['dm_like']=='DM1_like'].copy()
print(f"\n=== Within DM1 (n={len(dm1)}, events={int(dm1['OS_event'].sum())}) ===")
print(f"Fusion classes: {dm1['fusion_class'].value_counts().to_dict()}")

# Multivariate logrank across fusion classes
res = multivariate_logrank_test(dm1['OS_time'], dm1['fusion_class'], dm1['OS_event'])
print(f"Multivariate logrank across fusion classes: chi2={res.test_statistic:.2f}, p={res.p_value:.4g}")

# DM1 fusion+ vs DM1 fusion-
fp = dm1[dm1['fusion_any']==1]
fn = dm1[dm1['fusion_any']==0]
lr = logrank_test(fp['OS_time'], fn['OS_time'], fp['OS_event'], fn['OS_event'])
print(f"\nDM1 fusion+ (n={len(fp)}, events={int(fp['OS_event'].sum())}) vs fusion- (n={len(fn)}, events={int(fn['OS_event'].sum())}): logrank p={lr.p_value:.4g}")

# Cox: within DM1 — fusion + TERT
try:
    cox = CoxPHFitter()
    df = dm1[['OS_time','OS_event','fusion_any','TERT_pos']].copy()
    cox.fit(df, 'OS_time', 'OS_event')
    print("\nDM1 Cox (fusion + TERT):")
    print(cox.summary[['exp(coef)','exp(coef) lower 95%','exp(coef) upper 95%','p']].round(4))
except Exception as e:
    print(f"Cox error: {e}")

# Whole cohort: DM1 + TERT joint
try:
    cox2 = CoxPHFitter()
    m_all = m.copy(); m_all['DM1'] = (m_all['dm_like']=='DM1_like').astype(int)
    df2 = m_all[['OS_time','OS_event','DM1','TERT_pos']]
    cox2.fit(df2, 'OS_time', 'OS_event')
    print("\nWhole cohort Cox (DM1 + TERT joint):")
    print(cox2.summary[['exp(coef)','exp(coef) lower 95%','exp(coef) upper 95%','p']].round(4))
except Exception as e:
    print(f"Cox2 error: {e}")

# Save
out_rows = []
out_rows.append({'analysis':'DM1_multi_fusion_logrank','stat':float(res.test_statistic),'p':float(res.p_value)})
out_rows.append({'analysis':'DM1_fusion_pos_vs_neg_logrank','stat':float(lr.test_statistic),'p':float(lr.p_value),'n_pos':len(fp),'n_neg':len(fn)})
pd.DataFrame(out_rows).to_csv(OUT/"v7_dm1_OS_logrank.tsv", sep="\t", index=False)
m.to_csv(OUT/"v7_merge_table.tsv", sep="\t", index=False)
print("\n=== DONE v7 ===")
