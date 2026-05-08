"""v10: External RNA validation of MAPK→8-gene-suppression — TCGA + Lee + GSE76039."""
from pathlib import Path
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from scipy import stats

OUT = Path("/home/seungho/personal/THCA_data_analysis/project/results/p_deconv_2026_05_08")

MAPK_GENES = ['DUSP4','DUSP5','DUSP6','SPRY2','SPRY4','ETV4','ETV5','PHLDA1','CCND1']
PANEL_GENES = ['DIO1','FOXE1','NKX2-1','PAX8','SLC5A5','TG','TPO','TSHR']

def load_z(path, name):
    z = pd.read_csv(path, sep="\t", index_col=0)
    print(f"{name}: {z.shape}, sample sample: {z.columns[0]}")
    return z

def score_corr(z, name):
    available_mapk = [g for g in MAPK_GENES if g in z.index]
    available_panel = [g for g in PANEL_GENES if g in z.index]
    print(f"  {name}: MAPK genes={len(available_mapk)}/{len(MAPK_GENES)}, Panel genes={len(available_panel)}/{len(PANEL_GENES)}")
    if len(available_mapk)<5 or len(available_panel)<5: return None
    mapk_score = z.loc[available_mapk].mean(axis=0)
    panel_score = z.loc[available_panel].mean(axis=0)
    rho, p = stats.spearmanr(mapk_score, panel_score)
    pearson_r, pp = stats.pearsonr(mapk_score, panel_score)
    return {'cohort':name,'n_samples':z.shape[1],'mapk_n':len(available_mapk),
            'panel_n':len(available_panel),'spearman_rho':rho,'spearman_p':p,
            'pearson_r':pearson_r,'pearson_p':pp,
            'mapk_mean':mapk_score.mean(),'panel_mean':panel_score.mean()}

print("=== TCGA-THCA ===")
tcga = load_z("/data/thca/data_processed/bulk_rnaseq/TCGA-THCA_rnaseq_expression_zscore.tsv","TCGA-THCA")
print("\n=== Lee GSE213647 ===")
lee = load_z("/data/thca/data_processed/bulk_rnaseq/GSE213647_rnaseq_expression_zscore.tsv","Lee")
print("\n=== GSE126698 (Korean PTC+HT) ===")
gse126 = load_z("/data/thca/data_processed/bulk_rnaseq/GSE126698_rnaseq_expression_zscore.tsv","GSE126698")

results = []
for z, n in [(tcga,'TCGA-THCA'), (lee,'Lee_GSE213647'), (gse126,'GSE126698_Korean_PTC_HT')]:
    r = score_corr(z, n)
    if r: results.append(r)

df = pd.DataFrame(results)
df.to_csv(OUT/"v10_mapk_panel_correlation.tsv", sep="\t", index=False)
print(f"\n=== MAPK output × 8-gene panel correlation (cross-cohort) ===")
print(df[['cohort','n_samples','mapk_n','panel_n','spearman_rho','spearman_p']].round(4).to_string(index=False))

# Within TCGA: per-driver-class MAPK score
print("\n=== Within TCGA: per-driver-class MAPK score ===")
clin = pd.read_csv("/data/thca/repo_results/dark_matter_phase1/tcga_dm_master_with_pfi.tsv", sep="\t")
clin['patient'] = clin['tcga_short']
def patient_id(s): return "-".join(s.split("-")[:3])
cbio = pd.read_csv("/data/thca/repo_results/audit_2026_04_30/round3/cbio_sv_thca.tsv", sep="\t")
cbio['patient'] = cbio['sampleId'].apply(patient_id)
def fc(s):
    s = str(s) if pd.notna(s) else ''
    if 'RET' in s.upper(): return 'RET_fusion'
    if 'NTRK' in s.upper(): return 'NTRK_fusion'
    return None
cbio['fusion_class'] = cbio['eventInfo'].apply(fc)
fusion = cbio.dropna(subset=['fusion_class']).groupby('patient')['fusion_class'].first().reset_index()

available_mapk = [g for g in MAPK_GENES if g in tcga.index]
available_panel = [g for g in PANEL_GENES if g in tcga.index]
mapk_score = tcga.loc[available_mapk].mean(axis=0)
panel_score = tcga.loc[available_panel].mean(axis=0)
sample_df = pd.DataFrame({'sample':mapk_score.index, 'mapk':mapk_score.values, 'panel':panel_score.values})
sample_df['patient'] = sample_df['sample'].apply(patient_id)
sample_df = sample_df.merge(clin[['patient','has_braf_v600e','has_ras_mut']], on='patient', how='inner')
sample_df = sample_df.merge(fusion, on='patient', how='left')

def driver(r):
    if r['has_braf_v600e']: return 'BRAF_V600E'
    if r['has_ras_mut']: return 'RAS_mut'
    if r['fusion_class']=='RET_fusion': return 'RET_fusion'
    if r['fusion_class']=='NTRK_fusion': return 'NTRK_fusion'
    return 'driver_neg'
sample_df['driver'] = sample_df.apply(driver, axis=1)
print(f"  n={len(sample_df)}, {sample_df['driver'].value_counts().to_dict()}")

print("\nMAPK score by driver class (mean ± std):")
g = sample_df.groupby('driver').agg({'mapk':['mean','std','count'], 'panel':['mean','std']})
print(g.round(3))
g.to_csv(OUT/"v10_driver_mapk_panel_means.tsv", sep="\t")

# Cohen's d MAPK
def cd(a,b):
    a,b = np.asarray(a),np.asarray(b)
    if len(a)<2 or len(b)<2: return 0
    p = np.sqrt(((len(a)-1)*a.var(ddof=1)+(len(b)-1)*b.var(ddof=1))/(len(a)+len(b)-2))
    return (a.mean()-b.mean())/p if p>0 else 0

print("\nMAPK score Cohen's d:")
contrasts = [('BRAF_V600E','RAS_mut'),('RET_fusion','RAS_mut'),('BRAF_V600E','RET_fusion'),
             ('BRAF_V600E','driver_neg'),('RAS_mut','driver_neg'),('RET_fusion','driver_neg')]
rows = []
for ca, cb in contrasts:
    a = sample_df.loc[sample_df['driver']==ca,'mapk'].values
    b = sample_df.loc[sample_df['driver']==cb,'mapk'].values
    if len(a)<3 or len(b)<3: continue
    d_mapk = cd(a, b)
    a2 = sample_df.loc[sample_df['driver']==ca,'panel'].values
    b2 = sample_df.loc[sample_df['driver']==cb,'panel'].values
    d_panel = cd(a2, b2)
    print(f"  {ca} vs {cb}: MAPK d={d_mapk:+.3f}, Panel d={d_panel:+.3f}")
    rows.append({'contrast':f"{ca}_vs_{cb}",'mapk_d':d_mapk,'panel_d':d_panel,'n_a':len(a),'n_b':len(b)})
pd.DataFrame(rows).to_csv(OUT/"v10_driver_mapk_cohens_d.tsv", sep="\t", index=False)

print("\n=== DONE v10 ===")
