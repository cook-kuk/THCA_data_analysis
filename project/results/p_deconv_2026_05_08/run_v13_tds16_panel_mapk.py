"""v13: TDS-16 × 8-gene panel × MAPK — canonical-panel reproduction of MAPK→silencing.

Question: is the MAPK-driven silencing finding (v9-v12) specific to the deployable
8-gene compact readout, or is it a property of the full canonical Yoo 2014 TDS-16
thyroid differentiation score? If both panels behave identically, the 8-gene paper's
"compact readout, not cherry-pick" claim is bulletproof for Reviewer Q.

Outputs (all under p_deconv_2026_05_08/):
  v13_cross_cohort_corr.tsv           — MAPK × {Panel, TDS16, TDS8only} Spearman ρ for TCGA + Lee
  v13_per_driver_class_means.tsv      — TCGA per-driver-class score means (MAPK / Panel / TDS16)
  v13_per_driver_class_d.tsv          — pairwise Cohen's d (BRAF/RAS/RET/driver-neg) on the 3 scores
  v13_subAB_three_panels.tsv          — sub-A vs sub-B d(A-B) for MAPK / HT / Panel / TDS16 / TDS8only
  v13_per_gene_mapk_corr.tsv          — TDS-16 per-gene × MAPK_score Spearman ρ in TCGA + Lee
  v13_decile_trajectory.tsv           — MAPK-decile mean Panel/TDS16/TDS8only in TCGA + Lee
  v13_auc_panel_vs_tds16.tsv          — ROC-AUC of (MAPK-high vs MAPK-low) using Panel z vs TDS-16 z
"""
from pathlib import Path
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score, roc_curve

OUT = Path("/home/seungho/personal/THCA_data_analysis/project/results/p_deconv_2026_05_08")

PANEL_8 = ['DIO1','FOXE1','NKX2-1','PAX8','SLC5A5','TG','TPO','TSHR']
TDS_16  = ['DIO1','DIO2','DUOX1','DUOX2','FOXE1','GLIS3','NKX2-1','PAX8',
           'SLC26A4','SLC5A5','SLC5A8','TG','THRA','THRB','TPO','TSHR']
TDS_8ONLY = sorted(set(TDS_16) - set(PANEL_8))   # DIO2/DUOX1/DUOX2/GLIS3/SLC26A4/SLC5A8/THRA/THRB
MAPK = ['DUSP4','DUSP5','DUSP6','SPRY2','SPRY4','ETV4','ETV5','PHLDA1','CCND1']
HT   = ['HLA-DRA','HLA-DRB1','HLA-DPA1','HLA-DPB1','HLA-DQA1','HLA-DQB1',
        'CD79A','CD79B','MS4A1','AICDA','IGHM','IGKC','CXCL13','CCR6','IFNG']

def patient_id(s): return "-".join(s.split("-")[:3])

def cohen_d(a, b):
    a, b = np.asarray(a), np.asarray(b)
    if len(a) < 2 or len(b) < 2: return np.nan
    p = np.sqrt(((len(a)-1)*a.var(ddof=1)+(len(b)-1)*b.var(ddof=1))/(len(a)+len(b)-2))
    return (a.mean() - b.mean()) / p if p > 0 else np.nan

# =========================================================================
# 1. Load TCGA + Lee z-score, compute scores
# =========================================================================
print("=== Loading TCGA-THCA z-score ===")
tcga = pd.read_csv("/data/thca/data_processed/bulk_rnaseq/TCGA-THCA_rnaseq_expression_zscore.tsv",
                   sep="\t", index_col=0)
print(f"  shape: {tcga.shape}")

def score(z, genes):
    g = [x for x in genes if x in z.index]
    return z.loc[g].mean(axis=0), len(g)

t_mapk, n_t_mapk = score(tcga, MAPK)
t_pan,  n_t_pan  = score(tcga, PANEL_8)
t_tds,  n_t_tds  = score(tcga, TDS_16)
t_8o,   n_t_8o   = score(tcga, TDS_8ONLY)
t_ht,   n_t_ht   = score(tcga, HT)
print(f"  TCGA scored: MAPK={n_t_mapk}/9 PANEL={n_t_pan}/8 TDS16={n_t_tds}/16 TDS8only={n_t_8o}/8 HT={n_t_ht}/15")

print("\n=== Loading Lee/GSE213647 (Ensembl→symbol) ===")
lee = pd.read_csv("/data/thca/data_processed/bulk_rnaseq/GSE213647_rnaseq_expression_zscore.tsv",
                  sep="\t", index_col=0)
gmap = pd.read_csv("/data/thca/repo_results/v17p3/tables/F1_gene_recovery_mapping.tsv", sep="\t")
ens_to_sym = dict(zip(gmap['ensembl'], gmap['symbol']))
target = sorted(set(PANEL_8 + TDS_16 + MAPK + HT))
ens_target = [e for e, s in ens_to_sym.items() if s in target]
lee_sub = lee.loc[lee.index.intersection(ens_target)].copy()
lee_sub.index = [ens_to_sym[e] for e in lee_sub.index]
lee_sub = lee_sub.groupby(lee_sub.index).mean()
print(f"  Lee scored shape after symbol collapse: {lee_sub.shape}")
l_mapk, n_l_mapk = score(lee_sub, MAPK)
l_pan,  n_l_pan  = score(lee_sub, PANEL_8)
l_tds,  n_l_tds  = score(lee_sub, TDS_16)
l_8o,   n_l_8o   = score(lee_sub, TDS_8ONLY)
l_ht,   n_l_ht   = score(lee_sub, HT)
print(f"  Lee scored: MAPK={n_l_mapk}/9 PANEL={n_l_pan}/8 TDS16={n_l_tds}/16 TDS8only={n_l_8o}/8 HT={n_l_ht}/15")

# =========================================================================
# 2. Cross-cohort MAPK × {Panel, TDS16, TDS8only} Spearman ρ
# =========================================================================
print("\n=== A. Cross-cohort MAPK × thyroid-panel Spearman ρ ===")
rows = []
for cohort, mapk_s, pan_s, tds_s, t8o_s in [
    ("TCGA-THCA",        t_mapk, t_pan, t_tds, t_8o),
    ("Lee_GSE213647",    l_mapk, l_pan, l_tds, l_8o),
]:
    for label, score_s in [("Panel_8", pan_s), ("TDS_16", tds_s), ("TDS_8only", t8o_s)]:
        rho, p = stats.spearmanr(mapk_s, score_s)
        pr,  pp = stats.pearsonr(mapk_s, score_s)
        rows.append({'cohort':cohort,'panel':label,'n':len(mapk_s),
                     'spearman_rho':rho,'spearman_p':p,
                     'pearson_r':pr,'pearson_p':pp})
        print(f"  {cohort:14s}  MAPK × {label:10s}  ρ={rho:+.3f}  p={p:.2e}  (n={len(mapk_s)})")
pd.DataFrame(rows).to_csv(OUT/"v13_cross_cohort_corr.tsv", sep="\t", index=False)

# =========================================================================
# 3. Per-driver-class TCGA: MAPK / Panel / TDS16
# =========================================================================
print("\n=== B. Per-driver-class TCGA score means ===")
clin = pd.read_csv("/data/thca/repo_results/dark_matter_phase1/tcga_dm_master_with_pfi.tsv", sep="\t")
clin['patient'] = clin['tcga_short']
cbio = pd.read_csv("/data/thca/repo_results/audit_2026_04_30/round3/cbio_sv_thca.tsv", sep="\t")
cbio['patient'] = cbio['sampleId'].apply(patient_id)
def fc(s):
    s = str(s) if pd.notna(s) else ''
    if 'RET' in s.upper(): return 'RET_fusion'
    if 'NTRK' in s.upper(): return 'NTRK_fusion'
    if 'ALK' in s.upper(): return 'ALK_fusion'
    return None
cbio['fusion_class'] = cbio['eventInfo'].apply(fc)
fusion = cbio.dropna(subset=['fusion_class']).groupby('patient')['fusion_class'].first().reset_index()

samp = pd.DataFrame({
    'sample':t_mapk.index,
    'mapk':t_mapk.values, 'panel':t_pan.values,
    'tds16':t_tds.values, 'tds8only':t_8o.values,
    'ht':t_ht.values,
})
samp['patient'] = samp['sample'].apply(patient_id)
m = samp.merge(clin[['patient','has_braf_v600e','has_ras_mut']], on='patient', how='inner')
m = m.merge(fusion, on='patient', how='left')
def build_class(r):
    if r['has_braf_v600e']: return 'BRAF_V600E'
    if r['has_ras_mut']:    return 'RAS_mut'
    if r['fusion_class'] == 'RET_fusion':  return 'RET_fusion'
    if r['fusion_class'] == 'NTRK_fusion': return 'NTRK_fusion'
    return 'driver_neg'
m['driver'] = m.apply(build_class, axis=1)
counts = m['driver'].value_counts().to_dict()
print(f"  Class composition (n={len(m)}): {counts}")

means_rows = []
for cls in ['BRAF_V600E','RAS_mut','RET_fusion','NTRK_fusion','driver_neg']:
    sub = m[m['driver']==cls]
    if len(sub) < 3: continue
    means_rows.append({'driver':cls,'n':len(sub),
                       'mapk_mean':sub['mapk'].mean(),
                       'panel_mean':sub['panel'].mean(),
                       'tds16_mean':sub['tds16'].mean(),
                       'tds8only_mean':sub['tds8only'].mean()})
means_df = pd.DataFrame(means_rows)
print(means_df.round(3).to_string(index=False))
means_df.to_csv(OUT/"v13_per_driver_class_means.tsv", sep="\t", index=False)

# pairwise d on each score
print("\n  Pairwise Cohen's d (each score):")
contrasts = [('BRAF_V600E','driver_neg'),('RAS_mut','driver_neg'),('RET_fusion','driver_neg'),
             ('BRAF_V600E','RAS_mut'),('BRAF_V600E','RET_fusion')]
d_rows = []
for ca, cb in contrasts:
    if (m['driver']==ca).sum()<3 or (m['driver']==cb).sum()<3: continue
    for col in ['mapk','panel','tds16','tds8only']:
        a = m.loc[m['driver']==ca, col].values
        b = m.loc[m['driver']==cb, col].values
        d = cohen_d(a, b)
        _, p = stats.mannwhitneyu(a, b, alternative='two-sided')
        d_rows.append({'contrast':f"{ca}_vs_{cb}",'score':col,
                       'cohen_d':d,'p':p,'n_a':len(a),'n_b':len(b)})
d_df = pd.DataFrame(d_rows)
print(d_df.pivot_table(index='contrast', columns='score', values='cohen_d').round(3))
d_df.to_csv(OUT/"v13_per_driver_class_d.tsv", sep="\t", index=False)

# =========================================================================
# 4. sub-A vs sub-B convergence on TDS-16
# =========================================================================
print("\n=== C. sub-A vs sub-B d(A−B) on MAPK / HT / Panel / TDS16 / TDS8only ===")
sub_lab = pd.read_csv("/data/thca/repo_results/d6p7_dm1_subcluster/dm1_subcluster_labels.tsv", sep="\t")
sub_lab.columns = ['sample_id','sub_cluster']
sub_lab['patient'] = sub_lab['sample_id'].apply(patient_id)
ms = samp.merge(sub_lab[['patient','sub_cluster']], on='patient', how='inner')
print(f"  sub × TCGA RNA: {len(ms)} rows; {ms['sub_cluster'].value_counts().to_dict()}")
A = ms[ms['sub_cluster']=='sub_A']
B = ms[ms['sub_cluster']=='sub_B']

ab_rows = []
for col in ['mapk','ht','panel','tds16','tds8only']:
    d = cohen_d(A[col].values, B[col].values)
    _, p = stats.mannwhitneyu(A[col], B[col], alternative='two-sided')
    ab_rows.append({'metric':col,'d_A_minus_B':d,'p':p,'n_A':len(A),'n_B':len(B),
                    'A_mean':A[col].mean(),'B_mean':B[col].mean()})
ab_df = pd.DataFrame(ab_rows)
print(ab_df.round(3).to_string(index=False))
ab_df.to_csv(OUT/"v13_subAB_three_panels.tsv", sep="\t", index=False)

# =========================================================================
# 5. Per-gene MAPK ρ heatmap (TDS-16 × {TCGA, Lee})
# =========================================================================
print("\n=== D. Per-gene MAPK × TDS-16 Spearman ρ ===")
pg_rows = []
for g in TDS_16:
    if g in tcga.index:
        rho_t, p_t = stats.spearmanr(t_mapk, tcga.loc[g])
    else:
        rho_t, p_t = (np.nan, np.nan)
    if g in lee_sub.index:
        rho_l, p_l = stats.spearmanr(l_mapk, lee_sub.loc[g])
    else:
        rho_l, p_l = (np.nan, np.nan)
    pg_rows.append({'gene':g,'in_panel_8':g in PANEL_8,
                    'rho_TCGA':rho_t,'p_TCGA':p_t,
                    'rho_Lee':rho_l,'p_Lee':p_l})
pg_df = pd.DataFrame(pg_rows)
print(pg_df.round(3).to_string(index=False))
pg_df.to_csv(OUT/"v13_per_gene_mapk_corr.tsv", sep="\t", index=False)

# =========================================================================
# 6. Decile pseudotime: rank by MAPK, mean Panel/TDS16/TDS8only per decile
# =========================================================================
print("\n=== E. MAPK-decile pseudotime (TCGA + Lee) ===")
dec_rows = []
for cohort, mapk_s, pan_s, tds_s, t8o_s in [
    ("TCGA-THCA",     t_mapk, t_pan, t_tds, t_8o),
    ("Lee_GSE213647", l_mapk, l_pan, l_tds, l_8o),
]:
    df = pd.DataFrame({'mapk':mapk_s,'panel':pan_s,'tds16':tds_s,'tds8only':t8o_s})
    df['decile'] = pd.qcut(df['mapk'].rank(method='first'), 10, labels=False)
    g = df.groupby('decile').agg(mapk_mean=('mapk','mean'),
                                 panel_mean=('panel','mean'),
                                 tds16_mean=('tds16','mean'),
                                 tds8only_mean=('tds8only','mean'),
                                 n=('mapk','size')).reset_index()
    g['cohort'] = cohort
    dec_rows.append(g)
    rho_pan, _ = stats.spearmanr(g['decile'], g['panel_mean'])
    rho_tds, _ = stats.spearmanr(g['decile'], g['tds16_mean'])
    rho_8o,  _ = stats.spearmanr(g['decile'], g['tds8only_mean'])
    print(f"  {cohort:14s}  decile-rank ρ:  Panel={rho_pan:+.3f}  TDS16={rho_tds:+.3f}  TDS_8only={rho_8o:+.3f}")
dec_df = pd.concat(dec_rows, ignore_index=True)
dec_df.to_csv(OUT/"v13_decile_trajectory.tsv", sep="\t", index=False)

# =========================================================================
# 7. AUC: classify MAPK-high vs MAPK-low via Panel z vs TDS-16 z (median split)
# =========================================================================
print("\n=== F. ROC-AUC: MAPK-high vs MAPK-low using Panel vs TDS-16 score ===")
auc_rows = []
for cohort, mapk_s, pan_s, tds_s, t8o_s in [
    ("TCGA-THCA",     t_mapk, t_pan, t_tds, t_8o),
    ("Lee_GSE213647", l_mapk, l_pan, l_tds, l_8o),
]:
    y = (mapk_s.values > np.median(mapk_s.values)).astype(int)
    # high-MAPK should be low-thyroid → invert sign so AUC reflects "low score = high MAPK"
    auc_pan = roc_auc_score(y, -pan_s.values)
    auc_tds = roc_auc_score(y, -tds_s.values)
    auc_8o  = roc_auc_score(y, -t8o_s.values)
    auc_rows.append({'cohort':cohort,'AUC_Panel':auc_pan,'AUC_TDS16':auc_tds,
                     'AUC_TDS8only':auc_8o,'delta_TDS16_vs_Panel':auc_tds-auc_pan,'n':len(y)})
    print(f"  {cohort:14s}  AUC: Panel={auc_pan:.3f}  TDS16={auc_tds:.3f}  TDS_8only={auc_8o:.3f}"
          f"  Δ(TDS16−Panel)={auc_tds-auc_pan:+.3f}")
pd.DataFrame(auc_rows).to_csv(OUT/"v13_auc_panel_vs_tds16.tsv", sep="\t", index=False)

print("\n=== DONE v13 — outputs in p_deconv_2026_05_08/v13_*.tsv ===")
