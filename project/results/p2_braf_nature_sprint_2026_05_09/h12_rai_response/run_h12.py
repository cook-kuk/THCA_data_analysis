"""H12 — RAI response x HT/DM1 axis in BRAF-cPTC.

Two cohorts:
  (A) GSE151179 (Colombo 2020) Clariom D — pre/post-RAI primaries + LN-mets, with patient
      response ('Avid' = remission/responder, 'Refractory' = RAI-refractory).
  (B) TCGA-THCA primary tumor cohort, RAI proxies = (i) `additional_radiation_therapy`
      = YES (post-surgery RT after RAI given), (ii) `new_tumor_event_after_initial_treatment`,
      and (iii) i_131 dose strata.

For each cohort we:
  1. Compute HT-13, FA-12, MAPK-9 panel scores (within-cohort z-mean) using genes available.
  2. Stratify by RAI-response / refractory status, by DM1/DM2, by driver and stage.
  3. Cohen's d, Mann-Whitney p, Fisher exact for DM x refractory rates, and
     Cox PFI for time-to-recurrence x DM x HT-quartile in BRAF-cPTC.

Output: h12_rai_results.tsv (long-form) + h12_summary.json + figure
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

OUT = Path('/home/seungho/personal/THCA_data_analysis/project/results/p2_braf_nature_sprint_2026_05_09/h12_rai_response')
OUT.mkdir(parents=True, exist_ok=True)

# ----------------------------------------------------------------------
# Panel definitions (matched to H4 conventions)
# ----------------------------------------------------------------------
HT13 = ['HLA-DRA','HLA-DRB1','HLA-DPA1','HLA-DPB1','HLA-DQA1','HLA-DQB1',
        'CD79A','CD79B','MS4A1','AICDA','CXCL13','CCR6','IFNG']
FA12 = ['CPT1A','CPT1B','ACADM','ACADVL','HADHA','HADHB','ECHS1','HMGCS2',
        'ACAA2','ACOX1','PPARA','PPARGC1A']
MAPK9 = ['DUSP4','DUSP5','DUSP6','SPRY1','SPRY2','SPRY4','ETV4','ETV5','PHLDA1']
PANELS = {'HT13': HT13, 'FA12': FA12, 'MAPK9': MAPK9}

# 8-gene DM panel (paper-1 standard)
DM8 = ['SLC5A5','TPO','TG','TSHR','PAX8','DIO1','DIO2','FOXE1']  # canonical RAI/lineage


def panel_zscore_mean(expr_long: pd.DataFrame, genes: list[str]) -> pd.Series:
    """expr_long: rows=samples, cols=genes (already z-scored within cohort)."""
    keep = [g for g in genes if g in expr_long.columns]
    if not keep:
        return pd.Series(np.nan, index=expr_long.index)
    return expr_long[keep].mean(axis=1)


def cohens_d(a: np.ndarray, b: np.ndarray) -> float:
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    if len(a) < 2 or len(b) < 2:
        return np.nan
    sa, sb = np.var(a, ddof=1), np.var(b, ddof=1)
    pooled = np.sqrt(((len(a)-1)*sa + (len(b)-1)*sb) / (len(a)+len(b)-2))
    if not np.isfinite(pooled) or pooled == 0:
        return np.nan
    return (np.mean(a) - np.mean(b)) / pooled


def ci95_d(d: float, n1: int, n2: int) -> tuple[float, float]:
    if not np.isfinite(d) or n1 < 2 or n2 < 2:
        return (np.nan, np.nan)
    se = np.sqrt((n1 + n2) / (n1 * n2) + d * d / (2 * (n1 + n2)))
    return (d - 1.96 * se, d + 1.96 * se)


def mannwhitney_p(a, b) -> float:
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    if len(a) < 2 or len(b) < 2:
        return np.nan
    try:
        return stats.mannwhitneyu(a, b, alternative='two-sided').pvalue
    except Exception:
        return np.nan


# ----------------------------------------------------------------------
# Cohort A: GSE151179
# ----------------------------------------------------------------------
def run_gse151179():
    """The aggressive_sprint_2026_05_06 file already gives us per-sample expression
    of 8 RAI/lineage genes + module annotations. We additionally pull the
    HT-13 module score from the precomputed paper3 file.
    """
    rows = []

    rai = pd.read_csv('/data/thca/repo_results/aggressive_sprint_2026_05_06/gse151179_rai_scores.tsv', sep='\t')
    rai = rai.rename(columns={'Unnamed: 0': 'sample_id'})

    mod = pd.read_csv('/data/thca/repo_results/paper3_ici_track_b_lite/scores_per_cohort/GSE151179_module_scores.tsv', sep='\t')

    df = rai.merge(mod, on='sample_id', how='left')

    # The HLA-class-II module here = subset of HT-13. We treat HLA_class_II as the
    # HT-immune proxy (pre-validated in H4 at d=+0.62 post-vs-pre).
    df['HT_proxy'] = df['HLA_class_II']  # primary HT proxy on Clariom D
    # Alternate composite: average of HLA-I, HLA-II, TLS, IFNG
    df['HT_composite'] = df[['HLA_class_I','HLA_class_II','TLS_CXCL13_like','IFNG_T_cell_inflamed']].mean(axis=1)
    df['rai_lineage'] = df['rai6_z_mean']

    # Tumor-only subset (drop non-neoplastic thyroid)
    tumors = df[df['is_tumor'] == True].copy()

    # ----- Contrast 1: post-RAI vs pre-RAI (all tumors) -----
    pre = tumors[tumors['is_pre_rai'] == True]
    post = tumors[tumors['is_pre_rai'] == False]
    for col, name in [('HT_proxy','HT_HLA-II'), ('HT_composite','HT_composite'),
                      ('rai_lineage','RAI_lineage'),
                      ('HLA_class_II','HLA_class_II'),
                      ('IFNG_T_cell_inflamed','IFNG_T_inflamed'),
                      ('TLS_CXCL13_like','TLS_CXCL13'),
                      ('myeloid_suppressive','myeloid_supp'),
                      ('thyroid_differentiation','thyroid_diff')]:
        d = cohens_d(post[col].dropna(), pre[col].dropna())
        ci = ci95_d(d, post[col].notna().sum(), pre[col].notna().sum())
        p = mannwhitney_p(post[col], pre[col])
        rows.append({'cohort':'GSE151179','contrast':'post_vs_pre_RAI_tumor',
                     'group_high':'post-RAI','group_low':'pre-RAI','score':name,
                     'n_high':int(post[col].notna().sum()),'n_low':int(pre[col].notna().sum()),
                     'cohens_d':d,'ci_lo':ci[0],'ci_hi':ci[1],'mw_p':p})

    # ----- Contrast 2: pre-RAI primaries Refractory vs Avid (response prediction) -----
    pre_prim = tumors[(tumors['is_pre_rai'] == True) & (tumors['is_primary'] == True)].copy()
    refr = pre_prim[pre_prim['patient_rai_responce'] == 'Refractory']
    avid = pre_prim[pre_prim['patient_rai_responce'] == 'Avid']
    for col, name in [('HT_proxy','HT_HLA-II'), ('HT_composite','HT_composite'),
                      ('rai_lineage','RAI_lineage'),
                      ('thyroid_differentiation','thyroid_diff'),
                      ('myeloid_suppressive','myeloid_supp')]:
        d = cohens_d(refr[col].dropna(), avid[col].dropna())
        ci = ci95_d(d, refr[col].notna().sum(), avid[col].notna().sum())
        p = mannwhitney_p(refr[col], avid[col])
        rows.append({'cohort':'GSE151179','contrast':'pre_RAI_primary_Refractory_vs_Avid',
                     'group_high':'Refractory','group_low':'Avid','score':name,
                     'n_high':int(refr[col].notna().sum()),'n_low':int(avid[col].notna().sum()),
                     'cohens_d':d,'ci_lo':ci[0],'ci_hi':ci[1],'mw_p':p})

    # ----- Contrast 3: BRAF-only pre-RAI primaries refractory vs avid -----
    braf_pre = pre_prim[pre_prim['lesion_class'].isin(['BRAFV600E','BRAFV600E+pTERT'])]
    refr_b = braf_pre[braf_pre['patient_rai_responce'] == 'Refractory']
    avid_b = braf_pre[braf_pre['patient_rai_responce'] == 'Avid']
    for col, name in [('HT_proxy','HT_HLA-II'), ('HT_composite','HT_composite'),
                      ('rai_lineage','RAI_lineage'),
                      ('thyroid_differentiation','thyroid_diff')]:
        d = cohens_d(refr_b[col].dropna(), avid_b[col].dropna())
        ci = ci95_d(d, refr_b[col].notna().sum(), avid_b[col].notna().sum())
        p = mannwhitney_p(refr_b[col], avid_b[col])
        rows.append({'cohort':'GSE151179','contrast':'BRAFV600E_pre_RAI_Refractory_vs_Avid',
                     'group_high':'Refractory','group_low':'Avid','score':name,
                     'n_high':int(refr_b[col].notna().sum()),'n_low':int(avid_b[col].notna().sum()),
                     'cohens_d':d,'ci_lo':ci[0],'ci_hi':ci[1],'mw_p':p})

    # ----- Contrast 4: RAI uptake at met site Yes vs No (post-RAI LN-mets) -----
    post_lnmet = tumors[tumors['is_pre_rai'] == False].copy()
    yes = post_lnmet[post_lnmet['rai_uptake_at_the_metastatic_site'] == 'Yes']
    no = post_lnmet[post_lnmet['rai_uptake_at_the_metastatic_site'] == 'No']
    for col, name in [('HT_proxy','HT_HLA-II'), ('HT_composite','HT_composite'),
                      ('rai_lineage','RAI_lineage'),
                      ('thyroid_differentiation','thyroid_diff')]:
        d = cohens_d(no[col].dropna(), yes[col].dropna())   # high=No (refractory) / low=Yes (avid)
        ci = ci95_d(d, no[col].notna().sum(), yes[col].notna().sum())
        p = mannwhitney_p(no[col], yes[col])
        rows.append({'cohort':'GSE151179','contrast':'post_RAI_LNmet_NoUptake_vs_YesUptake',
                     'group_high':'No_uptake','group_low':'Yes_uptake','score':name,
                     'n_high':int(no[col].notna().sum()),'n_low':int(yes[col].notna().sum()),
                     'cohens_d':d,'ci_lo':ci[0],'ci_hi':ci[1],'mw_p':p})

    # ----- Logistic: HT_proxy predicts Refractory (pre-RAI primaries) -----
    logit_summary = {}
    if len(pre_prim) >= 12 and pre_prim['patient_rai_responce'].nunique() == 2:
        from sklearn.linear_model import LogisticRegression
        from sklearn.preprocessing import StandardScaler
        sub = pre_prim.dropna(subset=['HT_proxy','rai_lineage','thyroid_differentiation']).copy()
        sub['y'] = (sub['patient_rai_responce'] == 'Refractory').astype(int)
        for predictor in ['HT_proxy','rai_lineage','thyroid_differentiation']:
            X = StandardScaler().fit_transform(sub[[predictor]].values)
            try:
                lr = LogisticRegression(penalty=None, max_iter=1000).fit(X, sub['y'].values)
                from scipy.stats import norm
                # bootstrap p
                rng = np.random.default_rng(0)
                coefs = []
                for _ in range(2000):
                    idx = rng.integers(0, len(sub), len(sub))
                    if len(np.unique(sub['y'].values[idx])) < 2:
                        continue
                    try:
                        coefs.append(LogisticRegression(penalty=None, max_iter=500).fit(X[idx], sub['y'].values[idx]).coef_[0,0])
                    except Exception:
                        pass
                coefs = np.array(coefs)
                if len(coefs) >= 100:
                    or_ = np.exp(lr.coef_[0,0])
                    or_lo, or_hi = np.exp(np.percentile(coefs, 2.5)), np.exp(np.percentile(coefs, 97.5))
                    p = 2 * min((coefs <= 0).mean(), (coefs >= 0).mean())
                    logit_summary[predictor] = {'OR':or_,'OR_lo':or_lo,'OR_hi':or_hi,'p_boot':p,
                                                'n':len(sub),'n_refractory':int(sub['y'].sum())}
            except Exception:
                pass

    return rows, logit_summary, df


# ----------------------------------------------------------------------
# Cohort B: TCGA-THCA RAI annotation
# ----------------------------------------------------------------------
def parse_dose(x):
    """Parse i_131_total_administered_dose (mixed mCi / millicures)."""
    if pd.isna(x):
        return np.nan
    s = str(x).lower().replace('millicures','').replace('mci','').strip()
    try:
        return float(s)
    except Exception:
        return np.nan


def run_tcga():
    rows = []
    m = pd.read_csv('/data/thca/repo_results/audit_2026_04_30/round9/r9_1_master_with_rai.tsv', sep='\t', low_memory=False)
    expr = pd.read_csv('/data/thca/data_processed/bulk_rnaseq/TCGA-THCA_rnaseq_expression_zscore.tsv', sep='\t', index_col=0)
    expr = expr.T  # rows = samples
    # match sample_id -> expr columns
    common = m['sample_id'].isin(expr.index)
    m = m[common].copy()
    expr = expr.loc[m['sample_id'].values]
    # compute panel scores within full TCGA primary (z-scores already applied)
    for name, genes in PANELS.items():
        m[f'panel_{name}'] = panel_zscore_mean(expr, genes).values
    m['panel_RAI8'] = panel_zscore_mean(expr, DM8).values

    # Parse RAI dose
    m['i131_dose_num'] = m['i_131_total_administered_dose'].apply(parse_dose)
    # RAI-refractory proxies
    m['rai_refractory_v1'] = (m['additional_radiation_therapy'] == 'YES').astype(int)
    m['rai_refractory_v2'] = (m['new_tumor_event_after_initial_treatment'] == 'YES').astype(int)
    m['rai_refractory_combined'] = ((m['rai_refractory_v1']==1) | (m['rai_refractory_v2']==1)).astype(int)

    # BRAF-cPTC stratum (per H6 convention)
    s1 = m[(m['molecular_subtype']=='BRAF_like') & (m['histology_subtype'].str.lower().str.contains('classical', na=False))].copy()
    if len(s1) < 20:  # fallback if histology empty -> use all BRAF_like
        s1 = m[m['molecular_subtype']=='BRAF_like'].copy()
    s1_label = 'BRAF_like (cPTC if labeled)' if 'classical' in str(m['histology_subtype'].dropna().iloc[0]).lower() else 'BRAF_like'

    # ---- Stratum-level contrasts ----
    strata = {
        'BRAF_like_all': m[m['molecular_subtype']=='BRAF_like'],
        'BRAF_like_cPTC': s1,
        'all_primary': m,
    }
    panels = ['panel_HT13','panel_FA12','panel_MAPK9','panel_RAI8']

    for sname, ss in strata.items():
        # 1) refractory vs not (combined proxy)
        for pcol in panels:
            for ref_col, label in [('rai_refractory_combined','refractory_combined'),
                                   ('rai_refractory_v1','postRT_post_RAI'),
                                   ('rai_refractory_v2','new_tumor_event')]:
                a = ss[ss[ref_col]==1][pcol]
                b = ss[ss[ref_col]==0][pcol]
                d = cohens_d(a.dropna(), b.dropna())
                ci = ci95_d(d, a.notna().sum(), b.notna().sum())
                p = mannwhitney_p(a, b)
                rows.append({'cohort':'TCGA-THCA','stratum':sname,'contrast':f'{label}_yes_vs_no',
                             'group_high':label,'group_low':'no','score':pcol.replace('panel_',''),
                             'n_high':int(a.notna().sum()),'n_low':int(b.notna().sum()),
                             'cohens_d':d,'ci_lo':ci[0],'ci_hi':ci[1],'mw_p':p})

        # 2) DM1 vs DM2 refractory rates (Fisher) within stratum
        for ref_col, label in [('rai_refractory_combined','refractory_combined'),
                               ('rai_refractory_v1','postRT'),
                               ('rai_refractory_v2','new_tumor_event')]:
            t = pd.crosstab(ss['dm'], ss[ref_col])
            for d1 in ['DM1','DM2','not_DM']:
                if d1 not in t.index or 0 not in t.columns or 1 not in t.columns:
                    continue
            if all(x in t.index for x in ['DM1','DM2']):
                table = [[t.loc['DM1',1] if 1 in t.columns else 0, t.loc['DM1',0] if 0 in t.columns else 0],
                         [t.loc['DM2',1] if 1 in t.columns else 0, t.loc['DM2',0] if 0 in t.columns else 0]]
                try:
                    odds, p_fish = stats.fisher_exact(table)
                except Exception:
                    odds, p_fish = np.nan, np.nan
                rate_dm1 = table[0][0] / max(1, sum(table[0]))
                rate_dm2 = table[1][0] / max(1, sum(table[1]))
                rows.append({'cohort':'TCGA-THCA','stratum':sname,'contrast':f'DM1_vs_DM2_{label}_rate',
                             'group_high':'DM1','group_low':'DM2','score':label,
                             'n_high':sum(table[0]),'n_low':sum(table[1]),
                             'cohens_d':rate_dm1 - rate_dm2,'ci_lo':np.nan,'ci_hi':np.nan,'mw_p':p_fish,
                             'extra':f'rate_DM1={rate_dm1:.3f}; rate_DM2={rate_dm2:.3f}; OR={odds:.3f}'})

        # 3) HT-quartile refractory rates (within BRAF-cPTC primarily)
        if sname in ['BRAF_like_cPTC','BRAF_like_all']:
            ss = ss.copy()
            ss['HT_q'] = pd.qcut(ss['panel_HT13'], 4, labels=['Q1','Q2','Q3','Q4'], duplicates='drop')
            t = pd.crosstab(ss['HT_q'], ss['rai_refractory_combined'])
            if 'Q1' in t.index and 'Q4' in t.index and 1 in t.columns:
                table = [[t.loc['Q4',1], t.loc['Q4',0] if 0 in t.columns else 0],
                         [t.loc['Q1',1], t.loc['Q1',0] if 0 in t.columns else 0]]
                try:
                    odds, p_fish = stats.fisher_exact(table)
                except Exception:
                    odds, p_fish = np.nan, np.nan
                rate_q4 = table[0][0] / max(1, sum(table[0]))
                rate_q1 = table[1][0] / max(1, sum(table[1]))
                rows.append({'cohort':'TCGA-THCA','stratum':sname,'contrast':'HT_Q4_vs_Q1_refractory_rate',
                             'group_high':'HT_Q4','group_low':'HT_Q1','score':'rai_refractory_combined',
                             'n_high':sum(table[0]),'n_low':sum(table[1]),
                             'cohens_d':rate_q4 - rate_q1,'ci_lo':np.nan,'ci_hi':np.nan,'mw_p':p_fish,
                             'extra':f'rate_Q4={rate_q4:.3f}; rate_Q1={rate_q1:.3f}; OR={odds:.3f}'})

    # 4) Cox: time-to-new-tumor-event (PFI proxy) by HT score & DM, BRAF-cPTC
    cox_rows = []
    try:
        from lifelines import CoxPHFitter
        sb = strata['BRAF_like_cPTC'].copy()
        # crude PFI = OS_days as time, new_tumor_event_after_initial_treatment as event
        sb['time'] = pd.to_numeric(sb['os_days'], errors='coerce')
        sb['event'] = (sb['new_tumor_event_after_initial_treatment']=='YES').astype(int)
        sb['age_n'] = pd.to_numeric(sb['age_at_diagnosis'], errors='coerce')
        sb['stage_n'] = sb['stage'].astype(str).str.extract(r'(IV|III|II|I)', expand=False).map({'I':1,'II':2,'III':3,'IV':4})
        sb['sex_n'] = (sb['sex'].astype(str).str.lower()=='male').astype(int)

        # HT score Cox in BRAF-cPTC
        for predictor, label in [('panel_HT13','HT13_per_SD'),
                                 ('panel_FA12','FA12_per_SD'),
                                 ('panel_MAPK9','MAPK9_per_SD'),
                                 ('panel_RAI8','RAI8_per_SD')]:
            sub = sb[['time','event',predictor,'age_n','stage_n','sex_n']].dropna()
            if sub['event'].sum() < 5 or len(sub) < 30:
                continue
            sub[predictor] = (sub[predictor] - sub[predictor].mean())/sub[predictor].std()
            try:
                cph = CoxPHFitter(penalizer=0.01).fit(sub, 'time', 'event')
                hr = float(np.exp(cph.params_[predictor]))
                lo = float(np.exp(cph.confidence_intervals_.loc[predictor].iloc[0]))
                hi = float(np.exp(cph.confidence_intervals_.loc[predictor].iloc[1]))
                p = float(cph.summary.loc[predictor,'p'])
                rows.append({'cohort':'TCGA-THCA','stratum':'BRAF_like_cPTC',
                             'contrast':'cox_PFI_per_SD',
                             'group_high':predictor,'group_low':'unit_SD','score':label,
                             'n_high':len(sub),'n_low':int(sub['event'].sum()),
                             'cohens_d':np.log(hr),'ci_lo':np.log(lo),'ci_hi':np.log(hi),'mw_p':p,
                             'extra':f'HR={hr:.3f} [{lo:.2f},{hi:.2f}]'})
            except Exception as e:
                rows.append({'cohort':'TCGA-THCA','stratum':'BRAF_like_cPTC',
                             'contrast':'cox_PFI_per_SD','group_high':predictor,'group_low':'fit_fail',
                             'score':label,'n_high':np.nan,'n_low':np.nan,
                             'cohens_d':np.nan,'ci_lo':np.nan,'ci_hi':np.nan,'mw_p':np.nan,
                             'extra':str(e)[:80]})

        # DM1/DM2 Cox in BRAF-cPTC
        sb['is_DM1'] = (sb['dm']=='DM1').astype(int)
        sb['is_DM2'] = (sb['dm']=='DM2').astype(int)
        sub = sb[['time','event','is_DM1','is_DM2','age_n','stage_n','sex_n']].dropna()
        if sub['event'].sum() >= 5 and len(sub) >= 30:
            try:
                cph = CoxPHFitter(penalizer=0.01).fit(sub, 'time', 'event')
                for v in ['is_DM1','is_DM2']:
                    hr = float(np.exp(cph.params_[v])); lo = float(np.exp(cph.confidence_intervals_.loc[v].iloc[0]))
                    hi = float(np.exp(cph.confidence_intervals_.loc[v].iloc[1])); p = float(cph.summary.loc[v,'p'])
                    rows.append({'cohort':'TCGA-THCA','stratum':'BRAF_like_cPTC',
                                 'contrast':'cox_PFI_DM_vs_others',
                                 'group_high':v,'group_low':'others','score':v,
                                 'n_high':int(sub[v].sum()),'n_low':int((sub[v]==0).sum()),
                                 'cohens_d':np.log(hr),'ci_lo':np.log(lo),'ci_hi':np.log(hi),'mw_p':p,
                                 'extra':f'HR={hr:.3f} [{lo:.2f},{hi:.2f}]'})
            except Exception:
                pass
    except ImportError:
        pass

    return rows, m


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------
def main():
    r1, logit_a, gse_df = run_gse151179()
    r2, tcga_m = run_tcga()
    rows = r1 + r2
    df = pd.DataFrame(rows)
    df.to_csv(OUT/'h12_rai_results.tsv', sep='\t', index=False)

    # also save tidy stratum tables for reference
    gse_df.to_csv(OUT/'h12_gse151179_per_sample.tsv', sep='\t', index=False)
    keep_cols = ['sample_id','dataset','molecular_subtype','histology_subtype','dm','driver_anchor',
                 'tert_status','panel_HT13','panel_FA12','panel_MAPK9','panel_RAI8',
                 'rai_refractory_v1','rai_refractory_v2','rai_refractory_combined',
                 'new_tumor_event_after_initial_treatment','additional_radiation_therapy',
                 'i131_dose_num','os_days','stage','age_at_diagnosis','sex']
    keep_cols = [c for c in keep_cols if c in tcga_m.columns]
    tcga_m[keep_cols].to_csv(OUT/'h12_tcga_per_sample.tsv', sep='\t', index=False)

    summary = {
        'gse151179_logistic_predictors_pre_RAI_primary_Refractory_vs_Avid': logit_a,
        'tcga_braf_cptc_n': int((tcga_m['molecular_subtype']=='BRAF_like').sum()),
        'tcga_braf_cptc_refractory_combined_n': int(((tcga_m['molecular_subtype']=='BRAF_like') & (tcga_m['rai_refractory_combined']==1)).sum()),
        'rows_in_h12_results': len(df),
    }
    with open(OUT/'h12_summary.json','w') as f:
        json.dump(summary, f, indent=2, default=str)

    print(df.to_string(index=False))
    print('\nSummary:', json.dumps(summary, indent=2, default=str))


if __name__ == '__main__':
    main()
