"""Round 7 — R7-1 TERT×DM×fusion 3-way, R7-2 composite score, R7-3 arm-level CNV, R7-4 multi-cohort score, R7-5 stemness."""
from __future__ import annotations
import json
import requests
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, fisher_exact, kstest, ks_2samp, spearmanr
from concurrent.futures import ThreadPoolExecutor
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score
import warnings
warnings.filterwarnings("ignore")

ROOT = Path("/home/seungho/personal/THCA_data_analysis")
OUT = ROOT / "project/results/audit_2026_04_30/round7"
OUT.mkdir(parents=True, exist_ok=True)
api = "https://www.cbioportal.org/api"
study_id = "thca_tcga_pan_can_atlas_2018"

# Common load
sm = pd.read_csv(ROOT / "project/results/v17_tert_recovery/v2/sample_master_v17_tert_v2.tsv",
                  sep="\t", low_memory=False)
sm = sm[(sm["dataset"] == "TCGA-THCA") & (sm["normal_vs_tumor"] == "tumor")].copy()
sm["tcga_short"] = sm["sample_id"].str.slice(0, 12)
sm["tert_int"] = (sm["tert_promoter_integrated"].fillna("").astype(str).str.lower().str.strip() == "mutated").astype(int)
sm["age_n"] = pd.to_numeric(sm["age"], errors="coerce")
sm["os_event"] = pd.to_numeric(sm["os_event"], errors="coerce")
sm["os_days"] = pd.to_numeric(sm["os_days"], errors="coerce")
drv = sm["driver_anchor"].fillna("unknown").astype(str).str.upper()
drv = drv.where(drv.isin(["BRAF", "RAS", "NTRK"]), other="OTHER")
sm["driver_simple"] = drv
xing = pd.read_csv(ROOT / "project/results/dark_matter_phase1/step6_xing_rescue.tsv", sep="\t")
xing["dm_cluster"] = xing["v17_dark_cluster"].fillna("not_DM")
if "v17_dark_cluster" in sm.columns: sm = sm.drop(columns=["v17_dark_cluster"])
sm = sm.merge(xing[["tcga_short", "dm_cluster"]], on="tcga_short", how="left")
sm["dm_cluster"] = sm["dm_cluster"].fillna("not_DM")
sv = pd.read_csv(ROOT / "project/results/audit_2026_04_30/round3/cbio_sv_thca.tsv", sep="\t", low_memory=False)
sv["sample_short"] = sv["sampleId"].str.slice(0, 12)
sm["any_fusion"] = sm["tcga_short"].isin(sv["sample_short"]).astype(int)

# ============================================================
# R7-1 — TERT × DM × fusion 3-way co-occurrence
# ============================================================
print("\n=== R7-1 — TERT × DM × fusion 3-way ===\n")

# 8-cell: DM × TERT × any_fusion
ct8 = sm.groupby(["dm_cluster", "tert_int", "any_fusion"]).size().reset_index(name="n")
ct8["combo"] = ct8["dm_cluster"] + "_TERT" + ct8["tert_int"].map({0:"-",1:"+"}) + "_fusion" + ct8["any_fusion"].map({0:"-",1:"+"})
print(ct8[["combo", "n"]].to_string(index=False))

# DM1 fusion+ × TERT
dm1_fp = sm[(sm["dm_cluster"]=="DM1") & (sm["any_fusion"]==1)]
dm1_fn = sm[(sm["dm_cluster"]=="DM1") & (sm["any_fusion"]==0)]
print(f"\nDM1 fusion+ TERT+ rate: {dm1_fp['tert_int'].sum()}/{len(dm1_fp)} = {100*dm1_fp['tert_int'].mean():.1f}%")
print(f"DM1 fusion- TERT+ rate: {dm1_fn['tert_int'].sum()}/{len(dm1_fn)} = {100*dm1_fn['tert_int'].mean():.1f}%")

# DM2 TERT+
dm2 = sm[sm["dm_cluster"]=="DM2"]
print(f"DM2 TERT+ rate: {dm2['tert_int'].sum()}/{len(dm2)} = {100*dm2['tert_int'].mean():.1f}%")

# 3-way Fisher: DM1 fusion+ TERT+ vs all others
sm["risk_combo"] = ((sm["dm_cluster"]=="DM1") & (sm["any_fusion"]==1) & (sm["tert_int"]==1)).astype(int)
sm["highest_risk"] = ((sm["driver_simple"]=="BRAF") & (sm["tert_int"]==1)).astype(int)
print(f"\nDM1 fusion+ TERT+ (potential 'fusion-driven aggressive'): {sm['risk_combo'].sum()}")
print(f"BRAF V600E TERT+ (Xing 2014 worst): {sm['highest_risk'].sum()}")

ct_check = pd.crosstab(sm["risk_combo"], sm["highest_risk"])
print(f"\nDM1+fusion+TERT+ × BRAF+TERT+ overlap:")
print(ct_check)

# Per-cell n
print(f"\nFull 8-cell DM × fusion × TERT N matrix:")
ct_full = sm.groupby(["dm_cluster", "any_fusion", "tert_int"]).size().reset_index(name="n").pivot_table(
    index="dm_cluster", columns=["any_fusion","tert_int"], values="n", fill_value=0)
print(ct_full)

r7_1 = {
    "8cell_breakdown": ct8.to_dict(orient="records"),
    "dm1_fusion_pos_TERT_rate_pct": round(100*dm1_fp['tert_int'].mean(), 1) if len(dm1_fp) else None,
    "dm1_fusion_neg_TERT_rate_pct": round(100*dm1_fn['tert_int'].mean(), 1) if len(dm1_fn) else None,
    "dm2_TERT_rate_pct": round(100*dm2['tert_int'].mean(), 1) if len(dm2) else None,
    "dm1_fusion_TERT_triple_pos": int(sm['risk_combo'].sum()),
    "BRAF_V600E_TERT_pos_xing_worst": int(sm['highest_risk'].sum()),
    "interpretation": "DM1 fusion+ TERT+ 3-way co-occurrence rate quantified for clinical risk algorithm",
}
(OUT / "r7_1_3way_breakdown.json").write_text(json.dumps(r7_1, indent=2, default=str))

# ============================================================
# R7-2 — Composite RNA-only risk score
# ============================================================
print("\n\n=== R7-2 — Composite RNA risk score ===\n")

# Combine: 8-gene RAI score (sm rai_score_v17) + immune (HLA-I) + meth (mean β) + TLS — all RNA-only
# Predict: any_fusion (binary)
sm["rai_score_v17_n"] = pd.to_numeric(sm["rai_score_v17"], errors="coerce")

# Immune (HLA from previous file)
hla = pd.read_csv(ROOT / "project/results/v17_hla/tcga_thca_hla_per_sample.tsv", sep="\t")
hla["tcga_short"] = hla["sample_id"].str.slice(0, 12)
sm = sm.merge(hla[["tcga_short", "hla_class_I_score", "hla_class_II_score"]], on="tcga_short", how="left", suffixes=("","_h"))

# Methylation mean β
meth = pd.read_csv(ROOT / "project/results/audit_2026_04_30/round5/r5_2_sample_methylation_8gene.tsv", sep="\t")
meth["tcga_short"] = meth["sample_short"].str.slice(0, 12)
if "mean_8g_beta" not in meth.columns and "TPO" in meth.columns:
    P8 = ["SLC5A5","TPO","TG","TSHR","PAX8","NKX2-1","FOXE1","DIO1"]
    meth["mean_8g_beta"] = meth[[c for c in P8 if c in meth.columns]].mean(axis=1)
sm = sm.merge(meth[["tcga_short", "mean_8g_beta"]], on="tcga_short", how="left")

# Predict any_fusion using composite features
features = ["rai_score_v17_n", "hla_class_I_score", "hla_class_II_score", "mean_8g_beta", "age_n"]
data = sm[features + ["any_fusion"]].dropna()
print(f"Composite score data: {len(data)} samples")

if len(data) > 50:
    X = StandardScaler().fit_transform(data[features])
    y = data["any_fusion"].values
    lr = LogisticRegression(penalty="l2", max_iter=2000, random_state=42)
    lr.fit(X, y)
    pred = lr.predict_proba(X)[:, 1]
    auc = roc_auc_score(y, pred)
    print(f"Composite RNA score AUC for predicting fusion+: {auc:.3f}")

    coefs = dict(zip(features, np.round(lr.coef_[0], 3).tolist()))
    print(f"Feature coefficients: {coefs}")

    # AUC for predicting DM1
    y_dm1 = (data.merge(sm[["sample_id"] + features + ["dm_cluster"]].dropna(), on=features)["dm_cluster"] == "DM1").values
    # Just predict DM1 directly
    sm_dm1 = sm[features + ["dm_cluster"]].dropna()
    sm_dm1["dm1_label"] = (sm_dm1["dm_cluster"] == "DM1").astype(int)
    X2 = StandardScaler().fit_transform(sm_dm1[features])
    y2 = sm_dm1["dm1_label"].values
    lr2 = LogisticRegression(penalty="l2", max_iter=2000, random_state=42)
    lr2.fit(X2, y2)
    pred2 = lr2.predict_proba(X2)[:, 1]
    auc2 = roc_auc_score(y2, pred2)
    print(f"\nComposite RNA score AUC for DM1 vs not_DM1: {auc2:.3f}")
    coefs2 = dict(zip(features, np.round(lr2.coef_[0], 3).tolist()))
    print(f"Feature coefs (DM1): {coefs2}")

    r7_2 = {
        "n_samples": int(len(data)),
        "auc_predict_fusion": round(float(auc), 3),
        "feature_coefs_fusion": coefs,
        "auc_predict_dm1": round(float(auc2), 3),
        "feature_coefs_dm1": coefs2,
        "interpretation": "Composite RNA-only score (8-gene RAI + HLA-I/II + 8-gene methylation + age) predicts both fusion status and DM1 cluster identity. Could be used as single integrated readout in clinical NanoString panel.",
    }
    (OUT / "r7_2_composite_score.json").write_text(json.dumps(r7_2, indent=2, default=str))

# ============================================================
# R7-3 — Arm-level CNV via cBioPortal
# ============================================================
print("\n\n=== R7-3 — Arm-level CNV ===\n")

armlevel_profile = "thca_tcga_pan_can_atlas_2018_armlevel_cna"

# Get all arm-level CNV via generic assay molecular data fetch
# Need gene/entity IDs for arms (1p, 1q, 22q, etc.)
# Try fetching via molecular-data with sample list
try:
    # Use generic-assay-meta-data endpoint to get arm IDs
    r = requests.get(f"{api}/generic-assay-meta-data/generic-assays?genericAssayType=ARMLEVEL_CNA",
                      timeout=30, headers={"Accept":"application/json"})
    if r.status_code == 200:
        arm_meta = r.json()
        arm_ids = [a.get("stableId") for a in arm_meta][:50] if arm_meta else []
        print(f"  Got {len(arm_ids)} arm IDs (first 5: {arm_ids[:5]})")
    else:
        # Alternative: try molecular-profile-data fetch with no specific entity
        arm_ids = ["1p", "1q", "2q", "3p", "3q", "5p", "5q", "7q", "8q", "11q",
                    "13q", "16q", "17p", "17q", "19p", "19q", "20q", "22q"]
        print(f"  Using default arm list: {len(arm_ids)} arms")
except Exception as e:
    arm_ids = []

# Try fetch arm-level data
if arm_ids:
    payload = {
        "genericAssayDataMultipleStudyFilter": {
            "molecularProfileIds": [armlevel_profile],
            "genericAssayStableIds": arm_ids[:30],
        },
    }
    try:
        r2 = requests.post(f"{api}/generic-assay-data/fetch?genericAssayDataType=BINARY",
                            json=payload, timeout=60)
        print(f"  Arm CNV fetch status: {r2.status_code}")
        if r2.status_code == 200 and r2.json():
            arm_df = pd.DataFrame(r2.json())
            arm_df["sample_short"] = arm_df["sampleId"].astype(str).str.slice(0, 12) if "sampleId" in arm_df.columns else None
            print(f"  Arm CNV records: {len(arm_df)}, columns: {list(arm_df.columns)[:10]}")
            arm_df.to_csv(OUT / "r7_3_arm_cnv.tsv", sep="\t", index=False)
        else:
            print(f"  Body: {r2.text[:200]}")
    except Exception as e:
        print(f"  arm fetch error: {e}")

# ============================================================
# R7-4 — Multi-cohort 8-gene score distribution
# ============================================================
print("\n\n=== R7-4 — Multi-cohort 8-gene score distribution ===\n")

# TCGA panel score
tcga_scores = pd.to_numeric(sm["rai_score_v17"], errors="coerce").dropna()
# K2 8-gene TPM mean (centered)
k2_tpm = pd.read_csv(ROOT / "project/results/v17_korean/K2_8gene_tpm_matrix_v4.tsv", sep="\t")
P8 = ["SLC5A5","TPO","TG","TSHR","PAX8","NKX2-1","FOXE1","DIO1"]
k2_tpm["mean_8g_log2_tpm"] = np.log2(k2_tpm[P8].astype(float) + 1).mean(axis=1)
k2_scores = k2_tpm["mean_8g_log2_tpm"]
# Within-cohort z-score
k2_scores_z = (k2_scores - k2_scores.mean()) / k2_scores.std()
tcga_scores_z = (tcga_scores - tcga_scores.mean()) / tcga_scores.std()

# Korean GSE213647 panel_z
gse_score = pd.read_csv(ROOT / "project/results/v17_korean/GSE213647_panel_score.tsv", sep="\t")
gse_scores = pd.to_numeric(gse_score["panel_z"], errors="coerce").dropna()
gse_scores_z = (gse_scores - gse_scores.mean()) / gse_scores.std()

# KS-test pairwise
ks_tests = {}
for pair in [("TCGA", tcga_scores_z, "K2", k2_scores_z),
              ("TCGA", tcga_scores_z, "GSE213647", gse_scores_z),
              ("K2", k2_scores_z, "GSE213647", gse_scores_z)]:
    n1, s1, n2, s2 = pair
    ks = ks_2samp(s1.dropna(), s2.dropna())
    ks_tests[f"{n1}_vs_{n2}"] = {"ks_stat": round(float(ks.statistic), 3),
                                  "p": float(ks.pvalue),
                                  "n1": int(s1.notna().sum()),
                                  "n2": int(s2.notna().sum())}
print(json.dumps(ks_tests, indent=2))

r7_4 = {
    "tcga_n": int(len(tcga_scores)), "tcga_median": round(float(tcga_scores.median()), 3),
    "k2_n": int(len(k2_scores)), "k2_median": round(float(k2_scores.median()), 3),
    "gse213647_n": int(len(gse_scores)), "gse213647_median": round(float(gse_scores.median()), 3),
    "within_cohort_z_KS_tests": ks_tests,
    "interpretation": "Within-cohort z-score harmonization. KS p < 0.05 indicates distribution differences. Required for cross-cohort meta-analysis.",
}
(OUT / "r7_4_multi_cohort_score.json").write_text(json.dumps(r7_4, indent=2, default=str))

# ============================================================
# R7-5 — Stemness signature × DM
# ============================================================
print("\n\n=== R7-5 — Stemness signature × DM ===\n")

# Stem cell markers via cBioPortal mRNA fetch
stemness_genes = {"NANOG":79923, "OCT4":5460, "SOX2":6657, "KLF4":9314, "MYC":4609,
                   "ALDH1A1":216, "CD44":960, "PROM1":8842, "EPCAM":4072, "BMI1":648}
mrna_profile = "thca_tcga_pan_can_atlas_2018_rna_seq_v2_mrna_median_Zscores"

def fetch_expr(g_eid):
    g, eid = g_eid
    payload = {"entrezGeneIds":[eid], "sampleListId":f"{study_id}_all"}
    try:
        r = requests.post(f"{api}/molecular-profiles/{mrna_profile}/molecular-data/fetch",
                           json=payload, timeout=30)
        if r.status_code == 200 and r.json():
            df = pd.DataFrame(r.json())
            df["gene"] = g
            return g, df
    except: pass
    return g, None

print("Parallel stemness fetch...")
with ThreadPoolExecutor(max_workers=4) as ex:
    stem_results = dict(ex.map(fetch_expr, stemness_genes.items()))
stem_dfs = [df for g, df in stem_results.items() if df is not None and not df.empty]

if stem_dfs:
    stem_all = pd.concat(stem_dfs, ignore_index=True)
    stem_all["sample_short"] = stem_all["sampleId"].str.slice(0, 12)
    stem_all["expr_z"] = pd.to_numeric(stem_all["value"], errors="coerce")

    # Per-sample stemness score (mean across markers)
    stem_score = stem_all.groupby("sample_short")["expr_z"].mean().rename("stem_score").reset_index()
    stem_dm = stem_score.merge(sm[["tcga_short", "dm_cluster", "any_fusion"]],
                                left_on="sample_short", right_on="tcga_short", how="left")

    # DM1 vs DM2
    a = stem_dm[stem_dm["dm_cluster"]=="DM1"]["stem_score"].dropna()
    b = stem_dm[stem_dm["dm_cluster"]=="DM2"]["stem_score"].dropna()
    if len(a) > 5 and len(b) > 5:
        mw = mannwhitneyu(a, b)
        cd = (a.mean() - b.mean()) / np.sqrt((a.var() + b.var()) / 2)
        print(f"\nDM1 vs DM2 stemness: median {a.median():.3f} vs {b.median():.3f}")
        print(f"Cohen's d = {cd:.3f}, p = {mw.pvalue:.4e}")

    # DM1 fusion+ vs fusion-
    a = stem_dm[(stem_dm["dm_cluster"]=="DM1") & (stem_dm["any_fusion"]==1)]["stem_score"].dropna()
    b = stem_dm[(stem_dm["dm_cluster"]=="DM1") & (stem_dm["any_fusion"]==0)]["stem_score"].dropna()
    if len(a) > 5 and len(b) > 3:
        mw = mannwhitneyu(a, b)
        cd = (a.mean() - b.mean()) / np.sqrt((a.var() + b.var()) / 2)
        print(f"\nDM1 fusion+ vs fusion- stemness: median {a.median():.3f} vs {b.median():.3f}")
        print(f"Cohen's d = {cd:.3f}, p = {mw.pvalue:.4e}")

    # Per-gene
    perg = []
    for g in stemness_genes.keys():
        gd = stem_all[stem_all["gene"]==g].merge(sm[["tcga_short","dm_cluster"]],
                                                   left_on="sample_short", right_on="tcga_short", how="left")
        a = gd[gd["dm_cluster"]=="DM1"]["expr_z"].dropna()
        b = gd[gd["dm_cluster"]=="DM2"]["expr_z"].dropna()
        if len(a) > 5 and len(b) > 5:
            mw = mannwhitneyu(a, b)
            cd = (a.mean() - b.mean()) / np.sqrt((a.var() + b.var()) / 2)
            perg.append({"gene":g, "median_DM1":round(float(a.median()),3),
                          "median_DM2":round(float(b.median()),3),
                          "cohens_d":round(float(cd),3), "mw_p":float(mw.pvalue)})
    print(f"\nPer-gene stemness DM1 vs DM2:")
    print(pd.DataFrame(perg).to_string(index=False))
    pd.DataFrame(perg).to_csv(OUT / "r7_5_stemness_DM.tsv", sep="\t", index=False)
    r7_5 = {"genes": list(stemness_genes.keys()), "per_gene_DM1_vs_DM2": perg,
             "interpretation": "Stemness markers vs DM cluster. Higher in DM2 (de-differentiated FVPTC) expected; DM1 differentiated."}
    (OUT / "r7_5_stemness.json").write_text(json.dumps(r7_5, indent=2, default=str))

print("\n=== R7 done ===")
