"""Audit 2026-04-30 — P1 (BRAF-/TERT- subset survival) + P7 (sample_master 6 angles) + P4 partial (DM1 sub-clustering).

Outputs to project/results/audit_2026_04_30/.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import multivariate_logrank_test, logrank_test
from scipy.stats import chi2_contingency, fisher_exact, mannwhitneyu

ROOT = Path("/home/seungho/personal/THCA_data_analysis")
OUT = ROOT / "project/results/audit_2026_04_30"
OUT.mkdir(parents=True, exist_ok=True)

# Load data
sm = pd.read_csv(ROOT / "project/results/v17_tert_recovery/v2/sample_master_v17_tert_v2.tsv",
                  sep="\t", low_memory=False)
sm = sm[(sm["dataset"] == "TCGA-THCA") & (sm["normal_vs_tumor"] == "tumor")].copy()
sm["tcga_short"] = sm["sample_id"].str.slice(0, 12)
sm["tert_int"] = (sm["tert_promoter_integrated"].fillna("").astype(str).str.lower().str.strip() == "mutated").astype(int)
sm["os_event"] = pd.to_numeric(sm["os_event"], errors="coerce")
sm["os_days"] = pd.to_numeric(sm["os_days"], errors="coerce")
sm["age_n"] = pd.to_numeric(sm["age"], errors="coerce")
sm["sex_male"] = (sm["sex"].astype(str).str.lower().str.strip() == "male").astype(int)
sm["stage_advanced"] = sm["clinical_stage"].astype(str).str.contains("III|IV", case=False, na=False).astype(int)
drv = sm["driver_anchor"].fillna("unknown").astype(str).str.upper()
drv = drv.where(drv.isin(["BRAF", "RAS", "NTRK"]), other="OTHER")
sm["driver_simple"] = drv

xing = pd.read_csv(ROOT / "project/results/dark_matter_phase1/step6_xing_rescue.tsv", sep="\t")
xing["dm_cluster"] = xing["v17_dark_cluster"].fillna("not_DM")
# sample_master also has v17_dark_cluster — drop to avoid suffix
if "v17_dark_cluster" in sm.columns:
    sm = sm.drop(columns=["v17_dark_cluster"])
sm = sm.merge(xing[["tcga_short", "xing_group", "dm_cluster"]], on="tcga_short", how="left")
sm["dm_cluster"] = sm["dm_cluster"].fillna("not_DM")
sm = sm.rename(columns={"dm_cluster": "v17_dark_cluster"})

# ==============================================================
# P1 — BRAF-/TERT- subset survival stratification
# ==============================================================
print("\n=== P1 — BRAF-/TERT- subset survival ===\n")
braf_tert_neg = sm[sm["xing_group"] == "BRAF-/TERT-"].copy()
braf_tert_neg["dm_3group"] = braf_tert_neg["v17_dark_cluster"].fillna("not_DM")
print("3-group breakdown:")
print(braf_tert_neg["dm_3group"].value_counts())

surv = braf_tert_neg.dropna(subset=["os_event", "os_days"])
print(f"\nWith survival: {len(surv)} / {len(braf_tert_neg)}")
print("Events per group:")
print(surv.groupby("dm_3group")["os_event"].agg(["sum", "count"]))

# Cox univariate per group (DM2 vs DM1)
p1_cox_rows = []
for ref, comp in [("DM1", "DM2"), ("DM1", "not_DM"), ("DM2", "not_DM")]:
    pair = surv[surv["dm_3group"].isin([ref, comp])].copy()
    pair["G"] = (pair["dm_3group"] == comp).astype(int)
    if pair["G"].nunique() < 2 or pair["os_event"].sum() < 2:
        p1_cox_rows.append({"ref": ref, "comp": comp, "n_ref": (pair["G"]==0).sum(), "n_comp": (pair["G"]==1).sum(),
                             "events_ref": pair[pair["G"]==0]["os_event"].sum(),
                             "events_comp": pair[pair["G"]==1]["os_event"].sum(),
                             "hr": None, "p": None, "note": "insufficient events"})
        continue
    try:
        c = CoxPHFitter(penalizer=0.05)
        c.fit(pair[["os_days", "os_event", "G"]].rename(columns={"os_days":"T", "os_event":"E"}),
              duration_col="T", event_col="E", show_progress=False)
        coef = c.params_.loc["G"]
        ci_lo = c.confidence_intervals_.loc["G", "95% lower-bound"]
        ci_hi = c.confidence_intervals_.loc["G", "95% upper-bound"]
        p = c.summary.loc["G", "p"]
        p1_cox_rows.append({"ref": ref, "comp": comp,
                             "n_ref": int((pair["G"]==0).sum()), "n_comp": int((pair["G"]==1).sum()),
                             "events_ref": int(pair[pair["G"]==0]["os_event"].sum()),
                             "events_comp": int(pair[pair["G"]==1]["os_event"].sum()),
                             "hr": round(np.exp(coef), 3),
                             "ci_lo": round(np.exp(ci_lo), 3), "ci_hi": round(np.exp(ci_hi), 3),
                             "p": round(p, 4)})
    except Exception as e:
        p1_cox_rows.append({"ref": ref, "comp": comp, "error": str(e)})

p1_cox_df = pd.DataFrame(p1_cox_rows)
p1_cox_df.to_csv(OUT / "p1_cox_hr_table.tsv", sep="\t", index=False)
print("\nCox HR (pairwise):")
print(p1_cox_df.to_string(index=False))

# 3-group omnibus logrank
try:
    lr = multivariate_logrank_test(surv["os_days"], surv["dm_3group"], surv["os_event"])
    omnibus_p = float(lr.p_value)
except Exception:
    omnibus_p = None
print(f"\n3-group omnibus logrank p = {omnibus_p}")

# Multivariate Cox: outcome ~ DM_cluster + age + stage + sex
mv_data = surv.dropna(subset=["age_n", "stage_advanced", "sex_male"])
mv_data["dm_DM2"] = (mv_data["dm_3group"] == "DM2").astype(int)
mv_data["dm_notDM"] = (mv_data["dm_3group"] == "not_DM").astype(int)
mv_cols = ["os_days", "os_event", "dm_DM2", "dm_notDM", "age_n", "sex_male", "stage_advanced"]
mv_cox_rows = []
try:
    c = CoxPHFitter(penalizer=0.1)  # heavier ridge for very small N
    c.fit(mv_data[mv_cols], duration_col="os_days", event_col="os_event", show_progress=False)
    for var in c.summary.index:
        coef = c.summary.loc[var, "coef"]
        ci_lo = c.summary.loc[var, "coef lower 95%"]
        ci_hi = c.summary.loc[var, "coef upper 95%"]
        p = c.summary.loc[var, "p"]
        mv_cox_rows.append({"var": var, "hr": round(np.exp(coef), 3),
                             "ci_lo": round(np.exp(ci_lo), 3), "ci_hi": round(np.exp(ci_hi), 3),
                             "p": round(p, 4)})
except Exception as e:
    mv_cox_rows.append({"error": str(e)})
mv_cox_df = pd.DataFrame(mv_cox_rows)
mv_cox_df.to_csv(OUT / "p1_multivariate_cox.tsv", sep="\t", index=False)
print(f"\nMultivariate Cox (n={len(mv_data)}, events={int(mv_data['os_event'].sum())}):")
print(mv_cox_df.to_string(index=False))

# Continuous 8-gene RAI score Cox (in BRAF-/TERT- subset)
sm_btn = surv.dropna(subset=["rai_score_v17"]).copy()
sm_btn["rai_score_v17"] = pd.to_numeric(sm_btn["rai_score_v17"], errors="coerce")
sm_btn = sm_btn.dropna(subset=["rai_score_v17"])
cont_cox_rows = []
try:
    c = CoxPHFitter(penalizer=0.05)
    c.fit(sm_btn[["os_days", "os_event", "rai_score_v17"]],
          duration_col="os_days", event_col="os_event", show_progress=False)
    coef = c.summary.loc["rai_score_v17", "coef"]
    ci_lo = c.summary.loc["rai_score_v17", "coef lower 95%"]
    ci_hi = c.summary.loc["rai_score_v17", "coef upper 95%"]
    p = c.summary.loc["rai_score_v17", "p"]
    cont_cox_rows.append({"score": "rai_score_v17 (per unit)", "hr_per_unit": round(np.exp(coef), 3),
                           "ci_lo": round(np.exp(ci_lo), 3), "ci_hi": round(np.exp(ci_hi), 3),
                           "p": round(p, 4),
                           "note": "lower score = more dedifferentiated; HR < 1 means score↑ → safer"})
except Exception as e:
    cont_cox_rows.append({"error": str(e)})
pd.DataFrame(cont_cox_rows).to_csv(OUT / "p1_continuous_cox.tsv", sep="\t", index=False)
print(f"\nContinuous score Cox (BRAF-/TERT- subset, n={len(sm_btn)}, events={int(sm_btn['os_event'].sum())}):")
print(pd.DataFrame(cont_cox_rows).to_string(index=False))

p1_summary = {
    "cohort": "TCGA-THCA BRAF-/TERT- (Xing dark matter)",
    "n_total": int(len(braf_tert_neg)),
    "n_with_survival": int(len(surv)),
    "events_total": int(surv["os_event"].sum()),
    "underpowered_warning": "YES — events < 10 per group" if surv["os_event"].sum() < 10 else "moderate",
    "dm_3group_n": braf_tert_neg["dm_3group"].value_counts().to_dict(),
    "omnibus_logrank_p": omnibus_p,
}
(OUT / "p1_summary.json").write_text(json.dumps(p1_summary, indent=2), encoding="utf-8")

# ==============================================================
# P7 — sample_master 6 추가 angles (TCGA + use whatever is available)
# ==============================================================
print("\n\n=== P7 — sample_master 6 추가 angles ===\n")
# Use 8-cell or DM cluster as stratifier — we have DM cluster from xing
sm["dm_cluster"] = sm["v17_dark_cluster"].fillna("not_DM")

p7_results = []

# Angle 1-4: Need multifocality / ETE / LN_mets / M1 — TCGA에는 없음. K2에서 활용
# K2 처리
k2_meta = pd.read_csv(ROOT / "project/results/dark_matter_phase2/k2_yoo2016_mutations_parsed.tsv",
                      sep="\t", low_memory=False)
print("K2 columns:", [c for c in k2_meta.columns if "Multif" in c or "extension" in c.lower() or "invasion" in c.lower() or "metastas" in c.lower()])

k2_clean = k2_meta[k2_meta["mol_subtype_label"].notna()].copy()
phen_map = {
    "Multifocality": "Multifocality",
    "ETE": "Extrathyroidal extension",
    "LymphInvasion": "Lymphatic invasion",
    "VascularInvasion": "Blood vassel invasion",
    "DistantMets": "Distant metastasis",
}
for short, col in phen_map.items():
    if col not in k2_clean.columns: continue
    valid = k2_clean[col].astype(str).str.match(r"^[01]$").fillna(False)
    sub = k2_clean.loc[valid].copy()
    sub["_v"] = sub[col].astype(int)
    if len(sub) < 30: continue
    ct = pd.crosstab(sub["mol_subtype_label"], sub["_v"])
    if ct.shape[1] < 2: continue
    chi2, p, _, _ = chi2_contingency(ct.values + 0.5)  # Yates-like correction
    n = ct.values.sum()
    cramers_v = float(np.sqrt(chi2 / (n * min(ct.shape) - 1))) if n > 0 else None
    p7_results.append({"cohort": "K2 (n=180)", "angle": short, "test": "chi-square",
                        "n": int(len(sub)), "chi2": round(chi2, 3),
                        "p": round(p, 4), "cramers_v": round(cramers_v, 3) if cramers_v else None,
                        "note": f"K2 mol_subtype_label crosstab: {dict(zip(*np.unique(sub['mol_subtype_label'], return_counts=True)))}"})

# Angle 5: Age × DM cluster (TCGA — sample_master)
age_data = sm.dropna(subset=["age_n", "dm_cluster"])
dm1_age = age_data[age_data["dm_cluster"] == "DM1"]["age_n"]
dm2_age = age_data[age_data["dm_cluster"] == "DM2"]["age_n"]
notdm_age = age_data[age_data["dm_cluster"] == "not_DM"]["age_n"]

if len(dm1_age) >= 5 and len(dm2_age) >= 5:
    mw_age = mannwhitneyu(dm1_age, dm2_age, alternative="two-sided")
    cohen_d_age = (dm1_age.mean() - dm2_age.mean()) / np.sqrt((dm1_age.var() + dm2_age.var()) / 2)
    p7_results.append({"cohort": "TCGA-THCA", "angle": "Age × DM cluster (DM1 vs DM2)", "test": "Mann-Whitney + Cohen's d",
                        "n": int(len(dm1_age) + len(dm2_age)),
                        "p": round(float(mw_age.pvalue), 4),
                        "cohens_d": round(float(cohen_d_age), 3),
                        "median_dm1": round(float(dm1_age.median()), 1),
                        "median_dm2": round(float(dm2_age.median()), 1)})
# Young-onset (<45) enrichment
sm["young_onset"] = (sm["age_n"] < 45).astype(int)
ct_young = pd.crosstab(sm["dm_cluster"], sm["young_onset"])
ct_young_dm = ct_young.loc[["DM1", "DM2"]] if "DM1" in ct_young.index and "DM2" in ct_young.index else None
if ct_young_dm is not None and ct_young_dm.shape[1] == 2 and (ct_young_dm.values >= 5).all():
    chi2, p, _, _ = chi2_contingency(ct_young_dm.values + 0.5)
    p7_results.append({"cohort": "TCGA-THCA", "angle": "Young-onset (<45) × DM cluster",
                        "test": "chi-square",
                        "n": int(ct_young_dm.values.sum()),
                        "p": round(p, 4),
                        "young_pct_dm1": round(100*ct_young_dm.loc["DM1", 1] / ct_young_dm.loc["DM1"].sum(), 1),
                        "young_pct_dm2": round(100*ct_young_dm.loc["DM2", 1] / ct_young_dm.loc["DM2"].sum(), 1)})

# Angle 6: Sex × DM cluster (TCGA)
ct_sex = pd.crosstab(sm["dm_cluster"], sm["sex_male"])
if "DM1" in ct_sex.index and "DM2" in ct_sex.index:
    ct_sex_dm = ct_sex.loc[["DM1", "DM2"]]
    if ct_sex_dm.shape[1] == 2:
        chi2, p, _, _ = chi2_contingency(ct_sex_dm.values + 0.5)
        p7_results.append({"cohort": "TCGA-THCA", "angle": "Sex × DM cluster (M:F ratio)",
                            "test": "chi-square",
                            "n": int(ct_sex_dm.values.sum()),
                            "p": round(p, 4),
                            "male_pct_dm1": round(100*ct_sex_dm.loc["DM1", 1] / ct_sex_dm.loc["DM1"].sum(), 1),
                            "male_pct_dm2": round(100*ct_sex_dm.loc["DM2", 1] / ct_sex_dm.loc["DM2"].sum(), 1)})

# Aggressive flag × DM cluster
agg_data = sm.dropna(subset=["aggressive_flag", "dm_cluster"])
agg_data["agg_yes"] = (agg_data["aggressive_flag"] == "yes").astype(int)
ct_agg = pd.crosstab(agg_data["dm_cluster"], agg_data["agg_yes"])
if "DM1" in ct_agg.index and "DM2" in ct_agg.index:
    ct_agg_dm = ct_agg.loc[["DM1", "DM2"]]
    if ct_agg_dm.shape[1] == 2:
        chi2, p, _, _ = chi2_contingency(ct_agg_dm.values + 0.5)
        p7_results.append({"cohort": "TCGA-THCA", "angle": "Aggressive_flag × DM cluster",
                            "test": "chi-square",
                            "n": int(ct_agg_dm.values.sum()),
                            "p": round(p, 4),
                            "agg_pct_dm1": round(100*ct_agg_dm.loc["DM1", 1] / ct_agg_dm.loc["DM1"].sum(), 1),
                            "agg_pct_dm2": round(100*ct_agg_dm.loc["DM2", 1] / ct_agg_dm.loc["DM2"].sum(), 1)})

p7_df = pd.DataFrame(p7_results)
p7_df.to_csv(OUT / "p7_six_angles.tsv", sep="\t", index=False)
print(p7_df.to_string(index=False))

# ==============================================================
# P4 partial — DM1 sub-clustering (외부 fusion/methylation/CNV 없음)
# ==============================================================
print("\n\n=== P4 partial — DM1 sub-clustering (no fusion/meth/CNV data) ===\n")
dm1_subset = sm[sm["v17_dark_cluster"] == "DM1"].copy()
print(f"DM1 n={len(dm1_subset)}")

# Sub-clustering attempt on available scores
score_cols = ["tds_score", "rai_score_v17", "tds16_score_v17", "dedifferentiation_proxy_score"]
score_cols = [c for c in score_cols if c in dm1_subset.columns]
dm1_scores = dm1_subset[score_cols].apply(pd.to_numeric, errors="coerce").dropna()
print(f"DM1 with full scores: {len(dm1_scores)}")

# Standardize + KMeans K=2..5
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import silhouette_score
sc = StandardScaler().fit_transform(dm1_scores)
sub_results = []
for k in range(2, 6):
    if len(sc) < k: continue
    km = KMeans(n_clusters=k, random_state=42, n_init=10)
    labs = km.fit_predict(sc)
    sil = silhouette_score(sc, labs)
    sub_results.append({"k": k, "silhouette": round(float(sil), 3),
                         "cluster_sizes": np.bincount(labs).tolist(),
                         "interpretable": sil > 0.3})
print("DM1 sub-clustering attempts:")
print(pd.DataFrame(sub_results).to_string(index=False))
pd.DataFrame(sub_results).to_csv(OUT / "p4_dm1_subclustering.tsv", sep="\t", index=False)

# DM1 patient driver_anchor breakdown (since DM1 contains BRAF-neg primarily)
dm1_drivers = dm1_subset["driver_simple"].value_counts().to_dict()
dm1_summary = {
    "n_dm1": int(len(dm1_subset)),
    "n_with_full_scores": int(len(dm1_scores)),
    "dm1_driver_breakdown": dm1_drivers,
    "subclustering_attempts": sub_results,
    "best_k": int(max(sub_results, key=lambda x: x["silhouette"])["k"]) if sub_results else None,
    "best_silhouette": round(max(s["silhouette"] for s in sub_results), 3) if sub_results else None,
    "interpretable_subcluster": any(s.get("interpretable") for s in sub_results),
    "limitation": "Without TCGA Fusion DB, 450K methylation, GISTIC2 CNV, only score-space sub-clustering possible. Phase 2 P2-C full pipeline needs external download.",
    "next_steps": [
        "Download TCGA Fusion Database (Hu Y 2018 NAR — http://www.tumorfusions.org)",
        "Download TCGA-THCA 450K methylation via GDC API",
        "Download TCGA-THCA GISTIC2 CNV results",
        "Re-run with all 4 modalities for definitive DM1 mechanism characterization",
    ],
}
(OUT / "p4_dm1_partial_summary.json").write_text(json.dumps(dm1_summary, indent=2), encoding="utf-8")
print(f"\nDM1 summary: {dm1_summary['best_silhouette']} max silhouette (interpretable: {dm1_summary['interpretable_subcluster']})")

print("\n\n=== Done. Outputs in", OUT, "===")
