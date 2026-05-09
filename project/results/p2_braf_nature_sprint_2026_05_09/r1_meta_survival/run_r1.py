#!/usr/bin/env python3
"""
R1 — REINFORCEMENT meta-survival pooling for Paper 1+2 BRAF Nature sprint.

Goal: tighten the DM2-vs-not_DM PFI HR (H6 lead row, n=25/5 events, HR=5.91
[1.79, 19.5]) via random-effects + fixed-effects meta over all available
cohorts that carry both DM-status and survival.

Cohort-availability audit (full prose in R1_REPORT.md):
  TCGA-THCA          DM-status from HM450 8-gene methylation. PFI/OS/DSS.
                      The ONLY cohort where DM2-vs-not_DM is separable
                      (DM2 = methylation-classifier-positive secondary
                      class; needs both classifier + non-DM cases).
  Landa GSE76039     RNA-proxy DM1 only (microarray; no HM450). All
                      advanced (PDTC + ATC). DM1 vs not-DM1 ~ ATC vs PDTC
                      so DM2-vs-not_DM is NOT testable. OK for DM1xTERT.
  Lee 2024 GSE213647 No survival, no TERT. EXCLUDED.
  K2 / PRJEB11591    No TERT, no survival. EXCLUDED.
  GSE286332          n=18, no survival. EXCLUDED.
  Pozdeyev 2018      Mutations only, no RNA, no DM-status. EXCLUDED.
  cBio thyroid_mskcc_2016 = same Landa cohort.

For DM2-vs-not_DM (the #1 weak spot) the only available pooled grain
within TCGA is by **driver_anchor sub-cohort** (BRAF / RAS / NBNR), which
gives 3 disjoint subsets with events. We pool those (within-TCGA meta).
We additionally report TCGA full primary as "discovery" reference and
TCGA BRAF-cPTC stratum as "lead row".

For DM1xTERT we have 2 truly independent cohorts (TCGA + Landa).

Outputs:
  r1_per_cohort_hr.tsv  -- cohort x contrast x HR/CI/p/n/events/log_HR/SE
  r1_pooled_meta.tsv    -- pooled HR + I^2 + Q-test (FE + DL + REML)
  r1_leave_one_out.tsv  -- LOO sensitivity
  r1_forest_data.tsv    -- ready for forest plot
  R1_REPORT.md          -- under 400 words, lead w/ pooled DM2 PFI HR
"""

import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from lifelines.exceptions import ConvergenceError, ConvergenceWarning
from scipy import stats

warnings.filterwarnings("ignore", category=ConvergenceWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=UserWarning)

ROOT = Path("/home/seungho/personal/THCA_data_analysis/project/results/p2_braf_nature_sprint_2026_05_09")
OUT = ROOT / "r1_meta_survival"
OUT.mkdir(parents=True, exist_ok=True)

H6_MERGED = ROOT / "h6_survival" / "h6_merged_clinical.tsv"
LANDA_MERGED = ROOT / "h8_tert_replication" / "h8_landa_merged.tsv"

# =============================================================================
# Helpers
# =============================================================================

def cox_hr(df, time_col, event_col, primary, covariates=None, penalizer=0.01):
    covariates = covariates or []
    cols = [primary] + [c for c in covariates if c != primary]
    sub = df[cols + [time_col, event_col]].dropna().copy()
    n = len(sub)
    events = int(sub[event_col].sum())
    out = dict(n=n, events=events, HR=np.nan, ci_lo=np.nan, ci_hi=np.nan,
               p=np.nan, log_HR=np.nan, SE_logHR=np.nan, note="")
    if events < 3 or n < 8:
        out["note"] = f"too_few (n={n}, events={events})"
        return out
    keep = [c for c in cols if sub[c].nunique(dropna=True) > 1]
    if primary not in keep:
        out["note"] = "primary_no_variance"
        return out
    try:
        cph = CoxPHFitter(penalizer=penalizer)
        cph.fit(sub[keep + [time_col, event_col]],
                duration_col=time_col, event_col=event_col)
        s = cph.summary.loc[primary]
        out["HR"] = float(s["exp(coef)"])
        out["ci_lo"] = float(s["exp(coef) lower 95%"])
        out["ci_hi"] = float(s["exp(coef) upper 95%"])
        out["p"] = float(s["p"])
        out["log_HR"] = float(s["coef"])
        out["SE_logHR"] = float(s["se(coef)"])
    except (ConvergenceError, np.linalg.LinAlgError, ValueError) as e:
        out["note"] = f"cox_fail:{type(e).__name__}"
    return out


def meta_inverse_variance(log_hr, se_logHR, method="REML"):
    log_hr = np.asarray(log_hr, dtype=float)
    se = np.asarray(se_logHR, dtype=float)
    keep = np.isfinite(log_hr) & np.isfinite(se) & (se > 0)
    log_hr = log_hr[keep]
    se = se[keep]
    k = len(log_hr)
    if k < 2:
        return dict(method=method, k=k, pooled_logHR=np.nan, pooled_SE=np.nan,
                    pooled_HR=np.nan, ci_lo=np.nan, ci_hi=np.nan, p=np.nan,
                    Q=np.nan, Q_df=np.nan, Q_p=np.nan, I2=np.nan, tau2=np.nan,
                    note="meta needs k>=2")
    w_fe = 1.0 / (se ** 2)
    pooled_fe = np.sum(w_fe * log_hr) / np.sum(w_fe)
    Q = np.sum(w_fe * (log_hr - pooled_fe) ** 2)
    df = k - 1
    Q_p = 1 - stats.chi2.cdf(Q, df) if df > 0 else np.nan
    c = np.sum(w_fe) - np.sum(w_fe ** 2) / np.sum(w_fe)
    tau2_DL = max(0.0, (Q - df) / c) if c > 0 else 0.0
    if method == "REML":
        tau2 = tau2_DL
        for _ in range(200):
            w = 1.0 / (se ** 2 + tau2)
            mu = np.sum(w * log_hr) / np.sum(w)
            num = np.sum((w ** 2) * ((log_hr - mu) ** 2 - (se ** 2 + tau2))) \
                  + np.sum(w) - np.sum((w ** 2) / np.sum(w))
            denom = np.sum(w ** 2) - np.sum(w ** 2) / np.sum(w)
            if denom <= 0:
                break
            new_tau2 = max(0.0, tau2 + num / denom)
            if abs(new_tau2 - tau2) < 1e-8:
                tau2 = new_tau2
                break
            tau2 = new_tau2
    elif method == "DL":
        tau2 = tau2_DL
    elif method == "FE":
        tau2 = 0.0
    else:
        tau2 = tau2_DL
    w = 1.0 / (se ** 2 + tau2)
    pooled = np.sum(w * log_hr) / np.sum(w)
    se_pooled = np.sqrt(1.0 / np.sum(w))
    z = pooled / se_pooled
    p = 2 * (1 - stats.norm.cdf(abs(z)))
    ci_lo = pooled - 1.96 * se_pooled
    ci_hi = pooled + 1.96 * se_pooled
    I2 = max(0.0, (Q - df) / Q) * 100 if Q > 0 else 0.0
    return dict(
        method=method, k=k,
        pooled_logHR=pooled, pooled_SE=se_pooled,
        pooled_HR=float(np.exp(pooled)),
        ci_lo=float(np.exp(ci_lo)), ci_hi=float(np.exp(ci_hi)),
        p=float(p), Q=float(Q), Q_df=int(df), Q_p=float(Q_p),
        I2=float(I2), tau2=float(tau2),
    )


# =============================================================================
# DATA LOAD + ENCODE
# =============================================================================
tcga = pd.read_csv(H6_MERGED, sep="\t")
landa = pd.read_csv(LANDA_MERGED, sep="\t")
print(f"[TCGA] n={len(tcga)}; [Landa] n={len(landa)}")

def tert_flag(x):
    if pd.isna(x):
        return np.nan
    return 1 if str(x).lower() in ("mutated", "1", "true") else 0

tcga["TERT"] = tcga["tert_promoter_integrated"].apply(tert_flag)
tcga["DM1"] = (tcga["dm"] == "DM1").astype(int)
tcga["DM2"] = (tcga["dm"] == "DM2").astype(int)
tcga["comutant"] = ((tcga["DM1"] == 1) & (tcga["TERT"] == 1)).astype(int)
tcga["BRAF_like"] = (tcga["molecular_subtype"] == "BRAF_like").astype(int)
tcga["RAS_like"] = (tcga["molecular_subtype"] == "RAS_like").astype(int)

landa["comutant"] = ((landa["DM1_bin"] == 1) & (landa["TERT"] == 1)).astype(int)

cov_age_stage_sex = ["age_yr", "adv_stage", "sex_m"]

# =============================================================================
# Per-cohort Cox runs
# =============================================================================
per_cohort = []

# --- DM2 vs not_DM in disjoint TCGA driver sub-cohorts (BRAF / RAS / NBNR) ---
# These are mutually exclusive samples so they form a within-TCGA meta.
driver_strata = [
    ("TCGA_BRAF_like",  tcga[tcga["molecular_subtype"] == "BRAF_like"]),
    ("TCGA_RAS_like",   tcga[tcga["molecular_subtype"] == "RAS_like"]),
    ("TCGA_NBNR",       tcga[~tcga["molecular_subtype"].isin(["BRAF_like", "RAS_like"])]),
]
# H6 lead-row stratum (BRAF + cPTC) — kept as the published anchor row, not in meta
anchor_strata = [
    ("TCGA_BRAF_cPTC_anchor",
     tcga[(tcga["molecular_subtype"] == "BRAF_like")
          & (tcga["histology_subtype"] == "cPTC")]),
    ("TCGA_full_primary_S3",  tcga),
]

ENDPOINTS = [("PFI", "PFI", "PFI.time"),
             ("OS",  "OS",  "OS.time"),
             ("DSS", "DSS", "DSS.time")]

for name, df in driver_strata + anchor_strata:
    sub = df[df["dm"].isin(["DM2", "not_DM"])].copy()
    sub["DM2"] = (sub["dm"] == "DM2").astype(int)
    for ep_label, ev_col, time_col in ENDPOINTS:
        r = cox_hr(sub, time_col, ev_col, "DM2",
                    covariates=["DM2"] + cov_age_stage_sex)
        per_cohort.append({"cohort": name, "contrast": "DM2_vs_notDM",
                            "endpoint": ep_label,
                            "model": "Cox + age/stage/sex", **r})

# --- DM1 vs DM2 (protective) per driver stratum + anchor + Landa ---
for name, df in driver_strata + anchor_strata:
    sub = df[df["dm"].isin(["DM1", "DM2"])].copy()
    sub["DM1"] = (sub["dm"] == "DM1").astype(int)
    for ep_label, ev_col, time_col in ENDPOINTS:
        r = cox_hr(sub, time_col, ev_col, "DM1",
                    covariates=["DM1"] + cov_age_stage_sex)
        per_cohort.append({"cohort": name, "contrast": "DM1_vs_DM2",
                            "endpoint": ep_label,
                            "model": "Cox + age/stage/sex", **r})

# Landa DM1 vs not-DM1 (proxy)
sub_l = landa.dropna(subset=["OS_event", "OS_months"]).copy()
sub_l["DM1"] = sub_l["DM1_bin"]
r = cox_hr(sub_l, "OS_months", "OS_event", "DM1",
            covariates=["DM1", "AGE_num"])
per_cohort.append({"cohort": "Landa_GSE76039_PDTC_ATC",
                    "contrast": "DM1_vs_notDM1_proxy",
                    "endpoint": "OS",
                    "model": "Cox + age (RNA-proxy DM)", **r})

# --- DM1xTERT comutant: TCGA + Landa OS ---
# TCGA full primary (driver-adjusted)
sub_t = tcga.dropna(subset=["DM1", "TERT"]).copy()
for ep_label, ev_col, time_col in ENDPOINTS:
    sub_te = sub_t.dropna(subset=[ev_col, time_col]).copy()
    r = cox_hr(sub_te, time_col, ev_col, "comutant",
                covariates=["comutant", "BRAF_like", "RAS_like",
                            "age_yr", "adv_stage", "sex_m"])
    per_cohort.append({"cohort": "TCGA_full_primary_S3",
                        "contrast": "DM1xTERT_comutant_vs_others_ADJ",
                        "endpoint": ep_label,
                        "model": "Cox + driver+age+stage+sex", **r})
    r = cox_hr(sub_te, time_col, ev_col, "comutant",
                covariates=["comutant"])
    per_cohort.append({"cohort": "TCGA_full_primary_S3",
                        "contrast": "DM1xTERT_comutant_vs_others_UNADJ",
                        "endpoint": ep_label,
                        "model": "Cox unadjusted", **r})

# Landa OS (adjusted is_ATC + age, and unadjusted)
r = cox_hr(sub_l, "OS_months", "OS_event", "comutant",
            covariates=["comutant", "is_ATC", "AGE_num"])
per_cohort.append({"cohort": "Landa_GSE76039_PDTC_ATC",
                    "contrast": "DM1xTERT_comutant_vs_others_ADJ",
                    "endpoint": "OS",
                    "model": "Cox + is_ATC+age", **r})
r = cox_hr(sub_l, "OS_months", "OS_event", "comutant",
            covariates=["comutant"])
per_cohort.append({"cohort": "Landa_GSE76039_PDTC_ATC",
                    "contrast": "DM1xTERT_comutant_vs_others_UNADJ",
                    "endpoint": "OS", "model": "Cox unadjusted", **r})

# Save per-cohort
pc = pd.DataFrame(per_cohort)
pc.to_csv(OUT / "r1_per_cohort_hr.tsv", sep="\t", index=False)
print(f"\n[wrote] r1_per_cohort_hr.tsv ({len(pc)} rows)")
print("\n=== HEADLINE per-cohort ===")
print(pc[["cohort", "contrast", "endpoint", "n", "events", "HR", "ci_lo",
          "ci_hi", "p", "note"]].to_string(index=False, max_colwidth=42))

# =============================================================================
# META RUNS
# =============================================================================
def meta_block(contrast, endpoint, scope_label, cohorts_subset=None,
                cohorts_pred=None, source_df=pc):
    if cohorts_subset is not None:
        sub = source_df[
            (source_df["contrast"] == contrast)
            & (source_df["endpoint"] == endpoint)
            & (source_df["cohort"].isin(cohorts_subset))
            & (source_df["log_HR"].notna())
        ]
    else:
        sub = source_df[
            (source_df["contrast"] == contrast)
            & (source_df["endpoint"] == endpoint)
            & (cohorts_pred(source_df["cohort"]))
            & (source_df["log_HR"].notna())
        ]
    rows = []
    for method in ["FE", "DL", "REML"]:
        m = meta_inverse_variance(sub["log_HR"].values,
                                    sub["SE_logHR"].values, method=method)
        m_with = {
            "contrast": contrast, "endpoint": endpoint,
            "scope": scope_label, "cohorts_in_meta": ",".join(sub["cohort"].tolist()),
            **m,
        }
        rows.append(m_with)
    return rows, sub

meta_rows = []

# (A) DM2_vs_notDM PFI — within-TCGA driver-anchor sub-cohorts
for ep in ["PFI", "OS", "DSS"]:
    rows, sub = meta_block("DM2_vs_notDM", ep,
                            "TCGA driver subcohorts (BRAF + RAS + NBNR)",
                            cohorts_subset=["TCGA_BRAF_like", "TCGA_RAS_like",
                                            "TCGA_NBNR"])
    print(f"\n[meta DM2_vs_notDM {ep}] cohorts in meta:")
    print(sub[["cohort", "n", "events", "HR", "log_HR", "SE_logHR"]])
    meta_rows.extend(rows)

# (B) DM1_vs_DM2 PFI / OS within-TCGA driver-anchor sub-cohorts
for ep in ["PFI", "OS", "DSS"]:
    rows, sub = meta_block("DM1_vs_DM2", ep,
                            "TCGA driver subcohorts (BRAF + RAS + NBNR)",
                            cohorts_subset=["TCGA_BRAF_like", "TCGA_RAS_like",
                                            "TCGA_NBNR"])
    print(f"\n[meta DM1_vs_DM2 {ep}] cohorts in meta:")
    print(sub[["cohort", "n", "events", "HR", "log_HR", "SE_logHR"]])
    meta_rows.extend(rows)

# (C) DM1xTERT comutant OS — TCGA + Landa (unadjusted, same model spec)
rows, sub = meta_block(
    "DM1xTERT_comutant_vs_others_UNADJ", "OS",
    "TCGA_full_primary + Landa_GSE76039 (independent cohorts, OS unadj)",
    cohorts_subset=["TCGA_full_primary_S3", "Landa_GSE76039_PDTC_ATC"])
print(f"\n[meta DM1xTERT comutant OS UNADJ] cohorts in meta:")
print(sub[["cohort", "n", "events", "HR", "log_HR", "SE_logHR"]])
meta_rows.extend(rows)

# (D) DM1xTERT comutant OS — TCGA + Landa (adjusted within-cohort covariates)
rows, sub = meta_block(
    "DM1xTERT_comutant_vs_others_ADJ", "OS",
    "TCGA_full_primary(adj driver+age+stage+sex) + Landa(adj is_ATC+age)",
    cohorts_subset=["TCGA_full_primary_S3", "Landa_GSE76039_PDTC_ATC"])
print(f"\n[meta DM1xTERT comutant OS ADJ] cohorts in meta:")
print(sub[["cohort", "n", "events", "HR", "log_HR", "SE_logHR"]])
meta_rows.extend(rows)

# (E) DM1 vs (DM2-or-notDM1-proxy) — TCGA driver + Landa proxy
# DM1_vs_DM2 in TCGA + DM1_vs_notDM1_proxy in Landa = different reference, so NOT
# pooled together. We document but don't pool.

meta_df = pd.DataFrame(meta_rows)
meta_df.to_csv(OUT / "r1_pooled_meta.tsv", sep="\t", index=False)
print(f"\n[wrote] r1_pooled_meta.tsv ({len(meta_df)} rows)")
cols = ["contrast", "endpoint", "scope", "method", "k", "pooled_HR", "ci_lo",
        "ci_hi", "p", "I2", "Q", "Q_p", "tau2"]
print("\n=== POOLED META ===")
print(meta_df[cols].to_string(index=False, max_colwidth=44))

# =============================================================================
# Leave-one-out
# =============================================================================
loo_rows = []
def run_loo(contrast, endpoint, cohorts_subset, scope):
    sub = pc[
        (pc["contrast"] == contrast)
        & (pc["endpoint"] == endpoint)
        & (pc["cohort"].isin(cohorts_subset))
        & (pc["log_HR"].notna())
    ].reset_index(drop=True)
    rows = []
    for i in range(len(sub)):
        drop = sub.iloc[i]["cohort"]
        keep = sub.drop(i)
        if len(keep) == 1:
            r = keep.iloc[0]
            rows.append({"contrast": contrast, "endpoint": endpoint,
                         "scope": scope, "left_out": drop, "method": "single",
                         "k": 1, "pooled_logHR": r["log_HR"],
                         "pooled_SE": r["SE_logHR"], "pooled_HR": r["HR"],
                         "ci_lo": r["ci_lo"], "ci_hi": r["ci_hi"],
                         "p": r["p"], "I2": np.nan, "tau2": np.nan,
                         "Q": np.nan, "Q_df": 0, "Q_p": np.nan})
            continue
        if len(keep) >= 2:
            m = meta_inverse_variance(keep["log_HR"].values,
                                        keep["SE_logHR"].values, method="REML")
            rows.append({"contrast": contrast, "endpoint": endpoint,
                         "scope": scope, "left_out": drop, **m})
    return rows

loo_rows.extend(run_loo("DM2_vs_notDM", "PFI",
                         ["TCGA_BRAF_like", "TCGA_RAS_like", "TCGA_NBNR"],
                         "TCGA driver subcohorts"))
loo_rows.extend(run_loo("DM2_vs_notDM", "OS",
                         ["TCGA_BRAF_like", "TCGA_RAS_like", "TCGA_NBNR"],
                         "TCGA driver subcohorts"))
loo_rows.extend(run_loo("DM1xTERT_comutant_vs_others_UNADJ", "OS",
                         ["TCGA_full_primary_S3", "Landa_GSE76039_PDTC_ATC"],
                         "TCGA_full + Landa"))
loo_rows.extend(run_loo("DM1xTERT_comutant_vs_others_ADJ", "OS",
                         ["TCGA_full_primary_S3", "Landa_GSE76039_PDTC_ATC"],
                         "TCGA_full + Landa"))

loo_df = pd.DataFrame(loo_rows)
loo_df.to_csv(OUT / "r1_leave_one_out.tsv", sep="\t", index=False)
print(f"\n[wrote] r1_leave_one_out.tsv ({len(loo_df)} rows)")

# =============================================================================
# Forest data
# =============================================================================
forest = []
def add_panel(panel_label, contrast, endpoint, cohorts):
    for _, r in pc[
        (pc["contrast"] == contrast)
        & (pc["endpoint"] == endpoint)
        & (pc["cohort"].isin(cohorts))
    ].iterrows():
        forest.append({"forest_panel": panel_label, "label": r["cohort"],
                        "kind": "cohort", "n": r["n"], "events": r["events"],
                        "HR": r["HR"], "ci_lo": r["ci_lo"],
                        "ci_hi": r["ci_hi"], "p": r["p"], "I2": np.nan,
                        "Q_p": np.nan})
    for _, r in meta_df[
        (meta_df["contrast"] == contrast)
        & (meta_df["endpoint"] == endpoint)
    ].iterrows():
        forest.append({"forest_panel": panel_label,
                        "label": f"Pooled ({r['method']})",
                        "kind": "pooled", "n": np.nan, "events": np.nan,
                        "HR": r["pooled_HR"], "ci_lo": r["ci_lo"],
                        "ci_hi": r["ci_hi"], "p": r["p"], "I2": r["I2"],
                        "Q_p": r["Q_p"]})

add_panel("DM2_vs_notDM_PFI_TCGA_drivers", "DM2_vs_notDM", "PFI",
          ["TCGA_BRAF_like", "TCGA_RAS_like", "TCGA_NBNR",
           "TCGA_BRAF_cPTC_anchor", "TCGA_full_primary_S3"])
add_panel("DM2_vs_notDM_OS_TCGA_drivers", "DM2_vs_notDM", "OS",
          ["TCGA_BRAF_like", "TCGA_RAS_like", "TCGA_NBNR",
           "TCGA_BRAF_cPTC_anchor", "TCGA_full_primary_S3"])
add_panel("DM1xTERT_OS_pooled", "DM1xTERT_comutant_vs_others_UNADJ", "OS",
          ["TCGA_full_primary_S3", "Landa_GSE76039_PDTC_ATC"])
add_panel("DM1xTERT_OS_ADJ_pooled", "DM1xTERT_comutant_vs_others_ADJ", "OS",
          ["TCGA_full_primary_S3", "Landa_GSE76039_PDTC_ATC"])
add_panel("DM1_vs_DM2_PFI_TCGA_drivers", "DM1_vs_DM2", "PFI",
          ["TCGA_BRAF_like", "TCGA_RAS_like", "TCGA_NBNR",
           "TCGA_BRAF_cPTC_anchor", "TCGA_full_primary_S3"])

forest_df = pd.DataFrame(forest)
forest_df.to_csv(OUT / "r1_forest_data.tsv", sep="\t", index=False)
print(f"\n[wrote] r1_forest_data.tsv ({len(forest_df)} rows)")

# =============================================================================
# Headline JSON
# =============================================================================
def get_pooled(contrast, endpoint, scope_substr, method):
    m = meta_df[
        (meta_df["contrast"] == contrast)
        & (meta_df["endpoint"] == endpoint)
        & (meta_df["scope"].str.contains(scope_substr, regex=False))
        & (meta_df["method"] == method)
    ]
    if len(m) == 0:
        return None
    r = m.iloc[0]
    return {
        "method": method, "k": int(r["k"]),
        "HR": float(r["pooled_HR"]),
        "ci_lo": float(r["ci_lo"]), "ci_hi": float(r["ci_hi"]),
        "p": float(r["p"]), "I2": float(r["I2"]),
        "Q_p": float(r["Q_p"]), "tau2": float(r["tau2"]),
        "cohorts": str(r["cohorts_in_meta"]),
    }

headline = {
    "DM2_vs_notDM_PFI_TCGA_driver_meta_REML": get_pooled(
        "DM2_vs_notDM", "PFI", "TCGA driver subcohorts", "REML"),
    "DM2_vs_notDM_PFI_TCGA_driver_meta_FE": get_pooled(
        "DM2_vs_notDM", "PFI", "TCGA driver subcohorts", "FE"),
    "DM2_vs_notDM_OS_TCGA_driver_meta_REML": get_pooled(
        "DM2_vs_notDM", "OS", "TCGA driver subcohorts", "REML"),
    "DM2_vs_notDM_DSS_TCGA_driver_meta_REML": get_pooled(
        "DM2_vs_notDM", "DSS", "TCGA driver subcohorts", "REML"),
    "DM1_vs_DM2_PFI_TCGA_driver_meta_REML": get_pooled(
        "DM1_vs_DM2", "PFI", "TCGA driver subcohorts", "REML"),
    "DM1xTERT_comutant_OS_TCGA_Landa_REML_UNADJ": get_pooled(
        "DM1xTERT_comutant_vs_others_UNADJ", "OS",
        "TCGA_full_primary + Landa_GSE76039", "REML"),
    "DM1xTERT_comutant_OS_TCGA_Landa_FE_UNADJ": get_pooled(
        "DM1xTERT_comutant_vs_others_UNADJ", "OS",
        "TCGA_full_primary + Landa_GSE76039", "FE"),
    "DM1xTERT_comutant_OS_TCGA_Landa_REML_ADJ": get_pooled(
        "DM1xTERT_comutant_vs_others_ADJ", "OS",
        "TCGA_full_primary(adj driver+age+stage+sex)", "REML"),
}

with open(OUT / "r1_headline.json", "w") as f:
    json.dump(headline, f, indent=2, default=str)
print("\n=== HEADLINE JSON ===")
print(json.dumps(headline, indent=2, default=str))

# Also print anchor row for context
anchor = pc[(pc["cohort"] == "TCGA_BRAF_cPTC_anchor")
             & (pc["contrast"] == "DM2_vs_notDM")
             & (pc["endpoint"] == "PFI")]
print("\n=== H6 ANCHOR (TCGA BRAF-cPTC DM2-vs-notDM PFI) ===")
print(anchor[["cohort", "n", "events", "HR", "ci_lo", "ci_hi", "p"]].to_string(index=False))
