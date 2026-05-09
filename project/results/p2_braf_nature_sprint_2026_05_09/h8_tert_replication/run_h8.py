#!/usr/bin/env python3
"""
H8 — DM1 x TERT promoter mutation replication HARD GATE.

External cohorts attempted:
  1. Lee 2024 GSE213647 (Korean PTC RNA-seq n=632)
       --> NO TERT or BRAF mutation calls released by authors.
           series_matrix Sample_characteristics_ch1: "genotype: NA" for all 632.
           Lee2024_MOESM5.xlsx contains metabolomics + DEG only.
           ==> SKIPPED (documented). Cannot fabricate.
  2. Landa GSE76039 (PDTC/ATC microarray n=37)
       RNA: /data/thca/repo_results/p_landa_2016/scores.tsv (DM1_like_score)
       TERT: /data/thca/repo_results/v17_tert_recovery/v2/FINAL_promoter_records_all_sources.tsv
             (35 samples called via Pozdeyev/cBioPortal MSKCC 2016 dump)
       Survival: /data/thca/repo_results/v17_tert_recovery/v3/v3_thyroid_mskcc_2016_clinical_*.tsv
                 (OS_MONTHS + OS_STATUS, n=117 patients, 47 deceased)
  3. K2 / PRJEB11591 (Yoo 2016 Korean RNA-seq n=180)
       Has BRAF/RAS calls but TERT promoter NOT sequenced (TERT count=0).
       ==> NOT a candidate.
  4. TCGA-THCA full primary tumors (n=513) — re-anchor cohort.
       Per-task instruction: re-run H6's DM1xTERT NOT restricted to BRAF-cPTC,
       with driver_anchor as covariate.

Output:
  h8_external_replication.tsv  rows = cohort x contrast, cols = HR/CI/p/n/events
  h8_data_paths.txt            search audit
  H8_REPORT.md                 verdict
"""
import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.exceptions import ConvergenceError, ConvergenceWarning
from scipy import stats

warnings.filterwarnings("ignore", category=ConvergenceWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=UserWarning)

OUT = Path("/home/seungho/personal/THCA_data_analysis/project/results/p2_braf_nature_sprint_2026_05_09/h8_tert_replication")
OUT.mkdir(parents=True, exist_ok=True)

H6_MERGED = "/home/seungho/personal/THCA_data_analysis/project/results/p2_braf_nature_sprint_2026_05_09/h6_survival/h6_merged_clinical.tsv"
LANDA_SCORES = "/data/thca/repo_results/p_landa_2016/scores.tsv"
LANDA_META = "/data/thca/repo_results/p_landa_2016/sample_metadata.tsv"
MSKCC_CLIN_PT = "/data/thca/repo_results/v17_tert_recovery/v3/v3_thyroid_mskcc_2016_clinical_patient.tsv"
MSKCC_CLIN_SAMP = "/data/thca/repo_results/v17_tert_recovery/v3/v3_thyroid_mskcc_2016_clinical_sample.tsv"
TERT_RECORDS = "/data/thca/repo_results/v17_tert_recovery/v2/FINAL_promoter_records_all_sources.tsv"

results = []  # rows for h8_external_replication.tsv

def cox_run(df, time_col, event_col, covariates, label, cohort, contrast, endpoint, strata=None):
    """Fit Cox PH with given covariates. Returns dict row for the table."""
    sub = df[covariates + [time_col, event_col]].dropna().copy()
    n = len(sub)
    events = int(sub[event_col].sum())
    row = dict(
        cohort=cohort, contrast=contrast, endpoint=endpoint,
        n=n, events=events,
        HR=np.nan, ci_lo=np.nan, ci_hi=np.nan, p=np.nan,
        primary_var=label, note="",
    )
    if events < 3 or n < 8:
        row["note"] = f"too_few_events_or_n (events={events}, n={n})"
        return row
    # collinearity guard: drop covariates with zero variance in subset
    keep = [c for c in covariates if sub[c].nunique(dropna=True) > 1]
    if label not in keep:
        row["note"] = f"primary_var_no_variance"
        return row
    try:
        cph = CoxPHFitter(penalizer=0.01)
        cph.fit(sub[keep + [time_col, event_col]], duration_col=time_col,
                event_col=event_col, strata=strata)
        s = cph.summary.loc[label]
        row["HR"] = float(s["exp(coef)"])
        row["ci_lo"] = float(s["exp(coef) lower 95%"])
        row["ci_hi"] = float(s["exp(coef) upper 95%"])
        row["p"] = float(s["p"])
    except (ConvergenceError, np.linalg.LinAlgError, ValueError) as e:
        row["note"] = f"cox_failed: {type(e).__name__}"
    return row


# =============================================================================
# COHORT 1: Lee 2024 GSE213647 -> SKIP, document reason
# =============================================================================
results.append(dict(
    cohort="Lee2024_GSE213647_Korean_n632",
    contrast="DM1xTERT",
    endpoint="(none)",
    n=0, events=0,
    HR=np.nan, ci_lo=np.nan, ci_hi=np.nan, p=np.nan,
    primary_var="(skipped)",
    note="GSE213647 series_matrix Sample_characteristics_ch1 'genotype: NA' for all 632; "
         "Lee2024_MOESM5.xlsx supp = metabolomics+DEG only; no TERT/BRAF call file released. "
         "No survival in GEO.",
))

# =============================================================================
# COHORT 2: Landa GSE76039 + MSKCC 2016 clinical
# =============================================================================
landa_scores = pd.read_csv(LANDA_SCORES, sep="\t")
landa_meta = pd.read_csv(LANDA_META, sep="\t")
clin_samp = pd.read_csv(MSKCC_CLIN_SAMP, sep="\t")
clin_pt = pd.read_csv(MSKCC_CLIN_PT, sep="\t")

# Landa GEO sample_id (GSM*) -> title (s_JF_thy_xxx_P) -> MSKCC SAMPLE_ID
landa = landa_scores.merge(landa_meta[["sample_id", "title", "histology", "disease_group"]],
                            on="sample_id", how="left", suffixes=("", "_meta"))
landa = landa.merge(
    clin_samp[["PATIENT_ID", "SAMPLE_ID", "CANCER_TYPE_DETAILED"]].rename(
        columns={"SAMPLE_ID": "title"}),
    on="title", how="left",
)
landa = landa.merge(
    clin_pt[["PATIENT_ID", "OS_STATUS", "OS_MONTHS", "AGE", "PDTC_DEFINITION"]],
    on="PATIENT_ID", how="left",
)

# TERT promoter mutation: pull MSKCC-source records, restrict to TERT (chr5:1295228 promoter)
tert_rec = pd.read_csv(TERT_RECORDS, sep="\t", low_memory=False)
# Promoter mutation = chr5 + 1295228..1295250 region OR variant_class == 5'Flank/5'UTR
# These records are pre-filtered to promoter region by upstream pipeline.
mskcc_tert = tert_rec[tert_rec["__source"] == "S4_pozdeyev2024"].copy()
landa_tert_samples = set(mskcc_tert["sample_barcode"].unique())
landa["TERT"] = landa["title"].isin(landa_tert_samples).astype(int)

# DM1_like_score: positive = DM1-like (per Landa scores file convention)
# Define DM1 (binary) as DM1_like_score > 0 vs <= 0
landa["DM1_bin"] = (landa["DM1_like_score"] > 0).astype(int)

# OS event
def parse_os(x):
    if pd.isna(x):
        return np.nan
    return 1 if str(x).startswith("1") else 0

landa["OS_event"] = landa["OS_STATUS"].apply(parse_os)
landa["OS_months"] = pd.to_numeric(landa["OS_MONTHS"], errors="coerce")
landa["AGE_num"] = pd.to_numeric(landa["AGE"], errors="coerce")
# Encode disease_group as ATC=1 vs PDTC=0 (covariate)
landa["is_ATC"] = (landa["disease_group"] == "ATC").astype(int)

print(f"[Landa] merged n={len(landa)}, with OS={landa['OS_event'].notna().sum()}, "
      f"events={int(landa['OS_event'].fillna(0).sum())}")
print(f"[Landa] TERT mutated: {int(landa['TERT'].sum())} / {len(landa)}; "
      f"DM1_bin: {int(landa['DM1_bin'].sum())}")

# crosstab DM1 x TERT
print("\n[Landa] DM1_bin x TERT crosstab (with OS):")
sub_land = landa.dropna(subset=["OS_event", "OS_months"])
print(pd.crosstab(sub_land["DM1_bin"], sub_land["TERT"], margins=True))
print("\n[Landa] event rates (OS) per DM1xTERT cell:")
print(sub_land.groupby(["DM1_bin", "TERT"])["OS_event"].agg(["sum", "count"]))

landa_save = landa[["sample_id", "title", "PATIENT_ID", "histology", "disease_group",
                    "DM1_like_score", "DM1_bin", "TERT", "OS_event", "OS_months",
                    "AGE_num", "is_ATC"]].copy()
landa_save.to_csv(OUT / "h8_landa_merged.tsv", sep="\t", index=False)

# Cox: DM1 x TERT interaction (OS) - Landa
sub_l = landa.dropna(subset=["OS_event", "OS_months", "DM1_bin", "TERT"]).copy()
sub_l["DM1xTERT"] = sub_l["DM1_bin"] * sub_l["TERT"]

# Model 1: DM1 main effect alone
results.append({**cox_run(sub_l, "OS_months", "OS_event", ["DM1_bin"],
                          "DM1_bin", "Landa_GSE76039_PDTC+ATC", "DM1_bin alone",
                          "OS")})
# Model 2: TERT main effect alone
results.append({**cox_run(sub_l, "OS_months", "OS_event", ["TERT"],
                          "TERT", "Landa_GSE76039_PDTC+ATC", "TERT alone", "OS")})
# Model 3: DM1xTERT interaction (full model w/ both main effects)
results.append({**cox_run(sub_l, "OS_months", "OS_event",
                          ["DM1_bin", "TERT", "DM1xTERT"],
                          "DM1xTERT", "Landa_GSE76039_PDTC+ATC",
                          "DM1xTERT interaction (with main effects)", "OS")})
# Model 4: DM1xTERT adjusted for is_ATC + AGE
results.append({**cox_run(sub_l, "OS_months", "OS_event",
                          ["DM1_bin", "TERT", "DM1xTERT", "is_ATC", "AGE_num"],
                          "DM1xTERT", "Landa_GSE76039_PDTC+ATC",
                          "DM1xTERT adj is_ATC+AGE", "OS")})
# Model 5: 4-group categorical (DM1+/TERT+ vs others) - dichotomous comutant flag
sub_l["comutant"] = ((sub_l["DM1_bin"] == 1) & (sub_l["TERT"] == 1)).astype(int)
results.append({**cox_run(sub_l, "OS_months", "OS_event", ["comutant"],
                          "comutant", "Landa_GSE76039_PDTC+ATC",
                          "DM1+/TERT+ comutant vs others", "OS")})
results.append({**cox_run(sub_l, "OS_months", "OS_event",
                          ["comutant", "is_ATC", "AGE_num"],
                          "comutant", "Landa_GSE76039_PDTC+ATC",
                          "DM1+/TERT+ comutant adj is_ATC+AGE", "OS")})

# =============================================================================
# COHORT 3: TCGA-THCA FULL primary tumor cohort (re-anchor, all drivers)
# =============================================================================
tcga = pd.read_csv(H6_MERGED, sep="\t")
print(f"\n[TCGA] full cohort n={len(tcga)}")
print(tcga[["dm", "molecular_subtype", "histology_subtype",
            "tert_promoter_integrated"]].head())

# TERT mutation flag from integrated column
def tert_flag(x):
    if pd.isna(x):
        return np.nan
    s = str(x).lower()
    return 1 if s in ("mutated", "1", "true") else 0

tcga["TERT"] = tcga["tert_promoter_integrated"].apply(tert_flag)
# DM1 binary - per H6 convention: dm == "DM1" -> 1, else 0
tcga["DM1_bin"] = (tcga["dm"] == "DM1").astype(int)
tcga["DM2_bin"] = (tcga["dm"] == "DM2").astype(int)
tcga["DM1xTERT"] = tcga["DM1_bin"] * tcga["TERT"]

# Driver anchor: BRAF_like / RAS_like / NBNR (from molecular_subtype)
# Encode as dummies
tcga["BRAF_like"] = (tcga["molecular_subtype"] == "BRAF_like").astype(int)
tcga["RAS_like"] = (tcga["molecular_subtype"] == "RAS_like").astype(int)
# NBNR = baseline reference

print("\n[TCGA] DM1xTERT crosstab (full cohort):")
print(pd.crosstab(tcga["DM1_bin"], tcga["TERT"], margins=True))
print("\n[TCGA] PFI events per DM1xTERT cell (full cohort):")
sub_t = tcga.dropna(subset=["PFI", "PFI.time", "DM1_bin", "TERT"])
print(sub_t.groupby(["DM1_bin", "TERT"])["PFI"].agg(["sum", "count"]))

# Endpoint helpers
endpoints = [("PFI", "PFI", "PFI.time"), ("OS", "OS", "OS.time"),
             ("DSS", "DSS", "DSS.time")]

for ep_label, ev_col, time_col in endpoints:
    sub_te = tcga.dropna(subset=[ev_col, time_col, "DM1_bin", "TERT"]).copy()
    sub_te["DM1xTERT"] = sub_te["DM1_bin"] * sub_te["TERT"]
    # full cohort: DM1xTERT alone
    results.append({**cox_run(sub_te, time_col, ev_col,
                              ["DM1_bin", "TERT", "DM1xTERT"],
                              "DM1xTERT",
                              "TCGA-THCA_full_primary",
                              "DM1xTERT (full cohort, no driver covar)",
                              ep_label)})
    # adjusted for driver_anchor + age + stage + sex
    covs = ["DM1_bin", "TERT", "DM1xTERT", "BRAF_like", "RAS_like",
            "age_yr", "adv_stage", "sex_m"]
    results.append({**cox_run(sub_te, time_col, ev_col, covs,
                              "DM1xTERT", "TCGA-THCA_full_primary",
                              "DM1xTERT adj driver+age+stage+sex", ep_label)})
    # Comutant (DM1+/TERT+) vs all others, adjusted
    sub_te["comutant"] = ((sub_te["DM1_bin"] == 1) & (sub_te["TERT"] == 1)).astype(int)
    results.append({**cox_run(sub_te, time_col, ev_col,
                              ["comutant", "BRAF_like", "RAS_like",
                               "age_yr", "adv_stage", "sex_m"],
                              "comutant", "TCGA-THCA_full_primary",
                              "DM1+/TERT+ comutant adj driver+age+stage+sex",
                              ep_label)})

# Per-driver stratum: does DM1xTERT replicate within RAS-like / NBNR?
for stratum_name, mask in [
    ("BRAF_like_only", tcga["molecular_subtype"] == "BRAF_like"),
    ("RAS_like_only", tcga["molecular_subtype"] == "RAS_like"),
    ("NBNR_only", tcga["molecular_subtype"] == "NBNR"),
]:
    for ep_label, ev_col, time_col in endpoints:
        sub_s = tcga[mask].dropna(subset=[ev_col, time_col, "DM1_bin", "TERT"]).copy()
        sub_s["DM1xTERT"] = sub_s["DM1_bin"] * sub_s["TERT"]
        results.append({**cox_run(sub_s, time_col, ev_col,
                                  ["DM1_bin", "TERT", "DM1xTERT"],
                                  "DM1xTERT", f"TCGA-THCA_{stratum_name}",
                                  "DM1xTERT (per-driver stratum)", ep_label)})

# =============================================================================
# Save
# =============================================================================
res_df = pd.DataFrame(results)
res_df.to_csv(OUT / "h8_external_replication.tsv", sep="\t", index=False)
print(f"\n[done] wrote {OUT / 'h8_external_replication.tsv'} ({len(res_df)} rows)")

# Summary print
print("\n=== HEADLINE RESULTS ===")
hl_cols = ["cohort", "contrast", "endpoint", "n", "events", "HR", "ci_lo", "ci_hi", "p", "note"]
print(res_df[hl_cols].to_string(index=False, max_colwidth=60))
