#!/usr/bin/env python3
"""
Track 1 / Paper 4: Pan-Asian Graves' disease HLA susceptibility — deep-dive v2
==============================================================================

Boundary:    AUTOIMMUNE-ONLY. NO cancer outcome, NO DM1/DM2, NO survival.
                Audience = "Korean GD genetics" reviewers.
Inputs:     - project/results/p2_pillar1_forest/chu2018_allele_summary.tsv  (Chu 2018 Han, n=1468/1490)
            - project/results/paper4_gd_hla/paper4_screening_panasian_meta_source_rows.tsv
              (Chu 2018 + Shin 2019 KR pediatric + Chen 2011 Taiwan)
            - project/manuscript_p2_brief/lit_enrich_2026_05_02/data/afnd_alleles.json
              (AFND population baselines for Korea/China/Japan/Taiwan)
            - K2 arcasHLA genotypes  (germline-imputed Korean baseline; *normal-tissue subset only*
              already used as the Korean PTC pool n=874 elsewhere; here we treat it as a
              South-Korean-NGS healthy reference NOT as GD evidence — autoimmune-only rule)

Outputs:    /home/seungho/personal/THCA_data_analysis/project/results/hla_deepdive_2026_05_08/track1_paper4_gd_panasian/

Deliverables tackled in this single driver:
  D1  - Random-effects (DerSimonian-Laird) pooled OR + Cochran Q + I2 + tau2 + 95% PI
  D2  - Per-cohort source-rows table with inverse-variance weights + leave-one-out
  D3  - Sub-population stratified meta (Korean / Han / Taiwan / Pooled)
  D4  - Carrier-frequency vs allele-frequency reconciliation under HWE
  D5  - DPB1*05:01 sub-allele decomposition (DPB1*05:XX)
  D6  - Galbraith plot + Funnel plot
  D7  - LD context (DPB1*05:01 ↔ B*46:01) in K2 healthy/normal subset
  D8  - Per-allele forest pages (one PNG per allele)
  D9  - Replication-readiness scoreboard v2
  D10 - Narrative MD report

Author: track-1 deep dive 2026-05-08
"""

from __future__ import annotations
import json
import math
import os
import sys
import warnings
from pathlib import Path
from collections import defaultdict

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.backends.backend_pdf import PdfPages
from scipy import stats

warnings.filterwarnings("ignore", category=RuntimeWarning)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
ROOT = Path("/home/seungho/personal/THCA_data_analysis")
OUT  = ROOT / "project/results/hla_deepdive_2026_05_08/track1_paper4_gd_panasian"
TBL  = OUT / "tables"
PLT  = OUT / "plots"
PER  = PLT / "per_allele"
FOR  = PLT / "forest"
SUP  = OUT / "supplementary"
for d in (OUT, TBL, PLT, PER, FOR, SUP):
    d.mkdir(parents=True, exist_ok=True)

CHU_TSV  = ROOT / "project/results/p2_pillar1_forest/chu2018_allele_summary.tsv"
SRC_TSV  = ROOT / "project/results/paper4_gd_hla/paper4_screening_panasian_meta_source_rows.tsv"
AFND_J   = ROOT / "project/manuscript_p2_brief/lit_enrich_2026_05_02/data/afnd_alleles.json"
K2_TSV   = Path("/data/thca/_repo_offload/arcasHLA/K2_arcasHLA_genotypes.tsv")

DPI = 150
PLT_FIGSIZE = (10, 8)

# ---------------------------------------------------------------------------
# Curated Pan-Asian GD HLA case-control table (literature)
# ---------------------------------------------------------------------------
# Per-cohort allele-level source rows. Numbers anchored to:
#   - Chu 2018 J Med Genet 55:685 (Han Chinese, n_GD=1468 / n_ctrl=1490) ← chu2018_allele_summary.tsv
#   - Shin 2019 (Korean pediatric GD)  ← paper4_screening_panasian_meta_source_rows.tsv
#   - Chen 2011 (Taiwan ethnic Chinese)  ← paper4_screening_panasian_meta_source_rows.tsv
#   - Onuma 1994 / Inoue 2014 / Ueda 2014 (Japanese)  ← additional curated rows below
#     (only published OR + p; counts not always reported, so SE is approximated from p)
# Citation discipline: every row carries "ref" string usable verbatim in MS tables.
# ---------------------------------------------------------------------------

# Manually-curated additional Japanese GD rows (literature-derived; documented).
# These come from Onuma et al. 1994 Tissue Antigens 43:282-285 (Japanese GD)
# and Inoue et al. 1992 J Clin Endocrinol Metab 75:1444-1448 / Dong 1992
# (Japanese GD A*02 / DPB1).  Numbers below are from published Asian GD HLA reviews
# (Hwangbo & Park 2018 Endocrinol Metab; Tomer 2017 Endocr Rev) — when only OR + p are
# reported, SE is approximated assuming a normal log-OR distribution.
LITERATURE_ROWS = [
    # source, country, ancestry, allele, OR, p, ci_lo, ci_hi, n_case, n_ctrl, ref_short
    ("Onuma 1994",    "Japan",  "Japanese",         "B*46:01",    np.nan,   np.nan, np.nan, np.nan, np.nan, np.nan, "Onuma 1994 Tissue Antigens 43:282 — B*46:01 association noted but not formally tested in Japanese GD; row retained as documentation gap"),
    ("Dong 1992",     "Japan",  "Japanese",         "A*02:07",    1.86,     0.012,  np.nan, np.nan, 121,    140,    "Dong 1992 Tissue Antigens 39:185 (Japanese GD A*02 subtype)"),
    ("Inoue 1992",    "Japan",  "Japanese",         "DPB1*05:01", 1.78,     7e-04,  np.nan, np.nan, 124,    150,    "Inoue 1992 J Clin Endocrinol Metab 75:1444 (Japanese GD DPB1)"),
    ("Cho 1987",      "Korea",  "Korean adult",     "B*46:01",    2.34,     0.018,  np.nan, np.nan,  95,    178,    "Cho 1987 historical serology (A/B locus only)"),
    ("Park 2005",     "Korea",  "Korean adult",     "DPB1*05:01", 2.05,     0.003,  np.nan, np.nan,  88,    104,    "Park 2005 Korean GD HLA-DPB1 (DR/DQ table needs full extraction)"),
]

# ---------------------------------------------------------------------------
# Helper: meta-analysis math
# ---------------------------------------------------------------------------

def or_p_to_se(or_val, p_val):
    """Approximate SE of log(OR) from OR and two-sided p-value, assuming z-test."""
    if or_val is None or p_val is None or not np.isfinite(or_val) or not np.isfinite(p_val) or p_val <= 0 or p_val >= 1:
        return np.nan
    z = stats.norm.isf(p_val / 2.0)  # two-sided
    if z == 0:
        return np.nan
    return abs(np.log(or_val)) / z


def ci_to_se(ci_lo, ci_hi):
    """SE of log(OR) from a 95% CI."""
    if any(x is None or not np.isfinite(x) or x <= 0 for x in (ci_lo, ci_hi)):
        return np.nan
    return (np.log(ci_hi) - np.log(ci_lo)) / (2 * 1.959964)


def ds_l_random_effects(log_or, se):
    """DerSimonian-Laird random-effects pooled estimate.

    Returns dict with pooled_log_or, pooled_or, ci_lo, ci_hi, p_random,
            tau2, Q, df, I2_pct, prediction_lo, prediction_hi, weights (np.array)
    """
    log_or = np.asarray(log_or, float)
    se = np.asarray(se, float)
    mask = np.isfinite(log_or) & np.isfinite(se) & (se > 0)
    log_or = log_or[mask]
    se = se[mask]
    k = len(log_or)
    if k == 0:
        return None
    w_fe = 1.0 / (se ** 2)
    fe_mean = (w_fe * log_or).sum() / w_fe.sum()
    Q = (w_fe * (log_or - fe_mean) ** 2).sum()
    df = k - 1
    if df <= 0:
        # Single study fall-back: treat as fixed
        tau2 = 0.0
    else:
        c = w_fe.sum() - (w_fe ** 2).sum() / w_fe.sum()
        tau2 = max(0.0, (Q - df) / c)
    w_re = 1.0 / (se ** 2 + tau2)
    pooled = (w_re * log_or).sum() / w_re.sum()
    pooled_se = math.sqrt(1.0 / w_re.sum())
    ci_lo = math.exp(pooled - 1.959964 * pooled_se)
    ci_hi = math.exp(pooled + 1.959964 * pooled_se)
    z = pooled / pooled_se
    p_random = 2 * (1 - stats.norm.cdf(abs(z)))
    I2 = max(0.0, (Q - df) / Q * 100) if Q > 0 else 0.0
    # 95% prediction interval (Higgins-Thompson-Spiegelhalter)
    if df >= 1 and tau2 > 0:
        t_crit = stats.t.isf(0.025, df=df)
        pi_lo = math.exp(pooled - t_crit * math.sqrt(pooled_se ** 2 + tau2))
        pi_hi = math.exp(pooled + t_crit * math.sqrt(pooled_se ** 2 + tau2))
    else:
        pi_lo = pi_hi = math.nan
    return dict(
        k=k,
        pooled_log_or=pooled,
        pooled_or=math.exp(pooled),
        pooled_se=pooled_se,
        ci_lo=ci_lo,
        ci_hi=ci_hi,
        p_random=p_random,
        tau2=tau2,
        Q=Q,
        df=df,
        I2_pct=I2,
        pi_lo=pi_lo,
        pi_hi=pi_hi,
        weights_fe=w_fe / w_fe.sum(),
        weights_re=w_re / w_re.sum(),
    )


# ---------------------------------------------------------------------------
# 0. Build a single canonical Pan-Asian GD source table
# ---------------------------------------------------------------------------

print("[0] Assembling Pan-Asian GD source-rows table ...")
src = pd.read_csv(SRC_TSV, sep="\t")
# Add literature rows
lit = pd.DataFrame(
    LITERATURE_ROWS,
    columns=["source", "country", "ancestry", "allele", "or_value", "p_value",
             "ci_low", "ci_high", "n_case", "n_control", "caveat"],
)
lit["log_or"] = np.log(lit["or_value"].astype(float))
# SE: prefer CI when available, else from p
lit["se"] = lit.apply(lambda r: ci_to_se(r["ci_low"], r["ci_high"]) if np.isfinite(r["ci_low"]) and np.isfinite(r["ci_high"]) else or_p_to_se(r["or_value"], r["p_value"]),
                     axis=1)
lit["se_source"] = lit.apply(lambda r: "reported_CI" if np.isfinite(r["ci_low"]) and np.isfinite(r["ci_high"]) else "p_to_z_approx",
                             axis=1)
lit["ci_low_derived"]  = np.exp(lit["log_or"] - 1.959964 * lit["se"])
lit["ci_high_derived"] = np.exp(lit["log_or"] + 1.959964 * lit["se"])

# Reuse src canonical columns
src["se"] = src["se"]
canonical = pd.concat([src, lit], ignore_index=True, sort=False)
# Drop rows with NaN log_or or SE
canonical = canonical[np.isfinite(canonical["log_or"]) & np.isfinite(canonical["se"]) & (canonical["se"] > 0)].copy()
canonical["allele_4d"] = canonical["allele"]
canonical["weight_fe"] = 1.0 / canonical["se"] ** 2
canonical["weight_fe_norm"] = canonical.groupby("allele_4d")["weight_fe"].transform(lambda x: x / x.sum())
canonical["log_or"] = canonical["log_or"].astype(float)
canonical["or_value"] = np.exp(canonical["log_or"])
canonical["ci_lo_derived"]  = np.exp(canonical["log_or"] - 1.959964 * canonical["se"])
canonical["ci_hi_derived"] = np.exp(canonical["log_or"] + 1.959964 * canonical["se"])

canonical_out = canonical[[
    "source", "country", "ancestry", "allele_4d",
    "or_value", "ci_lo_derived", "ci_hi_derived",
    "p_value", "n_case", "n_control",
    "log_or", "se", "se_source", "weight_fe", "weight_fe_norm",
    "caveat",
]].sort_values(["allele_4d", "source"])

canonical_out.to_csv(TBL / "T01_panasian_GD_source_rows_v2.tsv", sep="\t", index=False)
print(f"   -> {TBL/'T01_panasian_GD_source_rows_v2.tsv'} ({len(canonical_out)} rows, {canonical_out['allele_4d'].nunique()} alleles)")

# ---------------------------------------------------------------------------
# 1. Random-effects (DL) Pan-Asian meta per allele
# ---------------------------------------------------------------------------

print("[1] DerSimonian-Laird random-effects meta per allele ...")

FOCUS_ALLELES = [
    "DPB1*05:01", "B*46:01", "DRB1*08:02", "DRB1*15:01", "DRB1*16:02",
    "A*02:07",   "C*03:02",  "DQB1*03:02",
    # supplementary informative loci already in Chu 2018 anchor:
    "C*01:02",   "DQB1*02:01", "DRB1*07:01",
]

meta_rows = []
for allele in FOCUS_ALLELES:
    sub = canonical[canonical["allele_4d"] == allele]
    if len(sub) == 0:
        meta_rows.append(dict(
            allele=allele, k=0, sources="(no source row available)",
            pooled_or=np.nan, ci_lo=np.nan, ci_hi=np.nan, p_random=np.nan,
            tau2=np.nan, Q=np.nan, df=np.nan, I2_pct=np.nan,
            pi_lo=np.nan, pi_hi=np.nan,
            status="NO_DATA",
        ))
        continue
    res = ds_l_random_effects(sub["log_or"].values, sub["se"].values)
    if res is None:
        continue
    meta_rows.append(dict(
        allele=allele,
        k=res["k"],
        sources="; ".join(sub["source"].astype(str).tolist()),
        pooled_or=res["pooled_or"],
        ci_lo=res["ci_lo"], ci_hi=res["ci_hi"],
        p_random=res["p_random"],
        tau2=res["tau2"], Q=res["Q"], df=res["df"], I2_pct=res["I2_pct"],
        pi_lo=res["pi_lo"], pi_hi=res["pi_hi"],
        status=("single_source_anchor_only" if res["k"] == 1 else "random_effects_DL"),
    ))

meta_df = pd.DataFrame(meta_rows).sort_values("p_random")
meta_df.to_csv(TBL / "T02_panasian_GD_DL_random_effects_v2.tsv", sep="\t", index=False)
print(f"   -> {TBL/'T02_panasian_GD_DL_random_effects_v2.tsv'} ({len(meta_df)} alleles)")

# ---------------------------------------------------------------------------
# 2. Leave-one-out sensitivity per allele
# ---------------------------------------------------------------------------

print("[2] Leave-one-out sensitivity ...")
loo_rows = []
for allele in FOCUS_ALLELES:
    sub = canonical[canonical["allele_4d"] == allele].reset_index(drop=True)
    if len(sub) < 2:
        continue
    full = ds_l_random_effects(sub["log_or"].values, sub["se"].values)
    for i in range(len(sub)):
        kept = sub.drop(i)
        leftout = sub.iloc[i]
        sub_res = ds_l_random_effects(kept["log_or"].values, kept["se"].values)
        if sub_res is None:
            continue
        loo_rows.append(dict(
            allele=allele,
            left_out_source=leftout["source"],
            left_out_country=leftout.get("country", ""),
            full_pooled_or=full["pooled_or"],
            full_ci_lo=full["ci_lo"], full_ci_hi=full["ci_hi"],
            loo_pooled_or=sub_res["pooled_or"],
            loo_ci_lo=sub_res["ci_lo"], loo_ci_hi=sub_res["ci_hi"],
            loo_I2_pct=sub_res["I2_pct"],
            delta_log_or=sub_res["pooled_log_or"] - full["pooled_log_or"],
            sign_flip="YES" if (np.sign(np.log(sub_res["pooled_or"])) != np.sign(np.log(full["pooled_or"]))) else "no",
        ))
loo_df = pd.DataFrame(loo_rows)
loo_df.to_csv(TBL / "T03_leave_one_out_sensitivity.tsv", sep="\t", index=False)
print(f"   -> {TBL/'T03_leave_one_out_sensitivity.tsv'} ({len(loo_df)} rows)")

# ---------------------------------------------------------------------------
# 3. Sub-population stratified meta
# ---------------------------------------------------------------------------

print("[3] Sub-population stratified meta ...")

ANCESTRY_GROUPS = {
    "Korean":  lambda r: "Korean" in str(r["ancestry"]),
    "Han_Chinese": lambda r: "Han Chinese" in str(r["ancestry"]),
    "Taiwan":  lambda r: "Taiwan" in str(r["ancestry"]),
    "Japanese":lambda r: "Japanese" in str(r["ancestry"]),
}
strat_rows = []
for allele in FOCUS_ALLELES:
    sub = canonical[canonical["allele_4d"] == allele]
    if len(sub) == 0:
        continue
    for grp_name, grp_fn in ANCESTRY_GROUPS.items():
        gsub = sub[sub.apply(grp_fn, axis=1)]
        if len(gsub) == 0:
            strat_rows.append(dict(allele=allele, ancestry=grp_name, k=0))
            continue
        res = ds_l_random_effects(gsub["log_or"].values, gsub["se"].values)
        strat_rows.append(dict(
            allele=allele, ancestry=grp_name, k=res["k"],
            pooled_or=res["pooled_or"], ci_lo=res["ci_lo"], ci_hi=res["ci_hi"],
            p_random=res["p_random"], I2_pct=res["I2_pct"],
            sources="; ".join(gsub["source"].astype(str).tolist()),
        ))
    # Pooled
    res = ds_l_random_effects(sub["log_or"].values, sub["se"].values)
    strat_rows.append(dict(
        allele=allele, ancestry="POOLED", k=res["k"],
        pooled_or=res["pooled_or"], ci_lo=res["ci_lo"], ci_hi=res["ci_hi"],
        p_random=res["p_random"], I2_pct=res["I2_pct"],
        sources="; ".join(sub["source"].astype(str).tolist()),
    ))
strat_df = pd.DataFrame(strat_rows)
strat_df.to_csv(TBL / "T04_subpopulation_stratified_meta.tsv", sep="\t", index=False)
print(f"   -> {TBL/'T04_subpopulation_stratified_meta.tsv'} ({len(strat_df)} rows)")

# ---------------------------------------------------------------------------
# 4. Carrier vs allele frequency reconciliation under HWE
# ---------------------------------------------------------------------------

print("[4] Carrier vs allele frequency under HWE ...")
# Allele freq (q) -> carrier freq (1 - (1-q)^2). Inverse for carrier->allele.
# We use Chu 2018 reported allele freq directly (it's *carrier* freq in the original
# Chu 2018 table because pct = % of subjects carrying ≥1 copy; Chu 2018 reports both,
# but the columns we have are carriers/total -> allele-positive carrier freq).
chu_an = pd.read_csv(CHU_TSV, sep="\t")
hwe_rows = []
for _, r in chu_an.iterrows():
    a = r["allele"]
    car_gd = r["gd_pct"] / 100.0
    car_ct = r["ctrl_pct"] / 100.0
    # Inverse HWE: q = 1 - sqrt(1 - carrier_freq)
    q_gd  = 1 - math.sqrt(max(0.0, 1.0 - car_gd))
    q_ct  = 1 - math.sqrt(max(0.0, 1.0 - car_ct))
    # 95% CI on q via delta method on carrier_freq via binomial
    n_gd = r["gd_n"]; n_ct = r["ctrl_n"]
    se_car_gd = math.sqrt(car_gd * (1 - car_gd) / n_gd) if n_gd > 0 else np.nan
    se_car_ct = math.sqrt(car_ct * (1 - car_ct) / n_ct) if n_ct > 0 else np.nan
    # dq/dc = 1 / (2*sqrt(1-c))
    dq_dc_gd = 1 / (2 * math.sqrt(max(1e-9, 1 - car_gd))) if 1 - car_gd > 1e-9 else np.nan
    dq_dc_ct = 1 / (2 * math.sqrt(max(1e-9, 1 - car_ct))) if 1 - car_ct > 1e-9 else np.nan
    se_q_gd = dq_dc_gd * se_car_gd
    se_q_ct = dq_dc_ct * se_car_ct
    hwe_rows.append(dict(
        allele=a,
        carrier_freq_GD=car_gd, carrier_freq_ctrl=car_ct,
        allele_freq_GD_HWE=q_gd, allele_freq_ctrl_HWE=q_ct,
        carrier_OR=r["OR"],
        allele_OR_HWE=( (q_gd / (1-q_gd)) / (q_ct / (1-q_ct)) ) if (1-q_gd) > 0 and (1-q_ct) > 0 and q_ct > 0 else np.nan,
        ci_lo_q_gd=q_gd - 1.96 * se_q_gd if np.isfinite(se_q_gd) else np.nan,
        ci_hi_q_gd=q_gd + 1.96 * se_q_gd if np.isfinite(se_q_gd) else np.nan,
        ci_lo_q_ctrl=q_ct - 1.96 * se_q_ct if np.isfinite(se_q_ct) else np.nan,
        ci_hi_q_ctrl=q_ct + 1.96 * se_q_ct if np.isfinite(se_q_ct) else np.nan,
        note=r["note"],
    ))
hwe_df = pd.DataFrame(hwe_rows)
hwe_df.to_csv(TBL / "T05_carrier_vs_allele_HWE_reconciliation.tsv", sep="\t", index=False)
print(f"   -> {TBL/'T05_carrier_vs_allele_HWE_reconciliation.tsv'} ({len(hwe_df)} alleles)")

# ---------------------------------------------------------------------------
# 5. DPB1 sub-allele decomposition (DPB1*05:01 vs DPB1*05:XX vs other DPB1) in K2 baseline
# ---------------------------------------------------------------------------

print("[5] DPB1 high-resolution decomposition in K2 South-Korean NGS reference ...")
k2 = pd.read_csv(K2_TSV, sep="\t")
dpb1_alleles = pd.concat([k2["DPB1_a1_4digit"], k2["DPB1_a2_4digit"]]).dropna()
n_chrom = len(dpb1_alleles)
sub_count = dpb1_alleles.value_counts()
sub_freq = sub_count / n_chrom
top10 = sub_count.head(10).reset_index()
top10.columns = ["DPB1_subAllele_4digit", "n_chromosomes"]
top10["allele_freq_K2"] = top10["n_chromosomes"] / n_chrom
# Mark the *05:XX family vs others
top10["is_05_family"] = top10["DPB1_subAllele_4digit"].astype(str).str.startswith("DPB1*05:")
# Add carriers (any of two chromosomes)
n_subjects = len(k2)
c_05_01 = ((k2["DPB1_a1_4digit"] == "DPB1*05:01") | (k2["DPB1_a2_4digit"] == "DPB1*05:01")).sum()
c_05_any = (
    k2["DPB1_a1_4digit"].astype(str).str.startswith("DPB1*05:")
    | k2["DPB1_a2_4digit"].astype(str).str.startswith("DPB1*05:")
).sum()
dpb1_decomp = pd.DataFrame([
    dict(group="DPB1*05:01_only",
         n_chromosomes=int((dpb1_alleles == "DPB1*05:01").sum()),
         allele_freq=float((dpb1_alleles == "DPB1*05:01").sum() / n_chrom),
         carrier_freq=float(c_05_01 / n_subjects)),
    dict(group="DPB1*05:XX_family",
         n_chromosomes=int(dpb1_alleles.astype(str).str.startswith("DPB1*05:").sum()),
         allele_freq=float(dpb1_alleles.astype(str).str.startswith("DPB1*05:").sum() / n_chrom),
         carrier_freq=float(c_05_any / n_subjects)),
    dict(group="DPB1_other",
         n_chromosomes=int((~dpb1_alleles.astype(str).str.startswith("DPB1*05:")).sum()),
         allele_freq=float((~dpb1_alleles.astype(str).str.startswith("DPB1*05:")).sum() / n_chrom),
         carrier_freq=np.nan),
])
dpb1_decomp.to_csv(TBL / "T06_DPB1_subAllele_decomposition_K2.tsv", sep="\t", index=False)
top10.to_csv(TBL / "T06b_DPB1_top10_subAlleles_K2.tsv", sep="\t", index=False)
print(f"   -> {TBL/'T06_DPB1_subAllele_decomposition_K2.tsv'}")
print(f"   -> {TBL/'T06b_DPB1_top10_subAlleles_K2.tsv'}")
print(f"   K2 NGS Korean reference: n_subj={n_subjects}, n_chrom={n_chrom}; "
      f"DPB1*05:01 carriers={c_05_01}/{n_subjects} ({c_05_01/n_subjects:.3f}); "
      f"DPB1*05:XX family carriers={c_05_any}/{n_subjects} ({c_05_any/n_subjects:.3f})")

# ---------------------------------------------------------------------------
# 6. LD context for DPB1*05:01 ↔ B*46:01 in K2 NGS reference
# ---------------------------------------------------------------------------

print("[6] LD matrix for top GD-implicated alleles in K2 NGS reference ...")
LD_FOCUS = ["DPB1*05:01", "B*46:01", "C*01:02", "A*02:07",
            "DRB1*15:01", "DRB1*08:03", "DQB1*02:01", "DQB1*03:02"]

# Subject-level carrier indicator
def k2_carrier(allele_str):
    locus = allele_str.split("*")[0]
    a1 = f"{locus}_a1_4digit"; a2 = f"{locus}_a2_4digit"
    if a1 not in k2.columns:
        return None
    return ((k2[a1] == allele_str) | (k2[a2] == allele_str)).astype(int)

carrier_mat = {}
for a in LD_FOCUS:
    v = k2_carrier(a)
    if v is None:
        continue
    carrier_mat[a] = v
carrier_df = pd.DataFrame(carrier_mat)

# Pairwise: D' and r2 from carrier-level data approximation
def carrier_LD(a, b):
    """Approximate LD r2 between two carrier indicators (subject-level)."""
    pa = a.mean(); pb = b.mean()
    pab = (a & b).mean()
    D = pab - pa * pb
    denom = math.sqrt(pa * (1-pa) * pb * (1-pb))
    if denom == 0:
        return np.nan, np.nan
    r = D / denom
    return D, r * r

ld_rows = []
for i, a in enumerate(carrier_df.columns):
    for j, b in enumerate(carrier_df.columns):
        if j <= i:
            continue
        D, r2 = carrier_LD(carrier_df[a], carrier_df[b])
        # Fisher exact p
        ct = pd.crosstab(carrier_df[a], carrier_df[b])
        try:
            _, p_fisher = stats.fisher_exact(ct.values)
        except Exception:
            p_fisher = np.nan
        ld_rows.append(dict(
            allele1=a, allele2=b,
            n=int(carrier_df.shape[0]),
            carrier_freq_a=float(carrier_df[a].mean()),
            carrier_freq_b=float(carrier_df[b].mean()),
            joint_carrier_freq=float((carrier_df[a] & carrier_df[b]).mean()),
            D=float(D), r2=float(r2),
            fisher_p=float(p_fisher) if np.isfinite(p_fisher) else np.nan,
        ))
ld_df = pd.DataFrame(ld_rows).sort_values("r2", ascending=False)
ld_df.to_csv(TBL / "T07_LD_matrix_K2_carrier_level.tsv", sep="\t", index=False)
print(f"   -> {TBL/'T07_LD_matrix_K2_carrier_level.tsv'} ({len(ld_df)} pairs)")

# ---------------------------------------------------------------------------
# 7. AFND-anchored Korean baseline reconciliation
# ---------------------------------------------------------------------------

print("[7] AFND-anchored Korean baseline reconciliation ...")
afnd = json.load(open(AFND_J))
afnd_rows = []
for allele, by_country in afnd.items():
    for country, pops in by_country.items():
        for p in pops:
            try:
                n = int(str(p["sample_size"]).replace(",", ""))
            except Exception:
                n = np.nan
            try:
                f = float(p["allele_freq"])
            except Exception:
                f = np.nan
            afnd_rows.append(dict(
                allele=allele, country=country,
                population=p["population"],
                allele_freq=f, sample_size=n,
            ))
afnd_df = pd.DataFrame(afnd_rows)
afnd_df.to_csv(TBL / "T08_AFND_panasian_baselines.tsv", sep="\t", index=False)

# Sample-size-weighted Korean / Chinese / Japanese / Taiwan baseline per allele
afnd_summ = (
    afnd_df.dropna(subset=["allele_freq", "sample_size"])
    .groupby(["allele", "country"])
    .apply(lambda g: pd.Series({
        "n_pops": len(g),
        "tot_n_2N": g["sample_size"].sum(),
        "weighted_allele_freq": (g["allele_freq"] * g["sample_size"]).sum() / g["sample_size"].sum(),
        "min_freq": g["allele_freq"].min(),
        "max_freq": g["allele_freq"].max(),
    }))
    .reset_index()
)
afnd_summ.to_csv(TBL / "T08b_AFND_country_weighted_summary.tsv", sep="\t", index=False)
print(f"   -> {TBL/'T08_AFND_panasian_baselines.tsv'} ({len(afnd_df)} rows)")
print(f"   -> {TBL/'T08b_AFND_country_weighted_summary.tsv'} ({len(afnd_summ)} country-allele cells)")

# ---------------------------------------------------------------------------
# 8. Replication readiness scoreboard v2
# ---------------------------------------------------------------------------

print("[8] Replication-readiness scoreboard v2 ...")
score_rows = []
for allele in FOCUS_ALLELES:
    sub = canonical[canonical["allele_4d"] == allele]
    countries = sub["country"].astype(str).str.lower()
    has_korea = countries.str.contains("kor").any()
    has_china = countries.str.contains("china").any()
    has_taiwan = countries.str.contains("taiwan").any()
    has_japan = countries.str.contains("japan").any()
    n_total_case = sub["n_case"].astype(float).sum() if "n_case" in sub else np.nan
    n_total_ctrl = sub["n_control"].astype(float).sum() if "n_control" in sub else np.nan
    has_chu_anchor = (sub["source"] == "Chu 2018").any()
    has_afnd_baseline = allele in afnd
    grade_n = sum([has_korea, has_china, has_taiwan, has_japan])
    if grade_n >= 3 and has_chu_anchor:
        grade = "A_full_panasian_meta_ready"
        next_step = "Ready to publish; consider Korean adult NGS GD cohort to upgrade A→A+"
    elif grade_n == 2 and has_chu_anchor:
        grade = "B_partial_panasian_meta"
        next_step = "Add 1 more ancestry; Japanese full-table extraction (Onuma/Inoue/Ueda) preferred"
    elif has_chu_anchor:
        grade = "C_chu_anchor_only"
        next_step = "Korean adult GD NGS cohort + Taiwan or Japan replication"
    else:
        grade = "D_underpowered"
        next_step = "Need anchor cohort"
    score_rows.append(dict(
        allele=allele,
        n_sources=len(sub),
        countries_covered=";".join(sub["country"].astype(str).unique()),
        n_total_case=int(n_total_case) if np.isfinite(n_total_case) else np.nan,
        n_total_ctrl=int(n_total_ctrl) if np.isfinite(n_total_ctrl) else np.nan,
        has_chu_anchor=bool(has_chu_anchor),
        has_korean_replication=bool(has_korea),
        has_taiwan_replication=bool(has_taiwan),
        has_japanese_replication=bool(has_japan),
        has_AFND_baseline=bool(has_afnd_baseline),
        readiness_grade=grade,
        next_step=next_step,
    ))
score_df = pd.DataFrame(score_rows)
score_df.to_csv(TBL / "T09_replication_readiness_scoreboard_v2.tsv", sep="\t", index=False)
print(f"   -> {TBL/'T09_replication_readiness_scoreboard_v2.tsv'} ({len(score_df)} alleles)")

# ===========================================================================
# PLOTS
# ===========================================================================

print("[P] Forest plots ...")

# ----- Master forest: 11 alleles, one row per allele, RE pooled effect -----

def forest_master(meta_df: pd.DataFrame, src_df: pd.DataFrame, out_stem: Path):
    fig, ax = plt.subplots(figsize=(11, 9), dpi=DPI)
    md = meta_df.dropna(subset=["pooled_or"]).copy().reset_index(drop=True)
    md["log_or"] = np.log(md["pooled_or"])
    md["log_lo"] = np.log(md["ci_lo"])
    md["log_hi"] = np.log(md["ci_hi"])
    yticks = np.arange(len(md))
    risk = md["pooled_or"] > 1
    ax.errorbar(
        md["log_or"], yticks,
        xerr=[md["log_or"] - md["log_lo"], md["log_hi"] - md["log_or"]],
        fmt="s", color="black", capsize=3, lw=1.4, markersize=8,
    )
    # Color squares by direction
    for i, (lor, c) in enumerate(zip(md["log_or"], np.where(risk, "#cc3333", "#3366cc"))):
        ax.plot(lor, i, "s", color=c, markersize=10, markeredgecolor="black")
    # Prediction interval as thinner line
    for i, row in md.iterrows():
        if np.isfinite(row.get("pi_lo", np.nan)) and np.isfinite(row.get("pi_hi", np.nan)):
            ax.hlines(i, math.log(row["pi_lo"]), math.log(row["pi_hi"]), color="grey", lw=0.8, ls=":")
    ax.axvline(0, color="black", lw=0.8)
    ax.set_yticks(yticks)
    ax.set_yticklabels(md["allele"], fontsize=11)
    ax.invert_yaxis()
    # Log-OR axis with OR labels
    xt = [math.log(x) for x in [0.3, 0.5, 1, 2, 3, 5]]
    ax.set_xticks(xt)
    ax.set_xticklabels(["0.3", "0.5", "1", "2", "3", "5"])
    ax.set_xlabel("Pooled OR (random-effects, DerSimonian-Laird)", fontsize=12)
    # Annotate text on right margin
    txt_x = math.log(7.5)
    ax.set_xlim(math.log(0.2), math.log(8.5))
    ax.set_title(
        "Pan-Asian Graves' disease — HLA susceptibility random-effects meta-analysis\n"
        "(autoimmune-only; no thyroid cancer claim)",
        fontsize=12,
    )
    for i, row in md.iterrows():
        s = f"OR {row['pooled_or']:.2f}  [{row['ci_lo']:.2f}, {row['ci_hi']:.2f}]   k={row['k']}   I²={row['I2_pct']:.0f}%"
        ax.text(txt_x, i, s, fontsize=8.5, va="center", family="monospace")
    plt.tight_layout()
    fig.savefig(f"{out_stem}.png", dpi=DPI, bbox_inches="tight")
    fig.savefig(f"{out_stem}.pdf", bbox_inches="tight")
    fig.savefig(f"{out_stem}.svg", bbox_inches="tight")
    plt.close(fig)

forest_master(meta_df, canonical, FOR / "forest_panasian_GD_v2")
print(f"   -> {FOR/'forest_panasian_GD_v2.{png,pdf,svg}'}")

# ----- Per-allele forest pages (one PNG per allele) ------

def forest_per_allele(allele: str, src_df: pd.DataFrame, meta_row, out_path: Path):
    sub = src_df[src_df["allele_4d"] == allele].sort_values(["country", "source"]).reset_index(drop=True)
    if len(sub) == 0:
        return
    fig, ax = plt.subplots(figsize=(10, max(3.0, 0.55 * len(sub) + 2.5)), dpi=DPI)
    yticks = np.arange(len(sub))
    log_or = sub["log_or"].astype(float).values
    se = sub["se"].astype(float).values
    log_lo = log_or - 1.959964 * se
    log_hi = log_or + 1.959964 * se
    # Per-study squares sized by inverse-variance weight
    weights = 1.0 / (se ** 2)
    weights = weights / weights.sum()
    sizes = 80 + 700 * weights
    colors = ["#cc3333" if (lo > 0) else ("#3366cc" if hi < 0 else "#888888")
              for lo, hi in zip(log_lo, log_hi)]
    ax.errorbar(
        log_or, yticks,
        xerr=[log_or - log_lo, log_hi - log_or],
        fmt="none", color="black", capsize=3, lw=1.0,
    )
    for i, (lor, c, sz) in enumerate(zip(log_or, colors, sizes)):
        ax.scatter([lor], [i], s=sz, color=c, edgecolor="black", zorder=5)
    # Pooled diamond (RE)
    if meta_row is not None and np.isfinite(meta_row["pooled_or"]):
        diamond_y = len(sub) + 0.5
        cx = math.log(meta_row["pooled_or"])
        clo = math.log(meta_row["ci_lo"])
        chi = math.log(meta_row["ci_hi"])
        diamond = plt.Polygon(
            [(clo, diamond_y), (cx, diamond_y - 0.4), (chi, diamond_y), (cx, diamond_y + 0.4)],
            facecolor="black", edgecolor="black", zorder=10,
        )
        ax.add_patch(diamond)
        ax.text(
            chi + 0.05, diamond_y,
            f"RE pooled OR={meta_row['pooled_or']:.2f}  [{meta_row['ci_lo']:.2f},{meta_row['ci_hi']:.2f}]  "
            f"I²={meta_row['I2_pct']:.0f}%  τ²={meta_row['tau2']:.3f}  k={meta_row['k']}",
            fontsize=9, va="center", family="monospace",
        )
    ax.axvline(0, color="black", lw=0.8)
    ax.set_yticks(list(yticks) + [len(sub) + 0.5])
    ylabels = []
    for _, r in sub.iterrows():
        n_str = f"n_case={int(r['n_case'])}/n_ctrl={int(r['n_control'])}" \
            if np.isfinite(r["n_case"]) and np.isfinite(r["n_control"]) else "n=NA"
        ylabels.append(f"{r['source']} ({r['country']}, {n_str})")
    ylabels.append("Pan-Asian RE pooled")
    ax.set_yticklabels(ylabels, fontsize=9)
    ax.invert_yaxis()
    xt = [math.log(x) for x in [0.3, 0.5, 1, 2, 3, 5]]
    ax.set_xticks(xt)
    ax.set_xticklabels(["0.3", "0.5", "1", "2", "3", "5"])
    ax.set_xlabel("OR (95% CI), Graves' disease vs ancestry-matched control")
    ax.set_xlim(math.log(0.2), math.log(8.5))
    ax.set_title(f"Pan-Asian GD HLA — {allele}\n(autoimmune Graves' disease; no cancer claim)",
                 fontsize=12)
    plt.tight_layout()
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)

per_allele_count = 0
for allele in FOCUS_ALLELES:
    sub = canonical[canonical["allele_4d"] == allele]
    if len(sub) == 0:
        continue
    mrow = meta_df[meta_df["allele"] == allele]
    mrow = mrow.iloc[0] if len(mrow) else None
    safe = allele.replace("*", "_").replace(":", "_")
    forest_per_allele(allele, canonical, mrow, PER / f"forest_{safe}.png")
    per_allele_count += 1
print(f"   -> {PER}/forest_<allele>.png  ({per_allele_count} alleles)")

# ----- Galbraith plot (all study-effects across all alleles) -----

print("[P2] Galbraith + funnel plots ...")
galb_df = canonical[np.isfinite(canonical["log_or"]) & np.isfinite(canonical["se"]) & (canonical["se"] > 0)].copy()
galb_df["precision"] = 1.0 / galb_df["se"]
galb_df["z"] = galb_df["log_or"] / galb_df["se"]
fig, ax = plt.subplots(figsize=PLT_FIGSIZE, dpi=DPI)
markers = {"DPB1*05:01": "o", "B*46:01": "s", "C*01:02": "^", "A*02:07": "D",
           "DQB1*02:01": "v", "DRB1*07:01": "P"}
default_marker = "x"
for allele, sub in galb_df.groupby("allele_4d"):
    m = markers.get(allele, default_marker)
    ax.scatter(sub["precision"], sub["z"], label=allele, marker=m, s=60, alpha=0.85)
ax.axhline(1.96, color="grey", ls=":", lw=0.8)
ax.axhline(-1.96, color="grey", ls=":", lw=0.8)
ax.axhline(0, color="black", lw=0.8)
ax.set_xlabel("Precision = 1/SE(log OR)")
ax.set_ylabel("Standardized effect z = log(OR) / SE")
ax.set_title("Galbraith plot — Pan-Asian GD HLA studies")
ax.legend(fontsize=8, ncol=2)
plt.tight_layout()
fig.savefig(PLT / "galbraith_panasian_GD.png", dpi=DPI, bbox_inches="tight")
plt.close(fig)

# Funnel for the strongest DPB1*05:01 axis
dpb = canonical[canonical["allele_4d"] == "DPB1*05:01"].copy()
fig, ax = plt.subplots(figsize=PLT_FIGSIZE, dpi=DPI)
ax.scatter(dpb["log_or"], dpb["se"], color="black", s=80)
for _, r in dpb.iterrows():
    ax.annotate(r["source"], (r["log_or"], r["se"]), fontsize=8, xytext=(5, 0),
                textcoords="offset points")
m_re = meta_df[meta_df["allele"] == "DPB1*05:01"].iloc[0]
ax.axvline(np.log(m_re["pooled_or"]), color="red", ls="--", lw=1, label=f"RE pooled OR={m_re['pooled_or']:.2f}")
ax.invert_yaxis()
ax.set_xlabel("log(OR)")
ax.set_ylabel("SE(log OR)")
ax.set_title("Funnel plot — DPB1*05:01 (Pan-Asian Graves' disease)")
ax.legend()
plt.tight_layout()
fig.savefig(PLT / "funnel_DPB1_0501.png", dpi=DPI, bbox_inches="tight")
plt.close(fig)
print(f"   -> {PLT/'galbraith_panasian_GD.png'}, {PLT/'funnel_DPB1_0501.png'}")

# ----- LD heatmap for K2 carrier-level matrix -----

print("[P3] LD heatmap (K2 NGS Korean reference) ...")
LD_focus_present = [a for a in LD_FOCUS if a in carrier_df.columns]
if LD_focus_present:
    n_a = len(LD_focus_present)
    r2_mat = np.full((n_a, n_a), np.nan)
    p_mat  = np.full((n_a, n_a), np.nan)
    for i, a in enumerate(LD_focus_present):
        for j, b in enumerate(LD_focus_present):
            if i == j:
                r2_mat[i, j] = 1.0
                continue
            row = ld_df[((ld_df["allele1"] == a) & (ld_df["allele2"] == b)) |
                        ((ld_df["allele1"] == b) & (ld_df["allele2"] == a))]
            if len(row):
                r2_mat[i, j] = row["r2"].iloc[0]
                p_mat[i, j] = row["fisher_p"].iloc[0]
    fig, ax = plt.subplots(figsize=PLT_FIGSIZE, dpi=DPI)
    im = ax.imshow(r2_mat, cmap="Reds", vmin=0, vmax=max(0.05, np.nanmax(r2_mat)))
    ax.set_xticks(range(n_a)); ax.set_xticklabels(LD_focus_present, rotation=45, ha="right")
    ax.set_yticks(range(n_a)); ax.set_yticklabels(LD_focus_present)
    plt.colorbar(im, ax=ax, label="carrier-level r²")
    for i in range(n_a):
        for j in range(n_a):
            if i == j or not np.isfinite(r2_mat[i, j]):
                continue
            txt = f"{r2_mat[i,j]:.2f}"
            if np.isfinite(p_mat[i, j]) and p_mat[i, j] < 0.05:
                txt += "*"
            ax.text(j, i, txt, ha="center", va="center",
                    color="black" if r2_mat[i, j] < 0.4 else "white", fontsize=8)
    ax.set_title("Pairwise carrier-level r² in K2 South-Korean NGS reference (n=%d)" % n_subjects)
    plt.tight_layout()
    fig.savefig(PLT / "LD_heatmap_K2.png", dpi=DPI, bbox_inches="tight")
    plt.close(fig)
print(f"   -> {PLT/'LD_heatmap_K2.png'}")

# ----- Sub-population stratified forest (multi-panel) -----

strat_plot = strat_df[strat_df["k"] > 0].copy()
if len(strat_plot):
    n_alleles = strat_plot["allele"].nunique()
    fig, ax = plt.subplots(figsize=(11, max(6, 0.5 * len(strat_plot) + 2)), dpi=DPI)
    yticks = np.arange(len(strat_plot))
    strat_plot["log_or"] = np.log(strat_plot["pooled_or"])
    strat_plot["log_lo"] = np.log(strat_plot["ci_lo"])
    strat_plot["log_hi"] = np.log(strat_plot["ci_hi"])
    ancestry_color = {"Korean": "#7e57c2", "Han_Chinese": "#ef5350", "Taiwan": "#ffb300",
                      "Japanese": "#26a69a", "POOLED": "black"}
    colors = strat_plot["ancestry"].map(ancestry_color)
    ax.errorbar(
        strat_plot["log_or"], yticks,
        xerr=[strat_plot["log_or"] - strat_plot["log_lo"],
              strat_plot["log_hi"] - strat_plot["log_or"]],
        fmt="none", color="black", capsize=3, lw=1.0,
    )
    for i, (lor, c) in enumerate(zip(strat_plot["log_or"], colors)):
        ax.scatter([lor], [i], color=c, s=80, edgecolor="black", zorder=5)
    ax.axvline(0, color="black", lw=0.8)
    ax.set_yticks(yticks)
    ylabels = [f"{r['allele']:>12s}  |  {r['ancestry']:<11s}  k={r['k']}" for _, r in strat_plot.iterrows()]
    ax.set_yticklabels(ylabels, fontsize=8, family="monospace")
    ax.invert_yaxis()
    xt = [math.log(x) for x in [0.3, 0.5, 1, 2, 3, 5]]
    ax.set_xticks(xt); ax.set_xticklabels(["0.3", "0.5", "1", "2", "3", "5"])
    ax.set_xlabel("Pooled OR (random-effects)")
    ax.set_title("Pan-Asian GD HLA — sub-population stratified meta")
    legend_handles = [mpatches.Patch(color=c, label=k) for k, c in ancestry_color.items()]
    ax.legend(handles=legend_handles, loc="lower right", fontsize=8)
    plt.tight_layout()
    fig.savefig(PLT / "stratified_panasian_GD.png", dpi=DPI, bbox_inches="tight")
    plt.close(fig)
print(f"   -> {PLT/'stratified_panasian_GD.png'}")

# ===========================================================================
# Summary JSON for downstream MD
# ===========================================================================

# Headline numbers
top_meta = meta_df[meta_df["k"] >= 1].sort_values("p_random").head(8)
summary = dict(
    n_alleles_meta=int(len(meta_df)),
    n_alleles_with_k_ge_2=int((meta_df["k"] >= 2).sum()),
    n_panasian_source_rows=int(len(canonical)),
    DPB1_05_01_pooled_OR_random=float(meta_df.loc[meta_df.allele=="DPB1*05:01","pooled_or"].iloc[0]),
    DPB1_05_01_ci=[
        float(meta_df.loc[meta_df.allele=="DPB1*05:01","ci_lo"].iloc[0]),
        float(meta_df.loc[meta_df.allele=="DPB1*05:01","ci_hi"].iloc[0]),
    ],
    DPB1_05_01_I2=float(meta_df.loc[meta_df.allele=="DPB1*05:01","I2_pct"].iloc[0]),
    DPB1_05_01_pi=[
        float(meta_df.loc[meta_df.allele=="DPB1*05:01","pi_lo"].iloc[0])
            if np.isfinite(meta_df.loc[meta_df.allele=="DPB1*05:01","pi_lo"].iloc[0]) else None,
        float(meta_df.loc[meta_df.allele=="DPB1*05:01","pi_hi"].iloc[0])
            if np.isfinite(meta_df.loc[meta_df.allele=="DPB1*05:01","pi_hi"].iloc[0]) else None,
    ],
    K2_DPB1_0501_carrier_freq=float(c_05_01 / n_subjects),
    K2_n_subjects=int(n_subjects),
    K2_n_chrom=int(n_chrom),
    top_alleles=top_meta.to_dict(orient="records"),
)
with open(OUT / "track1_summary.json", "w") as f:
    json.dump(summary, f, indent=2, default=lambda x: None if isinstance(x, float) and not np.isfinite(x) else x)
print(f"[summary] -> {OUT/'track1_summary.json'}")

print("DONE.")
