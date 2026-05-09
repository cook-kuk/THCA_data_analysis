#!/usr/bin/env python3
"""R10 — Paper 3 BRAF Nature sprint
Anti-PD-1 vs anti-PD-L1 regimen stratification of R5 predictors.

Input:  R5 per-sample scores (n=448, 5 cohorts, response/HT13/TIS/DM1_inflam/HLA1)
Output: r10_regimen_stratified_meta.tsv, r10_2x2_partition.tsv,
        r10_pdl1_ihc_concordance.tsv (IMvigor210 only),
        R10_REPORT.md
"""
from __future__ import annotations

import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")

ROOT = Path("/home/seungho/personal/THCA_data_analysis/project/results/p2_braf_nature_sprint_2026_05_09")
R5 = ROOT / "r5_ici_cohorts" / "r5_per_sample_scores.tsv"
OUT = ROOT / "r10_ici_regimen"
OUT.mkdir(parents=True, exist_ok=True)


# ----------------------------------------------------------------------
# 1. Load R5 + annotate regimen
# ----------------------------------------------------------------------
df = pd.read_csv(R5, sep="\t")

def regimen_of(row):
    """Map (cohort, therapy) to PD1 / PDL1 / OTHER."""
    coh = row["cohort_r5"]
    th = str(row["therapy"]).lower()
    if coh == "IMvigor210":
        return "PDL1"
    if coh == "Hugo_GSE78220":
        return "PD1"
    if coh == "riaz_GSE91061":
        return "PD1"
    if coh == "MGH_GSE115821":
        if "anti-pd-1" == th:
            return "PD1"
        # anti-ctla-4 alone or +ctla4 combos -> exclude
        return "OTHER"
    if coh == "GSE176307":
        # mixed-regimen UC cohort
        if th in ("atezolizumab", "avelumab", "durvalumab"):
            return "PDL1"
        if th in ("nivolumab", "pembrolizumab"):
            return "PD1"
        return "OTHER"
    return "OTHER"


df["regimen"] = df.apply(regimen_of, axis=1)

# Build a regimen-aware "cohort key" so PD-L1 vs PD-1 sub-groups inside
# GSE176307 are treated as separate strata for meta-analysis.
df["regimen_cohort"] = np.where(
    df["cohort_r5"] == "GSE176307",
    "GSE176307_" + df["regimen"],
    df["cohort_r5"],
)

# Restrict to samples with response_binary defined and regimen ∈ {PD1, PDL1}
df = df[df["regimen"].isin(["PD1", "PDL1"])].copy()
df["resp"] = pd.to_numeric(df["response_binary_CRPR_vs_SD_PD"], errors="coerce")
df = df.dropna(subset=["resp"]).copy()
df["resp"] = df["resp"].astype(int)


# ----------------------------------------------------------------------
# 2. Helpers
# ----------------------------------------------------------------------
def fisher_high_low(score, resp, q=0.5):
    """High vs low (median split) → 2x2 OR + Fisher p."""
    score = np.asarray(score, dtype=float)
    resp = np.asarray(resp, dtype=int)
    mask = np.isfinite(score) & np.isfinite(resp)
    score, resp = score[mask], resp[mask]
    if len(score) < 6:
        return np.nan, np.nan, np.nan, np.nan, len(score), int(resp.sum())
    thr = np.quantile(score, q)
    high = score > thr
    a = int(((high) & (resp == 1)).sum())
    b = int(((high) & (resp == 0)).sum())
    c = int(((~high) & (resp == 1)).sum())
    d = int(((~high) & (resp == 0)).sum())
    table = np.array([[a, b], [c, d]])
    if (table.min(axis=0) == 0).any():
        # Haldane correction for OR but Fisher exact for p
        a_, b_, c_, d_ = a + 0.5, b + 0.5, c + 0.5, d + 0.5
        or_ = (a_ * d_) / (b_ * c_)
        _, p = stats.fisher_exact(table)
        return or_, p, len(score), int(resp.sum()), a, b
    or_ = (a * d) / (b * c)
    _, p = stats.fisher_exact(table)
    return or_, p, len(score), int(resp.sum()), a, b


def per_cohort_panel_or(sub: pd.DataFrame, panels=("HT13", "TIS", "DM1_inflam")):
    """Per-cohort OR for each panel."""
    rows = []
    for coh, g in sub.groupby("regimen_cohort"):
        for panel in panels:
            score = g[panel].to_numpy(float)
            resp = g["resp"].to_numpy(int)
            or_, p, n, n_resp, a, b = fisher_high_low(score, resp)
            # build 2x2 for MH later
            mask = np.isfinite(score)
            score_m = score[mask]
            resp_m = resp[mask]
            if len(score_m) < 6:
                continue
            thr = np.quantile(score_m, 0.5)
            high = score_m > thr
            aa = int(((high) & (resp_m == 1)).sum())
            bb = int(((high) & (resp_m == 0)).sum())
            cc = int(((~high) & (resp_m == 1)).sum())
            dd = int(((~high) & (resp_m == 0)).sum())
            rows.append(dict(regimen_cohort=coh, panel=panel, OR=or_, p=p,
                             n=n, n_resp=n_resp,
                             a=aa, b=bb, c=cc, d=dd))
    return pd.DataFrame(rows)


def per_cohort_adaptive_or(sub: pd.DataFrame):
    """DM1-adaptive = HT13 high AND HLA1 high (within cohort medians).
       Compare adaptive vs non-adaptive responder rate per cohort.
    """
    rows = []
    for coh, g in sub.groupby("regimen_cohort"):
        g = g.dropna(subset=["HT13", "HLA1", "resp"])
        if len(g) < 6:
            continue
        ht_thr = g["HT13"].median()
        hla_thr = g["HLA1"].median()
        ht_high = g["HT13"] > ht_thr
        hla_high = g["HLA1"] > hla_thr
        adaptive = ht_high & hla_high
        a = int(((adaptive) & (g["resp"] == 1)).sum())
        b = int(((adaptive) & (g["resp"] == 0)).sum())
        c = int(((~adaptive) & (g["resp"] == 1)).sum())
        d = int(((~adaptive) & (g["resp"] == 0)).sum())
        table = np.array([[a, b], [c, d]])
        if (table.min(axis=0) == 0).any():
            a_, b_, c_, d_ = a + 0.5, b + 0.5, c + 0.5, d + 0.5
            or_ = (a_ * d_) / (b_ * c_)
        else:
            or_ = (a * d) / (b * c) if b * c > 0 else np.nan
        _, p = stats.fisher_exact(table)
        rows.append(dict(regimen_cohort=coh, panel="DM1_adaptive_HT13xHLA1",
                          OR=or_, p=p, n=int(len(g)), n_resp=int(g["resp"].sum()),
                          a=a, b=b, c=c, d=d,
                          adapt_resp_rate=a/(a+b) if (a+b)>0 else np.nan,
                          other_resp_rate=c/(c+d) if (c+d)>0 else np.nan))
    return pd.DataFrame(rows)


def mh_pool(tab: pd.DataFrame):
    """Mantel-Haenszel pooled OR with Robins-Breslow-Greenland SE.
       tab columns: a, b, c, d.
    """
    a = tab["a"].to_numpy(float); b = tab["b"].to_numpy(float)
    c = tab["c"].to_numpy(float); d = tab["d"].to_numpy(float)
    n = a + b + c + d
    num = (a * d / n).sum()
    den = (b * c / n).sum()
    if den <= 0 or num <= 0:
        return dict(mh_OR=np.nan, lo=np.nan, hi=np.nan, p=np.nan,
                    I2=np.nan, k=int(len(tab)))
    mh_or = num / den
    P = (a + d) / n; Q = (b + c) / n
    R = (a * d) / n; S = (b * c) / n
    sumR = R.sum(); sumS = S.sum()
    var_ln = ((P * R).sum() / (2 * sumR ** 2)
              + ((P * S + Q * R).sum()) / (2 * sumR * sumS)
              + (Q * S).sum() / (2 * sumS ** 2))
    se_ln = float(np.sqrt(var_ln))
    ln_or = float(np.log(mh_or))
    lo = float(np.exp(ln_or - 1.96 * se_ln))
    hi = float(np.exp(ln_or + 1.96 * se_ln))
    z = ln_or / se_ln
    p = float(2 * (1 - stats.norm.cdf(abs(z))))
    # heterogeneity I2 with Haldane correction so zero cells don't kill it
    a_h = a + 0.5; b_h = b + 0.5; c_h = c + 0.5; d_h = d + 0.5
    or_i = (a_h * d_h) / (b_h * c_h)
    ln_or_i = np.log(or_i)
    w_i = 1.0 / (1.0/a_h + 1.0/b_h + 1.0/c_h + 1.0/d_h)
    ln_pool = (w_i * ln_or_i).sum() / w_i.sum()
    Q_stat = (w_i * (ln_or_i - ln_pool) ** 2).sum()
    df_q = len(tab) - 1
    I2 = max(0.0, (Q_stat - df_q) / Q_stat) * 100 if Q_stat > 0 else 0.0
    return dict(mh_OR=float(mh_or), lo=lo, hi=hi, p=p, I2=float(I2),
                k=int(len(tab)))


# ----------------------------------------------------------------------
# 3. Per-regimen meta — HT13 / TIS / DM1_inflam / DM1_adaptive
# ----------------------------------------------------------------------
out_rows = []
per_cohort_records = []

for regimen, sub in df.groupby("regimen"):
    panel_tab = per_cohort_panel_or(sub)
    adapt_tab = per_cohort_adaptive_or(sub)
    tab_all = pd.concat([panel_tab, adapt_tab], ignore_index=True)
    tab_all["regimen"] = regimen
    per_cohort_records.append(tab_all)
    for panel, g in tab_all.groupby("panel"):
        if len(g) == 0:
            continue
        pool = mh_pool(g)
        rec = dict(regimen=regimen, panel=panel,
                    n_total=int(g["n"].sum()),
                    n_resp_total=int((g["a"] + g["c"]).sum()),
                    sign_concordance=int((g["OR"] > 1).sum()))
        rec.update(pool)
        out_rows.append(rec)

# Also overall (regimen-pooled) and per-cohort raw
per_cohort_records = pd.concat(per_cohort_records, ignore_index=True)

# Append total (PD1 + PDL1) for reference
all_panel_tab = per_cohort_panel_or(df)
all_adapt_tab = per_cohort_adaptive_or(df)
all_tab = pd.concat([all_panel_tab, all_adapt_tab], ignore_index=True)
for panel, g in all_tab.groupby("panel"):
    pool = mh_pool(g)
    rec = dict(regimen="POOLED_ALL", panel=panel,
                n_total=int(g["n"].sum()),
                n_resp_total=int((g["a"] + g["c"]).sum()),
                sign_concordance=int((g["OR"] > 1).sum()))
    rec.update(pool)
    out_rows.append(rec)

meta = pd.DataFrame(out_rows)
meta = meta[["regimen", "panel", "k", "n_total", "n_resp_total",
             "mh_OR", "lo", "hi", "p", "I2", "sign_concordance"]]
meta.to_csv(OUT / "r10_regimen_stratified_meta.tsv", sep="\t", index=False)
per_cohort_records.to_csv(OUT / "r10_per_cohort_within_regimen.tsv",
                           sep="\t", index=False)


# ----------------------------------------------------------------------
# 4. Cox OS per regimen
# ----------------------------------------------------------------------
try:
    from lifelines import CoxPHFitter
    have_lifelines = True
except Exception:
    have_lifelines = False

cox_rows = []
if have_lifelines:
    for regimen, sub in df.groupby("regimen"):
        for panel in ("HT13", "TIS", "DM1_inflam"):
            g = sub.dropna(subset=["os_time", "os_event", panel]).copy()
            g["os_time"] = pd.to_numeric(g["os_time"], errors="coerce")
            g["os_event"] = pd.to_numeric(g["os_event"], errors="coerce")
            g = g.dropna(subset=["os_time", "os_event"])
            g = g[g["os_time"] > 0]
            if len(g) < 20 or g["os_event"].sum() < 5:
                continue
            # standardise score within regimen pool
            g["z"] = (g[panel] - g[panel].mean()) / g[panel].std(ddof=0)
            try:
                cph = CoxPHFitter()
                # stratify by cohort to keep regimen-level pooling clean
                cph.fit(g[["os_time", "os_event", "z", "regimen_cohort"]],
                        duration_col="os_time", event_col="os_event",
                        strata=["regimen_cohort"])
                hr = float(np.exp(cph.params_["z"]))
                ci = cph.confidence_intervals_.loc["z"].to_numpy()
                lo = float(np.exp(ci[0])); hi = float(np.exp(ci[1]))
                p = float(cph.summary.loc["z", "p"])
                cox_rows.append(dict(regimen=regimen, panel=panel,
                                      n=len(g), n_event=int(g["os_event"].sum()),
                                      HR=hr, lo=lo, hi=hi, p=p,
                                      strata="regimen_cohort"))
            except Exception as e:
                cox_rows.append(dict(regimen=regimen, panel=panel,
                                      n=len(g), n_event=int(g["os_event"].sum()),
                                      HR=np.nan, lo=np.nan, hi=np.nan, p=np.nan,
                                      strata=str(e)[:60]))

cox = pd.DataFrame(cox_rows)
cox.to_csv(OUT / "r10_regimen_cox_os.tsv", sep="\t", index=False)


# ----------------------------------------------------------------------
# 5. 2x2 HLA × HT13 partition — DM1-adaptive vs DM2-classical-escape vs rare
# ----------------------------------------------------------------------
part_rows = []
for regimen, sub in df.groupby("regimen"):
    sub = sub.dropna(subset=["HT13", "HLA1", "resp"]).copy()
    if len(sub) < 6:
        continue
    # use within-regimen medians so partition is balanced
    ht_med = sub["HT13"].median()
    hla_med = sub["HLA1"].median()
    sub["HT_high"] = sub["HT13"] > ht_med
    sub["HLA_high"] = sub["HLA1"] > hla_med
    classes = {
        "DM1_adaptive (HT_high & HLA_high)":   ( True,  True),
        "DM1_classical_escape (HT_high & HLA_low)": ( True, False),
        "DM2_adaptive (HT_low & HLA_high)":    (False,  True),
        "DM2_classical_escape (HT_low & HLA_low)": (False, False),
    }
    for label, (ht, hla) in classes.items():
        m = (sub["HT_high"] == ht) & (sub["HLA_high"] == hla)
        n = int(m.sum())
        nr = int(sub.loc[m, "resp"].sum())
        part_rows.append(dict(regimen=regimen, class_=label,
                               n=n, n_resp=nr,
                               resp_rate=nr / n if n > 0 else np.nan))
    # OR for DM1-adaptive vs DM2-classical-escape
    a = int(((sub["HT_high"]) & (sub["HLA_high"]) & (sub["resp"] == 1)).sum())
    b = int(((sub["HT_high"]) & (sub["HLA_high"]) & (sub["resp"] == 0)).sum())
    c = int(((~sub["HT_high"]) & (~sub["HLA_high"]) & (sub["resp"] == 1)).sum())
    d = int(((~sub["HT_high"]) & (~sub["HLA_high"]) & (sub["resp"] == 0)).sum())
    table = np.array([[a, b], [c, d]])
    a_, b_, c_, d_ = a + 0.5, b + 0.5, c + 0.5, d + 0.5
    or_diag = (a_ * d_) / (b_ * c_)
    _, p_diag = stats.fisher_exact(table)
    part_rows.append(dict(regimen=regimen,
                           class_="DIAGONAL_OR (DM1_adaptive vs DM2_classical)",
                           n=int(a + b + c + d),
                           n_resp=int(a + c), resp_rate=or_diag, p=p_diag))

part = pd.DataFrame(part_rows)
part.to_csv(OUT / "r10_2x2_partition.tsv", sep="\t", index=False)


# ----------------------------------------------------------------------
# 6. PD-L1 IHC concordance — IMvigor210 only (Ventana SP142, IC0/IC1/IC2+)
# ----------------------------------------------------------------------
def parse_ic(notes):
    if not isinstance(notes, str):
        return np.nan
    for tok in notes.split(";"):
        tok = tok.strip()
        if tok.startswith("IC:"):
            v = tok.split(":", 1)[1].strip()
            return v
    return np.nan


pdl1_rows = []
imv = df[df["cohort_r5"] == "IMvigor210"].copy()
imv["IC"] = imv["notes"].apply(parse_ic)
imv = imv.dropna(subset=["IC"])
imv["IC_high"] = imv["IC"].isin(["IC2+"]).astype(int)
imv["HT_high"] = (imv["HT13"] > imv["HT13"].median()).astype(int)
imv["TIS_high"] = (imv["TIS"] > imv["TIS"].median()).astype(int)

# concordance crosstab
ct = pd.crosstab(imv["IC"], imv["HT_high"], margins=True)
ct.to_csv(OUT / "r10_pdl1_ihc_HT_crosstab.tsv", sep="\t")

# IC alone OR
def or_resp(score, resp):
    a = int(((score == 1) & (resp == 1)).sum())
    b = int(((score == 1) & (resp == 0)).sum())
    c = int(((score == 0) & (resp == 1)).sum())
    d = int(((score == 0) & (resp == 0)).sum())
    a_, b_, c_, d_ = a + 0.5, b + 0.5, c + 0.5, d + 0.5
    or_ = (a_ * d_) / (b_ * c_)
    _, p = stats.fisher_exact(np.array([[a, b], [c, d]]))
    return or_, p, a, b, c, d


for label, x in [("IC2+ vs IC0/IC1", imv["IC_high"]),
                  ("HT13_high (median split)", imv["HT_high"]),
                  ("TIS_high (median split)", imv["TIS_high"])]:
    or_, p, a, b, c, d = or_resp(x.to_numpy(int), imv["resp"].to_numpy(int))
    pdl1_rows.append(dict(predictor=label, OR=or_, p=p,
                           a=a, b=b, c=c, d=d, n=len(imv)))

# Multivariable: HT13_high adjusted for IC
import statsmodels.api as sm
imv_clean = imv.dropna(subset=["resp", "HT_high", "IC_high", "TIS_high"])
X = imv_clean[["HT_high", "IC_high"]].astype(float)
X = sm.add_constant(X)
y = imv_clean["resp"].astype(int)
try:
    logit = sm.Logit(y, X).fit(disp=False)
    for name in ("HT_high", "IC_high"):
        beta = float(logit.params[name])
        se = float(logit.bse[name])
        or_ = float(np.exp(beta))
        lo = float(np.exp(beta - 1.96 * se))
        hi = float(np.exp(beta + 1.96 * se))
        p = float(logit.pvalues[name])
        pdl1_rows.append(dict(predictor=f"{name} | adj IC2+",
                               OR=or_, p=p, lo=lo, hi=hi,
                               n=len(imv_clean)))
except Exception as e:
    pdl1_rows.append(dict(predictor=f"multivariable_failed:{e}", OR=np.nan, p=np.nan))

pdl1 = pd.DataFrame(pdl1_rows)
pdl1.to_csv(OUT / "r10_pdl1_ihc_concordance.tsv", sep="\t", index=False)


# ----------------------------------------------------------------------
# 7. Summary JSON
# ----------------------------------------------------------------------
headline = {}
for (regimen, panel), g in meta.groupby(["regimen", "panel"]):
    headline[f"{regimen}__{panel}"] = {
        "k": int(g["k"].iloc[0]),
        "n_total": int(g["n_total"].iloc[0]),
        "OR": float(g["mh_OR"].iloc[0]) if pd.notna(g["mh_OR"].iloc[0]) else None,
        "lo": float(g["lo"].iloc[0]) if pd.notna(g["lo"].iloc[0]) else None,
        "hi": float(g["hi"].iloc[0]) if pd.notna(g["hi"].iloc[0]) else None,
        "p": float(g["p"].iloc[0]) if pd.notna(g["p"].iloc[0]) else None,
        "I2": float(g["I2"].iloc[0]) if pd.notna(g["I2"].iloc[0]) else None,
        "sign_concordance": int(g["sign_concordance"].iloc[0]),
    }

with open(OUT / "r10_headline.json", "w") as f:
    json.dump(headline, f, indent=2)

print("=== R10 regimen-stratified meta (DM1_adaptive_HT13xHLA1) ===")
print(meta[meta["panel"] == "DM1_adaptive_HT13xHLA1"].to_string(index=False))
print()
print("=== R10 regimen-stratified meta (HT13 / TIS / DM1_inflam) ===")
print(meta[meta["panel"].isin(["HT13", "TIS", "DM1_inflam"])].to_string(index=False))
print()
print("=== Cox OS per regimen ===")
print(cox.to_string(index=False))
print()
print("=== 2x2 HLA × HT partition ===")
print(part.to_string(index=False))
print()
print("=== PD-L1 IHC concordance (IMvigor210) ===")
print(pdl1.to_string(index=False))
print()
print("Outputs in", OUT)
