#!/usr/bin/env python3
"""R5 — HT-13 panel × ICI response across 5 immunotherapy cohorts.

Tests whether the HT-13 Hashimoto-thyroiditis-like B/MHC-II/IFNG panel
predicts ICI response (CR/PR vs SD/PD) — extending the H27 candidacy
claim to actual treated patients.

Cohorts:
  - IMvigor210 (UC anti-PD-L1, Mariathasan 2018)
  - Riaz GSE91061 pre (mel anti-PD-1, Riaz 2017)
  - MGH GSE115821 (mel anti-PD-1/CTLA-4, Auslander 2018)
  - GSE176307 (UC anti-PD-L1, BACI real-world)
  - Hugo GSE78220 (mel anti-PD-1, Hugo 2016)

For each cohort:
  - HT-13 = mean within-cohort z-score of available HT-13 genes
  - TIS = mean z-score of 18-gene Tumour Inflammation Signature (Ayers 2017)
  - DM1_inflam = mean(HT-13, TIS)
  - Logistic regression  score → response_binary (CR/PR vs SD/PD)
       OR + 95% CI + p (Wald), AUC, n
  - Cox PFS / OS (continuous score) where time available
  - DM1-like-adaptive label = HT-13 high (>median) AND HLA-I high (>median)
       Fisher exact for that vs response_binary

Pooled meta:
  - Random-effects (DerSimonian-Laird) logOR across cohorts
  - I^2 heterogeneity
  - Forest plot data table

Outputs to project/results/p2_braf_nature_sprint_2026_05_09/r5_ici_cohorts/.
"""
from __future__ import annotations
import json
import warnings
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
import gzip
import re

warnings.filterwarnings("ignore")

OUT = Path("/home/seungho/personal/THCA_data_analysis/project/results/p2_braf_nature_sprint_2026_05_09/r5_ici_cohorts")
OUT.mkdir(parents=True, exist_ok=True)
PHASE_C = Path("/home/seungho/personal/THCA_data_analysis/project/results/paper11_pancancer/phase_C_ICI/processed")

HT13 = ["HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1",
        "CD79A", "CD79B", "MS4A1", "AICDA", "CXCL13", "CCR6", "IFNG"]

# Tumour Inflammation Signature (Ayers 2017) 18-gene
TIS18 = ["CD27", "CD274", "CD276", "CD8A", "CMKLR1", "CXCL9", "CXCR6",
         "HLA-DQA1", "HLA-DRB1", "HLA-E", "IDO1", "LAG3", "NKG7",
         "PDCD1LG2", "PSMB10", "STAT1", "TIGIT", "CCL5"]

# HLA class I + APM core (for DM1-like-adaptive labelling)
HLA1 = ["HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2"]


def zrows_within_cohort(expr_df: pd.DataFrame) -> pd.DataFrame:
    """Z-score each gene row across samples (cohort-internal)."""
    mu = expr_df.mean(axis=1)
    sd = expr_df.std(axis=1).replace(0, np.nan)
    z = expr_df.sub(mu, axis=0).div(sd, axis=0)
    return z


def panel_score(zexpr: pd.DataFrame, genes: list[str]) -> tuple[pd.Series, list[str], list[str]]:
    avail = [g for g in genes if g in zexpr.index]
    miss = [g for g in genes if g not in zexpr.index]
    if not avail:
        return pd.Series(dtype=float), [], miss
    return zexpr.loc[avail].mean(axis=0), avail, miss


def load_phase_c_cohort(name: str) -> tuple[pd.DataFrame, pd.DataFrame] | None:
    p = PHASE_C / name
    e = p / "expr_matrix.tsv"
    m = p / "metadata.tsv"
    if not (e.exists() and m.exists()):
        return None
    expr = pd.read_csv(e, sep="\t", index_col=0)
    meta = pd.read_csv(m, sep="\t")
    return expr, meta


def load_hugo() -> tuple[pd.DataFrame, pd.DataFrame]:
    src = Path("/home/seungho/personal/THCA_data_analysis/project/results/paper3_ici_public_data/GSE78220")
    expr = pd.read_excel(src / "GSE78220_PatientFPKM.xlsx").set_index("Gene")
    expr = np.log2(expr + 1.0)
    # Pre-treatment baseline columns are "PtX.baseline"; also keep any "PtX.OnTx"? Use baselines only.
    base_cols = [c for c in expr.columns if c.endswith(".baseline")]
    expr = expr[base_cols]
    expr.columns = [c.replace(".baseline", "") for c in expr.columns]

    # Parse response from series matrix
    resp_map = {}
    with gzip.open(src / "GSE78220_series_matrix.txt.gz", "rt") as f:
        text = f.read()
    title_line = next(l for l in text.split("\n") if l.startswith("!Sample_title"))
    titles = [t.strip().strip('"') for t in title_line.split("\t")[1:]]
    src_line = next(l for l in text.split("\n") if l.startswith("!Sample_source_name_ch1"))
    sources = [s.strip().strip('"') for s in src_line.split("\t")[1:]]
    for t, s in zip(titles, sources):
        if "Complete Response" in s or "Partial Response" in s:
            resp_map[t] = ("CR" if "Complete" in s else "PR", 1)
        elif "Progressive Disease" in s:
            resp_map[t] = ("PD", 0)
        else:
            resp_map[t] = (s, np.nan)
    # Match to expr columns; some titles like Pt27A/Pt27B
    rows = []
    for c in expr.columns:
        cc = c
        if cc not in resp_map:
            # try strip suffix letters e.g. Pt27A -> Pt27
            stripped = re.sub(r"[A-Za-z]+$", "", cc)
            if stripped in resp_map and stripped != cc:
                cc = stripped
        if cc in resp_map:
            rraw, rbin = resp_map[cc]
            rows.append({"sample_id": c, "response_raw": rraw,
                         "response_binary_CRPR_vs_SD_PD": rbin,
                         "cohort": "Hugo_GSE78220",
                         "cancer_type": "melanoma",
                         "therapy": "anti-PD-1",
                         "timepoint": "pre"})
        else:
            rows.append({"sample_id": c, "response_raw": "NA",
                         "response_binary_CRPR_vs_SD_PD": np.nan,
                         "cohort": "Hugo_GSE78220",
                         "cancer_type": "melanoma",
                         "therapy": "anti-PD-1",
                         "timepoint": "pre"})
    meta = pd.DataFrame(rows)
    return expr, meta


def filter_pre_responder(meta: pd.DataFrame) -> pd.DataFrame:
    m = meta.copy()
    if "timepoint" in m.columns:
        m = m[m["timepoint"].astype(str).str.lower().isin(["pre", "baseline", "nan", ""]) | m["timepoint"].isna()]
    m = m[~m["response_binary_CRPR_vs_SD_PD"].isna()]
    return m


def logistic_or(score: np.ndarray, y: np.ndarray) -> dict:
    """Wald logistic regression on standardized score → binary y."""
    from sklearn.linear_model import LogisticRegression
    score = np.asarray(score, float)
    y = np.asarray(y, int)
    mask = ~np.isnan(score) & ~np.isnan(y)
    score, y = score[mask], y[mask]
    if len(np.unique(y)) < 2 or len(y) < 8:
        return dict(OR=np.nan, lo=np.nan, hi=np.nan, p=np.nan, n=int(len(y)),
                    n_resp=int(y.sum()), AUC=np.nan, beta=np.nan, se=np.nan)
    z = (score - score.mean()) / (score.std(ddof=1) + 1e-12)
    X = z.reshape(-1, 1)
    lr = LogisticRegression(C=1e6, solver="lbfgs", max_iter=500).fit(X, y)
    beta = lr.coef_[0, 0]
    p_hat = lr.predict_proba(X)[:, 1]
    W = np.diag(p_hat * (1 - p_hat))
    Xd = np.hstack([np.ones((len(z), 1)), X])
    try:
        cov = np.linalg.inv(Xd.T @ W @ Xd)
        se_beta = float(np.sqrt(cov[1, 1]))
    except np.linalg.LinAlgError:
        se_beta = float("nan")
    wald = beta / se_beta if se_beta and not np.isnan(se_beta) else np.nan
    p = 2 * (1 - stats.norm.cdf(abs(wald))) if not np.isnan(wald) else np.nan
    OR = float(np.exp(beta))
    lo = float(np.exp(beta - 1.96 * se_beta)) if not np.isnan(se_beta) else np.nan
    hi = float(np.exp(beta + 1.96 * se_beta)) if not np.isnan(se_beta) else np.nan
    try:
        from sklearn.metrics import roc_auc_score
        auc = float(roc_auc_score(y, score))
    except Exception:
        auc = float("nan")
    return dict(OR=OR, lo=lo, hi=hi, p=float(p), n=int(len(y)),
                n_resp=int(y.sum()), AUC=auc, beta=float(beta), se=float(se_beta))


def cox_continuous(score: np.ndarray, time: np.ndarray, event: np.ndarray) -> dict:
    """Univariate Cox regression on standardized continuous score."""
    score = np.asarray(score, float)
    time = pd.to_numeric(time, errors="coerce").to_numpy() if hasattr(time, "to_numpy") else np.asarray(time, float)
    event = pd.to_numeric(event, errors="coerce").to_numpy() if hasattr(event, "to_numpy") else np.asarray(event, float)
    mask = ~np.isnan(score) & ~np.isnan(time) & ~np.isnan(event) & (time > 0)
    score, time, event = score[mask], time[mask], event[mask]
    if len(score) < 12 or event.sum() < 4:
        return dict(HR=np.nan, lo=np.nan, hi=np.nan, p=np.nan, n=int(len(score)), events=int(event.sum()))
    try:
        from lifelines import CoxPHFitter
        z = (score - score.mean()) / (score.std(ddof=1) + 1e-12)
        df = pd.DataFrame({"t": time, "e": event.astype(int), "z": z})
        cph = CoxPHFitter()
        cph.fit(df, duration_col="t", event_col="e")
        s = cph.summary.iloc[0]
        return dict(HR=float(s["exp(coef)"]),
                    lo=float(s["exp(coef) lower 95%"]),
                    hi=float(s["exp(coef) upper 95%"]),
                    p=float(s["p"]), n=int(len(score)), events=int(event.sum()))
    except Exception as ex:
        return dict(HR=np.nan, lo=np.nan, hi=np.nan, p=np.nan, n=int(len(score)),
                    events=int(event.sum()), error=str(ex)[:100])


def random_effects_meta(beta: np.ndarray, se: np.ndarray) -> dict:
    """DerSimonian-Laird random-effects pooled OR."""
    beta = np.asarray(beta, float); se = np.asarray(se, float)
    mask = ~np.isnan(beta) & ~np.isnan(se) & (se > 0)
    beta, se = beta[mask], se[mask]
    k = len(beta)
    if k < 2:
        return dict(pooled_OR=np.nan, lo=np.nan, hi=np.nan, p=np.nan, I2=np.nan, k=k, tau2=np.nan)
    w_fe = 1.0 / (se ** 2)
    mu_fe = (w_fe * beta).sum() / w_fe.sum()
    Q = (w_fe * (beta - mu_fe) ** 2).sum()
    df = k - 1
    C = w_fe.sum() - (w_fe ** 2).sum() / w_fe.sum()
    tau2 = max(0.0, (Q - df) / C) if C > 0 else 0.0
    w_re = 1.0 / (se ** 2 + tau2)
    mu_re = (w_re * beta).sum() / w_re.sum()
    se_re = float(np.sqrt(1.0 / w_re.sum()))
    z = mu_re / se_re
    p = 2 * (1 - stats.norm.cdf(abs(z)))
    I2 = max(0.0, (Q - df) / Q) * 100 if Q > 0 else 0.0
    return dict(pooled_OR=float(np.exp(mu_re)),
                lo=float(np.exp(mu_re - 1.96 * se_re)),
                hi=float(np.exp(mu_re + 1.96 * se_re)),
                p=float(p), I2=float(I2), k=int(k),
                tau2=float(tau2), Q=float(Q))


# ---------------- driver ----------------

cohorts_spec = [
    ("IMvigor210", "phase_c"),
    ("riaz_GSE91061", "phase_c"),
    ("MGH_GSE115821", "phase_c"),
    ("GSE176307", "phase_c"),
    ("Hugo_GSE78220", "hugo"),
]

per_cohort_rows = []
adaptive_rows = []
per_sample_blocks = []
beta_logor = {"HT13": [], "TIS": [], "DM1_inflam": []}
se_logor = {"HT13": [], "TIS": [], "DM1_inflam": []}
cohort_names = []

for cname, kind in cohorts_spec:
    print(f"\n=== {cname} ===")
    if kind == "phase_c":
        load = load_phase_c_cohort(cname)
        if load is None:
            print("  missing — skip")
            continue
        expr, meta = load
    else:
        expr, meta = load_hugo()

    # restrict meta to pre + with response
    m = filter_pre_responder(meta)
    if "patient_id" in m.columns:
        # IMvigor210 has only one timepoint per patient already; still use sample_id keys
        pass
    print(f"  n_samples_with_response={len(m)}, n_resp={int(m['response_binary_CRPR_vs_SD_PD'].sum())}")

    sids = [s for s in m["sample_id"].astype(str).tolist() if s in expr.columns]
    if len(sids) < 8:
        print(f"  too few matched samples ({len(sids)}) — skip")
        continue
    expr_s = expr[sids]
    # drop all-zero / NA gene rows; coerce numeric
    expr_s = expr_s.apply(pd.to_numeric, errors="coerce")
    expr_s = expr_s.dropna(how="all").fillna(0.0)
    z = zrows_within_cohort(expr_s)

    ht13_score, ht13_avail, ht13_miss = panel_score(z, HT13)
    tis_score, tis_avail, tis_miss = panel_score(z, TIS18)
    hla1_score, hla1_avail, _ = panel_score(z, HLA1)

    if ht13_score.empty or tis_score.empty:
        print("  no panel genes available — skip"); continue

    print(f"  HT13 genes avail={len(ht13_avail)}/13 missing={ht13_miss}")
    print(f"  TIS genes avail={len(tis_avail)}/18 missing={tis_miss}")
    print(f"  HLA-I genes avail={len(hla1_avail)}/{len(HLA1)}")

    dm1_inflam = (ht13_score + tis_score) / 2.0

    # align with metadata response
    mm = m.set_index("sample_id").loc[sids].copy()
    mm["HT13"] = ht13_score
    mm["TIS"] = tis_score
    mm["DM1_inflam"] = dm1_inflam
    mm["HLA1"] = hla1_score.reindex(sids) if not hla1_score.empty else np.nan
    mm["cohort_r5"] = cname
    mm["n_HT13_genes"] = len(ht13_avail)
    mm["n_TIS_genes"] = len(tis_avail)
    per_sample_blocks.append(mm.reset_index())

    y = mm["response_binary_CRPR_vs_SD_PD"].astype(float).to_numpy()

    for panel_name, scores in [("HT13", mm["HT13"].to_numpy()),
                               ("TIS", mm["TIS"].to_numpy()),
                               ("DM1_inflam", mm["DM1_inflam"].to_numpy())]:
        res = logistic_or(scores, y)
        per_cohort_rows.append(dict(cohort=cname, panel=panel_name, **res,
                                    n_HT13_genes=len(ht13_avail),
                                    n_TIS_genes=len(tis_avail)))
        if not np.isnan(res["beta"]) and not np.isnan(res["se"]):
            beta_logor[panel_name].append(res["beta"])
            se_logor[panel_name].append(res["se"])

    # Cox PFS / OS where present
    for surv_col, ev_col, lab in [("pfs_time", "pfs_event", "PFS"),
                                  ("os_time", "os_event", "OS")]:
        if surv_col in mm.columns and ev_col in mm.columns:
            t = pd.to_numeric(mm[surv_col], errors="coerce")
            e = pd.to_numeric(mm[ev_col], errors="coerce")
            for panel_name, scores in [("HT13", mm["HT13"].to_numpy()),
                                       ("DM1_inflam", mm["DM1_inflam"].to_numpy())]:
                cox = cox_continuous(scores, t, e)
                per_cohort_rows.append(dict(cohort=cname, panel=panel_name + f"_{lab}",
                                            OR=np.nan, lo=cox["lo"], hi=cox["hi"],
                                            p=cox["p"], n=cox["n"],
                                            n_resp=cox["events"],
                                            AUC=np.nan, beta=np.nan, se=np.nan,
                                            HR=cox["HR"], surv_metric=lab,
                                            n_HT13_genes=len(ht13_avail),
                                            n_TIS_genes=len(tis_avail)))

    # DM1-like-adaptive label = HT-13 high AND HLA-I high
    if not hla1_score.empty:
        ht_med = float(np.nanmedian(mm["HT13"]))
        hla_med = float(np.nanmedian(mm["HLA1"]))
        adaptive = (mm["HT13"] > ht_med) & (mm["HLA1"] > hla_med)
        # Fisher exact 2x2: adaptive vs response
        a = int(((adaptive == True) & (y == 1)).sum())
        b = int(((adaptive == True) & (y == 0)).sum())
        c = int(((adaptive == False) & (y == 1)).sum())
        d = int(((adaptive == False) & (y == 0)).sum())
        try:
            OR_f, p_f = stats.fisher_exact([[a, b], [c, d]])
        except Exception:
            OR_f, p_f = np.nan, np.nan
        # response rate split
        rr_adapt = a / (a + b) if (a + b) > 0 else np.nan
        rr_other = c / (c + d) if (c + d) > 0 else np.nan
        adaptive_rows.append(dict(cohort=cname,
                                  n_adaptive=int(adaptive.sum()),
                                  n_other=int((~adaptive).sum()),
                                  resp_rate_adaptive=float(rr_adapt) if rr_adapt == rr_adapt else np.nan,
                                  resp_rate_other=float(rr_other) if rr_other == rr_other else np.nan,
                                  fisher_OR=float(OR_f) if OR_f == OR_f else np.nan,
                                  fisher_p=float(p_f) if p_f == p_f else np.nan,
                                  a_adapt_resp=a, b_adapt_nonresp=b,
                                  c_other_resp=c, d_other_nonresp=d,
                                  hla1_avail=len(hla1_avail)))
    cohort_names.append(cname)


# ---------------- pooled meta ----------------

pooled = {}
forest_rows = []
for panel_name in ["HT13", "TIS", "DM1_inflam"]:
    res = random_effects_meta(np.array(beta_logor[panel_name]),
                              np.array(se_logor[panel_name]))
    pooled[panel_name] = res
    for cname in cohort_names:
        # find row
        sub = [r for r in per_cohort_rows
               if r["cohort"] == cname and r["panel"] == panel_name]
        if sub:
            r = sub[0]
            forest_rows.append(dict(panel=panel_name, cohort=cname,
                                    OR=r["OR"], lo=r["lo"], hi=r["hi"],
                                    beta=r["beta"], se=r["se"], p=r["p"],
                                    n=r["n"], n_resp=r["n_resp"], AUC=r["AUC"]))
    forest_rows.append(dict(panel=panel_name, cohort="POOLED_RE",
                            OR=res.get("pooled_OR"),
                            lo=res.get("lo"), hi=res.get("hi"),
                            beta=np.log(res["pooled_OR"]) if res.get("pooled_OR") and not np.isnan(res["pooled_OR"]) else np.nan,
                            se=np.nan, p=res.get("p"),
                            n=int(np.sum([r["n"] for r in per_cohort_rows
                                          if r["panel"] == panel_name and r["cohort"] != "POOLED_RE"])),
                            n_resp=int(np.sum([r["n_resp"] for r in per_cohort_rows
                                               if r["panel"] == panel_name and r["cohort"] != "POOLED_RE"])),
                            AUC=np.nan, I2=res.get("I2"), k=res.get("k")))


# ---------------- save ----------------

per_cohort_df = pd.DataFrame(per_cohort_rows)
per_cohort_df.to_csv(OUT / "r5_per_cohort_HT13_response.tsv", sep="\t", index=False)

forest_df = pd.DataFrame(forest_rows)
forest_df.to_csv(OUT / "r5_pooled_ici_meta.tsv", sep="\t", index=False)

if adaptive_rows:
    pd.DataFrame(adaptive_rows).to_csv(OUT / "r5_dm1_adaptive_responder.tsv",
                                       sep="\t", index=False)

if per_sample_blocks:
    pd.concat(per_sample_blocks, ignore_index=True).to_csv(
        OUT / "r5_per_sample_scores.tsv", sep="\t", index=False)

with open(OUT / "r5_headline.json", "w") as f:
    json.dump(dict(pooled=pooled,
                   cohorts=cohort_names,
                   k=len(cohort_names)), f, indent=2, default=str)

print("\n=== POOLED ===")
for k, v in pooled.items():
    if v.get("pooled_OR") and not np.isnan(v["pooled_OR"]):
        print(f"  {k}: pooled OR={v['pooled_OR']:.3f} [{v['lo']:.3f}, {v['hi']:.3f}] p={v['p']:.3g} I2={v['I2']:.1f}% k={v['k']}")
    else:
        print(f"  {k}: NaN ({v})")

print("\nWrote:")
for fp in OUT.iterdir():
    print(f"  {fp.name}")
