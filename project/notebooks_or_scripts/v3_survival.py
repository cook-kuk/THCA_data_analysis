"""v3 STEP 5 — Survival analysis (TCGA-THCA).

KM curves, log-rank pairwise, univariate + multivariate Cox (adjust age+stage+sex).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import multivariate_logrank_test, logrank_test

sys.path.insert(0, "/opt/thyroid-dash/project/notebooks_or_scripts")
from v3_common import (DATASET_FILES, FIGS, RESULTS_TABLES, get_tds16,
                        load_expr, load_sample_master, small_n_warning,
                        write_tsv)

PROJECT = Path("/opt/thyroid-dash/project")
CLINICAL = PROJECT / "results" / "tables" / "tcga_thca_clinical_extended.tsv"


def load_clinical_merged():
    clin = pd.read_csv(CLINICAL, sep="\t")
    # sample_id column uses -01A form; master uses same
    sm = load_sample_master()
    sm_t = sm[(sm.dataset == "TCGA-THCA") & (sm.normal_vs_tumor == "tumor")].copy()
    merged = sm_t.merge(clin, on="sample_id", how="left", suffixes=("", "_c"))
    # build duration and event
    merged["os_days"] = merged["os_days"].astype(float)
    merged["os_event"] = merged["os_event"].fillna(0).astype(int)
    return merged


def compute_dediff_proxy_score(merged):
    """If `dedifferentiation_proxy_score` is present in sample_master, use it.
    Otherwise compute a simple TDS16-based decile (low TDS -> high dediff).
    """
    if "dedifferentiation_proxy_score" in merged.columns and \
            merged["dedifferentiation_proxy_score"].notna().any():
        return merged["dedifferentiation_proxy_score"].astype(float)

    expr = load_expr(DATASET_FILES["TCGA-THCA"])
    tds = [g for g in get_tds16() if g in expr.index]
    z = ((expr.loc[tds] - expr.loc[tds].mean(axis=1).values.reshape(-1, 1))
         / expr.loc[tds].std(axis=1).values.reshape(-1, 1))
    score = -z.mean(axis=0)  # high -> low TDS -> dediff-proxy
    return merged["sample_id"].map(score.to_dict()).astype(float)


def tds16_tertile(merged):
    expr = load_expr(DATASET_FILES["TCGA-THCA"])
    tds = [g for g in get_tds16() if g in expr.index]
    samps = [s for s in merged["sample_id"] if s in expr.columns]
    vals = expr.loc[tds, samps].mean(axis=0)
    bins = pd.qcut(vals, q=3, labels=["low", "mid", "high"])
    return merged["sample_id"].map(bins.to_dict())


def km_plot(durations, events, groups, title, out_html):
    fig = go.Figure()
    kmfs = {}
    for g in groups.dropna().unique():
        idx = groups == g
        n_g = int(idx.sum())
        events_g = int(events[idx].sum())
        if n_g < 3:
            continue
        kmf = KaplanMeierFitter()
        kmf.fit(durations[idx], events[idx], label=f"{g} (n={n_g}, d={events_g})")
        kmfs[g] = kmf
        sf = kmf.survival_function_.reset_index()
        fig.add_scatter(x=sf["timeline"],
                         y=sf[kmf._label].values,
                         mode="lines",
                         name=f"{g} (n={n_g}, d={events_g}"
                              f"{small_n_warning(events_g, 10)})")

    # log-rank
    try:
        valid = groups.dropna()
        lr = multivariate_logrank_test(durations.loc[valid.index],
                                        valid,
                                        events.loc[valid.index])
        subtitle = f"logrank p={lr.p_value:.3g}"
    except Exception as e:
        subtitle = f"logrank failed: {e}"

    fig.update_layout(title=f"{title}<br><sub>{subtitle}</sub>",
                      xaxis_title="Days", yaxis_title="Survival probability",
                      template="plotly_dark", height=500,
                      yaxis_range=[0, 1.02])
    fig.write_html(out_html, include_plotlyjs="cdn")
    return {"logrank_p": getattr(lr, "p_value", np.nan) if "lr" in dir() else np.nan}


def run():
    merged = load_clinical_merged()

    # filter: have survival
    merged = merged[merged["os_days"].notna() & (merged["os_days"] > 0)].copy()
    durations = merged["os_days"]
    events = merged["os_event"]

    # (a) molecular_subtype stratification
    sub = merged[merged["molecular_subtype"].isin(
        ["BRAF_like", "RAS_like", "dedifferentiated"])].copy()
    groups_mol = sub["molecular_subtype"]
    km_plot(sub["os_days"], sub["os_event"], groups_mol,
            title="OS by molecular subtype (TCGA-THCA)",
            out_html=FIGS / "v3_km_subtype.html")

    # (b) TDS16 tertile
    merged["tds_tertile"] = tds16_tertile(merged)
    km_plot(durations, events, merged["tds_tertile"],
            title="OS by TDS16 tertile",
            out_html=FIGS / "v3_km_tds_tertile.html")

    # (c) dedifferentiation-proxy tertile
    merged["dediff_score"] = compute_dediff_proxy_score(merged)
    merged["dediff_tertile"] = pd.qcut(merged["dediff_score"],
                                       q=3, labels=["low", "mid", "high"])
    km_plot(durations, events, merged["dediff_tertile"],
            title="OS by dediff-proxy tertile",
            out_html=FIGS / "v3_km_dediff_tertile.html")

    # Cox (univariate and multivariate)
    cox_rows = []
    # Univariate
    def cox_df(df, cols):
        d = df[cols + ["os_days", "os_event"]].copy()
        d = d.dropna()
        d["os_days"] = d["os_days"].astype(float)
        d["os_event"] = d["os_event"].astype(int)
        return d

    # Prepare multivariate covariates
    # encode stage as numeric
    stage_map = {}
    for s in merged["stage"].dropna().unique():
        lvl = 0
        if "IV" in s:
            lvl = 4
        elif "III" in s:
            lvl = 3
        elif "II" in s:
            lvl = 2
        elif "I" in s:
            lvl = 1
        stage_map[s] = lvl
    merged["stage_num"] = merged["stage"].map(stage_map)
    merged["sex_num"] = merged["gender"].map(
        {"female": 0, "male": 1, "Female": 0, "Male": 1})

    # Univariate features
    feats_uni = ["tds_tertile", "dediff_tertile", "molecular_subtype",
                  "stage_num", "sex_num", "age_at_diagnosis"]
    for f in feats_uni:
        try:
            d = merged[[f, "os_days", "os_event"]].dropna()
            if d[f].dtype.name in ("object", "category"):
                # dummy
                dum = pd.get_dummies(d[f], prefix=f, drop_first=True, dtype=float)
                dd = pd.concat([dum, d[["os_days", "os_event"]]], axis=1)
            else:
                dd = d.copy()
            if len(dd) < 10 or dd["os_event"].sum() < 3:
                cox_rows.append(dict(covariate=f, mode="univariate", n=len(dd),
                                      hr=np.nan, lo=np.nan, hi=np.nan,
                                      p=np.nan, note="insufficient events"))
                continue
            cph = CoxPHFitter()
            cph.fit(dd, duration_col="os_days", event_col="os_event")
            for c in cph.summary.index:
                row = cph.summary.loc[c]
                cox_rows.append(dict(covariate=c, mode="univariate", n=len(dd),
                                      hr=float(row["exp(coef)"]),
                                      lo=float(row["exp(coef) lower 95%"]),
                                      hi=float(row["exp(coef) upper 95%"]),
                                      p=float(row["p"]), note=""))
        except Exception as e:
            cox_rows.append(dict(covariate=f, mode="univariate", n=0,
                                  hr=np.nan, lo=np.nan, hi=np.nan,
                                  p=np.nan, note=f"error: {e}"))

    # Multivariate: molecular_subtype + stage_num + age + sex_num
    try:
        d = merged[["molecular_subtype", "stage_num", "age_at_diagnosis",
                    "sex_num", "os_days", "os_event"]].dropna()
        d = d[d["molecular_subtype"].isin(["BRAF_like", "RAS_like",
                                             "dedifferentiated"])]
        dum = pd.get_dummies(d["molecular_subtype"], prefix="subtype",
                             drop_first=True, dtype=float)
        dd = pd.concat([dum,
                        d[["stage_num", "age_at_diagnosis", "sex_num",
                           "os_days", "os_event"]]], axis=1)
        if dd["os_event"].sum() >= 3 and len(dd) > 10:
            cph = CoxPHFitter()
            cph.fit(dd, duration_col="os_days", event_col="os_event")
            for c in cph.summary.index:
                row = cph.summary.loc[c]
                cox_rows.append(dict(covariate=c, mode="multivariate",
                                      n=len(dd),
                                      hr=float(row["exp(coef)"]),
                                      lo=float(row["exp(coef) lower 95%"]),
                                      hi=float(row["exp(coef) upper 95%"]),
                                      p=float(row["p"]),
                                      note=("⚠ events<10"
                                            if dd["os_event"].sum() < 10 else "")))
    except Exception as e:
        cox_rows.append(dict(covariate="multivariate", mode="multivariate",
                              n=0, hr=np.nan, lo=np.nan, hi=np.nan,
                              p=np.nan, note=f"error: {e}"))

    cox_df_out = pd.DataFrame(cox_rows)
    write_tsv(cox_df_out, RESULTS_TABLES / "v3_survival_cox.tsv")

    # Cox forest plot (multivariate only)
    mv = cox_df_out[cox_df_out["mode"] == "multivariate"].copy()
    if not mv.empty:
        fig = go.Figure()
        mv = mv.sort_values("covariate")
        fig.add_scatter(x=mv["hr"], y=mv["covariate"], mode="markers",
                         marker=dict(color="#5eead4", size=12),
                         error_x=dict(type="data", symmetric=False,
                                      array=mv["hi"] - mv["hr"],
                                      arrayminus=mv["hr"] - mv["lo"]),
                         name="HR (95% CI)")
        fig.add_vline(x=1.0, line=dict(dash="dash", color="#94a3b8"))
        total_events = int(merged["os_event"].sum())
        fig.update_layout(
            title=f"Multivariate Cox forest (TCGA-THCA, events={total_events}"
                  f"{small_n_warning(total_events, 10)})",
            xaxis_title="Hazard ratio (log)",
            xaxis_type="log",
            template="plotly_dark", height=480)
        fig.write_html(FIGS / "v3_cox_forest.html", include_plotlyjs="cdn")

    # Summary
    n_samples = int(len(merged))
    n_events = int(merged["os_event"].sum())
    return dict(n_samples=n_samples, n_events=n_events,
                 flag=small_n_warning(n_events, 10),
                 out=str(RESULTS_TABLES / "v3_survival_cox.tsv"))


def main():
    res = run()
    print(json.dumps(res, indent=2, default=str))


if __name__ == "__main__":
    main()
