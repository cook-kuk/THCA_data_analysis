#!/usr/bin/env python3
from __future__ import annotations

import numpy as np
import pandas as pd
import plotly.express as px
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test
from scipy import stats

from v17_common import ROOT
from v17p2_common import stage_to_ord, tcga12
from v17p35_common import FIG35, LOG35, TAB35, THYROID_DIFF, load_dark, log_line, score_expr, write_json


def fisher_binary(df: pd.DataFrame, endpoint: str) -> tuple[float, float, int]:
    tmp = df[["v17_dark_cluster", endpoint]].dropna().copy()
    if tmp.empty or tmp[endpoint].nunique() < 2:
        return np.nan, np.nan, 0
    tab = pd.crosstab(tmp["v17_dark_cluster"], tmp[endpoint])
    if tab.shape != (2, 2):
        return np.nan, np.nan, int(tmp.shape[0])
    odds, p = stats.fisher_exact(tab.to_numpy())
    return float(odds), float(p), int(tmp.shape[0])


def main() -> None:
    log_line(LOG35, "v17p35 FIX3 start")
    expr, meta = load_dark()
    expr.columns = expr.columns.astype(str).str.upper()
    clin = pd.read_csv(ROOT / "results" / "tables" / "tcga_thca_clinical_extended.tsv", sep="\t")
    clin["tcga12"] = clin["sample_id"].astype(str).map(tcga12)
    df = meta.copy()
    df["tcga12"] = df["sample_id"].astype(str).map(tcga12)
    df = df.merge(clin, on="tcga12", how="left", suffixes=("", "_clin"))
    df["cluster_bin"] = (df["v17_dark_cluster"] == "DM2").astype(int)
    df["stage_ord"] = df["ajcc_stage_group"].map(stage_to_ord)
    df["stage_high"] = df["stage_ord"] >= 3
    df["sex_bin"] = df["sex"].astype(str).str.lower().map({"male": 1, "female": 0})
    df["histology_fvptc"] = (df["histology_subtype"] == "FVPTC").astype(int)
    df["rai_score_recalc"] = score_expr(expr, [g.upper() for g in THYROID_DIFF]).values
    df["os_event"] = pd.to_numeric(df["os_event"], errors="coerce")
    df["os_days"] = pd.to_numeric(df["os_days"], errors="coerce")
    df["tumor_size_mm"] = pd.to_numeric(df["tumor_size_mm"], errors="coerce")
    df.to_csv(TAB35 / "FIX3_cdr_clinical_full.tsv", sep="\t", index=False)

    rows = []
    # continuous endpoints
    for endpoint in ["age", "stage_ord", "tumor_size_mm", "rai_score_recalc", "tds16_score_v17"]:
        tmp = df[["v17_dark_cluster", endpoint]].dropna()
        if tmp.empty:
            continue
        a = pd.to_numeric(tmp.loc[tmp["v17_dark_cluster"] == "DM1", endpoint], errors="coerce")
        b = pd.to_numeric(tmp.loc[tmp["v17_dark_cluster"] == "DM2", endpoint], errors="coerce")
        if a.notna().sum() < 3 or b.notna().sum() < 3:
            continue
        t = stats.ttest_ind(a, b, equal_var=False, nan_policy="omit")
        rows.append(
            {
                "endpoint": endpoint,
                "endpoint_type": "continuous",
                "n_total": int(tmp.shape[0]),
                "n_events_or_positive": np.nan,
                "test_method": "Welch_t",
                "statistic": float(t.statistic),
                "p_value": float(t.pvalue),
                "effect_size": float(b.mean() - a.mean()),
                "effect_size_CI": np.nan,
            }
        )
    # binary endpoints
    for endpoint in ["stage_high", "os_event", "histology_fvptc"]:
        odds, p, n = fisher_binary(df, endpoint)
        rows.append(
            {
                "endpoint": endpoint,
                "endpoint_type": "binary",
                "n_total": n,
                "n_events_or_positive": int(df[endpoint].fillna(False).astype(bool).sum()) if endpoint in df.columns else np.nan,
                "test_method": "Fisher_exact",
                "statistic": odds,
                "p_value": p,
                "effect_size": odds,
                "effect_size_CI": np.nan,
            }
        )
    # survival log-rank
    tmp = df[["v17_dark_cluster", "os_days", "os_event"]].dropna()
    if not tmp.empty and tmp["os_event"].sum() > 0:
        d1 = tmp[tmp["v17_dark_cluster"] == "DM1"]
        d2 = tmp[tmp["v17_dark_cluster"] == "DM2"]
        lr = logrank_test(d1["os_days"], d2["os_days"], event_observed_A=d1["os_event"], event_observed_B=d2["os_event"])
        rows.append(
            {
                "endpoint": "OS_logrank",
                "endpoint_type": "survival",
                "n_total": int(tmp.shape[0]),
                "n_events_or_positive": int(tmp["os_event"].sum()),
                "test_method": "logrank",
                "statistic": float(lr.test_statistic),
                "p_value": float(lr.p_value),
                "effect_size": np.nan,
                "effect_size_CI": np.nan,
            }
        )
    alt = pd.DataFrame(rows)
    alt.to_csv(TAB35 / "FIX3_alternative_endpoints_full.tsv", sep="\t", index=False)

    # Cox
    cox_candidates = [
        ["os_days", "os_event", "cluster_bin", "age", "stage_ord"],
        ["os_days", "os_event", "cluster_bin", "age", "sex_bin"],
        ["os_days", "os_event", "cluster_bin", "age"],
        ["os_days", "os_event", "cluster_bin"],
    ]
    cox_out = pd.DataFrame()
    for cox_cols in cox_candidates:
        cox_df = df[cox_cols].dropna().copy()
        if not cox_df.empty and cox_df["os_event"].sum() > 1:
            try:
                cph = CoxPHFitter()
                cph.fit(cox_df, duration_col="os_days", event_col="os_event")
                summ = cph.summary.reset_index().rename(columns={"covariate": "feature"})
                cox_out = summ[["feature", "coef", "exp(coef)", "se(coef)", "p", "coef lower 95%", "coef upper 95%", "exp(coef) lower 95%", "exp(coef) upper 95%"]].copy()
                break
            except Exception:
                continue
    if not cox_out.empty:
        cox_out.to_csv(TAB35 / "FIX3_cox_multivariable.tsv", sep="\t", index=False)
    else:
        pd.DataFrame(columns=["feature", "coef", "exp(coef)", "se(coef)", "p"]).to_csv(TAB35 / "FIX3_cox_multivariable.tsv", sep="\t", index=False)

    km_rows = []
    if not tmp.empty:
        km_rows.append({"endpoint": "OS", "n_total": int(tmp.shape[0]), "n_events": int(tmp["os_event"].sum()), "p_value": float(rows[-1]["p_value"]) if rows and rows[-1]["endpoint"] == "OS_logrank" else np.nan})
    pd.DataFrame(km_rows).to_csv(TAB35 / "FIX3_km_logrank_per_endpoint.tsv", sep="\t", index=False)

    # figures
    if not alt.empty:
        volc = alt.copy()
        volc["mlog10p"] = -np.log10(np.clip(pd.to_numeric(volc["p_value"], errors="coerce").fillna(1.0), 1e-300, 1.0))
        px.scatter(volc, x="effect_size", y="mlog10p", color="endpoint_type", hover_name="endpoint", title="Clinical endpoint volcano").write_html(
            FIG35 / "FIX3_endpoint_volcano.html", include_plotlyjs="cdn"
        )
    if not cox_out.empty:
        px.scatter(cox_out, x="exp(coef)", y="feature", error_x=cox_out["exp(coef) upper 95%"] - cox_out["exp(coef)"], title="Cox multivariable forest").write_html(
            FIG35 / "FIX3_cox_forest.html", include_plotlyjs="cdn"
        )
    if not tmp.empty:
        kmf = KaplanMeierFitter()
        plot_rows = []
        for grp in ["DM1", "DM2"]:
            sub = tmp[tmp["v17_dark_cluster"] == grp]
            kmf.fit(sub["os_days"], sub["os_event"], label=grp)
            surv = kmf.survival_function_.reset_index()
            surv.columns = ["timeline", "survival"]
            surv["cluster"] = grp
            plot_rows.append(surv)
        if plot_rows:
            plot_df = pd.concat(plot_rows, ignore_index=True)
            px.line(plot_df, x="timeline", y="survival", color="cluster", title="OS Kaplan-Meier").write_html(
                FIG35 / "FIX3_km_panel.html", include_plotlyjs="cdn"
            )

    payload = {
        "n_endpoints_p_lt_0_05": int((pd.to_numeric(alt["p_value"], errors="coerce") < 0.05).sum()) if not alt.empty else 0,
        "cox_cluster_hr": float(cox_out.loc[cox_out["feature"] == "cluster_bin", "exp(coef)"].iloc[0]) if not cox_out.empty and (cox_out["feature"] == "cluster_bin").any() else None,
    }
    write_json(TAB35 / "FIX3_summary.json", payload)
    log_line(LOG35, f"v17p35 FIX3 done {payload}")


if __name__ == "__main__":
    main()
