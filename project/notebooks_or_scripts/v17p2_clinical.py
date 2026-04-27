#!/usr/bin/env python3
from __future__ import annotations

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from lifelines import CoxPHFitter, KaplanMeierFitter
from scipy import stats

from v17p2_common import FIG, LOG, TAB, dm_binary_labels, fisher_2x2, load_clinical, load_dark_cohort, log_line, stage_to_ord, tcga12, write_json


def main() -> None:
    log_line(LOG, "v17p2 layer2 start")
    _, meta = load_dark_cohort()
    clin = load_clinical()
    meta = meta.copy()
    meta["tcga12"] = meta["sample_id"].astype(str).map(tcga12)
    dm = meta.merge(clin, on="tcga12", how="left", suffixes=("", "_clin"))
    dm["cluster_bin"] = (dm["v17_dark_cluster"] == "DM2").astype(int)
    dm["stage_ord"] = dm["ajcc_stage_group"].map(stage_to_ord)
    dm["stage_high"] = dm["stage_ord"].ge(3).fillna(False).astype(int)
    dm["os_event_use"] = dm["os_event"].fillna(0).astype(int)
    dm["os_days_use"] = dm["os_days"].fillna(dm["days_to_last_follow_up"]).fillna(dm["days_to_death"]).astype(float)
    dm["sex_male"] = dm["sex"].astype(str).str.lower().eq("male").astype(int)

    age1 = dm.loc[dm["v17_dark_cluster"] == "DM1", "age"].astype(float)
    age2 = dm.loc[dm["v17_dark_cluster"] == "DM2", "age"].astype(float)
    age_welch = stats.ttest_ind(age1, age2, equal_var=False, nan_policy="omit")
    age_mw = stats.mannwhitneyu(age1, age2, alternative="two-sided")

    s1 = int(dm.loc[dm["v17_dark_cluster"] == "DM1", "stage_high"].sum())
    s2 = int((dm["v17_dark_cluster"] == "DM1").sum() - s1)
    s3 = int(dm.loc[dm["v17_dark_cluster"] == "DM2", "stage_high"].sum())
    s4 = int((dm["v17_dark_cluster"] == "DM2").sum() - s3)
    stage_p = fisher_2x2(s1, s2, s3, s4)
    rho, trend_p = stats.spearmanr(dm["cluster_bin"], dm["stage_ord"], nan_policy="omit")

    hist_tab = pd.crosstab(dm["v17_dark_cluster"], dm["histology_subtype"])
    hist_p = np.nan
    if hist_tab.shape == (2, 2):
        hist_p = fisher_2x2(int(hist_tab.iloc[0, 0]), int(hist_tab.iloc[0, 1]), int(hist_tab.iloc[1, 0]), int(hist_tab.iloc[1, 1]))
    else:
        try:
            hist_p = float(stats.chi2_contingency(hist_tab)[1])
        except Exception:
            pass

    tsize_p = float(stats.ttest_ind(
        dm.loc[dm["v17_dark_cluster"] == "DM1", "tumor_size_mm"].astype(float),
        dm.loc[dm["v17_dark_cluster"] == "DM2", "tumor_size_mm"].astype(float),
        equal_var=False,
        nan_policy="omit",
    ).pvalue)

    summary_rows = [
        {"analysis": "age_welch_t", "estimate": float(age2.mean() - age1.mean()), "pvalue": float(age_welch.pvalue), "detail": "DM2-DM1 mean age"},
        {"analysis": "age_mannwhitney", "estimate": float(age_mw.statistic), "pvalue": float(age_mw.pvalue), "detail": "rank-sum"},
        {"analysis": "stage_high_fisher", "estimate": float((s3 / max(1, s3 + s4)) - (s1 / max(1, s1 + s2))), "pvalue": stage_p, "detail": "Stage III+ fraction DM2-DM1"},
        {"analysis": "stage_trend_spearman", "estimate": float(rho), "pvalue": float(trend_p), "detail": "cluster vs stage ordinal"},
        {"analysis": "histology_enrichment", "estimate": float(hist_tab.to_numpy().max()), "pvalue": hist_p, "detail": "histology subtype by cluster"},
        {"analysis": "tumor_size_mm", "estimate": float(dm.loc[dm['v17_dark_cluster']=='DM2','tumor_size_mm'].mean() - dm.loc[dm['v17_dark_cluster']=='DM1','tumor_size_mm'].mean()), "pvalue": tsize_p, "detail": "DM2-DM1 tumor size"},
    ]
    pd.DataFrame(summary_rows).to_csv(TAB / "clinical_stratification_summary.tsv", sep="\t", index=False)

    age_adj = []
    cph = None
    cox_df = dm[["os_days_use", "os_event_use", "cluster_bin", "age", "sex_male", "stage_high"]].dropna().copy()
    if len(cox_df) >= 30 and cox_df["os_event_use"].sum() >= 5:
        cph = CoxPHFitter()
        cph.fit(cox_df, duration_col="os_days_use", event_col="os_event_use")
        cox_out = cph.summary.reset_index().rename(columns={"index": "covariate"})
        cox_out.to_csv(TAB / "cox_regression_results.tsv", sep="\t", index=False)
        if "cluster_bin" in cph.summary.index:
            s = cph.summary.loc["cluster_bin"]
            age_adj.append({"model": "Cox_OS", "coef": float(s["coef"]), "exp_coef": float(s["exp(coef)"]), "pvalue": float(s["p"])})
    else:
        pd.DataFrame([{"covariate": "NA", "note": "insufficient OS events"}]).to_csv(TAB / "cox_regression_results.tsv", sep="\t", index=False)

    reg_df = dm[["tds16_score_v17", "cluster_bin", "age"]].dropna().copy()
    if len(reg_df) >= 30:
        slope, intercept, r, p, se = stats.linregress(reg_df["age"], reg_df["tds16_score_v17"])
        resid = reg_df["tds16_score_v17"] - (intercept + slope * reg_df["age"])
        p_resid = stats.ttest_ind(resid[reg_df["cluster_bin"] == 0], resid[reg_df["cluster_bin"] == 1], equal_var=False).pvalue
        age_adj.append({"model": "TDS_residual_after_age", "coef": float(resid[reg_df["cluster_bin"] == 1].mean() - resid[reg_df["cluster_bin"] == 0].mean()), "exp_coef": np.nan, "pvalue": float(p_resid)})
    pd.DataFrame(age_adj).to_csv(TAB / "age_adjusted_cluster_effect.tsv", sep="\t", index=False)

    km = KaplanMeierFitter()
    fig = go.Figure()
    for cluster, color in [("DM1", "#22c55e"), ("DM2", "#f97316")]:
        sub = dm[dm["v17_dark_cluster"] == cluster]
        if sub["os_days_use"].notna().sum() == 0:
            continue
        km.fit(sub["os_days_use"], sub["os_event_use"], label=cluster)
        sf = km.survival_function_.reset_index()
        fig.add_trace(go.Scatter(x=sf.iloc[:, 0], y=sf.iloc[:, 1], mode="lines", name=cluster, line=dict(color=color)))
    fig.update_layout(title="DM1/DM2 overall survival", xaxis_title="Days", yaxis_title="Survival")
    fig.write_html(FIG / "km_curves_panel.html", include_plotlyjs="cdn")

    if cph is not None:
        csum = cph.summary.reset_index().rename(columns={"index": "covariate"})
        csum["lower"] = csum["exp(coef) lower 95%"]
        csum["upper"] = csum["exp(coef) upper 95%"]
        px.scatter(csum, x="exp(coef)", y="covariate", error_x=csum["upper"] - csum["exp(coef)"], error_x_minus=csum["exp(coef)"] - csum["lower"], color="p").write_html(
            FIG / "cox_forest_plot.html", include_plotlyjs="cdn"
        )
    else:
        go.Figure().update_layout(title="Cox forest plot unavailable: insufficient events").write_html(FIG / "cox_forest_plot.html", include_plotlyjs="cdn")

    px.violin(dm, x="v17_dark_cluster", y="age", box=True, points="all", color="v17_dark_cluster", title="Age by DM cluster").write_html(
        FIG / "age_violin_with_pvalues.html", include_plotlyjs="cdn"
    )
    stage_counts = dm.groupby(["v17_dark_cluster", "ajcc_stage_group"]).size().reset_index(name="n")
    px.bar(stage_counts, x="ajcc_stage_group", y="n", color="v17_dark_cluster", barmode="group", title="AJCC stage mosaic proxy").write_html(
        FIG / "stage_mosaic.html", include_plotlyjs="cdn"
    )

    summary = {
        "dm1_n": int((dm["v17_dark_cluster"] == "DM1").sum()),
        "dm2_n": int((dm["v17_dark_cluster"] == "DM2").sum()),
        "age_p_welch": float(age_welch.pvalue),
        "stage_high_p": float(stage_p) if pd.notna(stage_p) else None,
        "os_events": int(dm["os_event_use"].sum()),
    }
    write_json(TAB / "clinical_summary.json", summary)
    log_line(LOG, f"v17p2 layer2 done {summary}")


if __name__ == "__main__":
    main()

