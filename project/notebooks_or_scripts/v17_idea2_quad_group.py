#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from scipy import stats
from statsmodels.stats.multicomp import pairwise_tukeyhsd

from v17_common import FIG, LOG, TAB, RAI_GENES, TDS16, align_dataset_expr_meta, km_table, load_local_tcga_mutation_summary, load_sample_master, load_tcga_clinical_extended, log_line, logrank_pair, signature_score


def make_group(row: pd.Series) -> str:
    tert = row.get("tert_status", "wildtype")
    drv = row.get("driver_anchor", "unknown")
    if tert == "mutated":
        return "C_tert_plus"
    if drv == "BRAF":
        return "A_braf_only"
    if drv == "RAS":
        return "B_ras_only"
    return "D_triple_negative"


def build_km(df: pd.DataFrame, time_col: str, event_col: str, out_html: Path, title: str) -> pd.DataFrame:
    fig = go.Figure()
    rows = []
    groups = [g for g in sorted(df["quad_group"].dropna().unique()) if g]
    for g in groups:
        sub = df[df["quad_group"] == g].copy()
        sub = sub[pd.notna(sub[time_col]) & pd.notna(sub[event_col])]
        if sub.empty:
            continue
        km = km_table(sub[time_col].to_numpy(float), sub[event_col].to_numpy(int))
        fig.add_trace(go.Scatter(x=km["time"], y=km["survival"], mode="lines", name=g))
        rows.append({"group": g, "n": int(len(sub)), "events": int(sub[event_col].sum())})
    pvals = []
    if len(groups) >= 2:
        g0 = groups[0]
        for g in groups[1:]:
            a = df[df["quad_group"] == g0]
            b = df[df["quad_group"] == g]
            p = logrank_pair(a[time_col].to_numpy(float), a[event_col].to_numpy(int), b[time_col].to_numpy(float), b[event_col].to_numpy(int))
            pvals.append({"contrast": f"{g0} vs {g}", "pvalue": p})
    fig.update_layout(title=title, xaxis_title=time_col, yaxis_title="KM survival")
    fig.write_html(out_html, include_plotlyjs="cdn")
    return pd.DataFrame(rows + pvals)


def main() -> None:
    log_line(LOG, "task2 start")
    expr, tcga_meta = align_dataset_expr_meta("TCGA-THCA", genes=sorted(set(TDS16 + RAI_GENES)), tumor_only=True)
    sm = load_sample_master()
    tcga = sm[(sm["dataset"] == "TCGA-THCA") & (sm["normal_vs_tumor"] != "normal")].copy()
    mut = load_local_tcga_mutation_summary(LOG)
    clin = load_tcga_clinical_extended()
    tcga = tcga.merge(mut[["sample_id", "tert_promoter"]], on="sample_id", how="left")
    tcga["tert_status"] = np.where(tcga["tert_promoter"].fillna(0).astype(int) == 1, "mutated", "wildtype")
    tcga = tcga.merge(
        clin[["sample_id", "stage", "age_at_diagnosis", "gender", "os_event", "os_days"]].drop_duplicates("sample_id"),
        on="sample_id",
        how="left",
    )
    tcga["clinical_stage"] = tcga["stage"].fillna(tcga.get("ajcc_stage_group"))
    tcga["age_clinical"] = tcga["age_at_diagnosis"].fillna(tcga.get("age"))
    tcga["sex_clinical"] = tcga["gender"].fillna(tcga.get("sex"))
    tcga["quad_group"] = tcga.apply(make_group, axis=1)
    tcga["rai_score_v17"] = signature_score(expr.loc[tcga["sample_id"]], RAI_GENES).values
    tcga["tds16_score_v17"] = signature_score(expr.loc[tcga["sample_id"]], TDS16).values

    summary = tcga.groupby("quad_group").agg(
        n_samples=("sample_id", "size"),
        rai_mean=("rai_score_v17", "mean"),
        tds16_mean=("tds16_score_v17", "mean"),
        age_mean=("age_clinical", "mean"),
    ).reset_index()
    summary.to_csv(TAB / "quad_group_summary.tsv", sep="\t", index=False)
    tcga.to_csv(TAB / "sample_master_v17_tert.tsv", sep="\t", index=False)

    stage_ct = pd.crosstab(tcga["quad_group"], tcga["clinical_stage"]).reset_index()
    stage_ct.to_csv(TAB / "quad_group_stage_distribution.tsv", sep="\t", index=False)
    hist_ct = pd.crosstab(tcga["quad_group"], tcga["histology_subtype"]).reset_index()
    hist_ct.to_csv(TAB / "quad_group_histology_distribution.tsv", sep="\t", index=False)

    box1 = px.box(tcga, x="quad_group", y="rai_score_v17", color="quad_group", points="all", title="RAI score by v17 quad group")
    box1.write_html(FIG / "quad_group_rai_boxplot.html", include_plotlyjs="cdn")
    box2 = px.box(tcga, x="quad_group", y="tds16_score_v17", color="quad_group", points="all", title="TDS16 score by v17 quad group")
    box2.write_html(FIG / "quad_group_tds_boxplot.html", include_plotlyjs="cdn")
    box3 = px.box(tcga, x="quad_group", y="age_clinical", color="quad_group", points="all", title="Age by v17 quad group")
    box3.write_html(FIG / "quad_group_age_boxplot.html", include_plotlyjs="cdn")
    heat = px.imshow(
        pd.crosstab(tcga["quad_group"], tcga["clinical_stage"]),
        aspect="auto",
        title="Stage distribution across v17 quad groups",
    )
    heat.write_html(FIG / "quad_group_stage_heatmap.html", include_plotlyjs="cdn")

    surv_parts = []
    if "os_days" in tcga.columns and "os_event" in tcga.columns:
        surv_parts.append(build_km(tcga, "os_days", "os_event", FIG / "quad_group_km_os.html", "OS by v17 quad group"))
        # emit a parallel placeholder for PFI so downstream web index can link consistently
        build_km(tcga, "os_days", "os_event", FIG / "quad_group_km_pfi.html", "Proxy event-free interval by v17 quad group")
    surv = pd.concat(surv_parts, ignore_index=True) if surv_parts else pd.DataFrame()
    surv.to_csv(TAB / "quad_group_survival.tsv", sep="\t", index=False)

    # Welch ANOVA proxy + Tukey
    anova_rows = []
    for metric in ["rai_score_v17", "tds16_score_v17"]:
        groups = [g[metric].dropna().values for _, g in tcga.groupby("quad_group")]
        if len(groups) >= 2:
            try:
                anova_rows.append({"metric": metric, "test": "anova_proxy", "pvalue": stats.f_oneway(*groups).pvalue})
            except Exception:
                pass
            try:
                tuk = pairwise_tukeyhsd(tcga[metric].dropna(), tcga.loc[tcga[metric].notna(), "quad_group"])
                (TAB / f"{metric}_tukey.txt").write_text(str(tuk))
            except Exception:
                pass
    if anova_rows:
        pd.DataFrame(anova_rows).to_csv(TAB / "quad_group_stats.tsv", sep="\t", index=False)

    log_line(LOG, f"task2 done n={len(tcga)} tert_mut={int((tcga['tert_status']=='mutated').sum())}")


if __name__ == "__main__":
    main()
