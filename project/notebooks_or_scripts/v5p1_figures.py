#!/usr/bin/env python3
"""v5.1 Phase 5 — Figures."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px

sys.path.insert(0, "/opt/thyroid-dash/project/notebooks_or_scripts")
from v5p1_common import RESULTS_V5, FIGS_V5, LOGS, log_line

LOG = LOGS / "v5p1_figures.log"
AMBER = "#F5A623"

TEMPL = "plotly_dark"


def load_results():
    d = pd.read_csv(RESULTS_V5 / "v5p1_dial_all_cancers.tsv", sep="\t")
    lg_path = RESULTS_V5 / "v5p1_linearity_gap.tsv"
    lg = pd.read_csv(lg_path, sep="\t") if lg_path.exists() else pd.DataFrame()
    coh_path = RESULTS_V5 / "v5p1_cohort_availability.tsv"
    coh = pd.read_csv(coh_path, sep="\t") if coh_path.exists() else pd.DataFrame()
    return d, lg, coh


def fig1_dial_heatmap(d):
    piv = d.pivot_table(index="cancer", columns="classifier", values="dial", aggfunc="mean")
    fig = go.Figure(data=go.Heatmap(
        z=piv.values, x=piv.columns, y=piv.index,
        colorscale=[[0.0, "#102030"], [0.5, "#333333"], [1.0, AMBER]],
        colorbar=dict(title="DIAL"),
        zmin=0, zmax=max(0.5, float(np.nanmax(piv.values))),
        text=np.round(piv.values, 3), texttemplate="%{text}"))
    fig.update_layout(title="v5.1 DIAL: Cancer x Classifier (REAL data)",
                      template=TEMPL, height=420)
    fig.write_html(FIGS_V5 / "v5p1_fig1_dial_heatmap.html", include_plotlyjs="cdn")


def fig2_auc_pre_vs_post(d):
    fig = go.Figure()
    for cancer, g in d.groupby("cancer"):
        fig.add_trace(go.Scatter(x=g["auc_pre"], y=g["auc_post"], mode="markers+text",
                                 text=g["classifier"], textposition="top center",
                                 name=cancer, marker=dict(size=12)))
    xs = np.linspace(0, 1, 50)
    fig.add_trace(go.Scatter(x=xs, y=xs, mode="lines", name="y=x",
                             line=dict(dash="dash", color="grey")))
    fig.add_trace(go.Scatter(x=xs, y=1-xs, mode="lines", name="y=1-x (flip)",
                             line=dict(dash="dot", color=AMBER)))
    fig.update_layout(title="v5.1 AUC pre vs post-ComBat (REAL data)",
                      xaxis_title="AUC pre", yaxis_title="AUC post",
                      template=TEMPL, height=520)
    fig.write_html(FIGS_V5 / "v5p1_fig2_auc_pre_vs_post_scatter.html", include_plotlyjs="cdn")


def fig3_interpretation(d):
    piv = d.groupby(["cancer", "interpretation"]).size().reset_index(name="n")
    total = piv.groupby("cancer")["n"].transform("sum")
    piv["frac"] = piv["n"] / total
    fig = px.bar(piv, x="cancer", y="n", color="interpretation", barmode="stack",
                 title="v5.1 DIAL interpretation breakdown per cancer",
                 color_discrete_map={"batch_entangled": AMBER,
                                     "true_biology": "#4CC38A",
                                     "no_signal": "#888888",
                                     "ambiguous": "#6EC0FF"})
    fig.update_layout(template=TEMPL, height=420)
    fig.write_html(FIGS_V5 / "v5p1_fig3_interpretation_breakdown.html", include_plotlyjs="cdn")


def fig4_flip_preservation(d):
    fig = go.Figure()
    for cancer, g in d.groupby("cancer"):
        fig.add_trace(go.Scatter(x=g["auc_flip_pre"], y=g["auc_flip_post"],
                                 mode="markers+text", text=g["classifier"],
                                 textposition="top center", name=cancer,
                                 marker=dict(size=12)))
    xs = np.linspace(0.5, 1.0, 50)
    fig.add_trace(go.Scatter(x=xs, y=xs, mode="lines", name="y=x",
                             line=dict(dash="dash", color="grey")))
    try:
        x = d["auc_flip_pre"].values
        y = d["auc_flip_post"].values
        ok = ~(np.isnan(x) | np.isnan(y))
        r = np.corrcoef(x[ok], y[ok])[0, 1] if ok.sum() > 2 else float("nan")
    except Exception:
        r = float("nan")
    fig.update_layout(title=f"v5.1 flip_pre vs flip_post (r={r:.3f})",
                      xaxis_title="flip_pre", yaxis_title="flip_post",
                      template=TEMPL, height=460)
    fig.write_html(FIGS_V5 / "v5p1_fig4_flip_preservation.html", include_plotlyjs="cdn")


def fig5_identifiability_vs_dial(d):
    fig = go.Figure()
    for cancer, g in d.groupby("cancer"):
        fig.add_trace(go.Scatter(x=g["batch_identifiability_post"], y=g["dial"],
                                 mode="markers+text", text=g["classifier"],
                                 textposition="top center", name=cancer,
                                 marker=dict(size=12)))
    fig.add_hline(y=0.3, line_dash="dash", line_color=AMBER,
                  annotation_text="DIAL=0.3 (batch-entangled)",
                  annotation_position="top left")
    fig.add_vline(x=0.7, line_dash="dash", line_color="grey",
                  annotation_text="ident=0.7", annotation_position="top right")
    fig.update_layout(title="v5.1 Batch identifiability vs DIAL (quadrants)",
                      xaxis_title="Batch identifiability (post)", yaxis_title="DIAL",
                      template=TEMPL, height=480)
    fig.write_html(FIGS_V5 / "v5p1_fig5_identifiability_vs_dial.html", include_plotlyjs="cdn")


def fig6_linearity_gap(lg):
    if lg.empty:
        return
    fig = go.Figure()
    fig.add_trace(go.Bar(x=lg["cancer"], y=lg["dial_linear_mean"],
                         name="linear mean (LogReg)", marker_color=AMBER))
    fig.add_trace(go.Bar(x=lg["cancer"], y=lg["dial_nonlinear_mean"],
                         name="nonlinear mean (RF/GB/XGB)", marker_color="#4CC38A"))
    fig.add_trace(go.Scatter(x=lg["cancer"], y=lg["dial_linearity_gap"],
                             mode="markers+text+lines",
                             text=[f"gap={x:.2f} ({l})" for x, l in zip(lg["dial_linearity_gap"], lg["label"])],
                             textposition="top center",
                             name="gap (linear - nonlinear)",
                             marker=dict(color="#FF6F61", size=14)))
    fig.update_layout(title="v5.1 NEW: Linearity gap analysis",
                      xaxis_title="cancer", yaxis_title="DIAL",
                      template=TEMPL, height=520, barmode="group")
    fig.write_html(FIGS_V5 / "v5p1_fig6_linearity_gap.html", include_plotlyjs="cdn")


def fig7_summary_table(d, lg, coh):
    # one-row-per-cancer summary
    rows = []
    for cancer, g in d.groupby("cancer"):
        row = {
            "cancer": cancer,
            "n_samples": int(g["n_samples"].iloc[0]),
            "n_genes": int(g["n_genes"].iloc[0]),
            "dial_mean": round(g["dial"].mean(), 3),
            "dial_linear": round(g[g["classifier"].isin(["LogReg_l2","LogReg_elasticnet"])]["dial"].mean(), 3),
            "dial_nonlinear": round(g[g["classifier"].isin(["RandomForest","GradientBoosting","XGBoost","HistGB"])]["dial"].mean(), 3),
            "auc_post_mean": round(g["auc_post"].mean(), 3),
            "batch_ident_post": round(g["batch_identifiability_post"].mean(), 3),
        }
        if not lg.empty:
            m = lg[lg["cancer"] == cancer]
            if len(m):
                row["linearity_label"] = m.iloc[0]["label"]
                row["gap"] = round(m.iloc[0]["dial_linearity_gap"], 3)
        rows.append(row)
    df = pd.DataFrame(rows)
    fig = go.Figure(data=[go.Table(
        header=dict(values=list(df.columns), fill_color="#222", font_color=AMBER,
                    align="left"),
        cells=dict(values=[df[c] for c in df.columns], fill_color="#111",
                   font_color="white", align="left"))])
    fig.update_layout(title="v5.1 Cross-cancer DIAL summary (REAL data)",
                      template=TEMPL, height=60 + 40 * len(df))
    fig.write_html(FIGS_V5 / "v5p1_fig7_summary_table.html", include_plotlyjs="cdn")


def main():
    log_line(LOG, "=== v5p1 Phase 5 figures START ===")
    d, lg, coh = load_results()
    if d.empty:
        log_line(LOG, "no DIAL results — skipping")
        return
    fig1_dial_heatmap(d)
    fig2_auc_pre_vs_post(d)
    fig3_interpretation(d)
    fig4_flip_preservation(d)
    fig5_identifiability_vs_dial(d)
    fig6_linearity_gap(lg)
    fig7_summary_table(d, lg, coh)
    log_line(LOG, "=== v5p1 Phase 5 DONE ===")


if __name__ == "__main__":
    main()
