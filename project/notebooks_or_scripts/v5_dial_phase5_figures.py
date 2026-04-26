#!/usr/bin/env python3
"""v5 DIAL Phase 5 — Build 6 interactive Plotly figures.

All figures use template="plotly_dark" with transparent paper/plot bg
and amber accent #F5A623.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px

sys.path.insert(0, str(Path(__file__).parent))
from v5_dial_common import FIGS_V5, RESULTS_V5, LOGS, log_line

LOGFILE = LOGS / "v5_dial_run.log"

ACCENT = "#F5A623"
ACCENT_LIGHT = "#FFD78A"
BG = "rgba(0,0,0,0)"


def _style(fig: go.Figure, title: str, height: int = 520) -> go.Figure:
    fig.update_layout(
        template="plotly_dark",
        paper_bgcolor=BG,
        plot_bgcolor=BG,
        title=dict(text=title, x=0.02, y=0.97,
                   font=dict(size=18, color="#e6edf3", family="Inter, sans-serif")),
        margin=dict(l=60, r=30, t=70, b=60),
        height=height,
        font=dict(family="Inter, JetBrains Mono, monospace",
                  size=12, color="#c9d1d9"),
        hoverlabel=dict(bgcolor="#101418", font_size=12, font_family="Inter"),
    )
    return fig


def fig1_heatmap(df: pd.DataFrame) -> go.Figure:
    pivot = df.pivot_table(index="cancer", columns="classifier", values="dial",
                           aggfunc="mean")
    order_c = ["THCA", "SKCM", "LGG", "LUAD", "COAD"]
    pivot = pivot.reindex([c for c in order_c if c in pivot.index])
    fig = go.Figure(data=go.Heatmap(
        z=pivot.values, x=pivot.columns.tolist(), y=pivot.index.tolist(),
        colorscale=[[0, "#1a1e24"], [0.25, "#3a2a18"], [0.6, "#a66820"], [1, ACCENT]],
        zmin=0, zmax=0.5,
        colorbar=dict(title=dict(text="DIAL", side="right"),
                      tickfont=dict(color="#c9d1d9")),
        hovertemplate="<b>%{y}</b> × %{x}<br>DIAL=%{z:.3f}<extra></extra>",
        showscale=True,
    ))
    # add DIAL value annotations
    for i, r in enumerate(pivot.index):
        for j, c in enumerate(pivot.columns):
            v = pivot.values[i, j]
            if not np.isnan(v):
                col = "#0b0e12" if v > 0.25 else "#e6edf3"
                fig.add_annotation(x=c, y=r, text=f"{v:.2f}",
                                   showarrow=False, font=dict(color=col, size=12))
    return _style(fig, "Figure 1 — DIAL heatmap across 5 cancers × 5 classifiers", 460)


def fig2_pre_vs_post(df: pd.DataFrame) -> go.Figure:
    fig = go.Figure()
    for cancer in sorted(df["cancer"].unique()):
        sub = df[df["cancer"] == cancer]
        fig.add_trace(go.Scatter(
            x=sub["auc_pre"], y=sub["auc_post"],
            mode="markers+text",
            name=cancer,
            text=sub["classifier"],
            textposition="top center",
            textfont=dict(size=9, color="#c9d1d9"),
            marker=dict(size=11, line=dict(width=1, color="#0b0e12")),
            hovertemplate=f"<b>{cancer}</b><br>%{{text}}<br>AUC_pre=%{{x:.3f}}<br>AUC_post=%{{y:.3f}}<extra></extra>",
        ))
    # diagonals
    fig.add_shape(type="line", x0=0, y0=0, x1=1, y1=1,
                  line=dict(color="#555", dash="dot", width=1))
    fig.add_shape(type="line", x0=0, y0=1, x1=1, y1=0,
                  line=dict(color=ACCENT, dash="dash", width=1.5))
    fig.add_annotation(x=0.03, y=0.97, showarrow=False,
                       text="<i>y = 1-x:</i> perfect label-flip", font=dict(color=ACCENT, size=11))
    fig.update_xaxes(title="AUC pre-correction", range=[0, 1.02],
                     gridcolor="#22262d", zeroline=False)
    fig.update_yaxes(title="AUC post-correction", range=[0, 1.02],
                     gridcolor="#22262d", zeroline=False)
    return _style(fig, "Figure 2 — AUC pre vs AUC post (anti-diagonal = label flip)")


def fig3_interpretation_breakdown(df: pd.DataFrame) -> go.Figure:
    cats = ["batch_entangled", "ambiguous", "true_biology", "no_signal", "error"]
    palette = {
        "batch_entangled": ACCENT,
        "ambiguous": "#8c7853",
        "true_biology": "#5fb878",
        "no_signal": "#4a5263",
        "error": "#c24c4c",
    }
    cancers = ["THCA", "SKCM", "LGG", "LUAD", "COAD"]
    counts = {c: {k: 0 for k in cats} for c in cancers}
    for _, r in df.iterrows():
        c, k = r["cancer"], r["interpretation"]
        if c in counts and k in counts[c]:
            counts[c][k] += 1
    fig = go.Figure()
    for cat in cats:
        fig.add_trace(go.Bar(
            x=cancers, y=[counts[c][cat] for c in cancers],
            name=cat, marker_color=palette[cat],
            hovertemplate=f"<b>%{{x}}</b><br>{cat}: %{{y}}<extra></extra>",
        ))
    fig.update_layout(barmode="stack")
    fig.update_xaxes(gridcolor="#22262d")
    fig.update_yaxes(title="# classifier × cancer pairs", gridcolor="#22262d")
    return _style(fig, "Figure 3 — Interpretation breakdown per cancer")


def fig4_flip_preservation(df: pd.DataFrame) -> go.Figure:
    fig = go.Figure()
    for cancer in sorted(df["cancer"].unique()):
        sub = df[df["cancer"] == cancer]
        fig.add_trace(go.Scatter(
            x=sub["auc_flip_pre"], y=sub["auc_flip_post"],
            mode="markers+text",
            name=cancer,
            text=sub["classifier"],
            textposition="top center",
            textfont=dict(size=9, color="#c9d1d9"),
            marker=dict(size=11, line=dict(width=1, color="#0b0e12")),
        ))
    fig.add_shape(type="line", x0=0.4, y0=0.4, x1=1.02, y1=1.02,
                  line=dict(color=ACCENT, dash="dash", width=1.5))
    # correlation annotation
    valid = df.dropna(subset=["auc_flip_pre", "auc_flip_post"])
    if len(valid) > 3:
        r = float(np.corrcoef(valid["auc_flip_pre"], valid["auc_flip_post"])[0, 1])
        fig.add_annotation(x=0.45, y=1.0, showarrow=False,
                           text=f"<b>r = {r:.2f}</b>", font=dict(color=ACCENT, size=14))
    fig.update_xaxes(title="auc_flip pre", range=[0.4, 1.02], gridcolor="#22262d")
    fig.update_yaxes(title="auc_flip post", range=[0.4, 1.02], gridcolor="#22262d")
    return _style(fig, "Figure 4 — Direction-invariant leakage is preserved under ComBat")


def fig5_quadrant(df: pd.DataFrame) -> go.Figure:
    fig = go.Figure()
    for cancer in sorted(df["cancer"].unique()):
        sub = df[df["cancer"] == cancer]
        fig.add_trace(go.Scatter(
            x=sub["dial"], y=sub["batch_identifiability_post"],
            mode="markers+text",
            name=cancer, text=sub["classifier"],
            textposition="top center", textfont=dict(size=9, color="#c9d1d9"),
            marker=dict(size=12, line=dict(width=1, color="#0b0e12")),
            hovertemplate=f"<b>{cancer}</b><br>%{{text}}<br>DIAL=%{{x:.3f}}<br>Ident=%{{y:.3f}}<extra></extra>",
        ))
    # threshold lines
    fig.add_shape(type="line", x0=0.3, y0=0, x1=0.3, y1=1,
                  line=dict(color="#555", dash="dot", width=1))
    fig.add_shape(type="line", x0=0, y0=0.7, x1=0.55, y1=0.7,
                  line=dict(color="#555", dash="dot", width=1))
    # quadrant labels
    fig.add_annotation(x=0.45, y=0.35, showarrow=False,
                       text="<b>BATCH-ENTANGLED</b><br><i>(high DIAL, low ident)</i>",
                       font=dict(color=ACCENT, size=12))
    fig.add_annotation(x=0.05, y=0.95, showarrow=False,
                       text="true biology regime<br>(low DIAL, high ident)",
                       font=dict(color="#5fb878", size=11))
    fig.add_annotation(x=0.45, y=0.92, showarrow=False,
                       text="residual batch info<br>(rare)",
                       font=dict(color="#8c7853", size=11))
    fig.add_annotation(x=0.05, y=0.35, showarrow=False,
                       text="no signal<br>(low DIAL, low ident)",
                       font=dict(color="#4a5263", size=11))
    fig.update_xaxes(title="DIAL", range=[-0.02, 0.55], gridcolor="#22262d")
    fig.update_yaxes(title="post-correction batch identifiability (macro-OvR AUC)",
                     range=[0, 1.02], gridcolor="#22262d")
    return _style(fig, "Figure 5 — Quadrant diagnostic: DIAL vs batch identifiability")


def fig6_summary_table(df: pd.DataFrame) -> go.Figure:
    cols = ["cancer", "classifier", "auc_pre", "auc_post", "auc_flip_post",
            "dial", "batch_identifiability_post", "interpretation"]
    df2 = df[cols].copy()
    for col in ["auc_pre", "auc_post", "auc_flip_post", "dial", "batch_identifiability_post"]:
        df2[col] = df2[col].apply(lambda v: f"{v:.3f}" if pd.notna(v) else "—")

    # color rows by interpretation
    color_map = {"batch_entangled": "#3a2a18", "ambiguous": "#2a261c",
                 "true_biology": "#1a2e1f", "no_signal": "#1a1d24",
                 "error": "#2c1a1a"}
    row_colors = [color_map.get(c, "#101418") for c in df2["interpretation"]]
    n_rows = len(df2)

    fig = go.Figure(data=[go.Table(
        header=dict(
            values=[f"<b>{c}</b>" for c in df2.columns],
            fill_color="#1d2127", font=dict(color="#F5A623", size=12, family="Inter"),
            align="left", line_color="#2a2f36", height=36,
        ),
        cells=dict(
            values=[df2[c].tolist() for c in df2.columns],
            fill_color=[row_colors] * len(df2.columns),
            font=dict(color="#e6edf3", size=11, family="JetBrains Mono"),
            align="left", line_color="#2a2f36", height=28,
        ),
    )])
    return _style(fig, "Figure 6 — Summary table (sortable in place)", height=min(720, 120 + n_rows * 30))


def main():
    log_line(LOGFILE, "PHASE5 start — build 6 figures")
    df = pd.read_csv(RESULTS_V5 / "v5_dial_all_cancers.tsv", sep="\t")
    log_line(LOGFILE, f"PHASE5 loaded {len(df)} rows")

    figs = {
        "v5_dial_fig1_heatmap.html":         fig1_heatmap(df),
        "v5_dial_fig2_pre_vs_post_scatter.html": fig2_pre_vs_post(df),
        "v5_dial_fig3_interpretation_breakdown.html": fig3_interpretation_breakdown(df),
        "v5_dial_fig4_flip_preservation_scatter.html": fig4_flip_preservation(df),
        "v5_dial_fig5_identifiability_vs_dial.html": fig5_quadrant(df),
        "v5_dial_fig6_summary_table.html":   fig6_summary_table(df),
    }
    for name, fig in figs.items():
        fig.write_html(FIGS_V5 / name, include_plotlyjs="cdn", full_html=True)
        log_line(LOGFILE, f"PHASE5 wrote {name}")
    log_line(LOGFILE, "PHASE5 done")


if __name__ == "__main__":
    main()
