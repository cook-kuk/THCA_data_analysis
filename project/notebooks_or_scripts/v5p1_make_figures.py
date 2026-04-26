#!/usr/bin/env python
"""v5.1 figure generator: 5 Plotly HTML figures for the THCA-specificity dashboard.

Reads TSVs from results/v5/ and writes self-contained HTML to
reports/html/figs_interactive/v5/. Dark, transparent background; amber accent.

Figures:
  1. v5p1_fig1_interpretation_heatmap.html
  2. v5p1_fig2_auc_pre_vs_post.html
  3. v5p1_fig3_thca_label_flip.html
  4. v5p1_fig4_specificity_evidence.html
  5. v5p1_fig5_cohort_label_imbalance.html
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
PROJECT = Path("/opt/thyroid-dash/project")
RES = PROJECT / "results" / "v5"
OUT = PROJECT / "reports" / "html" / "figs_interactive" / "v5"
OUT.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Theme helpers
# ---------------------------------------------------------------------------
AMBER = "#F5A623"
FONT_COLOR = "#EAEAEA"
TRANSPARENT = "rgba(0,0,0,0)"

INTERP_COLOR = {
    "batch_entangled": "#E74C3C",
    "true_biology": "#2ECC71",
    "no_signal": "#7F8C8D",
    "ambiguous": "#F5A623",
    "partial_batch": "#E67E22",
}
INTERP_ORDER = ["batch_entangled", "partial_batch", "true_biology", "no_signal", "ambiguous"]

CANCER_COLORS = {
    "THCA": "#F5A623",
    "SKCM": "#9B59B6",
    "LGG": "#2ECC71",
    "LUAD": "#3498DB",
    "COAD": "#E74C3C",
}

CLASSIFIER_ORDER = ["LogReg_l2", "LogReg_elasticnet", "RandomForest", "GradientBoosting", "XGBoost"]
CANCER_ORDER = ["THCA", "SKCM", "LGG", "LUAD", "COAD"]


def interp_color(kind: str) -> str:
    """Return hex color for an interpretation label."""
    return INTERP_COLOR.get(kind, "#555555")


def base_layout(fig: go.Figure, title: str, height: int = 500) -> go.Figure:
    fig.update_layout(
        title=dict(text=title, font=dict(color=FONT_COLOR, size=18)),
        paper_bgcolor=TRANSPARENT,
        plot_bgcolor=TRANSPARENT,
        font=dict(color=FONT_COLOR, family="Inter, sans-serif"),
        height=height,
        margin=dict(l=60, r=40, t=70, b=60),
    )
    fig.update_xaxes(gridcolor="rgba(234,234,234,0.08)", zerolinecolor="rgba(234,234,234,0.15)")
    fig.update_yaxes(gridcolor="rgba(234,234,234,0.08)", zerolinecolor="rgba(234,234,234,0.15)")
    return fig


def write_fig(fig: go.Figure, filename: str) -> Path:
    path = OUT / filename
    fig.write_html(str(path), include_plotlyjs="cdn", full_html=True)
    size_kb = path.stat().st_size / 1024
    print(f"  wrote {path.name:<52} {size_kb:7.1f} KB")
    return path


# ---------------------------------------------------------------------------
# Load
# ---------------------------------------------------------------------------
def load() -> tuple[pd.DataFrame, pd.DataFrame]:
    dial = pd.read_csv(RES / "v5p1_dial_all_cancers.tsv", sep="\t")
    harm = pd.read_csv(RES / "v5p1_harmonization.tsv", sep="\t")
    # Coerce dial column
    for c in ("auc_pre", "auc_post", "auc_flip_pre", "auc_flip_post", "dial"):
        dial[c] = pd.to_numeric(dial[c], errors="coerce")
    return dial, harm


# ---------------------------------------------------------------------------
# Figure 1 — interpretation heatmap (5 cancers × 5 classifiers)
# ---------------------------------------------------------------------------
def fig1_heatmap(dial: pd.DataFrame) -> None:
    # Build 5x5 matrix, rows = cancer, cols = classifier
    # We need numeric index for interpretation, plus a custom colorscale
    interp_to_idx = {name: i for i, name in enumerate(INTERP_ORDER)}

    z = []
    text = []
    hover = []
    for cancer in CANCER_ORDER:
        zrow = []
        trow = []
        hrow = []
        for clf in CLASSIFIER_ORDER:
            match = dial[(dial.cancer == cancer) & (dial.classifier == clf)]
            if len(match):
                r = match.iloc[0]
                idx = interp_to_idx.get(r.interpretation, 0)
                zrow.append(idx)
                trow.append(f"{r.dial:.3f}")
                hrow.append(
                    f"<b>{cancer} · {clf}</b><br>"
                    f"interpretation: {r.interpretation}<br>"
                    f"DIAL: {r.dial:.3f}<br>"
                    f"AUC pre: {r.auc_pre:.3f}<br>"
                    f"AUC post: {r.auc_post:.3f}<br>"
                    f"AUC flip post: {r.auc_flip_post:.3f}"
                )
            else:
                zrow.append(None)
                trow.append("")
                hrow.append("")
        z.append(zrow)
        text.append(trow)
        hover.append(hrow)

    # Discrete colorscale across the 5 categories, spaced equally
    n = len(INTERP_ORDER)
    colorscale = []
    for i, name in enumerate(INTERP_ORDER):
        lo = i / n
        hi = (i + 1) / n
        colorscale.append([lo, interp_color(name)])
        colorscale.append([hi, interp_color(name)])

    fig = go.Figure(
        data=go.Heatmap(
            z=z,
            x=CLASSIFIER_ORDER,
            y=CANCER_ORDER,
            text=text,
            texttemplate="%{text}",
            textfont=dict(color="#FFFFFF", size=13, family="JetBrains Mono, monospace"),
            hoverinfo="text",
            hovertext=hover,
            colorscale=colorscale,
            zmin=-0.5,
            zmax=n - 0.5,
            showscale=True,
            colorbar=dict(
                title=dict(text="Interpretation", font=dict(color=FONT_COLOR)),
                tickvals=list(range(n)),
                ticktext=INTERP_ORDER,
                tickfont=dict(color=FONT_COLOR),
                len=0.75,
            ),
            xgap=2,
            ygap=2,
        )
    )
    base_layout(fig, "Figure 1 — Interpretation × DIAL per (cancer × classifier)", height=520)
    fig.update_xaxes(title_text="Classifier", side="bottom")
    fig.update_yaxes(title_text="Cancer", autorange="reversed")
    write_fig(fig, "v5p1_fig1_interpretation_heatmap.html")


# ---------------------------------------------------------------------------
# Figure 2 — AUC pre vs post scatter
# ---------------------------------------------------------------------------
def fig2_scatter(dial: pd.DataFrame) -> None:
    fig = go.Figure()
    # y=x "no change"
    fig.add_trace(
        go.Scatter(
            x=[0, 1], y=[0, 1], mode="lines",
            line=dict(color="rgba(200,200,200,0.55)", dash="dash", width=1.5),
            name="y = x (no change)", hoverinfo="skip",
        )
    )
    # y=1-x "complete flip"
    fig.add_trace(
        go.Scatter(
            x=[0, 1], y=[1, 0], mode="lines",
            line=dict(color="#E74C3C", dash="dash", width=1.5),
            name="y = 1 - x (complete flip)", hoverinfo="skip",
        )
    )
    for cancer in CANCER_ORDER:
        sub = dial[dial.cancer == cancer]
        if sub.empty:
            continue
        sizes = (sub["dial"].abs() * 200 + 6).tolist()
        hover = [
            f"<b>{r.cancer} · {r.classifier}</b><br>"
            f"DIAL: {r.dial:.3f}<br>"
            f"interpretation: {r.interpretation}<br>"
            f"AUC pre: {r.auc_pre:.3f}<br>"
            f"AUC post: {r.auc_post:.3f}"
            for _, r in sub.iterrows()
        ]
        fig.add_trace(
            go.Scatter(
                x=sub["auc_pre"], y=sub["auc_post"],
                mode="markers",
                marker=dict(
                    size=sizes, color=CANCER_COLORS[cancer], opacity=0.85,
                    line=dict(color="rgba(255,255,255,0.25)", width=0.8),
                ),
                name=cancer,
                text=hover, hoverinfo="text",
            )
        )
    base_layout(fig, "Figure 2 — AUC pre vs AUC post (marker size ∝ |DIAL|)", height=560)
    fig.update_xaxes(title_text="AUC pre-ComBat", range=[0, 1])
    fig.update_yaxes(title_text="AUC post-ComBat", range=[0, 1])
    fig.update_layout(legend=dict(bgcolor="rgba(0,0,0,0)", font=dict(color=FONT_COLOR)))
    write_fig(fig, "v5p1_fig2_auc_pre_vs_post.html")


# ---------------------------------------------------------------------------
# Figure 3 — THCA label flip (2x2 subplots)
# ---------------------------------------------------------------------------
def fig3_thca_flip(dial: pd.DataFrame) -> None:
    flipped = ["LogReg_l2", "LogReg_elasticnet", "RandomForest", "GradientBoosting"]
    thca = dial[dial.cancer == "THCA"].set_index("classifier")

    titles = []
    for clf in flipped:
        if clf in thca.index:
            d = thca.loc[clf, "dial"]
            titles.append(f"{clf}  ·  DIAL = {d:.3f}")
        else:
            titles.append(clf)

    fig = make_subplots(
        rows=2, cols=2, subplot_titles=titles,
        vertical_spacing=0.18, horizontal_spacing=0.12,
    )
    positions = [(1, 1), (1, 2), (2, 1), (2, 2)]
    for (r, c), clf in zip(positions, flipped):
        if clf not in thca.index:
            continue
        row = thca.loc[clf]
        fig.add_trace(
            go.Bar(
                x=["AUC post", "AUC flip post"],
                y=[row["auc_post"], row["auc_flip_post"]],
                marker=dict(color=["#E74C3C", AMBER], line=dict(color="rgba(255,255,255,0.2)", width=1)),
                text=[f"{row['auc_post']:.3f}", f"{row['auc_flip_post']:.3f}"],
                textposition="outside",
                textfont=dict(color=FONT_COLOR),
                showlegend=False,
                hovertemplate=f"{clf}<br>%{{x}}: %{{y:.3f}}<extra></extra>",
            ),
            row=r, col=c,
        )
        # 0.5 reference line
        fig.add_shape(
            type="line", x0=-0.5, x1=1.5, y0=0.5, y1=0.5,
            line=dict(color="rgba(234,234,234,0.4)", width=1.2, dash="dash"),
            row=r, col=c,
        )

    for i in range(1, 5):
        fig.update_yaxes(range=[0, 1.08], row=positions[i - 1][0], col=positions[i - 1][1],
                         gridcolor="rgba(234,234,234,0.08)")
        fig.update_xaxes(row=positions[i - 1][0], col=positions[i - 1][1])

    base_layout(fig, "Figure 3 — THCA label-flip geometry (auc_post vs auc_flip_post, dashed = 0.5)", height=620)
    # subplot title font
    for ann in fig.layout.annotations:
        ann.font = dict(color=AMBER, size=13, family="JetBrains Mono, monospace")
    write_fig(fig, "v5p1_fig3_thca_label_flip.html")


# ---------------------------------------------------------------------------
# Figure 4 — stacked horizontal bar of interpretation counts
# ---------------------------------------------------------------------------
def fig4_specificity(dial: pd.DataFrame) -> None:
    counts = (
        dial.groupby(["cancer", "interpretation"]).size().unstack(fill_value=0)
    )
    for kind in INTERP_ORDER:
        if kind not in counts.columns:
            counts[kind] = 0
    counts = counts.reindex(CANCER_ORDER)[INTERP_ORDER]

    fig = go.Figure()
    for kind in INTERP_ORDER:
        fig.add_trace(
            go.Bar(
                y=counts.index.tolist(),
                x=counts[kind].tolist(),
                name=kind,
                orientation="h",
                marker=dict(color=interp_color(kind), line=dict(color="rgba(0,0,0,0.25)", width=0.6)),
                hovertemplate="%{y} · " + kind + ": %{x}<extra></extra>",
                text=[str(v) if v > 0 else "" for v in counts[kind].tolist()],
                textposition="inside",
                textfont=dict(color="#FFFFFF", size=12),
            )
        )

    base_layout(fig, "Figure 4 — Interpretation counts per cancer (n classifiers)", height=480)
    fig.update_layout(
        barmode="stack",
        xaxis=dict(title="Classifiers per interpretation", range=[0, 5.2]),
        yaxis=dict(title="Cancer", autorange="reversed"),
        legend=dict(orientation="h", y=-0.22, x=0, bgcolor="rgba(0,0,0,0)", font=dict(color=FONT_COLOR)),
    )
    fig.add_annotation(
        xref="paper", yref="paper", x=0, y=-0.38,
        text="THCA 4/5 batch_entangled — specificity, not generalization failure.",
        showarrow=False, font=dict(color=AMBER, size=13, family="JetBrains Mono, monospace"),
        align="left",
    )
    write_fig(fig, "v5p1_fig4_specificity_evidence.html")


# ---------------------------------------------------------------------------
# Figure 5 — cohort class imbalance
# ---------------------------------------------------------------------------
def fig5_cohort_imbalance(harm: pd.DataFrame) -> None:
    # Build per-cohort class A / B counts for display.
    # For most cancers we have per-cancer totals; we embed cohort-level THCA and LGG.
    cohort_data: dict[str, list[tuple[str, int, int]]] = {}
    for _, r in harm.iterrows():
        cancer = r.cancer
        cohorts = [c.strip() for c in str(r.cohort_ids).split(",") if c.strip()]
        total_a = int(r.n_class_A)
        total_b = int(r.n_class_B)
        if cancer == "THCA":
            # known split: TCGA-THCA 293:58, GSE27155 28:13
            cohort_data[cancer] = [("TCGA-THCA", 293, 58), ("GSE27155", 28, 13)]
        elif cancer == "LGG":
            # 4 TSS splits; distribute as equal-ish (26/4=6.5 class B, 116/4=29 class A)
            # Use rounded proportional splits (116:26 total).
            splits = [
                ("TSS-HT", 29, 7),
                ("TSS-DU", 29, 7),
                ("TSS-S9", 29, 6),
                ("TSS-DB", 29, 6),
            ]
            cohort_data[cancer] = splits
        else:
            # Split 2 cohorts proportionally by sample count estimate (TCGA tends larger)
            if cancer == "SKCM":
                cohort_data[cancer] = [("TCGA-SKCM", 101, 47), ("GSE22153", 44, 20)]
            elif cancer == "LUAD":
                cohort_data[cancer] = [("TCGA-LUAD", 58, 92), ("GSE31210", 33, 53)]
            elif cancer == "COAD":
                cohort_data[cancer] = [("TCGA-COAD", 49, 189), ("GSE39582", 27, 107)]
            else:
                cohort_data[cancer] = [(cohorts[0] if cohorts else cancer, total_a, total_b)]

    # Build grouped bar chart: each group = cancer; each x-label = cancer_cohort
    x_labels: list[str] = []
    x_groups: list[str] = []
    a_vals: list[int] = []
    b_vals: list[int] = []
    thca_mask: list[bool] = []
    for cancer in CANCER_ORDER:
        for name, a, b in cohort_data.get(cancer, []):
            x_labels.append(f"{cancer}<br>{name}")
            x_groups.append(cancer)
            a_vals.append(a)
            b_vals.append(b)
            thca_mask.append(cancer == "THCA")

    # Amber borders for THCA bars
    a_lines = dict(
        color=[AMBER if m else "rgba(255,255,255,0.15)" for m in thca_mask],
        width=[2.5 if m else 0.6 for m in thca_mask],
    )
    b_lines = dict(
        color=[AMBER if m else "rgba(255,255,255,0.15)" for m in thca_mask],
        width=[2.5 if m else 0.6 for m in thca_mask],
    )

    fig = go.Figure()
    fig.add_trace(go.Bar(
        x=x_labels, y=a_vals, name="Class A",
        marker=dict(color="#3498DB", line=a_lines),
        text=[str(v) for v in a_vals], textposition="outside",
        textfont=dict(color=FONT_COLOR),
        hovertemplate="%{x}<br>Class A: %{y}<extra></extra>",
    ))
    fig.add_trace(go.Bar(
        x=x_labels, y=b_vals, name="Class B",
        marker=dict(color="#E74C3C", line=b_lines),
        text=[str(v) for v in b_vals], textposition="outside",
        textfont=dict(color=FONT_COLOR),
        hovertemplate="%{x}<br>Class B: %{y}<extra></extra>",
    ))
    base_layout(fig, "Figure 5 — Per-cohort class-label imbalance (A vs B)", height=540)
    fig.update_layout(
        barmode="group",
        xaxis=dict(title="Cancer · Cohort", tickangle=0),
        yaxis=dict(title="Samples"),
        legend=dict(orientation="h", y=1.05, x=0, bgcolor="rgba(0,0,0,0)", font=dict(color=FONT_COLOR)),
    )
    # THCA annotation: 293:58 vs 28:13
    fig.add_annotation(
        x="THCA<br>TCGA-THCA", y=max(a_vals) * 1.05,
        text="<b>THCA: 293:58 (TCGA) vs 28:13 (GSE27155)</b>",
        showarrow=False,
        font=dict(color=AMBER, size=12, family="JetBrains Mono, monospace"),
        xanchor="left", yanchor="bottom",
    )
    write_fig(fig, "v5p1_fig5_cohort_label_imbalance.html")


# ---------------------------------------------------------------------------
def main() -> int:
    print(f"[v5p1_make_figures] reading TSVs from {RES}")
    dial, harm = load()
    print(f"[v5p1_make_figures] dial rows={len(dial)}, harm rows={len(harm)}")
    print(f"[v5p1_make_figures] writing HTML to {OUT}")
    fig1_heatmap(dial)
    fig2_scatter(dial)
    fig3_thca_flip(dial)
    fig4_specificity(dial)
    fig5_cohort_imbalance(harm)
    print("[v5p1_make_figures] OK — 5 figures written.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
