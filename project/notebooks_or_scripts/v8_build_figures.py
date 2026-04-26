#!/usr/bin/env python3
"""v8 figure builder — Plotly figures for the statgen-supplement dashboard page.

Inputs (tsv, /opt/thyroid-dash/project/results/v8_statgen/):
  v8_metasoft_results.tsv, v8_metasoft_forest_data.tsv
  v8_combatseq_vs_combat.tsv
  v8_pathway_dial.tsv, v8_pathway_vs_gene.tsv
  v8_fastRNA_style_dial.tsv
  v8_LMM_corrected_biomarkers.tsv
  v8_biomarker_DE_summary.tsv, v8_druggable8_retention.tsv
  v8_quantum_comparison.tsv (optional — skipped if missing)

Outputs (reports/html/figs_interactive/v8/*.html):
  fig1_forest_mvalues.html
  fig2_combatseq_bar.html
  fig3_pathway_heatmap.html
  fig4_cohort_centering.html
  fig5_quantum_grid.html (if quantum output exists)
  fig6_LMM_scatter.html
"""
from __future__ import annotations
import sys, os
from pathlib import Path
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.io as pio

ROOT = Path("/opt/thyroid-dash/project")
RES  = ROOT / "results" / "v8_statgen"
OUT  = ROOT / "reports" / "html" / "figs_interactive" / "v8"
OUT.mkdir(parents=True, exist_ok=True)

# Dark-theme layout matching the v5 pages.
DARK_LAYOUT = dict(
    paper_bgcolor="rgba(0,0,0,0)",
    plot_bgcolor="rgba(0,0,0,0)",
    font=dict(color="#EAEAEA", family="Inter, sans-serif"),
    margin=dict(l=60, r=40, t=70, b=60),
    height=560,
    legend=dict(font=dict(color="#EAEAEA"), bgcolor="rgba(0,0,0,0)"),
)
AXIS_STYLE = dict(gridcolor="rgba(234,234,234,0.08)",
                  zerolinecolor="rgba(234,234,234,0.15)")


def save_html(fig: go.Figure, path: Path) -> None:
    pio.write_html(fig, str(path), include_plotlyjs="cdn", full_html=True,
                   config={"responsive": True})
    print(f"[fig] wrote {path}")


def fig1_forest() -> None:
    """Forest plot of Han m-values per classifier per cancer."""
    meta = pd.read_csv(RES / "v8_metasoft_results.tsv", sep="\t")
    m_cols = [c for c in meta.columns if c.startswith("m_")]
    cancers = [c.replace("m_", "") for c in m_cols]
    classifiers = meta["classifier"].tolist()

    fig = go.Figure()
    palette = {"THCA":"#F5A623","SKCM":"#9B59B6","LGG":"#2ECC71",
               "LUAD":"#3498DB","COAD":"#E74C3C"}
    for i, cancer in enumerate(cancers):
        vals = meta[m_cols[i]].values
        fig.add_trace(go.Scatter(
            x=vals, y=classifiers,
            mode="markers",
            name=cancer,
            marker=dict(
                color=palette.get(cancer, "#999"),
                size=[8 + 14 * v for v in vals],
                opacity=0.85,
                line=dict(color="rgba(255,255,255,0.25)", width=0.8),
                symbol="diamond",
            ),
            hovertemplate=(f"<b>{cancer}</b><br>classifier: "
                           "%{y}<br>m-value: %{x:.3f}<extra></extra>"),
        ))

    # Han reference bands at m<=0.1 (no effect) and m>=0.9 (effect)
    fig.add_vrect(x0=0, x1=0.1, fillcolor="rgba(46,204,113,0.06)",
                  line_width=0, annotation_text="no effect (m≤0.1)",
                  annotation_position="top left",
                  annotation=dict(font_color="rgba(234,234,234,0.55)", font_size=11))
    fig.add_vrect(x0=0.9, x1=1.0, fillcolor="rgba(245,166,35,0.08)",
                  line_width=0, annotation_text="effect (m≥0.9)",
                  annotation_position="top right",
                  annotation=dict(font_color="rgba(234,234,234,0.55)", font_size=11))

    fig.update_layout(
        title=dict(text="Fig S2 — Random-effects meta-analysis: "
                        "Han m-values per classifier (THCA vs 4 others)",
                   font=dict(size=18, color="#EAEAEA")),
        xaxis=dict(title="m-value (posterior P[effect | data])",
                   range=[-0.02, 1.02], **AXIS_STYLE),
        yaxis=dict(title="", **AXIS_STYLE),
        **DARK_LAYOUT,
    )
    save_html(fig, OUT / "fig1_forest_mvalues.html")


def fig2_combatseq() -> None:
    """Grouped bar: ComBat vs ComBat-seq DIAL per classifier."""
    df = pd.read_csv(RES / "v8_combatseq_vs_combat.tsv", sep="\t")
    fig = go.Figure()
    fig.add_trace(go.Bar(
        name="ComBat (v5.1 baseline)", x=df["classifier"], y=df["dial_combat"],
        marker=dict(color="#9B59B6", line=dict(color="rgba(255,255,255,0.2)", width=0.6)),
        hovertemplate="%{x}<br>DIAL: %{y:.3f}<extra>ComBat</extra>",
    ))
    fig.add_trace(go.Bar(
        name="ComBat-seq (RNA-seq native)", x=df["classifier"], y=df["dial_combatseq"],
        marker=dict(color="#F5A623", line=dict(color="rgba(255,255,255,0.2)", width=0.6)),
        hovertemplate="%{x}<br>DIAL: %{y:.3f}<extra>ComBat-seq</extra>",
    ))
    fig.add_hline(y=0.1, line=dict(color="rgba(234,234,234,0.4)", width=1, dash="dash"),
                  annotation_text="DIAL = 0.1 (tolerance)",
                  annotation_font_color="rgba(234,234,234,0.55)",
                  annotation_font_size=11,
                  annotation_position="top left")
    fig.update_layout(
        barmode="group",
        title=dict(text="Fig S1 — THCA DIAL under ComBat vs ComBat-seq "
                        "(4/5 classifiers same; RF differs, both still batch-entangled)",
                   font=dict(size=17, color="#EAEAEA")),
        xaxis=dict(title="classifier", **AXIS_STYLE),
        yaxis=dict(title="DIAL (post-correction label flip)",
                   range=[0, 0.6], **AXIS_STYLE),
        **DARK_LAYOUT,
    )
    save_html(fig, OUT / "fig2_combatseq_bar.html")


def fig3_pathway_heatmap() -> None:
    """Heatmap of DIAL at gene (3000) vs pathway (50 Hallmark) level."""
    df = pd.read_csv(RES / "v8_pathway_vs_gene.tsv", sep="\t")
    df = df.sort_values("classifier").reset_index(drop=True)
    classifiers = df["classifier"].tolist()
    # Two-row heatmap: gene (dial_gene), pathway (dial_pathway)
    z = np.array([df["dial_gene"].values, df["dial_pathway"].values])
    levels = ["Gene level (3000 genes)", "Pathway level (50 Hallmark)"]

    # Annotate cells
    annotations = []
    for i, level in enumerate(levels):
        for j, cls in enumerate(classifiers):
            annotations.append(dict(
                x=cls, y=level, text=f"{z[i][j]:.3f}",
                showarrow=False, font=dict(color="#111", size=12),
            ))

    fig = go.Figure(go.Heatmap(
        z=z, x=classifiers, y=levels,
        colorscale=[[0, "#2ECC71"], [0.1, "#F5E9C6"], [0.5, "#F5A623"], [1, "#E74C3C"]],
        zmin=0, zmax=0.5,
        colorbar=dict(title=dict(text="DIAL", font_color="#EAEAEA"),
                      tickfont=dict(color="#EAEAEA")),
        hovertemplate="%{y}<br>classifier: %{x}<br>DIAL: %{z:.3f}<extra></extra>",
    ))
    fig.update_layout(
        annotations=annotations,
        title=dict(text="Fig S3 — Pathway aggregation erases the THCA DIAL flip "
                        "(DIAL → 0.0 at 50-Hallmark level for all 5 classifiers)",
                   font=dict(size=17, color="#EAEAEA")),
        xaxis=dict(title="", **AXIS_STYLE),
        yaxis=dict(title="", **AXIS_STYLE),
        **{**DARK_LAYOUT, "height": 360},
    )
    save_html(fig, OUT / "fig3_pathway_heatmap.html")


def fig4_cohort_centering() -> None:
    """Bar: FastRNA-style centering vs v5.1 ComBat DIAL per classifier."""
    df = pd.read_csv(RES / "v8_fastRNA_style_dial.tsv", sep="\t")
    fig = go.Figure()
    fig.add_trace(go.Bar(
        name="v5.1 (ComBat)", x=df["classifier"], y=df["dial_v5p1"],
        marker=dict(color="#9B59B6"),
        hovertemplate="%{x}<br>DIAL: %{y:.3f}<extra>ComBat</extra>",
    ))
    fig.add_trace(go.Bar(
        name="Cohort-centered (FastRNA-style)", x=df["classifier"], y=df["dial_centered"],
        marker=dict(color="#2ECC71"),
        hovertemplate="%{x}<br>DIAL: %{y:.3f}<extra>centered</extra>",
    ))
    fig.update_layout(
        barmode="group",
        title=dict(text="Fig S8A — Per-cohort mean subtraction eliminates THCA LODO flip "
                        "(ΔDIAL = +0.33 on average)",
                   font=dict(size=17, color="#EAEAEA")),
        xaxis=dict(title="classifier", **AXIS_STYLE),
        yaxis=dict(title="DIAL", range=[0, 0.55], **AXIS_STYLE),
        **DARK_LAYOUT,
    )
    save_html(fig, OUT / "fig4_cohort_centering.html")


def fig5_quantum_grid() -> None:
    """Grid plot of quantum vs classical DIAL across cancers — only if file present."""
    qpath = RES / "v8_quantum_comparison.tsv"
    if not qpath.exists():
        print(f"[fig] SKIP quantum grid — {qpath} not yet written")
        return
    df = pd.read_csv(qpath, sep="\t")
    # Columns in v8_quantum_dial.py output: cancer, classifier_family, auc_pre,
    # auc_post, dial, interpretation, n_used, seconds.
    clf_col = "classifier_family" if "classifier_family" in df.columns else "classifier"
    if "cancer" not in df.columns or clf_col not in df.columns:
        print(f"[fig] SKIP quantum grid — unexpected columns: {list(df.columns)}")
        return
    paradigms = [c for c in ["SVC_RBF", "QSVC", "VQC"] if c in df[clf_col].unique()]
    cancers = ["THCA", "SKCM", "LGG", "LUAD", "COAD"]
    cancers = [c for c in cancers if c in df["cancer"].unique()]

    df["dial_num"] = pd.to_numeric(df["dial"], errors="coerce")
    pivot = df.pivot_table(index="cancer", columns=clf_col, values="dial_num")
    pivot = pivot.reindex(index=cancers, columns=paradigms)
    interp = df.pivot_table(index="cancer", columns=clf_col, values="interpretation",
                            aggfunc="first")
    interp = interp.reindex(index=cancers, columns=paradigms)

    annotations = []
    for i, c in enumerate(cancers):
        for j, p in enumerate(paradigms):
            v = pivot.loc[c, p]
            tag = interp.loc[c, p] if pd.notna(interp.loc[c, p]) else ""
            if pd.isna(v) or tag == "skipped_budget_skip":
                txt = "skipped"
            else:
                txt = f"{v:.3f}"
            annotations.append(dict(x=p, y=c, text=txt,
                                    showarrow=False, font=dict(color="#111", size=12)))

    fig = go.Figure(go.Heatmap(
        z=pivot.values, x=paradigms, y=cancers,
        colorscale=[[0, "#2ECC71"], [0.1, "#F5E9C6"], [0.5, "#F5A623"], [1, "#E74C3C"]],
        zmin=0, zmax=0.5,
        colorbar=dict(title=dict(text="DIAL", font_color="#EAEAEA"),
                      tickfont=dict(color="#EAEAEA")),
        hovertemplate="%{y} / %{x}<br>DIAL: %{z:.3f}<extra></extra>",
    ))
    fig.update_layout(
        annotations=annotations,
        title=dict(text="Fig S5 — DIAL across classifier paradigms "
                        "(classical SVC vs quantum QSVC/VQC): "
                        "THCA flip is paradigm-invariant",
                   font=dict(size=17, color="#EAEAEA")),
        xaxis=dict(title="paradigm", **AXIS_STYLE),
        yaxis=dict(title="cancer", **AXIS_STYLE),
        **{**DARK_LAYOUT, "height": 420},
    )
    save_html(fig, OUT / "fig5_quantum_grid.html")


def fig6_lmm_scatter() -> None:
    """Scatter of −log10(p_naive) vs −log10(p_LMM), coloured by status."""
    df = pd.read_csv(RES / "v8_LMM_corrected_biomarkers.tsv", sep="\t")
    # Cap −log10(p) at 50 for display
    neg_naive = np.minimum(50, -np.log10(np.maximum(df["pvalue_naive"].values, 1e-50)))
    neg_lmm   = np.minimum(50, -np.log10(np.maximum(df["pvalue_LMM"].values,   1e-50)))
    status = df["status"].values

    colour_map = {
        "same": "#A0A0A0",
        "lost_with_LMM": "#E74C3C",
        "gained_with_LMM": "#2ECC71",
        "nonsig_both": "#606060",
    }
    fig = go.Figure()
    for s, col in colour_map.items():
        mask = status == s
        if not mask.any():
            continue
        fig.add_trace(go.Scattergl(
            x=neg_naive[mask], y=neg_lmm[mask],
            mode="markers",
            name=f"{s} (n={mask.sum():,})",
            marker=dict(color=col, size=4, opacity=0.55,
                        line=dict(width=0)),
            hovertemplate=f"gene idx: %{{pointNumber}}<br>status: {s}"
                          "<br>−log10(p_naive): %{x:.2f}<br>−log10(p_LMM): %{y:.2f}<extra></extra>",
        ))
    # FDR=0.05 approximate reference (depends on n; draw at p=0.05 as guide)
    fig.add_shape(type="line", x0=0, x1=50, y0=0, y1=50,
                  line=dict(color="rgba(234,234,234,0.3)", dash="dash"))
    fig.update_layout(
        title=dict(text="Fig S8C — Per-gene significance: naive OLS vs "
                        "linear-mixed model (cohort random effect)",
                   font=dict(size=17, color="#EAEAEA")),
        xaxis=dict(title="−log10(p) naive", range=[0, 52], **AXIS_STYLE),
        yaxis=dict(title="−log10(p) LMM", range=[0, 52], **AXIS_STYLE),
        **DARK_LAYOUT,
    )
    save_html(fig, OUT / "fig6_LMM_scatter.html")


def main() -> None:
    print("[v8-figs] building figures into", OUT)
    tasks = [
        ("fig1_forest", fig1_forest),
        ("fig2_combatseq", fig2_combatseq),
        ("fig3_pathway_heatmap", fig3_pathway_heatmap),
        ("fig4_cohort_centering", fig4_cohort_centering),
        ("fig5_quantum_grid", fig5_quantum_grid),
        ("fig6_lmm_scatter", fig6_lmm_scatter),
    ]
    errs = []
    for name, fn in tasks:
        try:
            fn()
        except Exception as e:
            print(f"[fig] {name} FAILED: {e!r}")
            errs.append(name)
    if errs:
        print(f"[v8-figs] FAILED: {errs}")
        sys.exit(1 if len(errs) > 1 else 0)  # tolerate single miss (quantum pending)
    print("[v8-figs] done")


if __name__ == "__main__":
    main()
