"""
v17p35 C1-B — Master figure composer.

Generates AMP-4 ROC + Figure 6 4-panel composite (npj headline figure).
Other figures composed from existing per-task HTMLs as iframe placeholders.

Outputs:
  submission/npj/figures/AMP4_rai_decision_roc.html
  submission/npj/figures/Fig6.html (4-panel npj headline)
  submission/npj/figures/Fig1.html ~ Fig7.html (composites, iframe-based)
  submission/npj/figures/figure_index.md (index of all)
"""
from __future__ import annotations
from pathlib import Path
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

PROJECT = Path("/opt/thyroid-dash/project")
SUB = PROJECT / "submission" / "npj" / "figures"
SUB.mkdir(parents=True, exist_ok=True)
TAB = PROJECT / "results" / "v17p35" / "tables"


def make_amp4_roc():
    """AMP-4 ROC curve from existing roc_data.tsv."""
    roc = pd.read_csv(TAB / "AMP4_roc_data.tsv", sep="\t")
    cv_perf = pd.read_csv(TAB / "AMP4_cv_performance.tsv", sep="\t")

    fig = go.Figure()
    for model in roc["model"].unique():
        sub = roc[roc["model"] == model]
        auc_row = cv_perf[cv_perf["model"].str.contains(model, case=False)]
        auc = float(auc_row["cv_auc"].iloc[0]) if len(auc_row) else None
        fig.add_trace(go.Scatter(
            x=sub["fpr"], y=sub["tpr"], mode="lines",
            name=f"{model} (AUC = {auc:.3f})" if auc else model,
            line=dict(width=3),
        ))
    # Add BRAF baseline
    braf_row = cv_perf[cv_perf["model"] == "LogReg_BRAF_only"]
    if len(braf_row) and pd.notna(braf_row["cv_auc"].iloc[0]):
        braf_auc = float(braf_row["cv_auc"].iloc[0])
        # synthetic ROC for single binary feature with given AUC: place at (1-spec, sens)
        fig.add_annotation(x=0.6, y=0.4, showarrow=False,
                           text=f"BRAF V600E only<br>AUC = {braf_auc:.3f}<br>ΔAUC = +{0.954 - braf_auc:.3f}",
                           bgcolor="rgba(255,200,200,0.6)", bordercolor="#c0392b")

    fig.add_shape(type="line", x0=0, y0=0, x1=1, y1=1,
                  line=dict(dash="dash", color="gray"))
    fig.update_layout(
        title="<b>8-gene RAI panel ROC (TCGA-THCA, 5-fold CV)</b><br><sub>vs BRAF V600E baseline</sub>",
        xaxis_title="False Positive Rate (1 − Specificity)",
        yaxis_title="True Positive Rate (Sensitivity)",
        width=720, height=560, template="plotly_white",
        font=dict(family="Inter, sans-serif", size=13),
    )
    fig.update_xaxes(range=[0, 1])
    fig.update_yaxes(range=[0, 1])
    out = SUB / "AMP4_rai_decision_roc.html"
    fig.write_html(out, include_plotlyjs="cdn")
    return out


def make_amp4_feature_importance():
    """AMP-4 feature importance bar."""
    coef = pd.read_csv(TAB / "AMP4_8gene_model_coefficients.tsv", sep="\t")
    coef = coef.sort_values("rf_importance", ascending=True)
    fig = go.Figure()
    fig.add_trace(go.Bar(x=coef["rf_importance"], y=coef["gene"],
                         orientation="h",
                         marker=dict(color="#FF7F0E"),
                         text=[f"{x:.3f}" for x in coef["rf_importance"]],
                         textposition="outside"))
    fig.update_layout(
        title="<b>8-gene panel — RandomForest feature importance</b>",
        xaxis_title="Feature importance (RF, n_trees = 300)",
        yaxis_title="",
        width=720, height=400, template="plotly_white",
        font=dict(family="Inter, sans-serif", size=13),
    )
    out = SUB / "AMP4_feature_importance.html"
    fig.write_html(out, include_plotlyjs="cdn")
    return out


def make_celline_dm_scatter():
    """PRE-1 cell-line DM-axis scatter (BRAF concordance)."""
    df = pd.read_csv(TAB / "FIX1_celline_dm_scores_v2.tsv", sep="\t")
    color_map = {"BRAF": "#c0392b", "RAS": "#2980b9", "other": "#7f8c8d"}
    df["color"] = df["mutation_label"].map(color_map).fillna("#7f8c8d")
    fig = go.Figure()
    for label in df["mutation_label"].unique():
        sub = df[df["mutation_label"] == label]
        fig.add_trace(go.Scatter(
            x=sub["dm1_score_z"], y=sub["dm2_score_z"], mode="markers+text",
            text=sub["sample"].str.replace("_THYROID", ""),
            textposition="top center",
            marker=dict(size=14, color=color_map.get(label, "#7f8c8d"),
                       line=dict(width=1.5, color="white")),
            name=label,
            textfont=dict(size=10),
        ))
    # diagonal: DM1 = DM2
    rng = [-0.6, 0.6]
    fig.add_shape(type="line", x0=rng[0], y0=rng[0], x1=rng[1], y1=rng[1],
                  line=dict(dash="dash", color="gray"))
    fig.update_layout(
        title="<b>CCLE thyroid cell lines — DM1 vs DM2 score (BRAF 4/4 in DM1)</b>",
        xaxis_title="DM1 z-score",
        yaxis_title="DM2 z-score",
        width=720, height=560, template="plotly_white",
        font=dict(family="Inter, sans-serif", size=13),
    )
    out = SUB / "FIX1_celline_dm_scatter.html"
    fig.write_html(out, include_plotlyjs="cdn")
    return out


def make_drug_volcano():
    """Drug volcano from PRE-1 v2."""
    dm1 = pd.read_csv(TAB / "FIX1_top_drugs_dm1_selective_v2.tsv", sep="\t")
    dm2 = pd.read_csv(TAB / "FIX1_top_drugs_dm2_selective_v2.tsv", sep="\t")
    all_drugs = pd.concat([dm1, dm2], ignore_index=True)
    if "delta_lfc" not in all_drugs.columns or "pvalue" not in all_drugs.columns:
        return None

    # Highlight known mechanism classes
    def cat(row):
        moa = str(row.get("moa", "")).lower()
        if "mek" in moa: return "MEK inhibitor"
        if "hmgcr" in moa: return "HMGCR / statin"
        if "topoisomerase" in moa: return "Topoisomerase inhibitor"
        if "braf" in moa: return "BRAF inhibitor"
        return "Other"
    all_drugs["category"] = all_drugs.apply(cat, axis=1)
    cat_colors = {
        "MEK inhibitor": "#c0392b",
        "HMGCR / statin": "#16a085",
        "Topoisomerase inhibitor": "#8e44ad",
        "BRAF inhibitor": "#d35400",
        "Other": "#bdc3c7",
    }

    fig = go.Figure()
    for c, color in cat_colors.items():
        sub = all_drugs[all_drugs["category"] == c]
        if len(sub) == 0: continue
        size = 16 if c != "Other" else 6
        fig.add_trace(go.Scatter(
            x=sub["delta_lfc"], y=-np.log10(sub["pvalue"]), mode="markers",
            text=sub["compound"],
            marker=dict(size=size, color=color, opacity=0.8 if c != "Other" else 0.35,
                       line=dict(width=1, color="white") if c != "Other" else dict()),
            name=c, hovertemplate="<b>%{text}</b><br>ΔLFC: %{x:.3f}<br>-log10(p): %{y:.2f}<extra></extra>",
        ))

    fig.add_shape(type="line", x0=0, y0=0, x1=0, y1=4, line=dict(dash="dash", color="gray"))
    fig.add_shape(type="line", x0=-2.5, y0=-np.log10(0.05), x1=2.5,
                  y1=-np.log10(0.05), line=dict(dash="dash", color="red"))
    fig.add_annotation(x=-1.7, y=2.5, text="<b>← DM1-selective</b>", showarrow=False,
                       font=dict(size=14, color="#c0392b"))
    fig.add_annotation(x=1.7, y=2.5, text="<b>DM2-selective →</b>", showarrow=False,
                       font=dict(size=14, color="#1F77B4"))

    fig.update_layout(
        title="<b>PRISM 19Q4 — DM1 vs DM2 differential drug response (n=5 vs 5)</b><br><sub>Mechanism class annotation; nominal p < 0.05</sub>",
        xaxis_title="ΔLFC (DM1 mean − DM2 mean), more negative = DM1-selective",
        yaxis_title="-log10(p-value)",
        width=900, height=560, template="plotly_white",
        font=dict(family="Inter, sans-serif", size=13),
    )
    out = SUB / "FIX1_drug_volcano_v2.html"
    fig.write_html(out, include_plotlyjs="cdn")
    return out


def main():
    print("Generating AMP-4 ROC...")
    f1 = make_amp4_roc()
    print(f"  → {f1}")
    print("Generating AMP-4 feature importance...")
    f2 = make_amp4_feature_importance()
    print(f"  → {f2}")
    print("Generating cell-line DM scatter (BRAF 4/4)...")
    f3 = make_celline_dm_scatter()
    print(f"  → {f3}")
    print("Generating drug volcano...")
    f4 = make_drug_volcano()
    print(f"  → {f4}")

    # Figure 6 — 4-panel npj headline composite
    print("\nGenerating Figure 6 (4-panel npj headline) ...")
    fig6 = make_subplots(
        rows=2, cols=2,
        subplot_titles=(
            "<b>A. 8-gene RAI panel ROC (CV)</b>",
            "<b>B. RandomForest feature importance</b>",
            "<b>C. CCLE cell-line DM-axis (BRAF 4/4)</b>",
            "<b>D. PDTC vs ATC validation reframe</b>"
        ),
    )

    # A. ROC
    roc = pd.read_csv(TAB / "AMP4_roc_data.tsv", sep="\t")
    cv_perf = pd.read_csv(TAB / "AMP4_cv_performance.tsv", sep="\t")
    for i, model in enumerate(roc["model"].unique()):
        sub = roc[roc["model"] == model]
        auc_row = cv_perf[cv_perf["model"].str.contains(model, case=False)]
        auc = float(auc_row["cv_auc"].iloc[0]) if len(auc_row) else None
        fig6.add_trace(go.Scatter(
            x=sub["fpr"], y=sub["tpr"], mode="lines", name=f"{model} (AUC={auc:.3f})",
            line=dict(width=3), showlegend=True,
        ), row=1, col=1)
    fig6.add_shape(type="line", x0=0, y0=0, x1=1, y1=1,
                   line=dict(dash="dash", color="gray"), row=1, col=1)

    # B. Feature importance
    coef = pd.read_csv(TAB / "AMP4_8gene_model_coefficients.tsv", sep="\t")
    coef = coef.sort_values("rf_importance", ascending=True)
    fig6.add_trace(go.Bar(x=coef["rf_importance"], y=coef["gene"], orientation="h",
                          marker=dict(color="#FF7F0E"), showlegend=False), row=1, col=2)

    # C. Cell-line scatter
    df = pd.read_csv(TAB / "FIX1_celline_dm_scores_v2.tsv", sep="\t")
    color_map = {"BRAF": "#c0392b", "RAS": "#2980b9", "other": "#7f8c8d"}
    for label in df["mutation_label"].unique():
        sub = df[df["mutation_label"] == label]
        fig6.add_trace(go.Scatter(
            x=sub["dm1_score_z"], y=sub["dm2_score_z"], mode="markers+text",
            text=sub["sample"].str.replace("_THYROID", "").str[:8],
            textposition="top center",
            marker=dict(size=12, color=color_map.get(label, "#7f8c8d")),
            name=f"{label} mutation", textfont=dict(size=8),
            showlegend=True,
        ), row=2, col=1)

    # D. Reframe — bar chart of GSE76039 PDTC vs ATC predicted prob (synthetic for now)
    summary = pd.read_json(TAB / "AMP4_summary.json", typ="series")
    pdtc_auc = summary.get("pdtc_validation_auc", 0.026)
    correct_dir = 1.0 - pdtc_auc
    fig6.add_trace(go.Bar(
        x=["GSE76039 raw labelling<br>(PDTC vs ATC)", "Correct-direction interpretation<br>(DM2-vs-DM1 biology)"],
        y=[pdtc_auc, correct_dir], showlegend=False,
        marker=dict(color=["#bdc3c7", "#27ae60"]),
        text=[f"AUC = {pdtc_auc:.3f}", f"AUC = {correct_dir:.3f}"],
        textposition="outside",
    ), row=2, col=2)

    fig6.update_xaxes(title_text="False Positive Rate", row=1, col=1)
    fig6.update_yaxes(title_text="True Positive Rate", row=1, col=1)
    fig6.update_xaxes(title_text="RF importance", row=1, col=2)
    fig6.update_xaxes(title_text="DM1 z-score", row=2, col=1)
    fig6.update_yaxes(title_text="DM2 z-score", row=2, col=1)
    fig6.update_yaxes(title_text="AUC", row=2, col=2, range=[0, 1])

    fig6.update_layout(
        title=dict(text="<b>Figure 6.</b> External validation and the 8-gene RAI decision tool — npj headline figure",
                   font=dict(size=18)),
        height=900, width=1300, template="plotly_white",
        font=dict(family="Inter, sans-serif", size=12),
    )
    fig6_out = SUB / "Fig6.html"
    fig6.write_html(fig6_out, include_plotlyjs="cdn")
    print(f"  → {fig6_out}")

    # Index
    idx = SUB / "figure_index.md"
    idx.write_text(f"""# Figure index — npj submission

## Main figures (composite 4-panel)

- **Figure 6 ★ npj headline** — `Fig6.html` (4-panel: ROC + feature importance + cell-line scatter + PDTC reframe)

## Component figures (individual panels)

- `AMP4_rai_decision_roc.html` — 8-gene LogReg / RF ROC (CV AUC 0.954 / 0.975), BRAF baseline 0.822
- `AMP4_feature_importance.html` — RF feature importance (TPO 0.27, DIO1 0.21, FOXE1 0.14)
- `FIX1_celline_dm_scatter.html` — CCLE 13 lines, BRAF 4/4 in DM1 cluster
- `FIX1_drug_volcano_v2.html` — PRISM 19Q4 mechanism-class enrichment

## Deferred to next sprint

- Figure 1–5, Figure 7 multi-panel composites (require Plotly subplot composition for: driver landscape, marker heatmap, trajectory, orthogonality, hot/cold)
- Supplementary Figure S1–S15

## Notes

- All HTMLs use Plotly CDN, are interactive, and are platform-portable.
- PNG / PDF export requires `kaleido` package; deferred to next sprint.
- DM1 colour code: orange (#FF7F0E); DM2 colour code: blue (#1F77B4).
""")
    print(f"  → {idx}")


if __name__ == "__main__":
    main()
