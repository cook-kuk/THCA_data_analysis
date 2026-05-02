"""
Fig 9 candidate — composite RNA score visualization.
Panel A: TCGA train ROC (5-feat DM1, 5-feat fusion, 4-feat DM1)
Panel B: GSE213647 composite-DM1 probability density by histology
Panel C: Coefficient bar (5-feat DM1)
"""
from __future__ import annotations
import json
from pathlib import Path
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_curve, roc_auc_score
from sklearn.preprocessing import StandardScaler

ROOT = Path("/home/seungho/personal/THCA_data_analysis")
R8 = ROOT / "project/results/audit_2026_04_30/round8"
FIG = ROOT / "project/reports/html/figs_interactive/v17"
FIG.mkdir(parents=True, exist_ok=True)


def main() -> None:
    # rebuild composite to get ROC curves
    master = pd.read_csv(
        ROOT / "project/results/v17_tert_recovery/v2/sample_master_v17_tert_v2.tsv", sep="\t"
    )
    master["sample_short"] = master["sample_id"].str.slice(0, 15)
    master["patient12"] = master["sample_short"].str.slice(0, 12)
    hla = pd.read_csv(ROOT / "project/results/v17_hla/tcga_thca_hla_per_sample.tsv", sep="\t")
    hla["patient12"] = hla["sample_id"].str.slice(0, 12)
    meth = pd.read_csv(R8.parent / "round5/r5_2_sample_methylation_8gene.tsv", sep="\t")
    meth["patient12"] = meth["sample_short"].str.slice(0, 12)
    sv = pd.read_csv(R8.parent / "round3/cbio_sv_thca.tsv", sep="\t")
    fusion_set = set(sv["sampleId"].unique())
    master["any_fusion"] = master["sample_short"].isin(fusion_set).astype(int)
    df = (
        master[["patient12", "rai_score_v17", "age", "v17_dark_cluster", "any_fusion"]]
        .merge(hla[["patient12", "hla_class_I_score", "hla_class_II_score"]], on="patient12", how="inner")
        .merge(meth[["patient12", "mean_8g_beta"]], on="patient12", how="inner")
        .dropna(subset=["rai_score_v17", "age", "hla_class_I_score", "hla_class_II_score", "mean_8g_beta"])
    )
    df["dm1"] = (df["v17_dark_cluster"] == "DM1").astype(int)

    feats5 = ["rai_score_v17", "hla_class_I_score", "hla_class_II_score", "mean_8g_beta", "age"]
    feats4 = ["rai_score_v17", "hla_class_I_score", "hla_class_II_score", "age"]

    sc5 = StandardScaler().fit(df[feats5].to_numpy())
    X5 = sc5.transform(df[feats5].to_numpy())
    m_dm1 = LogisticRegression(max_iter=2000).fit(X5, df["dm1"])
    m_fus = LogisticRegression(max_iter=2000).fit(X5, df["any_fusion"])

    sc4 = StandardScaler().fit(df[feats4].to_numpy())
    X4 = sc4.transform(df[feats4].to_numpy())
    m_dm1_4 = LogisticRegression(max_iter=2000).fit(X4, df["dm1"])

    # ROC curves
    p_dm1 = m_dm1.predict_proba(X5)[:, 1]
    p_fus = m_fus.predict_proba(X5)[:, 1]
    p_dm1_4 = m_dm1_4.predict_proba(X4)[:, 1]
    fpr1, tpr1, _ = roc_curve(df["dm1"], p_dm1)
    fpr2, tpr2, _ = roc_curve(df["any_fusion"], p_fus)
    fpr3, tpr3, _ = roc_curve(df["dm1"], p_dm1_4)
    auc1 = roc_auc_score(df["dm1"], p_dm1)
    auc2 = roc_auc_score(df["any_fusion"], p_fus)
    auc3 = roc_auc_score(df["dm1"], p_dm1_4)

    # GSE213647 composite densities by histology
    g3 = pd.read_csv(R8 / "r8_4_gse213647_composite_pred.tsv", sep="\t")

    # build figure
    fig = make_subplots(
        rows=1,
        cols=3,
        subplot_titles=(
            "A. TCGA composite ROC",
            "B. GSE213647 composite-DM1 by histology",
            "C. 5-feature DM1 coefficients (z-scored)",
        ),
        column_widths=[0.33, 0.34, 0.33],
    )

    fig.add_trace(
        go.Scatter(
            x=fpr1, y=tpr1, mode="lines", name=f"DM1 (5-feat) AUC={auc1:.3f}",
            line=dict(color="#d62728", width=3),
        ),
        row=1, col=1,
    )
    fig.add_trace(
        go.Scatter(
            x=fpr2, y=tpr2, mode="lines", name=f"any-fusion (5-feat) AUC={auc2:.3f}",
            line=dict(color="#1f77b4", width=2),
        ),
        row=1, col=1,
    )
    fig.add_trace(
        go.Scatter(
            x=fpr3, y=tpr3, mode="lines", name=f"DM1 (4-feat no-meth) AUC={auc3:.3f}",
            line=dict(color="#ff7f0e", width=2, dash="dash"),
        ),
        row=1, col=1,
    )
    fig.add_trace(
        go.Scatter(x=[0, 1], y=[0, 1], mode="lines", line=dict(color="grey", dash="dot"), showlegend=False),
        row=1, col=1,
    )
    fig.update_xaxes(title="FPR", row=1, col=1)
    fig.update_yaxes(title="TPR", row=1, col=1)

    # B
    color_map = {"Normal": "#2ca02c", "PDFP": "#bcbd22", "PTC": "#d62728", "UTC/ATC": "#000000"}
    for h in ["Normal", "PDFP", "PTC", "UTC/ATC"]:
        sub = g3[g3["histology"] == h]
        if len(sub) < 2:
            continue
        fig.add_trace(
            go.Violin(
                y=sub["composite_dm1_p"], name=f"{h} (n={len(sub)})", box_visible=True,
                meanline_visible=True, line_color=color_map.get(h, "#888"),
                fillcolor=color_map.get(h, "#888"), opacity=0.45,
            ),
            row=1, col=2,
        )
    fig.update_yaxes(title="composite DM1 prob.", row=1, col=2)

    # C: coefficients
    coefs = m_dm1.coef_.ravel()
    fig.add_trace(
        go.Bar(
            x=feats5, y=coefs,
            marker_color=["#d62728" if c > 0 else "#1f77b4" for c in coefs],
            text=[f"{c:.2f}" for c in coefs], textposition="auto",
            showlegend=False,
        ),
        row=1, col=3,
    )
    fig.update_yaxes(title="LogReg coef (z-scored)", row=1, col=3)

    fig.update_layout(
        height=520, width=1500,
        title=f"Fig 9 (candidate). RNA-only composite score for DM1 / fusion+ identification — TCGA train (n={len(df)}) + GSE213647 cross-cohort apply (n={len(g3)})",
        template="plotly_white", margin=dict(l=60, r=20, t=80, b=60),
    )
    out_html = FIG / "v17_fig9_composite_audit.html"
    fig.write_html(str(out_html), include_plotlyjs="cdn")
    print(f"wrote {out_html}")
    summary = dict(
        auc_dm1_5feat=round(auc1, 3),
        auc_fusion_5feat=round(auc2, 3),
        auc_dm1_4feat=round(auc3, 3),
        n_tcga=len(df),
        n_gse213647=len(g3),
        coefs_5feat={k: round(float(v), 3) for k, v in zip(feats5, coefs)},
    )
    (R8 / "fig9_composite_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
