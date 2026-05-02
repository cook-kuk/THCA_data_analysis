"""Export Plotly Fig 10/11/12 (npj v7) to PDF + PNG using kaleido."""
from __future__ import annotations
import importlib.util
from pathlib import Path
import sys

ROOT = Path("/home/seungho/personal/THCA_data_analysis")

# We rebuild figs in-process to use write_image directly (the saved HTML is
# JS-rendered, not directly convertible).
def export_fig9():
    spec = importlib.util.spec_from_file_location(
        "fig9", ROOT / "project/notebooks_or_scripts/v17_audit_fig9_composite.py"
    )
    mod = importlib.util.module_from_spec(spec)
    sys.modules["fig9"] = mod
    spec.loader.exec_module(mod)


def export_fig10():
    spec = importlib.util.spec_from_file_location(
        "fig10", ROOT / "project/notebooks_or_scripts/v17_audit_fig10_chr7_forest.py"
    )
    mod = importlib.util.module_from_spec(spec)
    sys.modules["fig10"] = mod
    spec.loader.exec_module(mod)


def export_fig11():
    spec = importlib.util.spec_from_file_location(
        "fig11", ROOT / "project/notebooks_or_scripts/v17_audit_fig11.py"
    )
    mod = importlib.util.module_from_spec(spec)
    sys.modules["fig11"] = mod
    spec.loader.exec_module(mod)


# Cleaner: write a small re-render that builds the figs and saves PDF/PNG
import json
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

OUT = ROOT / "project/submission/npj/figures"
OUT.mkdir(parents=True, exist_ok=True)

# ------------------------------- Fig 10 (composite ROC + GSE213647 + coef)
def render_fig10():
    R8 = ROOT / "project/results/audit_2026_04_30/round8"
    summary = json.loads((R8 / "fig9_composite_summary.json").read_text())
    g3 = pd.read_csv(R8 / "r8_4_gse213647_composite_pred.tsv", sep="\t")

    # ROC curves cannot be reconstructed without raw — so re-train minimal
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import roc_curve, roc_auc_score
    from sklearn.preprocessing import StandardScaler
    import numpy as np

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
    X5 = StandardScaler().fit_transform(df[feats5].to_numpy())
    X4 = StandardScaler().fit_transform(df[feats4].to_numpy())
    m_dm1 = LogisticRegression(max_iter=2000).fit(X5, df["dm1"])
    m_fus = LogisticRegression(max_iter=2000).fit(X5, df["any_fusion"])
    m_dm1_4 = LogisticRegression(max_iter=2000).fit(X4, df["dm1"])
    fpr1, tpr1, _ = roc_curve(df["dm1"], m_dm1.predict_proba(X5)[:, 1])
    fpr2, tpr2, _ = roc_curve(df["any_fusion"], m_fus.predict_proba(X5)[:, 1])
    fpr3, tpr3, _ = roc_curve(df["dm1"], m_dm1_4.predict_proba(X4)[:, 1])
    auc1 = roc_auc_score(df["dm1"], m_dm1.predict_proba(X5)[:, 1])
    auc2 = roc_auc_score(df["any_fusion"], m_fus.predict_proba(X5)[:, 1])
    auc3 = roc_auc_score(df["dm1"], m_dm1_4.predict_proba(X4)[:, 1])
    coefs = m_dm1.coef_.ravel()

    fig = make_subplots(rows=1, cols=3, subplot_titles=(
        "A. TCGA composite ROC", "B. GSE213647 by histology", "C. 5-feat DM1 coefs"
    ))
    fig.add_trace(go.Scatter(x=fpr1, y=tpr1, mode="lines", name=f"DM1 5-feat AUC {auc1:.3f}",
                             line=dict(color="#d62728", width=3)), row=1, col=1)
    fig.add_trace(go.Scatter(x=fpr2, y=tpr2, mode="lines", name=f"fusion 5-feat AUC {auc2:.3f}",
                             line=dict(color="#1f77b4", width=2)), row=1, col=1)
    fig.add_trace(go.Scatter(x=fpr3, y=tpr3, mode="lines", name=f"DM1 4-feat AUC {auc3:.3f}",
                             line=dict(color="#ff7f0e", dash="dash", width=2)), row=1, col=1)
    fig.add_trace(go.Scatter(x=[0, 1], y=[0, 1], mode="lines", line=dict(color="grey", dash="dot"),
                             showlegend=False), row=1, col=1)
    color_map = {"Normal": "#2ca02c", "PDFP": "#bcbd22", "PTC": "#d62728", "UTC/ATC": "#000000"}
    for h in ["Normal", "PDFP", "PTC", "UTC/ATC"]:
        sub = g3[g3["histology"] == h]
        if len(sub) < 2:
            continue
        fig.add_trace(go.Violin(y=sub["composite_dm1_p"], name=f"{h} (n={len(sub)})",
                                box_visible=True, line_color=color_map[h], opacity=0.5),
                      row=1, col=2)
    fig.add_trace(go.Bar(x=feats5, y=coefs, marker_color=["#d62728" if c > 0 else "#1f77b4" for c in coefs],
                          text=[f"{c:.2f}" for c in coefs], textposition="auto",
                          showlegend=False), row=1, col=3)
    fig.update_layout(height=520, width=1500, template="plotly_white",
                      title="Fig 10. Composite RNA-only score (R8-4)")
    fig.write_image(str(OUT / "Fig10.pdf"), width=1500, height=520)
    fig.write_image(str(OUT / "Fig10.png"), width=1500, height=520, scale=3)
    print("Fig 10 ->", OUT)


def render_fig11():
    R9 = ROOT / "project/results/audit_2026_04_30/round9"
    boot = json.loads((R9 / "r9_3_summary.json").read_text())["bootstrap"]
    i131 = pd.DataFrame(json.loads((R9 / "r9_1_summary.json").read_text())["i131_by_dm"])
    k2_z = pd.read_csv(R9 / "r9_4_k2_within_z.tsv", sep="\t")
    rtk = pd.read_csv(R9 / "r9_2_chr7_rtk_summary.tsv", sep="\t")

    fig = make_subplots(rows=1, cols=4, subplot_titles=(
        "A. chr arm gain Δ DM1−not_DM (boot CI)", "B. chr7 RTK Cohen's d",
        "C. I-131 dose by DM", "D. K2 within-z calibration"
    ))
    arms = list(boot.keys())
    means = [boot[a]["mean_diff"] for a in arms]
    lows = [m - boot[a]["ci_low"] for m, a in zip(means, arms)]
    highs = [boot[a]["ci_high"] - m for m, a in zip(means, arms)]
    fig.add_trace(go.Scatter(x=means, y=arms, mode="markers",
                             error_x=dict(type="data", array=highs, arrayminus=lows),
                             marker=dict(size=14, color=["#d62728" if a in ("7p","7q") else "#888" for a in arms]),
                             showlegend=False), row=1, col=1)
    fig.add_vline(x=0, line=dict(color="grey", dash="dot"), row=1, col=1)
    rtk_s = rtk.sort_values("cohens_d_dm1_vs_notDM")
    fig.add_trace(go.Bar(y=rtk_s["gene"], x=rtk_s["cohens_d_dm1_vs_notDM"], orientation="h",
                          marker_color=["#d62728" if g in ("BRAF","EGFR","MET") else ("#1f77b4" if g in ("RET","RAF1") else "#888") for g in rtk_s["gene"]],
                          text=[f"d={d:.2f}" for d in rtk_s["cohens_d_dm1_vs_notDM"]], textposition="auto",
                          showlegend=False), row=1, col=2)
    fig.add_vline(x=0, line=dict(color="grey", dash="dot"), row=1, col=2)
    fig.add_trace(go.Bar(x=i131["cluster"], y=i131["pct_received_i131"],
                          marker_color=["#d62728","#ff7f0e","#888"],
                          text=[f"{r['n_with_i131']}/{r['n_total']} ({r['pct_received_i131']}%)" for _, r in i131.iterrows()],
                          textposition="auto", showlegend=False), row=1, col=3)
    fig.add_trace(go.Scatter(x=k2_z["rai_within_z_top4"], y=k2_z["p_DM2"], mode="markers",
                             marker=dict(color=k2_z["rai_within_z_top4"], colorscale="RdBu_r", size=6),
                             showlegend=False), row=1, col=4)
    fig.update_layout(height=560, width=1700, template="plotly_white",
                      title="Fig 11. chr7 + I-131 + K2 within-z (R8/R9)")
    fig.write_image(str(OUT / "Fig11.pdf"), width=1700, height=560)
    fig.write_image(str(OUT / "Fig11.png"), width=1700, height=560, scale=3)
    print("Fig 11 ->", OUT)


def render_fig12():
    R10 = ROOT / "project/results/audit_2026_04_30/round10"
    surv = pd.read_csv(R10 / "r10_2_3way_combo_summary.tsv", sep="\t")
    chr7 = pd.read_csv(R10 / "r10_3_chr7_x_fusion_4cell.tsv", sep="\t")
    msk = json.loads((R10 / "r10_1_msk_chr7_replication.json").read_text())
    tls_per = pd.read_csv(R10 / "r10_5_tcga_tls_per_sample.tsv", sep="\t")
    master = pd.read_csv(
        ROOT / "project/results/v17_tert_recovery/v2/sample_master_v17_tert_v2.tsv", sep="\t"
    )
    master["sample_short"] = master["sample_id"].str.slice(0, 15)
    master["dm"] = master["v17_dark_cluster"].fillna("not_DM").replace("", "not_DM")
    merged = tls_per.merge(master[["sample_short", "dm"]], on="sample_short", how="inner")

    fig = make_subplots(rows=1, cols=4, subplot_titles=(
        "A. TCGA Cabrita TLS by DM", "B. 3-way median OS",
        "C. chr7 × fus within DM1", "D. MSK chr7 RTK %Gain"
    ))
    color_map = {"DM1": "#d62728", "DM2": "#ff7f0e", "not_DM": "#888"}
    for c in ["DM1", "DM2", "not_DM"]:
        sub = merged[merged["dm"] == c]["TLS_score"].dropna()
        fig.add_trace(go.Violin(y=sub, name=f"{c} (n={len(sub)})", box_visible=True,
                                line_color=color_map[c], opacity=0.5), row=1, col=1)
    surv_s = surv.sort_values("median_os")
    cols = []
    for c in surv_s["combo"]:
        if c.startswith("DM1_fus1"): cols.append("#d62728")
        elif c.startswith("DM1"): cols.append("#ff7f0e")
        else: cols.append("#888")
    fig.add_trace(go.Bar(y=surv_s["combo"], x=surv_s["median_os"], orientation="h",
                          marker_color=cols, text=[f"n={n} ({int(e)} ev)" for n, e in zip(surv_s["n"], surv_s["events"])],
                          textposition="auto", showlegend=False), row=1, col=2)
    dm1_only = chr7[chr7["cohort"] == "DM1"].copy()
    dm1_only["label"] = dm1_only.apply(lambda r: f"chr7={r['chr7_gain']}|fus={r['fusion']}", axis=1)
    fig.add_trace(go.Bar(y=dm1_only["label"], x=dm1_only["n"], orientation="h",
                          marker_color=["#d62728" if r["chr7_gain"]==1 and r["fusion"]==1 else "#888" for _, r in dm1_only.iterrows()],
                          text=[f"n={r['n']}" for _, r in dm1_only.iterrows()], textposition="auto",
                          showlegend=False), row=1, col=3)
    msk_rows = []
    for st, p in msk.items():
        if "per_gene" not in p:
            continue
        for g in p["per_gene"]:
            msk_rows.append({"study": st, **g})
    msk_df = pd.DataFrame(msk_rows)
    if not msk_df.empty:
        for st in msk_df["study"].unique():
            sub = msk_df[msk_df["study"] == st]
            fig.add_trace(go.Bar(x=sub["gene"], y=sub["gain_amp_pct"], name=st), row=1, col=4)
    fig.update_layout(height=560, width=1800, template="plotly_white", barmode="group",
                      title="Fig 12. TLS replication + 3-way + chr7×fus + MSK (R10)")
    fig.write_image(str(OUT / "Fig12.pdf"), width=1800, height=560)
    fig.write_image(str(OUT / "Fig12.png"), width=1800, height=560, scale=3)
    print("Fig 12 ->", OUT)


if __name__ == "__main__":
    render_fig10()
    render_fig11()
    render_fig12()
