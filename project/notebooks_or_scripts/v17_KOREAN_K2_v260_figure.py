#!/usr/bin/env python3
"""K2 Korean n=260 figure — replaces v17_KOREAN_K2_figure.py (n=9 pilot)."""
from __future__ import annotations
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import plotly.graph_objects as go
from plotly.subplots import make_subplots

FIG_DIR = Path("/opt/thyroid-dash/project/submission/npj/figures")
HTML_DIR = Path("/opt/thyroid-dash/project/reports/html/figs_interactive/v17")
FIG_DIR.mkdir(parents=True, exist_ok=True)
HTML_DIR.mkdir(parents=True, exist_ok=True)

GENE_8 = ["SLC5A5", "TPO", "TG", "TSHR", "PAX8", "NKX2-1", "FOXE1", "DIO1"]
DM1_C = "#1f77b4"
DM2_C = "#d62728"


def load_tcga_centered():
    tpm = pd.read_csv(
        "/data/thca/data_processed/bulk_rnaseq/TCGA-THCA_rnaseq_expression_log2.tsv",
        sep="\t", index_col=0,
    )
    lbl = pd.read_csv(
        "/opt/thyroid-dash/project/results/v17_realfix/R1A_cluster_labels.tsv", sep="\t"
    )
    g_have = [g for g in GENE_8 if g in tpm.index]
    samples = [s for s in lbl["sample_id"] if s in tpm.columns]
    X = tpm.loc[g_have, samples].T.values
    y = lbl.set_index("sample_id").loc[samples, "cluster"].str.startswith("DM2").astype(int).values
    Xc = X - X.mean(axis=1, keepdims=True)
    return g_have, Xc, y


def load_korean():
    mat = pd.read_csv(
        "/opt/thyroid-dash/project/results/v17_korean/K2_8gene_tpm_matrix_v4.tsv",
        sep="\t", index_col=0,
    )
    pred = pd.read_csv(
        "/opt/thyroid-dash/project/results/v17_korean/K2_korean_predictions_v4.tsv",
        sep="\t",
    )
    summary = json.loads(
        Path("/opt/thyroid-dash/project/results/v17_korean/K2_korean_summary_v4.json").read_text()
    )
    return mat, pred, summary


def build_png(g_have, Xc_tcga, y_tcga, mat, pred, summary):
    Xk = np.log2(mat[g_have].values + 1.0)
    Xkc = Xk - Xk.mean(axis=1, keepdims=True)
    pred_sorted = pred.sort_values("p_DM2", ascending=False).reset_index(drop=True)
    korean_idx_map = {r: i for i, r in enumerate(mat.index)}
    korean_idx_sorted = [korean_idx_map[r] for r in pred_sorted["run"]]
    Xkc_sorted = Xkc[korean_idx_sorted]

    dm2_centroid = Xc_tcga[y_tcga == 1].mean(0)
    dm1_centroid = Xc_tcga[y_tcga == 0].mean(0)
    n_dm1, n_dm2 = int((y_tcga == 0).sum()), int((y_tcga == 1).sum())
    n_kor = len(pred_sorted)

    top_n = 15
    show_idx = list(range(top_n)) + list(range(n_kor - top_n, n_kor))
    rows_label = (
        [f"TCGA DM2 centroid (n={n_dm2})", f"TCGA DM1 centroid (n={n_dm1})"]
        + [f"{pred_sorted.iloc[i]['run']}  (p={pred_sorted.iloc[i]['p_DM2']:.2f})" for i in show_idx]
    )
    z = np.vstack([dm2_centroid, dm1_centroid, Xkc_sorted[show_idx]])

    mpl.rcParams.update({
        "font.family": ["DejaVu Sans", "Arial Unicode MS"],
        "font.size": 10,
        "axes.titlesize": 11.5,
    })

    fig = plt.figure(figsize=(16, 11), dpi=120)
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.15], width_ratios=[1.0, 1.4],
                          hspace=0.32, wspace=0.28)

    ax_h = fig.add_subplot(gs[0, 0])
    p_vals = pred_sorted["p_DM2"].values
    bins = np.linspace(0, 1, 21)
    counts, edges, patches = ax_h.hist(p_vals, bins=bins, edgecolor="black", linewidth=0.4)
    for c, patch in zip(edges[:-1], patches):
        patch.set_facecolor(DM1_C if c < 0.5 else DM2_C)
    ax_h.axvline(0.5, color="#444", linestyle="--", linewidth=1.0)
    ax_h.set_xlabel("predicted p(DM2)")
    ax_h.set_ylabel("Korean samples (count)")
    ax_h.set_title(f"a   p(DM2) distribution — n = {n_kor}", loc="left", fontweight="bold")
    n_dm1_kor = int((p_vals < 0.5).sum())
    n_dm2_kor = int((p_vals >= 0.5).sum())
    ax_h.text(0.02, 0.95, f"DM1 calls: {n_dm1_kor}\nDM2 calls: {n_dm2_kor}\nmean p = {p_vals.mean():.3f}\nmedian p = {np.median(p_vals):.3f}",
              transform=ax_h.transAxes, va="top", fontsize=9.5,
              bbox=dict(facecolor="white", edgecolor="#999", boxstyle="round,pad=0.4"))
    ax_h.spines["top"].set_visible(False); ax_h.spines["right"].set_visible(False)

    ax_v = fig.add_subplot(gs[0, 1])
    box_data = []
    box_pos = []
    box_colors = []
    for j, g in enumerate(g_have):
        tcga_dm1 = Xc_tcga[y_tcga == 0, j]
        tcga_dm2 = Xc_tcga[y_tcga == 1, j]
        kor = Xkc[:, j]
        box_data.extend([tcga_dm1, tcga_dm2, kor])
        base = j * 4
        box_pos.extend([base, base + 1, base + 2])
        box_colors.extend([DM1_C, DM2_C, "#10b981"])
    bp = ax_v.boxplot(box_data, positions=box_pos, widths=0.85, patch_artist=True,
                      medianprops=dict(color="black", linewidth=1.2), showfliers=False)
    for patch, c in zip(bp["boxes"], box_colors):
        patch.set_facecolor(c); patch.set_alpha(0.65); patch.set_edgecolor("black"); patch.set_linewidth(0.4)
    ax_v.set_xticks([j * 4 + 1 for j in range(len(g_have))])
    ax_v.set_xticklabels(g_have, rotation=18, ha="right")
    ax_v.set_ylabel("centered log2(TPM+1)")
    ax_v.set_title(f"b   per-gene centered profile — TCGA DM1 (blue, n={n_dm1}) · DM2 (red, n={n_dm2}) · Korean (green, n={n_kor})", loc="left", fontweight="bold", fontsize=10.5)
    ax_v.axhline(0, color="#999", linewidth=0.6, linestyle="--")
    ax_v.spines["top"].set_visible(False); ax_v.spines["right"].set_visible(False)

    ax_hm = fig.add_subplot(gs[1, :])
    vmax = float(np.nanmax(np.abs(z)))
    im = ax_hm.imshow(z, cmap="RdBu_r", aspect="auto", vmin=-vmax, vmax=vmax)
    ax_hm.set_xticks(range(len(g_have)))
    ax_hm.set_xticklabels(g_have, rotation=18, ha="right")
    ax_hm.set_yticks(range(len(rows_label)))
    ax_hm.set_yticklabels(rows_label, fontsize=8.5)
    ax_hm.axhline(1.5, color="black", linewidth=1.4)
    ax_hm.axhline(1.5 + top_n, color="#444", linewidth=0.8, linestyle=":")
    ax_hm.set_title(f"c   centered profile heatmap — TCGA centroids + top {top_n} highest p(DM2) + bottom {top_n} lowest p(DM2) Korean samples", loc="left", fontweight="bold", fontsize=10.5)
    cbar = fig.colorbar(im, ax=ax_hm, fraction=0.025, pad=0.014)
    cbar.set_label("centered log2(TPM+1)", fontsize=9)
    for i in range(z.shape[0]):
        for j in range(z.shape[1]):
            ax_hm.text(j, i, f"{z[i, j]:+.1f}", ha="center", va="center",
                       fontsize=7.0, color="black" if abs(z[i, j]) < vmax * 0.55 else "white")

    auc_cv = summary["tcga_5fold_cv_auc_mean"]
    auc_std = summary["tcga_5fold_cv_auc_std"]
    fig.suptitle(
        f"Fig K2.  PRJEB11591 Korean PTC validation (n = {n_kor}) — within-sample-centered 8-gene profile vs TCGA DM1/DM2  ·  TCGA 5-fold CV AUC = {auc_cv:.3f} ± {auc_std:.3f}  ·  Korean call: DM1 {n_dm1_kor} / DM2 {n_dm2_kor}",
        fontsize=12.5, fontweight="bold", y=0.995,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.965))
    out_png = FIG_DIR / "Korean_K2_v260.png"
    fig.savefig(out_png, dpi=180, bbox_inches="tight")
    fig.savefig(FIG_DIR / "Korean_K2_v260.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"  ✓ {out_png}")
    return n_dm1_kor, n_dm2_kor


def build_html(g_have, Xc_tcga, y_tcga, mat, pred, summary):
    Xk = np.log2(mat[g_have].values + 1.0)
    Xkc = Xk - Xk.mean(axis=1, keepdims=True)
    pred_sorted = pred.sort_values("p_DM2", ascending=False).reset_index(drop=True)
    korean_idx_map = {r: i for i, r in enumerate(mat.index)}
    Xkc_sorted = Xkc[[korean_idx_map[r] for r in pred_sorted["run"]]]
    dm2_centroid = Xc_tcga[y_tcga == 1].mean(0)
    dm1_centroid = Xc_tcga[y_tcga == 0].mean(0)
    n_dm1, n_dm2 = int((y_tcga == 0).sum()), int((y_tcga == 1).sum())
    n_kor = len(pred_sorted)

    top_n = 20
    show = list(range(top_n)) + list(range(n_kor - top_n, n_kor))
    rows_label = (
        [f"TCGA DM2 centroid (n={n_dm2})", f"TCGA DM1 centroid (n={n_dm1})"]
        + [f"{pred_sorted.iloc[i]['run']} (p={pred_sorted.iloc[i]['p_DM2']:.2f})" for i in show]
    )
    z = np.vstack([dm2_centroid, dm1_centroid, Xkc_sorted[show]])

    fig = make_subplots(rows=1, cols=2, column_widths=[0.42, 0.58], horizontal_spacing=0.10,
                        subplot_titles=(f"<b>p(DM2) distribution — n={n_kor}</b>",
                                        f"<b>centered profile heatmap — TCGA centroids + 40 representative Korean</b>"))

    p_vals = pred_sorted["p_DM2"].values
    bins = np.linspace(0, 1, 21)
    centers = (bins[:-1] + bins[1:]) / 2
    counts, _ = np.histogram(p_vals, bins=bins)
    colors_h = [DM1_C if c < 0.5 else DM2_C for c in centers]
    fig.add_trace(
        go.Bar(x=centers, y=counts, marker_color=colors_h, marker_line_color="black",
               marker_line_width=0.5, width=0.045,
               hovertemplate="p(DM2)≈%{x:.2f}<br>count=%{y}<extra></extra>",
               showlegend=False),
        row=1, col=1,
    )
    fig.add_vline(x=0.5, line_dash="dash", line_color="#444", row=1, col=1)
    fig.add_annotation(x=0.05, y=max(counts) * 0.92,
                       text=f"DM1: {int((p_vals<0.5).sum())}<br>DM2: {int((p_vals>=0.5).sum())}<br>mean={p_vals.mean():.3f}<br>median={np.median(p_vals):.3f}",
                       showarrow=False, align="left",
                       bgcolor="white", bordercolor="#999", borderwidth=1, font=dict(size=11),
                       row=1, col=1)
    fig.update_xaxes(title_text="predicted p(DM2)", row=1, col=1)
    fig.update_yaxes(title_text="Korean samples (count)", row=1, col=1)

    vmax = float(np.nanmax(np.abs(z)))
    fig.add_trace(
        go.Heatmap(z=z, x=g_have, y=rows_label, colorscale="RdBu_r", zmin=-vmax, zmax=vmax,
                   colorbar=dict(title="centered<br>log2(TPM+1)", thickness=12, len=0.85, x=1.02),
                   hovertemplate="%{y}<br>%{x}: %{z:+.2f}<extra></extra>"),
        row=1, col=2,
    )
    fig.add_hline(y=1.5, line_color="black", line_width=1.5, row=1, col=2)

    auc_cv = summary["tcga_5fold_cv_auc_mean"]
    auc_std = summary["tcga_5fold_cv_auc_std"]
    fig.update_layout(
        title=dict(
            text=f"<b>Fig K2 (n={n_kor})</b> — PRJEB11591 Korean PTC validation · within-sample-centered 8-gene profile vs TCGA DM1/DM2 · TCGA 5-fold CV AUC = {auc_cv:.3f} ± {auc_std:.3f}",
            font=dict(size=14), x=0.5, xanchor="center",
        ),
        font=dict(family="-apple-system,Segoe UI,Noto Sans KR,sans-serif", size=11),
        height=720, width=1500, margin=dict(l=80, r=80, t=80, b=60),
        paper_bgcolor="white", plot_bgcolor="white",
    )
    out_html = HTML_DIR / "v17_K2_PRJEB11591_v260.html"
    fig.write_html(out_html, include_plotlyjs="cdn", config={"displaylogo": False, "toImageButtonOptions": {"format": "png", "scale": 2}})
    print(f"  ✓ {out_html}")
    return out_html


def main():
    g_have, Xc_tcga, y_tcga = load_tcga_centered()
    mat, pred, summary = load_korean()
    n_dm1_k, n_dm2_k = build_png(g_have, Xc_tcga, y_tcga, mat, pred, summary)
    out_html = build_html(g_have, Xc_tcga, y_tcga, mat, pred, summary)
    print(f"\nDone — Korean n={len(pred)}, DM1:{n_dm1_k} DM2:{n_dm2_k}")
    print(f"  PNG: {FIG_DIR / 'Korean_K2_v260.png'}")
    print(f"  HTML: {out_html}")


if __name__ == "__main__":
    main()
