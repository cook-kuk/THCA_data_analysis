#!/usr/bin/env python3
"""K2 Korean figure (matplotlib PNG, embeddable in dashboard).
  a) 8-gene within-sample-centered profile heatmap: TCGA DM2 + TCGA DM1 + 9 Korean
  b) p_DM2 horizontal bar chart per Korean sample
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

FIG_DIR = Path("/opt/thyroid-dash/project/submission/npj/figures")
FIG_DIR.mkdir(parents=True, exist_ok=True)
GENE_8 = ["SLC5A5", "TPO", "TG", "TSHR", "PAX8", "NKX2-1", "FOXE1", "DIO1"]
DM1_COLOR = "#1f77b4"
DM2_COLOR = "#d62728"


def main():
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
    dm1_centroid = Xc[y == 0].mean(0)
    dm2_centroid = Xc[y == 1].mean(0)
    n_dm1, n_dm2 = int((y == 0).sum()), int((y == 1).sum())

    mat = pd.read_csv(
        "/opt/thyroid-dash/project/results/v17_korean/K2_8gene_tpm_matrix_v4.tsv",
        sep="\t", index_col=0,
    )
    pred = pd.read_csv(
        "/opt/thyroid-dash/project/results/v17_korean/K2_korean_predictions_v4.tsv",
        sep="\t",
    )
    Xk = np.log2(mat[g_have].values + 1.0)
    Xkc = Xk - Xk.mean(axis=1, keepdims=True)

    pred_sorted = pred.sort_values("p_DM2", ascending=False).reset_index(drop=True)
    korean_order = list(pred_sorted["run"])
    korean_idx = [list(mat.index).index(r) for r in korean_order]

    rows_label = (
        [f"TCGA DM2 centroid (n={n_dm2})", f"TCGA DM1 centroid (n={n_dm1})"]
        + [f"{r}  (p={p:.2f})" for r, p in zip(pred_sorted["run"], pred_sorted["p_DM2"])]
    )
    z = np.vstack([dm2_centroid, dm1_centroid, Xkc[korean_idx]])

    mpl.rcParams.update({
        "font.family": ["DejaVu Sans", "Arial Unicode MS", "WenQuanYi Zen Hei"],
        "font.size": 10,
        "axes.titlesize": 11,
        "axes.labelsize": 10,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
    })

    fig = plt.figure(figsize=(15.5, 6.0), dpi=120)
    gs = fig.add_gridspec(1, 2, width_ratios=[1.55, 1.0], wspace=0.32)

    ax0 = fig.add_subplot(gs[0])
    vmax = float(np.nanmax(np.abs(z)))
    im = ax0.imshow(z, cmap="RdBu_r", aspect="auto", vmin=-vmax, vmax=vmax)
    ax0.set_xticks(range(len(g_have)))
    ax0.set_xticklabels(g_have, rotation=20, ha="right")
    ax0.set_yticks(range(len(rows_label)))
    ax0.set_yticklabels(rows_label)
    ax0.axhline(1.5, color="black", linewidth=1.4)
    ax0.set_title("a  8-gene within-sample-centered profile (TCGA centroids vs 9 Korean samples)", loc="left", fontweight="bold")
    cbar = fig.colorbar(im, ax=ax0, fraction=0.04, pad=0.018)
    cbar.set_label("centered log2(TPM+1)", fontsize=9)
    for i in range(z.shape[0]):
        for j in range(z.shape[1]):
            ax0.text(j, i, f"{z[i, j]:+.1f}", ha="center", va="center",
                     fontsize=7.5, color="black" if abs(z[i, j]) < vmax * 0.55 else "white")

    ax1 = fig.add_subplot(gs[1])
    bar_colors = [DM2_COLOR if p >= 0.5 else DM1_COLOR for p in pred_sorted["p_DM2"]]
    y_pos = np.arange(len(pred_sorted))
    ax1.barh(y_pos, pred_sorted["p_DM2"], color=bar_colors, edgecolor="black", linewidth=0.5)
    ax1.set_yticks(y_pos)
    ax1.set_yticklabels(pred_sorted["run"])
    ax1.invert_yaxis()
    ax1.axvline(0.5, color="#555", linestyle="--", linewidth=1.0)
    ax1.text(0.51, len(pred_sorted) - 0.5, "DM1 / DM2 cut", fontsize=8, color="#555", va="center")
    ax1.set_xlim(0, 1.05)
    ax1.set_xlabel("p(DM2) — TCGA-trained scale-invariant LogReg")
    ax1.set_title("b  Korean p(DM2) per sample (sorted)", loc="left", fontweight="bold")
    for i, p in enumerate(pred_sorted["p_DM2"]):
        ax1.text(p + 0.012, i, f"{p:.3f}", va="center", fontsize=8.5)
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)

    fig.suptitle(
        "Fig K2.  PRJEB11591 Korean PTC pilot (n = 9) — within-sample-centered 8-gene profile vs TCGA DM1/DM2 centroids;  9/9 → DM2 (mean p = 0.898)",
        fontsize=12.5, fontweight="bold", y=0.998,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.965))

    out_png = FIG_DIR / "Korean_K2.png"
    out_pdf = FIG_DIR / "Korean_K2.pdf"
    fig.savefig(out_png, dpi=180, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✓ {out_png}")
    print(f"  ✓ {out_pdf}")


if __name__ == "__main__":
    main()
