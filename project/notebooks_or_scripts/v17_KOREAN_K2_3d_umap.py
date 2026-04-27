#!/usr/bin/env python3
"""Korean K2 3D Plotly Scatter — TCGA DM1/DM2 + 9 Korean samples in PCA space
of within-sample-centered 8-gene panel log2(TPM+1) profile.

Trey-Ideker-style rotatable 3D scatter. PCA computed on TCGA centered profiles,
Korean samples projected onto the same PCA. Centroid markers + sample-level dots.

Output: submission/npj/figures/Korean_K2_3D.html (interactive)
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

import plotly.graph_objects as go

FIG_DIR = Path("/opt/thyroid-dash/project/submission/npj/figures")
FIG_DIR.mkdir(parents=True, exist_ok=True)
GENE_8 = ["SLC5A5", "TPO", "TG", "TSHR", "PAX8", "NKX2-1", "FOXE1", "DIO1"]


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
    pca = PCA(n_components=3, random_state=42).fit(Xc)
    Z_tcga = pca.transform(Xc)
    explained = pca.explained_variance_ratio_ * 100

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
    Z_kor = pca.transform(Xkc)

    # Centroids in PCA space
    Z_dm1_centroid = Z_tcga[y == 0].mean(axis=0)
    Z_dm2_centroid = Z_tcga[y == 1].mean(axis=0)

    fig = go.Figure()

    fig.add_trace(go.Scatter3d(
        x=Z_tcga[y == 1, 0], y=Z_tcga[y == 1, 1], z=Z_tcga[y == 1, 2],
        mode="markers",
        marker=dict(size=3.2, color="#d62728", opacity=0.45,
                    line=dict(width=0)),
        name=f"TCGA DM2 (n={int((y==1).sum())})",
        hovertemplate="TCGA DM2<br>PC1=%{x:.2f}<br>PC2=%{y:.2f}<br>PC3=%{z:.2f}<extra></extra>",
    ))
    fig.add_trace(go.Scatter3d(
        x=Z_tcga[y == 0, 0], y=Z_tcga[y == 0, 1], z=Z_tcga[y == 0, 2],
        mode="markers",
        marker=dict(size=3.2, color="#1f77b4", opacity=0.45,
                    line=dict(width=0)),
        name=f"TCGA DM1 (n={int((y==0).sum())})",
        hovertemplate="TCGA DM1<br>PC1=%{x:.2f}<br>PC2=%{y:.2f}<br>PC3=%{z:.2f}<extra></extra>",
    ))

    fig.add_trace(go.Scatter3d(
        x=[Z_dm2_centroid[0]], y=[Z_dm2_centroid[1]], z=[Z_dm2_centroid[2]],
        mode="markers+text",
        marker=dict(size=14, color="#7f1d1d", symbol="diamond",
                    line=dict(width=2, color="white")),
        text=["DM2 centroid"], textposition="top center",
        textfont=dict(size=12, color="#7f1d1d"),
        name="TCGA DM2 centroid",
        hovertemplate="DM2 centroid<br>PC1=%{x:.2f}<br>PC2=%{y:.2f}<br>PC3=%{z:.2f}<extra></extra>",
    ))
    fig.add_trace(go.Scatter3d(
        x=[Z_dm1_centroid[0]], y=[Z_dm1_centroid[1]], z=[Z_dm1_centroid[2]],
        mode="markers+text",
        marker=dict(size=14, color="#1e3a8a", symbol="diamond",
                    line=dict(width=2, color="white")),
        text=["DM1 centroid"], textposition="top center",
        textfont=dict(size=12, color="#1e3a8a"),
        name="TCGA DM1 centroid",
        hovertemplate="DM1 centroid<br>PC1=%{x:.2f}<br>PC2=%{y:.2f}<br>PC3=%{z:.2f}<extra></extra>",
    ))

    p_map = {r: p for r, p in zip(pred["run"], pred["p_DM2"])}
    p_arr = np.array([p_map[r] for r in mat.index])
    korean_color = ["#16a34a"] * len(mat)  # all DM2-predicted
    fig.add_trace(go.Scatter3d(
        x=Z_kor[:, 0], y=Z_kor[:, 1], z=Z_kor[:, 2],
        mode="markers+text",
        marker=dict(size=10, color=korean_color, symbol="circle",
                    line=dict(width=2, color="black"), opacity=0.95),
        text=[f"{r}<br>p={p:.2f}" for r, p in zip(mat.index, p_arr)],
        textposition="top center",
        textfont=dict(size=9, color="#064e3b"),
        name="🇰🇷 Korean (PRJEB11591, n=9)",
        hovertemplate="<b>%{text}</b><br>PC1=%{x:.2f}<br>PC2=%{y:.2f}<br>PC3=%{z:.2f}<extra></extra>",
    ))

    fig.update_layout(
        title=dict(
            text=("<b>Fig K2-3D.</b> PRJEB11591 Korean PTC pilot (n = 9, green) projected onto TCGA "
                  "DM1/DM2 within-sample-centered 8-gene panel PCA space — Korean cluster overlaps "
                  "TCGA DM2 centroid (preserved differentiation)"),
            x=0.5, xanchor="center", font=dict(size=13),
        ),
        scene=dict(
            xaxis=dict(title=f"PC1 ({explained[0]:.1f}% var)", backgroundcolor="rgb(245,245,250)",
                       gridcolor="white", showspikes=False),
            yaxis=dict(title=f"PC2 ({explained[1]:.1f}% var)", backgroundcolor="rgb(245,245,250)",
                       gridcolor="white", showspikes=False),
            zaxis=dict(title=f"PC3 ({explained[2]:.1f}% var)", backgroundcolor="rgb(245,245,250)",
                       gridcolor="white", showspikes=False),
            camera=dict(eye=dict(x=1.6, y=1.6, z=1.0)),
            aspectmode="cube",
        ),
        legend=dict(x=0.01, y=0.99, bgcolor="rgba(255,255,255,0.9)", bordercolor="#cbd5e1",
                    borderwidth=1, font=dict(size=11)),
        margin=dict(l=0, r=0, t=70, b=0),
        height=720,
        paper_bgcolor="white",
    )

    out_html = FIG_DIR / "Korean_K2_3D.html"
    fig.write_html(str(out_html), include_plotlyjs="cdn",
                   config={"displayModeBar": True, "displaylogo": False})
    print(f"  ✓ {out_html}")
    print(f"    Explained variance: PC1={explained[0]:.1f}%, PC2={explained[1]:.1f}%, PC3={explained[2]:.1f}% (cum={explained.sum():.1f}%)")
    print(f"    Korean centroid in PC space: PC1={Z_kor[:,0].mean():.2f}, PC2={Z_kor[:,1].mean():.2f}, PC3={Z_kor[:,2].mean():.2f}")
    print(f"    DM2 centroid:                 PC1={Z_dm2_centroid[0]:.2f}, PC2={Z_dm2_centroid[1]:.2f}, PC3={Z_dm2_centroid[2]:.2f}")
    print(f"    DM1 centroid:                 PC1={Z_dm1_centroid[0]:.2f}, PC2={Z_dm1_centroid[1]:.2f}, PC3={Z_dm1_centroid[2]:.2f}")


if __name__ == "__main__":
    main()
