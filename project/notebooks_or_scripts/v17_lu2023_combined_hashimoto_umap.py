#!/usr/bin/env python3
"""v17 — Combined Lu 2023 (GSE193581) UMAP + Hashimoto-like overlay.

Reads:
  - results/v17_lu2023/GSE193581_hvg_adata.h5ad         (Q14 UMAP + DM_score)
  - results/v17_lu2023_hashimoto/GSE193581_cells_DM_HLA.tsv (Q15 hashi_like flag)

Produces 4 plotly_dark interactive figures under
  reports/html/figs_interactive/v17/v17_lu2023_combined_*.html
plus a 2x2 dashboard.
"""
from __future__ import annotations
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import scanpy as sc
from plotly.subplots import make_subplots

H5AD = Path("/opt/thyroid-dash/project/results/v17_lu2023/GSE193581_hvg_adata.h5ad")
HTSV = Path("/opt/thyroid-dash/project/results/v17_lu2023_hashimoto/GSE193581_cells_DM_HLA.tsv")
FIG = Path("/opt/thyroid-dash/project/reports/html/figs_interactive/v17")
FIG.mkdir(parents=True, exist_ok=True)
LOG = Path("/opt/thyroid-dash/project/logs/v17_lu2023_combined_hashimoto_umap.log")
LOG.parent.mkdir(parents=True, exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.FileHandler(LOG, mode="w"), logging.StreamHandler(sys.stdout)],
)
log = logging.getLogger(__name__)

SEED = 42
np.random.seed(SEED)

BG, INK = "#0b0e12", "#F2F2F2"
DARK = dict(template="plotly_dark", paper_bgcolor=BG, plot_bgcolor=BG,
            font=dict(color=INK, size=13))
COLOR_HISTO = {"NORM": "#2ECC71", "PTC": "#F5A623", "ATC": "#c24c4c"}
TEAL = "#7ccfcd"
GRAY = "rgba(150,150,150,0.18)"
HASHI_COLOR = TEAL
OTHER_PTC_COLOR = "rgba(220,220,220,0.45)"

# ---------------- LOAD ----------------
log.info("Loading h5ad …")
adata = sc.read_h5ad(H5AD)
log.info(f"adata: {adata.shape}, obsm: {list(adata.obsm.keys())}")

um = pd.DataFrame(adata.obsm["X_umap"], index=adata.obs_names, columns=["UMAP1", "UMAP2"])
um["sample"] = adata.obs["sample"].astype(str).values
um["histology"] = adata.obs["histology"].astype(str).values
um["author_celltype"] = adata.obs["author_celltype"].astype(str).values
um["DM_score"] = adata.obs["DM_score"].values

log.info("Loading hashimoto TSV …")
hashi = pd.read_csv(HTSV, sep="\t", index_col=0)
hashi.index.name = "cell"
log.info(f"hashi: {hashi.shape}, hashi_like True={int(hashi['hashi_like'].sum())}")

# Merge
um["HLA_I_score"] = hashi["HLA_I_score"].reindex(um.index)
um["HLA_II_score"] = hashi["HLA_II_score"].reindex(um.index)
um["DM_class"] = hashi["DM_class"].reindex(um.index)
um["hashi_like"] = hashi["hashi_like"].reindex(um.index)
# overwrite DM_score with the hashi version where present (consistent with Q15)
um.loc[hashi.index.intersection(um.index), "DM_score"] = hashi.loc[
    hashi.index.intersection(um.index), "DM_score"
].values

n_hashi = int(um["hashi_like"].fillna(False).sum())
n_ptc_mal = int(um["hashi_like"].notna().sum())
log.info(f"on UMAP: hashi_like cells={n_hashi}, PTC malignant cells={n_ptc_mal}")


def _hover_text(df: pd.DataFrame) -> list[str]:
    out = []
    for r in df.itertuples():
        dm = r.DM_score
        h2 = r.HLA_II_score
        dm_str = f"{dm:+.2f}" if pd.notna(dm) else "—"
        h2_str = f"{h2:+.2f}" if pd.notna(h2) else "—"
        out.append(
            f"<b>{r.Index}</b><br>sample: {r.sample}<br>"
            f"histology: {r.histology}<br>celltype: {r.author_celltype}<br>"
            f"DM_score: {dm_str}<br>HLA-II z: {h2_str}"
        )
    return out


# =========================================================
# FIG 1 — UMAP with Hashimoto-like cells highlighted
# =========================================================
log.info("Fig 1 — UMAP × hashimoto overlay")
mask_ptc_mal = um["hashi_like"].notna()
mask_hashi = um["hashi_like"].fillna(False).astype(bool)
mask_other_ptc_mal = mask_ptc_mal & ~mask_hashi
mask_other = ~mask_ptc_mal

fig = go.Figure()
sub = um[mask_other]
fig.add_trace(go.Scattergl(
    x=sub["UMAP1"], y=sub["UMAP2"], mode="markers",
    marker=dict(size=2, color=GRAY),
    name=f"non-PTC-malignant (n={len(sub):,})",
    hoverinfo="skip",
))
sub = um[mask_other_ptc_mal]
fig.add_trace(go.Scattergl(
    x=sub["UMAP1"], y=sub["UMAP2"], mode="markers",
    marker=dict(size=3.5, color=OTHER_PTC_COLOR, line=dict(width=0)),
    name=f"PTC malignant — other (n={len(sub):,})",
    text=_hover_text(sub), hoverinfo="text",
))
sub = um[mask_hashi]
fig.add_trace(go.Scattergl(
    x=sub["UMAP1"], y=sub["UMAP2"], mode="markers",
    marker=dict(size=7, color=HASHI_COLOR,
                line=dict(width=1, color="#0b0e12"), opacity=0.95),
    name=f"PTC Hashimoto-like (n={len(sub):,})",
    text=_hover_text(sub), hoverinfo="text",
))
fig.update_layout(
    title=dict(text=(
        "<b>Lu 2023 GSE193581 — UMAP × Hashimoto-like overlay</b><br>"
        f"<sub style='color:#7ccfcd'>teal = PTC malignant cell with HLA-II z>1 "
        f"(n={int(mask_hashi.sum()):,}); white = other PTC malignant "
        f"(n={int(mask_other_ptc_mal.sum()):,}); gray = other cells</sub>"
    )),
    xaxis=dict(title="UMAP1", showgrid=False, zeroline=False),
    yaxis=dict(title="UMAP2", showgrid=False, zeroline=False),
    height=640, width=980, margin=dict(l=60, r=20, t=90, b=60), **DARK,
    legend=dict(itemsizing="constant", bgcolor="rgba(0,0,0,0.3)"),
)
out1 = FIG / "v17_lu2023_combined_umap_hashimoto.html"
fig.write_html(out1, include_plotlyjs="cdn", full_html=True)
log.info(f"wrote {out1}")


# =========================================================
# FIG 2 — UMAP colored by HLA-II score (malignant only)
# =========================================================
log.info("Fig 2 — UMAP × HLA-II score gradient")
um_mal = um[mask_ptc_mal]
um_mal_atc = um[(um["histology"] == "ATC") & (um["author_celltype"] == "Malignant cell")]
um_other = um[~mask_ptc_mal]

fig = go.Figure()
fig.add_trace(go.Scattergl(
    x=um_other["UMAP1"], y=um_other["UMAP2"], mode="markers",
    marker=dict(size=2, color=GRAY),
    name=f"non-PTC-malignant (n={len(um_other):,})",
    hoverinfo="skip",
))
fig.add_trace(go.Scattergl(
    x=um_mal["UMAP1"], y=um_mal["UMAP2"], mode="markers",
    marker=dict(
        size=4,
        color=um_mal["HLA_II_score"],
        colorscale="Viridis",
        cmin=float(np.nanpercentile(um_mal["HLA_II_score"], 2)),
        cmax=float(np.nanpercentile(um_mal["HLA_II_score"], 98)),
        showscale=True,
        colorbar=dict(
            title=dict(text="HLA-II z", font=dict(color=INK)),
            tickfont=dict(color=INK), thickness=14, len=0.7,
        ),
        line=dict(width=0),
    ),
    name=f"PTC malignant (n={len(um_mal):,})",
    text=_hover_text(um_mal), hoverinfo="text",
))
fig.update_layout(
    title=dict(text=(
        "<b>Lu 2023 GSE193581 — UMAP colored by HLA-II score (PTC malignant cells only)</b><br>"
        "<sub style='color:#7ccfcd'>HLA-II = mean z of HLA-DRA/DRB1/DPA1/DPB1/DQA1/DQB1/CIITA; "
        f"hashi-like cutoff z>1 ⇒ n={int(mask_hashi.sum()):,} cells</sub>"
    )),
    xaxis=dict(title="UMAP1", showgrid=False, zeroline=False),
    yaxis=dict(title="UMAP2", showgrid=False, zeroline=False),
    height=640, width=980, margin=dict(l=60, r=20, t=90, b=60), **DARK,
    legend=dict(itemsizing="constant", bgcolor="rgba(0,0,0,0.3)"),
)
out2 = FIG / "v17_lu2023_combined_umap_hlaII.html"
fig.write_html(out2, include_plotlyjs="cdn", full_html=True)
log.info(f"wrote {out2}")


# =========================================================
# FIG 3 — Zoom into PTC malignant cluster, colored by DM_score with hashi rings
# =========================================================
log.info("Fig 3 — UMAP zoom × DM_score with hashi rings")
# bbox = PTC-malignant cell bounding box w/ small padding
xy = um_mal[["UMAP1", "UMAP2"]].dropna()
pad_x = 0.05 * (xy["UMAP1"].max() - xy["UMAP1"].min())
pad_y = 0.05 * (xy["UMAP2"].max() - xy["UMAP2"].min())
xlim = (xy["UMAP1"].min() - pad_x, xy["UMAP1"].max() + pad_x)
ylim = (xy["UMAP2"].min() - pad_y, xy["UMAP2"].max() + pad_y)

# context cells inside bbox (gray)
ctx = um_other[
    (um_other["UMAP1"] >= xlim[0]) & (um_other["UMAP1"] <= xlim[1]) &
    (um_other["UMAP2"] >= ylim[0]) & (um_other["UMAP2"] <= ylim[1])
]
fig = go.Figure()
fig.add_trace(go.Scattergl(
    x=ctx["UMAP1"], y=ctx["UMAP2"], mode="markers",
    marker=dict(size=2, color=GRAY),
    name=f"context — non-PTC-malignant (n={len(ctx):,})",
    hoverinfo="skip",
))
# all PTC malignant colored by DM_score
fig.add_trace(go.Scattergl(
    x=um_mal["UMAP1"], y=um_mal["UMAP2"], mode="markers",
    marker=dict(
        size=5,
        color=um_mal["DM_score"],
        colorscale="RdBu",
        cmid=0,
        cmin=float(np.nanpercentile(um_mal["DM_score"], 2)),
        cmax=float(np.nanpercentile(um_mal["DM_score"], 98)),
        showscale=True,
        colorbar=dict(
            title=dict(text="DM_score", font=dict(color=INK)),
            tickfont=dict(color=INK), thickness=14, len=0.7, x=1.02,
        ),
        line=dict(width=0),
    ),
    name=f"PTC malignant (n={len(um_mal):,})",
    text=_hover_text(um_mal), hoverinfo="text",
))
# Hashi-like rings overlay
hashi_sub = um[mask_hashi]
fig.add_trace(go.Scattergl(
    x=hashi_sub["UMAP1"], y=hashi_sub["UMAP2"], mode="markers",
    marker=dict(size=10, color="rgba(0,0,0,0)",
                line=dict(width=1.6, color=HASHI_COLOR), opacity=0.95),
    name=f"Hashimoto-like ring (n={len(hashi_sub):,})",
    text=_hover_text(hashi_sub), hoverinfo="text",
))
fig.update_layout(
    title=dict(text=(
        "<b>Lu 2023 GSE193581 — PTC malignant zoom · DM_score (heat) + Hashimoto rings</b><br>"
        "<sub style='color:#7ccfcd'>Hashimoto-like teal rings sit at slightly LOWER DM_score within "
        "the PTC malignant cluster (median DM 0.192 vs 0.298, MW p=6.4e-12)</sub>"
    )),
    xaxis=dict(title="UMAP1", range=xlim, showgrid=False, zeroline=False),
    yaxis=dict(title="UMAP2", range=ylim, showgrid=False, zeroline=False),
    height=680, width=1000, margin=dict(l=60, r=20, t=90, b=60), **DARK,
    legend=dict(itemsizing="constant", bgcolor="rgba(0,0,0,0.3)"),
)
out3 = FIG / "v17_lu2023_combined_umap_zoom_dmhashi.html"
fig.write_html(out3, include_plotlyjs="cdn", full_html=True)
log.info(f"wrote {out3}")


# =========================================================
# FIG 4 — Combined dashboard (2×2)
# =========================================================
log.info("Fig 4 — Combined 2x2 dashboard")
fig = make_subplots(
    rows=2, cols=2,
    subplot_titles=(
        "(1) UMAP × histology",
        "(2) UMAP × Hashimoto-like overlay",
        "(3) DM_score × HLA-II (PTC malignant)",
        "(4) Per-PTC-sample Hashimoto-like %",
    ),
    horizontal_spacing=0.09, vertical_spacing=0.13,
    specs=[[{"type": "scattergl"}, {"type": "scattergl"}],
           [{"type": "scattergl"}, {"type": "bar"}]],
)

# (1) UMAP × histology
for h in ["NORM", "PTC", "ATC"]:
    sub = um[um["histology"] == h]
    fig.add_trace(go.Scattergl(
        x=sub["UMAP1"], y=sub["UMAP2"], mode="markers",
        marker=dict(size=2.2, color=COLOR_HISTO[h], opacity=0.55),
        name=f"{h} (n={len(sub):,})",
        legendgroup="histology", hoverinfo="skip",
        showlegend=True,
    ), row=1, col=1)

# (2) UMAP × hashimoto
sub = um[mask_other]
fig.add_trace(go.Scattergl(
    x=sub["UMAP1"], y=sub["UMAP2"], mode="markers",
    marker=dict(size=2, color=GRAY),
    name="non-PTC-mal", legendgroup="hashi",
    hoverinfo="skip", showlegend=True,
), row=1, col=2)
sub = um[mask_other_ptc_mal]
fig.add_trace(go.Scattergl(
    x=sub["UMAP1"], y=sub["UMAP2"], mode="markers",
    marker=dict(size=3, color=OTHER_PTC_COLOR),
    name=f"PTC mal — other ({len(sub):,})",
    legendgroup="hashi", hoverinfo="skip", showlegend=True,
), row=1, col=2)
sub = um[mask_hashi]
fig.add_trace(go.Scattergl(
    x=sub["UMAP1"], y=sub["UMAP2"], mode="markers",
    marker=dict(size=6, color=HASHI_COLOR,
                line=dict(width=1, color="#0b0e12")),
    name=f"Hashi-like ({len(sub):,})",
    legendgroup="hashi", hoverinfo="skip", showlegend=True,
), row=1, col=2)

# (3) DM_score × HLA-II 2D
ptc_mal_df = um[mask_ptc_mal].dropna(subset=["DM_score", "HLA_II_score"])
fig.add_trace(go.Scattergl(
    x=ptc_mal_df.loc[~ptc_mal_df["hashi_like"].astype(bool), "DM_score"],
    y=ptc_mal_df.loc[~ptc_mal_df["hashi_like"].astype(bool), "HLA_II_score"],
    mode="markers",
    marker=dict(size=3.5, color=OTHER_PTC_COLOR),
    name=f"PTC mal — other ({int((~ptc_mal_df['hashi_like'].astype(bool)).sum()):,})",
    legendgroup="scatter", hoverinfo="skip", showlegend=True,
), row=2, col=1)
fig.add_trace(go.Scattergl(
    x=ptc_mal_df.loc[ptc_mal_df["hashi_like"].astype(bool), "DM_score"],
    y=ptc_mal_df.loc[ptc_mal_df["hashi_like"].astype(bool), "HLA_II_score"],
    mode="markers",
    marker=dict(size=5, color=HASHI_COLOR,
                line=dict(width=0.6, color="#0b0e12"), opacity=0.95),
    name=f"Hashi-like ({int(ptc_mal_df['hashi_like'].astype(bool).sum()):,})",
    legendgroup="scatter", hoverinfo="skip", showlegend=True,
), row=2, col=1)
fig.add_hline(y=1.0, line=dict(color="#FFFFFF", dash="dot", width=1),
              annotation_text="HLA-II z=1 (Hashi cutoff)",
              annotation_font=dict(color="#bbb", size=10),
              annotation_position="top right",
              row=2, col=1)

# (4) Per-PTC-sample Hashimoto-like %
ptc_only = hashi[hashi["histology"] == "PTC"]
per_sample = (
    ptc_only.groupby("sample")
    .agg(n_cells=("hashi_like", "size"), n_hashi=("hashi_like", "sum"))
)
per_sample["pct_hashi"] = per_sample["n_hashi"] / per_sample["n_cells"] * 100.0
per_sample = per_sample.sort_values("pct_hashi", ascending=False)
bar_colors = [HASHI_COLOR if v >= 20 else "rgba(150,150,150,0.65)"
              for v in per_sample["pct_hashi"]]
fig.add_trace(go.Bar(
    x=per_sample.index.tolist(),
    y=per_sample["pct_hashi"].values,
    marker=dict(color=bar_colors, line=dict(color="#0b0e12", width=0.5)),
    text=[f"{v:.1f}%" for v in per_sample["pct_hashi"]],
    textposition="outside",
    textfont=dict(color=INK, size=11),
    name="Hashi-like %",
    showlegend=False,
    customdata=np.stack([per_sample["n_cells"].values,
                         per_sample["n_hashi"].values], axis=-1),
    hovertemplate=("sample=%{x}<br>%{y:.1f}% Hashi-like<br>"
                   "n_cells=%{customdata[0]}<br>n_hashi=%{customdata[1]}"
                   "<extra></extra>"),
), row=2, col=2)
# 20% threshold line
fig.add_hline(y=20.0, line=dict(color="#FFFFFF", dash="dot", width=1),
              annotation_text="20% threshold",
              annotation_font=dict(color="#bbb", size=10),
              annotation_position="top left",
              row=2, col=2)

fig.update_xaxes(title_text="UMAP1", showgrid=False, zeroline=False, row=1, col=1)
fig.update_yaxes(title_text="UMAP2", showgrid=False, zeroline=False, row=1, col=1)
fig.update_xaxes(title_text="UMAP1", showgrid=False, zeroline=False, row=1, col=2)
fig.update_yaxes(title_text="UMAP2", showgrid=False, zeroline=False, row=1, col=2)
fig.update_xaxes(title_text="DM_score", showgrid=False, zeroline=False, row=2, col=1)
fig.update_yaxes(title_text="HLA-II z", showgrid=False, zeroline=False, row=2, col=1)
fig.update_xaxes(title_text="PTC sample", showgrid=False, zeroline=False, row=2, col=2)
fig.update_yaxes(title_text="% Hashimoto-like cells",
                 showgrid=False, zeroline=False, row=2, col=2)

fig.update_layout(
    title=dict(text=(
        "<b>Lu 2023 GSE193581 — Hashimoto-like single-cell story (combined dashboard)</b><br>"
        f"<sub style='color:#7ccfcd'>67,678 cells / 23 samples; PTC malignant {n_ptc_mal:,} cells, "
        f"Hashimoto-like {n_hashi:,} (HLA-II z>1); 1/6 PTC samples >20% threshold ≈ Q12 Korean 17%</sub>"
    )),
    height=900, width=1400, **DARK,
    margin=dict(l=60, r=30, t=110, b=60),
    legend=dict(itemsizing="constant", bgcolor="rgba(0,0,0,0.3)",
                orientation="v", x=1.02, y=1.0),
)
out4 = FIG / "v17_lu2023_combined_dashboard.html"
fig.write_html(out4, include_plotlyjs="cdn", full_html=True)
log.info(f"wrote {out4}")

log.info("=== combined hashimoto UMAP figures DONE ===")
print("\n".join(str(p) for p in [out1, out2, out3, out4]))
