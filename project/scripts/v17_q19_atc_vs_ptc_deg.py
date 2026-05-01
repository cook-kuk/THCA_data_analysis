"""v17 Q19: ATC vs PTC malignant-cell DEG (Lu 2023 GSE193581).

Identifies non-panel dedifferentiation markers driving ATC transformation.
Author: Seungho Cook
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import scanpy as sc

PROJ = Path("/opt/thyroid-dash/project")
H5AD = PROJ / "results/v17_lu2023/GSE193581_hvg_adata.h5ad"
PANEL_TSV = PROJ / "results/v17_lu2023/GSE193581_cell_panel_score.tsv"
OUT = PROJ / "results/v17_lu2023_DEG"
OUT.mkdir(parents=True, exist_ok=True)

BG = "#0b0e12"
PANEL8 = ["TG", "TPO", "SLC5A5", "TSHR", "PAX8", "FOXE1", "NKX2-1", "DIO2"]
# common alt symbols in some panels
PANEL8_ALT = {"NKX2-1": ["NKX2-1", "TTF1", "TITF1"], "SLC5A5": ["SLC5A5", "NIS"]}

# curated novel-marker watchlist (literature-driven)
WATCHLIST = {
    "EMT/ECM":    ["VIM", "ZEB1", "ZEB2", "SNAI1", "SNAI2", "TWIST1", "FN1", "CDH2", "S100A4", "MMP9", "MMP2"],
    "Immune evasion": ["CD274", "PDCD1LG2", "CTLA4", "LGALS9", "IDO1"],
    "Stemness":   ["POU5F1", "NANOG", "SOX2", "KLF4", "ALDH1A3", "PROM1", "CD44"],
    "ATC-known":  ["CDKN2A", "MKI67", "TOP2A", "EZH2", "BIRC5", "AURKB", "CCNB1", "CDK1"],
    "Stress/inflam": ["S100A8", "S100A9", "LCN2", "SPP1", "IL6", "IL1B"],
}
WATCH_FLAT = {g: cat for cat, gs in WATCHLIST.items() for g in gs}


def main() -> None:
    print("[load] h5ad")
    adata = sc.read_h5ad(H5AD)
    panel_df = pd.read_csv(PANEL_TSV, sep="\t")
    malignant_cells = set(panel_df["cell"])

    # restrict to malignant ATC + PTC
    mask = adata.obs_names.isin(malignant_cells) & adata.obs["histology"].isin(["ATC", "PTC", "NORM"])
    adata_all = adata[mask].copy()
    print("[filter] malignant subset:", adata_all.shape, adata_all.obs["histology"].value_counts().to_dict())

    # --- DEG: ATC vs PTC (malignant only) ---
    adata_pa = adata_all[adata_all.obs["histology"].isin(["ATC", "PTC"])].copy()
    adata_pa.obs["histology"] = adata_pa.obs["histology"].astype("category")
    # NOTE: adata.X is z-scored/scaled (range ~ -2.8 .. 10), so scanpy's log2FC=NaN.
    # We replace log2FC with a "delta mean (scaled)" effect size: mean(ATC) - mean(PTC) on scaled X.
    # This is monotone-equivalent to a standardised mean diff, scientifically interpretable as
    # an effect size in z-units across HVGs.
    is_atc = adata_pa.obs["histology"].values == "ATC"
    is_ptc = adata_pa.obs["histology"].values == "PTC"
    Xpa = np.asarray(adata_pa.X)
    mean_atc = Xpa[is_atc].mean(axis=0)
    mean_ptc = Xpa[is_ptc].mean(axis=0)
    delta_atc_minus_ptc = mean_atc - mean_ptc
    delta_df = pd.DataFrame({"gene": adata_pa.var_names,
                             "delta_mean_scaled": delta_atc_minus_ptc,
                             "mean_ATC_scaled": mean_atc,
                             "mean_PTC_scaled": mean_ptc})

    print("[deg] rank_genes_groups ATC vs PTC, n_genes=200")
    sc.tl.rank_genes_groups(
        adata_pa, groupby="histology", groups=["ATC"], reference="PTC",
        method="wilcoxon", n_genes=200, use_raw=False,
    )
    atc_up = sc.get.rank_genes_groups_df(adata_pa, group="ATC")
    atc_up = atc_up.rename(columns={"names": "gene", "logfoldchanges": "log2FC_invalid",
                                    "pvals": "pval", "pvals_adj": "padj", "scores": "wilcoxon_z"})
    atc_up = atc_up.drop(columns=["log2FC_invalid"]).merge(delta_df, on="gene", how="left")

    # reverse: PTC vs ATC
    print("[deg] rank_genes_groups PTC vs ATC, n_genes=200")
    sc.tl.rank_genes_groups(
        adata_pa, groupby="histology", groups=["PTC"], reference="ATC",
        method="wilcoxon", n_genes=200, use_raw=False,
    )
    ptc_up = sc.get.rank_genes_groups_df(adata_pa, group="PTC")
    ptc_up = ptc_up.rename(columns={"names": "gene", "logfoldchanges": "log2FC_invalid",
                                    "pvals": "pval", "pvals_adj": "padj", "scores": "wilcoxon_z"})
    ptc_up = ptc_up.drop(columns=["log2FC_invalid"]).merge(delta_df, on="gene", how="left")

    # pct expressing per group (use_raw=False; "expressing" = X > 0)
    print("[pct] computing pct_in_ATC / pct_in_PTC")
    Xa = adata_pa.X
    if hasattr(Xa, "toarray"):
        # sparse handling
        idx_atc = np.where(adata_pa.obs["histology"].values == "ATC")[0]
        idx_ptc = np.where(adata_pa.obs["histology"].values == "PTC")[0]
        pct_atc = np.asarray((Xa[idx_atc] > 0).mean(axis=0)).ravel()
        pct_ptc = np.asarray((Xa[idx_ptc] > 0).mean(axis=0)).ravel()
    else:
        is_atc = adata_pa.obs["histology"].values == "ATC"
        is_ptc = adata_pa.obs["histology"].values == "PTC"
        pct_atc = (Xa[is_atc] > 0).mean(axis=0)
        pct_ptc = (Xa[is_ptc] > 0).mean(axis=0)
    pct_df = pd.DataFrame({"gene": adata_pa.var_names, "pct_in_ATC": pct_atc, "pct_in_PTC": pct_ptc})

    atc_up_full = atc_up.merge(pct_df, on="gene", how="left")
    ptc_up_full = ptc_up.merge(pct_df, on="gene", how="left")

    # save unified top200 ATC-up
    atc_up_full.to_csv(OUT / "atc_vs_ptc_top200.tsv", sep="\t", index=False)
    ptc_up_full.to_csv(OUT / "ptc_vs_atc_top200.tsv", sep="\t", index=False)
    print("[save] atc_vs_ptc_top200.tsv, ptc_vs_atc_top200.tsv")

    # --- panel rank within "PTC > ATC" (= differentiation) list ---
    panel_resolved = []
    for sym in PANEL8:
        for alt in PANEL8_ALT.get(sym, [sym]):
            if alt in adata_pa.var_names:
                panel_resolved.append((sym, alt))
                break
        else:
            panel_resolved.append((sym, None))
    print("[panel] resolved:", panel_resolved)

    # full ranking PTC vs ATC (all genes) for panel rank
    sc.tl.rank_genes_groups(
        adata_pa, groupby="histology", groups=["PTC"], reference="ATC",
        method="wilcoxon", n_genes=adata_pa.n_vars, use_raw=False,
    )
    ptc_full_rank = sc.get.rank_genes_groups_df(adata_pa, group="PTC")
    ptc_full_rank["rank"] = np.arange(1, len(ptc_full_rank) + 1)
    panel_rank = []
    for sym, present in panel_resolved:
        if present is None:
            panel_rank.append({"panel_gene": sym, "resolved_symbol": None, "rank_PTC_up": None,
                               "delta_mean_scaled": None, "padj": None, "in_HVG": False})
            continue
        row = ptc_full_rank[ptc_full_rank["names"] == present]
        delta_row = delta_df[delta_df["gene"] == present]
        if len(row) == 0:
            panel_rank.append({"panel_gene": sym, "resolved_symbol": present, "rank_PTC_up": None,
                               "delta_mean_scaled": None, "padj": None, "in_HVG": True})
        else:
            r = row.iloc[0]
            panel_rank.append({
                "panel_gene": sym, "resolved_symbol": present,
                "rank_PTC_up": int(r["rank"]),
                "delta_mean_scaled": float(delta_row["delta_mean_scaled"].iloc[0]) if len(delta_row) else None,
                "padj": float(r["pvals_adj"]),
                "in_HVG": True,
            })
    panel_rank_df = pd.DataFrame(panel_rank)
    print("[panel rank in PTC>ATC]")
    print(panel_rank_df)

    # --- novel marker candidates: top ATC-up not in panel, flagged with category ---
    panel_set = {p[1] for p in panel_resolved if p[1] is not None}
    novel = atc_up_full.copy()
    novel["category"] = novel["gene"].map(WATCH_FLAT).fillna("other")
    novel["is_panel"] = novel["gene"].isin(panel_set)
    novel = novel[~novel["is_panel"]].head(50)
    novel_high = novel[novel["category"] != "other"].head(20)
    print("[novel] watchlist hits in top200 ATC-up:")
    print(novel_high[["gene", "wilcoxon_z", "delta_mean_scaled", "padj", "pct_in_ATC", "pct_in_PTC", "category"]])

    # --- Volcano plot (use the ATC-vs-PTC ranking on all HVGs) ---
    sc.tl.rank_genes_groups(
        adata_pa, groupby="histology", groups=["ATC"], reference="PTC",
        method="wilcoxon", n_genes=adata_pa.n_vars, use_raw=False,
    )
    atc_full = sc.get.rank_genes_groups_df(adata_pa, group="ATC")
    atc_full = atc_full.rename(columns={"names": "gene", "pvals_adj": "padj"})
    atc_full = atc_full.merge(delta_df, on="gene", how="left")
    atc_full["neglog10_padj"] = -np.log10(np.clip(atc_full["padj"].values, 1e-300, 1.0))
    atc_full["category"] = atc_full["gene"].map(WATCH_FLAT).fillna("other")
    atc_full["is_panel"] = atc_full["gene"].isin(panel_set)

    def label_color(row):
        if row["is_panel"]:
            return "panel-8"
        if row["category"] != "other":
            return row["category"]
        if row["padj"] < 1e-10 and abs(row["delta_mean_scaled"]) > 0.5:
            return "sig"
        return "ns"

    atc_full["color_class"] = atc_full.apply(label_color, axis=1)
    color_map = {
        "panel-8": "#fde047",
        "EMT/ECM": "#f87171",
        "Immune evasion": "#a78bfa",
        "Stemness": "#4ade80",
        "ATC-known": "#fb923c",
        "Stress/inflam": "#22d3ee",
        "other": "#475569",
        "sig": "#94a3b8",
        "ns": "#334155",
    }
    fig = px.scatter(
        atc_full, x="delta_mean_scaled", y="neglog10_padj", color="color_class",
        color_discrete_map=color_map, hover_name="gene", template="plotly_dark",
        title="ATC vs PTC malignant DEG (Lu 2023, n=14,655 cells)",
        labels={"delta_mean_scaled": "Δ mean expression (ATC − PTC, scaled units)",
                "neglog10_padj": "-log10(adj. p)"},
        opacity=0.7,
    )
    # add gene labels for panel + top watchlist
    label_genes = list(panel_set) + novel_high["gene"].tolist()
    for g in set(label_genes):
        sub = atc_full[atc_full["gene"] == g]
        if len(sub):
            fig.add_annotation(x=float(sub["delta_mean_scaled"].iloc[0]),
                               y=float(sub["neglog10_padj"].iloc[0]),
                               text=g, showarrow=True, arrowhead=1, font=dict(color="#e2e8f0", size=10),
                               arrowcolor="#94a3b8")
    fig.update_layout(paper_bgcolor=BG, plot_bgcolor=BG, font_color="#e2e8f0",
                      legend=dict(bgcolor=BG, bordercolor="#334155"))
    fig.write_html(OUT / "v17_q19_volcano_ATC_vs_PTC.html", include_plotlyjs="cdn")

    # --- Top 20 ATC-up dot plot across NORM/PTC/ATC ---
    top20_atc = atc_up_full.head(20)["gene"].tolist()
    groups = ["NORM", "PTC", "ATC"]
    rows = []
    for grp in groups:
        sub = adata_all[adata_all.obs["histology"] == grp]
        Xs = sub.X
        for g in top20_atc:
            j = adata_all.var_names.get_loc(g)
            col = Xs[:, j]
            mean_expr = float(np.asarray(col).mean())
            pct_expr = float(np.asarray(col > 0).mean())
            rows.append({"gene": g, "histology": grp, "mean_expr": mean_expr, "pct_expr": pct_expr * 100})
    dot_df = pd.DataFrame(rows)
    fig2 = px.scatter(
        dot_df, x="histology", y="gene", size="pct_expr", color="mean_expr",
        color_continuous_scale="plasma", template="plotly_dark",
        title="Top-20 ATC-up genes: mean expression and pct expressing",
        labels={"mean_expr": "mean (log-norm)", "pct_expr": "% expressing"},
        size_max=22, category_orders={"histology": groups, "gene": top20_atc[::-1]},
    )
    fig2.update_layout(paper_bgcolor=BG, plot_bgcolor=BG, font_color="#e2e8f0",
                       coloraxis_colorbar=dict(bgcolor=BG))
    fig2.write_html(OUT / "v17_q19_top20_genes_dotplot.html", include_plotlyjs="cdn")

    # --- Per-cell heatmap: 8-panel + top 10 novel ---
    novel_top10 = novel.head(10)["gene"].tolist()
    heat_genes = [p for p in panel_set] + novel_top10
    # sample 1500 cells across NORM/PTC/ATC, sorted by histology
    rng = np.random.default_rng(42)
    sampled_idx = []
    for grp in groups:
        idx = np.where(adata_all.obs["histology"].values == grp)[0]
        n = min(500, len(idx))
        sampled_idx.append(rng.choice(idx, size=n, replace=False))
    sampled_idx = np.concatenate(sampled_idx)
    sub = adata_all[sampled_idx, :]
    gene_idx = [adata_all.var_names.get_loc(g) for g in heat_genes]
    M = np.asarray(sub.X[:, gene_idx])
    # z-score per gene for visualization
    M_z = (M - M.mean(axis=0)) / (M.std(axis=0) + 1e-6)
    M_z = np.clip(M_z.T, -3, 3)  # genes x cells
    histo_order = sub.obs["histology"].values
    # color bar above
    fig3 = go.Figure(data=go.Heatmap(
        z=M_z, x=np.arange(M_z.shape[1]), y=heat_genes,
        colorscale="RdBu_r", zmin=-3, zmax=3, colorbar=dict(title="z-score"),
    ))
    # histology shading via shapes
    boundaries = []
    prev = None
    for i, h in enumerate(histo_order):
        if h != prev:
            boundaries.append((i, h))
            prev = h
    boundaries.append((len(histo_order), None))
    band_color = {"NORM": "#22d3ee", "PTC": "#fde047", "ATC": "#f87171"}
    shapes = []
    annots = []
    for k in range(len(boundaries) - 1):
        x0, h = boundaries[k]
        x1, _ = boundaries[k + 1]
        shapes.append(dict(type="rect", xref="x", yref="paper",
                           x0=x0 - 0.5, x1=x1 - 0.5, y0=1.01, y1=1.04,
                           fillcolor=band_color.get(h, "#888"), line=dict(width=0)))
        annots.append(dict(x=(x0 + x1) / 2, y=1.07, xref="x", yref="paper",
                           text=h, showarrow=False, font=dict(color=band_color.get(h, "#fff"), size=11)))
    fig3.update_layout(
        title="Per-cell expression: 8-panel + top-10 novel ATC-up genes (sampled, sorted by histology)",
        template="plotly_dark", paper_bgcolor=BG, plot_bgcolor=BG, font_color="#e2e8f0",
        shapes=shapes, annotations=annots, height=620,
        xaxis=dict(title="cells (n=1500)", showticklabels=False),
        yaxis=dict(title="gene"),
    )
    fig3.write_html(OUT / "v17_q19_panel_vs_novel_heatmap.html", include_plotlyjs="cdn")

    # --- summary JSON ---
    summary = {
        "dataset": "GSE193581 Lu 2023 (malignant cells only)",
        "n_cells": {"ATC": int((adata_pa.obs['histology']=='ATC').sum()),
                    "PTC": int((adata_pa.obs['histology']=='PTC').sum()),
                    "NORM": int((adata_all.obs['histology']=='NORM').sum())},
        "method": "scanpy.rank_genes_groups Wilcoxon, use_raw=False, HVG=2000",
        "effect_size_note": "adata.X is scaled (z-units); 'log2FC' is undefined. We report wilcoxon_z and delta_mean_scaled (ATC mean − PTC mean on scaled X) instead.",
        "top20_ATC_up": atc_up_full.head(20)[["gene","wilcoxon_z","delta_mean_scaled","padj","pct_in_ATC","pct_in_PTC"]].to_dict(orient="records"),
        "top20_PTC_up": ptc_up_full.head(20)[["gene","wilcoxon_z","delta_mean_scaled","padj","pct_in_ATC","pct_in_PTC"]].to_dict(orient="records"),
        "panel8_rank_in_PTC_up_list": panel_rank_df.to_dict(orient="records"),
        "novel_marker_candidates_top": novel.head(20)[["gene","wilcoxon_z","delta_mean_scaled","padj","pct_in_ATC","pct_in_PTC","category"]].to_dict(orient="records"),
        "novel_marker_watchlist_hits_top": novel_high[["gene","wilcoxon_z","delta_mean_scaled","padj","pct_in_ATC","pct_in_PTC","category"]].to_dict(orient="records"),
        "figures": [
            "v17_q19_volcano_ATC_vs_PTC.html",
            "v17_q19_top20_genes_dotplot.html",
            "v17_q19_panel_vs_novel_heatmap.html",
        ],
        "files": ["atc_vs_ptc_top200.tsv", "ptc_vs_atc_top200.tsv"],
        "author": "Seungho Cook",
    }
    with open(OUT / "q19_summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=float)
    print("[done] outputs in", OUT)


if __name__ == "__main__":
    main()
