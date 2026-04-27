#!/usr/bin/env python3
from __future__ import annotations

import json

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import scanpy as sc
from inmoose.pycombat import pycombat_norm
from scipy import stats

from v17_common import FIG, LOG, TAB, TDS16, TIERA67, align_dataset_expr_meta, bh_fdr, log_line, signature_score


def main() -> None:
    log_line(LOG, "task5 start")
    sm = pd.read_csv(TAB / "sample_master_v17_full.tsv", sep="\t")
    dark_path = TAB / "dark_matter_cohort.tsv"
    if dark_path.exists() and "v17_dark_cluster" not in sm.columns:
        sm = sm.merge(pd.read_csv(dark_path, sep="\t")[["sample_id", "v17_dark_cluster"]], on="sample_id", how="left")
    tumor = sm[sm["normal_vs_tumor"] != "normal"].copy()
    tumor = tumor.sort_values(["sample_id", "dataset"]).drop_duplicates("sample_id", keep="first").reset_index(drop=True)
    datasets = ["TCGA-THCA", "GSE27155", "GSE126698", "GSE76039", "GSE213647"]
    expr_parts = []
    meta_parts = []
    common_genes = None
    for ds in datasets:
        expr, meta = align_dataset_expr_meta(ds, genes=TIERA67, tumor_only=True)
        if expr.empty or meta.empty:
            log_line(LOG, f"task5 skip dataset={ds} empty alignment")
            continue
        meta = meta.merge(tumor[["sample_id", "driver_anchor_v17"] + (["v17_dark_cluster"] if "v17_dark_cluster" in tumor.columns else [])], on="sample_id", how="left")
        common = expr.index.intersection(meta["sample_id"])
        if len(common) == 0:
            log_line(LOG, f"task5 skip dataset={ds} no common samples")
            continue
        meta_use = meta.drop_duplicates("sample_id").set_index("sample_id").loc[common].reset_index().rename(columns={"index": "sample_id"})
        expr_use = expr.loc[common]
        common_genes = set(expr_use.columns) if common_genes is None else (common_genes & set(expr_use.columns))
        expr_parts.append(expr_use)
        meta_parts.append(meta_use)
    common_genes = sorted(common_genes or [])
    if not expr_parts or not common_genes:
        raise RuntimeError("task5 could not build a non-empty multi-cohort expression matrix")
    expr_parts = [x.loc[:, common_genes].copy() for x in expr_parts]
    expr_all = pd.concat(expr_parts, axis=0)
    meta_all = pd.concat(meta_parts, axis=0).reset_index(drop=True)
    try:
        corrected = pycombat_norm(expr_all.T, batch=meta_all["dataset"].tolist(), covar_mod=None)
        corrected = corrected.T if isinstance(corrected, pd.DataFrame) else pd.DataFrame(np.asarray(corrected).T, index=expr_all.index, columns=expr_all.columns)
        log_line(LOG, f"task5 combat success genes={len(common_genes)} samples={len(expr_all)}")
    except Exception as e:
        log_line(LOG, f"task5 combat failed, fallback to uncorrected matrix: {e}")
        corrected = expr_all.copy()

    adata = sc.AnnData(corrected.values)
    adata.obs_names = corrected.index
    adata.var_names = corrected.columns
    adata.obs = meta_all.set_index("sample_id").loc[corrected.index]
    adata.obs["tds_score_v17"] = signature_score(corrected, TDS16).values
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=min(20, adata.n_vars - 1, adata.n_obs - 1), svd_solver="arpack")
    sc.pp.neighbors(adata, n_neighbors=min(15, adata.n_obs - 1), random_state=42)
    sc.tl.umap(adata, random_state=42, n_components=3)
    sc.tl.diffmap(adata)
    root_idx = int(np.nanargmax(adata.obs["tds_score_v17"].to_numpy(float)))
    adata.uns["iroot"] = root_idx
    sc.tl.dpt(adata)

    out = adata.obs.reset_index().rename(columns={"index": "sample_id", "dpt_pseudotime": "pseudotime"})
    out.to_csv(TAB / "trajectory_pseudotime.tsv", sep="\t", index=False)

    cors = []
    pt = out.set_index("sample_id").loc[corrected.index, "pseudotime"].to_numpy(float)
    for g in corrected.columns:
        rho, p = stats.spearmanr(corrected[g].to_numpy(float), pt, nan_policy="omit")
        cors.append({"gene": g, "spearman_rho": rho, "pvalue": p})
    dyn = pd.DataFrame(cors)
    dyn["fdr"] = bh_fdr(dyn["pvalue"].to_numpy(float))
    dyn = dyn.sort_values("spearman_rho", key=lambda s: s.abs(), ascending=False)
    dyn.to_csv(TAB / "trajectory_dynamic_genes.tsv", sep="\t", index=False)

    umap3 = pd.DataFrame(adata.obsm["X_umap"], columns=["UMAP1", "UMAP2", "UMAP3"], index=corrected.index)
    plot_df = out.set_index("sample_id").join(umap3).reset_index()
    fig = px.scatter_3d(plot_df, x="UMAP1", y="UMAP2", z="UMAP3", color="pseudotime", hover_data=["driver_anchor_v17", "histology_subtype", "dataset"])
    fig.update_layout(title="v17 trajectory UMAP 3D")
    fig.write_html(FIG / "trajectory_umap_3d.html", include_plotlyjs="cdn")

    grad = px.scatter(plot_df, x="pseudotime", y="tds_score_v17", color="driver_anchor_v17", title="TDS gradient by pseudotime")
    grad.write_html(FIG / "trajectory_tds_gradient.html", include_plotlyjs="cdn")

    high = plot_df.copy()
    cutoff = high["pseudotime"].quantile(0.9)
    high["high_risk_like"] = np.where(high["pseudotime"] >= cutoff, "high", "other")
    fig2 = px.scatter(high, x="UMAP1", y="UMAP2", color="high_risk_like", symbol="driver_anchor_v17", title="High-risk-like trajectory tail")
    fig2.write_html(FIG / "trajectory_high_risk_cluster.html", include_plotlyjs="cdn")
    cluster_fig = px.scatter(
        plot_df[pd.notna(plot_df.get("v17_dark_cluster"))],
        x="UMAP1",
        y="UMAP2",
        color="v17_dark_cluster",
        hover_data=["driver_anchor_v17", "histology_subtype", "dataset"],
        title="Dark matter clusters projected onto global trajectory",
    )
    cluster_fig.write_html(FIG / "trajectory_dark_cluster_projection.html", include_plotlyjs="cdn")

    log_line(LOG, f"task5 done samples={len(out)}")


if __name__ == "__main__":
    main()
