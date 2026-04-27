#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import scanpy as sc
from scipy import stats
from sklearn.cluster import AgglomerativeClustering
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

from v17_common import FIG, LOG, TAB, BRS71, RAI_GENES, TDS16, TIERA67, align_dataset_expr_meta, bh_fdr, load_sample_master, log_line, signature_score


OUT_COHORT = TAB / "dark_matter_cohort.tsv"
OUT_MARKERS = TAB / "dark_matter_cluster_markers.tsv"
OUT_CLIN = TAB / "dark_matter_cluster_clinical.tsv"


def leiden_target_k(X: np.ndarray, n_clusters: int, seed: int) -> np.ndarray:
    adata = sc.AnnData(X)
    sc.pp.neighbors(adata, n_neighbors=min(15, max(5, X.shape[0] - 1)), use_rep="X", random_state=seed)
    lo, hi = 0.05, 3.0
    best = None
    for _ in range(12):
        mid = (lo + hi) / 2
        sc.tl.leiden(adata, resolution=mid, key_added="leiden_tmp", random_state=seed)
        labels = adata.obs["leiden_tmp"].astype(str).to_numpy()
        k = len(np.unique(labels))
        best = labels
        if k == n_clusters:
            break
        if k < n_clusters:
            lo = mid
        else:
            hi = mid
    return pd.factorize(best)[0]


def cluster_grid(X: np.ndarray) -> tuple[np.ndarray, pd.DataFrame, np.ndarray]:
    seeds = list(range(42, 72))
    best = None
    best_score = -np.inf
    best_cons = None
    rows = []
    for k in range(2, 11):
        mats = []
        label_runs = []
        for seed in seeds:
            try:
                labels = leiden_target_k(X, k, seed)
            except Exception:
                labels = AgglomerativeClustering(n_clusters=k, linkage="ward").fit_predict(X)
            label_runs.append(labels)
            mats.append((labels[:, None] == labels[None, :]).astype(float))
        consensus = np.mean(mats, axis=0)
        scores = []
        for labels in label_runs:
            for c in np.unique(labels):
                idx = np.where(labels == c)[0]
                if len(idx) <= 1:
                    continue
                scores.append(consensus[np.ix_(idx, idx)].mean())
        stability = float(np.mean(scores)) if scores else np.nan
        rows.append({"k": k, "stability": stability})
        if stability > best_score:
            best_score = stability
            best = label_runs[0]
            best_cons = consensus
    return best, pd.DataFrame(rows), best_cons


def marker_table(expr: pd.DataFrame, labels: np.ndarray) -> pd.DataFrame:
    rows = []
    for c in sorted(np.unique(labels)):
        in_mask = labels == c
        out_mask = ~in_mask
        pvals = []
        lfc = []
        genes = expr.columns.tolist()
        for g in genes:
            a = expr.loc[in_mask, g].astype(float).values
            b = expr.loc[out_mask, g].astype(float).values
            lfc.append(float(np.nanmean(a) - np.nanmean(b)))
            try:
                pvals.append(stats.ttest_ind(a, b, equal_var=False, nan_policy="omit").pvalue)
            except Exception:
                pvals.append(np.nan)
        fdr = bh_fdr(np.asarray(pvals))
        sub = pd.DataFrame({"cluster": c, "gene": genes, "log2fc": lfc, "pvalue": pvals, "fdr": fdr})
        sub = sub.sort_values(["fdr", "log2fc"], ascending=[True, False]).head(20)
        rows.append(sub)
    return pd.concat(rows, ignore_index=True)


def main() -> None:
    log_line(LOG, "task1 start")
    expr, meta = align_dataset_expr_meta("TCGA-THCA", genes=TIERA67, tumor_only=True)
    mask = ~meta["driver_anchor"].isin(["BRAF", "RAS"])
    expr = expr.loc[mask.values].copy()
    meta = meta.loc[mask.values].copy()
    if expr.empty:
        raise RuntimeError("no TCGA dark matter cohort found")

    X = StandardScaler().fit_transform(expr.fillna(expr.mean()).values)
    pca = PCA(n_components=min(10, X.shape[1], X.shape[0] - 1), random_state=42)
    pcs = pca.fit_transform(X)
    labels, stability_df, consensus = cluster_grid(pcs)
    meta["v17_dark_cluster"] = [f"DM{int(x)+1}" for x in labels]
    meta["tds16_score_v17"] = signature_score(expr, TDS16).values
    meta["brs71_score_v17"] = signature_score(expr, BRS71).values
    meta["rai_score_v17"] = signature_score(expr, RAI_GENES).values

    markers = marker_table(expr, labels)

    clin = meta.groupby("v17_dark_cluster").agg(
        n_samples=("sample_id", "size"),
        tds16_mean=("tds16_score_v17", "mean"),
        brs71_mean=("brs71_score_v17", "mean"),
        rai_mean=("rai_score_v17", "mean"),
        age_mean=("age", "mean"),
    ).reset_index()

    contingency = pd.crosstab(meta["v17_dark_cluster"], meta["molecular_subtype"]).reset_index()
    clin = clin.merge(contingency, on="v17_dark_cluster", how="left")

    meta.to_csv(OUT_COHORT, sep="\t", index=False)
    markers.to_csv(OUT_MARKERS, sep="\t", index=False)
    clin.to_csv(OUT_CLIN, sep="\t", index=False)
    stability_df.to_csv(TAB / "dark_matter_cluster_stability.tsv", sep="\t", index=False)

    adata = sc.AnnData(X)
    sc.pp.neighbors(adata, n_neighbors=min(15, len(meta) - 1), random_state=42)
    sc.tl.umap(adata, random_state=42)
    umap = pd.DataFrame(adata.obsm["X_umap"], columns=["UMAP1", "UMAP2"], index=meta.index)
    fig = px.scatter(
        pd.concat([meta.reset_index(drop=True), umap.reset_index(drop=True)], axis=1),
        x="UMAP1", y="UMAP2", color="v17_dark_cluster",
        hover_data=["sample_id", "histology_subtype", "driver_anchor", "rai_score_v17"],
        title="v17 Dark Matter UMAP"
    )
    fig.write_html(FIG / "dark_matter_umap.html", include_plotlyjs="cdn")

    heat = go.Figure(data=go.Heatmap(z=consensus, colorscale="Viridis"))
    heat.update_layout(title="Dark Matter Consensus Heatmap")
    heat.write_html(FIG / "dark_matter_consensus_heatmap.html", include_plotlyjs="cdn")

    top_genes = markers.groupby("cluster")["gene"].apply(list).explode().drop_duplicates().tolist()[:40]
    hm_mat = expr[top_genes].copy()
    hm_mat["cluster"] = meta["v17_dark_cluster"].values
    hm_mat = hm_mat.sort_values("cluster")
    heat2 = go.Figure(data=go.Heatmap(z=hm_mat[top_genes].T.values, x=hm_mat.index.astype(str), y=top_genes, colorscale="RdBu"))
    heat2.update_layout(title="Dark Matter Marker Heatmap")
    heat2.write_html(FIG / "dark_matter_marker_heatmap.html", include_plotlyjs="cdn")

    box = px.box(meta, x="v17_dark_cluster", y="rai_score_v17", points="all", color="v17_dark_cluster", title="RAI uptake score by dark cluster")
    box.write_html(FIG / "dark_matter_rai_boxplot.html", include_plotlyjs="cdn")

    summary = {
        "n_samples": int(len(meta)),
        "cluster_counts": meta["v17_dark_cluster"].value_counts().sort_index().to_dict(),
        "best_k": int(len(np.unique(labels))),
        "stability": float(stability_df.sort_values("stability", ascending=False).iloc[0]["stability"]),
    }
    (TAB / "dark_matter_summary.json").write_text(json.dumps(summary, indent=2))
    log_line(LOG, f"task1 done {json.dumps(summary)}")


if __name__ == "__main__":
    main()
