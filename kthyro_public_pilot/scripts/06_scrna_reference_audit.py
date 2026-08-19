#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from scipy import sparse

from pilot_utils import clean_gene_symbol, compute_module_scores, config, ensure_standard_dirs, gene_sets, pilot_root, savefig, setup_logging, warn, write_not_run_table


MARKERS = {
    "thyrocyte_malignant": ["EPCAM", "KRT8", "KRT18", "KRT19", "PAX8", "TG", "TPO"],
    "t_cell": ["CD3D", "CD3E", "CD8A"],
    "b_cell": ["MS4A1", "CD79A"],
    "myeloid": ["LST1", "CD68", "CD163", "AIF1"],
    "fibroblast_caf": ["COL1A1", "ACTA2", "FAP", "PDGFRB"],
    "endothelial": ["PECAM1", "VWF", "KDR"],
}


def expr_subset(ad, genes: list[str]) -> pd.DataFrame:
    mapping = {}
    for i, g in enumerate(ad.var_names):
        cg = clean_gene_symbol(g)
        if cg not in mapping:
            mapping[cg] = i
    present = [g for g in genes if g in mapping]
    if not present:
        return pd.DataFrame(index=ad.obs_names)
    X = ad.X[:, [mapping[g] for g in present]]
    if sparse.issparse(X):
        X = X.toarray()
    return pd.DataFrame(X, index=ad.obs_names, columns=present)


def rough_label(expr: pd.DataFrame) -> pd.Series:
    module_means = {}
    for name, genes in MARKERS.items():
        present = [g for g in genes if g in expr]
        module_means[name] = expr[present].mean(axis=1) if present else pd.Series(-np.inf, index=expr.index)
    mat = pd.DataFrame(module_means)
    return mat.idxmax(axis=1).replace(
        {
            "thyrocyte_malignant": "malignant/thyrocyte-like",
            "t_cell": "T cell",
            "b_cell": "B cell",
            "myeloid": "myeloid",
            "fibroblast_caf": "fibroblast/CAF",
            "endothelial": "endothelial",
        }
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="Audit scRNA reference marker/module expression.")
    parser.parse_args()
    ensure_standard_dirs()
    logger = setup_logging("06_scrna_reference_audit")
    path = Path(config().get("local_sources", {}).get("scrna_lu2023_h5ad", ""))
    tables = pilot_root() / "results" / "tables"
    fig_dir = pilot_root() / "results" / "figures"
    if not path.exists():
        reason = "No local scRNA h5ad found. Proposed fallback: Cell2location/RCTD/SPOTlight using thyroid scRNA reference."
        write_not_run_table(tables / "scrna_celltype_module_scores.tsv", reason)
        (pilot_root() / "results" / "reports" / "scrna_deconvolution_plan.md").write_text(
            "# scRNA Reference Not Run\n\n" + reason + "\n\nSpatial HLA-II/CD74 may originate from APCs; CAF/myeloid modules need mIHC/GeoMx validation.\n"
        )
        logger.warning(reason)
        return
    ad = sc.read_h5ad(path)
    logger.info("Loaded scRNA reference %s: %d cells x %d genes", path, ad.n_obs, ad.n_vars)
    all_genes = []
    for payload in gene_sets().values():
        all_genes.extend(payload.get("core", []) or [])
        all_genes.extend(payload.get("barrier_negative", []) or [])
    for genes in MARKERS.values():
        all_genes.extend(genes)
    all_genes = list(dict.fromkeys(clean_gene_symbol(g) for g in all_genes))
    expr = expr_subset(ad, all_genes)
    if "author_celltype" in ad.obs:
        celltype = ad.obs["author_celltype"].astype(str)
    else:
        celltype = rough_label(expr)
    scores, coverage = compute_module_scores(expr, gene_sets())
    scores["cell_type"] = celltype.values
    scores["sample"] = ad.obs["sample"].astype(str).values if "sample" in ad.obs else "unknown"
    score_cols = [c for c in scores.columns if c.endswith("_score")]
    summary = scores.groupby("cell_type")[score_cols].agg(["mean", "median", "count"])
    summary.columns = ["_".join(c) for c in summary.columns]
    summary = summary.reset_index()
    summary.to_csv(tables / "scrna_celltype_module_scores.tsv", sep="\t", index=False)
    coverage.to_csv(tables / "scrna_gene_set_coverage.tsv", sep="\t", index=False)

    marker_genes = [g for genes in MARKERS.values() for g in genes if g in expr.columns]
    if marker_genes:
        means = expr[marker_genes].assign(cell_type=celltype.values).groupby("cell_type").mean()
        means = means.loc[means.index[:30]]
        fig, ax = plt.subplots(figsize=(max(8, len(marker_genes) * 0.35), max(4, means.shape[0] * 0.25)))
        sns.heatmap(means, cmap="viridis", ax=ax)
        ax.set_title("scRNA marker expression by annotated cell type")
        savefig(fig, fig_dir / "scrna_marker_dotplot.png")
        plt.close(fig)

    plot_scores = [c for c in ["hla_i_apm_score", "hla_ii_apc_score", "cytotoxic_t_score", "myeloid_tam_score", "caf_ecm_tgfb_score", "tumor_epithelial_score"] if c in scores]
    if plot_scores:
        long = scores.sample(min(20000, len(scores)), random_state=1).melt(id_vars="cell_type", value_vars=plot_scores, var_name="module", value_name="score")
        fig, ax = plt.subplots(figsize=(12, 6))
        sns.boxplot(data=long, x="module", y="score", hue="cell_type", fliersize=0, ax=ax)
        ax.tick_params(axis="x", rotation=35)
        ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False, fontsize=7)
        ax.set_title("scRNA module score sanity check by cell type")
        savefig(fig, fig_dir / "scrna_module_score_by_celltype.png")
        plt.close(fig)
    logger.info("Wrote scRNA cell-type module summary for %d cell types.", summary.shape[0])


if __name__ == "__main__":
    main()
