#!/usr/bin/env python3
from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from pilot_utils import clean_gene_symbol, compute_module_scores, config, ensure_standard_dirs, fdr_bh, gene_sets, pilot_root, savefig, setup_logging, warn, write_not_run_table


MECH = {
    "RAI redifferentiation / MAPK": ["BRAF", "RAF", "MEK", "MAPK", "ERK", "RET", "NTRK", "TRAMETINIB", "SELUMETINIB", "DABRAFENIB", "VEMURAFENIB"],
    "HLA/APM restoration / IFN-epigenetic": ["JAK", "STAT", "IFN", "HDAC", "DNMT", "AZA", "DECITABINE", "EPIGEN"],
    "Proliferation stress": ["CDK", "AURK", "MYC", "ATR", "CHEK", "BCL", "MCL"],
    "Myeloid/CAF barrier validation only": ["CSF1", "TGFB", "CXCR4", "FAP"],
}


def needed_genes() -> list[str]:
    genes = []
    for payload in gene_sets().values():
        genes.extend(payload.get("core", []) or [])
        genes.extend(payload.get("barrier_negative", []) or [])
    return list(dict.fromkeys(clean_gene_symbol(g) for g in genes if g))


def load_depmap_expression(path: Path, logger) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    header = pd.read_csv(path, nrows=0).columns.tolist()
    meta = [c for c in ["ModelID", "CellLineName", "OncotreeLineage", "OncotreePrimaryDisease", "OncotreeSubtype"] if c in header]
    genes = set(needed_genes())
    gene_cols = []
    rename = {}
    for c in header:
        symbol = clean_gene_symbol(re.sub(r"\s+\(\d+\)$", "", c))
        if symbol in genes:
            gene_cols.append(c)
            rename[c] = symbol
    usecols = list(dict.fromkeys(meta + gene_cols))
    logger.info("Reading DepMap expression subset: %d columns from %s", len(usecols), path)
    df = pd.read_csv(path, usecols=usecols)
    df = df.rename(columns=rename)
    gene_df = df[[c for c in df.columns if c in genes]].copy()
    meta_df = df[[c for c in meta if c in df.columns]].copy()
    scores, cov = compute_module_scores(gene_df, gene_sets())
    out = pd.concat([meta_df, scores], axis=1)
    cov.to_csv(pilot_root() / "results" / "tables" / "drug_pilot_depmap_expression_signature_coverage.tsv", sep="\t", index=False)
    return out


def mechanism(row: pd.Series) -> str:
    text = f"{row.get('target','')} {row.get('drug_label','')}".upper()
    for label, pats in MECH.items():
        if any(p in text for p in pats):
            return label
    return "Other exploratory"


def evidence_category(row: pd.Series) -> str:
    n_t = row.get("n_thyroid_models", 0)
    p_t = row.get("p_lineage_sensitivity_thyroid", np.nan)
    p_pan = row.get("p_lineage_sensitivity", np.nan)
    mech = row.get("mechanism_axis", "")
    if n_t >= 10 and pd.notna(p_t) and p_t < 0.10 and mech != "Other exploratory":
        return "A: thyroid line sensitivity + mechanistic relevance"
    if pd.notna(p_pan) and p_pan < 0.05 and mech != "Other exploratory":
        return "B: pan-cancer association + thyroid expression support"
    if mech != "Other exploratory":
        return "C: literature/plausibility only"
    return "D: insufficient evidence"


def main() -> None:
    parser = argparse.ArgumentParser(description="Public DepMap/PRISM drug vulnerability pilot.")
    parser.parse_args()
    ensure_standard_dirs()
    logger = setup_logging("07_depm_prism_drug_pilot")
    cfg = config().get("local_sources", {})
    tables = pilot_root() / "results" / "tables"
    fig_dir = pilot_root() / "results" / "figures"

    expr_scores = load_depmap_expression(Path(cfg.get("depmap_expression", "")), logger)
    lineage_path = Path(cfg.get("depmap_lineage_scores", ""))
    if lineage_path.exists():
        meta = pd.read_csv(lineage_path, sep="\t")
        thyroid = meta[meta["OncotreeLineage"].astype(str).str.contains("Thyroid", case=False, na=False)].copy()
        if not expr_scores.empty and "ModelID" in expr_scores.columns:
            thyroid = thyroid.merge(expr_scores.drop(columns=[c for c in expr_scores.columns if c in thyroid.columns and c != "ModelID"], errors="ignore"), on="ModelID", how="left")
        thyroid.to_csv(tables / "drug_pilot_thyroid_cell_lines.tsv", sep="\t", index=False)
    else:
        thyroid = pd.DataFrame()
        write_not_run_table(tables / "drug_pilot_thyroid_cell_lines.tsv", "DepMap model metadata unavailable")

    assoc_path = Path(cfg.get("depmap_drug_associations", ""))
    if assoc_path.exists():
        assoc = pd.read_csv(assoc_path, sep="\t")
        assoc["mechanism_axis"] = assoc.apply(mechanism, axis=1)
        assoc["evidence_category"] = assoc.apply(evidence_category, axis=1)
        assoc["fdr_pan_cancer"] = fdr_bh(assoc["p_lineage_sensitivity"])
        assoc["fdr_thyroid"] = fdr_bh(assoc["p_lineage_sensitivity_thyroid"])
        assoc["rank_score"] = (
            assoc["spearman_lineage_sensitivity_thyroid"].abs().fillna(0) * np.log1p(assoc["n_thyroid_models"].fillna(0))
            + assoc["spearman_lineage_sensitivity_metric"].abs().fillna(0)
        )
        cand = assoc.sort_values(["evidence_category", "rank_score"], ascending=[True, False]).copy()
        cand.to_csv(tables / "drug_pilot_candidate_rankings.tsv", sep="\t", index=False)
        assoc.to_csv(tables / "drug_pilot_signature_drug_correlations.tsv", sep="\t", index=False)
    else:
        cand = pd.DataFrame()
        write_not_run_table(tables / "drug_pilot_candidate_rankings.tsv", "Drug association table unavailable")
        write_not_run_table(tables / "drug_pilot_signature_drug_correlations.tsv", "Drug association table unavailable")

    dep_path = Path(cfg.get("depmap_dependency_associations", ""))
    if dep_path.exists():
        dep = pd.read_csv(dep_path, sep="\t")
        dep["fdr_lineage_dependency"] = fdr_bh(dep["p_lineage_dependency"])
        dep.to_csv(tables / "drug_pilot_dependency_associations.tsv", sep="\t", index=False)

    if not cand.empty:
        top = cand[cand["evidence_category"].ne("D: insufficient evidence")].head(25)
        if top.empty:
            top = cand.head(25)
        fig, ax = plt.subplots(figsize=(11, max(4, top.shape[0] * 0.25)))
        heat = top.set_index("drug_label")[["spearman_lineage_sensitivity_metric", "spearman_lineage_sensitivity_thyroid"]].fillna(0)
        sns.heatmap(heat, cmap="vlag", center=0, ax=ax, cbar_kws={"label": "Spearman"})
        ax.set_title("Drug/perturbation candidate pilot")
        savefig(fig, fig_dir / "drug_pilot_candidate_heatmap.png")
        plt.close(fig)

        ex = top.dropna(subset=["spearman_lineage_sensitivity_metric", "spearman_lineage_sensitivity_thyroid"]).head(30)
        fig, ax = plt.subplots(figsize=(7, 6))
        sns.scatterplot(data=ex, x="spearman_lineage_sensitivity_metric", y="spearman_lineage_sensitivity_thyroid", hue="mechanism_axis", size="n_thyroid_models", ax=ax)
        ax.axhline(0, color="0.5", lw=1)
        ax.axvline(0, color="0.5", lw=1)
        ax.set_title("Pan-cancer vs thyroid-line drug association")
        ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8, frameon=False)
        savefig(fig, fig_dir / "drug_pilot_signature_vs_drug_scatter_examples.png")
        plt.close(fig)
    logger.info("Drug pilot complete: %d thyroid cell lines, %d candidate rows.", len(thyroid), len(cand))


if __name__ == "__main__":
    main()
