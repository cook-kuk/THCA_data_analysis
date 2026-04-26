#!/usr/bin/env python3
"""Build gene explorer data cache.

Produces two JSON files consumed by assets/js/gene-explorer.js:
  - gene_universe.json     : flat list of HGNC symbols available for search.
  - gene_expression.json   : {gene: {dataset: [values, ...]}} for a curated
                             subset (panels + top-ranked novel biomarkers),
                             so violins render client-side without huge
                             payload.

Only processed log2 matrices already on disk are read; no heavy recompute.
"""
from __future__ import annotations
import json
from pathlib import Path

import pandas as pd

PROJECT = Path("/opt/thyroid-dash/project")
OUT = PROJECT / "reports" / "html" / "assets" / "data"
OUT.mkdir(parents=True, exist_ok=True)

EXPR_SOURCES = {
    "TCGA-THCA":  PROJECT / "data_processed/bulk_rnaseq/TCGA-THCA_rnaseq_expression_log2.tsv",
    "GSE126698":  PROJECT / "data_processed/bulk_rnaseq/GSE126698_rnaseq_expression_log2.tsv",
    "GSE213647":  PROJECT / "data_processed/bulk_rnaseq/GSE213647_rnaseq_expression_log2.tsv",
    "GSE27155":   PROJECT / "data_processed/microarray/GSE27155_microarray_expression_log2.tsv",
    "GSE76039":   PROJECT / "data_processed/microarray/GSE76039_microarray_expression_log2.tsv",
}

PANELS_JSON = PROJECT / "reports/html/assets/data/gene_panels.json"
BIOMARKER_VALIDATED = PROJECT / "results/tables/biomarker_validated.tsv"


def normalize(sym: str) -> str:
    return str(sym).strip().upper().replace(" ", "")


def read_expr(path: Path) -> pd.DataFrame | None:
    if not path.exists():
        return None
    df = pd.read_csv(path, sep="\t")
    gene_col = df.columns[0]
    df = df.rename(columns={gene_col: "gene_symbol"}).set_index("gene_symbol")
    df.index = df.index.map(normalize)
    df = df[~df.index.duplicated(keep="first")]
    return df


def main() -> None:
    # --- load all matrices, build universe from intersection of ≥1 matrix
    matrices = {name: read_expr(p) for name, p in EXPR_SOURCES.items()}
    matrices = {k: v for k, v in matrices.items() if v is not None and not v.empty}
    print(f"loaded {len(matrices)} expression matrices")

    all_genes: set[str] = set()
    for name, df in matrices.items():
        all_genes |= set(df.index)
        print(f"  {name:<12s} {df.shape[0]:>6d} genes x {df.shape[1]:>4d} samples")

    # --- curated subset = panels ∪ top novel biomarkers ∪ key driver genes
    panels = json.loads(PANELS_JSON.read_text(encoding="utf-8"))
    panel_genes: set[str] = set()
    for v in panels.values():
        if isinstance(v, list):
            panel_genes.update(normalize(g) for g in v)
        elif isinstance(v, dict):
            for sub in v.values():
                if isinstance(sub, list):
                    panel_genes.update(normalize(g) for g in sub)

    novel_top: list[str] = []
    if BIOMARKER_VALIDATED.exists():
        bio = pd.read_csv(BIOMARKER_VALIDATED, sep="\t")
        gene_col = "gene" if "gene" in bio.columns else bio.columns[0]
        score_col = next((c for c in ["novelty_score", "cohens_d", "abs_cohens_d"] if c in bio.columns), None)
        if score_col:
            bio["_abs"] = bio[score_col].abs() if score_col == "cohens_d" else bio[score_col]
            bio = bio.sort_values("_abs", ascending=False)
        novel_top = [normalize(g) for g in bio[gene_col].head(150).dropna().tolist()]

    driver_and_extras = [
        "BRAF","NRAS","KRAS","HRAS","RET","NTRK1","NTRK3","ALK","PAX8","PPARG","EIF1AX","TERT",
        "TP53","CDKN2A","CDKN2B","PIK3CA","PTEN","AKT1","MET","CTNNB1","MYC","VIM","CDH1","CDH2",
        "SNAI1","SNAI2","TWIST1","ZEB1","ZEB2","KRT19","CITED1","LGALS3","CALCA","CEACAM5","TACSTD2",
        "PLEKHA6","CYP1B1","TMPRSS4","LDLR","DUSP4","DUSP5","DUSP6","SPRY1","SPRY2","SPRY4",
        "ETV4","ETV5","FOSL1","PHLDA1","SERPINA1","TIMP1","CD274","CD8A","FOXP3","IFNG","HLA-DRA",
        "KLK10","KLK6","KLK7","MMP9","LOX","CDK6","CCND1","STAT1","IRF1","CXCL10","CXCL9","IDO1",
    ]

    curated = (panel_genes | set(novel_top) | {normalize(g) for g in driver_and_extras}) & all_genes
    print(f"curated gene universe for violins: {len(curated)} genes")

    # --- export universe (full autocomplete list = every gene present in ≥1 matrix)
    universe_sorted = sorted(all_genes)
    (OUT / "gene_universe.json").write_text(
        json.dumps({"symbols": universe_sorted, "n": len(universe_sorted),
                    "curated": sorted(curated), "n_curated": len(curated)},
                   separators=(",", ":")),
        encoding="utf-8",
    )
    print(f"wrote gene_universe.json ({len(universe_sorted)} symbols, {len(curated)} curated)")

    # --- export expression cache for curated subset (samples kept to reduce size)
    cache: dict[str, dict[str, list[float]]] = {}
    for gene in sorted(curated):
        per_dataset: dict[str, list[float]] = {}
        for name, df in matrices.items():
            if gene in df.index:
                vals = df.loc[gene].astype(float).values
                # filter NaN
                vals = [float(v) for v in vals if v == v]
                if vals:
                    # Subsample large cohorts to 400 points max (TCGA 500+), keep small cohorts whole
                    if len(vals) > 400:
                        import random
                        random.seed(42)
                        vals = random.sample(vals, 400)
                    per_dataset[name] = [round(v, 3) for v in vals]
        if per_dataset:
            cache[gene] = per_dataset

    cache_path = OUT / "gene_expression.json"
    cache_path.write_text(json.dumps(cache, separators=(",", ":")), encoding="utf-8")
    size_mb = cache_path.stat().st_size / 1024 / 1024
    print(f"wrote gene_expression.json ({len(cache)} genes, {size_mb:.2f} MB)")


if __name__ == "__main__":
    main()
