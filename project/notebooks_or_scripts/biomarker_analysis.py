#!/usr/bin/env python3
"""
THCA Biomarker Analysis: Known vs Novel BRAF_like vs RAS_like candidates.

Training cohort: TCGA-THCA RNA-seq (log2).
External cohorts: GSE27155 microarray (log2), GSE126698 RNA-seq (log2).

Outputs (all NEW files; does not touch any existing script/template):
  - results/tables/biomarker_known_vs_novel.tsv
  - results/tables/biomarker_de_full.tsv
  - results/tables/biomarker_validated.tsv
  - reports/html/assets/data/biomarker_payload.json
  - reports/html/figs_interactive/volcano_braf_vs_ras.html
  - reports/html/figs_interactive/biomarker_heatmap_known_vs_novel.html
  - reports/html/figs_interactive/biomarker_effect_vs_coverage.html
  - reports/html/figs_interactive/biomarker_replication_scatter.html
  - reports/biomarker_analysis.md
"""

from __future__ import annotations

import json
import math
import os
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

import plotly.graph_objects as go
from plotly.subplots import make_subplots


# ----------------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------------
PROJECT_ROOT = Path("/home/seungho/personal/THCA_data_analysis/project")
DATA_TCGA = PROJECT_ROOT / "data_processed/bulk_rnaseq/TCGA-THCA_rnaseq_expression_log2.tsv"
DATA_GSE27155 = PROJECT_ROOT / "data_processed/microarray/GSE27155_microarray_expression_log2.tsv"
DATA_GSE126698 = PROJECT_ROOT / "data_processed/bulk_rnaseq/GSE126698_rnaseq_expression_log2.tsv"
SAMPLE_MASTER = PROJECT_ROOT / "metadata/sample_master.tsv"
PANEL_TDS16 = PROJECT_ROOT / "metadata/tds16_genes.txt"
PANEL_TIERA67 = PROJECT_ROOT / "metadata/tierA67_genes.txt"
PANEL_BRS71 = PROJECT_ROOT / "metadata/brs71_genes.txt"

OUT_TABLES = PROJECT_ROOT / "results/tables"
OUT_FIGS_INTERACTIVE = PROJECT_ROOT / "reports/html/figs_interactive"
OUT_DATA_ASSETS = PROJECT_ROOT / "reports/html/assets/data"
OUT_REPORTS = PROJECT_ROOT / "reports"

OUT_TABLES.mkdir(parents=True, exist_ok=True)
OUT_FIGS_INTERACTIVE.mkdir(parents=True, exist_ok=True)
OUT_DATA_ASSETS.mkdir(parents=True, exist_ok=True)
OUT_REPORTS.mkdir(parents=True, exist_ok=True)


# ----------------------------------------------------------------------
# Step 1: Known biomarker reference set (hardcoded with citations)
# ----------------------------------------------------------------------
KNOWN_BIOMARKERS: Dict[str, Dict[str, str]] = {}


def _add_known(genes: List[str], category: str, citation: str) -> None:
    for g in genes:
        if g not in KNOWN_BIOMARKERS:
            KNOWN_BIOMARKERS[g] = {"category": category, "citation": citation}
        else:
            # append category if not already present
            cats = KNOWN_BIOMARKERS[g]["category"].split("; ")
            if category not in cats:
                KNOWN_BIOMARKERS[g]["category"] = "; ".join(cats + [category])


# Thyroid-lineage / TDS core (Cell 2014 PMID 25417114)
_add_known(
    ["TG", "TPO", "TSHR", "SLC5A5", "SLC26A4", "PAX8", "NKX2-1",
     "FOXE1", "DIO1", "DIO2", "DUOX1", "DUOX2", "GLIS3",
     "THRA", "THRB", "SLC5A8", "IYD"],
    category="Thyroid_lineage_TDS",
    citation="TCGA PTC Cell 2014 (PMID 25417114)",
)
# Driver mutations/fusions (expression proxies)
_add_known(
    ["BRAF", "NRAS", "HRAS", "KRAS", "RET", "NTRK1", "NTRK3", "ALK",
     "PPARG", "EIF1AX", "TERT"],
    category="Driver_mutation_fusion_proxy",
    citation="TCGA PTC Cell 2014 (PMID 25417114)",
)
# MAPK output (Pratilas 2009 PMID 19470482)
_add_known(
    ["DUSP4", "DUSP5", "DUSP6", "SPRY1", "SPRY2", "SPRY4",
     "ETV4", "ETV5", "PHLDA1", "FOSL1"],
    category="MAPK_output",
    citation="Pratilas 2009 (PMID 19470482)",
)
# Aggressiveness / dedifferentiation / EMT
_add_known(
    ["TP53", "CDKN2A", "CDKN2B", "PIK3CA", "PTEN", "AKT1", "MET",
     "VIM", "CDH1", "CDH2", "SNAI1", "SNAI2", "ZEB1", "ZEB2",
     "TWIST1", "MMP9", "LOX"],
    category="Aggressiveness_EMT",
    citation="Landa 2016 (PMID 27525637); Xing 2013 (PMID 23550066)",
)
# Clinical markers
_add_known(
    ["CALCA", "CEACAM5", "CITED1", "LGALS3", "KRT19", "MSLN"],
    category="Clinical_IHC",
    citation="de Matos 2013 (PMID 23956041); Nikiforova IHC practice",
)
# BRS-adjacent literature signature (BRAF-like/RAS-like)
_add_known(
    ["ERBB3", "LRP4", "TIMP1", "SERPINA1", "ANGPTL4", "HMGA2"],
    category="BRS_literature",
    citation="Yoo 2016 (PMID 27601545); Chakravarty 2011 (PMID 22105173)",
)


def load_panel(path: Path) -> List[str]:
    genes = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#") or line.startswith("["):
            continue
        genes.append(line)
    return genes


PANEL_SETS: Dict[str, List[str]] = {
    "TDS16": load_panel(PANEL_TDS16),
    "TierA67": load_panel(PANEL_TIERA67),
    "BRS71": load_panel(PANEL_BRS71),
}

PANEL_MEMBERSHIP: Dict[str, List[str]] = {}
for panel_name, panel_genes in PANEL_SETS.items():
    for g in panel_genes:
        PANEL_MEMBERSHIP.setdefault(g, []).append(panel_name)
for g in KNOWN_BIOMARKERS:
    PANEL_MEMBERSHIP.setdefault(g, []).append("HARDCODED_REF")


def is_known(gene: str) -> bool:
    return gene in KNOWN_BIOMARKERS or any(
        gene in panel_genes for panel_genes in PANEL_SETS.values()
    )


# ----------------------------------------------------------------------
# Step 2/3 helpers
# ----------------------------------------------------------------------
def load_expression(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", index_col=0)
    # Collapse duplicate gene symbols by max mean expression
    if df.index.has_duplicates:
        df["_mean"] = df.mean(axis=1, numeric_only=True)
        df = df.sort_values("_mean", ascending=False)
        df = df[~df.index.duplicated(keep="first")]
        df = df.drop(columns="_mean")
    return df


def load_sample_master() -> pd.DataFrame:
    df = pd.read_csv(SAMPLE_MASTER, sep="\t", low_memory=False)
    return df


def derive_126698_labels(columns: List[str]) -> Dict[str, str]:
    """GSE126698 column codes: A*=ATC, F*=FTC (RAS_like), N*=normal, P*=PTC (BRAF_like).
    Only BRAF_like / RAS_like retained for the subtype DE.
    """
    mapping = {}
    for c in columns:
        prefix = c[0].upper() if c else ""
        if prefix == "P":
            mapping[c] = "BRAF_like"
        elif prefix == "F":
            mapping[c] = "RAS_like"
        # A (ATC) and N (normal) are excluded from BRAF vs RAS test
    return mapping


def split_samples(
    sample_master: pd.DataFrame, dataset: str, expr_columns: List[str]
) -> Tuple[List[str], List[str]]:
    """Return (braf_like_samples, ras_like_samples) for a given dataset,
    restricted to tumor samples present in the expression matrix.
    """
    sub = sample_master[sample_master["dataset"] == dataset].copy()
    # restrict to tumor where available (TCGA-THCA has Primary Tumor; others vary)
    if "normal_vs_tumor" in sub.columns:
        tumor_rows = sub[sub["normal_vs_tumor"] == "tumor"]
        if len(tumor_rows) > 0:
            sub = tumor_rows
    present = set(expr_columns)
    sub = sub[sub["sample_id"].isin(present)]
    braf = sub.loc[sub["molecular_subtype"] == "BRAF_like", "sample_id"].tolist()
    ras = sub.loc[sub["molecular_subtype"] == "RAS_like", "sample_id"].tolist()
    return braf, ras


def welch_de(expr: pd.DataFrame, braf_cols: List[str], ras_cols: List[str]) -> pd.DataFrame:
    """Per-gene Welch t-test between BRAF_like (group A) and RAS_like (group B).
    Returns DataFrame indexed by gene with log2FC, cohens_d, pval, and means.
    log2FC = mean(braf) - mean(ras)  (since input already log2).
    Positive => higher in BRAF_like.
    """
    if len(braf_cols) < 2 or len(ras_cols) < 2:
        return pd.DataFrame()
    A = expr[braf_cols].to_numpy(dtype=float)
    B = expr[ras_cols].to_numpy(dtype=float)
    meanA = np.nanmean(A, axis=1)
    meanB = np.nanmean(B, axis=1)
    varA = np.nanvar(A, axis=1, ddof=1)
    varB = np.nanvar(B, axis=1, ddof=1)
    nA = np.sum(~np.isnan(A), axis=1)
    nB = np.sum(~np.isnan(B), axis=1)
    # Welch t-test
    with np.errstate(divide="ignore", invalid="ignore"):
        se = np.sqrt(varA / np.maximum(nA, 1) + varB / np.maximum(nB, 1))
        t = (meanA - meanB) / se
        # degrees of freedom
        num = (varA / np.maximum(nA, 1) + varB / np.maximum(nB, 1)) ** 2
        den = (varA / np.maximum(nA, 1)) ** 2 / np.maximum(nA - 1, 1) \
              + (varB / np.maximum(nB, 1)) ** 2 / np.maximum(nB - 1, 1)
        df_t = num / den
    pvals = 2 * stats.t.sf(np.abs(t), df_t)
    pvals = np.where(np.isfinite(pvals), pvals, 1.0)
    # Cohen's d (pooled SD)
    pooled_sd = np.sqrt(((nA - 1) * varA + (nB - 1) * varB) / np.maximum(nA + nB - 2, 1))
    cohens_d = np.where(pooled_sd > 0, (meanA - meanB) / pooled_sd, 0.0)
    out = pd.DataFrame({
        "gene": expr.index,
        "log2FC": meanA - meanB,
        "cohens_d": cohens_d,
        "pval": pvals,
        "braf_mean": meanA,
        "ras_mean": meanB,
        "n_braf": nA,
        "n_ras": nB,
    }).set_index("gene")
    return out


def add_fdr(df: pd.DataFrame, col: str = "pval") -> pd.DataFrame:
    if len(df) == 0:
        df["fdr"] = []
        return df
    mask = df[col].notna()
    fdr = np.ones(len(df))
    if mask.sum() > 0:
        _, q, _, _ = multipletests(df.loc[mask, col].values, method="fdr_bh")
        fdr[mask.values] = q
    df["fdr"] = fdr
    return df


# ----------------------------------------------------------------------
# Main pipeline
# ----------------------------------------------------------------------
def main() -> int:
    print("[biomarker] loading sample master ...")
    sm = load_sample_master()

    print("[biomarker] loading TCGA-THCA expression ...")
    tcga = load_expression(DATA_TCGA)
    print(f"         shape = {tcga.shape}")

    print("[biomarker] loading GSE27155 microarray expression ...")
    gse27 = load_expression(DATA_GSE27155)
    print(f"         shape = {gse27.shape}")

    print("[biomarker] loading GSE126698 RNA-seq expression ...")
    gse126 = load_expression(DATA_GSE126698)
    print(f"         shape = {gse126.shape}")

    # ---- Step 2: TCGA DE BRAF_like vs RAS_like ----
    tcga_braf, tcga_ras = split_samples(sm, "TCGA-THCA", tcga.columns.tolist())
    print(f"[biomarker] TCGA BRAF_like n={len(tcga_braf)}  RAS_like n={len(tcga_ras)}")
    tcga_de = welch_de(tcga, tcga_braf, tcga_ras)
    tcga_de = add_fdr(tcga_de)
    print(f"[biomarker] TCGA DE genes computed: {len(tcga_de)}")

    # ---- Step 3: external replication ----
    # GSE27155: use sample_master molecular_subtype labels
    e27_braf, e27_ras = split_samples(sm, "GSE27155", gse27.columns.tolist())
    print(f"[biomarker] GSE27155 BRAF_like n={len(e27_braf)}  RAS_like n={len(e27_ras)}")
    de_27 = welch_de(gse27, e27_braf, e27_ras)
    de_27 = add_fdr(de_27) if len(de_27) else de_27

    # GSE126698: derive labels from column prefix (A/F/N/P)
    label_126 = derive_126698_labels(gse126.columns.tolist())
    c126_braf = [c for c, l in label_126.items() if l == "BRAF_like"]
    c126_ras = [c for c, l in label_126.items() if l == "RAS_like"]
    print(f"[biomarker] GSE126698 BRAF_like n={len(c126_braf)}  RAS_like n={len(c126_ras)}")
    de_126 = welch_de(gse126, c126_braf, c126_ras)
    de_126 = add_fdr(de_126) if len(de_126) else de_126

    # ---- Merge per-gene results ----
    de_all = tcga_de.rename(columns={
        "log2FC": "log2FC_tcga", "cohens_d": "cohens_d_tcga",
        "pval": "pval_tcga", "fdr": "fdr_tcga",
        "braf_mean": "braf_mean_tcga", "ras_mean": "ras_mean_tcga",
        "n_braf": "n_braf_tcga", "n_ras": "n_ras_tcga",
    })

    if len(de_27):
        de_all = de_all.join(de_27[["log2FC", "cohens_d", "pval", "fdr"]].rename(
            columns={"log2FC": "log2FC_27155", "cohens_d": "cohens_d_27155",
                     "pval": "pval_27155", "fdr": "fdr_27155"}), how="left")
    else:
        for c in ("log2FC_27155", "cohens_d_27155", "pval_27155", "fdr_27155"):
            de_all[c] = np.nan
    if len(de_126):
        de_all = de_all.join(de_126[["log2FC", "cohens_d", "pval", "fdr"]].rename(
            columns={"log2FC": "log2FC_126698", "cohens_d": "cohens_d_126698",
                     "pval": "pval_126698", "fdr": "fdr_126698"}), how="left")
    else:
        for c in ("log2FC_126698", "cohens_d_126698", "pval_126698", "fdr_126698"):
            de_all[c] = np.nan

    # Replication flags
    def _rep_flag(row, suffix):
        l2 = row.get(f"log2FC_{suffix}")
        p = row.get(f"pval_{suffix}")
        l_tcga = row.get("log2FC_tcga")
        if pd.isna(l2) or pd.isna(p) or pd.isna(l_tcga):
            return 0
        same_sign = np.sign(l2) == np.sign(l_tcga) and l_tcga != 0
        return int(bool(same_sign) and p < 0.05)

    de_all["replicated_27155"] = de_all.apply(lambda r: _rep_flag(r, "27155"), axis=1)
    de_all["replicated_126698"] = de_all.apply(lambda r: _rep_flag(r, "126698"), axis=1)
    de_all["replication_rate"] = de_all["replicated_27155"] + de_all["replicated_126698"]

    # Presence across cohorts (gene exists in matrix at all)
    present_tcga = set(tcga.index)
    present_27 = set(gse27.index)
    present_126 = set(gse126.index)
    de_all["present_tcga"] = [g in present_tcga for g in de_all.index]
    de_all["present_27155"] = [g in present_27 for g in de_all.index]
    de_all["present_126698"] = [g in present_126 for g in de_all.index]
    de_all["n_cohorts_present"] = (
        de_all["present_tcga"].astype(int)
        + de_all["present_27155"].astype(int)
        + de_all["present_126698"].astype(int)
    )

    # Median expression in TCGA (for noise-floor filter)
    de_all["median_expr_tcga"] = tcga.loc[de_all.index.intersection(tcga.index)].median(
        axis=1, numeric_only=True
    )
    de_all["median_expr_tcga"] = de_all["median_expr_tcga"].reindex(de_all.index)

    # Detection-rate coverage metric (fraction of cohorts where median>1)
    def _above_thr(df, gene, thr=1.0):
        if gene not in df.index:
            return False
        v = df.loc[gene].median()
        return bool(v > thr)

    cov = []
    for g in de_all.index:
        c = (
            int(_above_thr(tcga, g, 1.0))
            + int(_above_thr(gse27, g, 1.0))
            + int(_above_thr(gse126, g, 1.0))
        )
        cov.append(c / 3.0)
    de_all["coverage_fraction"] = cov

    # Known / novel classification
    de_all["is_known"] = [is_known(g) for g in de_all.index]
    de_all["panel_memberships"] = [
        ";".join(PANEL_MEMBERSHIP.get(g, [])) for g in de_all.index
    ]
    de_all["known_category"] = [
        KNOWN_BIOMARKERS.get(g, {}).get("category", "") for g in de_all.index
    ]
    de_all["known_citation"] = [
        KNOWN_BIOMARKERS.get(g, {}).get("citation", "") for g in de_all.index
    ]

    # Novel candidate filter
    novel_mask = (
        (~de_all["is_known"])
        & (de_all["fdr_tcga"] < 0.05)
        & (de_all["cohens_d_tcga"].abs() > 0.5)
        & (de_all["replication_rate"] >= 1)
        & (de_all["median_expr_tcga"] > 1)
        & (de_all["n_cohorts_present"] >= 2)
    )
    de_all["is_novel_validated"] = novel_mask.astype(bool)

    # Novelty score (for any gene; meaningful mostly for novels)
    with np.errstate(divide="ignore", invalid="ignore"):
        nlog = -np.log10(de_all["fdr_tcga"].clip(lower=1e-300))
    de_all["novelty_score"] = (
        de_all["cohens_d_tcga"].abs() * de_all["replication_rate"] * nlog
    ).fillna(0.0)

    # "Known-replicated" bookkeeping
    de_all["known_replicated"] = (
        de_all["is_known"]
        & (de_all["replication_rate"] >= 1)
        & (de_all["fdr_tcga"] < 0.05)
    )

    # ---- Write full DE table ----
    de_full_path = OUT_TABLES / "biomarker_de_full.tsv"
    de_all.reset_index().rename(columns={"index": "gene"}).to_csv(
        de_full_path, sep="\t", index=False
    )
    print(f"[biomarker] wrote {de_full_path}")

    # ---- known_vs_novel table (just the labels + top fields) ----
    knvn_cols = [
        "is_known", "is_novel_validated", "known_category", "panel_memberships",
        "log2FC_tcga", "cohens_d_tcga", "fdr_tcga",
        "log2FC_27155", "pval_27155", "log2FC_126698", "pval_126698",
        "replicated_27155", "replicated_126698", "replication_rate",
        "median_expr_tcga", "n_cohorts_present", "coverage_fraction",
        "novelty_score", "known_replicated",
    ]
    knvn = de_all[knvn_cols].copy()
    knvn_path = OUT_TABLES / "biomarker_known_vs_novel.tsv"
    knvn.reset_index().rename(columns={"index": "gene"}).to_csv(
        knvn_path, sep="\t", index=False
    )
    print(f"[biomarker] wrote {knvn_path}")

    # ---- validated table (known-replicated OR novel-validated) ----
    validated = de_all[
        de_all["known_replicated"] | de_all["is_novel_validated"]
    ].copy()
    validated = validated.sort_values("novelty_score", ascending=False)
    val_path = OUT_TABLES / "biomarker_validated.tsv"
    validated.reset_index().rename(columns={"index": "gene"}).to_csv(
        val_path, sep="\t", index=False
    )
    print(f"[biomarker] wrote {val_path} ({len(validated)} rows)")

    # ---- Step 5: Compact JSON payload ----
    def _safe(v):
        if pd.isna(v):
            return None
        if isinstance(v, (np.floating, np.integer)):
            v = float(v)
            if math.isnan(v) or math.isinf(v):
                return None
            return round(v, 4)
        if isinstance(v, (np.bool_, bool)):
            return bool(v)
        return v

    top_known = de_all[de_all["is_known"]].copy()
    top_known = top_known.sort_values("cohens_d_tcga", key=lambda s: s.abs(), ascending=False).head(200)
    top_novel = de_all[~de_all["is_known"]].copy()
    top_novel = top_novel.sort_values("novelty_score", ascending=False).head(200)

    def _rowdict(g, row):
        return {
            "gene": g,
            "log2FC": _safe(row["log2FC_tcga"]),
            "cohens_d": _safe(row["cohens_d_tcga"]),
            "fdr": _safe(row["fdr_tcga"]),
            "replicated": [int(row["replicated_27155"]), int(row["replicated_126698"])],
            "panel_memberships": PANEL_MEMBERSHIP.get(g, []),
            "is_novel": bool(row["is_novel_validated"]),
            "is_known": bool(row["is_known"]),
            "novelty_score": _safe(row["novelty_score"]),
            "per_cohort_means": [
                _safe(row["braf_mean_tcga"]),
                _safe(row["ras_mean_tcga"]),
            ],
            "external_log2FC": [_safe(row["log2FC_27155"]), _safe(row["log2FC_126698"])],
            "known_category": row["known_category"] or None,
        }

    payload = {
        "generated_by": "notebooks_or_scripts/biomarker_analysis.py",
        "cohorts": {
            "TCGA-THCA": {"n_braf": int(len(tcga_braf)), "n_ras": int(len(tcga_ras))},
            "GSE27155": {"n_braf": int(len(e27_braf)), "n_ras": int(len(e27_ras))},
            "GSE126698": {"n_braf": int(len(c126_braf)), "n_ras": int(len(c126_ras))},
        },
        "summary": {
            "n_genes_tested_tcga": int(len(tcga_de)),
            "n_known_total": int(de_all["is_known"].sum()),
            "n_known_replicated": int(de_all["known_replicated"].sum()),
            "n_novel_candidates": int(
                ((~de_all["is_known"]) & (de_all["fdr_tcga"] < 0.05)
                 & (de_all["cohens_d_tcga"].abs() > 0.5)).sum()
            ),
            "n_novel_validated": int(de_all["is_novel_validated"].sum()),
            "n_genes_in_2plus_cohorts": int((de_all["n_cohorts_present"] >= 2).sum()),
            "n_genes_in_all_3_cohorts": int((de_all["n_cohorts_present"] == 3).sum()),
        },
        "top_known": [_rowdict(g, r) for g, r in top_known.iterrows()],
        "top_novel": [_rowdict(g, r) for g, r in top_novel.iterrows()],
    }
    payload_path = OUT_DATA_ASSETS / "biomarker_payload.json"
    with payload_path.open("w") as f:
        json.dump(payload, f, separators=(",", ":"))
    size_kb = payload_path.stat().st_size / 1024
    print(f"[biomarker] wrote {payload_path} ({size_kb:.1f} KB)")

    # ---- Step 6: Plotly figures ----
    # 1. Volcano
    vol_df = de_all.copy()
    vol_df["status"] = "not_significant"
    sig_mask = (vol_df["fdr_tcga"] < 0.05) & (vol_df["cohens_d_tcga"].abs() > 0.5)
    vol_df.loc[sig_mask & vol_df["is_known"] & (vol_df["replication_rate"] >= 1), "status"] = "known_replicated"
    vol_df.loc[sig_mask & vol_df["is_known"] & (vol_df["replication_rate"] == 0), "status"] = "known_unreplicated"
    vol_df.loc[sig_mask & (~vol_df["is_known"]) & (vol_df["replication_rate"] >= 1), "status"] = "novel_validated"
    vol_df.loc[sig_mask & (~vol_df["is_known"]) & (vol_df["replication_rate"] == 0), "status"] = "novel_unvalidated"

    status_colors = {
        "known_replicated": "#1f77b4",
        "known_unreplicated": "#aec7e8",
        "novel_validated": "#d62728",
        "novel_unvalidated": "#ff9896",
        "not_significant": "#cccccc",
    }

    fig = go.Figure()
    vol_df["neglog10_fdr"] = -np.log10(vol_df["fdr_tcga"].clip(lower=1e-300))
    for status, color in status_colors.items():
        sub = vol_df[vol_df["status"] == status]
        fig.add_trace(go.Scattergl(
            x=sub["log2FC_tcga"], y=sub["neglog10_fdr"],
            mode="markers",
            name=f"{status} (n={len(sub)})",
            marker=dict(size=6, color=color, opacity=0.8 if status != "not_significant" else 0.3),
            hovertext=[
                f"{g}<br>log2FC={r['log2FC_tcga']:.2f}<br>d={r['cohens_d_tcga']:.2f}<br>"
                f"FDR={r['fdr_tcga']:.2e}<br>repl={int(r['replication_rate'])}/2<br>"
                f"panels={PANEL_MEMBERSHIP.get(g, [])}"
                for g, r in sub.iterrows()
            ],
            hoverinfo="text",
        ))
    # threshold lines
    fig.add_hline(y=-np.log10(0.05), line_dash="dash", line_color="grey",
                  annotation_text="FDR=0.05")
    fig.add_vline(x=1, line_dash="dash", line_color="grey")
    fig.add_vline(x=-1, line_dash="dash", line_color="grey")
    fig.update_layout(
        title="Volcano: BRAF_like vs RAS_like (TCGA-THCA, log2FC = BRAF - RAS)",
        xaxis_title="log2 fold change (BRAF_like − RAS_like)",
        yaxis_title="-log10(FDR)",
        template="plotly_white",
        hovermode="closest",
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
    )
    volcano_path = OUT_FIGS_INTERACTIVE / "volcano_braf_vs_ras.html"
    fig.write_html(volcano_path, include_plotlyjs="cdn")
    print(f"[biomarker] wrote {volcano_path}")

    # 2. Replication scatter: TCGA vs mean external
    ext_mean_logfc = vol_df[["log2FC_27155", "log2FC_126698"]].mean(axis=1, skipna=True)
    rep_df = vol_df[sig_mask].copy()
    rep_df["log2FC_external_mean"] = ext_mean_logfc.reindex(rep_df.index)
    rep_df = rep_df.dropna(subset=["log2FC_external_mean"])
    fig2 = go.Figure()
    for is_novel, color, label in [(True, "#d62728", "novel"), (False, "#1f77b4", "known")]:
        sub = rep_df[rep_df["is_novel_validated"] == is_novel] if is_novel else rep_df[rep_df["is_known"]]
        fig2.add_trace(go.Scattergl(
            x=sub["log2FC_tcga"],
            y=sub["log2FC_external_mean"],
            mode="markers",
            name=label + f" (n={len(sub)})",
            marker=dict(
                size=np.clip(sub["novelty_score"].values / max(1e-6, rep_df["novelty_score"].max()) * 14 + 4, 4, 20),
                color=color, opacity=0.7,
            ),
            hovertext=[
                f"{g}<br>TCGA log2FC={r['log2FC_tcga']:.2f}<br>"
                f"ext mean log2FC={r['log2FC_external_mean']:.2f}<br>"
                f"panels={PANEL_MEMBERSHIP.get(g, [])}<br>"
                f"novelty={r['novelty_score']:.2f}"
                for g, r in sub.iterrows()
            ],
            hoverinfo="text",
        ))
    # y=x line
    lo, hi = -6, 6
    if len(rep_df) > 0:
        lo = float(min(rep_df["log2FC_tcga"].min(), rep_df["log2FC_external_mean"].min()))
        hi = float(max(rep_df["log2FC_tcga"].max(), rep_df["log2FC_external_mean"].max()))
    fig2.add_trace(go.Scatter(
        x=[lo, hi], y=[lo, hi], mode="lines",
        line=dict(color="grey", dash="dash"), name="y=x",
    ))
    fig2.update_layout(
        title="Replication: TCGA effect vs mean external effect (GSE27155, GSE126698)",
        xaxis_title="log2FC (TCGA, BRAF−RAS)",
        yaxis_title="mean log2FC (GSE27155 + GSE126698)",
        template="plotly_white",
    )
    rep_path = OUT_FIGS_INTERACTIVE / "biomarker_replication_scatter.html"
    fig2.write_html(rep_path, include_plotlyjs="cdn")
    print(f"[biomarker] wrote {rep_path}")

    # 3. Heatmap: top 30 known + top 30 novel
    top30_known = (
        de_all[de_all["is_known"] & (de_all["fdr_tcga"] < 0.05)]
        .sort_values("cohens_d_tcga", key=lambda s: s.abs(), ascending=False)
        .head(30).index.tolist()
    )
    top30_novel = (
        de_all[de_all["is_novel_validated"]]
        .sort_values("novelty_score", ascending=False)
        .head(30).index.tolist()
    )
    heat_genes = [g for g in (top30_known + top30_novel) if g in tcga.index]

    # Use TCGA samples annotated as BRAF_like or RAS_like
    keep_cols = tcga_braf + tcga_ras
    keep_cols = [c for c in keep_cols if c in tcga.columns]
    hm_mat = tcga.loc[heat_genes, keep_cols].copy()
    # z-score within TCGA cohort (row-wise)
    mu = hm_mat.mean(axis=1)
    sd = hm_mat.std(axis=1).replace(0, 1)
    hm_z = hm_mat.sub(mu, axis=0).div(sd, axis=0)

    col_labels = ["BRAF" if c in set(tcga_braf) else "RAS" for c in keep_cols]
    # sort columns so BRAF first then RAS
    order = np.argsort([0 if l == "BRAF" else 1 for l in col_labels])
    hm_z = hm_z.iloc[:, order]
    col_labels = [col_labels[i] for i in order]
    row_labels = [f"{g} [{'known' if g in top30_known else 'novel'}]" for g in heat_genes]

    fig3 = go.Figure(data=go.Heatmap(
        z=hm_z.values,
        x=[f"{s}|{l}" for s, l in zip(hm_z.columns, col_labels)],
        y=row_labels,
        colorscale="RdBu_r",
        zmin=-3, zmax=3,
        colorbar=dict(title="z-score"),
    ))
    fig3.update_layout(
        title="Top 30 known + top 30 novel BRAF/RAS markers (TCGA; row z-score)",
        xaxis=dict(title="samples (BRAF_like → RAS_like)", tickangle=-90, showticklabels=False),
        yaxis=dict(title="gene [category]"),
        template="plotly_white",
        height=900,
    )
    heat_path = OUT_FIGS_INTERACTIVE / "biomarker_heatmap_known_vs_novel.html"
    fig3.write_html(heat_path, include_plotlyjs="cdn")
    print(f"[biomarker] wrote {heat_path} (z-scored row-wise within TCGA cohort — single-cohort view is cleanest because external cohort platforms/sample sizes differ)")

    # 4. Effect vs coverage
    eff_df = de_all.copy()
    eff_df["abs_d"] = eff_df["cohens_d_tcga"].abs()
    fig4 = go.Figure()
    for is_novel, color, label in [(True, "#d62728", "novel_validated"),
                                   (False, "#1f77b4", "known/other")]:
        sub = eff_df[eff_df["is_novel_validated"] == is_novel] if is_novel else eff_df[~eff_df["is_novel_validated"]]
        fig4.add_trace(go.Scattergl(
            x=sub["coverage_fraction"] + np.random.default_rng(7).normal(0, 0.01, len(sub)),
            y=sub["abs_d"],
            mode="markers",
            name=f"{label} (n={len(sub)})",
            marker=dict(size=5, color=color,
                         opacity=0.75 if is_novel else 0.25),
            hovertext=[
                f"{g}<br>|d|={r['abs_d']:.2f}<br>cov={r['coverage_fraction']:.2f}<br>"
                f"panels={PANEL_MEMBERSHIP.get(g, [])}"
                for g, r in sub.iterrows()
            ],
            hoverinfo="text",
        ))
    fig4.update_layout(
        title="Effect vs coverage: |Cohen's d| (TCGA) vs fraction of cohorts where median log2>1",
        xaxis_title="coverage fraction (0, 1/3, 2/3, 1) + small jitter",
        yaxis_title="|Cohen's d| (TCGA BRAF vs RAS)",
        template="plotly_white",
    )
    cov_path = OUT_FIGS_INTERACTIVE / "biomarker_effect_vs_coverage.html"
    fig4.write_html(cov_path, include_plotlyjs="cdn")
    print(f"[biomarker] wrote {cov_path}")

    # ---- Step 7: Markdown summary ----
    known_all = de_all[de_all["is_known"]].copy()
    known_all = known_all.sort_values("cohens_d_tcga", key=lambda s: s.abs(), ascending=False)
    top_known_md = known_all.head(20)

    novel_sorted = de_all[de_all["is_novel_validated"]].sort_values(
        "novelty_score", ascending=False
    )
    top_novel_md = novel_sorted.head(20)
    top5_novel = novel_sorted.head(5)
    top5_known_md = known_all[known_all["known_replicated"]].head(5)

    def _fmt_row_known(g, r):
        return (
            f"| {g} | {r['known_category'] or '-'} | {r['log2FC_tcga']:+.2f} | "
            f"{r['cohens_d_tcga']:+.2f} | {r['fdr_tcga']:.2e} | "
            f"{int(r['replicated_27155'])} | {int(r['replicated_126698'])} | "
            f"{int(r['replication_rate'])} |"
        )

    def _fmt_row_novel(g, r):
        return (
            f"| {g} | {r['log2FC_tcga']:+.2f} | {r['cohens_d_tcga']:+.2f} | "
            f"{r['fdr_tcga']:.2e} | {int(r['replicated_27155'])} | "
            f"{int(r['replicated_126698'])} | {r['novelty_score']:.2f} | "
            f"{r['median_expr_tcga']:.2f} |"
        )

    md_lines = []
    md_lines.append("# THCA Biomarker Analysis: Known vs Novel\n")
    md_lines.append("## TL;DR\n")

    top3_novel_short = ", ".join(
        f"**{g}** (d={r['cohens_d_tcga']:+.2f}, replicated {int(r['replication_rate'])}/2)"
        for g, r in novel_sorted.head(3).iterrows()
    )

    tldr = (
        f"Tested {len(tcga_de):,} genes on TCGA-THCA BRAF_like (n={len(tcga_braf)}) "
        f"vs RAS_like (n={len(tcga_ras)}). "
        f"{int(de_all['is_known'].sum())} of the {len(KNOWN_BIOMARKERS) + sum(len(v) for v in PANEL_SETS.values())} "
        f"hardcoded + panel reference genes were present; "
        f"{int(de_all['known_replicated'].sum())} replicated (same direction, external p<0.05) "
        f"in ≥1 external cohort. "
        f"{int(de_all['is_novel_validated'].sum())} genes passed the novel-validated filter "
        f"(FDR<0.05, |d|>0.5, replication≥1, median log2>1, present in ≥2 cohorts). "
        f"Top novel candidates: {top3_novel_short}. "
        f"Genes with large TCGA effects that fail to replicate in BOTH externals are **not** called biomarkers."
    )
    md_lines.append(tldr + "\n")

    md_lines.append("## Methods summary\n")
    md_lines.append(
        "- **Discovery**: Welch t-test on log2 expression between BRAF_like and RAS_like TCGA samples, BH-FDR across all tested genes. "
        "log2FC is defined as mean(BRAF) − mean(RAS).\n"
        "- **Replication**: same test re-run on GSE27155 microarray and GSE126698 RNA-seq. "
        "A gene replicates in a cohort if its external log2FC has the same sign as TCGA and the external p-value is <0.05 (per-cohort, uncorrected).\n"
        "- **External labels**: GSE27155 uses `molecular_subtype` from `sample_master.tsv` (cPTC/tall-cell→BRAF_like, FVPTC/FTC→RAS_like as already encoded by rerun_v2). "
        "GSE126698 sample codes are mapped by prefix: P*→BRAF_like (PTC), F*→RAS_like (FTC); A* (ATC) and N* (normal) are excluded.\n"
        "- **Known set**: hardcoded reference table (below) plus genes in `tds16_genes.txt`, `tierA67_genes.txt`, `brs71_genes.txt`.\n"
        "- **Novel filter**: NOT in any known list AND TCGA FDR<0.05 AND |d|>0.5 AND replication_rate≥1 AND median TCGA log2>1 AND present in ≥2 cohorts.\n"
        "- **Novelty score** = |d|_TCGA × replication_rate × −log10(FDR) (used for ranking, not as a hard threshold).\n"
    )

    md_lines.append("## Known biomarkers: what replicated, what didn't\n")
    md_lines.append(
        "| gene | category | log2FC_TCGA | d | FDR_TCGA | repl_GSE27155 | repl_GSE126698 | rep_rate |"
    )
    md_lines.append("|------|----------|-------------|---|----------|-----|-----|---|")
    for g, r in known_all.iterrows():
        md_lines.append(_fmt_row_known(g, r))
    md_lines.append("")

    md_lines.append("## Novel candidate biomarkers\n")
    if len(top_novel_md) == 0:
        md_lines.append("_No gene passed the novel-validated filter. See Caveats section._\n")
    else:
        md_lines.append(
            "| gene | log2FC_TCGA | d | FDR_TCGA | repl_GSE27155 | repl_GSE126698 | novelty_score | median_log2_TCGA |"
        )
        md_lines.append("|------|-------------|---|----------|-----|-----|---|---|")
        for g, r in top_novel_md.iterrows():
            md_lines.append(_fmt_row_novel(g, r))
        md_lines.append("")

        md_lines.append("### Top 5 novel candidates — short context\n")
        for g, r in top5_novel.iterrows():
            repstr = (
                "both externals" if r["replication_rate"] == 2
                else ("GSE27155 only" if r["replicated_27155"] else "GSE126698 only")
            )
            direction = "up in BRAF_like" if r["log2FC_tcga"] > 0 else "up in RAS_like"
            md_lines.append(
                f"- **{g}** — TCGA d={r['cohens_d_tcga']:+.2f}, log2FC={r['log2FC_tcga']:+.2f}, "
                f"FDR={r['fdr_tcga']:.2e}; replicates in {repstr}; direction: {direction}. "
                f"Gene symbol only, no literature claim in this pipeline."
            )
        md_lines.append("")

    md_lines.append("### Top 5 known markers with clearest signal (replicated)\n")
    if len(top5_known_md) == 0:
        md_lines.append("_No known marker passed replication + FDR<0.05._\n")
    else:
        for g, r in top5_known_md.iterrows():
            repstr = (
                "both externals" if r["replication_rate"] == 2
                else ("GSE27155 only" if r["replicated_27155"] else "GSE126698 only")
            )
            md_lines.append(
                f"- **{g}** ({r['known_category']}) — d={r['cohens_d_tcga']:+.2f}, "
                f"FDR={r['fdr_tcga']:.2e}, replicates in {repstr}."
            )
        md_lines.append("")

    md_lines.append("## Caveats\n")
    md_lines.append(
        "- **Label proxies in external cohorts**: GSE27155 subtype labels come from histology-derived "
        "mappings already baked into `sample_master.tsv`; GSE126698 labels are inferred from sample-ID prefix "
        "(P→PTC→BRAF_like, F→FTC→RAS_like). These are not molecularly-defined subtypes.\n"
        "- **Small N in GSE126698**: only ~6 BRAF-like and ~6 RAS-like tumors. Per-cohort p-values are low-power; "
        "replicate-in-126698 flags should be read as directional support, not independent statistical validation.\n"
        "- **Microarray coverage**: GSE27155 carries ~13k probes mapped to genes; ~60% of TCGA-tested genes are simply not present there, "
        "which mechanically caps `replication_rate` at 1 for many otherwise-valid candidates.\n"
        "- **No cross-cohort joint test**: FDR is computed per cohort. There is no joint meta-analysis FDR, so the "
        "novel-validated set is a conjunction of per-cohort filters, not a joint-model discovery.\n"
        "- **Novelty is operational, not literature-exhaustive**: 'novel' here strictly means 'not in any hardcoded reference panel or gene list "
        "used by this pipeline'. Many genes flagged novel likely have prior thyroid-cancer literature that simply is not encoded here.\n"
        "- A gene with |d|≈2 in TCGA that fails to replicate in BOTH externals is NOT reported as a biomarker; the novel-validated filter "
        "explicitly requires replication_rate≥1.\n"
    )

    md_path = OUT_REPORTS / "biomarker_analysis.md"
    md_path.write_text("\n".join(md_lines))
    print(f"[biomarker] wrote {md_path}")

    # ---- Small runtime summary ----
    print("\n[biomarker] === SUMMARY ===")
    print(f"  genes tested (TCGA)       : {len(tcga_de):,}")
    print(f"  known genes present       : {int(de_all['is_known'].sum())}")
    print(f"  known replicated          : {int(de_all['known_replicated'].sum())}")
    print(f"  novel candidates (pre-rep): "
          f"{int(((~de_all['is_known']) & (de_all['fdr_tcga']<0.05) & (de_all['cohens_d_tcga'].abs()>0.5)).sum())}")
    print(f"  novel validated           : {int(de_all['is_novel_validated'].sum())}")

    if len(top5_novel):
        print("  top-5 novel:")
        for g, r in top5_novel.iterrows():
            print(f"    - {g:<12} d={r['cohens_d_tcga']:+.2f}  fdr={r['fdr_tcga']:.2e}  "
                  f"rep={int(r['replication_rate'])}/2  logFC={r['log2FC_tcga']:+.2f}")
    if len(top5_known_md):
        print("  top-5 known (replicated):")
        for g, r in top5_known_md.iterrows():
            print(f"    - {g:<12} d={r['cohens_d_tcga']:+.2f}  fdr={r['fdr_tcga']:.2e}  "
                  f"rep={int(r['replication_rate'])}/2  cat={r['known_category']}")

    # ---- Verify outputs ----
    expected = [
        OUT_TABLES / "biomarker_known_vs_novel.tsv",
        OUT_TABLES / "biomarker_de_full.tsv",
        OUT_TABLES / "biomarker_validated.tsv",
        OUT_DATA_ASSETS / "biomarker_payload.json",
        OUT_FIGS_INTERACTIVE / "volcano_braf_vs_ras.html",
        OUT_FIGS_INTERACTIVE / "biomarker_heatmap_known_vs_novel.html",
        OUT_FIGS_INTERACTIVE / "biomarker_effect_vs_coverage.html",
        OUT_FIGS_INTERACTIVE / "biomarker_replication_scatter.html",
        OUT_REPORTS / "biomarker_analysis.md",
    ]
    print("\n[biomarker] output files:")
    for p in expected:
        ok = p.exists()
        size = p.stat().st_size if ok else 0
        print(f"    [{'OK' if ok else 'MISS'}] {p}  ({size/1024:.1f} KB)")

    return 0


if __name__ == "__main__":
    sys.exit(main())
