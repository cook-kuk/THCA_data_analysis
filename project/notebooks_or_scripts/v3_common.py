"""Shared utilities for v3_* scripts.

Loads cohort expression matrices on tumor-only samples with BRAF_like/RAS_like
labels (and dedifferentiated where available), reusing metadata/sample_master.tsv.
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

PROJECT = Path("/opt/thyroid-dash/project")
DATA_PROC = PROJECT / "data_processed"
META = PROJECT / "metadata"
RESULTS = PROJECT / "results"
RESULTS_TABLES = RESULTS / "tables"
REPORTS = PROJECT / "reports"
HTML = REPORTS / "html"
FIGS = HTML / "figs_interactive"
PAGES = HTML / "pages"
LOGS = PROJECT / "logs"

for d in (RESULTS_TABLES, FIGS, PAGES, LOGS, META):
    d.mkdir(parents=True, exist_ok=True)

# Ensure we can import rerun_v2
sys.path.insert(0, str(PROJECT / "notebooks_or_scripts"))

TCGA_RNA = DATA_PROC / "bulk_rnaseq" / "TCGA-THCA_rnaseq_expression_log2.tsv"
GSE126698_RNA = DATA_PROC / "bulk_rnaseq" / "GSE126698_rnaseq_expression_log2.tsv"
GSE213647_RNA = DATA_PROC / "bulk_rnaseq" / "GSE213647_rnaseq_expression_log2.tsv"
GSE27155_MA = DATA_PROC / "microarray" / "GSE27155_microarray_expression_log2.tsv"
GSE76039_MA = DATA_PROC / "microarray" / "GSE76039_microarray_expression_log2.tsv"
GSE97466_METH = DATA_PROC / "methylation" / "GSE97466_beta_top5000.tsv"

DATASET_FILES: Dict[str, Path] = {
    "TCGA-THCA": TCGA_RNA,
    "GSE27155": GSE27155_MA,
    "GSE126698": GSE126698_RNA,
    "GSE213647": GSE213647_RNA,
    "GSE76039": GSE76039_MA,
}

MAPK_OUTPUT = ["DUSP4", "DUSP5", "DUSP6", "SPRY1", "SPRY2", "SPRY4",
               "ETV4", "ETV5", "FOSL1", "PHLDA1"]


def get_tierA67_clean() -> List[str]:
    from rerun_v2 import TIERA67_CLEAN_UNIQUE
    return list(TIERA67_CLEAN_UNIQUE)


def get_tds16() -> List[str]:
    from rerun_v2 import TDS16
    return list(TDS16)


def load_expr(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    gene_col = df.columns[0]
    df = df.rename(columns={gene_col: "gene_symbol"}).set_index("gene_symbol")
    df = df[~df.index.duplicated(keep="first")]
    return df


def load_sample_master() -> pd.DataFrame:
    return pd.read_csv(META / "sample_master.tsv", sep="\t")


def _expr_sample_id_column(dataset: str, sm: pd.DataFrame) -> pd.Series:
    """Return the column of sample_master to use when matching to expression
    matrix columns.  GSE126698 stores the TSV alias (e.g. 'A3') in
    clinical_subtype_tag as '<alias> [RNA-Seq]'."""
    if dataset == "GSE126698":
        return sm["clinical_subtype_tag"].astype(str).str.replace(
            " [RNA-Seq]", "", regex=False).str.strip()
    return sm["sample_id"]


def build_feature_matrix(dataset: str, gene_list: List[str],
                         keep_tumor_only: bool = True,
                         labels: Tuple[str, ...] = ("BRAF_like", "RAS_like"),
                         ) -> Tuple[pd.DataFrame, pd.Series]:
    """Return (X[samples x genes], y) where y=1 for first label, 0 for second.

    Missing genes are filled with 0.0 (so cross-platform feature union is safe).
    """
    path = DATASET_FILES[dataset]
    expr = load_expr(path)
    sm = load_sample_master()
    sm = sm[sm["dataset"] == dataset].copy()
    if keep_tumor_only:
        sm = sm[sm["normal_vs_tumor"] == "tumor"]
    sm = sm[sm["molecular_subtype"].isin(labels)]
    sm["_expr_id"] = _expr_sample_id_column(dataset, sm)
    sm = sm[sm["_expr_id"].isin(expr.columns)].copy()
    y = sm["molecular_subtype"].map({labels[0]: 1, labels[1]: 0}).astype(int)
    y.index = sm["_expr_id"].values
    X = expr.reindex(gene_list).fillna(0.0).T.loc[sm["_expr_id"]]
    X.columns = list(gene_list)
    return X, y


def build_feature_matrix_multi(dataset: str, gene_list: List[str],
                               label_col: str = "molecular_subtype",
                               allowed_labels: Tuple[str, ...] = (
                                   "BRAF_like", "RAS_like", "dedifferentiated"),
                               keep_tumor_only: bool = True
                               ) -> Tuple[pd.DataFrame, pd.Series]:
    path = DATASET_FILES[dataset]
    expr = load_expr(path)
    sm = load_sample_master()
    sm = sm[sm["dataset"] == dataset].copy()
    if keep_tumor_only:
        sm = sm[sm["normal_vs_tumor"] == "tumor"]
    sm = sm[sm[label_col].isin(allowed_labels)]
    sm["_expr_id"] = _expr_sample_id_column(dataset, sm)
    sm = sm[sm["_expr_id"].isin(expr.columns)].copy()
    y = sm[label_col].astype(str)
    y.index = sm["_expr_id"].values
    X = expr.reindex(gene_list).fillna(0.0).T.loc[sm["_expr_id"]]
    X.columns = list(gene_list)
    return X, y


def write_tsv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)


def write_skipped(path: Path, reason: str, extra: dict | None = None) -> None:
    payload = {"status": "skipped_no_raw", "reason": reason}
    if extra:
        payload.update(extra)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        json.dump(payload, f, indent=2)


def bootstrap_ci(values: np.ndarray, n: int = 1000, alpha: float = 0.05,
                 rng_seed: int = 42) -> Tuple[float, float]:
    rng = np.random.default_rng(rng_seed)
    vals = np.asarray(values, dtype=float)
    if vals.size == 0 or np.all(np.isnan(vals)):
        return (float("nan"), float("nan"))
    boot = np.empty(n, dtype=float)
    for i in range(n):
        idx = rng.integers(0, vals.size, vals.size)
        boot[i] = np.nanmean(vals[idx])
    lo = float(np.nanpercentile(boot, 100 * alpha / 2))
    hi = float(np.nanpercentile(boot, 100 * (1 - alpha / 2)))
    return lo, hi


def small_n_warning(n: int, threshold: int = 20) -> str:
    return " ⚠ small n" if n < threshold else ""


def tumor_label_map(allow_dediff: bool = False):
    if allow_dediff:
        return {"BRAF_like": "BRAF_like", "RAS_like": "RAS_like",
                "dedifferentiated": "dedifferentiated"}
    return {"BRAF_like": 1, "RAS_like": 0}
