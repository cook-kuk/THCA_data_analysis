#!/usr/bin/env python3
"""Shared helpers for v5 sprint scripts.

All v5 scripts import labels, cohort loaders, and basic utilities
from here. Prefixed with `v5_` per sprint rules; reads v4/v3 sample
master to stay compatible with existing labeling.
"""
from __future__ import annotations

import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

ROOT = Path("/opt/thyroid-dash/project")
META = ROOT / "metadata"
BULK_ORIG = ROOT / "data_processed" / "bulk_rnaseq"
MICRO_ORIG = ROOT / "data_processed" / "microarray"
MICRO_V3 = ROOT / "data_processed" / "microarray_v3"
BULK_V3 = ROOT / "data_processed" / "bulk_rnaseq_v3"
RESULTS_ML = ROOT / "results" / "ml"
RESULTS_TABLES = ROOT / "results" / "tables"
REPORTS_HTML = ROOT / "reports" / "html"
FIGS_INTERACTIVE = REPORTS_HTML / "figs_interactive"
PAGES = REPORTS_HTML / "pages"
ASSETS_DATA = REPORTS_HTML / "assets" / "data"
REPORTS = ROOT / "reports"
LOGS = ROOT / "logs"
DATA_RAW = ROOT / "data_raw"

for p in [RESULTS_ML, RESULTS_TABLES, FIGS_INTERACTIVE, PAGES, ASSETS_DATA, LOGS, REPORTS]:
    p.mkdir(parents=True, exist_ok=True)


COHORTS = {
    "TCGA-THCA": BULK_ORIG / "TCGA-THCA_rnaseq_expression_log2.tsv",
    "GSE213647": BULK_ORIG / "GSE213647_rnaseq_expression_log2.tsv",
    "GSE126698": BULK_ORIG / "GSE126698_rnaseq_expression_log2.tsv",
    "GSE27155":  MICRO_ORIG / "GSE27155_microarray_expression_log2.tsv",
    "GSE33630":  MICRO_V3 / "GSE33630_v3_log2.tsv",
    "GSE76039":  MICRO_ORIG / "GSE76039_microarray_expression_log2.tsv",
}


def load_sample_master() -> pd.DataFrame:
    sm = pd.read_csv(META / "sample_master_v3_merged.tsv", sep="\t")
    sm["sample_id"] = sm["sample_id"].astype(str)
    return sm


def binary_label(row: pd.Series):
    nvt = str(row.get("normal_vs_tumor", "")).strip().lower()
    if nvt == "normal":
        return 0
    if nvt == "tumor":
        return 1
    return None


def three_class(row: pd.Series):
    nvt = str(row.get("normal_vs_tumor", "")).strip().lower()
    if nvt == "normal":
        return "normal"
    if nvt != "tumor":
        return None
    hist = str(row.get("histology_subtype", "")).strip()
    if hist in {"ATC", "PDTC"}:
        return "aggressive"
    if str(row.get("aggressive_flag", "")).lower() == "yes":
        return "aggressive"
    return "indolent"


def bio_label_bvr(row: pd.Series):
    """Biomarker label: 1=BRAF_like, 0=RAS_like, None otherwise.

    Uses molecular_subtype column if present (BRAF_like, RAS_like),
    falls back to driver_anchor for TCGA.
    """
    ms = str(row.get("molecular_subtype", "")).strip().lower()
    if "braf" in ms:
        return 1
    if "ras" in ms:
        return 0
    da = str(row.get("driver_anchor", "")).strip().lower()
    if "braf" in da:
        return 1
    if "ras" in da:
        return 0
    return None


def load_cohort(name: str, path: Path):
    if not path.exists():
        return None
    df = pd.read_csv(path, sep="\t", index_col=0)
    df.index = df.index.astype(str).str.upper()
    df = df.loc[~df.index.duplicated(keep="first")]
    ensg_frac = float(sum(str(x).startswith("ENSG") for x in df.index)) / max(len(df.index), 1)
    if ensg_frac > 0.5:
        return None
    return df


def build_pooled_matrix(min_samples_per_cohort: int = 5, top_var_genes: int = 3000,
                        label_fn=binary_label):
    """Build a shared-gene matrix across COHORTS with a label vector.

    Returns: (expr_df [genes x samples], meta_df [samples -> cohort,label])
    """
    sm_idx = load_sample_master().set_index("sample_id")

    per_cohort = {}
    for name, path in COHORTS.items():
        df = load_cohort(name, path)
        if df is None:
            continue
        per_cohort[name] = df

    if not per_cohort:
        raise RuntimeError("no cohorts loaded")

    shared = None
    for df in per_cohort.values():
        shared = set(df.index) if shared is None else shared & set(df.index)
    shared = sorted(shared)

    frames = []
    meta_rows = []
    for name, df in per_cohort.items():
        sub = df.loc[shared]
        keep_cols = [c for c in sub.columns if c in sm_idx.index]
        sub = sub[keep_cols]
        for c in keep_cols:
            row = sm_idx.loc[c]
            if isinstance(row, pd.DataFrame):
                row = row.iloc[0]
            lbl = label_fn(row)
            if lbl is None:
                continue
            meta_rows.append({"sample_id": c, "cohort": name, "label": lbl})
        frames.append(sub)

    meta_df = pd.DataFrame(meta_rows).drop_duplicates("sample_id").set_index("sample_id")
    big = pd.concat(frames, axis=1)
    big = big.loc[:, ~big.columns.duplicated()]
    keep = [s for s in meta_df.index if s in big.columns]
    big = big[keep]
    meta_df = meta_df.loc[keep]

    # drop cohorts with < min_samples
    for c in list(meta_df["cohort"].unique()):
        n = int((meta_df["cohort"] == c).sum())
        if n < min_samples_per_cohort:
            keep_rows = meta_df["cohort"] != c
            meta_df = meta_df.loc[keep_rows]
            big = big.loc[:, meta_df.index]

    # filter low-variance genes and NaN
    gene_var = big.var(axis=1)
    big = big.loc[gene_var > 1e-6]
    row_means = big.mean(axis=1)
    big = big.T.fillna(row_means).T.dropna(axis=0)

    if top_var_genes is not None and top_var_genes < big.shape[0]:
        top = big.var(axis=1).sort_values(ascending=False).head(top_var_genes).index
        big = big.loc[top]

    return big, meta_df


def annotate_small_n(ax, counts, threshold=20):
    """Add warning text to a matplotlib axis if any class has n<threshold."""
    import matplotlib.pyplot as plt  # noqa
    small = [k for k, v in dict(counts).items() if v < threshold]
    if small:
        ax.text(0.02, 0.98, f"! small n: {', '.join(map(str, small))}",
                transform=ax.transAxes, fontsize=8, va="top",
                color="#c0392b", alpha=0.9)


def safe_write_json(path: Path, obj) -> None:
    path.write_text(json.dumps(obj, indent=2, default=str))
