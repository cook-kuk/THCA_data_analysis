#!/usr/bin/env python3
from __future__ import annotations

import json
import os
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse, stats
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from v17_common import (
    BRS71,
    CACHE,
    EXPR_PATHS,
    FIG as V17FIG,
    LOG as V17LOG,
    META,
    RAI_GENES,
    ROOT,
    TAB as V17TAB,
    TDS16,
    TIERA67,
    align_dataset_expr_meta,
    bh_fdr,
    log_line,
    read_expr,
    signature_score,
)

RES = ROOT / "results" / "v17p3"
TAB = RES / "tables"
FIG = RES / "figs"
RPT = ROOT / "reports" / "v17p3"
LOG = RES / "orchestrator.log"
P2TAB = ROOT / "results" / "v17p2" / "tables"
P2FIG = ROOT / "results" / "v17p2" / "figs"
SCRNA = Path("/data/thca/scrna/raw/scrna_raw.h5ad")
MSIG_H = ROOT / "data_raw" / "msigdb_cache" / "MSigDB_Hallmark_2020.json"
MSIG_R = ROOT / "data_raw" / "msigdb_cache" / "Reactome_2022.txt"
for d in (TAB, FIG, RPT):
    d.mkdir(parents=True, exist_ok=True)

SEED = 42

MAPK_PANEL = ["DUSP4", "DUSP5", "DUSP6", "ETV4", "ETV5", "FOSL1", "SPRY1", "SPRY2", "SPRY4", "MET", "PHLDA1"]
THYROID_DIFF = ["PAX8", "NKX2-1", "FOXE1", "SLC5A5", "TPO", "TG", "TSHR", "DIO1", "DIO2", "GLIS3", "SLC26A4", "IYD"]
CYTOLYTIC = ["GZMA", "PRF1"]
IMMUNE_EVASION = ["CD274", "PDCD1", "CTLA4", "IDO1", "HLA-A", "HLA-B", "HLA-C", "B2M"]
TF_TARGETS = {
    "PAX8": ["TG", "TPO", "TSHR", "SLC5A5", "SLC26A4", "DIO1", "DIO2", "IYD"],
    "NKX2-1": ["TG", "TPO", "TSHR", "SLC5A5", "PAX8", "FOXE1", "DUOX1", "DUOX2"],
    "FOXE1": ["TG", "TPO", "PAX8", "NKX2-1", "SLC26A4", "IYD"],
    "HHEX": ["FOXE1", "PAX8", "NKX2-1", "TG", "TSHR"],
}


def write_json(path: Path, payload: dict) -> None:
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")


def load_dark() -> tuple[pd.DataFrame, pd.DataFrame]:
    meta = pd.read_csv(V17TAB / "dark_matter_cohort.tsv", sep="\t")
    expr = read_expr("TCGA-THCA", genes=None)
    expr = expr.loc[meta["sample_id"]].copy()
    expr.columns = expr.columns.astype(str)
    return expr, meta


def load_dark_tiera() -> tuple[pd.DataFrame, pd.DataFrame]:
    meta = pd.read_csv(V17TAB / "dark_matter_cohort.tsv", sep="\t")
    expr, _ = align_dataset_expr_meta("TCGA-THCA", genes=TIERA67, tumor_only=True)
    expr = expr.loc[meta["sample_id"]].copy()
    return expr, meta


def dm_labels(meta: pd.DataFrame) -> np.ndarray:
    return (meta["v17_dark_cluster"].astype(str) == "DM2").astype(int).to_numpy()


def preprocess_expr(expr: pd.DataFrame, genes: list[str]) -> pd.DataFrame:
    use = [g for g in genes if g in expr.columns]
    out = expr[use].apply(pd.to_numeric, errors="coerce")
    return out


def fit_dm_classifier(train_expr: pd.DataFrame, meta: pd.DataFrame, genes: list[str]) -> tuple[object, list[str]]:
    use = [g for g in genes if g in train_expr.columns]
    clf = make_pipeline(StandardScaler(), LogisticRegression(max_iter=4000, random_state=SEED))
    clf.fit(train_expr[use], dm_labels(meta))
    return clf, use


def dia_auc(y_true: np.ndarray, prob: np.ndarray) -> float:
    auc = roc_auc_score(y_true, prob)
    return float(max(auc, 1.0 - auc))


def score_expr(expr: pd.DataFrame, genes: list[str]) -> pd.Series:
    keep = [g for g in genes if g in expr.columns]
    if not keep:
        return pd.Series(np.nan, index=expr.index)
    return expr[keep].mean(axis=1)


def read_cross_cancer(cancer: str) -> tuple[np.ndarray, pd.Series, pd.Series, list[str]]:
    base = Path("/data/thca/data_processed/v5_cross_cancer") / cancer
    X = np.load(base / "X_combined.npz")["X"]
    y = pd.read_csv(base / "Y.tsv", sep="\t", header=None).iloc[:, 0].astype(str)
    b = pd.read_csv(base / "B.tsv", sep="\t", header=None).iloc[:, 0].astype(str)
    genes = [x.strip() for x in (base / "shared_genes.txt").read_text().splitlines() if x.strip()]
    n = min(X.shape[0], len(y), len(b))
    return X[:n], y.iloc[:n].reset_index(drop=True), b.iloc[:n].reset_index(drop=True), genes


def parse_reactome_txt(path: Path) -> dict[str, list[str]]:
    out: dict[str, list[str]] = {}
    for line in path.read_text().splitlines():
        parts = [p for p in line.split("\t") if p]
        if len(parts) >= 3:
            out[parts[0]] = parts[2:]
    return out


def parse_enrichr_text(text: str) -> dict[str, list[str]]:
    out = {}
    for line in text.splitlines():
        parts = [p for p in line.split("\t") if p]
        if len(parts) >= 2:
            out[parts[0]] = parts[1:]
    return out


def html_table(df: pd.DataFrame, table_id: str) -> str:
    head = "".join(f"<th>{c}</th>" for c in df.columns)
    body = []
    for _, row in df.iterrows():
        body.append("<tr>" + "".join(f"<td>{v}</td>" for v in row.tolist()) + "</tr>")
    return f"<table id='{table_id}' class='display compact stripe'><thead><tr>{head}</tr></thead><tbody>{''.join(body)}</tbody></table>"

