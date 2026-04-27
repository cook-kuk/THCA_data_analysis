#!/usr/bin/env python3
from __future__ import annotations

import json
import math
import os
import time
from pathlib import Path
from typing import Iterable

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
    EXPR_PATHS,
    EXTRA_DRIVER_GENES,
    META,
    RAI_GENES,
    ROOT,
    TDS16,
    TIERA67,
    TCGA_CLINICAL_EXT,
    align_dataset_expr_meta,
    bh_fdr,
    load_sample_master,
    log_line,
    read_expr,
    signature_score,
)

RES = ROOT / "results" / "v17p2"
TAB = RES / "tables"
FIG = RES / "figs"
RPT = ROOT / "reports" / "v17p2"
LOG = RES / "orchestrator.log"
PHASE1_TAB = ROOT / "results" / "v17" / "tables"
PHASE1_FIG = ROOT / "results" / "v17" / "figs"
SCRNA_H5AD = Path("/data/thca/scrna/raw/scrna_raw.h5ad")
V5_CC = Path("/data/thca/data_processed/v5_cross_cancer")
for _d in (TAB, FIG, RPT):
    _d.mkdir(parents=True, exist_ok=True)

SEED = 42


def getenv_int(name: str, default: int) -> int:
    try:
        return int(os.environ.get(name, default))
    except Exception:
        return default


def zscore_pca(expr: pd.DataFrame, n_components: int = 10) -> tuple[np.ndarray, StandardScaler, PCA]:
    X = expr.fillna(expr.mean()).to_numpy(float)
    scaler = StandardScaler()
    Xz = scaler.fit_transform(X)
    pca = PCA(n_components=min(n_components, Xz.shape[1], Xz.shape[0] - 1), random_state=SEED)
    pcs = pca.fit_transform(Xz)
    return pcs, scaler, pca


def load_dark_cohort() -> tuple[pd.DataFrame, pd.DataFrame]:
    meta = pd.read_csv(PHASE1_TAB / "dark_matter_cohort.tsv", sep="\t")
    expr, _ = align_dataset_expr_meta("TCGA-THCA", genes=TIERA67, tumor_only=True)
    expr = expr.loc[meta["sample_id"]].copy()
    return expr, meta


def load_dark_transcriptome() -> tuple[pd.DataFrame, pd.DataFrame]:
    meta = pd.read_csv(PHASE1_TAB / "dark_matter_cohort.tsv", sep="\t")
    expr = read_expr("TCGA-THCA", genes=None)
    expr = expr.loc[meta["sample_id"]].copy()
    return expr, meta


def dm_binary_labels(meta: pd.DataFrame) -> np.ndarray:
    return (meta["v17_dark_cluster"].astype(str) == "DM2").astype(int).to_numpy()


def tcga12(sample_id: str) -> str:
    sample_id = str(sample_id)
    return sample_id[:12] if sample_id.startswith("TCGA-") else sample_id


def load_clinical() -> pd.DataFrame:
    df = pd.read_csv(TCGA_CLINICAL_EXT, sep="\t")
    df["tcga12"] = df["sample_id"].astype(str).map(tcga12)
    return df


def stage_to_ord(x: str) -> float:
    if pd.isna(x):
        return np.nan
    x = str(x).strip().upper().replace("STAGE ", "")
    mapping = {"0A": 0.0, "I": 1.0, "II": 2.0, "III": 3.0, "IVA": 4.0, "IVB": 4.5, "IVC": 5.0, "IV": 4.0}
    return mapping.get(x, np.nan)


def fisher_2x2(a: int, b: int, c: int, d: int) -> float:
    try:
        return float(stats.fisher_exact([[a, b], [c, d]])[1])
    except Exception:
        return np.nan


def dia_auc(y_true: np.ndarray, prob: np.ndarray) -> float:
    auc = roc_auc_score(y_true, prob)
    return float(max(auc, 1.0 - auc))


def identifiability_auc(expr: pd.DataFrame, batches: np.ndarray) -> float:
    y = pd.Series(batches).astype(str).factorize()[0]
    if len(np.unique(y)) < 2:
        return np.nan
    clf = make_pipeline(StandardScaler(), LogisticRegression(max_iter=4000, random_state=SEED))
    cv = StratifiedKFold(5, shuffle=True, random_state=SEED)
    prob = cross_val_predict(clf, expr.to_numpy(float), y, cv=cv, method="predict_proba")[:, 1]
    return float(roc_auc_score(y, prob))


def ensure_overlap(tcga_expr: pd.DataFrame, ext_expr: pd.DataFrame, genes: Iterable[str] | None = None) -> list[str]:
    if genes is None:
        genes = tcga_expr.columns
    return [g for g in genes if g in tcga_expr.columns and g in ext_expr.columns]


def top_marker_sets() -> dict[str, list[str]]:
    markers = pd.read_csv(PHASE1_TAB / "dark_matter_cluster_markers.tsv", sep="\t")
    out: dict[str, list[str]] = {}
    for cluster in sorted(markers["cluster"].unique()):
        key = f"DM{int(cluster) + 1}"
        sub = markers[markers["cluster"] == cluster].sort_values(["fdr", "log2fc"], ascending=[True, False])
        pos = sub[sub["log2fc"] > 0]["gene"].tolist()[:20]
        out[key] = pos
    return out


def sparse_gene_score(X, idx: list[int]) -> np.ndarray:
    if not idx:
        return np.full(X.shape[0], np.nan)
    sub = X[:, idx]
    if sparse.issparse(sub):
        return np.asarray(sub.mean(axis=1)).ravel()
    return np.asarray(sub).mean(axis=1)


def write_json(path: Path, payload: dict) -> None:
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")


def html_table(df: pd.DataFrame, table_id: str) -> str:
    cols = "".join(f"<th>{c}</th>" for c in df.columns)
    rows = []
    for _, rec in df.iterrows():
        rows.append("<tr>" + "".join(f"<td>{v}</td>" for v in rec.tolist()) + "</tr>")
    return f"<table id='{table_id}' class='display compact stripe'><thead><tr>{cols}</tr></thead><tbody>{''.join(rows)}</tbody></table>"

