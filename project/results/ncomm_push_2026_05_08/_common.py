"""Shared utilities for the NComm-push public-data analyses (2026-05-08).

8-gene RAI panel + scoring helpers identical to Paper 1 v8 STAR Methods:
  z-score within sample, average across the 8 genes.
"""
from __future__ import annotations

import json
import os
import time
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

# 8-gene panel — Paper 1 v8 canonical
PANEL_8 = ["SLC5A5", "TPO", "TG", "TSHR", "PAX8", "NKX2-1", "FOXE1", "DIO1"]
PANEL_8_ENTREZ = {
    "SLC5A5": 6528,
    "TPO": 7173,
    "TG": 7038,
    "TSHR": 7253,
    "PAX8": 7849,
    "NKX2-1": 7080,
    "FOXE1": 2304,
    "DIO1": 1733,
}
# Some studies use NKX2-1 as TITF1
ALIASES = {"NKX2-1": ["NKX2-1", "TITF1", "NKX2.1", "NKX2_1"]}

ROOT = Path(__file__).resolve().parent
LOGS = ROOT / "logs"
LOGS.mkdir(exist_ok=True)


def log(name: str, msg: str) -> None:
    """Append a timestamped line to the worker's log."""
    ts = time.strftime("%Y-%m-%d %H:%M:%S")
    line = f"[{ts}] {msg}\n"
    print(line, end="")
    with (LOGS / f"{name}.log").open("a", encoding="utf-8") as fh:
        fh.write(line)


def write_status(worker: str, status: str, **kw) -> None:
    """Drop a JSON status file the launcher can grep for done/error."""
    payload = {"worker": worker, "status": status, "timestamp": time.time(), **kw}
    (ROOT / worker / "_status.json").write_text(json.dumps(payload, indent=2))


def resolve_panel(gene_index: Iterable[str]) -> dict[str, str]:
    """Map canonical PANEL_8 names to whatever symbol the study uses."""
    avail = {g.upper(): g for g in gene_index}
    out: dict[str, str] = {}
    for canonical in PANEL_8:
        candidates = ALIASES.get(canonical, [canonical])
        for c in candidates:
            if c.upper() in avail:
                out[canonical] = avail[c.upper()]
                break
    return out


def score_panel(expr: pd.DataFrame, gene_map: dict[str, str]) -> pd.DataFrame:
    """expr: rows=samples, cols=genes (after gene_map applied).

    Returns DataFrame with columns RAI_8 (= mean z-score), DM1_like (= -RAI_8).
    Within-sample z-score is canonical (cf. STAR Methods 'within-sample-centered profile').
    """
    sub = expr.loc[:, list(gene_map.values())].copy()
    # within-sample center then panel-z (per-gene cohort z, then mean across panel)
    cohort_z = (sub - sub.mean(axis=0)) / sub.std(axis=0).replace(0, 1)
    rai = cohort_z.mean(axis=1)
    out = pd.DataFrame({"RAI_8": rai, "DM1_like": -rai}, index=expr.index)
    return out


def kmeans_dm(scores: pd.DataFrame, seed: int = 42) -> pd.Series:
    """KMeans k=2 on the panel z-scores. DM1 = lower-RAI cluster."""
    from sklearn.cluster import KMeans

    km = KMeans(n_clusters=2, n_init=10, random_state=seed).fit(scores[["RAI_8"]])
    labels = km.labels_
    # DM1 = the cluster with the LOWER mean RAI_8
    centers = km.cluster_centers_.flatten()
    dm1_label = int(np.argmin(centers))
    out = pd.Series(
        np.where(labels == dm1_label, "DM1", "DM2"),
        index=scores.index,
        name="DM_call",
    )
    return out


def set_threads(n: int) -> None:
    """Limit BLAS / OMP threads so 4 workers in parallel don't fight."""
    for var in (
        "OMP_NUM_THREADS",
        "MKL_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        os.environ[var] = str(n)
