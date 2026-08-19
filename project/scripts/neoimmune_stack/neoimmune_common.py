#!/usr/bin/env python3
"""Shared utilities for NeoImmune-Stack / CLEAN-Neo++.

The scripts in this folder intentionally integrate existing result artifacts.
They do not train a new foundation model and they keep clean-science and
production-stack tracks separated.
"""

from __future__ import annotations

import gzip
import json
import math
import os
import re
from dataclasses import dataclass
from datetime import date
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[3]
TODAY = date.today().isoformat().replace("-", "_")
DEFAULT_RUN_NAME = f"neoimmune_stack_{TODAY}"
RESULTS_ROOT = REPO_ROOT / "project" / "results"
HUB_ROOT = Path("/data/neoantigen_vaccine_hub")


DIRS = [
    "audit",
    "registry",
    "data",
    "external_scores",
    "local_scores",
    "clean_track",
    "production_track",
    "metrics",
    "predictions",
    "figures",
    "reports",
]


CANONICAL_COLUMNS = [
    "candidate_id",
    "patient_id",
    "sample_id",
    "dataset_source",
    "tumor_type",
    "gene",
    "mutation_id",
    "source_type",
    "peptide_mut",
    "peptide_wt",
    "peptide_length",
    "flank_left",
    "flank_right",
    "hla_allele",
    "hla_class",
    "expression_tpm",
    "vaf",
    "clonality",
    "hla_loh_status",
    "apm_score",
    "b2m_status",
    "tap1_expr",
    "tap2_expr",
    "label_presentation",
    "label_immunogenicity",
    "assay_type",
    "assay_result",
    "validation_level",
    "train_test_group",
    "source_study",
]


EXTERNAL_SCORE_COLUMNS = [
    "candidate_id",
    "model_name",
    "score_type",
    "raw_score",
    "normalized_score",
    "rank_percentile",
    "higher_is_better",
    "runnable_status",
    "error_message",
    "runtime_sec",
]


LOCAL_SCORE_COLUMNS = [
    "candidate_id",
    "local_model_name",
    "local_score",
    "normalized_score",
    "fold_id",
    "split_type",
    "training_scope",
    "notes",
]


EXTERNAL_MODELS = [
    "NetMHCpan_4.1_EL",
    "NetMHCpan_4.1_BA",
    "MHCflurry_2.0_presentation",
    "MHCflurry_2.0_affinity",
    "BigMHC_EL",
    "BigMHC_IM",
    "PRIME",
    "PRIME2.1",
    "HLApollo",
    "pVACtools_parser",
    "NeoDisc_parser",
    "Topiary_pending",
    "MixMHCpred_pending",
    "NetMHCIIpan_pending",
]


LOCAL_MODELS = [
    "Structure_LR",
    "Wave8_TCR_SelfSim_full",
    "ESM2_Bayesian",
    "GP_quantum",
    "VQC",
    "W7A_QK_only",
    "W7A_full",
    "W7B_stacked",
    "quantum_kernel_no_anchor_gamma1.0",
    "quantum_kernel_compact_gamma0.5",
    "Stack_mean",
    "Stack_median",
    "Stack_LR",
    "BAR-Neo",
    "BAR-Neo-BMA",
    "BAR-Neo-X",
    "stress_guarded_clean_stack",
]


def ensure_run_dir(run_dir: str | Path | None = None) -> Path:
    if run_dir is None:
        run_dir = RESULTS_ROOT / DEFAULT_RUN_NAME
    run_dir = Path(run_dir)
    for d in DIRS:
        (run_dir / d).mkdir(parents=True, exist_ok=True)
    return run_dir


def rel(path: Path | str) -> str:
    p = Path(path)
    try:
        return str(p.relative_to(REPO_ROOT))
    except Exception:
        return str(p)


def safe_read_table(path: Path | str, nrows: int | None = None) -> pd.DataFrame:
    path = Path(path)
    sep = "\t" if path.suffix in {".tsv", ".gz"} and ".tsv" in path.name else ","
    return pd.read_csv(path, sep=sep, nrows=nrows, low_memory=False)


def write_tsv_gz(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False, compression="gzip")


def write_tsv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)


def write_md(text: str, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def norm_series(s: pd.Series, higher_is_better: bool = True) -> pd.Series:
    x = pd.to_numeric(s, errors="coerce")
    if x.notna().sum() == 0:
        return pd.Series(np.nan, index=s.index)
    lo = x.min(skipna=True)
    hi = x.max(skipna=True)
    if not np.isfinite(lo) or not np.isfinite(hi) or math.isclose(float(lo), float(hi)):
        out = pd.Series(0.5, index=s.index)
    else:
        out = (x - lo) / (hi - lo)
    if not higher_is_better:
        out = 1 - out
    return out.clip(0, 1)


def rank_percentile(score: pd.Series, higher_is_better: bool = True) -> pd.Series:
    x = pd.to_numeric(score, errors="coerce")
    return x.rank(pct=True, ascending=not higher_is_better)


def first_existing(paths: Iterable[Path | str]) -> Path | None:
    for p in paths:
        pp = Path(p)
        if pp.exists():
            return pp
    return None


def find_files(patterns: Iterable[str], roots: Iterable[Path | str]) -> list[Path]:
    out: list[Path] = []
    for root in roots:
        root = Path(root)
        if not root.exists():
            continue
        for pat in patterns:
            out.extend(root.glob(pat))
    return sorted(set(out))


def clean_hla(x: object) -> str:
    if pd.isna(x):
        return ""
    s = str(x).strip()
    if not s:
        return ""
    s = s.replace("HLA-", "")
    if s.startswith(("A*", "B*", "C*", "D")):
        return "HLA-" + s
    return s


def make_key(peptide: object, hla: object, source: object | None = None) -> str:
    pep = "" if pd.isna(peptide) else str(peptide).strip()
    allele = clean_hla(hla)
    src = "" if source is None or pd.isna(source) else str(source).strip()
    return f"{pep}|{allele}|{src}"


def build_candidate_lookup(canon: pd.DataFrame) -> dict[str, str]:
    lookup: dict[str, str] = {}
    for _, r in canon.iterrows():
        cid = str(r.get("candidate_id", ""))
        keys = [
            make_key(r.get("peptide_mut"), r.get("hla_allele"), r.get("dataset_source")),
            make_key(r.get("peptide_mut"), r.get("hla_allele"), r.get("source_study")),
            make_key(r.get("peptide_mut"), r.get("hla_allele"), None),
        ]
        for k in keys:
            lookup.setdefault(k, cid)
    return lookup


def metric_binary(y_true: pd.Series, y_score: pd.Series) -> dict[str, float]:
    y = pd.to_numeric(y_true, errors="coerce")
    s = pd.to_numeric(y_score, errors="coerce")
    m = y.notna() & s.notna()
    y = y[m].astype(int)
    s = s[m]
    out: dict[str, float] = {"n": float(len(y)), "positives": float(y.sum())}
    if len(y) < 2 or y.nunique() < 2:
        out.update({"AUROC": np.nan, "AUPRC": np.nan, "calibration_brier": np.nan})
        return out
    try:
        from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score

        out["AUROC"] = float(roc_auc_score(y, s))
        out["AUPRC"] = float(average_precision_score(y, s))
        out["calibration_brier"] = float(brier_score_loss(y, s.clip(0, 1)))
    except Exception:
        out.update({"AUROC": np.nan, "AUPRC": np.nan, "calibration_brier": np.nan})
    return out


def precision_recall_at_k(df: pd.DataFrame, score_col: str, label_col: str, k: int) -> dict[str, float]:
    if df.empty or score_col not in df or label_col not in df:
        return {f"Precision@{k}": np.nan, f"Recall@{k}": np.nan}
    d = df.copy()
    d[score_col] = pd.to_numeric(d[score_col], errors="coerce")
    d = d[d[score_col].notna()]
    if d.empty:
        return {f"Precision@{k}": np.nan, f"Recall@{k}": np.nan}
    d[label_col] = pd.to_numeric(d[label_col], errors="coerce").fillna(0).astype(int)
    total_pos = int(d[label_col].sum())
    d = d.sort_values(score_col, ascending=False).head(k)
    hits = int(d[label_col].sum())
    return {
        f"Precision@{k}": float(hits / max(len(d), 1)),
        f"Recall@{k}": float(hits / total_pos) if total_pos else np.nan,
    }


def patient_topn_metrics(df: pd.DataFrame, score_col: str, label_col: str, n: int) -> dict[str, float]:
    if df.empty or "patient_id" not in df:
        return {f"patient_hit_rate@{n}": np.nan, f"patient_recall@{n}": np.nan}
    d = df.copy()
    d[score_col] = pd.to_numeric(d[score_col], errors="coerce")
    d = d[d[score_col].notna()]
    if d.empty:
        return {f"patient_hit_rate@{n}": np.nan, f"patient_recall@{n}": np.nan, "patients_evaluated": 0.0}
    d[label_col] = pd.to_numeric(d[label_col], errors="coerce").fillna(0).astype(int)
    d["patient_id"] = d["patient_id"].fillna("unknown_patient")
    hits = []
    recalls = []
    for _, g in d.groupby("patient_id", dropna=False):
        g = g.sort_values(score_col, ascending=False)
        top = g.head(n)
        total_pos = int(g[label_col].sum())
        hit = int(top[label_col].sum() > 0)
        hits.append(hit)
        if total_pos:
            recalls.append(float(top[label_col].sum() / total_pos))
    return {
        f"patient_hit_rate@{n}": float(np.mean(hits)) if hits else np.nan,
        f"patient_recall@{n}": float(np.mean(recalls)) if recalls else np.nan,
        "patients_evaluated": float(len(hits)),
    }


def choose_label(df: pd.DataFrame) -> str | None:
    for c in ["label_immunogenicity", "label", "label_binary", "y", "immunogenic_positive"]:
        if c in df.columns and pd.to_numeric(df[c], errors="coerce").notna().sum() > 0:
            return c
    return None


def peptide_features(df: pd.DataFrame) -> pd.DataFrame:
    pep = df.get("peptide_mut", pd.Series("", index=df.index)).fillna("").astype(str)
    aa_sets = {
        "hydrophobic_frac": set("AILMFWYV"),
        "charged_frac": set("DEKRH"),
        "aromatic_frac": set("FWY"),
        "polar_frac": set("STNQCY"),
    }
    out = pd.DataFrame(index=df.index)
    out["peptide_length_num"] = pep.str.len().replace(0, np.nan)
    for name, aset in aa_sets.items():
        out[name] = pep.map(lambda x, aset=aset: sum(ch in aset for ch in x) / len(x) if x else np.nan)
    wt = df.get("peptide_wt", pd.Series("", index=df.index)).fillna("").astype(str)
    out["mut_wt_hamming"] = [
        sum(a != b for a, b in zip(m, w)) + abs(len(m) - len(w)) if m and w else np.nan
        for m, w in zip(pep, wt)
    ]
    for c in ["expression_tpm", "vaf", "apm_score", "tap1_expr", "tap2_expr"]:
        if c in df.columns:
            out[c] = pd.to_numeric(df[c], errors="coerce")
    return out


def simple_train_predict(
    df: pd.DataFrame,
    feature_cols: list[str],
    label_col: str,
    split_col: str | None,
) -> pd.Series:
    """Small, dependency-light estimator used only as a stack/ranker wrapper."""
    X = df[feature_cols].apply(pd.to_numeric, errors="coerce")
    X = X.fillna(X.median(numeric_only=True)).fillna(0.0)
    y = pd.to_numeric(df[label_col], errors="coerce").fillna(0).astype(int)
    if len(feature_cols) == 0 or y.nunique() < 2:
        return pd.Series(0.5, index=df.index)
    try:
        from sklearn.ensemble import GradientBoostingClassifier
        from sklearn.linear_model import LogisticRegression
        from sklearn.model_selection import GroupKFold, StratifiedKFold
        from sklearn.pipeline import make_pipeline
        from sklearn.preprocessing import StandardScaler

        pred = pd.Series(np.nan, index=df.index, dtype=float)
        groups = None
        if split_col and split_col in df.columns:
            groups = df[split_col].fillna("missing").astype(str)
        if groups is not None and groups.nunique() >= 3:
            splitter = GroupKFold(n_splits=min(5, groups.nunique()))
            splits = splitter.split(X, y, groups)
        else:
            splitter = StratifiedKFold(n_splits=min(5, int(y.value_counts().min())), shuffle=True, random_state=17)
            splits = splitter.split(X, y)
        for tr, te in splits:
            if y.iloc[tr].nunique() < 2:
                pred.iloc[te] = float(y.iloc[tr].mean())
                continue
            clf = make_pipeline(
                StandardScaler(),
                LogisticRegression(max_iter=1000, class_weight="balanced", random_state=17),
            )
            clf.fit(X.iloc[tr], y.iloc[tr])
            pred.iloc[te] = clf.predict_proba(X.iloc[te])[:, 1]
        if pred.isna().any():
            clf = GradientBoostingClassifier(random_state=17)
            clf.fit(X, y)
            pred[pred.isna()] = clf.predict_proba(X.loc[pred.isna()])[:, 1]
        return pred.clip(0, 1)
    except Exception:
        z = X.apply(norm_series).mean(axis=1)
        return norm_series(z)


def json_load(path: Path | str) -> object | None:
    try:
        return json.loads(Path(path).read_text(encoding="utf-8"))
    except Exception:
        return None
