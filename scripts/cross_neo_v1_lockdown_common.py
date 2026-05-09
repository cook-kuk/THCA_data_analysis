#!/usr/bin/env python3
"""Shared helpers for CROSS-Neo v1 lockdown analyses.

The lockdown layer is intentionally conservative. It reuses v0 locked OOF
predictions, records every split/fold decision, and keeps public predictor
scores out of the feature set.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score


SEED = 20260509
REPO = Path(__file__).resolve().parents[1]
V0 = REPO / "project/results/cross_neo_v0"
OUT = REPO / "project/results/cross_neo_v1_lockdown"
FIG = OUT / "figures"
INPUT = REPO / "project/results/p_neo_bayesian_2026_05_09"

PRIMARY_SPLITS = [
    "repeated_stratified_5x5_internal",
    "exact_peptide_hla_holdout",
    "near_peptide_cluster_holdout",
    "hla_stratified_group_5fold",
    "hla_supertype_heldout",
]
SOURCE_STUDIES = ["CEDAR", "NEPdb", "TESLA_mmc4", "TESLA_mmc7_validation"]

ANCHOR_MAP = {
    ("C_counterfactual", "rf_secondary"): "anchor_rf",
    ("C_counterfactual", "elastic_net_lr"): "anchor_lr",
}
QK_BRANCHES = {
    "qk_no_anchor_gamma1": "qk_no_anchor",
    "qk_quantum_only_gamma1": "qk_quantum_only",
}


def ensure_dirs() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)


def sigmoid(x: np.ndarray) -> np.ndarray:
    return 1.0 / (1.0 + np.exp(-np.clip(x, -40, 40)))


def safe_float(x: object, default: float = np.nan) -> float:
    try:
        return float(x)
    except Exception:
        return default


def ece_score(y: np.ndarray, score: np.ndarray, n_bins: int = 10) -> float:
    y = np.asarray(y, dtype=int)
    s = np.clip(np.asarray(score, dtype=float), 0, 1)
    if len(y) == 0:
        return np.nan
    bins = np.linspace(0, 1, n_bins + 1)
    ece = 0.0
    for lo, hi in zip(bins[:-1], bins[1:]):
        mask = (s >= lo) & (s < hi if hi < 1 else s <= hi)
        if mask.sum() == 0:
            continue
        ece += float(mask.mean() * abs(s[mask].mean() - y[mask].mean()))
    return ece


def metric_row(y: np.ndarray, score: np.ndarray, k_values: tuple[int, ...] = (5, 10, 20)) -> dict[str, float]:
    y = np.asarray(y, dtype=int)
    score = np.asarray(score, dtype=float)
    mask = np.isfinite(score)
    y = y[mask]
    score = score[mask]
    out: dict[str, float] = {
        "n": int(len(y)),
        "n_pos": int(y.sum()) if len(y) else 0,
        "prevalence": float(y.mean()) if len(y) else np.nan,
    }
    if len(y) == 0 or len(np.unique(y)) < 2:
        out.update({"AUPRC": np.nan, "AUROC": np.nan, "Brier": np.nan, "ECE": np.nan})
    else:
        clipped = np.clip(score, 0, 1)
        out.update(
            {
                "AUPRC": float(average_precision_score(y, score)),
                "AUROC": float(roc_auc_score(y, score)),
                "Brier": float(brier_score_loss(y, clipped)),
                "ECE": float(ece_score(y, clipped)),
            }
        )
    order = np.argsort(-score)
    for k in k_values:
        kk = min(k, len(y))
        if kk == 0:
            prec = rec = enrich = np.nan
        else:
            top = y[order[:kk]]
            prec = float(top.mean())
            rec = float(top.sum() / max(1, y.sum()))
            enrich = float(prec / out["prevalence"]) if out["prevalence"] and out["prevalence"] > 0 else np.nan
        out[f"top{k}_precision"] = prec
        out[f"recall_at_{k}"] = rec
        out[f"enrichment_at_{k}"] = enrich
    return out


def add_ranks(df: pd.DataFrame, group_cols: list[str], score_col: str = "score") -> pd.DataFrame:
    out = df.copy()
    out["rank"] = out.groupby(group_cols)[score_col].rank(ascending=False, method="first").astype(int)
    out["rank_pct"] = out.groupby(group_cols)[score_col].rank(ascending=False, pct=True, method="average")
    out["n_in_rank_group"] = out.groupby(group_cols)[score_col].transform("size").astype(int)
    return out


def summarize_predictions(pred: pd.DataFrame, method_col: str = "method") -> pd.DataFrame:
    rows = []
    for keys, g in pred.groupby(["split_name", method_col], dropna=False):
        split_name, method = keys
        mm = metric_row(g["label"].to_numpy(int), g["score"].to_numpy(float))
        rows.append({"split_name": split_name, method_col: method, **mm})
    if not rows:
        return pd.DataFrame()
    return pd.DataFrame(rows).sort_values(["split_name", "AUPRC", "top10_precision"], ascending=[True, False, False])


def load_master() -> pd.DataFrame:
    return pd.read_csv(V0 / "master_table.tsv", sep="\t")


def load_folds() -> pd.DataFrame:
    return pd.read_csv(V0 / "folds.tsv", sep="\t")


def load_retrieval() -> pd.DataFrame:
    path = V0 / "retrieval_features_by_fold.tsv"
    return pd.read_csv(path, sep="\t") if path.exists() else pd.DataFrame()


def load_structure() -> pd.DataFrame:
    path = V0 / "structure_geometry_features.tsv"
    return pd.read_csv(path, sep="\t") if path.exists() else pd.DataFrame()


def source_split_name(study: str) -> str:
    return f"source_heldout_{study}"


def load_locked_base_predictions() -> pd.DataFrame:
    """Return anchor/QK predictions in one long locked table."""
    rows = []
    oof = pd.read_csv(V0 / "oof_predictions.tsv", sep="\t")
    for (feature_group, model), method in ANCHOR_MAP.items():
        sub = oof[(oof["feature_group"].eq(feature_group)) & (oof["model"].eq(model))].copy()
        if sub.empty:
            continue
        sub = sub[["split_name", "fold_id", "sample_id", "label", "score"]]
        sub["method"] = method
        sub["method_family"] = "clean_anchor"
        rows.append(sub)

    qk = pd.read_csv(V0 / "qk_fallback_predictions.tsv", sep="\t")
    for branch, method in QK_BRANCHES.items():
        sub = qk[qk["branch"].eq(branch)].copy()
        if sub.empty:
            continue
        sub = sub[["split_name", "fold_id", "sample_id", "label", "score"]]
        sub["method"] = method
        sub["method_family"] = "qk_fallback"
        rows.append(sub)

    source_path = V0 / "source_heldout_predictions.tsv"
    if source_path.exists():
        src = pd.read_csv(source_path, sep="\t")
        src_map = {
            "sourceheld_counterfactual_rf": ("anchor_rf", "source_clean_anchor"),
            "sourceheld_qk_compact_gamma1": ("source_qk_compact_gamma1", "source_qk_fallback"),
        }
        for src_method, (method, family) in src_map.items():
            sub = src[src["method"].eq(src_method)].copy()
            if sub.empty:
                continue
            sub["split_name"] = sub["heldout_study"].map(source_split_name)
            sub["fold_id"] = sub["heldout_study"]
            sub = sub[["split_name", "fold_id", "sample_id", "label", "score"]]
            sub["method"] = method
            sub["method_family"] = family
            rows.append(sub)

    if not rows:
        return pd.DataFrame(columns=["split_name", "fold_id", "sample_id", "label", "score", "method", "method_family"])
    pred = pd.concat(rows, ignore_index=True)
    return add_ranks(pred, ["split_name", "fold_id", "method"])


def pivot_methods(pred: pd.DataFrame, methods: list[str]) -> pd.DataFrame:
    sub = pred[pred["method"].isin(methods)].copy()
    if sub.empty:
        return pd.DataFrame()
    wide = sub.pivot_table(
        index=["split_name", "fold_id", "sample_id", "label"],
        columns="method",
        values="score",
        aggfunc="mean",
    ).reset_index()
    wide.columns.name = None
    return wide


def strict_train_ids_for_outer_fold(split_name: str, fold_id: str, folds: pd.DataFrame, master: pd.DataFrame) -> set[str]:
    strict_ids = set(master.loc[master["strict_set_flag"].astype(bool), "sample_id"])
    test_ids = set(folds.loc[(folds["split_name"].eq(split_name)) & (folds["fold_id"].eq(fold_id)), "sample_id"])
    return strict_ids - test_ids


def select_weight_from_inner(
    wide: pd.DataFrame,
    split_name: str,
    fold_id: str,
    test_ids: set[str],
    anchor_col: str,
    qk_col: str,
    candidate_anchor_weights: list[float],
) -> tuple[float, float, int, str]:
    inner = wide[
        (wide["split_name"].eq(split_name))
        & (~wide["sample_id"].isin(test_ids))
        & wide[anchor_col].notna()
        & wide[qk_col].notna()
    ].copy()
    if inner["label"].nunique() < 2 or len(inner) < 10:
        return 0.5, np.nan, int(len(inner)), "default_no_inner_signal"
    scored = []
    for aw in candidate_anchor_weights:
        pred = aw * inner[anchor_col].to_numpy(float) + (1.0 - aw) * inner[qk_col].to_numpy(float)
        mm = metric_row(inner["label"].to_numpy(int), pred)
        scored.append((aw, safe_float(mm["AUPRC"], -1), safe_float(mm["top10_precision"], -1)))
    scored.sort(key=lambda x: (x[1], x[2]), reverse=True)
    return float(scored[0][0]), float(scored[0][1]), int(len(inner)), "inner_oof_selected"


def bootstrap_ci(y: np.ndarray, score: np.ndarray, metric: str, n_boot: int = 300) -> tuple[float, float]:
    y = np.asarray(y, dtype=int)
    score = np.asarray(score, dtype=float)
    if len(y) < 5:
        return np.nan, np.nan
    rng = np.random.default_rng(SEED)
    vals: list[float] = []
    for _ in range(n_boot):
        idx = rng.integers(0, len(y), len(y))
        yy = y[idx]
        ss = score[idx]
        if len(np.unique(yy)) < 2:
            continue
        if metric == "AUPRC":
            vals.append(float(average_precision_score(yy, ss)))
        elif metric == "top10_precision":
            vals.append(metric_row(yy, ss)["top10_precision"])
    if not vals:
        return np.nan, np.nan
    return float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


def paired_bootstrap_delta(
    df: pd.DataFrame,
    score_col: str,
    ref_score_col: str,
    metric: str,
    n_boot: int = 300,
) -> tuple[float, float, float]:
    d = df[["label", score_col, ref_score_col]].dropna()
    if len(d) < 5 or d["label"].nunique() < 2:
        return np.nan, np.nan, np.nan
    y = d["label"].to_numpy(int)
    s = d[score_col].to_numpy(float)
    r = d[ref_score_col].to_numpy(float)
    def val(yy: np.ndarray, ss: np.ndarray) -> float:
        if metric == "AUPRC":
            return float(average_precision_score(yy, ss))
        if metric == "top10_precision":
            return metric_row(yy, ss)["top10_precision"]
        raise ValueError(metric)
    delta = val(y, s) - val(y, r)
    rng = np.random.default_rng(SEED + 17)
    vals = []
    for _ in range(n_boot):
        idx = rng.integers(0, len(y), len(y))
        yy = y[idx]
        if len(np.unique(yy)) < 2:
            continue
        vals.append(val(yy, s[idx]) - val(yy, r[idx]))
    if not vals:
        return delta, np.nan, np.nan
    return float(delta), float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


def peptide_cluster(seq: str) -> str:
    s = str(seq or "")
    if len(s) < 4:
        return f"L{len(s)}_{s}"
    return f"L{len(s)}_{s[:2]}_{s[-2:]}"


def hla_supertype_from_hla(hla: str) -> str:
    h = str(hla or "")
    if "*" not in h:
        return ""
    locus = h.split("*", 1)[0].replace("HLA-", "")
    fam = h.split("*", 1)[1].split(":", 1)[0]
    return f"{locus}{fam}"
