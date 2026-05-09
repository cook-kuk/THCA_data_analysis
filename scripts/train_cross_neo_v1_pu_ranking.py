#!/usr/bin/env python3
"""PU-aware ranking models for CROSS-Neo v1."""

from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

from cross_neo_v1_common import V1, build_base_feature_table, ensure_v1_dirs, fill_by_train_median, get_fold_ids, load_folds, load_master, metrics


def feature_cols(base: pd.DataFrame) -> list[str]:
    cols = [c for c in base.columns if c.startswith("cf_")]
    cols += [c for c in ["GP_quantum", "VQC", "W7A_QK_only", "W7A_full", "mean_pLDDT_peptide", "structure_missing", "structure_low_confidence", "peptide_len"] if c in base.columns]
    return cols


def pu_logistic(train: pd.DataFrame, test: pd.DataFrame, cols: list[str]) -> np.ndarray:
    xtr, xte = fill_by_train_median(train, test, cols)
    y = train["label"].to_numpy(int)
    w = np.where(y == 1, 1.5, 0.30)
    clf = LogisticRegression(max_iter=3000, class_weight=None, C=0.25)
    scaler = StandardScaler()
    xtr = scaler.fit_transform(xtr)
    xte = scaler.transform(xte)
    clf.fit(xtr, y, sample_weight=w)
    return clf.predict_proba(xte)[:, 1]


def bagging_pu_rf(train: pd.DataFrame, test: pd.DataFrame, cols: list[str], n_bags: int = 25) -> np.ndarray:
    xtr, xte = fill_by_train_median(train, test, cols)
    y = train["label"].to_numpy(int)
    pos = np.where(y == 1)[0]
    unl = np.where(y == 0)[0]
    if len(pos) == 0 or len(unl) == 0:
        return np.full(len(test), y.mean() if len(y) else 0.5)
    rng = np.random.default_rng(20260509)
    preds = []
    for b in range(n_bags):
        neg = rng.choice(unl, size=min(len(unl), max(len(pos) * 2, 8)), replace=len(unl) < max(len(pos) * 2, 8))
        keep = np.concatenate([pos, neg])
        clf = RandomForestClassifier(n_estimators=80, max_depth=3, min_samples_leaf=4, class_weight="balanced", random_state=20260509 + b, n_jobs=-1)
        clf.fit(xtr[keep], y[keep])
        preds.append(clf.predict_proba(xte)[:, 1])
    return np.mean(preds, axis=0)


def positive_centroid(train: pd.DataFrame, test: pd.DataFrame, cols: list[str]) -> np.ndarray:
    xtr, xte = fill_by_train_median(train, test, cols)
    scaler = StandardScaler()
    xtr = scaler.fit_transform(xtr)
    xte = scaler.transform(xte)
    pos = xtr[train["label"].to_numpy(int) == 1]
    if len(pos) == 0:
        return np.full(len(test), 0.5)
    centroid = pos.mean(axis=0)
    d = np.linalg.norm(xte - centroid, axis=1)
    return 1.0 / (1.0 + d)


def pairwise_ranker(train: pd.DataFrame, test: pd.DataFrame, cols: list[str]) -> np.ndarray:
    xtr, xte = fill_by_train_median(train, test, cols)
    scaler = StandardScaler()
    xtr = scaler.fit_transform(xtr)
    xte = scaler.transform(xte)
    y = train["label"].to_numpy(int)
    rng = np.random.default_rng(20260509)
    diffs, labels = [], []
    pos_idx = np.where(y == 1)[0]
    neg_idx = np.where(y == 0)[0]
    if len(pos_idx) < 2 or len(neg_idx) < 2:
        return pu_logistic(train, test, cols)
    for p in pos_idx:
        choices = rng.choice(neg_idx, size=min(8, len(neg_idx)), replace=False)
        for n in choices:
            diffs.append(xtr[p] - xtr[n]); labels.append(1)
            diffs.append(xtr[n] - xtr[p]); labels.append(0)
    clf = LogisticRegression(max_iter=2000, C=0.5)
    clf.fit(np.asarray(diffs), np.asarray(labels))
    return 1.0 / (1.0 + np.exp(-np.clip(xte @ clf.coef_.ravel(), -30, 30)))


def eval_locked(base: pd.DataFrame, master: pd.DataFrame, folds: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    rows = []
    models = {
        "pu_logistic_weighted": pu_logistic,
        "bagging_pu_rf": bagging_pu_rf,
        "positive_only_centroid": positive_centroid,
        "pairwise_positive_over_unlabeled": pairwise_ranker,
    }
    for (split, fold), _ in folds.groupby(["split_name", "fold_id"], sort=False):
        train_ids, test_ids = get_fold_ids(master, folds, split, fold)
        train = base[base["sample_id"].isin(train_ids)].copy()
        test = base[base["sample_id"].isin(test_ids)].copy()
        if len(train) < 12 or len(test) < 2 or train["label"].nunique() < 2 or test["label"].nunique() < 2:
            continue
        for name, func in models.items():
            score = func(train, test, cols)
            for sid, y, p in zip(test["sample_id"], test["label"], score):
                rows.append({"split_name": split, "fold_id": fold, "model": name, "sample_id": sid, "label": int(y), "score": float(p)})
    return pd.DataFrame(rows)


def eval_tesla_source(base: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    rows = []
    train_pool = base[base["study"].isin(["CEDAR", "NEPdb", "TESLA_mmc4", "TESLA_mmc7_validation"])].copy()
    models = {"pu_logistic_weighted": pu_logistic, "bagging_pu_rf": bagging_pu_rf, "positive_only_centroid": positive_centroid, "pairwise_positive_over_unlabeled": pairwise_ranker}
    for heldout in ["TESLA_mmc4", "TESLA_mmc7_validation", "NEPdb"]:
        test = train_pool[train_pool["study"] == heldout].copy()
        train = train_pool[train_pool["study"] != heldout].copy()
        if len(test) < 10 or test["label"].nunique() < 2 or train["label"].nunique() < 2:
            continue
        for name, func in models.items():
            score = func(train, test, cols)
            for sid, y, p in zip(test["sample_id"], test["label"], score):
                rows.append({"split_name": f"source_heldout_{heldout}", "fold_id": heldout, "model": name, "sample_id": sid, "label": int(y), "score": float(p)})
    return pd.DataFrame(rows)


def main() -> None:
    ensure_v1_dirs()
    master = load_master()
    folds = load_folds()
    base = build_base_feature_table(master)
    cols = feature_cols(base)
    pred = pd.concat([eval_locked(base, master, folds, cols), eval_tesla_source(base, cols)], ignore_index=True)
    pred.to_csv(V1 / "pu_ranking_predictions.tsv", sep="\t", index=False)
    rows = []
    for (split, model), sub in pred.groupby(["split_name", "model"]):
        rows.append({"split_name": split, "model": model, **metrics(sub["label"].to_numpy(), sub["score"].to_numpy())})
    met = pd.DataFrame(rows).sort_values(["split_name", "AUPRC"], ascending=[True, False])
    met.to_csv(V1 / "pu_ranking_metrics.tsv", sep="\t", index=False)
    lines = [
        "# CROSS-Neo v1 PU Ranking Report",
        "",
        "Negatives are treated as unlabeled/ambiguous for weighted and bagged PU variants. No public predictor scores are used.",
        "",
        "## Top PU Rows",
        met.sort_values("AUPRC", ascending=False).head(18).to_markdown(index=False),
        "",
        "PU models are promotion candidates only if they improve top-k precision, not merely AUROC.",
    ]
    (V1 / "pu_ranking_report.md").write_text("\n".join(lines) + "\n")
    print(f"[v1-pu] predictions={len(pred)} metrics={len(met)}")


if __name__ == "__main__":
    main()
