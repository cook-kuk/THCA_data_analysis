#!/usr/bin/env python3
"""Stress tests for the top quantum-kernel fallback candidate."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import GroupKFold, RepeatedStratifiedKFold, StratifiedGroupKFold, StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC

import run_fallback_multiprocess_sweep as sweep


ROOT = Path(__file__).resolve().parent


FEATURE_SETS = {
    "quantum_structure_tcr": [
        *sweep.QUANTUM_SCORES,
        *sweep.STRUCTURE_FEATS,
        "tcr_motif_score",
    ],
    "quantum_structure_no_tcr": [
        *sweep.QUANTUM_SCORES,
        *sweep.STRUCTURE_FEATS,
    ],
    "quantum_only": sweep.QUANTUM_SCORES,
    "structure_only": sweep.STRUCTURE_FEATS,
    "quantum_tcr_no_structure": [*sweep.QUANTUM_SCORES, "tcr_motif_score"],
    "best_single_qk": ["W7A_QK_only"],
}


def qkernel(a: np.ndarray, b: np.ndarray, gamma: float = 1.0) -> np.ndarray:
    diff = a[:, None, :] - b[None, :, :]
    log_k = np.log(np.clip(np.cos(gamma * diff) ** 2, 1e-8, 1.0)).mean(axis=2)
    return np.exp(log_k)


def score(y: np.ndarray, pred: np.ndarray) -> tuple[float, float]:
    return float(roc_auc_score(y, pred)), float(average_precision_score(y, pred))


def cv_predict(x: np.ndarray, y: np.ndarray, splitter, groups=None) -> tuple[np.ndarray, int]:
    pred = np.full(len(y), np.nan, dtype=float)
    counts = np.zeros(len(y), dtype=float)
    n_used = 0
    splits = splitter.split(x, y, groups) if groups is not None else splitter.split(x, y)
    for train_idx, test_idx in splits:
        if len(np.unique(y[train_idx])) < 2 or len(np.unique(y[test_idx])) < 2:
            continue
        scaler = StandardScaler()
        xtr = scaler.fit_transform(x[train_idx])
        xte = scaler.transform(x[test_idx])
        clf = SVC(kernel="precomputed", class_weight="balanced", C=1.0)
        clf.fit(qkernel(xtr, xtr), y[train_idx])
        pred[test_idx] = np.nan_to_num(pred[test_idx], nan=0.0) + clf.decision_function(qkernel(xte, xtr))
        counts[test_idx] += 1
        n_used += 1
    mask = counts > 0
    pred[mask] = pred[mask] / counts[mask]
    return pred, n_used


def eval_set(df: pd.DataFrame, name: str, features: list[str]) -> list[dict]:
    cols = [c for c in features if c in df.columns]
    x = df[cols].to_numpy(float)
    y = df["label"].to_numpy(int)
    rows = []

    splitters = [
        ("repeated_stratified_5x5", RepeatedStratifiedKFold(n_splits=5, n_repeats=5, random_state=23), None),
        ("single_stratified_5fold", StratifiedKFold(n_splits=5, shuffle=True, random_state=23), None),
        ("hla_stratified_group_5fold", StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=23), df["hla"].to_numpy()),
    ]
    length_groups = df["peptide"].astype(str).str.len().to_numpy()
    if len(np.unique(length_groups)) >= 3:
        splitters.append(("length_group_3fold", GroupKFold(n_splits=3), length_groups))

    for split_name, splitter, groups in splitters:
        pred, n_used = cv_predict(x, y, splitter, groups)
        mask = np.isfinite(pred)
        if mask.sum() == 0 or len(np.unique(y[mask])) < 2:
            auc = auprc = np.nan
        else:
            auc, auprc = score(y[mask], pred[mask])
        rows.append(
            {
                "feature_set": name,
                "split": split_name,
                "n": int(mask.sum()),
                "n_pos": int(y[mask].sum()),
                "n_splits_used": int(n_used),
                "features": ",".join(cols),
                "auroc": auc,
                "auprc": auprc,
            }
        )
    return rows


def main() -> None:
    df = sweep.load_data()
    rows = []
    for name, feats in FEATURE_SETS.items():
        rows.extend(eval_set(df, name, feats))
    out = pd.DataFrame(rows).sort_values(["split", "auroc"], ascending=[True, False])
    out.to_csv(ROOT / "top_quantum_stress_tests.tsv", sep="\t", index=False)

    display = out[out["split"].isin(["repeated_stratified_5x5", "hla_stratified_group_5fold"])]
    lines = [
        "# Top Quantum Candidate Stress Tests",
        "",
        "These are still internal strict-set tests. They are meant to identify failure modes before external benchmarking.",
        "",
        "| feature_set | split | n / pos | AUROC | AUPRC |",
        "|---|---|---:|---:|---:|",
    ]
    for _, r in display.iterrows():
        lines.append(
            f"| {r['feature_set']} | {r['split']} | {int(r['n'])} / {int(r['n_pos'])} | {r['auroc']:.3f} | {r['auprc']:.3f} |"
        )
    lines.extend(
        [
            "",
            "Paper boundary: promote only if the signal survives fixed-feature, fixed-gamma, HLA/study/time-held-out external splits.",
        ]
    )
    (ROOT / "TOP_QUANTUM_STRESS_TESTS.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
