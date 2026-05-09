#!/usr/bin/env python3
"""Fixed-gamma fold-safe quantum-kernel fallback predictions."""

from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC

from cross_neo_v0_common import INPUT, OUT, ensure_dirs, load_master, load_or_build_folds


Q_COLS = ["GP_quantum", "VQC", "W7A_QK_only", "W7A_full"]
GEOM_COLS = [
    "mean_pLDDT_peptide",
    "min_pLDDT_peptide",
    "interface_contacts_8A",
    "interface_contacts_10A",
    "radius_of_gyration_peptide",
    "peptide_helicity_proxy",
]
GAMMA = 1.0


def qkernel(a: np.ndarray, b: np.ndarray, gamma: float = GAMMA) -> np.ndarray:
    diff = a[:, None, :] - b[None, :, :]
    log_k = np.log(np.clip(np.cos(gamma * diff) ** 2, 1e-8, 1.0)).mean(axis=2)
    return np.exp(log_k)


def fit_predict(train_x: np.ndarray, train_y: np.ndarray, test_x: np.ndarray) -> np.ndarray:
    scaler = StandardScaler()
    xtr = scaler.fit_transform(train_x)
    xte = scaler.transform(test_x)
    clf = SVC(kernel="precomputed", C=1.0, class_weight="balanced")
    clf.fit(qkernel(xtr, xtr), train_y)
    raw = clf.decision_function(qkernel(xte, xtr))
    return 1.0 / (1.0 + np.exp(-raw))


def main() -> None:
    ensure_dirs()
    master = load_master()
    folds = load_or_build_folds(master)
    wide = pd.read_csv(INPUT / "curation_2026_05_09/curated_predictions_itsndb_wide.tsv", sep="\t")
    q = master[["sample_id", "peptide_mut", "hla", "label", "strict_set_flag", "tcr_motif_score"]].copy()
    q = q.merge(
        wide[["peptide", "hla", "label", *Q_COLS]],
        left_on=["peptide_mut", "hla", "label"],
        right_on=["peptide", "hla", "label"],
        how="left",
    )
    geom = pd.read_csv(OUT / "structure_geometry_features.tsv", sep="\t")
    q = q.merge(geom[["sample_id", *GEOM_COLS]], on="sample_id", how="left")
    for c in Q_COLS + GEOM_COLS + ["tcr_motif_score"]:
        q[c] = q[c].fillna(q[c].median() if q[c].notna().any() else 0.0)

    strict_ids = set(master.loc[master["strict_set_flag"].astype(bool), "sample_id"])
    rows = []
    for (split_name, fold_id), fold in folds.groupby(["split_name", "fold_id"], sort=False):
        test_ids = set(fold["sample_id"])
        train_ids = list(strict_ids - test_ids)
        test_ids_list = list(test_ids)
        train = q[q["sample_id"].isin(train_ids)].copy()
        test = q[q["sample_id"].isin(test_ids_list)].copy()
        if len(train) == 0 or len(test) == 0 or train["label"].nunique() < 2:
            continue
        feature_sets = {
            "qk_no_anchor_gamma1": Q_COLS + GEOM_COLS + ["tcr_motif_score"],
            "qk_clean_no_tcr_gamma1": Q_COLS + GEOM_COLS,
            "qk_quantum_only_gamma1": Q_COLS,
        }
        for name, cols in feature_sets.items():
            pred = fit_predict(train[cols].to_numpy(float), train["label"].to_numpy(int), test[cols].to_numpy(float))
            for sid, y, p in zip(test["sample_id"], test["label"], pred):
                rows.append(
                    {
                        "split_name": split_name,
                        "fold_id": fold_id,
                        "sample_id": sid,
                        "branch": name,
                        "label": int(y),
                        "score": float(p),
                        "features": ",".join(cols),
                        "gamma": GAMMA,
                    }
                )
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "qk_fallback_predictions.tsv", sep="\t", index=False)
    print(f"[qk] wrote rows={len(out)}")


if __name__ == "__main__":
    main()
