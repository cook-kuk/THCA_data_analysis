#!/usr/bin/env python3
"""Train shallow split-safe CROSS-Neo v0 models."""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
from sklearn.calibration import CalibratedClassifierCV
from sklearn.ensemble import RandomForestClassifier
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import LinearSVC

from cross_neo_v0_common import INPUT, OUT, ece_score, ensure_dirs, load_master, load_or_build_folds, metrics

warnings.filterwarnings("ignore", category=ConvergenceWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

try:
    from xgboost import XGBClassifier
except Exception:  # pragma: no cover - optional dependency
    XGBClassifier = None


FORBIDDEN = {"MHCflurry", "NetMHCpan_4.1", "NetMHCstabpan", "BigMHC_IM", "PRIME", "MixMHCpred"}


def load_feature_matrix(master: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, list[str]]]:
    base = master[["sample_id", "label", "strict_set_flag"]].copy()

    struct = pd.read_csv(INPUT / "wave2/structure_features.tsv", sep="\t")
    struct_cols = [c for c in struct.columns if c not in ("peptide", "HLA_norm", "label", "split", "in_master")]
    base = base.merge(
        master[["sample_id", "peptide_mut", "hla"]].merge(
            struct[["peptide", "HLA_norm", *struct_cols]],
            left_on=["peptide_mut", "hla"],
            right_on=["peptide", "HLA_norm"],
            how="left",
        ).drop(columns=["peptide_mut", "hla", "peptide", "HLA_norm"]),
        on="sample_id",
        how="left",
    )

    geom = pd.read_csv(OUT / "structure_geometry_features.tsv", sep="\t")
    geom_cols = [
        c
        for c in geom.columns
        if c
        not in {
            "sample_id",
            "peptide_mut",
            "hla",
            "label",
            "strict_set_flag",
            "mutant_wt_surface_delta",
            "ensemble_variance",
        }
    ]
    base = base.merge(geom[["sample_id", *geom_cols]], on="sample_id", how="left")

    cf = np.load(OUT / "counterfactual_embeddings.npy")
    names = pd.read_csv(OUT / "counterfactual_feature_names.tsv", sep="\t")["feature"].tolist()
    idx = pd.read_csv(OUT / "counterfactual_feature_index.tsv", sep="\t")
    cf_df = pd.DataFrame(cf, columns=[f"cf__{n}" for n in names])
    cf_df.insert(0, "sample_id", idx["sample_id"].values)
    # Keep dimensions manageable for n=89: use first 80 deterministic features.
    cf_cols = [c for c in cf_df.columns if c != "sample_id"][:80]
    base = base.merge(cf_df[["sample_id", *cf_cols]], on="sample_id", how="left")

    wide = pd.read_csv(INPUT / "curation_2026_05_09/curated_predictions_itsndb_wide.tsv", sep="\t")
    quantum_cols = ["GP_quantum", "VQC", "W7A_QK_only", "W7A_full"]
    q_df = master[["sample_id", "peptide_mut", "hla", "label"]].merge(
        wide[["peptide", "hla", "label", *quantum_cols]],
        left_on=["peptide_mut", "hla", "label"],
        right_on=["peptide", "hla", "label"],
        how="left",
    )[["sample_id", *quantum_cols]]
    base = base.merge(q_df, on="sample_id", how="left")

    # Fold-specific retrieval features are merged inside each fold.
    groups = {
        "A_structure_baseline_features": struct_cols,
        "C_counterfactual": cf_cols,
        "D_structure_geometry": geom_cols,
        "E_quantum_fixed": quantum_cols,
    }
    for c in base.columns:
        if c not in ("sample_id", "label", "strict_set_flag"):
            base[c] = base[c].astype(float)
            base[c] = base[c].fillna(base[c].median() if base[c].notna().any() else 0.0)
    return base, groups


def make_model(model_name: str):
    if model_name == "elastic_net_lr":
        return make_pipeline(
            StandardScaler(),
            LogisticRegression(
                penalty="elasticnet",
                solver="saga",
                l1_ratio=0.25,
                C=0.5,
                class_weight="balanced",
                max_iter=5000,
                random_state=20260509,
            ),
        )
    if model_name == "calibrated_linear_svm":
        return make_pipeline(
            StandardScaler(),
            CalibratedClassifierCV(
                LinearSVC(C=0.5, class_weight="balanced", random_state=20260509, max_iter=5000),
                cv=3,
                method="sigmoid",
            ),
        )
    if model_name == "rf_secondary":
        return RandomForestClassifier(
            n_estimators=120,
            max_depth=3,
            min_samples_leaf=5,
            class_weight="balanced",
            random_state=20260509,
            n_jobs=-1,
        )
    if model_name == "xgboost_depth2":
        if XGBClassifier is None:
            raise ValueError("xgboost is not installed")
        return XGBClassifier(
            n_estimators=60,
            max_depth=2,
            learning_rate=0.05,
            subsample=0.85,
            colsample_bytree=0.85,
            objective="binary:logistic",
            eval_metric="logloss",
            random_state=20260509,
            n_jobs=1,
            reg_lambda=2.0,
            min_child_weight=2,
        )
    raise ValueError(model_name)


def predict_proba(model, x):
    if hasattr(model, "predict_proba"):
        return model.predict_proba(x)[:, 1]
    return model.decision_function(x)


def available_model_names() -> list[str]:
    names = ["elastic_net_lr", "calibrated_linear_svm"]
    if XGBClassifier is not None:
        names.append("xgboost_depth2")
    names.append("rf_secondary")
    return names


def main() -> None:
    ensure_dirs()
    master = load_master()
    folds = load_or_build_folds(master)
    feats, groups = load_feature_matrix(master)
    retr = pd.read_csv(OUT / "retrieval_features_by_fold.tsv", sep="\t")
    retr_cols = [
        "exact_peptide_hla_hit_train",
        "exact_peptide_hit_train",
        "near_peptide_similarity_train",
        "near_positive_peptide_similarity_train",
        "same_hla_near_peptide_similarity_train",
        "wt_self_similarity",
        "same_hla_train_density",
        "same_hla_positive_rate_train",
        "train_fold_prevalence",
    ]
    groups["B_retrieval_only_clean_flags"] = retr_cols
    groups["F_CROSS_all"] = (
        groups["A_structure_baseline_features"]
        + retr_cols
        + groups["C_counterfactual"]
        + groups["D_structure_geometry"]
        + groups["E_quantum_fixed"]
    )

    strict_ids = set(master.loc[master["strict_set_flag"].astype(bool), "sample_id"])
    pred_rows = []
    metric_rows = []
    model_names = available_model_names()
    for (split_name, fold_id), fold in folds.groupby(["split_name", "fold_id"], sort=False):
        test_ids = set(fold["sample_id"])
        train_ids = strict_ids - test_ids
        rfold = retr[(retr["split_name"] == split_name) & (retr["fold_id"] == fold_id)]
        fold_feat = feats.merge(rfold[["sample_id", *retr_cols]], on="sample_id", how="left")
        for c in retr_cols:
            fold_feat[c] = fold_feat[c].fillna(0.0).astype(float)
        train = fold_feat[fold_feat["sample_id"].isin(train_ids)].copy()
        test = fold_feat[fold_feat["sample_id"].isin(test_ids)].copy()
        if len(train) < 10 or len(test) < 2 or train["label"].nunique() < 2 or test["label"].nunique() < 2:
            continue
        for group_name, cols in groups.items():
            cols = [c for c in cols if c in fold_feat.columns and c not in FORBIDDEN]
            if not cols:
                continue
            xtr = train[cols].to_numpy(float)
            ytr = train["label"].to_numpy(int)
            xte = test[cols].to_numpy(float)
            yte = test["label"].to_numpy(int)
            for model_name in model_names:
                model = make_model(model_name)
                model.fit(xtr, ytr)
                pred = np.clip(predict_proba(model, xte), 0, 1)
                for sid, y, p in zip(test["sample_id"], yte, pred):
                    pred_rows.append(
                        {
                            "split_name": split_name,
                            "fold_id": fold_id,
                            "sample_id": sid,
                            "feature_group": group_name,
                            "model": model_name,
                            "label": int(y),
                            "score": float(p),
                            "n_features": len(cols),
                        }
                    )
    pred_df = pd.DataFrame(pred_rows)
    pred_df.to_csv(OUT / "oof_predictions.tsv", sep="\t", index=False)
    for keys, sub in pred_df.groupby(["split_name", "feature_group", "model"], sort=False):
        m = metrics(sub["label"].to_numpy(), sub["score"].to_numpy())
        m["ECE"] = ece_score(sub["label"].to_numpy(), sub["score"].to_numpy())
        metric_rows.append({"split_name": keys[0], "feature_group": keys[1], "model": keys[2], **m})
    metrics_df = pd.DataFrame(metric_rows).sort_values(["split_name", "AUPRC"], ascending=[True, False])
    metrics_df.to_csv(OUT / "metrics_by_split.tsv", sep="\t", index=False)

    ablation = (
        metrics_df[metrics_df["split_name"] == "repeated_stratified_5x5_internal"]
        .sort_values("AUPRC", ascending=False)
        .copy()
    )
    ablation.to_csv(OUT / "ablation_table.tsv", sep="\t", index=False)
    print(f"[train] predictions={len(pred_df)} metrics={len(metrics_df)}")


if __name__ == "__main__":
    main()
