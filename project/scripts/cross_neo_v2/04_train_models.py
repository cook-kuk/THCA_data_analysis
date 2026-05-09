#!/usr/bin/env python3
"""Train CROSS-Neo 2.0 first-pass locked models and import v1 baselines."""

from __future__ import annotations

import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import GradientBoostingClassifier, HistGradientBoostingClassifier, RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from common import OUT, SEED, V1, ensure_dirs, safe_parquet

warnings.filterwarnings("ignore")


PRIMARY_SPLITS = {
    "exact_peptide_hla_holdout",
    "near_peptide_cluster_holdout",
    "hla_stratified_group_5fold",
    "hla_supertype_heldout",
}


def read_feature(path: Path) -> pd.DataFrame:
    if path.exists():
        try:
            return pd.read_parquet(path)
        except Exception:
            pass
    pkl = path.with_suffix(path.suffix + ".pkl")
    if pkl.exists():
        return pd.read_pickle(pkl)
    tsv = path.with_suffix(".tsv")
    if tsv.exists():
        return pd.read_csv(tsv, sep="\t")
    raise FileNotFoundError(path)


def numeric_features(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    for c in out.columns:
        if c != "row_id":
            out[c] = pd.to_numeric(out[c], errors="coerce")
    return out.fillna(0)


def merge_features(frames: list[pd.DataFrame]) -> pd.DataFrame:
    x = frames[0].copy()
    for f in frames[1:]:
        x = x.merge(f, on="row_id", how="left", suffixes=("", "_dup"))
    dup = [c for c in x.columns if c.endswith("_dup")]
    if dup:
        x = x.drop(columns=dup)
    return numeric_features(x).fillna(0)


def make_estimator(kind: str):
    if kind == "lr":
        return Pipeline([
            ("imputer", SimpleImputer(strategy="median")),
            ("scaler", StandardScaler()),
            ("clf", LogisticRegression(C=0.5, class_weight="balanced", max_iter=3000, solver="liblinear", random_state=SEED)),
        ])
    if kind == "rf":
        return Pipeline([
            ("imputer", SimpleImputer(strategy="median")),
            ("clf", RandomForestClassifier(
                n_estimators=500,
                max_depth=3,
                min_samples_leaf=3,
                class_weight="balanced_subsample",
                random_state=SEED,
                n_jobs=-1,
            )),
        ])
    if kind == "hgb":
        return Pipeline([
            ("imputer", SimpleImputer(strategy="median")),
            ("clf", HistGradientBoostingClassifier(
                learning_rate=0.04,
                max_iter=90,
                max_leaf_nodes=7,
                l2_regularization=1.0,
                random_state=SEED,
            )),
        ])
    return Pipeline([
        ("imputer", SimpleImputer(strategy="median")),
        ("clf", GradientBoostingClassifier(n_estimators=80, learning_rate=0.04, max_depth=2, random_state=SEED)),
    ])


def predict_score(model, x: pd.DataFrame) -> np.ndarray:
    if hasattr(model, "predict_proba"):
        return model.predict_proba(x)[:, 1]
    if hasattr(model, "decision_function"):
        raw = np.asarray(model.decision_function(x), dtype=float)
    else:
        raw = np.asarray(model.predict(x), dtype=float)
    lo, hi = np.nanmin(raw), np.nanmax(raw)
    return (raw - lo) / (hi - lo + 1e-9)


def fit_model(model, x: pd.DataFrame, y: np.ndarray, weights: np.ndarray | None):
    try:
        if weights is not None:
            return model.fit(x, y, clf__sample_weight=weights)
    except Exception:
        pass
    return model.fit(x, y)


def train_weights(train_meta: pd.DataFrame, mode: str) -> np.ndarray | None:
    if mode == "none":
        return None
    y = train_meta["label_binary"].astype(int)
    w = pd.Series(1.0, index=train_meta.index)
    if mode in {"source_balanced", "groupdro_proxy"}:
        source_counts = train_meta["source_dataset"].fillna("unknown").value_counts()
        w *= train_meta["source_dataset"].fillna("unknown").map(lambda s: 1.0 / max(1, source_counts.get(s, 1))).astype(float)
        w *= len(train_meta) / max(1e-9, w.sum())
    if mode in {"class_balanced", "groupdro_proxy"}:
        cls_counts = y.value_counts()
        w *= y.map(lambda yy: len(y) / (2.0 * max(1, cls_counts.get(yy, 1)))).astype(float)
    if mode == "groupdro_proxy":
        hla_counts = train_meta["hla_supertype"].fillna("unknown").value_counts()
        hla_w = train_meta["hla_supertype"].fillna("unknown").map(lambda s: 1.0 / np.sqrt(max(1, hla_counts.get(s, 1)))).astype(float)
        w *= hla_w
    w = np.asarray(w, dtype=float)
    return np.clip(w / max(1e-9, np.mean(w)), 0.2, 8.0)


def train_percentile(train_scores: np.ndarray, train_groups: pd.Series, test_scores: np.ndarray, test_groups: pd.Series) -> np.ndarray:
    global_scores = np.asarray(train_scores, dtype=float)
    out = []
    for score, group in zip(test_scores, test_groups):
        ref = global_scores
        m = train_groups.astype(str).eq(str(group)).values
        if m.sum() >= 8:
            ref = global_scores[m]
        out.append(float((ref <= score).mean()))
    return np.asarray(out, dtype=float)


def import_v1_predictions(reg_ids: set[str]) -> pd.DataFrame:
    rows = []
    anchor = V1 / "locked_anchor_predictions.tsv"
    if anchor.exists():
        a = pd.read_csv(anchor, sep="\t")
        a = a.rename(columns={"sample_id": "row_id", "method": "model_name", "method_family": "model_family"})
        a = a[a["row_id"].astype(str).isin(reg_ids)].copy()
        a["claim_status"] = np.where(a["model_name"].isin(["anchor_rf", "anchor_lr"]), "reviewer_safe_internal_locked", "diagnostic")
        a["source_script"] = "v1_lockdown_import"
        rows.append(a[["split_name", "fold_id", "row_id", "label", "score", "model_name", "model_family", "claim_status", "source_script"]])
    fusion = V1 / "foldsafe_fusion_predictions.tsv"
    if fusion.exists():
        f = pd.read_csv(fusion, sep="\t")
        f = f.rename(columns={"sample_id": "row_id", "method": "model_name"})
        f = f[f["row_id"].astype(str).isin(reg_ids)].copy()
        f["model_family"] = "v1_foldsafe_fusion"
        safe = f["model_name"].astype(str).str.contains("w0.5|nested|rule_gate", regex=True)
        f["claim_status"] = np.where(safe, "reviewer_safe_internal_locked", "exploratory_descriptive_only")
        f["source_script"] = "v1_lockdown_import"
        rows.append(f[["split_name", "fold_id", "row_id", "label", "score", "model_name", "model_family", "claim_status", "source_script"]])
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def add_ranks(pred: pd.DataFrame) -> pd.DataFrame:
    pred = pred.copy()
    pred["rank"] = pred.groupby(["split_name", "fold_id", "model_name"])["score"].rank(method="first", ascending=False)
    sizes = pred.groupby(["split_name", "fold_id", "model_name"])["row_id"].transform("size")
    pred["rank_pct"] = pred["rank"] / sizes.clip(lower=1)
    pred["n_in_rank_group"] = sizes
    return pred


def main() -> None:
    ensure_dirs()
    rng = np.random.default_rng(SEED)
    reg = pd.read_csv(OUT / "canonical_registry.tsv", sep="\t")
    reg["label_binary"] = pd.to_numeric(reg["label_binary"], errors="coerce").fillna(0).astype(int)
    reg_ids = set(reg["row_id"].astype(str))

    cf = read_feature(OUT / "features/counterfactual_features.parquet")
    plm = read_feature(OUT / "features/plm_embeddings.parquet")
    struct = read_feature(OUT / "features/structure_features.parquet")
    proc = read_feature(OUT / "features/processing_context_features.parquet")
    feature_sets = {
        "counterfactual": numeric_features(cf),
        "plm_frozen_hash": numeric_features(plm),
        "counterfactual_plm": merge_features([cf, plm]),
        "multimodal": merge_features([cf, plm, struct, proc]),
    }
    esm2_path = OUT / "features/esm2_35m_features.parquet"
    if esm2_path.exists() or esm2_path.with_suffix(esm2_path.suffix + ".pkl").exists() or (OUT / "features/esm2_35m_features.tsv").exists():
        esm2 = read_feature(esm2_path)
        feature_sets["esm2_35m"] = numeric_features(esm2)
        feature_sets["counterfactual_esm2_35m"] = merge_features([cf, esm2])
        feature_sets["multimodal_esm2_35m"] = merge_features([cf, esm2, struct, proc])

    specs = [
        ("v2_counterfactual_lr", "classical_counterfactual", "counterfactual", "lr", "none", "reviewer_safe_internal_locked"),
        ("v2_counterfactual_rf", "classical_counterfactual", "counterfactual", "rf", "none", "reviewer_safe_internal_locked"),
        ("v2_counterfactual_hgb", "classical_counterfactual", "counterfactual", "hgb", "none", "reviewer_safe_internal_locked"),
        ("v2_plm_lr_frozen_hash_pilot", "frozen_plm_pilot", "plm_frozen_hash", "lr", "none", "diagnostic"),
        ("v2_plm_rf_frozen_hash_pilot", "frozen_plm_pilot", "plm_frozen_hash", "rf", "none", "diagnostic"),
        ("v2_cf_plm_lr", "counterfactual_plm", "counterfactual_plm", "lr", "none", "reviewer_safe_internal_locked"),
        ("v2_cf_plm_rf", "counterfactual_plm", "counterfactual_plm", "rf", "none", "reviewer_safe_internal_locked"),
        ("v2_multimodal_lr", "multimodal_structure_processing", "multimodal", "lr", "none", "reviewer_safe_internal_locked"),
        ("v2_multimodal_rf", "multimodal_structure_processing", "multimodal", "rf", "none", "reviewer_safe_internal_locked"),
        ("v2_source_balanced_cf_rf", "source_robust", "counterfactual", "rf", "source_balanced", "reviewer_safe_internal_locked"),
        ("v2_groupdro_proxy_cf_lr", "source_robust", "counterfactual", "lr", "groupdro_proxy", "reviewer_safe_internal_locked"),
        ("v2_class_balanced_cf_plm_hgb", "ranking_proxy", "counterfactual_plm", "hgb", "class_balanced", "reviewer_safe_internal_locked"),
    ]
    if "esm2_35m" in feature_sets:
        specs.extend([
            ("v2_esm2_35m_lr", "real_frozen_plm", "esm2_35m", "lr", "none", "reviewer_safe_internal_locked"),
            ("v2_cf_esm2_35m_lr", "counterfactual_real_plm", "counterfactual_esm2_35m", "lr", "none", "reviewer_safe_internal_locked"),
            ("v2_multimodal_esm2_35m_lr", "multimodal_real_plm", "multimodal_esm2_35m", "lr", "none", "reviewer_safe_internal_locked"),
            ("v2_source_balanced_cf_esm2_35m_lr", "source_robust_real_plm", "counterfactual_esm2_35m", "lr", "source_balanced", "reviewer_safe_internal_locked"),
            ("v2_groupdro_proxy_cf_esm2_35m_lr", "source_robust_real_plm", "counterfactual_esm2_35m", "lr", "groupdro_proxy", "reviewer_safe_internal_locked"),
        ])

    all_pred = [import_v1_predictions(reg_ids)]
    model_errors = []
    split_files = sorted((OUT / "splits").glob("*.tsv"))
    for split_file in split_files:
        split_df = pd.read_csv(split_file, sep="\t")
        split_name = split_file.stem
        for fold_id, fold in split_df.groupby("fold_id"):
            train_ids = fold.loc[fold["role"].eq("train"), "row_id"].astype(str).tolist()
            test_ids = fold.loc[fold["role"].eq("test"), "row_id"].astype(str).tolist()
            train_meta = reg[reg["row_id"].astype(str).isin(train_ids)].copy()
            test_meta = reg[reg["row_id"].astype(str).isin(test_ids)].copy()
            if train_meta["label_binary"].nunique() < 2 or test_meta.empty:
                continue
            for model_name, family, fset, kind, weight_mode, status in specs:
                feats = feature_sets[fset]
                train = train_meta[["row_id", "label_binary", "source_dataset", "hla_supertype"]].merge(feats, on="row_id", how="left").fillna(0)
                test = test_meta[["row_id", "label_binary", "source_dataset", "hla_supertype"]].merge(feats, on="row_id", how="left").fillna(0)
                x_cols = [c for c in train.columns if c not in {"row_id", "label_binary", "source_dataset", "hla_supertype"}]
                model = make_estimator(kind)
                weights = train_weights(train_meta, weight_mode)
                try:
                    fit_model(model, train[x_cols], train["label_binary"].values, weights)
                    score = predict_score(model, test[x_cols])
                    out = pd.DataFrame({
                        "split_name": split_name,
                        "fold_id": fold_id,
                        "row_id": test["row_id"].astype(str).values,
                        "label": test["label_binary"].astype(int).values,
                        "score": score,
                        "model_name": model_name,
                        "model_family": family,
                        "claim_status": status,
                        "source_script": "cross_neo_v2_train_models",
                    })
                    all_pred.append(out)
                    if model_name in {"v2_counterfactual_lr", "v2_cf_plm_rf"}:
                        train_score = predict_score(model, train[x_cols])
                        rn = train_percentile(train_score, train_meta["hla_supertype"], score, test_meta["hla_supertype"])
                        all_pred.append(out.assign(
                            score=rn,
                            model_name=model_name + "_hla_ranknorm",
                            model_family="source_robust_rank_normalization",
                        ))
                except Exception as exc:
                    model_errors.append({"split_name": split_name, "fold_id": fold_id, "model_name": model_name, "error": repr(exc)})

    pred = pd.concat([x for x in all_pred if x is not None and not x.empty], ignore_index=True)
    # Deterministic score-level MoE/gate: no label access, no test-fold tuning.
    moe_rows = []
    base_names = ["v2_counterfactual_rf", "v2_plm_lr_frozen_hash_pilot", "v2_multimodal_lr"]
    for (split, fold), g in pred[pred["model_name"].isin(base_names)].groupby(["split_name", "fold_id"]):
        wide = g.pivot_table(index=["row_id", "label"], columns="model_name", values="score", aggfunc="mean").reset_index()
        if set(base_names).issubset(wide.columns):
            score = 0.55 * wide["v2_counterfactual_rf"] + 0.25 * wide["v2_multimodal_lr"] + 0.20 * wide["v2_plm_lr_frozen_hash_pilot"]
            moe_rows.append(pd.DataFrame({
                "split_name": split,
                "fold_id": fold,
                "row_id": wide["row_id"].astype(str),
                "label": wide["label"].astype(int),
                "score": np.asarray(score, dtype=float),
                "model_name": "v2_rule_moe_cf_plm_structure",
                "model_family": "mixture_of_experts_prespecified",
                "claim_status": "reviewer_safe_internal_locked",
                "source_script": "cross_neo_v2_train_models",
            }))
    if moe_rows:
        pred = pd.concat([pred] + moe_rows, ignore_index=True)
    pred["score"] = pd.to_numeric(pred["score"], errors="coerce").fillna(0).clip(0, 1)
    pred = add_ranks(pred)
    safe_parquet(pred, OUT / "predictions/all_predictions.parquet")
    pred.to_csv(OUT / "predictions/all_predictions.tsv", sep="\t", index=False, na_rep="NA")
    pd.DataFrame(model_errors).to_csv(OUT / "logs/model_training_errors.tsv", sep="\t", index=False, na_rep="NA")
    (OUT / "models").mkdir(exist_ok=True)
    (OUT / "models/model_manifest.json").write_text(json.dumps({
        "seed": SEED,
        "n_prediction_rows": int(len(pred)),
        "n_models": int(pred["model_name"].nunique()),
        "n_splits": int(pred["split_name"].nunique()),
        "note": "v2 first pass: frozen hashed PLM pilot; public predictor scores excluded; QK imported only as v1 diagnostic/fusion baseline.",
    }, indent=2) + "\n")
    print(f"[v2-train] predictions={len(pred)} models={pred['model_name'].nunique()} splits={pred['split_name'].nunique()} errors={len(model_errors)}")


if __name__ == "__main__":
    main()
