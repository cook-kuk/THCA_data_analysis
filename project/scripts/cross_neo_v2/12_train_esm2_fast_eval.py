#!/usr/bin/env python3
"""Fast RunPod evaluation for ESM2-35M CROSS-Neo v2 candidates."""

from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from common import OUT, SEED, ensure_dirs, metrics, read_table


def read_feature(name: str) -> pd.DataFrame:
    path = OUT / "features" / name
    try:
        return pd.read_parquet(path)
    except Exception:
        pkl = path.with_suffix(path.suffix + ".pkl")
        if pkl.exists():
            return pd.read_pickle(pkl)
    return pd.read_csv(path.with_suffix(".tsv"), sep="\t")


def merge_features(frames: list[pd.DataFrame]) -> pd.DataFrame:
    x = frames[0].copy()
    for f in frames[1:]:
        x = x.merge(f, on="row_id", how="left")
    for c in x.columns:
        if c != "row_id":
            x[c] = pd.to_numeric(x[c], errors="coerce")
    return x.fillna(0)


def source_weights(meta: pd.DataFrame, mode: str) -> np.ndarray | None:
    if mode == "none":
        return None
    y = meta["label_binary"].astype(int)
    w = pd.Series(1.0, index=meta.index)
    if mode in {"source", "groupdro"}:
        vc = meta["source_dataset"].fillna("unknown").value_counts()
        w *= meta["source_dataset"].fillna("unknown").map(lambda s: 1 / max(1, vc.get(s, 1))).astype(float)
        w *= len(w) / max(1e-9, w.sum())
    if mode in {"class", "groupdro"}:
        vc = y.value_counts()
        w *= y.map(lambda yy: len(y) / (2 * max(1, vc.get(yy, 1)))).astype(float)
    return np.clip(np.asarray(w / max(1e-9, w.mean()), dtype=float), 0.2, 8.0)


def fit_lr(x, y, weights=None):
    model = Pipeline([
        ("imputer", SimpleImputer(strategy="median")),
        ("scaler", StandardScaler()),
        ("clf", LogisticRegression(C=0.5, class_weight="balanced", max_iter=3000, solver="liblinear", random_state=SEED)),
    ])
    if weights is not None:
        try:
            model.fit(x, y, clf__sample_weight=weights)
            return model
        except Exception:
            pass
    model.fit(x, y)
    return model


def add_ranks(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["rank"] = df.groupby(["split_name", "fold_id", "model_name"])["score"].rank(method="first", ascending=False)
    n = df.groupby(["split_name", "fold_id", "model_name"])["row_id"].transform("size")
    df["rank_pct"] = df["rank"] / n.clip(lower=1)
    return df


def main() -> None:
    ensure_dirs()
    tag = os.environ.get("CROSS_NEO_ESM2_TAG", "35m").lower()
    reg = pd.read_csv(OUT / "canonical_registry.tsv", sep="\t")
    reg["label_binary"] = pd.to_numeric(reg["label_binary"], errors="coerce").fillna(0).astype(int)
    cf = read_feature("counterfactual_features.parquet")
    esm2 = read_feature(f"esm2_{tag}_features.parquet")
    struct = read_feature("structure_features.parquet")
    proc = read_feature("processing_context_features.parquet")
    sets = {
        "esm2": esm2,
        "cf_esm2": merge_features([cf, esm2]),
        "multi_esm2": merge_features([cf, esm2, struct, proc]),
    }
    specs = [
        (f"v2_esm2_{tag}_lr_fast", "real_frozen_plm_fast", "esm2", "none"),
        (f"v2_cf_esm2_{tag}_lr_fast", "counterfactual_real_plm_fast", "cf_esm2", "none"),
        (f"v2_multimodal_esm2_{tag}_lr_fast", "multimodal_real_plm_fast", "multi_esm2", "none"),
        (f"v2_source_balanced_cf_esm2_{tag}_lr_fast", "source_robust_real_plm_fast", "cf_esm2", "source"),
        (f"v2_groupdro_proxy_cf_esm2_{tag}_lr_fast", "source_robust_real_plm_fast", "cf_esm2", "groupdro"),
    ]
    preds = []
    for split_file in sorted((OUT / "splits").glob("*.tsv")):
        split = split_file.stem
        sdf = pd.read_csv(split_file, sep="\t")
        for fold, fdf in sdf.groupby("fold_id"):
            train_ids = set(fdf.loc[fdf["role"].eq("train"), "row_id"].astype(str))
            test_ids = set(fdf.loc[fdf["role"].eq("test"), "row_id"].astype(str))
            train_meta = reg[reg["row_id"].astype(str).isin(train_ids)].copy()
            test_meta = reg[reg["row_id"].astype(str).isin(test_ids)].copy()
            if train_meta["label_binary"].nunique() < 2 or test_meta.empty:
                continue
            for model_name, fam, set_name, wmode in specs:
                feat = sets[set_name]
                tr = train_meta[["row_id", "label_binary"]].merge(feat, on="row_id", how="left").fillna(0)
                te = test_meta[["row_id", "label_binary"]].merge(feat, on="row_id", how="left").fillna(0)
                xcols = [c for c in tr.columns if c not in {"row_id", "label_binary"}]
                model = fit_lr(tr[xcols], tr["label_binary"].values, source_weights(train_meta, wmode))
                score = model.predict_proba(te[xcols])[:, 1]
                preds.append(pd.DataFrame({
                    "split_name": split,
                    "fold_id": fold,
                    "row_id": te["row_id"].astype(str).values,
                    "label": te["label_binary"].astype(int).values,
                    "score": score,
                    "model_name": model_name,
                    "model_family": fam,
                    "claim_status": "reviewer_safe_internal_locked",
                    "source_script": "cross_neo_v2_esm2_fast_eval",
                }))
    pred = add_ranks(pd.concat(preds, ignore_index=True))
    pred.to_csv(OUT / f"predictions/esm2_{tag}_fast_predictions.tsv", sep="\t", index=False)
    rows = []
    for (split, model), g in pred.groupby(["split_name", "model_name"]):
        rows.append({"split_name": split, "model_name": model, "model_family": g["model_family"].iloc[0], "claim_status": g["claim_status"].iloc[0], **metrics(g["label"].values, g["score"].values)})
    met = pd.DataFrame(rows)
    met.to_csv(OUT / f"metrics/esm2_{tag}_fast_metrics.tsv", sep="\t", index=False)
    primary = ["exact_peptide_hla_holdout", "near_peptide_cluster_holdout", "hla_stratified_group_5fold", "hla_supertype_heldout"]
    best = met[met["split_name"].isin(primary)].sort_values(["split_name", "AUPRC", "top10_precision"], ascending=[True, False, False]).groupby("split_name").head(3)
    report = [
        f"# ESM2-{tag} Fast RunPod Evaluation",
        "",
        best[["split_name", "model_name", "n", "n_pos", "prevalence", "AUPRC", "AUROC", "top10_precision", "top20_precision"]].to_markdown(index=False),
        "",
        f"Claim boundary: frozen ESM2-{tag} features, train-fold-only scaling/model fitting, no public predictor scores.",
    ]
    (OUT / f"ESM2_{tag}_fast_runpod_report.md").write_text("\n".join(report) + "\n")
    print(f"[v2-esm2-{tag}-fast] pred={len(pred)} metric_rows={len(met)}")
    print(best[["split_name", "model_name", "AUPRC", "top10_precision"]].to_string(index=False))


if __name__ == "__main__":
    main()
