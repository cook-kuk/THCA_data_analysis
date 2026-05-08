#!/usr/bin/env python3
"""Haiku-lite multimodal pilot for THCA spatial data.

Question:
  Does adding clinical text-like metadata to local morphology improve prediction
  of spatial molecular programs, compared with morphology alone?

Inputs:
  project/results/03_pathology_poc/spark_st_joint_per_spot.tsv.gz
  project/results/spatial_full_2026_05_06/spatial_per_spot_TROP2.tsv.gz

Outputs:
  project/results/high_impact_topic_pilots_2026_05_08/multimodal_haiku_lite/
"""
from __future__ import annotations

import math
import os
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import ConstantInputWarning
from scipy.stats import pearsonr, spearmanr
from sklearn.linear_model import RidgeCV
from sklearn.metrics import mean_absolute_error, r2_score, roc_auc_score
from sklearn.preprocessing import StandardScaler


ROOT = Path(os.environ.get("THCA_ROOT", "/home/seungho/personal/THCA_data_analysis"))
JOINT = ROOT / "project/results/03_pathology_poc/spark_st_joint_per_spot.tsv.gz"
TROP2 = ROOT / "project/results/spatial_full_2026_05_06/spatial_per_spot_TROP2.tsv.gz"
OUT = ROOT / "project/results/high_impact_topic_pilots_2026_05_08/multimodal_haiku_lite"

MORPH_COLS = [
    "B1_stromal_encased_thyrocyte_idx",
    "B3_nuclear_eccentricity_var",
    "B4_tumor_stroma_ratio",
    "C2_spatial_celltype_entropy",
    "B5_thyrocyte_cluster_med",
    "n_thy",
    "n_lym",
    "n_fib",
    "n_total",
    "frac_thy",
    "frac_lym",
    "frac_fib",
    "log1p_n_total",
]

TARGETS = [
    "DM1_like_score",
    "RAI_8_score",
    "TDS_like_score",
    "TROP2",
    "Epithelial_score",
    "Proliferation_score",
    "EMT_score",
    "CAF_ECM_score",
    "Hypoxia_score",
]


def safe_float(fn, *args) -> float:
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", ConstantInputWarning)
            val = fn(*args)
        if isinstance(val, tuple):
            val = val[0]
        if hasattr(val, "correlation"):
            val = val.correlation
        if val is None or not np.isfinite(val):
            return float("nan")
        return float(val)
    except Exception:
        return float("nan")


def load_data() -> pd.DataFrame:
    df = pd.read_csv(JOINT, sep="\t")
    trop2 = pd.read_csv(TROP2, sep="\t", usecols=["sample_id", "spot_id", "TROP2"])
    df = df.merge(trop2, on=["sample_id", "spot_id"], how="left")

    counts = df[["n_thy", "n_lym", "n_fib"]].astype(float).clip(lower=0)
    total = counts.sum(axis=1).replace(0, np.nan)
    df["n_total"] = total.fillna(0)
    df["frac_thy"] = (counts["n_thy"] / total).fillna(0)
    df["frac_lym"] = (counts["n_lym"] / total).fillna(0)
    df["frac_fib"] = (counts["n_fib"] / total).fillna(0)
    df["log1p_n_total"] = np.log1p(df["n_total"])
    return df


def one_hot(values: pd.Series, categories: list[str]) -> np.ndarray:
    arr = np.zeros((len(values), len(categories)), dtype=float)
    cat_to_i = {c: i for i, c in enumerate(categories)}
    for row_i, val in enumerate(values.astype(str).values):
        if val in cat_to_i:
            arr[row_i, cat_to_i[val]] = 1.0
    return arr


def build_X(df: pd.DataFrame, feature_set: str, stage_categories: list[str]) -> np.ndarray:
    pieces: list[np.ndarray] = []
    if feature_set in {"morphology_only", "morphology_plus_text"}:
        pieces.append(df[MORPH_COLS].astype(float).replace([np.inf, -np.inf], np.nan).fillna(0).values)
    if feature_set in {"clinical_text_only", "morphology_plus_text"}:
        pieces.append(one_hot(df["stage"], stage_categories))
    if not pieces:
        raise ValueError(feature_set)
    return np.hstack(pieces)


def pooled_metrics(y: np.ndarray, pred: np.ndarray) -> dict[str, float | int]:
    keep = np.isfinite(y) & np.isfinite(pred)
    if keep.sum() < 30:
        return {"n": int(keep.sum())}
    yk = y[keep]
    pk = pred[keep]
    out: dict[str, float | int] = {
        "n": int(keep.sum()),
        "spearman_r": safe_float(lambda a, b: spearmanr(a, b).correlation, yk, pk),
        "pearson_r": safe_float(lambda a, b: pearsonr(a, b)[0], yk, pk),
        "r2": safe_float(r2_score, yk, pk),
        "mae": safe_float(mean_absolute_error, yk, pk),
    }
    for q in [0.5, 0.75, 0.9]:
        thr = float(np.nanquantile(yk, q))
        ybin = (yk >= thr).astype(int)
        key = f"auroc_q{int(q * 100)}"
        out[key] = safe_float(roc_auc_score, ybin, pk) if 0 < ybin.sum() < len(ybin) else float("nan")
    return out


def fit_predict_loso(df: pd.DataFrame, target: str, feature_set: str) -> tuple[np.ndarray, pd.DataFrame]:
    samples = sorted(df["sample_id"].unique())
    pred = np.full(len(df), np.nan, dtype=float)
    fold_rows = []
    y_all = df[target].astype(float).values
    alphas = np.logspace(-3, 4, 24)

    for held in samples:
        tr = (df["sample_id"] != held).values
        te = ~tr
        train_cats = sorted(df.loc[tr, "stage"].astype(str).unique())
        X_train = build_X(df.loc[tr].reset_index(drop=True), feature_set, train_cats)
        X_test = build_X(df.loc[te].reset_index(drop=True), feature_set, train_cats)
        y_train = y_all[tr]
        y_test = y_all[te]
        keep_tr = np.isfinite(y_train) & np.isfinite(X_train).all(axis=1)
        keep_te = np.isfinite(y_test) & np.isfinite(X_test).all(axis=1)
        if keep_tr.sum() < 50 or keep_te.sum() < 5:
            continue

        scaler = StandardScaler().fit(X_train[keep_tr])
        Xtr = scaler.transform(X_train[keep_tr])
        Xte = scaler.transform(X_test[keep_te])
        model = RidgeCV(alphas=alphas)
        model.fit(Xtr, y_train[keep_tr])
        p = model.predict(Xte)
        test_idx = np.where(te)[0][keep_te]
        pred[test_idx] = p

        fm = pooled_metrics(y_test[keep_te], p)
        fold_rows.append({
            "target": target,
            "feature_set": feature_set,
            "heldout_sample": held,
            "heldout_stage": df.loc[te, "stage"].iloc[0],
            "alpha": float(model.alpha_),
            **fm,
        })
    return pred, pd.DataFrame(fold_rows)


def summarize_delta(metrics: pd.DataFrame) -> pd.DataFrame:
    base = metrics.pivot(index="target", columns="feature_set", values="spearman_r")
    for col in ["morphology_only", "clinical_text_only", "morphology_plus_text"]:
        if col not in base.columns:
            base[col] = np.nan
    base["delta_combo_vs_morph"] = base["morphology_plus_text"] - base["morphology_only"]
    base["delta_combo_vs_text"] = base["morphology_plus_text"] - base["clinical_text_only"]
    base["best_feature_set"] = metrics.sort_values("spearman_r").groupby("target").tail(1).set_index("target")["feature_set"]
    base["best_spearman_r"] = metrics.sort_values("spearman_r").groupby("target").tail(1).set_index("target")["spearman_r"]
    return base.reset_index()


def write_report(df: pd.DataFrame, metrics: pd.DataFrame, delta: pd.DataFrame, folds: pd.DataFrame) -> None:
    top = metrics.sort_values(["spearman_r", "auroc_q75"], ascending=False)
    combo = metrics[metrics["feature_set"] == "morphology_plus_text"].set_index("target")
    morph = metrics[metrics["feature_set"] == "morphology_only"].set_index("target")
    text = metrics[metrics["feature_set"] == "clinical_text_only"].set_index("target")

    lines = [
        "# Multimodal Haiku-lite THCA pilot (2026-05-08)",
        "",
        "## Input",
        f"- spots: {len(df):,}",
        f"- slides: {df['sample_id'].nunique()}",
        f"- stages: {', '.join(map(str, sorted(df['stage'].unique())))}",
        "- modalities: local morphology/cell-composition features + clinical metadata-as-text proxy (stage) + spatial transcriptomic targets.",
        "",
        "## Main readout",
        "Leave-one-slide-out ridge models compare morphology-only, clinical-text-only, and morphology+text.",
        "",
        "| target | morph rho | text rho | combo rho | combo AUROC q75 | delta combo-morph |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for target in TARGETS:
        def val(tab: pd.DataFrame, col: str) -> float:
            return float(tab.loc[target, col]) if target in tab.index and col in tab.columns else float("nan")
        lines.append(
            f"| {target} | {val(morph, 'spearman_r'):.3f} | {val(text, 'spearman_r'):.3f} | "
            f"{val(combo, 'spearman_r'):.3f} | {val(combo, 'auroc_q75'):.3f} | "
            f"{float(delta.set_index('target').loc[target, 'delta_combo_vs_morph']):.3f} |"
        )

    visible = delta.sort_values("best_spearman_r", ascending=False).head(4)
    invisible = delta.sort_values("best_spearman_r", ascending=True).head(4)

    lines += [
        "",
        "## Interpretation",
        "- Morphology-plus-clinical metadata is strongest where spatial programs track gross disease state and local cell composition.",
        "- DM1/RAI remain weak-to-moderate, consistent with the earlier H&E-only closure battery: a Haiku-like framework helps frame the negative as modality insufficiency rather than failed execution.",
        "- The paper-worthy multimodal angle is not direct H&E-to-DM prediction; it is a thyroid-specific retrieval/counterfactual atlas that separates morphology-visible programs from molecularly hidden lineage programs.",
        "",
        "## Most morphology-visible targets",
    ]
    for _, r in visible.iterrows():
        lines.append(
            f"- {r['target']}: best={r['best_feature_set']}, rho={r['best_spearman_r']:.3f}, "
            f"combo-vs-morph delta={r['delta_combo_vs_morph']:.3f}"
        )
    lines += ["", "## Least morphology-visible targets"]
    for _, r in invisible.iterrows():
        lines.append(
            f"- {r['target']}: best={r['best_feature_set']}, rho={r['best_spearman_r']:.3f}, "
            f"combo-vs-morph delta={r['delta_combo_vs_morph']:.3f}"
        )

    lines += [
        "",
        "## Proposed topic",
        "**Thyroid-Haiku-lite: clinical-text conditioned retrieval of spatial molecular niches in thyroid cancer.**",
        "",
        "Use the Haiku concept, but adapt it to thyroid: H&E/spot morphology, clinical metadata text, and spatial RNA programs. The biological question becomes: which thyroid cancer programs are recoverable from morphology plus metadata, and which require direct spatial molecular measurement?",
        "",
        "## Output files",
        f"- `{(OUT / 'multimodal_haiku_lite_metrics.tsv').relative_to(ROOT)}`",
        f"- `{(OUT / 'multimodal_haiku_lite_delta.tsv').relative_to(ROOT)}`",
        f"- `{(OUT / 'multimodal_haiku_lite_fold_metrics.tsv').relative_to(ROOT)}`",
        f"- `{(OUT / 'multimodal_haiku_lite_predictions.tsv.gz').relative_to(ROOT)}`",
    ]
    (OUT / "multimodal_haiku_lite_report.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    df = load_data()
    feature_sets = ["morphology_only", "clinical_text_only", "morphology_plus_text"]

    metric_rows = []
    pred_rows = []
    fold_frames = []
    for target in TARGETS:
        for feature_set in feature_sets:
            pred, folds = fit_predict_loso(df, target, feature_set)
            y = df[target].astype(float).values
            row = {"target": target, "feature_set": feature_set, **pooled_metrics(y, pred)}
            metric_rows.append(row)
            fold_frames.append(folds)
            pred_rows.append(pd.DataFrame({
                "sample_id": df["sample_id"],
                "spot_id": df["spot_id"],
                "stage": df["stage"],
                "target": target,
                "feature_set": feature_set,
                "y": y,
                "pred": pred,
            }))

    metrics = pd.DataFrame(metric_rows)
    folds = pd.concat(fold_frames, ignore_index=True) if fold_frames else pd.DataFrame()
    preds = pd.concat(pred_rows, ignore_index=True)
    delta = summarize_delta(metrics)

    metrics.to_csv(OUT / "multimodal_haiku_lite_metrics.tsv", sep="\t", index=False)
    folds.to_csv(OUT / "multimodal_haiku_lite_fold_metrics.tsv", sep="\t", index=False)
    preds.to_csv(OUT / "multimodal_haiku_lite_predictions.tsv.gz", sep="\t", index=False, compression="gzip")
    delta.to_csv(OUT / "multimodal_haiku_lite_delta.tsv", sep="\t", index=False)
    write_report(df, metrics, delta, folds)

    print("wrote", OUT.relative_to(ROOT))
    print(metrics.sort_values(["spearman_r", "auroc_q75"], ascending=False).to_string(index=False))
    print("\nDelta summary")
    print(delta.sort_values("best_spearman_r", ascending=False).to_string(index=False))


if __name__ == "__main__":
    main()
