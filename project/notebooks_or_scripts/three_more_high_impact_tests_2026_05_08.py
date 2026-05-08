#!/usr/bin/env python3
"""Three additional high-impact stress tests for T00/T02/T12.

Test 1. T00 CNV robustness:
  Does the cBioPortal ARMDRIVER signature still identify DM2 after stratifying
  by morphology/clinical covariates, and does it add cross-validated signal
  over those covariates?

Test 2. T12 retrieval:
  Does evidence-grounded KNN retrieval recover spatial targets from morphology
  and clinical text, without target-specific regression?

Test 3. T12 counterfactual text perturbation:
  Holding morphology fixed, does changing only the stage text alter retrieved
  spatial programs in the expected direction?
"""
from __future__ import annotations

import math
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, pearsonr, spearmanr
from sklearn.compose import ColumnTransformer
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import mean_absolute_error, r2_score, roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.neighbors import NearestNeighbors
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler


ROOT = Path(os.environ.get("THCA_ROOT", "/home/seungho/personal/THCA_data_analysis"))
BASE = ROOT / "project/results/high_impact_topic_pilots_2026_05_08"
DATA = BASE / "data_search"
OUT = BASE / "three_more_tests_2026_05_08"

CLASS6 = DATA / "cbio_class6_merged_sample_clinical.tsv"
JOINT = ROOT / "project/results/03_pathology_poc/spark_st_joint_per_spot.tsv.gz"
TROP2 = ROOT / "project/results/spatial_full_2026_05_06/spatial_per_spot_TROP2.tsv.gz"

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
    "TROP2",
    "DM1_like_score",
    "RAI_8_score",
    "TDS_like_score",
    "CAF_ECM_score",
    "Hypoxia_score",
]

STAGES = ["PT", "PTC", "LPTC", "ATC"]


def fmt(x: float) -> str:
    if x is None or not np.isfinite(x):
        return "NA"
    if abs(x) >= 1000 or (abs(x) < 0.001 and x != 0):
        return f"{x:.3e}"
    return f"{x:.3f}"


def safe_corr(fn, a, b) -> float:
    try:
        r = fn(a, b)
        if isinstance(r, tuple):
            r = r[0]
        if hasattr(r, "correlation"):
            r = r.correlation
        return float(r) if np.isfinite(r) else float("nan")
    except Exception:
        return float("nan")


def pooled_metrics(y: np.ndarray, pred: np.ndarray) -> dict[str, float | int]:
    keep = np.isfinite(y) & np.isfinite(pred)
    out: dict[str, float | int] = {"n": int(keep.sum())}
    if keep.sum() < 30:
        return out
    yk = y[keep]
    pk = pred[keep]
    out.update({
        "spearman_r": safe_corr(lambda x, z: spearmanr(x, z).correlation, yk, pk),
        "pearson_r": safe_corr(lambda x, z: pearsonr(x, z)[0], yk, pk),
        "r2": safe_corr(r2_score, yk, pk),
        "mae": safe_corr(mean_absolute_error, yk, pk),
    })
    for q in [0.5, 0.75, 0.9]:
        thr = float(np.nanquantile(yk, q))
        ybin = (yk >= thr).astype(int)
        out[f"auroc_q{int(q * 100)}"] = (
            safe_corr(roc_auc_score, ybin, pk)
            if 0 < ybin.sum() < len(ybin)
            else float("nan")
        )
    return out


def load_spatial() -> pd.DataFrame:
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


def one_hot(values: pd.Series | list[str], categories: list[str]) -> np.ndarray:
    vals = pd.Series(values).astype(str).values
    arr = np.zeros((len(vals), len(categories)), dtype=float)
    idx = {c: i for i, c in enumerate(categories)}
    for i, val in enumerate(vals):
        if val in idx:
            arr[i, idx[val]] = 1.0
    return arr


def row_l2(x: np.ndarray) -> np.ndarray:
    norm = np.linalg.norm(x, axis=1, keepdims=True)
    return x / np.maximum(norm, 1e-8)


def retrieval_features(
    train: pd.DataFrame,
    query: pd.DataFrame,
    mode: str,
    alpha: float = 0.6,
    query_stage_override: str | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    stage_cats = [s for s in STAGES if s in set(train["stage"].astype(str))]
    if not stage_cats:
        stage_cats = sorted(train["stage"].astype(str).unique())

    morph_scaler = StandardScaler().fit(train[MORPH_COLS].astype(float).fillna(0).values)
    m_train = row_l2(morph_scaler.transform(train[MORPH_COLS].astype(float).fillna(0).values))
    m_query = row_l2(morph_scaler.transform(query[MORPH_COLS].astype(float).fillna(0).values))

    train_text = one_hot(train["stage"], stage_cats)
    q_stage = [query_stage_override] * len(query) if query_stage_override else query["stage"]
    query_text = one_hot(q_stage, stage_cats)
    text_scaler = StandardScaler().fit(train_text)
    t_train = row_l2(text_scaler.transform(train_text))
    t_query = row_l2(text_scaler.transform(query_text))

    if mode == "morphology_only":
        return m_train, m_query
    if mode == "clinical_text_only":
        return t_train, t_query
    if mode == "fusion_a06":
        return (
            np.hstack([math.sqrt(alpha) * m_train, math.sqrt(1 - alpha) * t_train]),
            np.hstack([math.sqrt(alpha) * m_query, math.sqrt(1 - alpha) * t_query]),
        )
    if mode == "fusion_a08":
        alpha = 0.8
        return (
            np.hstack([math.sqrt(alpha) * m_train, math.sqrt(1 - alpha) * t_train]),
            np.hstack([math.sqrt(alpha) * m_query, math.sqrt(1 - alpha) * t_query]),
        )
    raise ValueError(mode)


def knn_predict(X_train: np.ndarray, y_train: np.ndarray, X_query: np.ndarray, k: int = 25) -> np.ndarray:
    keep = np.isfinite(y_train) & np.isfinite(X_train).all(axis=1)
    Xg = X_train[keep]
    yg = y_train[keep]
    kk = min(k, len(yg))
    nn = NearestNeighbors(n_neighbors=kk, metric="euclidean")
    nn.fit(Xg)
    dist, ind = nn.kneighbors(X_query)
    weights = 1.0 / np.maximum(dist, 1e-3)
    vals = yg[ind]
    return (weights * vals).sum(axis=1) / np.maximum(weights.sum(axis=1), 1e-8)


def fisher_row(df: pd.DataFrame, label: str, flag_col: str = "cbio_ARMDRIVER_signature_any") -> dict:
    sub = df[df["v17_dark_cluster"].isin(["DM1", "DM2"])].copy()
    dm2 = sub[sub["v17_dark_cluster"] == "DM2"][flag_col].astype(bool)
    dm1 = sub[sub["v17_dark_cluster"] == "DM1"][flag_col].astype(bool)
    if len(dm2) < 3 or len(dm1) < 3:
        return {
            "test": label,
            "n_DM2": len(dm2),
            "n_DM1": len(dm1),
            "DM2_rate": np.nan,
            "DM1_rate": np.nan,
            "OR": np.nan,
            "p": np.nan,
        }
    table = [[int(dm2.sum()), int((~dm2).sum())], [int(dm1.sum()), int((~dm1).sum())]]
    OR, p = fisher_exact(table)
    return {
        "test": label,
        "n_DM2": len(dm2),
        "n_DM1": len(dm1),
        "DM2_rate": float(dm2.mean()),
        "DM1_rate": float(dm1.mean()),
        "OR": float(OR),
        "p": float(p),
    }


def test1_cnv_robustness() -> tuple[pd.DataFrame, pd.DataFrame]:
    df = pd.read_csv(CLASS6, sep="\t")
    df = df[df["v17_dark_cluster"].isin(["DM1", "DM2"])].copy()
    df["is_DM2"] = df["v17_dark_cluster"].eq("DM2").astype(int)
    df["cbio_ARMDRIVER_signature_any"] = df["cbio_ARMDRIVER_signature_any"].astype(bool)
    df["follicular_pct_num"] = pd.to_numeric(df["follicular_pct_num"], errors="coerce")
    df["PURITY"] = pd.to_numeric(df["PURITY"], errors="coerce")
    df["age"] = pd.to_numeric(df["age"], errors="coerce")

    subset_rows = [fisher_row(df, "full_Class6_DM2_vs_DM1")]
    for col in ["histology_subtype", "BRAFV600E_RAS", "stage", "ARM_SCNA_CLUSTER"]:
        for val in sorted(df[col].dropna().astype(str).unique()):
            sub = df[df[col].astype(str).eq(val)]
            if (sub["v17_dark_cluster"].eq("DM2").sum() >= 5 and
                    sub["v17_dark_cluster"].eq("DM1").sum() >= 5):
                subset_rows.append(fisher_row(sub, f"{col}={val}"))
    for col, label in [("follicular_pct_num", "follicular_pct"), ("PURITY", "purity")]:
        med = float(df[col].median())
        subset_rows.append(fisher_row(df[df[col] >= med], f"{label}>=median({med:.3f})"))
        subset_rows.append(fisher_row(df[df[col] < med], f"{label}<median({med:.3f})"))

    subset = pd.DataFrame(subset_rows).sort_values(["p", "test"], na_position="last")

    cat_cols = ["histology_subtype", "stage", "BRAFV600E_RAS"]
    num_cols = ["follicular_pct_num", "PURITY", "age"]
    y = df["is_DM2"].values
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=13)

    def model_auc(use_sig: bool) -> list[float]:
        cols = cat_cols + num_cols + (["cbio_ARMDRIVER_signature_any"] if use_sig else [])
        X = df[cols].copy()
        if use_sig:
            X["cbio_ARMDRIVER_signature_any"] = X["cbio_ARMDRIVER_signature_any"].astype(int)
        num = num_cols + (["cbio_ARMDRIVER_signature_any"] if use_sig else [])
        pre = ColumnTransformer([
            ("num", Pipeline([("impute", SimpleImputer(strategy="median")), ("scale", StandardScaler())]), num),
            ("cat", Pipeline([("impute", SimpleImputer(strategy="most_frequent")), ("oh", OneHotEncoder(handle_unknown="ignore"))]), cat_cols),
        ])
        clf = Pipeline([
            ("pre", pre),
            ("lr", LogisticRegression(max_iter=3000, class_weight="balanced", C=1.0)),
        ])
        aucs = []
        for tr, te in cv.split(X, y):
            clf.fit(X.iloc[tr], y[tr])
            proba = clf.predict_proba(X.iloc[te])[:, 1]
            aucs.append(float(roc_auc_score(y[te], proba)))
        return aucs

    covar_auc = model_auc(use_sig=False)
    sig_auc = model_auc(use_sig=True)
    auc = pd.DataFrame({
        "fold": np.arange(1, len(covar_auc) + 1),
        "covariates_only_auc": covar_auc,
        "covariates_plus_signature_auc": sig_auc,
    })
    auc["delta_auc"] = auc["covariates_plus_signature_auc"] - auc["covariates_only_auc"]
    return subset, auc


def test2_retrieval(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    modes = ["morphology_only", "clinical_text_only", "fusion_a06", "fusion_a08"]
    pred_frames = []
    metric_rows = []
    for target in TARGETS:
        y_all = df[target].astype(float).values
        for mode in modes:
            pred = np.full(len(df), np.nan)
            for held in sorted(df["sample_id"].unique()):
                tr = df["sample_id"].ne(held).values
                te = ~tr
                Xtr, Xte = retrieval_features(df.loc[tr], df.loc[te], mode)
                pred[np.where(te)[0]] = knn_predict(Xtr, y_all[tr], Xte, k=25)
            metric_rows.append({"target": target, "retrieval_mode": mode, **pooled_metrics(y_all, pred)})
            pred_frames.append(pd.DataFrame({
                "sample_id": df["sample_id"],
                "spot_id": df["spot_id"],
                "stage": df["stage"],
                "target": target,
                "retrieval_mode": mode,
                "y": y_all,
                "pred": pred,
            }))
    return pd.DataFrame(metric_rows), pd.concat(pred_frames, ignore_index=True)


def test3_counterfactual(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    row_frames = []
    for target in TARGETS:
        y_all = df[target].astype(float).values
        for held in sorted(df["sample_id"].unique()):
            tr = df["sample_id"].ne(held).values
            te = ~tr
            train = df.loc[tr]
            query = df.loc[te]
            for cf_stage in STAGES:
                Xtr, Xq = retrieval_features(train, query, "fusion_a06", alpha=0.6, query_stage_override=cf_stage)
                pred = knn_predict(Xtr, y_all[tr], Xq, k=25)
                row_frames.append(pd.DataFrame({
                    "sample_id": query["sample_id"].values,
                    "spot_id": query["spot_id"].values,
                    "actual_stage": query["stage"].values,
                    "counterfactual_stage": cf_stage,
                    "target": target,
                    "pred": pred,
                }))
    cf = pd.concat(row_frames, ignore_index=True)
    summary_rows = []
    for target, sub in cf.groupby("target"):
        piv = sub.pivot_table(index=["sample_id", "spot_id", "actual_stage"], columns="counterfactual_stage", values="pred")
        for stage in STAGES:
            summary_rows.append({
                "target": target,
                "comparison": f"mean_pred_if_{stage}",
                "mean": float(piv[stage].mean()),
                "median": float(piv[stage].median()),
                "n_spots": int(piv[stage].notna().sum()),
            })
        for a, b in [("PTC", "PT"), ("LPTC", "PT"), ("ATC", "PT"), ("LPTC", "PTC"), ("ATC", "PTC")]:
            diff = piv[a] - piv[b]
            summary_rows.append({
                "target": target,
                "comparison": f"{a}_minus_{b}",
                "mean": float(diff.mean()),
                "median": float(diff.median()),
                "n_spots": int(diff.notna().sum()),
            })
    return pd.DataFrame(summary_rows), cf


def write_report(test1: pd.DataFrame, auc: pd.DataFrame, retrieval: pd.DataFrame, cf_sum: pd.DataFrame) -> None:
    lines = [
        "# Three more high-impact tests (2026-05-08)",
        "",
        "## Test 1 — T00 CNV robustness",
        "",
        "Question: does the cBioPortal ARMDRIVER signature survive simple morphology/purity/clinical stratification, and does it add signal over covariates?",
        "",
        "| subset | n DM2 | n DM1 | DM2 rate | DM1 rate | OR | p |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for row in test1.head(12).itertuples(index=False):
        lines.append(
            f"| {row.test} | {row.n_DM2} | {row.n_DM1} | {fmt(row.DM2_rate)} | "
            f"{fmt(row.DM1_rate)} | {fmt(row.OR)} | {fmt(row.p)} |"
        )
    lines += [
        "",
        f"Cross-validated covariates-only AUC: {auc['covariates_only_auc'].mean():.3f}",
        f"Cross-validated covariates+CNV-signature AUC: {auc['covariates_plus_signature_auc'].mean():.3f}",
        f"Mean delta AUC: {auc['delta_auc'].mean():.3f}",
        "",
        "## Test 2 — Evidence-grounded retrieval",
        "",
        "Leave-one-slide-out KNN retrieval; no target-specific regression. `fusion_a06` uses 60% morphology and 40% stage text.",
        "",
        "| target | best mode | best rho | best AUROC q75 | morph rho | text rho | fusion_a06 rho |",
        "|---|---|---:|---:|---:|---:|---:|",
    ]
    for target, sub in retrieval.groupby("target"):
        best = sub.sort_values("spearman_r", ascending=False).iloc[0]
        pivot = sub.set_index("retrieval_mode")
        lines.append(
            f"| {target} | {best.retrieval_mode} | {fmt(best.spearman_r)} | {fmt(best.auroc_q75)} | "
            f"{fmt(pivot.loc['morphology_only', 'spearman_r'])} | "
            f"{fmt(pivot.loc['clinical_text_only', 'spearman_r'])} | "
            f"{fmt(pivot.loc['fusion_a06', 'spearman_r'])} |"
        )
    lines += [
        "",
        "## Test 3 — Haiku-style counterfactual stage text",
        "",
        "Holding morphology fixed, change only the stage text in fusion retrieval and summarize predicted target shifts.",
        "",
        "| target | PTC-PT | LPTC-PT | ATC-PT | LPTC-PTC | ATC-PTC |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for target, sub in cf_sum.groupby("target"):
        vals = sub.set_index("comparison")["mean"]
        lines.append(
            f"| {target} | {fmt(vals.get('PTC_minus_PT', np.nan))} | "
            f"{fmt(vals.get('LPTC_minus_PT', np.nan))} | {fmt(vals.get('ATC_minus_PT', np.nan))} | "
            f"{fmt(vals.get('LPTC_minus_PTC', np.nan))} | {fmt(vals.get('ATC_minus_PTC', np.nan))} |"
        )
    lines += [
        "",
        "## Bottom line",
        "",
        "- T00 remains the main biology flagship if the signature improves DM2 classification beyond covariates and remains enriched in key strata.",
        "- T12 is strengthened only if retrieval, not just regression, recovers TROP2 with morphology+text.",
        "- Counterfactual shifts should be presented as hypothesis-generating, matching Haiku's own caution.",
        "",
        "## Files",
        f"- `{(OUT / 'test1_cnv_stratified_robustness.tsv').relative_to(ROOT)}`",
        f"- `{(OUT / 'test1_cnv_covariate_auc.tsv').relative_to(ROOT)}`",
        f"- `{(OUT / 'test2_knn_retrieval_metrics.tsv').relative_to(ROOT)}`",
        f"- `{(OUT / 'test2_knn_retrieval_predictions.tsv.gz').relative_to(ROOT)}`",
        f"- `{(OUT / 'test3_counterfactual_stage_summary.tsv').relative_to(ROOT)}`",
        f"- `{(OUT / 'test3_counterfactual_stage_predictions.tsv.gz').relative_to(ROOT)}`",
    ]
    (OUT / "three_more_tests_report.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    t1, auc = test1_cnv_robustness()
    spatial = load_spatial()
    retrieval, retrieval_preds = test2_retrieval(spatial)
    cf_sum, cf_preds = test3_counterfactual(spatial)

    t1.to_csv(OUT / "test1_cnv_stratified_robustness.tsv", sep="\t", index=False)
    auc.to_csv(OUT / "test1_cnv_covariate_auc.tsv", sep="\t", index=False)
    retrieval.to_csv(OUT / "test2_knn_retrieval_metrics.tsv", sep="\t", index=False)
    retrieval_preds.to_csv(OUT / "test2_knn_retrieval_predictions.tsv.gz", sep="\t", index=False, compression="gzip")
    cf_sum.to_csv(OUT / "test3_counterfactual_stage_summary.tsv", sep="\t", index=False)
    cf_preds.to_csv(OUT / "test3_counterfactual_stage_predictions.tsv.gz", sep="\t", index=False, compression="gzip")
    write_report(t1, auc, retrieval, cf_sum)

    print("wrote", OUT.relative_to(ROOT))
    print("\nTEST1 top robustness")
    print(t1.head(12).to_string(index=False))
    print("\nTEST1 AUC")
    print(auc.describe().loc[["mean", "std"]].to_string())
    print("\nTEST2 retrieval")
    print(retrieval.sort_values(["target", "spearman_r"], ascending=[True, False]).to_string(index=False))
    print("\nTEST3 counterfactual summary")
    print(cf_sum[cf_sum["comparison"].str.contains("minus")].to_string(index=False))


if __name__ == "__main__":
    main()
