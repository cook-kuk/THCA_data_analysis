#!/usr/bin/env python3
"""Clean-CV scoring and lightweight fine-tuning for CROSS-Neo.

The goal is to produce real scores without inflating them with overlap-heavy rows.
We tune a small regularized stacker only on low-leakage clean rows, score it with
source/HLA held-out OOF protocols, then refit the selected configuration to score
all candidates. This is an analysis/fine-tuning layer, not a clinical model.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score
from sklearn.model_selection import GroupKFold, StratifiedGroupKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


ROOT = Path(__file__).resolve().parents[1]
BAR_DIR = ROOT / "project/results/clean_neobench_barneo_2026_05_09"
P0_DIR = ROOT / "project/results/cross_neo_kaggle_winner_playbook_2026_05_10/merged_p0_experiments"
OUT_DIR = ROOT / "project/results/cross_neo_kaggle_winner_playbook_2026_05_10/score_finetune_2026_05_10"

METHOD_SCORES = BAR_DIR / "clean_neobench_method_scores.tsv"
CANDIDATES = BAR_DIR / "barneo_stress_guarded_candidate_scores.tsv"
OVERLAP = BAR_DIR / "clean_neobench_overlap_flags.tsv"
SPLIT_METRICS = BAR_DIR / "clean_neobench_split_metrics.tsv"
METHOD_V2 = P0_DIR / "barneo_bma_v2_method_weights.tsv"
CANDIDATE_V2 = P0_DIR / "barneo_bma_v2_diversity_validity_scores.tsv"
WETLAB_V2 = P0_DIR / "cross_neo_wetlab_plate_v2.tsv"


@dataclass(frozen=True)
class ConfigResult:
    feature_set: str
    c_value: float
    protocol: str
    n_features: int
    n_rows: int
    n_pos: int
    auprc: float
    auroc: float
    brier: float
    ece: float
    top10_precision: float
    top20_precision: float
    oof_available: int


def read_tsv(path: Path, required: bool = False) -> pd.DataFrame:
    if not path.exists():
        if required:
            raise FileNotFoundError(path)
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t")


def write_tsv(df: pd.DataFrame, name: str) -> Path:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    path = OUT_DIR / name
    df.to_csv(path, sep="\t", index=False)
    return path


def numeric(df: pd.DataFrame, col: str, default: float = 0.0) -> pd.Series:
    if col not in df.columns:
        return pd.Series(default, index=df.index, dtype="float64")
    return pd.to_numeric(df[col], errors="coerce").fillna(default)


def text(df: pd.DataFrame, col: str, default: str = "") -> pd.Series:
    if col not in df.columns:
        return pd.Series(default, index=df.index, dtype="object")
    return df[col].fillna(default).astype(str)


def truthy(df: pd.DataFrame, col: str, default: bool = False) -> pd.Series:
    if col not in df.columns:
        return pd.Series(default, index=df.index, dtype="bool")
    raw = df[col]
    if raw.dtype == bool:
        return raw.fillna(default)
    filled = raw.astype("object").where(raw.notna(), default)
    return filled.astype(str).str.lower().str.strip().isin({"1", "true", "t", "yes", "y"})


def ece_score(y_true: np.ndarray, y_prob: np.ndarray, n_bins: int = 10) -> float:
    y_true = np.asarray(y_true).astype(float)
    y_prob = np.asarray(y_prob).astype(float)
    mask = np.isfinite(y_true) & np.isfinite(y_prob)
    y_true = y_true[mask]
    y_prob = y_prob[mask]
    if len(y_true) == 0:
        return float("nan")
    bins = np.linspace(0, 1, n_bins + 1)
    ece = 0.0
    for lo, hi in zip(bins[:-1], bins[1:]):
        bin_mask = (y_prob >= lo) & (y_prob < hi if hi < 1 else y_prob <= hi)
        if not bin_mask.any():
            continue
        ece += bin_mask.mean() * abs(y_true[bin_mask].mean() - y_prob[bin_mask].mean())
    return float(ece)


def topk_precision(y_true: np.ndarray, y_prob: np.ndarray, k: int) -> float:
    y_true = np.asarray(y_true).astype(float)
    y_prob = np.asarray(y_prob).astype(float)
    mask = np.isfinite(y_true) & np.isfinite(y_prob)
    y_true = y_true[mask]
    y_prob = y_prob[mask]
    if len(y_true) == 0:
        return float("nan")
    k = int(min(k, len(y_true)))
    order = np.argsort(-y_prob)[:k]
    return float(y_true[order].mean())


def metric_row(y_true: np.ndarray, y_prob: np.ndarray) -> dict[str, float]:
    y_true = np.asarray(y_true).astype(int)
    y_prob = np.asarray(y_prob).astype(float)
    mask = np.isfinite(y_prob)
    y_true = y_true[mask]
    y_prob = y_prob[mask]
    out = {
        "n_rows": int(len(y_true)),
        "n_pos": int(y_true.sum()) if len(y_true) else 0,
        "auprc": float("nan"),
        "auroc": float("nan"),
        "brier": float("nan"),
        "ece": float("nan"),
        "top10_precision": float("nan"),
        "top20_precision": float("nan"),
    }
    if len(y_true) == 0 or len(np.unique(y_true)) < 2:
        return out
    out["auprc"] = float(average_precision_score(y_true, y_prob))
    out["auroc"] = float(roc_auc_score(y_true, y_prob))
    out["brier"] = float(brier_score_loss(y_true, np.clip(y_prob, 0, 1)))
    out["ece"] = ece_score(y_true, np.clip(y_prob, 0, 1))
    out["top10_precision"] = topk_precision(y_true, y_prob, 10)
    out["top20_precision"] = topk_precision(y_true, y_prob, 20)
    return out


def load_metadata() -> pd.DataFrame:
    cand = read_tsv(CANDIDATES, required=True)
    overlap = read_tsv(OVERLAP)
    if not overlap.empty:
        cand = cand.merge(overlap, on="candidate_id", how="left", suffixes=("", "_overlap"))
        if "leakage_risk_level_overlap" in cand.columns:
            cand["leakage_risk_level"] = cand["leakage_risk_level"].fillna(cand["leakage_risk_level_overlap"])

    v2 = read_tsv(CANDIDATE_V2)
    if not v2.empty:
        v2_cols = [
            "candidate_id",
            "bma_v2_claim_safe_score",
            "bma_v2_discovery_score",
            "bma_v2_action",
            "primary_claim_blocker",
            "validity_dag_cap",
            "overlap_clean_cap",
            "patient_context_gate",
            "tcr_or_unknown_gate",
            "md_structure_gate",
        ]
        cand = cand.merge(v2[[c for c in v2_cols if c in v2.columns]], on="candidate_id", how="left")

    cand["label_binary"] = numeric(cand, "label", np.nan)
    leakage = text(cand, "leakage_risk_level").str.lower()
    cand["clean_low_leakage_trainable"] = (
        cand["label_binary"].isin([0, 1])
        & leakage.ne("high")
        & ~truthy(cand, "exact_peptide_train_overlap")
        & ~truthy(cand, "exact_peptide_hla_train_overlap")
        & text(cand, "mhc_class").eq("I")
    )
    cand["strict_no_near_overlap"] = cand["clean_low_leakage_trainable"] & ~truthy(cand, "near_peptide_train_overlap") & ~truthy(cand, "source_protein_window_train_overlap")
    return cand


def load_feature_matrix(clean_ids: set[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    ms = read_tsv(METHOD_SCORES, required=True)
    ms = ms.copy()
    ms["score_calibrated"] = pd.to_numeric(ms["score_calibrated"], errors="coerce")
    ms = ms[np.isfinite(ms["score_calibrated"])]

    method_meta = (
        ms.sort_values(["method_name"])
        .groupby("method_name", dropna=False)
        .agg(
            method_family=("method_family", "first"),
            method_role=("method_role", "first"),
            uses_public_pretraining=("uses_public_pretraining", "max"),
            clean_comparator_allowed=("clean_comparator_allowed", "min"),
            n_scored_rows=("candidate_id", "nunique"),
        )
        .reset_index()
    )

    pair_scores = (
        ms.groupby(["candidate_id", "method_name"], dropna=False)
        .agg(
            method_score_median=("score_calibrated", "median"),
            n_score_contexts=("score_context", "nunique"),
        )
        .reset_index()
    )
    matrix = pair_scores.pivot(index="candidate_id", columns="method_name", values="method_score_median")

    if clean_ids:
        coverage = matrix.loc[matrix.index.intersection(clean_ids)].notna().mean(axis=0)
    else:
        coverage = matrix.notna().mean(axis=0)
    method_meta = method_meta.merge(
        coverage.rename("clean_train_coverage").reset_index().rename(columns={"index": "method_name"}),
        on="method_name",
        how="left",
    )
    return matrix, method_meta


def build_feature_sets(matrix: pd.DataFrame, method_meta: pd.DataFrame, method_v2: pd.DataFrame) -> dict[str, list[str]]:
    meta = method_meta.copy()
    if not method_v2.empty:
        add_cols = [
            "method_name",
            "bma_v2_weight",
            "bma_v2_utility_norm",
            "bma_v2_use_class",
            "stress_min_floor",
            "stress_delta_floor",
            "public_or_clean_comparator_penalty",
            "fallback_role_cap",
        ]
        meta = meta.merge(method_v2[[c for c in add_cols if c in method_v2.columns]], on="method_name", how="left")

    meta["is_clean_internal"] = ~truthy(meta, "uses_public_pretraining") & truthy(meta, "clean_comparator_allowed", True)
    meta["available"] = meta["method_name"].isin(matrix.columns)
    meta["coverage_ok"] = numeric(meta, "clean_train_coverage", 0.0) >= 0.65
    meta = meta[meta["available"] & meta["coverage_ok"]]
    clean = meta[meta["is_clean_internal"]].copy()
    clean = clean.sort_values(["bma_v2_weight", "clean_train_coverage"], ascending=False)

    feature_sets: dict[str, list[str]] = {}
    feature_sets["all_clean_internal"] = clean["method_name"].tolist()
    feature_sets["top20_bma_v2_clean"] = clean.head(20)["method_name"].tolist()
    feature_sets["top40_bma_v2_clean"] = clean.head(40)["method_name"].tolist()

    stress_threshold = numeric(clean, "stress_min_floor", 0.0).median() if len(clean) else 0
    feature_sets["stress_floor_above_median"] = clean[numeric(clean, "stress_min_floor", 0.0) >= stress_threshold]["method_name"].tolist()

    families = text(clean, "method_family").str.lower()
    roles = text(clean, "method_role").str.lower()
    feature_sets["classical_stack_mix"] = clean[
        families.str.contains("classical|biophysical|stack|ensemble", regex=True)
        | roles.str.contains("internal|anchor", regex=True)
    ]["method_name"].head(45).tolist()

    return {k: v for k, v in feature_sets.items() if len(v) >= 3}


def make_pipeline(c_value: float) -> Pipeline:
    return Pipeline(
        steps=[
            ("imputer", SimpleImputer(strategy="median")),
            ("scaler", StandardScaler()),
            (
                "model",
                LogisticRegression(
                    C=c_value,
                    class_weight="balanced",
                    penalty="l2",
                    solver="liblinear",
                    max_iter=1000,
                    random_state=13,
                ),
            ),
        ]
    )


def get_splits(y: np.ndarray, groups: np.ndarray, protocol: str) -> list[tuple[np.ndarray, np.ndarray]]:
    y = np.asarray(y).astype(int)
    groups = np.asarray(groups).astype(str)
    n_groups = len(pd.unique(groups))
    if n_groups < 2:
        return []

    if protocol == "source_heldout_oof":
        splitter = GroupKFold(n_splits=min(n_groups, 5))
        return list(splitter.split(np.zeros_like(y), y, groups))

    n_splits = min(5, n_groups)
    splitter = StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=13)
    return list(splitter.split(np.zeros_like(y), y, groups))


def run_oof(
    x: pd.DataFrame,
    y: np.ndarray,
    groups: np.ndarray,
    protocol: str,
    c_value: float,
) -> tuple[np.ndarray, int]:
    oof = np.full(len(y), np.nan)
    splits = get_splits(y, groups, protocol)
    used = 0
    for train_idx, test_idx in splits:
        if len(np.unique(y[train_idx])) < 2:
            continue
        model = make_pipeline(c_value)
        model.fit(x.iloc[train_idx], y[train_idx])
        oof[test_idx] = model.predict_proba(x.iloc[test_idx])[:, 1]
        used += 1
    return oof, used


def evaluate_configs(
    matrix: pd.DataFrame,
    meta: pd.DataFrame,
    feature_sets: dict[str, list[str]],
) -> tuple[pd.DataFrame, dict[str, np.ndarray]]:
    train = meta[meta["clean_low_leakage_trainable"] & meta["candidate_id"].isin(matrix.index)].copy()
    train = train.sort_values("candidate_id")
    y = train["label_binary"].astype(int).to_numpy()

    c_values = [0.03, 0.10, 0.30, 1.00, 3.00]
    rows: list[ConfigResult] = []
    oof_store: dict[str, np.ndarray] = {}

    for feature_set, features in feature_sets.items():
        x = matrix.reindex(train["candidate_id"])[features]
        for c_value in c_values:
            for protocol, group_col in [
                ("source_heldout_oof", "source_name"),
                ("hla_heldout_oof", "hla_allele_4digit"),
            ]:
                groups = text(train, group_col, "unknown").to_numpy()
                oof, n_folds = run_oof(x, y, groups, protocol, c_value)
                metrics = metric_row(y, oof)
                key = f"{feature_set}__C{c_value:g}__{protocol}"
                oof_store[key] = oof
                rows.append(
                    ConfigResult(
                        feature_set=feature_set,
                        c_value=c_value,
                        protocol=protocol,
                        n_features=len(features),
                        n_rows=metrics["n_rows"],
                        n_pos=metrics["n_pos"],
                        auprc=metrics["auprc"],
                        auroc=metrics["auroc"],
                        brier=metrics["brier"],
                        ece=metrics["ece"],
                        top10_precision=metrics["top10_precision"],
                        top20_precision=metrics["top20_precision"],
                        oof_available=int(np.isfinite(oof).sum()),
                    )
                )

    metrics_df = pd.DataFrame([r.__dict__ for r in rows])
    return metrics_df, oof_store


def evaluate_reference_scores(meta: pd.DataFrame) -> pd.DataFrame:
    train = meta[meta["clean_low_leakage_trainable"]].copy()
    y = train["label_binary"].astype(int).to_numpy()
    refs = [
        "stress_guarded_final_review_score",
        "stress_guarded_clean_score",
        "stress_guarded_discovery_score",
        "bma_v2_discovery_score",
        "bma_v2_claim_safe_score",
        "validity_dag_cap",
    ]
    rows = []
    for col in refs:
        if col not in train.columns:
            continue
        m = metric_row(y, numeric(train, col, np.nan).to_numpy())
        row = {
            "feature_set": f"reference_{col}",
            "c_value": np.nan,
            "protocol": "clean_subset_reference_not_refit",
            "n_features": 1,
            "oof_available": int(np.isfinite(numeric(train, col, np.nan).to_numpy()).sum()),
        }
        row.update(m)
        rows.append(row)
    return pd.DataFrame(rows)


def select_best_config(metrics: pd.DataFrame) -> dict:
    if metrics.empty:
        raise ValueError("No finetune metrics available.")
    m = metrics.copy()
    m = m[np.isfinite(m["auprc"]) & m["protocol"].isin(["source_heldout_oof", "hla_heldout_oof"])]
    if m.empty:
        raise ValueError("No valid OOF metrics available.")
    grouped = (
        m.groupby(["feature_set", "c_value"], dropna=False)
        .agg(
            mean_auprc=("auprc", "mean"),
            min_auprc=("auprc", "min"),
            mean_auroc=("auroc", "mean"),
            mean_ece=("ece", "mean"),
            mean_top10_precision=("top10_precision", "mean"),
            n_protocols=("protocol", "nunique"),
            n_features=("n_features", "first"),
        )
        .reset_index()
    )
    grouped["selection_score"] = (
        grouped["mean_auprc"]
        + 0.15 * grouped["min_auprc"]
        + 0.05 * grouped["mean_top10_precision"]
        - 0.20 * grouped["mean_ece"]
        - 0.0008 * grouped["n_features"]
    )
    grouped = grouped.sort_values(["selection_score", "min_auprc"], ascending=False)
    return grouped.iloc[0].to_dict()


def fit_final_and_score(
    matrix: pd.DataFrame,
    meta: pd.DataFrame,
    features: list[str],
    c_value: float,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    train = meta[meta["clean_low_leakage_trainable"] & meta["candidate_id"].isin(matrix.index)].copy()
    y = train["label_binary"].astype(int).to_numpy()
    x_train = matrix.reindex(train["candidate_id"])[features]
    model = make_pipeline(c_value)
    model.fit(x_train, y)

    all_ids = pd.Index(sorted(set(meta["candidate_id"]).intersection(matrix.index)))
    all_x = matrix.reindex(all_ids)[features]
    scores = model.predict_proba(all_x)[:, 1]
    scored = pd.DataFrame({"candidate_id": all_ids, "finetuned_clean_stack_score": scores})
    scored = scored.merge(meta, on="candidate_id", how="left")
    scored["finetuned_rank_global"] = scored["finetuned_clean_stack_score"].rank(ascending=False, method="first").astype(int)
    scored["finetune_vs_bma_discovery_delta"] = scored["finetuned_clean_stack_score"] - numeric(scored, "bma_v2_discovery_score", 0.0)
    scored["finetune_action"] = np.select(
        [
            (scored["finetuned_clean_stack_score"] >= 0.70) & (numeric(scored, "overlap_clean_cap", 0.0) <= 0.42),
            (scored["finetuned_clean_stack_score"] >= 0.70) & (numeric(scored, "patient_context_gate", 0.0) <= 0.30),
            (scored["finetuned_clean_stack_score"] >= 0.70),
            (scored["finetuned_clean_stack_score"] <= 0.20) & (numeric(scored, "bma_v2_discovery_score", 0.0) >= 0.55),
            scored["finetune_vs_bma_discovery_delta"] >= 0.20,
        ],
        [
            "high_finetune_score_but_overlap_blocked",
            "high_finetune_score_but_patient_context_blocked",
            "high_finetune_score_manual_candidate_audit",
            "demote_bma_high_but_clean_finetune_low",
            "finetune_rescue_review",
        ],
        default="score_support_only",
    )

    coef = model.named_steps["model"].coef_[0]
    coef_df = pd.DataFrame({"method_name": features, "coefficient": coef})
    coef_df["abs_coefficient"] = coef_df["coefficient"].abs()
    coef_df = coef_df.sort_values("abs_coefficient", ascending=False)
    return scored.sort_values("finetuned_clean_stack_score", ascending=False), coef_df


def method_underperformance_actions(
    method_meta: pd.DataFrame,
    method_v2: pd.DataFrame,
    split_metrics: pd.DataFrame,
    selected_features: Iterable[str],
) -> pd.DataFrame:
    meta = method_meta.copy()
    if not method_v2.empty:
        cols = [
            "method_name",
            "bma_v2_weight",
            "bma_v2_utility_norm",
            "bma_v2_use_class",
            "stress_min_floor",
            "stress_delta_floor",
            "public_or_clean_comparator_penalty",
            "fallback_role_cap",
        ]
        meta = meta.merge(method_v2[[c for c in cols if c in method_v2.columns]], on="method_name", how="left")

    if not split_metrics.empty:
        sm = split_metrics.copy()
        sm["axis"] = "other"
        raw = (text(sm, "split_contract") + " " + text(sm, "split_group")).str.lower()
        sm.loc[raw.str.contains("source"), "axis"] = "source"
        sm.loc[raw.str.contains("hla"), "axis"] = "hla"
        sm.loc[raw.str.contains("low_prevalence|low prevalence", regex=True), "axis"] = "low_prevalence"
        agg = sm.groupby("method_name").agg(
            split_rows=("split_contract", "size"),
            mean_AUPRC=("AUPRC", "mean"),
            min_AUPRC=("AUPRC", "min"),
            source_min_AUPRC=("AUPRC", lambda s: float(s[sm.loc[s.index, "axis"].eq("source")].min()) if sm.loc[s.index, "axis"].eq("source").any() else np.nan),
            hla_min_AUPRC=("AUPRC", lambda s: float(s[sm.loc[s.index, "axis"].eq("hla")].min()) if sm.loc[s.index, "axis"].eq("hla").any() else np.nan),
            low_prevalence_AUPRC=("AUPRC", lambda s: float(s[sm.loc[s.index, "axis"].eq("low_prevalence")].median()) if sm.loc[s.index, "axis"].eq("low_prevalence").any() else np.nan),
            mean_ECE=("calibration_ece", "mean"),
            mean_top10_precision=("top10_precision", "mean"),
        ).reset_index()
        meta = meta.merge(agg, on="method_name", how="left")

    selected = set(selected_features)
    meta["selected_in_finetuned_stack"] = meta["method_name"].isin(selected)
    meta["finetune_component_action"] = np.select(
        [
            truthy(meta, "uses_public_pretraining") | ~truthy(meta, "clean_comparator_allowed", True),
            numeric(meta, "clean_train_coverage", 0.0) < 0.65,
            numeric(meta, "source_min_AUPRC", 1.0) < 0.25,
            numeric(meta, "hla_min_AUPRC", 1.0) < 0.15,
            numeric(meta, "mean_ECE", 0.0) > 0.35,
            meta["selected_in_finetuned_stack"],
        ],
        [
            "demote_to_ablation_or_support_only_public_overlap_risk",
            "drop_from_finetune_missing_clean_coverage",
            "source_shift_cap_or_retrain_needed",
            "hla_shift_cap_or_allele_reweight_needed",
            "calibration_repair_needed_before_claim_use",
            "kept_in_clean_finetuned_stack",
        ],
        default="reserve_support_only",
    )
    meta["underperformance_severity"] = np.select(
        [
            meta["finetune_component_action"].str.contains("demote|drop", regex=True),
            meta["finetune_component_action"].str.contains("cap|repair", regex=True),
            meta["selected_in_finetuned_stack"],
        ],
        ["high", "medium", "kept"],
        default="low",
    )
    front = [
        "method_name",
        "method_family",
        "method_role",
        "finetune_component_action",
        "underperformance_severity",
        "selected_in_finetuned_stack",
        "clean_train_coverage",
        "bma_v2_weight",
        "source_min_AUPRC",
        "hla_min_AUPRC",
        "low_prevalence_AUPRC",
        "mean_ECE",
        "mean_top10_precision",
        "uses_public_pretraining",
        "clean_comparator_allowed",
    ]
    front = [c for c in front if c in meta.columns]
    return meta.sort_values(["underperformance_severity", "bma_v2_weight"], ascending=[True, False])[front + [c for c in meta.columns if c not in front]]


def make_report(
    summary: dict,
    all_metrics: pd.DataFrame,
    selected: dict,
    scored: pd.DataFrame,
    coef: pd.DataFrame,
    under: pd.DataFrame,
) -> str:
    best_rows = all_metrics[
        all_metrics["feature_set"].eq(selected["feature_set"])
        & all_metrics["c_value"].eq(selected["c_value"])
    ][
        [
            "feature_set",
            "c_value",
            "protocol",
            "n_features",
            "n_rows",
            "n_pos",
            "auprc",
            "auroc",
            "brier",
            "ece",
            "top10_precision",
            "top20_precision",
        ]
    ]
    ref_rows = all_metrics[all_metrics["protocol"].eq("clean_subset_reference_not_refit")][
        ["feature_set", "protocol", "auprc", "auroc", "brier", "ece", "top10_precision", "top20_precision"]
    ].sort_values("auprc", ascending=False)
    top_scores = scored.head(15)[
        [
            "candidate_id",
            "peptide",
            "hla_allele_4digit",
            "label",
            "source_name",
            "leakage_risk_level",
            "finetuned_clean_stack_score",
            "bma_v2_discovery_score",
            "bma_v2_claim_safe_score",
            "finetune_action",
            "primary_claim_blocker",
        ]
    ]
    top_coef = coef.head(15)[["method_name", "coefficient", "abs_coefficient"]]
    weak = under[under["underperformance_severity"].isin(["high", "medium"])].head(15)
    weak_cols = [c for c in ["method_name", "finetune_component_action", "underperformance_severity", "clean_train_coverage", "source_min_AUPRC", "hla_min_AUPRC", "mean_ECE"] if c in weak.columns]

    def table(df: pd.DataFrame) -> str:
        if df.empty:
            return "_no rows_"
        return df.to_markdown(index=False, floatfmt=".3f")

    return f"""# CROSS-Neo clean-CV scoring + finetune report

Generated: {summary["generated_at"]}

## Result

Clean low-leakage subset에서 source/HLA held-out OOF로 lightweight stacker를 튜닝했다. Best config는 `{selected["feature_set"]}` / `C={selected["c_value"]}`이고, 선택 기준은 mean AUPRC, worst-axis AUPRC, top-10 precision, calibration penalty를 같이 본 것이다.

## Best fine-tuned CV score

{table(best_rows)}

## Reference scores on same clean subset

These are not refit OOF models; they are context/reference scores on the same clean rows.

{table(ref_rows)}

## Top rescored candidates

{table(top_scores)}

## Top positive coefficients

{table(top_coef)}

## Underperforming components to cap/drop/recalibrate

{table(weak[weak_cols])}

## Decision

1. 점수는 `finetuned_clean_stack_score`로 냈다. 이 점수는 clean-CV tuned score이고, overlap/patient gate를 통과하지 못한 후보는 여전히 claim-safe가 아니다.
2. 너무 안 나오는 축은 새 모델을 바로 키우기보다 `underperforming_component_finetune_actions.tsv`의 cap/drop/recalibrate 지시대로 줄인다.
3. P0 wetlab 후보는 fine-tune score가 높아도 `overlap_clean_cap`과 `patient_context_gate`가 잠긴 상태라 실험 우선순위와 논문 claim을 분리한다.

## Files

- `clean_cv_finetune_metrics.tsv`
- `fine_tuned_candidate_scores.tsv`
- `fine_tuned_stack_coefficients.tsv`
- `underperforming_component_finetune_actions.tsv`
- `clean_cv_finetune_summary.json`
"""


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    meta = load_metadata()
    clean_ids = set(meta.loc[meta["clean_low_leakage_trainable"], "candidate_id"].astype(str))
    matrix, method_meta = load_feature_matrix(clean_ids)
    method_v2 = read_tsv(METHOD_V2)
    split_metrics = read_tsv(SPLIT_METRICS)
    feature_sets = build_feature_sets(matrix, method_meta, method_v2)

    metrics, _ = evaluate_configs(matrix, meta, feature_sets)
    refs = evaluate_reference_scores(meta)
    all_metrics = pd.concat([metrics, refs], ignore_index=True, sort=False)
    selected = select_best_config(metrics)
    selected_features = feature_sets[selected["feature_set"]]

    scored, coef = fit_final_and_score(matrix, meta, selected_features, float(selected["c_value"]))
    under = method_underperformance_actions(method_meta, method_v2, split_metrics, selected_features)

    write_tsv(all_metrics.sort_values(["protocol", "auprc"], ascending=[True, False]), "clean_cv_finetune_metrics.tsv")
    write_tsv(scored, "fine_tuned_candidate_scores.tsv")
    write_tsv(coef, "fine_tuned_stack_coefficients.tsv")
    write_tsv(under, "underperforming_component_finetune_actions.tsv")

    summary = {
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "output_dir": str(OUT_DIR.relative_to(ROOT)),
        "n_candidates_metadata": int(len(meta)),
        "n_candidates_with_method_scores": int(matrix.shape[0]),
        "n_clean_trainable": int(meta["clean_low_leakage_trainable"].sum()),
        "n_clean_positive": int(meta.loc[meta["clean_low_leakage_trainable"], "label_binary"].sum()),
        "n_feature_sets": int(len(feature_sets)),
        "selected_config": selected,
        "selected_feature_count": int(len(selected_features)),
        "top_scored_candidates": scored.head(10)[["candidate_id", "peptide", "hla_allele_4digit", "finetuned_clean_stack_score", "finetune_action"]].to_dict("records"),
        "output_files": {
            "clean_cv_finetune_metrics": "clean_cv_finetune_metrics.tsv",
            "fine_tuned_candidate_scores": "fine_tuned_candidate_scores.tsv",
            "fine_tuned_stack_coefficients": "fine_tuned_stack_coefficients.tsv",
            "underperforming_component_finetune_actions": "underperforming_component_finetune_actions.tsv",
            "report": "CLEAN_CV_FINETUNE_REPORT_KR.md",
        },
    }
    (OUT_DIR / "clean_cv_finetune_summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n")
    report = make_report(summary, all_metrics, selected, scored, coef, under)
    (OUT_DIR / "CLEAN_CV_FINETUNE_REPORT_KR.md").write_text(report)
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
