#!/usr/bin/env python3
"""Candidate-level score booster for CROSS-Neo.

This follows the failed method-matrix fine-tune with a safer fallback:
fine-tune/calibrate the already strong BAR-Neo stress score family under
source-heldout and HLA-heldout clean CV, then score all candidates with separate
experiment-priority and claim-capped scores.
"""

from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import ExtraTreesClassifier, GradientBoostingClassifier, RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score
from sklearn.model_selection import GroupKFold, StratifiedGroupKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


ROOT = Path(__file__).resolve().parents[1]
BAR_DIR = ROOT / "project/results/clean_neobench_barneo_2026_05_09"
P0_DIR = ROOT / "project/results/cross_neo_kaggle_winner_playbook_2026_05_10/merged_p0_experiments"
OLD_FINETUNE_DIR = ROOT / "project/results/cross_neo_kaggle_winner_playbook_2026_05_10/score_finetune_2026_05_10"
OUT_DIR = ROOT / "project/results/cross_neo_kaggle_winner_playbook_2026_05_10/score_booster_finetune_2026_05_10"

CANDIDATES = BAR_DIR / "barneo_stress_guarded_candidate_scores.tsv"
OVERLAP = BAR_DIR / "clean_neobench_overlap_flags.tsv"
CANDIDATE_V2 = P0_DIR / "barneo_bma_v2_diversity_validity_scores.tsv"
WETLAB_V2 = P0_DIR / "cross_neo_wetlab_plate_v2.tsv"
MATRIX_METRICS = OLD_FINETUNE_DIR / "clean_cv_finetune_metrics.tsv"


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


def ece_score(y_true: np.ndarray, y_prob: np.ndarray, bins: int = 10) -> float:
    y_true = np.asarray(y_true).astype(float)
    y_prob = np.asarray(y_prob).astype(float)
    mask = np.isfinite(y_true) & np.isfinite(y_prob)
    y_true = y_true[mask]
    y_prob = np.clip(y_prob[mask], 0, 1)
    if len(y_true) == 0:
        return float("nan")
    edges = np.linspace(0, 1, bins + 1)
    out = 0.0
    for lo, hi in zip(edges[:-1], edges[1:]):
        in_bin = (y_prob >= lo) & (y_prob < hi if hi < 1 else y_prob <= hi)
        if in_bin.any():
            out += in_bin.mean() * abs(y_true[in_bin].mean() - y_prob[in_bin].mean())
    return float(out)


def topk_precision(y_true: np.ndarray, y_prob: np.ndarray, k: int) -> float:
    y_true = np.asarray(y_true).astype(float)
    y_prob = np.asarray(y_prob).astype(float)
    mask = np.isfinite(y_true) & np.isfinite(y_prob)
    y_true = y_true[mask]
    y_prob = y_prob[mask]
    if len(y_true) == 0:
        return float("nan")
    order = np.argsort(-y_prob)[: min(k, len(y_true))]
    return float(y_true[order].mean())


def metrics(y: np.ndarray, p: np.ndarray) -> dict:
    y = np.asarray(y).astype(int)
    p = np.asarray(p).astype(float)
    mask = np.isfinite(p)
    y = y[mask]
    p = p[mask]
    out = {
        "n_rows": int(len(y)),
        "n_pos": int(y.sum()) if len(y) else 0,
        "auprc": float("nan"),
        "auroc": float("nan"),
        "brier": float("nan"),
        "ece": float("nan"),
        "top10_precision": float("nan"),
        "top20_precision": float("nan"),
    }
    if len(y) == 0 or len(np.unique(y)) < 2:
        return out
    out["auprc"] = float(average_precision_score(y, p))
    out["auroc"] = float(roc_auc_score(y, p))
    out["brier"] = float(brier_score_loss(y, np.clip(p, 0, 1)))
    out["ece"] = ece_score(y, p)
    out["top10_precision"] = topk_precision(y, p, 10)
    out["top20_precision"] = topk_precision(y, p, 20)
    return out


def load_candidates() -> pd.DataFrame:
    c = read_tsv(CANDIDATES, required=True)
    overlap = read_tsv(OVERLAP)
    if not overlap.empty:
        c = c.merge(overlap, on="candidate_id", how="left", suffixes=("", "_overlap"))
        if "leakage_risk_level_overlap" in c.columns:
            c["leakage_risk_level"] = c["leakage_risk_level"].fillna(c["leakage_risk_level_overlap"])
    v2 = read_tsv(CANDIDATE_V2, required=True)
    v2_keep = [
        "candidate_id",
        "bma_v2_claim_safe_score",
        "bma_v2_discovery_score",
        "bma_v2_action",
        "primary_claim_blocker",
        "validity_dag_cap",
        "antigen_model_gate",
        "presentation_hla_gate",
        "tcr_or_unknown_gate",
        "md_structure_gate",
        "patient_context_gate",
        "overlap_clean_cap",
        "source_hla_generalization_cap",
    ]
    c = c.merge(v2[[col for col in v2_keep if col in v2.columns]], on="candidate_id", how="left")
    leakage = text(c, "leakage_risk_level").str.lower()
    c["clean_low_leakage_trainable"] = (
        numeric(c, "label", np.nan).isin([0, 1])
        & leakage.ne("high")
        & ~truthy(c, "exact_peptide_train_overlap")
        & ~truthy(c, "exact_peptide_hla_train_overlap")
        & text(c, "mhc_class").eq("I")
    )
    return c


def model_factory(model_name: str) -> Pipeline:
    if model_name.startswith("logistic_C"):
        c_value = float(model_name.replace("logistic_C", ""))
        return Pipeline(
            [
                ("imputer", SimpleImputer(strategy="median")),
                ("scaler", StandardScaler()),
                ("model", LogisticRegression(C=c_value, class_weight="balanced", solver="liblinear", max_iter=1000, random_state=13)),
            ]
        )
    if model_name == "rf_depth3":
        return Pipeline(
            [
                ("imputer", SimpleImputer(strategy="median")),
                (
                    "model",
                    RandomForestClassifier(
                        n_estimators=300,
                        max_depth=3,
                        min_samples_leaf=8,
                        class_weight="balanced_subsample",
                        random_state=13,
                        n_jobs=-1,
                    ),
                ),
            ]
        )
    if model_name == "extra_depth3":
        return Pipeline(
            [
                ("imputer", SimpleImputer(strategy="median")),
                (
                    "model",
                    ExtraTreesClassifier(
                        n_estimators=300,
                        max_depth=3,
                        min_samples_leaf=8,
                        class_weight="balanced",
                        random_state=13,
                        n_jobs=-1,
                    ),
                ),
            ]
        )
    if model_name == "gb_depth2_lr003":
        return Pipeline(
            [
                ("imputer", SimpleImputer(strategy="median")),
                ("model", GradientBoostingClassifier(n_estimators=80, learning_rate=0.03, max_depth=2, random_state=13)),
            ]
        )
    raise ValueError(model_name)


def feature_sets(c: pd.DataFrame) -> dict[str, list[str]]:
    sets = {
        "stress_only": ["stress_guarded_final_review_score"],
        "stress_plus_components": [
            "stress_guarded_final_review_score",
            "stress_guarded_clean_score",
            "stress_guarded_discovery_score",
            "stress_guarded_confidence",
            "stress_guarded_expert_disagreement",
        ],
        "bma_v2_gates": [
            "bma_v2_discovery_score",
            "bma_v2_claim_safe_score",
            "validity_dag_cap",
            "antigen_model_gate",
            "presentation_hla_gate",
            "tcr_or_unknown_gate",
            "md_structure_gate",
            "patient_context_gate",
            "overlap_clean_cap",
            "source_hla_generalization_cap",
        ],
        "stress_plus_bma_v2": [
            "stress_guarded_final_review_score",
            "stress_guarded_clean_score",
            "stress_guarded_discovery_score",
            "stress_guarded_confidence",
            "stress_guarded_expert_disagreement",
            "bma_v2_discovery_score",
            "bma_v2_claim_safe_score",
            "validity_dag_cap",
            "antigen_model_gate",
            "presentation_hla_gate",
            "tcr_or_unknown_gate",
            "md_structure_gate",
            "patient_context_gate",
            "overlap_clean_cap",
            "source_hla_generalization_cap",
        ],
    }
    return {name: [col for col in cols if col in c.columns] for name, cols in sets.items()}


def split_indices(y: np.ndarray, groups: np.ndarray, protocol: str) -> list[tuple[np.ndarray, np.ndarray]]:
    unique_groups = len(pd.unique(groups))
    if unique_groups < 2:
        return []
    if protocol == "source_heldout_oof":
        return list(GroupKFold(n_splits=min(5, unique_groups)).split(np.zeros_like(y), y, groups))
    return list(
        StratifiedGroupKFold(n_splits=min(5, unique_groups), shuffle=True, random_state=13).split(
            np.zeros_like(y), y, groups
        )
    )


def run_cv(c: pd.DataFrame, fsets: dict[str, list[str]]) -> pd.DataFrame:
    train = c[c["clean_low_leakage_trainable"]].copy().sort_values("candidate_id")
    y = numeric(train, "label", np.nan).astype(int).to_numpy()
    rows = []
    model_names = ["logistic_C0.03", "logistic_C0.1", "logistic_C0.3", "logistic_C1.0", "rf_depth3", "extra_depth3", "gb_depth2_lr003"]

    for fs_name, cols in fsets.items():
        x = train[cols]
        for model_name in model_names:
            for protocol, group_col in [("source_heldout_oof", "source_name"), ("hla_heldout_oof", "hla_allele_4digit")]:
                groups = text(train, group_col, "unknown").to_numpy()
                oof = np.full(len(train), np.nan)
                folds = 0
                for tr, te in split_indices(y, groups, protocol):
                    if len(np.unique(y[tr])) < 2:
                        continue
                    model = model_factory(model_name)
                    model.fit(x.iloc[tr], y[tr])
                    oof[te] = model.predict_proba(x.iloc[te])[:, 1]
                    folds += 1
                row = {
                    "feature_set": fs_name,
                    "model_name": model_name,
                    "protocol": protocol,
                    "n_features": len(cols),
                    "folds_used": folds,
                    "oof_available": int(np.isfinite(oof).sum()),
                }
                row.update(metrics(y, oof))
                rows.append(row)

    for score_col in [
        "stress_guarded_final_review_score",
        "stress_guarded_discovery_score",
        "stress_guarded_clean_score",
        "bma_v2_discovery_score",
        "bma_v2_claim_safe_score",
    ]:
        row = {
            "feature_set": f"reference_{score_col}",
            "model_name": "reference_not_refit",
            "protocol": "clean_subset_reference_not_refit",
            "n_features": 1,
            "folds_used": 0,
            "oof_available": int(numeric(train, score_col, np.nan).notna().sum()),
        }
        row.update(metrics(y, numeric(train, score_col, np.nan).to_numpy()))
        rows.append(row)
    return pd.DataFrame(rows)


def select_config(cv: pd.DataFrame) -> dict:
    oof = cv[cv["protocol"].isin(["source_heldout_oof", "hla_heldout_oof"]) & np.isfinite(cv["auprc"])].copy()
    grouped = (
        oof.groupby(["feature_set", "model_name"], dropna=False)
        .agg(
            mean_auprc=("auprc", "mean"),
            min_auprc=("auprc", "min"),
            mean_auroc=("auroc", "mean"),
            mean_brier=("brier", "mean"),
            mean_ece=("ece", "mean"),
            mean_top10_precision=("top10_precision", "mean"),
            mean_top20_precision=("top20_precision", "mean"),
            n_features=("n_features", "first"),
            protocols=("protocol", "nunique"),
        )
        .reset_index()
    )
    grouped["selection_score"] = (
        grouped["mean_auprc"]
        + 0.15 * grouped["min_auprc"]
        + 0.05 * grouped["mean_top10_precision"]
        - 0.15 * grouped["mean_ece"]
        - 0.001 * grouped["n_features"]
    )
    grouped = grouped.sort_values(["selection_score", "min_auprc"], ascending=False)
    return grouped.iloc[0].to_dict()


def fit_score_all(c: pd.DataFrame, fsets: dict[str, list[str]], selected: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    train = c[c["clean_low_leakage_trainable"]].copy().sort_values("candidate_id")
    cols = fsets[selected["feature_set"]]
    y = numeric(train, "label", np.nan).astype(int).to_numpy()
    model = model_factory(selected["model_name"])
    model.fit(train[cols], y)

    scored = c.copy()
    scored["finetuned_score_booster_prob"] = model.predict_proba(scored[cols])[:, 1]
    scored["finetuned_experiment_priority_score"] = np.clip(
        0.70 * scored["finetuned_score_booster_prob"] + 0.30 * numeric(scored, "bma_v2_discovery_score", 0.0),
        0,
        1,
    )
    cap_cols = [
        "validity_dag_cap",
        "overlap_clean_cap",
        "patient_context_gate",
        "source_hla_generalization_cap",
    ]
    available_caps = [col for col in cap_cols if col in scored.columns]
    scored["finetuned_claim_capped_score"] = np.minimum(
        scored["finetuned_score_booster_prob"],
        scored[available_caps].astype(float).min(axis=1),
    )
    scored["finetuned_score_action"] = np.select(
        [
            (scored["finetuned_score_booster_prob"] >= 0.70) & (numeric(scored, "overlap_clean_cap", 1.0) <= 0.42),
            (scored["finetuned_score_booster_prob"] >= 0.70) & (numeric(scored, "patient_context_gate", 1.0) <= 0.30),
            (scored["finetuned_claim_capped_score"] >= 0.55),
            (scored["finetuned_experiment_priority_score"] >= 0.65),
            (scored["finetuned_score_booster_prob"] <= 0.25) & (numeric(scored, "bma_v2_discovery_score", 0.0) >= 0.55),
        ],
        [
            "score_high_but_overlap_blocked",
            "score_high_but_patient_context_blocked",
            "claim_score_candidate_after_manual_audit",
            "experiment_priority_with_claim_boundary",
            "demote_previous_high_score",
        ],
        default="score_support_only",
    )

    scored["finetuned_prob_rank"] = scored["finetuned_score_booster_prob"].rank(ascending=False, method="first").astype(int)
    scored["finetuned_experiment_priority_rank"] = scored["finetuned_experiment_priority_score"].rank(ascending=False, method="first").astype(int)
    scored["finetuned_claim_capped_rank"] = scored["finetuned_claim_capped_score"].rank(ascending=False, method="first").astype(int)

    fitted = model.named_steps["model"]
    if hasattr(fitted, "coef_"):
        imp = pd.DataFrame({"feature": cols, "importance": fitted.coef_[0]})
    elif hasattr(fitted, "feature_importances_"):
        imp = pd.DataFrame({"feature": cols, "importance": fitted.feature_importances_})
    else:
        imp = pd.DataFrame({"feature": cols, "importance": np.nan})
    imp["abs_importance"] = imp["importance"].abs()
    imp = imp.sort_values("abs_importance", ascending=False)
    return scored.sort_values("finetuned_experiment_priority_score", ascending=False), imp


def make_report(summary: dict, cv: pd.DataFrame, selected: dict, scored: pd.DataFrame, imp: pd.DataFrame) -> str:
    best_rows = cv[
        cv["feature_set"].eq(selected["feature_set"])
        & cv["model_name"].eq(selected["model_name"])
    ][
        [
            "feature_set",
            "model_name",
            "protocol",
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
    refs = cv[cv["protocol"].eq("clean_subset_reference_not_refit")][
        ["feature_set", "auprc", "auroc", "brier", "ece", "top10_precision", "top20_precision"]
    ].sort_values("auprc", ascending=False)
    top_exp = scored.head(12)[
        [
            "candidate_id",
            "peptide",
            "hla_allele_4digit",
            "label",
            "source_name",
            "leakage_risk_level",
            "finetuned_score_booster_prob",
            "finetuned_experiment_priority_score",
            "finetuned_claim_capped_score",
            "finetuned_score_action",
            "primary_claim_blocker",
        ]
    ]
    top_claim = scored.sort_values("finetuned_claim_capped_score", ascending=False).head(12)[
        [
            "candidate_id",
            "peptide",
            "hla_allele_4digit",
            "label",
            "source_name",
            "leakage_risk_level",
            "finetuned_score_booster_prob",
            "finetuned_claim_capped_score",
            "finetuned_score_action",
            "primary_claim_blocker",
        ]
    ]

    def table(df: pd.DataFrame) -> str:
        if df.empty:
            return "_no rows_"
        return df.to_markdown(index=False, floatfmt=".3f")

    return f"""# CROSS-Neo score booster fine-tune report

Generated: {summary["generated_at"]}

## Result

Method-matrix stacker가 source-heldout에서 낮게 나와서, 기존 BAR-Neo stress score 계열을 clean source/HLA OOF로 calibration/booster fine-tune했다. Best config는 `{selected["feature_set"]}` + `{selected["model_name"]}`.

## Best clean-CV score

{table(best_rows)}

## Reference scores

{table(refs)}

## Top experiment-priority scores

{table(top_exp)}

## Top claim-capped scores

{table(top_claim)}

## Feature importance

{table(imp.head(10))}

## Decision

1. Final score column: `finetuned_score_booster_prob`.
2. Experiment ranking column: `finetuned_experiment_priority_score`.
3. Claim-safe ranking column: `finetuned_claim_capped_score`, which keeps overlap/patient/source gates active.
4. Low-performing method-matrix fine-tune is not used as the main score; it remains a diagnostic table for cap/drop/recalibration.
"""


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    candidates = load_candidates()
    fsets = feature_sets(candidates)
    cv = run_cv(candidates, fsets)
    selected = select_config(cv)
    scored, imp = fit_score_all(candidates, fsets, selected)

    matrix_metrics = read_tsv(MATRIX_METRICS)
    if not matrix_metrics.empty:
        matrix_best = matrix_metrics[
            matrix_metrics["protocol"].isin(["source_heldout_oof", "hla_heldout_oof"])
        ]["auprc"].max()
    else:
        matrix_best = np.nan

    write_tsv(cv.sort_values(["protocol", "auprc"], ascending=[True, False]), "score_booster_cv_metrics.tsv")
    write_tsv(scored, "score_booster_candidate_scores.tsv")
    write_tsv(imp, "score_booster_feature_importance.tsv")

    summary = {
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "output_dir": str(OUT_DIR.relative_to(ROOT)),
        "n_candidates": int(len(candidates)),
        "n_clean_trainable": int(candidates["clean_low_leakage_trainable"].sum()),
        "n_clean_positive": int(candidates.loc[candidates["clean_low_leakage_trainable"], "label"].sum()),
        "selected_config": selected,
        "method_matrix_best_oof_auprc_previous_run": float(matrix_best) if np.isfinite(matrix_best) else None,
        "top_experiment_priority": scored.head(10)[
            [
                "candidate_id",
                "peptide",
                "hla_allele_4digit",
                "finetuned_score_booster_prob",
                "finetuned_experiment_priority_score",
                "finetuned_claim_capped_score",
                "finetuned_score_action",
            ]
        ].to_dict("records"),
        "action_counts": scored["finetuned_score_action"].value_counts().to_dict(),
    }
    (OUT_DIR / "score_booster_summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n")
    (OUT_DIR / "SCORE_BOOSTER_FINETUNE_REPORT_KR.md").write_text(make_report(summary, cv, selected, scored, imp))
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
