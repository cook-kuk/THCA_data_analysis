#!/usr/bin/env python3
"""Train BAR-Neo, a benchmark-adaptive reliability ranking layer."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from common import (
    GradientBoostingClassifier,
    LogisticRegression,
    RandomForestClassifier,
    SKLEARN_AVAILABLE,
    ensure_dir,
    sequence_features,
    update_manifest,
    write_tsv,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", default=".", help="Repository root")
    parser.add_argument("--output-root", required=True, help="Output directory")
    return parser.parse_args()


def source_failure_lookup(repo: Path) -> dict[str, dict[str, float]]:
    path = repo / "project/results/cross_neo_v1_lockdown/source_collapse_diagnostics.tsv"
    if not path.exists():
        return {}
    df = pd.read_csv(path, sep="\t")
    out = {}
    for source, g in df.groupby("heldout_study"):
        out[source] = {
            "historical_source_top10": float(pd.to_numeric(g["top10_precision"], errors="coerce").fillna(0).max()),
            "historical_source_auprc": float(pd.to_numeric(g["AUPRC"], errors="coerce").fillna(0).max()),
        }
    return out


def build_candidate_features(master: pd.DataFrame, scores: pd.DataFrame, repo: Path) -> pd.DataFrame:
    base = master.copy()
    seq_df = pd.DataFrame([sequence_features(p) for p in base["peptide"]], index=base.index)
    base = pd.concat([base, seq_df], axis=1)

    score_pivot = scores.pivot_table(index="candidate_id", columns="method_name", values="score_calibrated", aggfunc="max")
    score_pivot.columns = [f"method_score__{c}" for c in score_pivot.columns]
    rank_pivot = scores.pivot_table(index="candidate_id", columns="method_name", values="rank_global", aggfunc="min")
    rank_pivot.columns = [f"method_rank__{c}" for c in rank_pivot.columns]
    role_counts = scores.pivot_table(index="candidate_id", columns="method_role", values="method_name", aggfunc="count", fill_value=0)
    role_counts.columns = [f"n_role__{c}" for c in role_counts.columns]
    method_n = scores.groupby("candidate_id")["method_name"].nunique().rename("n_methods_available")
    method_std = scores.groupby("candidate_id")["score_calibrated"].std().rename("method_disagreement_std")
    method_range = (scores.groupby("candidate_id")["score_calibrated"].max() - scores.groupby("candidate_id")["score_calibrated"].min()).rename("method_disagreement_range")
    feat = base.set_index("candidate_id").join([score_pivot, rank_pivot, role_counts, method_n, method_std, method_range])

    internal = scores[scores["method_role"].isin(["anchor", "internal_candidate", "bounded_fallback"])]
    public = scores[scores["method_role"].eq("caveated_public_comparator")]
    anchor = scores[scores["method_role"].eq("anchor")]
    feat["best_clean_internal_score"] = internal.groupby("candidate_id")["score_calibrated"].max()
    feat["best_caveated_public_score"] = public.groupby("candidate_id")["score_calibrated"].max()
    feat["best_anchor_score"] = anchor.groupby("candidate_id")["score_calibrated"].max()
    feat["public_vs_internal_disagreement"] = (feat["best_caveated_public_score"] - feat["best_clean_internal_score"]).abs()
    feat["anchor_vs_candidate_disagreement"] = (feat["best_anchor_score"] - feat["best_clean_internal_score"]).abs()
    feat["only_caveated_public_available"] = (
        feat.get("n_role__caveated_public_comparator", 0).fillna(0).gt(0)
        & feat.get("n_role__anchor", 0).fillna(0).eq(0)
        & feat.get("n_role__internal_candidate", 0).fillna(0).eq(0)
    )

    allele_counts = base["hla_allele_4digit"].value_counts().to_dict()
    source_prev = base.groupby("source_name")["label"].mean(numeric_only=True).to_dict()
    source_pos = base.groupby("source_name")["label"].sum(numeric_only=True).to_dict()
    failure = source_failure_lookup(repo)
    feat["hla_allele_support_count"] = feat["hla_allele_4digit"].map(allele_counts).fillna(0)
    feat["underrepresented_hla_allele"] = feat["hla_allele_support_count"].lt(10)
    feat["source_positive_prevalence"] = feat["source_name"].map(source_prev).fillna(0)
    feat["source_positive_count"] = feat["source_name"].map(source_pos).fillna(0)
    feat["low_prevalence_source"] = feat["source_positive_prevalence"].lt(0.10) | feat["source_positive_count"].lt(10)
    feat["historical_source_top10"] = feat["source_name"].map(lambda s: failure.get(s, {}).get("historical_source_top10", np.nan))
    feat["historical_source_auprc"] = feat["source_name"].map(lambda s: failure.get(s, {}).get("historical_source_auprc", np.nan))
    feat["source_ood_high"] = feat["historical_source_top10"].fillna(1).le(0.0)
    for col in [
        "exact_peptide_train_overlap",
        "exact_peptide_hla_train_overlap",
        "near_peptide_train_overlap",
        "study_train_overlap",
        "patient_train_overlap",
        "public_tool_training_overlap_any",
    ]:
        feat[col] = feat[col].fillna(False).astype(bool)
    return feat.reset_index()


def choose_model() -> object | None:
    if not SKLEARN_AVAILABLE:
        return None
    # A shallow forest is fast, robust to unscaled sparse reliability features,
    # and avoids long optimizer convergence in marathon-mode reruns.
    return RandomForestClassifier(
        n_estimators=80,
        max_depth=6,
        min_samples_leaf=5,
        class_weight="balanced_subsample",
        random_state=20260509,
        n_jobs=-1,
    )


def prepare_matrix(feat: pd.DataFrame) -> tuple[pd.DataFrame, list[str]]:
    exclude = {
        "candidate_id",
        "sample_id",
        "patient_id",
        "study_id",
        "source_name",
        "source_dataset",
        "peptide",
        "mut_peptide",
        "wt_peptide",
        "hla",
        "hla_allele_4digit",
        "label",
        "label_type",
        "public_tool_training_overlap_detail",
    }
    cols = []
    data = {}
    for col in feat.columns:
        if col in exclude:
            continue
        s = feat[col]
        if s.dtype == object:
            if col in {"hla_gene", "hla_supertype", "mhc_class", "leakage_risk_level", "disease_context", "cancer_type"}:
                d = pd.get_dummies(s.fillna("NA").astype(str), prefix=col)
                for dcol in d.columns:
                    data[dcol] = d[dcol].astype(float)
                    cols.append(dcol)
            continue
        data[col] = pd.to_numeric(s, errors="coerce").astype(float)
        cols.append(col)
    X = pd.DataFrame(data, index=feat.index).replace([np.inf, -np.inf], np.nan).fillna(0)
    return X, cols


def rule_score(feat: pd.DataFrame) -> pd.Series:
    score = feat["best_clean_internal_score"].fillna(feat["best_caveated_public_score"]).fillna(0.5).astype(float)
    penalty = (
        0.15 * feat["public_tool_training_overlap_any"].astype(float)
        + 0.12 * feat["underrepresented_hla_allele"].astype(float)
        + 0.15 * feat["source_ood_high"].astype(float)
        + 0.10 * feat["low_prevalence_source"].astype(float)
        + 0.10 * feat["method_disagreement_std"].fillna(0).clip(0, 1)
    )
    return (score - penalty).clip(0, 1)


def disease_gate(row: pd.Series) -> float:
    ctx = str(row.get("disease_context", "")).lower()
    cancer = str(row.get("cancer_type", "")).lower()
    if "low-risk" in ctx or "routine ptc" in ctx:
        return 0.3
    if any(tok in ctx for tok in ["resected", "mrd", "low-burden", "atc", "pdtc", "rai-refractory", "rr-dtc", "high-risk"]):
        return 1.0
    if "paad" in cancer or "thca" in cancer:
        return 0.8
    return 0.8


def presentation_gate(row: pd.Series) -> float:
    hla_loh = str(row.get("hla_loh", "")).lower()
    b2m = str(row.get("b2m_status", "")).lower()
    ap = str(row.get("antigen_processing_status", "")).lower()
    if "loss" in hla_loh or "loh" == hla_loh or "loss" in b2m or "deficient" in ap:
        return 0.0
    if not hla_loh and not b2m and not ap:
        return 0.8
    return 1.0


def antigen_gate(row: pd.Series) -> float:
    available = 0
    for col in ["expression_tpm", "mutant_expression", "vaf", "clonality"]:
        val = pd.to_numeric(pd.Series([row.get(col)]), errors="coerce").iloc[0]
        if pd.notna(val):
            available += 1
    return 1.0 if available >= 2 else 0.8


def immune_gate(row: pd.Series) -> float:
    vals = []
    for col in ["immune_context_score", "tls_score", "ifng_score", "cytolytic_score"]:
        val = pd.to_numeric(pd.Series([row.get(col)]), errors="coerce").iloc[0]
        if pd.notna(val):
            vals.append(val)
    if not vals:
        return 0.8
    return 1.0 if np.nanmean(vals) >= 0 else 0.7


def abstention_reasons(row: pd.Series) -> list[str]:
    reasons = []
    if str(row.get("leakage_risk_level", "")).lower() == "high":
        reasons.append("High leakage risk invalidates clean benchmark claim")
    if bool(row.get("underrepresented_hla_allele", False)):
        reasons.append("HLA allele underrepresented in benchmark")
    if bool(row.get("source_ood_high", False)):
        reasons.append("High source-shift risk: similar heldout source had poor top-k survival")
    if float(row.get("method_disagreement_std", 0) or 0) > 0.25 or float(row.get("method_disagreement_range", 0) or 0) > 0.60:
        reasons.append("High disagreement across predictors")
    if float(row.get("public_vs_internal_disagreement", 0) or 0) > 0.45:
        reasons.append("High disagreement between internal anchor and public pretrained predictors")
    if bool(row.get("only_caveated_public_available", False)):
        reasons.append("Only caveated public predictors available")
    if str(row.get("mhc_class", "I")).upper() != "I":
        reasons.append("MHC-II candidate excluded from class-I benchmark")
    if bool(row.get("low_prevalence_source", False)) and float(row.get("historical_source_top10", 1) or 0) <= 0:
        reasons.append("Low-prevalence source with poor historical top-k performance")
    n_methods = pd.to_numeric(pd.Series([row.get("n_methods_available", 0)]), errors="coerce").fillna(0).iloc[0]
    if int(n_methods) < 2:
        reasons.append("Low confidence due to sparse method coverage")
    if presentation_gate(row) == 0.0:
        reasons.append("Patient presentation hard gate failed")
    return reasons


def main() -> None:
    args = parse_args()
    repo = Path(args.repo_root).resolve()
    output_root = Path(args.output_root)
    ensure_dir(output_root)
    master = pd.read_csv(output_root / "clean_neobench_master.tsv", sep="\t")
    scores = pd.read_csv(output_root / "clean_neobench_method_scores.tsv", sep="\t")
    feat = build_candidate_features(master, scores, repo)
    X, feature_cols = prepare_matrix(feat)
    y = pd.to_numeric(feat["label"], errors="coerce")
    train_mask = y.notna() & y.isin([0, 1])
    model = choose_model()
    if model is not None and train_mask.sum() >= 20 and y[train_mask].nunique() == 2:
        try:
            model.fit(X.loc[train_mask], y.loc[train_mask].astype(int))
            raw_score = pd.Series(model.predict_proba(X)[:, 1], index=feat.index)
            model_name = type(model).__name__
        except Exception:
            raw_score = rule_score(feat)
            model_name = "rule_score_fallback_after_fit_error"
    else:
        raw_score = rule_score(feat)
        model_name = "rule_score_fallback"

    feat["barneo_score"] = raw_score.clip(0, 1)
    feat["disease_gate"] = feat.apply(disease_gate, axis=1)
    feat["presentation_gate"] = feat.apply(presentation_gate, axis=1)
    feat["antigen_gate"] = feat.apply(antigen_gate, axis=1)
    feat["immune_context_gate"] = feat.apply(immune_gate, axis=1)
    feat["model_confidence_gate"] = (1.0 - feat["method_disagreement_std"].fillna(0).clip(0, 0.7)).clip(0.3, 1.0)
    feat["patient_gated_score"] = (
        feat["barneo_score"]
        * feat["disease_gate"]
        * feat["presentation_gate"]
        * feat["antigen_gate"]
        * feat["immune_context_gate"]
        * feat["model_confidence_gate"]
    )
    feat["confidence_score"] = (
        0.55 * feat["barneo_score"]
        + 0.25 * feat["model_confidence_gate"]
        + 0.20 * (1.0 - feat["underrepresented_hla_allele"].astype(float))
    ).clip(0, 1)
    feat["confidence_bin"] = pd.cut(feat["confidence_score"], bins=[-0.01, 0.45, 0.70, 1.01], labels=["low", "medium", "high"])
    feat["abstention_reason_all"] = feat.apply(lambda r: "; ".join(abstention_reasons(r)), axis=1)
    feat["abstention_reason_primary"] = feat["abstention_reason_all"].map(lambda s: s.split("; ")[0] if s else "")
    feat["abstain"] = feat["abstention_reason_all"].astype(str).str.len().gt(0) | feat["confidence_bin"].astype(str).eq("low")
    feat["barneo_rank_global"] = feat["patient_gated_score"].rank(ascending=False, method="first").astype(int)
    feat["barneo_rank_within_patient"] = feat.groupby("patient_id")["patient_gated_score"].rank(ascending=False, method="first")

    out_cols = [
        "candidate_id",
        "barneo_score",
        "patient_gated_score",
        "barneo_rank_global",
        "barneo_rank_within_patient",
        "confidence_score",
        "confidence_bin",
        "abstain",
        "abstention_reason_primary",
        "abstention_reason_all",
        "disease_gate",
        "presentation_gate",
        "antigen_gate",
        "immune_context_gate",
        "model_confidence_gate",
    ]
    write_tsv(feat[out_cols].sort_values("barneo_rank_global"), output_root / "barneo_candidate_scores.tsv")
    reasons = feat[["candidate_id", "abstain", "confidence_bin", "abstention_reason_primary", "abstention_reason_all"]].copy()
    exploded = []
    for _, r in reasons.iterrows():
        parts = [p for p in str(r["abstention_reason_all"]).split("; ") if p]
        if not parts:
            parts = ["no_abstention_reason"]
        for p in parts:
            exploded.append({**r.to_dict(), "abstention_reason": p})
    write_tsv(pd.DataFrame(exploded), output_root / "barneo_abstention_reasons.tsv")
    update_manifest(
        output_root,
        "train_barneo",
        {
            "model": model_name,
            "n_candidates_scored": int(len(feat)),
            "n_training_labels": int(train_mask.sum()),
            "n_abstain": int(feat["abstain"].sum()),
            "confidence_bin_counts": feat["confidence_bin"].astype(str).value_counts().to_dict(),
            "feature_count": int(len(feature_cols)),
            "warnings": [],
        },
    )
    print(f"[bar-neo] candidates={len(feat)} abstain={int(feat['abstain'].sum())} model={model_name}")


if __name__ == "__main__":
    main()
