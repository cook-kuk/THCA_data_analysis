#!/usr/bin/env python3
"""Evaluate CLEAN-NeoBench methods across leakage-aware split contracts."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from common import KOREAN_HLA_ALLELES, binary_metrics, dataframe_to_markdown, ensure_dir, update_manifest, write_tsv


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", default=".", help="Repository root")
    parser.add_argument("--output-root", required=True, help="Output directory")
    return parser.parse_args()


def add_metric_row(rows: list[dict], split_contract: str, split_group: str, method: str, method_meta: pd.Series, df: pd.DataFrame) -> None:
    m = binary_metrics(df["label"], df["score_calibrated"])
    abstain_hint = (
        df["leakage_risk_level"].isin(["high"]).mean()
        if "leakage_risk_level" in df.columns and len(df)
        else np.nan
    )
    m["abstention_rate_recommended"] = float(abstain_hint) if pd.notna(abstain_hint) else np.nan
    rows.append(
        {
            "split_contract": split_contract,
            "split_group": split_group,
            "method_name": method,
            "method_family": method_meta.get("method_family", ""),
            "method_role": method_meta.get("method_role", ""),
            "uses_public_pretraining": bool(method_meta.get("uses_public_pretraining", False)),
            "training_overlap_audited": bool(method_meta.get("training_overlap_audited", False)),
            "clean_comparator_allowed": bool(method_meta.get("clean_comparator_allowed", False)),
            **m,
        }
    )


def evaluate_grouped(rows: list[dict], scored: pd.DataFrame, split_contract: str, group_col: str, min_n: int = 5) -> None:
    if group_col not in scored.columns:
        return
    for group, g0 in scored.groupby(group_col, dropna=False):
        group_name = str(group) if str(group) else "UNKNOWN"
        if len(g0[["candidate_id", "label"]].drop_duplicates()) < min_n:
            continue
        for method, g in g0.groupby("method_name"):
            meta = g.iloc[0]
            add_metric_row(rows, split_contract, group_name, method, meta, g)


def build_leaderboard(metrics: pd.DataFrame) -> pd.DataFrame:
    if metrics.empty:
        return metrics
    preferred = metrics[metrics["split_contract"].isin(["exact_peptide_hla_holdout", "near_peptide_cluster_holdout", "hla_heldout", "source_heldout"])]
    if preferred.empty:
        preferred = metrics
    agg = (
        preferred.groupby(["method_name", "method_family", "method_role", "uses_public_pretraining", "training_overlap_audited", "clean_comparator_allowed"], dropna=False)
        .agg(
            mean_AUPRC=("AUPRC", "mean"),
            median_AUPRC=("AUPRC", "median"),
            mean_AUROC=("AUROC", "mean"),
            mean_top10_precision=("top10_precision", "mean"),
            mean_top20_precision=("top20_precision", "mean"),
            mean_ECE=("calibration_ece", "mean"),
            mean_Brier=("calibration_brier", "mean"),
            mean_abstention_rate=("abstention_rate_recommended", "mean"),
            n_split_rows=("split_contract", "size"),
            n_total_scored=("n_total", "sum"),
        )
        .reset_index()
    )
    role_weight = agg["method_role"].map(
        {
            "anchor": 0.04,
            "internal_candidate": 0.02,
            "bounded_fallback": 0.0,
            "uncertainty_only": -0.04,
            "caveated_public_comparator": -0.03,
        }
    ).fillna(0)
    caveat_penalty = np.where(agg["uses_public_pretraining"] & ~agg["training_overlap_audited"], 0.05, 0.0)
    agg["reviewer_safe_score"] = (
        agg["mean_AUPRC"].fillna(0)
        + 0.20 * agg["mean_top10_precision"].fillna(0)
        - 0.10 * agg["mean_ECE"].fillna(0.5)
        + role_weight
        - caveat_penalty
    )
    return agg.sort_values(["reviewer_safe_score", "mean_AUPRC", "mean_top10_precision"], ascending=False)


def main() -> None:
    args = parse_args()
    output_root = Path(args.output_root)
    ensure_dir(output_root)
    master = pd.read_csv(output_root / "clean_neobench_master.tsv", sep="\t")
    scores = pd.read_csv(output_root / "clean_neobench_method_scores.tsv", sep="\t")
    scored = scores.merge(master, on="candidate_id", how="left")
    scored = scored[pd.to_numeric(scored["label"], errors="coerce").notna()].copy()
    scored["label"] = pd.to_numeric(scored["label"], errors="coerce").astype(int)
    rows: list[dict] = []

    for method, g in scored.groupby("method_name"):
        add_metric_row(rows, "overall_labeled", "all", method, g.iloc[0], g)

    exact = scored[~scored["exact_peptide_hla_train_overlap"].astype(bool)]
    for method, g in exact.groupby("method_name"):
        add_metric_row(rows, "exact_peptide_hla_holdout", "no_exact_phla_train_overlap", method, g.iloc[0], g)

    near = scored[~scored["near_peptide_train_overlap"].astype(bool)]
    for method, g in near.groupby("method_name"):
        add_metric_row(rows, "near_peptide_cluster_holdout", "no_near_peptide_train_overlap", method, g.iloc[0], g)

    source_window = scored[scored["source_protein_window"].fillna("").astype(str).str.len() > 0]
    if len(source_window):
        no_window_overlap = source_window[~source_window["source_protein_window_train_overlap"].astype(bool)]
        for method, g in no_window_overlap.groupby("method_name"):
            add_metric_row(rows, "source_protein_window_holdout", "no_source_window_train_overlap", method, g.iloc[0], g)

    evaluate_grouped(rows, scored, "source_heldout", "source_name", min_n=5)
    evaluate_grouped(rows, scored, "study_heldout", "study_id", min_n=5)
    evaluate_grouped(rows, scored, "hla_heldout", "hla_allele_4digit", min_n=5)
    evaluate_grouped(rows, scored, "supertype_heldout", "hla_supertype", min_n=5)
    if scored["patient_id"].fillna("").astype(str).str.len().gt(0).any():
        evaluate_grouped(rows, scored, "patient_heldout", "patient_id", min_n=3)

    low_prev = scored[scored["split_low_prevalence"].eq("low_prevalence")]
    for method, g in low_prev.groupby("method_name"):
        add_metric_row(rows, "low_prevalence_heldout", "low_prevalence_sources", method, g.iloc[0], g)

    korean = scored[scored["hla_allele_4digit"].isin(KOREAN_HLA_ALLELES)]
    for method, g in korean.groupby("method_name"):
        add_metric_row(rows, "korean_hla_focus", "korean_relevant_alleles", method, g.iloc[0], g)

    # Preserve existing locked split contexts when present in prediction files.
    if "score_context" in scored.columns:
        ctx = scored[scored["score_context"].fillna("").astype(str).str.len() > 0]
        for context, g0 in ctx.groupby("score_context"):
            if len(g0[["candidate_id", "label"]].drop_duplicates()) < 5:
                continue
            for method, g in g0.groupby("method_name"):
                add_metric_row(rows, "existing_prediction_context", str(context), method, g.iloc[0], g)

    metrics = pd.DataFrame(rows)
    write_tsv(metrics, output_root / "clean_neobench_split_metrics.tsv")
    leaderboard = build_leaderboard(metrics)
    write_tsv(leaderboard, output_root / "clean_neobench_leaderboard.tsv")
    update_manifest(
        output_root,
        "evaluate_clean_neobench",
        {
            "n_metric_rows": int(len(metrics)),
            "n_leaderboard_methods": int(leaderboard["method_name"].nunique()) if len(leaderboard) else 0,
            "split_contract_counts": metrics["split_contract"].value_counts().to_dict() if len(metrics) else {},
            "top_methods": leaderboard.head(10)["method_name"].tolist() if len(leaderboard) else [],
            "warnings": [],
        },
    )
    print(f"[clean-neobench-eval] metrics={len(metrics)} leaderboard_methods={leaderboard['method_name'].nunique() if len(leaderboard) else 0}")


if __name__ == "__main__":
    main()
