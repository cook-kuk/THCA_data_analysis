#!/usr/bin/env python3
"""Build a stress-guarded BAR-Neo ranking layer.

This is the first "make it practical" layer after benchmarking:
method weights are no longer global only. They are downweighted when the
source/HLA stress audit says that method collapses in the candidate's context.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from common import dataframe_to_markdown, ensure_dir, update_manifest, write_tsv


OUTPUTS = [
    "barneo_stress_guarded_candidate_scores.tsv",
    "barneo_stress_guarded_method_weights.tsv",
    "barneo_stress_guarded_expert_contributions.tsv",
    "barneo_stress_guarded_review_queue.tsv",
    "BAR_NEO_STRESS_GUARDED_REPORT.md",
    "BAR_NEO_STRESS_GUARDED_REPORT_KR.md",
]

KOREAN_RELEVANT_HLAS = {
    "HLA-A*24:02",
    "HLA-A*11:01",
    "HLA-A*02:01",
    "HLA-B*15:01",
    "HLA-B*40:01",
    "HLA-C*01:02",
    "HLA-C*03:03",
    "HLA-C*07:02",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", default=".", help="Repository root")
    parser.add_argument("--output-root", required=True, help="CLEAN-NeoBench output root")
    parser.add_argument("--top-contrib", type=int, default=7, help="Top expert contributions per candidate")
    return parser.parse_args()


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t") if path.exists() else pd.DataFrame()


def to_num(series: pd.Series, default: float = 0.0) -> pd.Series:
    return pd.to_numeric(series, errors="coerce").fillna(default)


def clamp(x: Any, lo: float = 0.0, hi: float = 1.0) -> float:
    try:
        y = float(x)
    except Exception:
        return lo
    if not math.isfinite(y):
        return lo
    return max(lo, min(hi, y))


def yes(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"true", "1", "yes"}


def split_flags(value: Any) -> set[str]:
    return {x.strip() for x in str(value or "").split(";") if x.strip()}


def sigmoid_like_delta(delta: Any) -> float:
    """Map benchmark delta roughly into [0, 1] without overfitting."""
    try:
        d = float(delta)
    except Exception:
        d = 0.0
    if not math.isfinite(d):
        d = 0.0
    return clamp(0.5 + d, 0.05, 0.95)


def role_prior(role: str) -> float:
    return {
        "anchor": 0.80,
        "internal_candidate": 1.00,
        "bounded_fallback": 0.38,
        "uncertainty_only": 0.22,
        "caveated_public_comparator": 0.16,
        "negative_baseline": 0.03,
    }.get(str(role), 0.45)


def build_method_weights(stress: pd.DataFrame, leaderboard: pd.DataFrame) -> pd.DataFrame:
    if stress.empty:
        return pd.DataFrame()
    out = stress.copy()
    if not leaderboard.empty and "method_name" in leaderboard.columns:
        extra = leaderboard[
            [
                c
                for c in [
                    "method_name",
                    "mean_AUPRC",
                    "mean_top10_precision",
                    "mean_ECE",
                    "reviewer_safe_score",
                    "n_split_rows",
                ]
                if c in leaderboard.columns
            ]
        ].drop_duplicates("method_name")
        out = out.merge(extra, on="method_name", how="left", suffixes=("", "_leaderboard"))

    for col in [
        "mean_AUPRC",
        "reviewer_safe_score",
        "source_heldout_median_delta_AUPRC",
        "hla_stress_median_delta_AUPRC",
        "korean_hla_median_delta_AUPRC",
        "low_prevalence_median_delta_AUPRC",
        "mean_ECE",
        "mean_top10_precision",
    ]:
        if col not in out.columns:
            out[col] = np.nan
        out[col] = pd.to_numeric(out[col], errors="coerce")

    out["source_stress_prior"] = out["source_heldout_median_delta_AUPRC"].map(sigmoid_like_delta)
    out["hla_stress_prior"] = out["hla_stress_median_delta_AUPRC"].map(sigmoid_like_delta)
    out["korean_hla_prior"] = out["korean_hla_median_delta_AUPRC"].map(sigmoid_like_delta)
    out["low_prevalence_prior"] = out["low_prevalence_median_delta_AUPRC"].map(sigmoid_like_delta)
    out["calibration_prior"] = 1.0 - out["mean_ECE"].fillna(0.25).clip(0, 0.5) * 1.5
    out["role_prior"] = out["method_role"].map(role_prior)

    out["stress_guarded_base_utility"] = (
        0.30 * out["mean_AUPRC"].fillna(0.50).clip(0, 1)
        + 0.24 * out["reviewer_safe_score"].fillna(0.45).clip(0, 1)
        + 0.16 * out["source_stress_prior"]
        + 0.12 * out["hla_stress_prior"]
        + 0.08 * out["korean_hla_prior"]
        + 0.05 * out["low_prevalence_prior"]
        + 0.05 * out["calibration_prior"].clip(0, 1)
    )
    out["stress_guarded_base_weight"] = out["stress_guarded_base_utility"] * out["role_prior"]
    out.loc[out["method_role"].eq("caveated_public_comparator"), "stress_guarded_base_weight"] *= 0.65
    out.loc[out["method_role"].eq("bounded_fallback"), "stress_guarded_base_weight"] *= 0.75
    out.loc[out["method_role"].eq("uncertainty_only"), "stress_guarded_base_weight"] *= 0.55
    out["stress_flags_set"] = out["stress_flags"].map(split_flags) if "stress_flags" in out.columns else [set() for _ in range(len(out))]
    out["stress_guarded_base_weight"] = out["stress_guarded_base_weight"].clip(0.005, 1.0)
    return out.sort_values("stress_guarded_base_weight", ascending=False)


def prepare_scores(scores: pd.DataFrame) -> pd.DataFrame:
    if scores.empty:
        return scores
    out = scores.copy()
    out["score_calibrated_num"] = pd.to_numeric(out.get("score_calibrated"), errors="coerce")
    out["score_raw_num"] = pd.to_numeric(out.get("score_raw"), errors="coerce")
    out["rank_global_num"] = pd.to_numeric(out.get("rank_global"), errors="coerce")
    out["method_n"] = out.groupby("method_name")["candidate_id"].transform("count").clip(lower=1)
    out["rank_percentile_score"] = 1.0 - ((out["rank_global_num"] - 1.0) / out["method_n"])
    out["score_for_fusion"] = out["score_calibrated_num"]
    missing = out["score_for_fusion"].isna()
    out.loc[missing, "score_for_fusion"] = out.loc[missing, "rank_percentile_score"]
    out["score_for_fusion"] = out["score_for_fusion"].clip(0, 1)
    out = out[out.get("runtime_status", "ok").astype(str).str.lower().eq("ok")].copy()
    return out


def candidate_context(master: pd.DataFrame) -> pd.DataFrame:
    keep = [
        "candidate_id",
        "source_name",
        "source_dataset",
        "study_id",
        "patient_id",
        "hla_allele_4digit",
        "hla_gene",
        "hla_supertype",
        "peptide",
        "label",
        "leakage_risk_level",
        "split_low_prevalence",
        "split_korean_hla_focus",
        "mhc_class",
    ]
    out = master[[c for c in keep if c in master.columns]].copy()
    out["label"] = pd.to_numeric(out.get("label"), errors="coerce")
    source_prev = out.groupby("source_name")["label"].mean().rename("source_positive_prevalence")
    hla_support = out.groupby("hla_allele_4digit")["candidate_id"].nunique().rename("hla_allele_support_count")
    out = out.merge(source_prev, on="source_name", how="left")
    out = out.merge(hla_support, on="hla_allele_4digit", how="left")
    lowprev = out.get("split_low_prevalence", pd.Series("", index=out.index)).astype(str).str.lower()
    out["is_low_prevalence_context"] = (
        lowprev.isin({"low_prevalence", "low_prevalence_heldout", "low_prevalence_focus"})
        | out["source_positive_prevalence"].fillna(1.0).lt(0.20)
    )
    out["is_korean_hla_focus"] = (
        out.get("split_korean_hla_focus", "").astype(str).eq("korean_hla_focus")
        | out.get("hla_allele_4digit", "").astype(str).isin(KOREAN_RELEVANT_HLAS)
    )
    out["is_underrepresented_hla"] = out["hla_allele_support_count"].fillna(0).lt(25)
    out["is_high_leakage"] = out.get("leakage_risk_level", "").astype(str).str.lower().eq("high")
    out["is_mhc_ii"] = out.get("mhc_class", "").astype(str).str.upper().eq("II")
    return out


def context_multiplier(row: pd.Series) -> tuple[float, list[str]]:
    mult = 1.0
    reasons: list[str] = []
    flags = row.get("stress_flags_set", set())
    role = str(row.get("method_role", ""))

    if "source_shift_collapse" in flags:
        mult *= 0.55
        reasons.append("method source-shift collapse")
        if yes(row.get("is_low_prevalence_context")):
            mult *= 0.55
            reasons.append("candidate low-prevalence/source-stress context")
    if "source_low_floor" in flags and yes(row.get("is_low_prevalence_context")):
        mult *= 0.55
        reasons.append("method has low source-heldout floor")
    if "hla_shift_collapse" in flags:
        mult *= 0.70
        reasons.append("method HLA-heldout collapse")
        if yes(row.get("is_underrepresented_hla")):
            mult *= 0.55
            reasons.append("candidate underrepresented HLA")
    if "korean_hla_collapse" in flags and yes(row.get("is_korean_hla_focus")):
        mult *= 0.45
        reasons.append("method Korean-HLA collapse")
    if "low_prevalence_collapse" in flags and yes(row.get("is_low_prevalence_context")):
        mult *= 0.45
        reasons.append("method low-prevalence collapse")
    if "calibration_risk" in flags:
        mult *= 0.92
        reasons.append("method calibration risk")
    if role == "caveated_public_comparator":
        mult *= 0.55
        reasons.append("public pretrained caveat")
        if yes(row.get("is_high_leakage")):
            mult *= 0.70
            reasons.append("public support suppressed under high leakage")
    if role == "bounded_fallback":
        mult *= 0.75
        reasons.append("bounded fallback cap")
    if role == "uncertainty_only":
        mult *= 0.60
        reasons.append("uncertainty branch not ranking headline")
    if yes(row.get("is_mhc_ii")):
        mult *= 0.25
        reasons.append("MHC-II excluded from class-I stress policy")
    return clamp(mult, 0.01, 1.0), reasons


def weighted_mean(values: pd.Series, weights: pd.Series) -> float:
    w = pd.to_numeric(weights, errors="coerce").fillna(0)
    v = pd.to_numeric(values, errors="coerce").fillna(0)
    denom = float(w.sum())
    if denom <= 0:
        return math.nan
    return float((v * w).sum() / denom)


def weighted_std(values: pd.Series, weights: pd.Series) -> float:
    mu = weighted_mean(values, weights)
    if not math.isfinite(mu):
        return math.nan
    w = pd.to_numeric(weights, errors="coerce").fillna(0)
    v = pd.to_numeric(values, errors="coerce").fillna(mu)
    denom = float(w.sum())
    if denom <= 0:
        return math.nan
    return float(np.sqrt(((w * (v - mu) ** 2).sum()) / denom))


def build_candidate_scores(master: pd.DataFrame, scores: pd.DataFrame, weights: pd.DataFrame, contextual: pd.DataFrame, xscore: pd.DataFrame, top_contrib: int) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    if master.empty or scores.empty or weights.empty:
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()

    ctx = candidate_context(master)
    sc = prepare_scores(scores)
    wcols = [
        "method_name",
        "method_role",
        "method_family",
        "stress_guarded_base_weight",
        "stress_guarded_base_utility",
        "source_heldout_median_delta_AUPRC",
        "hla_stress_median_delta_AUPRC",
        "korean_hla_median_delta_AUPRC",
        "low_prevalence_median_delta_AUPRC",
        "stress_flags",
        "stress_flags_set",
        "recommended_use",
    ]
    merged = sc.merge(weights[[c for c in wcols if c in weights.columns]], on=["method_name", "method_role", "method_family"], how="left")
    merged = merged.merge(ctx, on="candidate_id", how="left")
    merged["stress_guarded_base_weight"] = pd.to_numeric(merged["stress_guarded_base_weight"], errors="coerce").fillna(0.02)

    multipliers = merged.apply(context_multiplier, axis=1)
    merged["context_multiplier"] = [x[0] for x in multipliers]
    merged["context_adjustment_reasons"] = ["; ".join(x[1]) if x[1] else "no additional stress penalty" for x in multipliers]
    merged["candidate_method_weight_raw"] = merged["stress_guarded_base_weight"] * merged["context_multiplier"]

    # Cap caveated public and bounded fallback support. They may support triage but cannot dominate.
    merged.loc[merged["method_role"].eq("caveated_public_comparator"), "candidate_method_weight_raw"] *= 0.40
    merged.loc[merged["method_role"].eq("bounded_fallback"), "candidate_method_weight_raw"] *= 0.65
    merged.loc[merged["method_role"].eq("uncertainty_only"), "candidate_method_weight_raw"] *= 0.45
    merged["candidate_method_contribution_raw"] = merged["score_for_fusion"] * merged["candidate_method_weight_raw"]

    rows = []
    contrib_rows = []
    clean_roles = {"anchor", "internal_candidate"}
    support_roles = {"anchor", "internal_candidate", "bounded_fallback", "uncertainty_only", "caveated_public_comparator"}
    for candidate_id, g in merged.groupby("candidate_id", sort=False):
        g = g.copy()
        all_support = g[g["method_role"].isin(support_roles)]
        clean = g[g["method_role"].isin(clean_roles) & g.get("clean_comparator_allowed", pd.Series(True, index=g.index)).map(yes)]
        fallback = g[g["method_role"].eq("bounded_fallback")]
        public = g[g["method_role"].eq("caveated_public_comparator")]

        support_score = weighted_mean(all_support["score_for_fusion"], all_support["candidate_method_weight_raw"])
        clean_score = weighted_mean(clean["score_for_fusion"], clean["candidate_method_weight_raw"])
        fallback_score = weighted_mean(fallback["score_for_fusion"], fallback["candidate_method_weight_raw"])
        public_score = weighted_mean(public["score_for_fusion"], public["candidate_method_weight_raw"])
        disagreement = weighted_std(all_support["score_for_fusion"], all_support["candidate_method_weight_raw"])
        clean_n = int(clean["method_name"].nunique())
        support_n = int(all_support["method_name"].nunique())
        public_weight_frac = float(public["candidate_method_weight_raw"].sum() / all_support["candidate_method_weight_raw"].sum()) if all_support["candidate_method_weight_raw"].sum() > 0 else 0.0
        fallback_weight_frac = float(fallback["candidate_method_weight_raw"].sum() / all_support["candidate_method_weight_raw"].sum()) if all_support["candidate_method_weight_raw"].sum() > 0 else 0.0

        ctx_row = ctx[ctx["candidate_id"].eq(candidate_id)].iloc[0].to_dict()
        stress_penalty = 1.0
        reasons: list[str] = []
        if ctx_row.get("is_high_leakage"):
            stress_penalty *= 0.62
            reasons.append("High leakage risk blocks clean claim")
        if ctx_row.get("is_underrepresented_hla"):
            stress_penalty *= 0.82
            reasons.append("HLA allele underrepresented")
        if ctx_row.get("is_low_prevalence_context"):
            stress_penalty *= 0.80
            reasons.append("Low-prevalence/source-stress context")
        if ctx_row.get("is_mhc_ii"):
            stress_penalty *= 0.25
            reasons.append("MHC-II candidate excluded from class-I benchmark")
        if clean_n < 3:
            stress_penalty *= 0.72
            reasons.append("Low clean internal method availability")
        if public_weight_frac > 0.35:
            stress_penalty *= 0.80
            reasons.append("Too much caveated public support")
        if fallback_weight_frac > 0.25:
            stress_penalty *= 0.85
            reasons.append("Fallback branch contributes materially")
        if math.isfinite(disagreement) and disagreement > 0.25:
            stress_penalty *= 0.75
            reasons.append("High expert disagreement")

        clean_score_for_mix = clean_score if math.isfinite(clean_score) else 0.0
        support_score_for_mix = support_score if math.isfinite(support_score) else 0.0
        fallback_score_for_mix = fallback_score if math.isfinite(fallback_score) else 0.0
        public_score_for_mix = public_score if math.isfinite(public_score) else 0.0
        discovery = (
            0.58 * support_score_for_mix
            + 0.27 * clean_score_for_mix
            + 0.10 * fallback_score_for_mix
            + 0.05 * public_score_for_mix
        )
        confidence = clamp(
            0.28
            + 0.06 * min(clean_n, 7)
            + 0.03 * min(support_n, 10)
            - (0.45 * disagreement if math.isfinite(disagreement) else 0.20)
        )
        confidence = clamp(confidence * stress_penalty)
        claim_safe = clamp((0.78 * clean_score_for_mix + 0.22 * support_score_for_mix) * confidence)
        action = "priority_review_candidate"
        if confidence < 0.30:
            action = "stress_abstain"
        elif claim_safe >= 0.50 and not ctx_row.get("is_high_leakage") and clean_n >= 4:
            action = "claim_safe_candidate_after_manual_audit"
        elif discovery >= 0.60:
            action = "priority_review_candidate"
        elif clean_score_for_mix >= 0.55:
            action = "support_review_candidate"
        else:
            action = "watchlist"

        if not reasons:
            reasons.append("Stress gates passed at current metadata resolution")
        ordered = g.sort_values("candidate_method_contribution_raw", ascending=False).head(top_contrib)
        contrib_rows.extend(
            {
                "candidate_id": candidate_id,
                "method_name": r["method_name"],
                "method_role": r["method_role"],
                "method_family": r["method_family"],
                "score_for_fusion": r["score_for_fusion"],
                "candidate_method_weight": r["candidate_method_weight_raw"],
                "candidate_method_contribution": r["candidate_method_contribution_raw"],
                "context_adjustment_reasons": r["context_adjustment_reasons"],
            }
            for _, r in ordered.iterrows()
        )
        rows.append(
            {
                "candidate_id": candidate_id,
                "stress_guarded_discovery_score": discovery,
                "stress_guarded_claim_safe_score": claim_safe,
                "stress_guarded_clean_score": clean_score,
                "stress_guarded_support_score": support_score,
                "stress_guarded_fallback_score": fallback_score,
                "stress_guarded_public_score": public_score,
                "stress_guarded_confidence": confidence,
                "stress_guarded_expert_disagreement": disagreement,
                "stress_guarded_clean_method_count": clean_n,
                "stress_guarded_support_method_count": support_n,
                "public_weight_fraction": public_weight_frac,
                "fallback_weight_fraction": fallback_weight_frac,
                "stress_guarded_action": action,
                "stress_guarded_abstain": action == "stress_abstain" or ctx_row.get("is_high_leakage") or ctx_row.get("is_mhc_ii"),
                "stress_guarded_reason_primary": reasons[0],
                "stress_guarded_reason_all": "; ".join(reasons),
                **ctx_row,
            }
        )

    candidate_scores = pd.DataFrame(rows)
    missing_ctx = ctx[~ctx["candidate_id"].isin(candidate_scores.get("candidate_id", pd.Series(dtype=str)))].copy()
    if not missing_ctx.empty:
        missing_rows = []
        for _, r in missing_ctx.iterrows():
            ctx_row = r.to_dict()
            missing_rows.append(
                {
                    "candidate_id": ctx_row["candidate_id"],
                    "stress_guarded_discovery_score": 0.0,
                    "stress_guarded_claim_safe_score": 0.0,
                    "stress_guarded_clean_score": math.nan,
                    "stress_guarded_support_score": math.nan,
                    "stress_guarded_fallback_score": math.nan,
                    "stress_guarded_public_score": math.nan,
                    "stress_guarded_confidence": 0.0,
                    "stress_guarded_expert_disagreement": math.nan,
                    "stress_guarded_clean_method_count": 0,
                    "stress_guarded_support_method_count": 0,
                    "public_weight_fraction": 0.0,
                    "fallback_weight_fraction": 0.0,
                    "stress_guarded_action": "stress_abstain",
                    "stress_guarded_abstain": True,
                    "stress_guarded_reason_primary": "No usable method scores after runtime/status filtering",
                    "stress_guarded_reason_all": "No usable method scores after runtime/status filtering",
                    **ctx_row,
                }
            )
        candidate_scores = pd.concat([candidate_scores, pd.DataFrame(missing_rows)], ignore_index=True, sort=False)
    if not contextual.empty and "candidate_id" in contextual:
        candidate_scores = candidate_scores.merge(
            contextual[
                [
                    c
                    for c in [
                        "candidate_id",
                        "contextual_bma_score",
                        "clean_contextual_bma_score",
                        "contextual_confidence_score",
                        "contextual_abstain",
                    ]
                    if c in contextual.columns
                ]
            ],
            on="candidate_id",
            how="left",
        )
    if not xscore.empty and "candidate_id" in xscore:
        candidate_scores = candidate_scores.merge(
            xscore[
                [
                    c
                    for c in [
                        "candidate_id",
                        "barneo_x_discovery_score",
                        "barneo_x_claim_safe_score",
                        "barneo_x_primary_action",
                    ]
                    if c in xscore.columns
                ]
            ],
            on="candidate_id",
            how="left",
        )
    candidate_scores["stress_guarded_final_review_score"] = (
        0.55 * to_num(candidate_scores["stress_guarded_claim_safe_score"])
        + 0.25 * to_num(candidate_scores.get("barneo_x_claim_safe_score", pd.Series(index=candidate_scores.index)))
        + 0.20 * to_num(candidate_scores["stress_guarded_discovery_score"])
    ).clip(0, 1)
    candidate_scores = candidate_scores.sort_values(
        ["stress_guarded_final_review_score", "stress_guarded_claim_safe_score", "stress_guarded_confidence"],
        ascending=False,
    )
    candidate_scores["stress_guarded_rank_global"] = np.arange(1, len(candidate_scores) + 1)
    if "patient_id" in candidate_scores.columns:
        candidate_scores["stress_guarded_rank_within_patient"] = candidate_scores.groupby("patient_id")[
            "stress_guarded_final_review_score"
        ].rank(method="first", ascending=False)
    if "hla_allele_4digit" in candidate_scores.columns:
        candidate_scores["stress_guarded_rank_within_hla"] = candidate_scores.groupby("hla_allele_4digit")[
            "stress_guarded_final_review_score"
        ].rank(method="first", ascending=False)

    contrib = pd.DataFrame(contrib_rows)
    review_queue = candidate_scores[
        candidate_scores["stress_guarded_action"].isin(
            ["claim_safe_candidate_after_manual_audit", "priority_review_candidate", "support_review_candidate"]
        )
    ].head(350)
    return candidate_scores, contrib, review_queue


def write_reports(output_root: Path, candidate_scores: pd.DataFrame, method_weights: pd.DataFrame, review_queue: pd.DataFrame) -> None:
    action_counts = (
        candidate_scores["stress_guarded_action"].value_counts().rename_axis("action").reset_index(name="n")
        if not candidate_scores.empty and "stress_guarded_action" in candidate_scores
        else pd.DataFrame()
    )
    top_cols = [
        "candidate_id",
        "stress_guarded_final_review_score",
        "stress_guarded_claim_safe_score",
        "stress_guarded_discovery_score",
        "stress_guarded_confidence",
        "stress_guarded_action",
        "stress_guarded_reason_primary",
        "source_name",
        "hla_allele_4digit",
        "peptide",
        "label",
        "stress_guarded_clean_method_count",
        "public_weight_fraction",
        "fallback_weight_fraction",
    ]
    method_cols = [
        "method_name",
        "method_role",
        "method_family",
        "mean_AUPRC",
        "reviewer_safe_score",
        "source_heldout_median_delta_AUPRC",
        "hla_stress_median_delta_AUPRC",
        "korean_hla_median_delta_AUPRC",
        "stress_guarded_base_weight",
        "stress_flags",
        "recommended_use",
    ]
    text = f"""# BAR-Neo Stress-Guarded Ranker

## What Changed

This layer converts CLEAN-NeoBench source/HLA stress behavior into candidate-specific ensemble weights. Methods that collapse under source shift, HLA shift, Korean-HLA focus, low-prevalence settings, public-pretraining caveats, or fallback-only support are downweighted before candidate ranking.

## Action Counts

{dataframe_to_markdown(action_counts, max_rows=20)}

## Top Review Candidates

{dataframe_to_markdown(review_queue[[c for c in top_cols if c in review_queue.columns]].head(50), max_rows=50)}

## Method Weights

{dataframe_to_markdown(method_weights[[c for c in method_cols if c in method_weights.columns]].head(35), max_rows=35)}

## Claim Boundary

Allowed: stress-guarded research triage, benchmark-adaptive method weighting, priority manual review.

Forbidden: clinical vaccine selection, new SOTA predictor, quantum advantage, or clean public baselines without row-level training-corpus overlap audit.
"""
    (output_root / "BAR_NEO_STRESS_GUARDED_REPORT.md").write_text(text.strip() + "\n")

    kr = f"""# BAR-Neo Stress-Guarded Ranker KR

## 한 줄 결론

이제 가능성 표에서 멈추지 않고, source/HLA에서 실제로 깨지는 method를 candidate별로 자동 downweight하는 실전 ranking layer가 생겼다. 이건 SOTA claim이 아니라 reviewer-safe triage engine이다.

## Action counts

{dataframe_to_markdown(action_counts, max_rows=20)}

## Top stress-guarded review queue

{dataframe_to_markdown(review_queue[[c for c in top_cols if c in review_queue.columns]].head(40), max_rows=40)}

## 어떤 method가 실제로 올라가나

{dataframe_to_markdown(method_weights[[c for c in method_cols if c in method_weights.columns]].head(25), max_rows=25)}

## 실용 해석

- `W7B_stacked` / `W7A_full`처럼 source/HLA stress에서 버티는 내부 후보는 weight가 올라간다.
- source-balanced/PU RF처럼 overall은 강하지만 source-heldout collapse가 큰 method는 해당 context에서 downweight된다.
- public pretrained methods는 overlap audit 전까지 support signal 이상으로 못 올라간다.
- QK는 좋은 slice가 있어도 bounded fallback/fusion cap을 받는다.
- high leakage row는 priority review는 될 수 있어도 clean claim은 막힌다.
"""
    (output_root / "BAR_NEO_STRESS_GUARDED_REPORT_KR.md").write_text(kr.strip() + "\n")


def main() -> None:
    args = parse_args()
    repo_root = Path(args.repo_root).resolve()
    output_root = Path(args.output_root)
    if not output_root.is_absolute():
        output_root = repo_root / output_root
    ensure_dir(output_root)

    master = read_tsv(output_root / "clean_neobench_master.tsv")
    scores = read_tsv(output_root / "clean_neobench_method_scores.tsv")
    stress = read_tsv(output_root / "clean_neobench_source_hla_stress_method_summary.tsv")
    leaderboard = read_tsv(output_root / "clean_neobench_leaderboard.tsv")
    contextual = read_tsv(output_root / "barneo_contextual_bma_candidate_scores.tsv")
    xscore = read_tsv(output_root / "barneo_x_candidate_scores.tsv")

    method_weights = build_method_weights(stress, leaderboard)
    candidate_scores, contrib, review_queue = build_candidate_scores(
        master, scores, method_weights, contextual, xscore, args.top_contrib
    )

    write_tsv(method_weights.drop(columns=["stress_flags_set"], errors="ignore"), output_root / "barneo_stress_guarded_method_weights.tsv")
    write_tsv(candidate_scores, output_root / "barneo_stress_guarded_candidate_scores.tsv")
    write_tsv(contrib, output_root / "barneo_stress_guarded_expert_contributions.tsv")
    write_tsv(review_queue, output_root / "barneo_stress_guarded_review_queue.tsv")
    write_reports(output_root, candidate_scores, method_weights, review_queue)

    manifest_path = output_root / "run_manifest.json"
    manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
    manifest.setdefault("output_files", [])
    for name in OUTPUTS:
        if name not in manifest["output_files"]:
            manifest["output_files"].append(name)
    manifest.setdefault("summary", {})
    manifest["summary"].update(
        {
            "n_stress_guarded_candidates": int(len(candidate_scores)),
            "n_stress_guarded_review_queue": int(len(review_queue)),
            "stress_guarded_top_score": float(candidate_scores["stress_guarded_final_review_score"].max())
            if not candidate_scores.empty
            else None,
            "n_stress_guarded_claim_safe_after_manual_audit": int(
                candidate_scores["stress_guarded_action"].eq("claim_safe_candidate_after_manual_audit").sum()
            )
            if not candidate_scores.empty
            else 0,
            "n_stress_guarded_priority_review": int(
                candidate_scores["stress_guarded_action"].eq("priority_review_candidate").sum()
            )
            if not candidate_scores.empty
            else 0,
        }
    )
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    update_manifest(
        output_root,
        "barneo_stress_guarded_ranker",
        {
            "outputs": OUTPUTS,
            "n_candidates": int(len(candidate_scores)),
            "n_review_queue": int(len(review_queue)),
            "warnings": ["Stress-guarded ranking is research triage, not clinical vaccine selection or SOTA validation."],
        },
    )
    print(
        "[barneo-stress-guarded] "
        f"candidates={len(candidate_scores)} review_queue={len(review_queue)} "
        f"top={candidate_scores['stress_guarded_final_review_score'].max() if not candidate_scores.empty else 'NA'}"
    )


if __name__ == "__main__":
    main()
