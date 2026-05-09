#!/usr/bin/env python3
"""Apply failure-aware reliability adjustment to BAR-Neo-BMA scores.

This is a conservative post-hoc reliability layer. It does not train a new
predictor. It converts observed failure modes into candidate-level confidence
caps, score penalties, and human-readable abstention reasons.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from common import ensure_dir, update_manifest, write_tsv


LOW_PREVALENCE_CUTOFF = 0.10
RARE_HLA_SUPPORT = 20
LOW_HLA_SUPPORT = 50


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", default=".", help="Repository root")
    parser.add_argument("--output-root", required=True, help="CLEAN-NeoBench output directory")
    return parser.parse_args()


def load_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t") if path.exists() else pd.DataFrame()


def as_float(value: Any, default: float = math.nan) -> float:
    try:
        v = float(value)
        return v if math.isfinite(v) else default
    except Exception:
        return default


def truthy(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"true", "1", "yes", "y"}


def risk_maps(failure_dist: pd.DataFrame) -> dict[str, dict[str, float]]:
    maps: dict[str, dict[str, float]] = {}
    if failure_dist.empty:
        return maps
    fd = failure_dist.copy()
    fd["enrichment_vs_background"] = pd.to_numeric(fd["enrichment_vs_background"], errors="coerce")
    fd["n_error"] = pd.to_numeric(fd["n_error"], errors="coerce").fillna(0)
    for feature in ["source_name", "hla_allele_4digit", "hla_supertype", "peptide_length", "leakage_risk_level"]:
        sub = fd[fd["feature"].eq(feature)]
        if sub.empty:
            continue
        for failure_type, suffix in [
            ("high_ranked_negative_top10pct", "fp"),
            ("threshold_false_positive_score_ge_0.5", "fp_threshold"),
            ("missed_positive_bottom50pct", "fn"),
            ("threshold_false_negative_score_lt_0.5", "fn_threshold"),
        ]:
            s = sub[sub["primary_failure_type"].eq(failure_type)]
            if s.empty:
                continue
            key = f"{feature}:{suffix}"
            grouped = (
                s.groupby("value")
                .agg(max_enrichment=("enrichment_vs_background", "max"), n_error=("n_error", "sum"))
                .reset_index()
            )
            maps[key] = {
                str(r["value"]): float(r["max_enrichment"]) * math.log1p(float(r["n_error"]))
                for _, r in grouped.iterrows()
                if math.isfinite(float(r["max_enrichment"]))
            }
    return maps


def source_hla_stats(master: pd.DataFrame) -> tuple[dict[str, dict[str, float]], dict[str, dict[str, float]]]:
    src: dict[str, dict[str, float]] = {}
    hla: dict[str, dict[str, float]] = {}
    if master.empty:
        return src, hla
    for value, g in master.groupby("source_name", dropna=False):
        src[str(value)] = {
            "n": float(len(g)),
            "positive_prevalence": float(pd.to_numeric(g["label"], errors="coerce").mean()),
        }
    for value, g in master.groupby("hla_allele_4digit", dropna=False):
        hla[str(value)] = {
            "n": float(len(g)),
            "positive_prevalence": float(pd.to_numeric(g["label"], errors="coerce").mean()),
        }
    return src, hla


def bounded_conf_bin(confidence: float) -> str:
    if confidence >= 0.70:
        return "high"
    if confidence >= 0.45:
        return "medium"
    return "low"


def score_candidate(row: pd.Series, maps: dict[str, dict[str, float]], source_stats: dict[str, dict[str, float]], hla_stats: dict[str, dict[str, float]]) -> dict[str, Any]:
    source = str(row.get("source_name", ""))
    hla = str(row.get("hla_allele_4digit", ""))
    supertype = str(row.get("hla_supertype", ""))
    peptide_length = str(row.get("peptide_length", ""))
    leakage = str(row.get("leakage_risk_level", "")).lower()
    source_prev = source_stats.get(source, {}).get("positive_prevalence", math.nan)
    source_n = source_stats.get(source, {}).get("n", math.nan)
    hla_support = hla_stats.get(hla, {}).get("n", math.nan)
    hla_prev = hla_stats.get(hla, {}).get("positive_prevalence", math.nan)

    base_score = as_float(row.get("patient_gated_bma_score", row.get("barneo_bma_score", 0.0)), 0.0)
    base_conf = as_float(row.get("bma_confidence_score", 0.0), 0.0)
    selected_count = as_float(row.get("selected_method_count", 0), 0.0)

    hla_fp = max(maps.get("hla_allele_4digit:fp", {}).get(hla, 0.0), maps.get("hla_supertype:fp", {}).get(supertype, 0.0))
    hla_fn = max(maps.get("hla_allele_4digit:fn", {}).get(hla, 0.0), maps.get("hla_supertype:fn", {}).get(supertype, 0.0))
    source_fp = maps.get("source_name:fp", {}).get(source, 0.0)
    length_fp = maps.get("peptide_length:fp", {}).get(peptide_length, 0.0)

    score_multiplier = 1.0
    confidence_cap = 1.0
    risk_score = 0.0
    reasons: list[str] = []

    if leakage == "high":
        score_multiplier *= 0.82
        confidence_cap = min(confidence_cap, 0.40)
        risk_score += 0.25
        reasons.append("High leakage-risk region; clean benchmark claim unsafe")

    if math.isfinite(source_prev) and source_prev < LOW_PREVALENCE_CUTOFF:
        score_multiplier *= 0.86
        confidence_cap = min(confidence_cap, 0.45)
        risk_score += 0.18
        reasons.append("Low-prevalence source; false-positive pressure expected")

    if math.isfinite(hla_support) and hla_support < RARE_HLA_SUPPORT:
        score_multiplier *= 0.80
        confidence_cap = min(confidence_cap, 0.35)
        risk_score += 0.22
        reasons.append("Rare HLA allele support in benchmark")
    elif math.isfinite(hla_support) and hla_support < LOW_HLA_SUPPORT:
        score_multiplier *= 0.92
        confidence_cap = min(confidence_cap, 0.55)
        risk_score += 0.10
        reasons.append("Low HLA allele support in benchmark")

    if hla_fp >= 8:
        score_multiplier *= 0.88
        confidence_cap = min(confidence_cap, 0.45)
        risk_score += min(0.20, 0.02 * math.log1p(hla_fp))
        reasons.append("HLA group enriched among high-ranked negatives")
    if source_fp >= 8:
        score_multiplier *= 0.90
        confidence_cap = min(confidence_cap, 0.50)
        risk_score += min(0.15, 0.02 * math.log1p(source_fp))
        reasons.append("Source enriched among high-ranked negatives")
    if length_fp >= 8:
        score_multiplier *= 0.94
        risk_score += 0.05
        reasons.append("Peptide-length group enriched among high-ranked negatives")

    missed_positive_watch = hla_fn >= 8
    rescue_bonus = 0.0
    if missed_positive_watch:
        confidence_cap = min(confidence_cap, 0.55)
        risk_score += min(0.12, 0.015 * math.log1p(hla_fn))
        reasons.append("HLA group enriched among missed positives; manual low-score review warranted")
        if base_score < 0.25:
            rescue_bonus = min(0.05, 0.01 * math.log1p(hla_fn))

    if selected_count < 3:
        score_multiplier *= 0.90
        confidence_cap = min(confidence_cap, 0.40)
        risk_score += 0.12
        reasons.append("Sparse selected expert support")

    if truthy(row.get("bma_abstain", False)):
        confidence_cap = min(confidence_cap, 0.45)
        prior_reason = str(row.get("bma_abstention_reason_primary", "")).strip()
        if prior_reason:
            reasons.append(f"BMA abstention inherited: {prior_reason}")

    adjusted_score = min(1.0, max(0.0, base_score * score_multiplier + rescue_bonus))
    confidence = min(base_conf, confidence_cap)
    confidence = max(0.0, min(1.0, confidence - min(0.35, 0.30 * risk_score)))
    abstain = bool(reasons) or confidence < 0.45
    if not abstain and confidence >= 0.45:
        review_tier = "clean_review_prompt"
    elif adjusted_score >= 0.25 and base_score >= 0.35:
        review_tier = "manual_review_high_score_claim_blocked"
    elif missed_positive_watch and base_score < 0.25:
        review_tier = "manual_review_missed_positive_watchlist"
    elif adjusted_score >= 0.15:
        review_tier = "cautious_manual_review"
    else:
        review_tier = "deprioritize"
    return {
        "failure_aware_score": adjusted_score,
        "failure_aware_confidence": confidence,
        "failure_aware_confidence_bin": bounded_conf_bin(confidence),
        "failure_aware_abstain": abstain,
        "failure_aware_review_tier": review_tier,
        "failure_aware_reason_primary": reasons[0] if reasons else "",
        "failure_aware_reason_all": "; ".join(dict.fromkeys(reasons)),
        "base_patient_gated_bma_score": base_score,
        "base_bma_confidence": base_conf,
        "score_multiplier": score_multiplier,
        "rescue_bonus": rescue_bonus,
        "failure_risk_score": risk_score,
        "source_positive_prevalence": source_prev,
        "source_n": source_n,
        "hla_allele_support_count": hla_support,
        "hla_allele_positive_prevalence": hla_prev,
        "hla_false_positive_enrichment_score": hla_fp,
        "hla_missed_positive_enrichment_score": hla_fn,
        "source_false_positive_enrichment_score": source_fp,
        "peptide_length_false_positive_enrichment_score": length_fp,
    }


def apply_adjustment(output_root: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    master = load_tsv(output_root / "clean_neobench_master.tsv")
    bma = load_tsv(output_root / "barneo_bma_candidate_scores.tsv")
    failure_dist = load_tsv(output_root / "clean_neobench_failure_distribution_enrichment.tsv")
    if master.empty or bma.empty:
        return pd.DataFrame(), pd.DataFrame()

    source_stats, hla_stats = source_hla_stats(master)
    maps = risk_maps(failure_dist)
    merged = bma.merge(
        master[
            [
                "candidate_id",
                "label",
                "source_name",
                "hla_allele_4digit",
                "hla_supertype",
                "peptide",
                "peptide_length",
                "leakage_risk_level",
                "split_low_prevalence",
                "exact_peptide_hla_train_overlap",
                "near_peptide_train_overlap",
            ]
        ],
        on="candidate_id",
        how="left",
    )
    scored = pd.DataFrame([score_candidate(r, maps, source_stats, hla_stats) for _, r in merged.iterrows()], index=merged.index)
    out = pd.concat([merged, scored], axis=1)
    out["failure_aware_rank_global"] = out["failure_aware_score"].rank(ascending=False, method="first").astype(int)
    out = out.sort_values("failure_aware_rank_global")

    summary_rows = []
    for group_col in ["source_name", "hla_allele_4digit", "leakage_risk_level", "failure_aware_reason_primary", "failure_aware_review_tier"]:
        if group_col not in out.columns:
            continue
        for value, g in out.groupby(group_col, dropna=False):
            summary_rows.append(
                {
                    "group_type": group_col,
                    "group_value": value,
                    "n_candidates": int(len(g)),
                    "n_label_pos": int(pd.to_numeric(g["label"], errors="coerce").fillna(0).sum()) if "label" in g else 0,
                    "mean_base_score": float(g["base_patient_gated_bma_score"].mean()),
                    "mean_failure_aware_score": float(g["failure_aware_score"].mean()),
                    "mean_confidence": float(g["failure_aware_confidence"].mean()),
                    "abstention_rate": float(g["failure_aware_abstain"].mean()),
                }
            )
    summary = pd.DataFrame(summary_rows).sort_values(["group_type", "n_candidates"], ascending=[True, False])
    return out, summary


def md_table(df: pd.DataFrame, columns: list[str], n: int = 20) -> str:
    if df.empty:
        return "No rows available."
    cols = [c for c in columns if c in df.columns]
    if not cols:
        return "No requested columns available."
    try:
        return df[cols].head(n).to_markdown(index=False)
    except Exception:
        return df[cols].head(n).to_csv(sep="\t", index=False)


def write_report(output_root: Path, scores: pd.DataFrame, summary: pd.DataFrame) -> None:
    top_review = (
        scores[scores["failure_aware_review_tier"].isin(["clean_review_prompt", "manual_review_high_score_claim_blocked", "manual_review_missed_positive_watchlist", "cautious_manual_review"])].copy()
        if not scores.empty and "failure_aware_review_tier" in scores.columns
        else pd.DataFrame()
    )
    top_review = top_review.sort_values(["failure_aware_review_tier", "failure_aware_score"], ascending=[True, False]) if not top_review.empty else top_review
    biggest_down = scores.assign(delta=scores["failure_aware_score"] - scores["base_patient_gated_bma_score"]).sort_values("delta") if not scores.empty else pd.DataFrame()
    if not scores.empty:
        reason_counts = (
            scores["failure_aware_reason_primary"]
            .replace("", "no_failure_aware_reason")
            .value_counts()
            .rename_axis("reason")
            .reset_index(name="n")
        )
    else:
        reason_counts = pd.DataFrame()
    tier_counts = (
        scores["failure_aware_review_tier"].value_counts().reset_index().rename(columns={"index": "review_tier", "failure_aware_review_tier": "n"})
        if not scores.empty and "failure_aware_review_tier" in scores.columns
        else pd.DataFrame()
    )
    text = f"""# BAR-Neo Failure-Aware Reliability Adjustment

## Purpose

This layer applies the observed failure-mode analysis to BAR-Neo-BMA candidate scores. It does not train a new predictor. It reduces confidence and score where source prevalence, HLA support, leakage risk, or enriched false-positive groups suggest unstable behavior.

## Summary

- Candidates adjusted: {len(scores)}
- Non-abstain after failure-aware adjustment: {int((~scores["failure_aware_abstain"].astype(bool)).sum()) if len(scores) else 0}
- Abstain after failure-aware adjustment: {int(scores["failure_aware_abstain"].astype(bool).sum()) if len(scores) else 0}

## Non-Abstain Review Prompts

There may be zero clean non-abstain rows. The `failure_aware_review_tier` column separates clean claim eligibility from practical manual review priority.

{md_table(top_review, ["candidate_id", "label", "source_name", "hla_allele_4digit", "base_patient_gated_bma_score", "failure_aware_score", "failure_aware_confidence", "failure_aware_review_tier", "failure_aware_rank_global", "selected_method_count"], 40)}

## Biggest Down-Weighted Rows

{md_table(biggest_down, ["candidate_id", "label", "source_name", "hla_allele_4digit", "base_patient_gated_bma_score", "failure_aware_score", "delta", "failure_aware_reason_primary", "failure_aware_reason_all"], 30)}

## Primary Reason Counts

{md_table(reason_counts, ["reason", "n"], 30)}

## Review Tier Counts

{md_table(tier_counts, ["review_tier", "n"], 20)}

## Group Summary

{md_table(summary, ["group_type", "group_value", "n_candidates", "n_label_pos", "mean_base_score", "mean_failure_aware_score", "mean_confidence", "abstention_rate"], 60)}

## Claim Boundary

Failure-aware adjustment improves practical triage discipline. It is not clinical vaccine selection and not a new SOTA predictor claim.
"""
    (output_root / "BAR_NEO_FAILURE_AWARE_RELIABILITY_REPORT.md").write_text(text.rstrip() + "\n")


def main() -> None:
    args = parse_args()
    output_root = Path(args.output_root)
    ensure_dir(output_root)
    scores, summary = apply_adjustment(output_root)
    write_tsv(scores, output_root / "barneo_failure_aware_candidate_scores.tsv")
    write_tsv(summary, output_root / "barneo_failure_aware_group_summary.tsv")
    write_report(output_root, scores, summary)

    manifest_path = output_root / "run_manifest.json"
    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text())
    else:
        manifest = {}
    manifest.setdefault("output_files", [])
    for name in [
        "barneo_failure_aware_candidate_scores.tsv",
        "barneo_failure_aware_group_summary.tsv",
        "BAR_NEO_FAILURE_AWARE_RELIABILITY_REPORT.md",
    ]:
        if name not in manifest["output_files"]:
            manifest["output_files"].append(name)
    manifest.setdefault("summary", {})
    manifest["summary"].update(
        {
            "n_failure_aware_candidates": int(len(scores)),
            "n_failure_aware_nonabstain": int((~scores["failure_aware_abstain"].astype(bool)).sum()) if len(scores) else 0,
            "n_failure_aware_abstain": int(scores["failure_aware_abstain"].astype(bool).sum()) if len(scores) else 0,
        }
    )
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    update_manifest(
        output_root,
        "failure_aware_reliability_adjustment",
        {
            "n_candidates": int(len(scores)),
            "n_nonabstain": int((~scores["failure_aware_abstain"].astype(bool)).sum()) if len(scores) else 0,
            "warnings": [
                "Failure-aware adjustment is a conservative post-hoc triage layer, not a new predictor."
            ],
        },
    )
    print(f"[barneo-failure-aware] candidates={len(scores)} nonabstain={int((~scores['failure_aware_abstain'].astype(bool)).sum()) if len(scores) else 0}")


if __name__ == "__main__":
    main()
