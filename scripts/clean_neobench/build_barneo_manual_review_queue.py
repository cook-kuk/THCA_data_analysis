#!/usr/bin/env python3
"""Build BAR-Neo manual review queue from failure-aware candidate scores.

The queue is intentionally separated from clean benchmark claims. It identifies
rows worth human review even when leakage/source/HLA risk blocks confident
algorithmic promotion.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from common import ensure_dir, update_manifest, write_tsv


SUPPORT_METHODS = [
    "source_balanced_rf_train_prior_calibrated",
    "pu_weighted_rf_train_prior_calibrated",
    "source_balanced_plus_pu_rf_train_prior_calibrated",
    "Structure_LR",
    "W7A_full",
    "W7A_QK_only",
    "W7B_stacked",
    "ESM2_Bayesian",
    "BigMHC_IM",
    "MHCflurry",
    "PRIME",
    "NetMHCpan_4.1",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", default=".", help="Repository root")
    parser.add_argument("--output-root", required=True, help="CLEAN-NeoBench output directory")
    parser.add_argument("--max-per-tier", type=int, default=250, help="Maximum candidates retained per review tier")
    return parser.parse_args()


def load_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t") if path.exists() else pd.DataFrame()


def to_float(value: Any, default: float = 0.0) -> float:
    try:
        v = float(value)
        return v if np.isfinite(v) else default
    except Exception:
        return default


def method_support_table(scores: pd.DataFrame) -> pd.DataFrame:
    if scores.empty:
        return pd.DataFrame()
    s = scores[scores["method_name"].isin(SUPPORT_METHODS)].copy()
    if s.empty:
        return pd.DataFrame()
    s["score_calibrated"] = pd.to_numeric(s["score_calibrated"], errors="coerce")
    pivot = s.pivot_table(index="candidate_id", columns="method_name", values="score_calibrated", aggfunc="max")
    pivot = pivot.rename(columns={m: f"support_{m}" for m in pivot.columns}).reset_index()
    clean_methods = [
        "support_source_balanced_rf_train_prior_calibrated",
        "support_pu_weighted_rf_train_prior_calibrated",
        "support_source_balanced_plus_pu_rf_train_prior_calibrated",
        "support_Structure_LR",
        "support_W7A_full",
        "support_W7B_stacked",
    ]
    public_methods = ["support_BigMHC_IM", "support_MHCflurry", "support_PRIME", "support_NetMHCpan_4.1"]
    qk_methods = ["support_W7A_QK_only"]
    present_clean = [c for c in clean_methods if c in pivot.columns]
    present_public = [c for c in public_methods if c in pivot.columns]
    present_qk = [c for c in qk_methods if c in pivot.columns]
    pivot["best_clean_internal_support"] = pivot[present_clean].max(axis=1, skipna=True) if present_clean else np.nan
    pivot["best_caveated_public_support"] = pivot[present_public].max(axis=1, skipna=True) if present_public else np.nan
    pivot["best_qk_support"] = pivot[present_qk].max(axis=1, skipna=True) if present_qk else np.nan
    pivot["public_minus_clean_support"] = pivot["best_caveated_public_support"] - pivot["best_clean_internal_support"]
    return pivot


def action_flags(row: pd.Series) -> str:
    actions: list[str] = []
    tier = str(row.get("failure_aware_review_tier", ""))
    reason = str(row.get("failure_aware_reason_all", ""))
    if "High leakage-risk" in reason:
        actions.append("audit exact/near peptide-HLA overlap before any clean claim")
    if "Low-prevalence source" in reason:
        actions.append("use source-calibrated top-k threshold; inspect low-prevalence false-positive risk")
    if "Rare HLA" in reason or "Low HLA" in reason or "HLA group enriched" in reason:
        actions.append("require HLA-support review and allele-specific calibration")
    if "missed_positive_watchlist" in tier or "missed positives" in reason:
        actions.append("manual rescue review: do not downrank solely by model score")
    if to_float(row.get("best_caveated_public_support"), np.nan) > to_float(row.get("best_clean_internal_support"), np.nan) + 0.25:
        actions.append("public-strong/internal-weak disagreement; keep public score caveated")
    if to_float(row.get("selected_method_count"), 0.0) < 3:
        actions.append("sparse expert support; rerun with more predictors or metadata")
    if not actions:
        actions.append("standard manual triage")
    return "; ".join(dict.fromkeys(actions))


def priority_score(row: pd.Series) -> float:
    tier = str(row.get("failure_aware_review_tier", ""))
    tier_bonus = {
        "manual_review_high_score_claim_blocked": 0.35,
        "manual_review_missed_positive_watchlist": 0.24,
        "cautious_manual_review": 0.18,
        "clean_review_prompt": 0.50,
        "deprioritize": -0.10,
    }.get(tier, 0.0)
    score = to_float(row.get("failure_aware_score")) * 0.40
    score += to_float(row.get("base_patient_gated_bma_score")) * 0.25
    score += to_float(row.get("best_clean_internal_support")) * 0.18
    score += min(1.0, to_float(row.get("selected_method_count")) / 5.0) * 0.08
    score += to_float(row.get("failure_aware_confidence")) * 0.09
    return score + tier_bonus


def priority_bin(score: float, tier: str) -> str:
    if tier == "manual_review_high_score_claim_blocked" and score >= 0.60:
        return "A_claim_blocked_high_score"
    if tier == "manual_review_missed_positive_watchlist" and score >= 0.45:
        return "B_rescue_watchlist"
    if score >= 0.42:
        return "C_cautious_review"
    return "D_deprioritize"


def build_queue(output_root: Path, max_per_tier: int) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    failure = load_tsv(output_root / "barneo_failure_aware_candidate_scores.tsv")
    scores = load_tsv(output_root / "clean_neobench_method_scores.tsv")
    if failure.empty:
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
    support = method_support_table(scores)
    q = failure.merge(support, on="candidate_id", how="left") if not support.empty else failure.copy()
    q["review_action_required"] = q.apply(action_flags, axis=1)
    q["manual_review_priority_score"] = q.apply(priority_score, axis=1)
    q["manual_review_priority_bin"] = q.apply(lambda r: priority_bin(float(r["manual_review_priority_score"]), str(r["failure_aware_review_tier"])), axis=1)

    keep_tiers = [
        "manual_review_high_score_claim_blocked",
        "manual_review_missed_positive_watchlist",
        "cautious_manual_review",
        "clean_review_prompt",
    ]
    queue = q[q["failure_aware_review_tier"].isin(keep_tiers)].copy()
    queue = (
        queue.sort_values(["failure_aware_review_tier", "manual_review_priority_score"], ascending=[True, False])
        .groupby("failure_aware_review_tier", group_keys=False)
        .head(max_per_tier)
        .sort_values(["manual_review_priority_bin", "manual_review_priority_score"], ascending=[True, False])
    )
    front_cols = [
        "candidate_id",
        "manual_review_priority_bin",
        "manual_review_priority_score",
        "failure_aware_review_tier",
        "label",
        "source_name",
        "hla_allele_4digit",
        "peptide",
        "base_patient_gated_bma_score",
        "failure_aware_score",
        "failure_aware_confidence",
        "failure_aware_reason_primary",
        "review_action_required",
        "best_clean_internal_support",
        "best_caveated_public_support",
        "best_qk_support",
        "public_minus_clean_support",
        "selected_method_count",
        "selected_methods",
    ]
    rest = [c for c in queue.columns if c not in front_cols]
    queue = queue[[c for c in front_cols if c in queue.columns] + rest]

    summary_rows = []
    for cols in [
        ["manual_review_priority_bin"],
        ["failure_aware_review_tier"],
        ["manual_review_priority_bin", "source_name"],
        ["manual_review_priority_bin", "hla_allele_4digit"],
        ["failure_aware_reason_primary"],
    ]:
        if not all(c in queue.columns for c in cols):
            continue
        for keys, g in queue.groupby(cols, dropna=False):
            if not isinstance(keys, tuple):
                keys = (keys,)
            row = {col: val for col, val in zip(cols, keys)}
            row.update(
                {
                    "group_by": "+".join(cols),
                    "n_candidates": int(len(g)),
                    "n_label_pos": int(pd.to_numeric(g.get("label", 0), errors="coerce").fillna(0).sum()),
                    "mean_priority_score": float(g["manual_review_priority_score"].mean()),
                    "mean_failure_aware_score": float(g["failure_aware_score"].mean()),
                    "mean_best_clean_internal_support": float(g["best_clean_internal_support"].mean()) if "best_clean_internal_support" in g else np.nan,
                }
            )
            summary_rows.append(row)
    summary = pd.DataFrame(summary_rows)
    challenge = build_distribution_challenge_queue(q)
    return queue, summary, challenge


def build_balanced_queue(queue: pd.DataFrame, per_source: int = 35, per_hla: int = 12) -> pd.DataFrame:
    if queue.empty:
        return queue
    parts = []
    if "source_name" in queue.columns:
        parts.append(
            queue.sort_values("manual_review_priority_score", ascending=False)
            .groupby("source_name", group_keys=False)
            .head(per_source)
            .assign(balance_reason="source_balanced")
        )
    if "hla_allele_4digit" in queue.columns:
        parts.append(
            queue.sort_values("manual_review_priority_score", ascending=False)
            .groupby("hla_allele_4digit", group_keys=False)
            .head(per_hla)
            .assign(balance_reason="hla_balanced")
        )
    if parts:
        out = pd.concat(parts, ignore_index=True, sort=False)
        out = out.sort_values("manual_review_priority_score", ascending=False).drop_duplicates("candidate_id")
    else:
        out = queue.copy()
        out["balance_reason"] = "unbalanced"
    return out.sort_values(["manual_review_priority_bin", "manual_review_priority_score"], ascending=[True, False])


def build_distribution_challenge_queue(all_rows: pd.DataFrame, per_group: int = 45) -> pd.DataFrame:
    """Select rows that stress known distribution-failure axes."""
    if all_rows.empty:
        return all_rows
    df = all_rows.copy()
    for col in ["source_positive_prevalence", "hla_allele_support_count", "public_minus_clean_support"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    parts = []
    if "source_positive_prevalence" in df.columns:
        low_prev = df[df["source_positive_prevalence"].lt(0.10)]
        if not low_prev.empty:
            parts.append(
                low_prev.sort_values("base_patient_gated_bma_score", ascending=False)
                .groupby("source_name", group_keys=False)
                .head(per_group)
                .assign(challenge_axis="low_prevalence_source")
            )
    if "hla_allele_support_count" in df.columns:
        rare_hla = df[df["hla_allele_support_count"].lt(20)]
        if not rare_hla.empty:
            parts.append(
                rare_hla.sort_values("base_patient_gated_bma_score", ascending=False)
                .groupby("hla_allele_4digit", group_keys=False)
                .head(max(5, per_group // 3))
                .assign(challenge_axis="rare_hla_support")
            )
    if "public_minus_clean_support" in df.columns:
        public_disagree = df[df["public_minus_clean_support"].gt(0.25)]
        if not public_disagree.empty:
            parts.append(
                public_disagree.sort_values("public_minus_clean_support", ascending=False)
                .head(per_group)
                .assign(challenge_axis="public_strong_internal_weak")
            )
    missed = df[df["failure_aware_review_tier"].eq("manual_review_missed_positive_watchlist")]
    if not missed.empty:
        parts.append(
            missed.sort_values("manual_review_priority_score", ascending=False)
            .groupby("source_name", group_keys=False)
            .head(per_group)
            .assign(challenge_axis="missed_positive_watchlist")
        )
    if not parts:
        return pd.DataFrame()
    out = pd.concat(parts, ignore_index=True, sort=False)
    out = out.sort_values(["challenge_axis", "manual_review_priority_score"], ascending=[True, False]).drop_duplicates(["candidate_id", "challenge_axis"])
    front = [
        "candidate_id",
        "challenge_axis",
        "manual_review_priority_bin",
        "manual_review_priority_score",
        "failure_aware_review_tier",
        "label",
        "source_name",
        "source_positive_prevalence",
        "hla_allele_4digit",
        "hla_allele_support_count",
        "peptide",
        "base_patient_gated_bma_score",
        "failure_aware_score",
        "best_clean_internal_support",
        "best_caveated_public_support",
        "public_minus_clean_support",
        "review_action_required",
    ]
    return out[[c for c in front if c in out.columns] + [c for c in out.columns if c not in front]]


def md_table(df: pd.DataFrame, columns: list[str], n: int = 30) -> str:
    if df.empty:
        return "No rows available."
    cols = [c for c in columns if c in df.columns]
    try:
        return df[cols].head(n).to_markdown(index=False)
    except Exception:
        return df[cols].head(n).to_csv(sep="\t", index=False)


def count_summary(series: pd.Series, name: str) -> pd.DataFrame:
    counts = (
        series.fillna("NA")
        .replace("", "NA")
        .value_counts()
        .rename_axis(name)
        .reset_index(name="n")
    )
    return counts


def write_report(output_root: Path, queue: pd.DataFrame, balanced: pd.DataFrame, challenge: pd.DataFrame, summary: pd.DataFrame) -> None:
    bins = count_summary(queue["manual_review_priority_bin"], "priority_bin") if not queue.empty else pd.DataFrame()
    tiers = count_summary(queue["failure_aware_review_tier"], "review_tier") if not queue.empty else pd.DataFrame()
    text = f"""# BAR-Neo Manual Review Queue

## Purpose

This queue is for practical human triage. It is not a clean benchmark claim list. Candidates enter the queue because they are high scoring but claim-blocked, likely missed positives, or cautious review cases after source/HLA failure-aware adjustment.

## Queue Counts

- Queue candidates: {len(queue)}

### Priority Bins

{md_table(bins, ["priority_bin", "n"], 20)}

### Review Tiers

{md_table(tiers, ["review_tier", "n"], 20)}

## Top Queue Rows

{md_table(queue, ["candidate_id", "manual_review_priority_bin", "manual_review_priority_score", "failure_aware_review_tier", "label", "source_name", "hla_allele_4digit", "peptide", "base_patient_gated_bma_score", "failure_aware_score", "failure_aware_confidence", "failure_aware_reason_primary", "review_action_required"], 50)}

## Balanced Queue Rows

This view caps each source/HLA so low-prevalence and underrepresented groups are not buried by CEDAR-heavy high-score rows.

{md_table(balanced, ["candidate_id", "balance_reason", "manual_review_priority_bin", "manual_review_priority_score", "failure_aware_review_tier", "label", "source_name", "hla_allele_4digit", "peptide", "failure_aware_score", "failure_aware_reason_primary", "review_action_required"], 80)}

## Distribution Challenge Queue

This queue intentionally samples low-prevalence sources, rare HLA alleles, public-vs-internal disagreement, and missed-positive watchlist rows. It is for stress testing and data curation.

{md_table(challenge, ["candidate_id", "challenge_axis", "manual_review_priority_bin", "manual_review_priority_score", "failure_aware_review_tier", "label", "source_name", "source_positive_prevalence", "hla_allele_4digit", "hla_allele_support_count", "peptide", "base_patient_gated_bma_score", "failure_aware_score", "best_clean_internal_support", "best_caveated_public_support", "public_minus_clean_support", "review_action_required"], 100)}

## Group Summary

{md_table(summary, ["group_by", "manual_review_priority_bin", "failure_aware_review_tier", "source_name", "hla_allele_4digit", "failure_aware_reason_primary", "n_candidates", "n_label_pos", "mean_priority_score", "mean_failure_aware_score"], 80)}

## Interpretation

Priority `A_claim_blocked_high_score` means the candidate is algorithmically interesting but not claim-eligible until leakage/source/HLA issues are audited. `B_rescue_watchlist` means the failure analysis suggests this group is prone to missed positives. These rows should drive data curation and validation design before any predictor claim.
"""
    (output_root / "BAR_NEO_MANUAL_REVIEW_QUEUE.md").write_text(text.rstrip() + "\n")


def main() -> None:
    args = parse_args()
    output_root = Path(args.output_root)
    ensure_dir(output_root)
    queue, summary, challenge = build_queue(output_root, args.max_per_tier)
    balanced = build_balanced_queue(queue)
    write_tsv(queue, output_root / "barneo_manual_review_queue.tsv")
    write_tsv(balanced, output_root / "barneo_manual_review_balanced_queue.tsv")
    write_tsv(challenge, output_root / "barneo_distribution_challenge_queue.tsv")
    write_tsv(summary, output_root / "barneo_manual_review_queue_summary.tsv")
    write_report(output_root, queue, balanced, challenge, summary)

    manifest_path = output_root / "run_manifest.json"
    manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
    manifest.setdefault("output_files", [])
    for name in [
        "barneo_manual_review_queue.tsv",
        "barneo_manual_review_balanced_queue.tsv",
        "barneo_distribution_challenge_queue.tsv",
        "barneo_manual_review_queue_summary.tsv",
        "BAR_NEO_MANUAL_REVIEW_QUEUE.md",
    ]:
        if name not in manifest["output_files"]:
            manifest["output_files"].append(name)
    manifest.setdefault("summary", {})
    manifest["summary"].update(
        {
            "n_manual_review_queue": int(len(queue)),
            "n_manual_review_balanced_queue": int(len(balanced)),
            "n_distribution_challenge_queue": int(len(challenge)),
            "n_manual_review_priority_A": int(queue["manual_review_priority_bin"].eq("A_claim_blocked_high_score").sum()) if not queue.empty else 0,
            "n_manual_review_priority_B": int(queue["manual_review_priority_bin"].eq("B_rescue_watchlist").sum()) if not queue.empty else 0,
        }
    )
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    update_manifest(
        output_root,
        "barneo_manual_review_queue",
        {
            "n_queue": int(len(queue)),
            "priority_counts": queue["manual_review_priority_bin"].value_counts().to_dict() if not queue.empty else {},
            "warnings": ["Manual review queue is not a clean benchmark claim list."],
        },
    )
    print(f"[barneo-review-queue] rows={len(queue)}")


if __name__ == "__main__":
    main()
