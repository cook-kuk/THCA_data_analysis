#!/usr/bin/env python3
"""Build a CLEAN-NeoBench challenge pack for next experiments.

The pack converts benchmark failures into concrete next-review queues:
low-prevalence stress, rare-HLA rescue, high-score claim-blocked rows,
external/holdout fragility, public-vs-internal disagreement, Korean-HLA focus,
and PAAD/THCA metadata blockers.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from common import KOREAN_HLA_ALLELES, dataframe_to_markdown, ensure_dir, update_manifest, write_tsv


OUTPUTS = [
    "clean_neobench_challenge_pack.tsv",
    "clean_neobench_challenge_axis_summary.tsv",
    "clean_neobench_next_experiment_plan.tsv",
    "CLEAN_NEOBENCH_CHALLENGE_PACK.md",
    "CLEAN_NEOBENCH_CHALLENGE_PACK_KR.md",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", default=".", help="Repository root")
    parser.add_argument("--output-root", required=True, help="CLEAN-NeoBench output directory")
    parser.add_argument("--per-axis-limit", type=int, default=120, help="Maximum rows per challenge axis")
    return parser.parse_args()


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t") if path.exists() else pd.DataFrame()


def num(value: Any, default: float = 0.0) -> float:
    try:
        x = float(value)
        return x if math.isfinite(x) else default
    except Exception:
        return default


def safe_bool(value: Any) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes", "y"}


def add_missing_columns(df: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    out = df.copy()
    for col in cols:
        if col not in out.columns:
            out[col] = np.nan
    return out


def base_table(output_root: Path) -> pd.DataFrame:
    master = read_tsv(output_root / "clean_neobench_master.tsv")
    contextual = read_tsv(output_root / "barneo_contextual_bma_candidate_scores.tsv")
    failure = read_tsv(output_root / "barneo_failure_aware_candidate_scores.tsv")
    manual = read_tsv(output_root / "barneo_manual_review_queue.tsv")
    xscore = read_tsv(output_root / "barneo_x_candidate_scores.tsv")
    patient = read_tsv(output_root / "patient_gated_clean_neo_candidate_queue.tsv")

    if master.empty:
        return pd.DataFrame()

    cols_master = [
        "candidate_id",
        "source_name",
        "source_dataset",
        "study_id",
        "patient_id",
        "cancer_type",
        "disease_context",
        "peptide",
        "peptide_length",
        "hla_allele_4digit",
        "hla_supertype",
        "label",
        "mhc_class",
        "leakage_risk_level",
        "exact_peptide_hla_train_overlap",
        "near_peptide_train_overlap",
        "split_low_prevalence",
        "split_korean_hla_focus",
    ]
    df = master[[c for c in cols_master if c in master.columns]].copy()

    keep_contextual = [
        "candidate_id",
        "context_labels",
        "contextual_bma_score",
        "clean_contextual_bma_score",
        "contextual_confidence_score",
        "contextual_confidence_bin",
        "contextual_abstain",
        "contextual_abstention_reason_primary",
        "contextual_abstention_reason_all",
        "contextual_clean_claim_allowed",
        "selected_contextual_methods",
        "selected_clean_contextual_methods",
        "caveated_public_expert_used",
        "source_positive_prevalence",
        "hla_allele_support_count",
        "contextual_bma_rank_global",
        "clean_contextual_bma_rank_global",
    ]
    keep_failure = [
        "candidate_id",
        "failure_aware_score",
        "failure_aware_confidence",
        "failure_aware_review_tier",
        "failure_aware_reason_primary",
        "failure_aware_reason_all",
        "failure_risk_score",
        "source_false_positive_enrichment_score",
        "hla_false_positive_enrichment_score",
        "hla_missed_positive_enrichment_score",
        "failure_aware_rank_global",
    ]
    keep_manual = [
        "candidate_id",
        "manual_review_priority_bin",
        "manual_review_priority_score",
        "review_action_required",
        "best_clean_internal_support",
        "best_caveated_public_support",
        "best_qk_support",
        "public_minus_clean_support",
    ]
    keep_x = [
        "candidate_id",
        "barneo_x_discovery_score",
        "barneo_x_claim_safe_score",
        "barneo_x_discovery_rank",
        "barneo_x_claim_safe_rank",
        "barneo_x_primary_action",
    ]
    for src, keep in [
        (contextual, keep_contextual),
        (failure, keep_failure),
        (manual, keep_manual),
        (xscore, keep_x),
    ]:
        if not src.empty:
            df = df.merge(src[[c for c in keep if c in src.columns]].drop_duplicates("candidate_id"), on="candidate_id", how="left")

    if not patient.empty:
        patient_status = (
            patient.groupby("candidate_id")
            .agg(
                patient_gate_demo_rows=("scenario_id", "size"),
                max_demo_patient_gated_score=("demo_patient_gated_score", "max"),
                patient_metadata_status=("patient_metadata_status", "first"),
                patient_gate_abstention_reason=("abstention_reason_primary", "first"),
            )
            .reset_index()
        )
        df = df.merge(patient_status, on="candidate_id", how="left")

    df["source_positive_prevalence"] = pd.to_numeric(df.get("source_positive_prevalence"), errors="coerce")
    if df["source_positive_prevalence"].isna().all():
        df["source_positive_prevalence"] = df["source_name"].map(df.groupby("source_name")["label"].mean())
    df["hla_allele_support_count"] = pd.to_numeric(df.get("hla_allele_support_count"), errors="coerce")
    if df["hla_allele_support_count"].isna().all():
        df["hla_allele_support_count"] = df["hla_allele_4digit"].map(df.groupby("hla_allele_4digit")["candidate_id"].count())
    return df


def challenge_row(row: pd.Series, axis: str, priority_score: float, action: str, split_contract: str, success_metric: str, reason: str) -> dict[str, Any]:
    return {
        "challenge_axis": axis,
        "challenge_priority_score": float(priority_score),
        "candidate_id": row.get("candidate_id", ""),
        "label": row.get("label", np.nan),
        "source_name": row.get("source_name", ""),
        "source_positive_prevalence": row.get("source_positive_prevalence", np.nan),
        "hla_allele_4digit": row.get("hla_allele_4digit", ""),
        "hla_supertype": row.get("hla_supertype", ""),
        "hla_allele_support_count": row.get("hla_allele_support_count", np.nan),
        "peptide": row.get("peptide", ""),
        "peptide_length": row.get("peptide_length", np.nan),
        "mhc_class": row.get("mhc_class", ""),
        "leakage_risk_level": row.get("leakage_risk_level", ""),
        "split_low_prevalence": row.get("split_low_prevalence", ""),
        "contextual_bma_score": row.get("contextual_bma_score", np.nan),
        "clean_contextual_bma_score": row.get("clean_contextual_bma_score", np.nan),
        "contextual_confidence_score": row.get("contextual_confidence_score", np.nan),
        "contextual_abstain": row.get("contextual_abstain", np.nan),
        "contextual_abstention_reason_primary": row.get("contextual_abstention_reason_primary", ""),
        "failure_aware_score": row.get("failure_aware_score", np.nan),
        "failure_aware_review_tier": row.get("failure_aware_review_tier", ""),
        "manual_review_priority_bin": row.get("manual_review_priority_bin", ""),
        "manual_review_priority_score": row.get("manual_review_priority_score", np.nan),
        "barneo_x_discovery_score": row.get("barneo_x_discovery_score", np.nan),
        "barneo_x_claim_safe_score": row.get("barneo_x_claim_safe_score", np.nan),
        "barneo_x_primary_action": row.get("barneo_x_primary_action", ""),
        "best_clean_internal_support": row.get("best_clean_internal_support", np.nan),
        "best_caveated_public_support": row.get("best_caveated_public_support", np.nan),
        "public_minus_clean_support": row.get("public_minus_clean_support", np.nan),
        "selected_contextual_methods": row.get("selected_contextual_methods", ""),
        "recommended_next_action": action,
        "recommended_split_contract": split_contract,
        "success_metric": success_metric,
        "challenge_reason": reason,
        "claim_boundary": "research benchmark/manual review only; not clinical vaccine selection",
    }


def top_rows(df: pd.DataFrame, mask: pd.Series, score_cols: list[str], n: int) -> pd.DataFrame:
    sub = df[mask].copy()
    if sub.empty:
        return sub
    score = pd.Series(0.0, index=sub.index)
    for i, col in enumerate(score_cols):
        if col in sub.columns:
            vals = pd.to_numeric(sub[col], errors="coerce").fillna(0)
            score += vals * (1.0 / (i + 1))
    sub["_axis_score"] = score
    return sub.sort_values("_axis_score", ascending=False).head(n)


def build_pack(df: pd.DataFrame, per_axis_limit: int) -> pd.DataFrame:
    if df.empty:
        return pd.DataFrame()
    d = add_missing_columns(
        df,
        [
            "contextual_bma_score",
            "clean_contextual_bma_score",
            "contextual_confidence_score",
            "failure_aware_score",
            "manual_review_priority_score",
            "barneo_x_discovery_score",
            "barneo_x_claim_safe_score",
            "source_positive_prevalence",
            "hla_allele_support_count",
            "public_minus_clean_support",
        ],
    )
    label = pd.to_numeric(d["label"], errors="coerce")
    contextual = pd.to_numeric(d["contextual_bma_score"], errors="coerce").fillna(0)
    clean_contextual = pd.to_numeric(d["clean_contextual_bma_score"], errors="coerce").fillna(0)
    failure = pd.to_numeric(d["failure_aware_score"], errors="coerce").fillna(0)
    manual = pd.to_numeric(d["manual_review_priority_score"], errors="coerce").fillna(0)
    discovery = pd.to_numeric(d["barneo_x_discovery_score"], errors="coerce").fillna(0)
    source_prev = pd.to_numeric(d["source_positive_prevalence"], errors="coerce")
    hla_support = pd.to_numeric(d["hla_allele_support_count"], errors="coerce")
    public_minus_clean = pd.to_numeric(d["public_minus_clean_support"], errors="coerce")
    leakage_high = d["leakage_risk_level"].fillna("").astype(str).str.lower().eq("high")

    rows: list[dict[str, Any]] = []
    axis_defs = [
        (
            "claim_safe_priority_review",
            d["barneo_x_primary_action"].fillna("").eq("priority_review_candidate"),
            0.55 * discovery + 0.45 * pd.to_numeric(d["barneo_x_claim_safe_score"], errors="coerce").fillna(0),
            "inspect as the only current claim-safe review prompt; verify no hidden overlap and metadata gaps",
            "exact_phla + source_heldout + HLA_heldout",
            "AUPRC, top10 precision, calibration ECE, zero unresolved overlap",
            "BAR-Neo-X claim-safe score crossed priority-review threshold",
        ),
        (
            "high_score_claim_blocked",
            (contextual >= 0.75) & leakage_high,
            contextual,
            "audit exact/near peptide-HLA leakage before any clean claim",
            "exact_phla_holdout + near_peptide_cluster_holdout",
            "drop in top10 precision after removing overlap-risk rows",
            "high contextual score but leakage-risk blocks claim",
        ),
        (
            "low_prevalence_false_positive_stress",
            ((source_prev < 0.10) | d["split_low_prevalence"].fillna("").astype(str).eq("low_prevalence")) & label.eq(0),
            contextual + failure,
            "stress-test top-k threshold on TESLA-like low-prevalence negative-heavy sources",
            "low_prevalence_heldout",
            "top10 precision, top20 precision, false-positive rate among negatives",
            "low-prevalence setting where high-ranked negatives are most damaging",
        ),
        (
            "missed_positive_rescue_watchlist",
            label.eq(1) & (
                d["failure_aware_review_tier"].fillna("").astype(str).str.contains("missed_positive", regex=False)
                | (contextual < 0.35)
                | (failure < 0.20)
            ),
            (1.0 - contextual.clip(0, 1)) + (1.0 - failure.clip(0, 1)) + label.fillna(0),
            "manual rescue review; do not downrank solely by current model score",
            "HLA_heldout + source_heldout",
            "recall@20, mean positive rank, missed-positive rate",
            "positive candidate is low-ranked or flagged as missed-positive watchlist",
        ),
        (
            "rare_hla_support_gap",
            (hla_support < 20) | d["failure_aware_reason_all"].fillna("").astype(str).str.contains("Rare HLA", regex=False),
            (1.0 / (hla_support.fillna(999) + 1.0)).clip(0, 1) + contextual,
            "require allele-specific calibration or additional support examples",
            "HLA_heldout + Korean_HLA_focus_if_applicable",
            "allele-specific AUPRC, top10 precision, calibration ECE",
            "rare or underrepresented HLA allele",
        ),
        (
            "korean_hla_focus_stress",
            d["hla_allele_4digit"].isin(KOREAN_HLA_ALLELES),
            clean_contextual + discovery,
            "review Korean-relevant HLA calibration and allele-specific ranks",
            "korean_hla_focus",
            "allele-level AUPRC, top10 precision, positive rank",
            "Korean-HLA focus allele",
        ),
        (
            "external_holdout_fragility",
            d["context_labels"].fillna("").astype(str).str.contains("external_or_holdout", regex=False),
            (1.0 - pd.to_numeric(d["contextual_confidence_score"], errors="coerce").fillna(0)) + contextual,
            "treat as external/holdout fragility probe; compare against source-heldout calibration",
            "source_heldout_external",
            "source-heldout AUPRC, top-k survival, confidence-risk curve",
            "external/holdout context with unstable confidence",
        ),
        (
            "public_internal_disagreement",
            public_minus_clean.abs() >= 0.25,
            public_minus_clean.abs().fillna(0) + contextual,
            "separate caveated-public support from clean internal support before feature use",
            "public_overlap_audit + clean_internal_only_ablation",
            "delta AUPRC/top10 between public-inclusive and clean-internal-only views",
            "large public-vs-internal support disagreement",
        ),
        (
            "patient_gate_metadata_blocker",
            d["patient_metadata_status"].fillna("").astype(str).str.contains("missing", regex=False)
            | d["cancer_type"].isna()
            | d["disease_context"].isna(),
            clean_contextual + pd.to_numeric(d["barneo_x_claim_safe_score"], errors="coerce").fillna(0),
            "collect disease timing, presentation, antigen, immune-context, and safety fields before PAAD/THCA confidence",
            "patient_gated_PAAD_THCA_demo",
            "metadata completion rate, gate pass/fail concordance, patient-gated rank stability",
            "patient-gated PAAD/THCA interpretation blocked by missing metadata",
        ),
    ]

    for axis, mask, priority, action, split, metric, reason in axis_defs:
        subset = d[mask.fillna(False)].copy()
        if subset.empty:
            continue
        subset["_priority"] = pd.to_numeric(priority.loc[subset.index] if isinstance(priority, pd.Series) else priority, errors="coerce").fillna(0)
        subset = subset.sort_values("_priority", ascending=False).head(per_axis_limit)
        for _, row in subset.iterrows():
            rows.append(challenge_row(row, axis, row["_priority"], action, split, metric, reason))
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out = out.sort_values(["challenge_axis", "challenge_priority_score"], ascending=[True, False])
    out["challenge_rank_within_axis"] = out.groupby("challenge_axis")["challenge_priority_score"].rank(ascending=False, method="first").astype(int)
    out["global_challenge_rank"] = out["challenge_priority_score"].rank(ascending=False, method="first").astype(int)
    return out.sort_values("global_challenge_rank")


def build_axis_summary(pack: pd.DataFrame) -> pd.DataFrame:
    if pack.empty:
        return pd.DataFrame()
    rows = []
    for axis, g in pack.groupby("challenge_axis"):
        label = pd.to_numeric(g["label"], errors="coerce")
        rows.append(
            {
                "challenge_axis": axis,
                "n_rows": int(len(g)),
                "n_unique_candidates": int(g["candidate_id"].nunique()),
                "n_pos": int(label.fillna(0).sum()),
                "positive_prevalence": float(label.mean()) if label.notna().any() else math.nan,
                "mean_priority_score": float(pd.to_numeric(g["challenge_priority_score"], errors="coerce").mean()),
                "mean_contextual_bma_score": float(pd.to_numeric(g["contextual_bma_score"], errors="coerce").mean()),
                "mean_clean_contextual_bma_score": float(pd.to_numeric(g["clean_contextual_bma_score"], errors="coerce").mean()),
                "mean_confidence": float(pd.to_numeric(g["contextual_confidence_score"], errors="coerce").mean()),
                "top_sources": ", ".join(g["source_name"].fillna("NA").astype(str).value_counts().head(5).index.tolist()),
                "top_hla": ", ".join(g["hla_allele_4digit"].fillna("NA").astype(str).value_counts().head(5).index.tolist()),
                "recommended_next_action": g["recommended_next_action"].iloc[0],
                "recommended_split_contract": g["recommended_split_contract"].iloc[0],
                "success_metric": g["success_metric"].iloc[0],
            }
        )
    return pd.DataFrame(rows).sort_values(["n_rows", "mean_priority_score"], ascending=[False, False])


def build_experiment_plan(summary: pd.DataFrame) -> pd.DataFrame:
    base = [
        {
            "priority": "P0",
            "experiment": "public_training_corpus_row_overlap_audit",
            "input_rows": "public tool corpus files",
            "why": "public pretrained tools remain caveated until row-level overlap is known",
            "success_metric": "zero unresolved public-training overlap for clean comparator status",
            "expected_output": "clean_neobench_public_training_row_overlap.tsv",
        },
        {
            "priority": "P0",
            "experiment": "exact_near_overlap_lockdown",
            "input_rows": "high_score_claim_blocked axis",
            "why": "high scores are currently claim-blocked by exact/near overlap risk",
            "success_metric": "top-k performance after exact and near peptide-HLA removal",
            "expected_output": "lockdown leaderboard with no exact/near leakage",
        },
        {
            "priority": "P1",
            "experiment": "low_prevalence_topk_stress",
            "input_rows": "low_prevalence_false_positive_stress axis",
            "why": "TESLA-like settings punish false positives more than aggregate AUPRC shows",
            "success_metric": "top10/top20 precision and false-positive pressure in low-prevalence sources",
            "expected_output": "low-prevalence heldout board",
        },
        {
            "priority": "P1",
            "experiment": "rare_hla_and_korean_hla_calibration",
            "input_rows": "rare_hla_support_gap and korean_hla_focus_stress axes",
            "why": "underrepresented alleles drive both missed positives and unstable confidence",
            "success_metric": "allele-specific calibration ECE, AUPRC, positive rank",
            "expected_output": "HLA allele support and Korean-HLA calibration board",
        },
        {
            "priority": "P1",
            "experiment": "contextual_bma_clean_internal_ablation",
            "input_rows": "public_internal_disagreement axis",
            "why": "public tool support cannot be a clean feature until overlap audit is complete",
            "success_metric": "delta between all-expert contextual BMA and clean-contextual BMA",
            "expected_output": "public-inclusive vs clean-internal ablation table",
        },
        {
            "priority": "P2",
            "experiment": "PAAD_THCA_patient_gate_live_demo",
            "input_rows": "patient_gate_metadata_blocker axis plus real patient metadata",
            "why": "patient-gated score cannot be interpreted without disease/presentation/safety fields",
            "success_metric": "metadata completion and gate-driven rank stability",
            "expected_output": "research triage demo with clinical_use=false",
        },
        {
            "priority": "P2",
            "experiment": "MHC_II_separate_benchmark",
            "input_rows": "MHC-II rows only",
            "why": "Class I and Class II must not be pooled as a single predictor",
            "success_metric": "separate Class-II contracts and leaderboard",
            "expected_output": "CLEAN-NeoBench-II prototype",
        },
    ]
    plan = pd.DataFrame(base)
    if not summary.empty:
        counts = summary.set_index("challenge_axis")["n_unique_candidates"].to_dict()
        plan["challenge_rows_available"] = plan["input_rows"].map(
            lambda x: "; ".join(f"{axis}={counts.get(axis, 0)}" for axis in counts if axis in str(x)) or "NA"
        )
    return plan


def write_reports(output_root: Path, pack: pd.DataFrame, summary: pd.DataFrame, plan: pd.DataFrame) -> None:
    top = pack.sort_values("global_challenge_rank").head(60) if not pack.empty else pack
    text = f"""# CLEAN-NeoBench Challenge Pack

## Purpose

This pack converts current benchmark failures into concrete next experiments and manual-review queues. It is not a clean validation set and not a clinical selection list.

## Axis Summary

{dataframe_to_markdown(summary, max_rows=30)}

## Top Challenge Rows

{dataframe_to_markdown(top[["global_challenge_rank", "challenge_axis", "candidate_id", "label", "source_name", "source_positive_prevalence", "hla_allele_4digit", "hla_allele_support_count", "contextual_bma_score", "clean_contextual_bma_score", "contextual_confidence_score", "challenge_reason", "recommended_next_action"]], max_rows=60) if not top.empty else "No challenge rows available."}

## Next Experiment Plan

{dataframe_to_markdown(plan, max_rows=20)}

## Claim Boundary

Allowed: challenge benchmark design, manual-review prioritization, distribution stress testing, reliability/abstention refinement.

Forbidden: clinical vaccine selection, new SOTA claim, public pretrained tools as clean baselines without training-corpus overlap audit, quantum advantage, and pooled Class I/Class II predictor claims.
"""
    (output_root / "CLEAN_NEOBENCH_CHALLENGE_PACK.md").write_text(text.strip() + "\n")

    kr = f"""# CLEAN-NeoBench Challenge Pack KR

## 한 줄 결론

이제 다음 할 일이 명확하다. aggregate leaderboard가 아니라 **실패하는 분포를 고정 challenge set으로 만들고**, 그 축에서 다시 검증해야 한다.

## Challenge axis 요약

{dataframe_to_markdown(summary[["challenge_axis", "n_unique_candidates", "positive_prevalence", "mean_contextual_bma_score", "mean_confidence", "recommended_split_contract"]].head(30) if not summary.empty else summary, max_rows=30)}

## 다음 실험

{dataframe_to_markdown(plan[["priority", "experiment", "why", "success_metric"]], max_rows=20)}

## 실용 판단

- high score라도 leakage high면 claim-blocked.
- low-prevalence source에서는 top-k false-positive stress가 우선.
- rare/Korean HLA는 allele-specific calibration이 필요.
- public pretrained agreement는 audit 전까지 caveated support다.
- PAAD/THCA patient gate는 실제 patient metadata 전까지 demo-only다.
"""
    (output_root / "CLEAN_NEOBENCH_CHALLENGE_PACK_KR.md").write_text(kr.strip() + "\n")


def main() -> None:
    args = parse_args()
    output_root = Path(args.output_root)
    ensure_dir(output_root)

    df = base_table(output_root)
    pack = build_pack(df, args.per_axis_limit)
    summary = build_axis_summary(pack)
    plan = build_experiment_plan(summary)

    write_tsv(pack, output_root / "clean_neobench_challenge_pack.tsv")
    write_tsv(summary, output_root / "clean_neobench_challenge_axis_summary.tsv")
    write_tsv(plan, output_root / "clean_neobench_next_experiment_plan.tsv")
    write_reports(output_root, pack, summary, plan)

    manifest_path = output_root / "run_manifest.json"
    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text())
    else:
        manifest = {}
    manifest.setdefault("output_files", [])
    for name in OUTPUTS:
        if name not in manifest["output_files"]:
            manifest["output_files"].append(name)
    manifest.setdefault("summary", {})
    manifest["summary"].update(
        {
            "n_challenge_pack_rows": int(len(pack)),
            "n_challenge_axes": int(pack["challenge_axis"].nunique()) if not pack.empty else 0,
            "n_next_experiment_plan_rows": int(len(plan)),
        }
    )
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    update_manifest(
        output_root,
        "clean_neobench_challenge_pack",
        {
            "outputs": OUTPUTS,
            "n_rows": int(len(pack)),
            "n_axes": int(pack["challenge_axis"].nunique()) if not pack.empty else 0,
            "warnings": [
                "Challenge pack is a stress-test/manual-review design artifact, not a clinical selection list."
            ],
        },
    )
    print(
        "[clean-neobench-challenge-pack] "
        f"rows={len(pack)} axes={pack['challenge_axis'].nunique() if not pack.empty else 0}"
    )


if __name__ == "__main__":
    main()
