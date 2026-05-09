#!/usr/bin/env python3
"""Build BAR-Neo-X why-not and rescue-lane audit tables."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


WHY_NOT_OUTPUTS = [
    "barneo_x_why_not_audit.tsv",
    "barneo_x_why_not_reason_summary.tsv",
    "barneo_x_primary_blocker_summary.tsv",
    "barneo_x_rescue_lane_summary.tsv",
    "BAR_NEO_X_WHY_NOT_REPORT.md",
    "BAR_NEO_X_WHY_NOT_SUMMARY.json",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-root", required=True, help="CLEAN-NeoBench BAR-Neo output directory")
    return parser.parse_args()


def load(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t") if path.exists() else pd.DataFrame()


def write_tsv(df: pd.DataFrame, path: Path) -> None:
    df.to_csv(path, sep="\t", index=False)


def truthy(value: Any) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes", "y"}


def empty(value: Any) -> bool:
    return str(value).strip().lower() in {"", "na", "nan", "none", "null", "n/a"}


def num(value: Any, default: float = 0.0) -> float:
    try:
        v = float(value)
        return v if np.isfinite(v) else default
    except Exception:
        return default


def merge_context(output_root: Path) -> pd.DataFrame:
    explain = load(output_root / "barneo_x_candidate_explanations.tsv")
    master = load(output_root / "clean_neobench_master.tsv")
    barneo = load(output_root / "barneo_candidate_scores.tsv")
    bma = load(output_root / "barneo_bma_candidate_scores.tsv")
    if explain.empty:
        return pd.DataFrame()

    master_cols = [
        "candidate_id",
        "source_name",
        "study_id",
        "patient_id",
        "cancer_type",
        "disease_context",
        "exact_peptide_train_overlap",
        "exact_peptide_hla_train_overlap",
        "near_peptide_train_overlap",
        "study_train_overlap",
        "patient_train_overlap",
        "public_tool_training_overlap_any",
        "split_low_prevalence",
        "split_korean_hla_focus",
    ]
    barneo_cols = ["candidate_id", "confidence_bin", "abstain", "abstention_reason_primary"]
    bma_cols = [
        "candidate_id",
        "posterior_uncertainty",
        "expert_disagreement_range",
        "bma_abstain",
        "bma_abstention_reason_primary",
        "selected_method_count",
        "selected_methods",
    ]
    out = explain.copy()
    if not master.empty:
        out = out.merge(master[[c for c in master_cols if c in master.columns]], on="candidate_id", how="left", suffixes=("", "_master"))
    if not barneo.empty:
        out = out.merge(barneo[[c for c in barneo_cols if c in barneo.columns]], on="candidate_id", how="left", suffixes=("", "_barneo"))
    if not bma.empty:
        out = out.merge(bma[[c for c in bma_cols if c in bma.columns]], on="candidate_id", how="left", suffixes=("", "_bma"))
    return out


def reasons_for(row: pd.Series) -> list[str]:
    reasons: list[str] = []
    leak = str(row.get("leakage_risk_level", "")).lower()
    if leak == "high":
        reasons.append("high_leakage_claim_blocker")
    elif leak == "medium":
        reasons.append("medium_leakage_caution")
    if truthy(row.get("exact_peptide_hla_train_overlap")):
        reasons.append("exact_peptide_hla_overlap_claim_blocker")
    elif truthy(row.get("exact_peptide_train_overlap")):
        reasons.append("exact_peptide_overlap_claim_caveat")
    if truthy(row.get("near_peptide_train_overlap")):
        reasons.append("near_peptide_overlap_caution")
    if truthy(row.get("study_train_overlap")):
        reasons.append("study_overlap_caution")
    if truthy(row.get("patient_train_overlap")):
        reasons.append("patient_overlap_claim_blocker")
    if truthy(row.get("public_tool_training_overlap_any")):
        reasons.append("public_training_overlap_caveat")
    if empty(row.get("cancer_type")) or empty(row.get("disease_context")) or num(row.get("barneo_x_metadata_cap"), 1.0) < 1.0:
        reasons.append("missing_patient_disease_metadata")
    if str(row.get("split_low_prevalence", "")) == "low_prevalence":
        reasons.append("low_prevalence_source")
    abstention_text = " ".join(
        str(row.get(c, ""))
        for c in ["abstention_reason_primary", "bma_abstention_reason_primary", "abstention_reason_primary_barneo", "bma_abstention_reason_primary_bma"]
    ).lower()
    if "underrepresented" in abstention_text or "hla allele" in abstention_text:
        reasons.append("underrepresented_hla_or_sparse_allele")
    if "sparse" in abstention_text or num(row.get("selected_method_count"), 99) < 3:
        reasons.append("sparse_expert_support")
    if "disagreement" in abstention_text or num(row.get("expert_disagreement_range"), 0) >= 0.30:
        reasons.append("expert_disagreement")
    if "uncertainty" in abstention_text or num(row.get("posterior_uncertainty"), 0) >= 0.15:
        reasons.append("posterior_uncertainty")
    if truthy(row.get("abstain")) or truthy(row.get("bma_abstain")):
        reasons.append("model_abstention")
    if num(row.get("barneo_x_claim_safe_score"), 0) < 0.30:
        reasons.append("low_claim_safe_score")
    if str(row.get("label", "")) in {"0", "0.0"}:
        reasons.append("benchmark_negative_label")
    return list(dict.fromkeys(reasons))


def primary_blocker(reasons: list[str], action: str) -> str:
    if action == "priority_review_candidate":
        return "none_priority_candidate"
    if any(r in reasons for r in ["high_leakage_claim_blocker", "exact_peptide_hla_overlap_claim_blocker", "patient_overlap_claim_blocker"]):
        return "claim_blocked_by_leakage_or_identity_overlap"
    if "missing_patient_disease_metadata" in reasons:
        return "patient_metadata_incomplete"
    if any(r in reasons for r in ["expert_disagreement", "posterior_uncertainty", "model_abstention"]):
        return "model_disagreement_or_uncertainty"
    if any(r in reasons for r in ["sparse_expert_support", "underrepresented_hla_or_sparse_allele", "low_prevalence_source"]):
        return "support_gap_or_sparse_context"
    if "benchmark_negative_label" in reasons:
        return "benchmark_negative_or_low_priority"
    return "score_below_priority_threshold"


def rescue_lane(reasons: list[str], blocker: str) -> tuple[str, str, str]:
    if blocker == "none_priority_candidate":
        return (
            "priority_review_now",
            "complete patient metadata and independent validation before translational claim",
            "reviewer-facing research triage only",
        )
    if blocker == "claim_blocked_by_leakage_or_identity_overlap":
        return (
            "clean_claim_not_rescuable_without_new_holdout_or_overlap_audit",
            "row-level training-corpus audit, de-duplicated split, independent external cohort",
            "do not use as clean top-k benchmark evidence",
        )
    if blocker == "patient_metadata_incomplete":
        return (
            "rescuable_by_patient_metadata_completion",
            "cancer type, disease context, tumor stage, expression, VAF/clonality, HLA-LOH, B2M, immune context",
            "can be review-prioritized after metadata cap is lifted",
        )
    if blocker == "model_disagreement_or_uncertainty":
        return (
            "rescuable_by_consensus_or_manual_biology_review",
            "additional orthogonal predictors, allele-specific calibration, manual peptide/HLA biology review",
            "do not overstate score until disagreement is resolved",
        )
    if blocker == "support_gap_or_sparse_context":
        return (
            "rescuable_by_more_support_for_sparse_context",
            "more allele/source-specific data, low-prevalence stress testing, rare-HLA calibration",
            "treat as hypothesis, not reviewer-facing claim",
        )
    return (
        "low_priority_or_negative_benchmark_row",
        "manual review only if biology or disease context is compelling",
        "not a current claim candidate",
    )


def build_audit(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, row in df.iterrows():
        reasons = reasons_for(row)
        action = str(row.get("barneo_x_primary_action", ""))
        blocker = primary_blocker(reasons, action)
        lane, required, boundary = rescue_lane(reasons, blocker)
        discovery = num(row.get("barneo_x_discovery_score"))
        claim_gate = num(row.get("barneo_x_claim_gate"), 1.0)
        metadata_cap = num(row.get("barneo_x_metadata_cap"), 1.0)
        score_if_metadata_complete = float(np.clip(discovery * claim_gate, 0, 1))
        score_if_leakage_gate_removed = float(np.clip(discovery * metadata_cap, 0, 1))
        rows.append(
            {
                "candidate_id": row.get("candidate_id"),
                "peptide": row.get("peptide"),
                "hla_allele_4digit": row.get("hla_allele_4digit"),
                "label": row.get("label"),
                "leakage_risk_level": row.get("leakage_risk_level"),
                "barneo_x_claim_safe_rank": row.get("barneo_x_claim_safe_rank"),
                "barneo_x_claim_safe_score": row.get("barneo_x_claim_safe_score"),
                "barneo_x_primary_action": action,
                "primary_blocker": blocker,
                "why_not_reasons": "; ".join(reasons) if reasons else "none",
                "rescue_lane": lane,
                "required_next_data": required,
                "claim_boundary": boundary,
                "score_if_metadata_complete": score_if_metadata_complete,
                "score_if_leakage_gate_removed_forbidden": score_if_leakage_gate_removed,
                "would_be_priority_if_metadata_complete": bool(score_if_metadata_complete >= 0.50 and "high_leakage_claim_blocker" not in reasons),
                "top_positive_factors": row.get("barneo_x_top_positive_factors"),
                "top_negative_factors": row.get("barneo_x_top_negative_factors"),
            }
        )
    return pd.DataFrame(rows).sort_values(["barneo_x_claim_safe_rank", "candidate_id"])


def reason_summary(audit: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for reason in sorted({r for raw in audit["why_not_reasons"].astype(str) for r in raw.split("; ") if r and r != "none"}):
        sub = audit[audit["why_not_reasons"].astype(str).str.contains(reason, regex=False)]
        labels = pd.to_numeric(sub.get("label"), errors="coerce")
        rows.append(
            {
                "reason": reason,
                "n_candidates": int(len(sub)),
                "n_positive_label": int(labels.eq(1).sum()),
                "median_claim_safe_score": float(pd.to_numeric(sub["barneo_x_claim_safe_score"], errors="coerce").median()),
                "top_rescue_lane": sub["rescue_lane"].mode().iloc[0] if not sub["rescue_lane"].mode().empty else "NA",
            }
        )
    return pd.DataFrame(rows).sort_values(["n_candidates", "reason"], ascending=[False, True])


def lane_summary(audit: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for lane, sub in audit.groupby("rescue_lane", dropna=False):
        labels = pd.to_numeric(sub.get("label"), errors="coerce")
        rows.append(
            {
                "rescue_lane": lane,
                "n_candidates": int(len(sub)),
                "n_positive_label": int(labels.eq(1).sum()),
                "median_claim_safe_score": float(pd.to_numeric(sub["barneo_x_claim_safe_score"], errors="coerce").median()),
                "example_required_next_data": sub["required_next_data"].iloc[0] if len(sub) else "NA",
                "claim_boundary": sub["claim_boundary"].iloc[0] if len(sub) else "NA",
            }
        )
    return pd.DataFrame(rows).sort_values(["n_candidates", "rescue_lane"], ascending=[False, True])


def blocker_summary(audit: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for blocker, sub in audit.groupby("primary_blocker", dropna=False):
        labels = pd.to_numeric(sub.get("label"), errors="coerce")
        rows.append(
            {
                "primary_blocker": blocker,
                "n_candidates": int(len(sub)),
                "n_positive_label": int(labels.eq(1).sum()),
                "median_claim_safe_score": float(pd.to_numeric(sub["barneo_x_claim_safe_score"], errors="coerce").median()),
                "n_would_be_priority_if_metadata_complete": int(sub["would_be_priority_if_metadata_complete"].sum()),
                "top_rescue_lane": sub["rescue_lane"].mode().iloc[0] if not sub["rescue_lane"].mode().empty else "NA",
                "example_required_next_data": sub["required_next_data"].iloc[0] if len(sub) else "NA",
            }
        )
    return pd.DataFrame(rows).sort_values(["n_candidates", "primary_blocker"], ascending=[False, True])


def headline_counts(audit: pd.DataFrame) -> dict[str, int]:
    return {
        "n_hard_claim_blocked": int(audit["primary_blocker"].eq("claim_blocked_by_leakage_or_identity_overlap").sum()),
        "n_metadata_rescuable": int(audit["rescue_lane"].eq("rescuable_by_patient_metadata_completion").sum()),
        "n_priority_now": int(audit["rescue_lane"].eq("priority_review_now").sum()),
        "n_would_be_priority_if_metadata_complete": int(audit["would_be_priority_if_metadata_complete"].sum()),
    }


def write_report(output_root: Path, audit: pd.DataFrame, reasons: pd.DataFrame, blockers: pd.DataFrame, lanes: pd.DataFrame) -> None:
    cols = [
        "candidate_id",
        "peptide",
        "hla_allele_4digit",
        "label",
        "leakage_risk_level",
        "barneo_x_claim_safe_score",
        "primary_blocker",
        "rescue_lane",
        "required_next_data",
    ]
    counts = headline_counts(audit)
    lines = [
        "# BAR-Neo-X Why-Not Audit",
        "",
        "## Position",
        "",
        "This audit explains why a candidate is priority, secondary, blocked, or metadata-required. It is a reviewer-facing triage explanation layer, not a clinical decision layer.",
        "",
        "## Headline Counts",
        "",
        f"- Hard claim-blocked by leakage or identity overlap: {counts['n_hard_claim_blocked']}",
        f"- Rescuable by patient/disease metadata completion: {counts['n_metadata_rescuable']}",
        f"- Priority review now: {counts['n_priority_now']}",
        f"- Would cross priority threshold if metadata were complete: {counts['n_would_be_priority_if_metadata_complete']}",
        "",
        "## Primary Blocker Summary",
        "",
        blockers.to_markdown(index=False) if not blockers.empty else "No primary-blocker rows.",
        "",
        "## Reason Summary",
        "",
        reasons.to_markdown(index=False) if not reasons.empty else "No reason rows.",
        "",
        "## Rescue Lane Summary",
        "",
        lanes.to_markdown(index=False) if not lanes.empty else "No rescue-lane rows.",
        "",
        "## Top Claim-Safe Rows With Why-Not Context",
        "",
        audit[[c for c in cols if c in audit.columns]].head(30).round(4).to_markdown(index=False),
        "",
        "## Claim Boundary",
        "",
        "High leakage, exact peptide-HLA overlap, or patient overlap rows are not clean benchmark claims without a new holdout or row-level overlap audit. Missing metadata rows can be rescued only by completing patient/disease/presentation context.",
    ]
    (output_root / "BAR_NEO_X_WHY_NOT_REPORT.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    args = parse_args()
    output_root = Path(args.output_root)
    df = merge_context(output_root)
    if df.empty:
        raise SystemExit("missing BAR-Neo-X explanation inputs")
    audit = build_audit(df)
    reasons = reason_summary(audit)
    blockers = blocker_summary(audit)
    lanes = lane_summary(audit)
    write_tsv(audit, output_root / "barneo_x_why_not_audit.tsv")
    write_tsv(reasons, output_root / "barneo_x_why_not_reason_summary.tsv")
    write_tsv(blockers, output_root / "barneo_x_primary_blocker_summary.tsv")
    write_tsv(lanes, output_root / "barneo_x_rescue_lane_summary.tsv")
    write_report(output_root, audit, reasons, blockers, lanes)
    counts = headline_counts(audit)
    summary = {
        "n_candidates": int(len(audit)),
        "n_reasons": int(len(reasons)),
        "n_primary_blockers": int(len(blockers)),
        "n_rescue_lanes": int(len(lanes)),
        **counts,
        "top_reason": reasons.iloc[0].to_dict() if not reasons.empty else {},
        "top_primary_blocker": blockers.iloc[0].to_dict() if not blockers.empty else {},
    }
    (output_root / "BAR_NEO_X_WHY_NOT_SUMMARY.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
