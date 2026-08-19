#!/usr/bin/env python3
"""Write a compact Korean current-status report for CLEAN-NeoBench/BAR-Neo."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import pandas as pd

from common import dataframe_to_markdown, update_manifest


OUTPUT = "CLEAN_NEOBENCH_CURRENT_STATUS_KR.md"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", default=".", help="Repository root")
    parser.add_argument("--output-root", required=True, help="CLEAN-NeoBench output root")
    return parser.parse_args()


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t") if path.exists() else pd.DataFrame()


def fmt(value: Any, digits: int = 3) -> str:
    try:
        if pd.isna(value):
            return "NA"
        if isinstance(value, bool):
            return str(value)
        x = float(value)
        if x.is_integer():
            return f"{int(x):,}"
        if abs(x) >= 1000:
            return f"{x:,.0f}"
        return f"{x:.{digits}f}"
    except Exception:
        return str(value)


def metric(summary: dict[str, Any], key: str) -> str:
    return fmt(summary.get(key, "NA"))


def subset_columns(df: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    if df.empty:
        return df
    keep = [col for col in columns if col in df.columns]
    return df[keep].copy() if keep else df.iloc[:, :0].copy()


def main() -> None:
    args = parse_args()
    repo_root = Path(args.repo_root).resolve()
    output_root = Path(args.output_root)
    if not output_root.is_absolute():
        output_root = repo_root / output_root

    manifest_path = output_root / "run_manifest.json"
    manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
    summary = manifest.get("summary", {}) if isinstance(manifest.get("summary", {}), dict) else {}

    winloss = read_tsv(output_root / "clean_neobench_winloss_method_summary.tsv")
    leaderboard = read_tsv(output_root / "clean_neobench_leaderboard.tsv")
    contextual = read_tsv(output_root / "barneo_contextual_bma_context_summary.tsv")
    challenge = read_tsv(output_root / "clean_neobench_challenge_axis_summary.tsv")
    stress = read_tsv(output_root / "clean_neobench_source_hla_stress_method_summary.tsv")
    stress_guarded = read_tsv(output_root / "barneo_stress_guarded_candidate_scores.tsv")
    stress_guarded_weights = read_tsv(output_root / "barneo_stress_guarded_method_weights.tsv")
    high_impact_leads = read_tsv(output_root / "barneo_high_impact_lead_candidates.tsv")
    high_impact_topk = read_tsv(output_root / "barneo_high_impact_topk_audit.tsv")
    reviewer_kill = read_tsv(output_root / "barneo_high_impact_reviewer_kill_audit.tsv")
    top_pass = read_tsv(output_root / "barneo_top_pass_reviewer_evidence.tsv")
    lab_handoff = read_tsv(output_root / "barneo_high_impact_lab_handoff_queue.tsv")
    t1_readiness = read_tsv(output_root / "barneo_t1_translational_readiness.tsv")
    t1_assay = read_tsv(output_root / "barneo_t1_assay_handoff.tsv")
    t1_impact_assay = read_tsv(output_root / "barneo_t1_assay_design_matrix.tsv")
    t1_intake = read_tsv(output_root / "barneo_t1_metadata_intake_template.tsv")
    nature_claims = read_tsv(output_root / "clean_neobench_nature_claim_ladder.tsv")
    nature_gaps = read_tsv(output_root / "clean_neobench_nature_validation_gap_table.tsv")
    nature_sprint = read_tsv(output_root / "clean_neobench_nature_30day_sprint.tsv")
    ng_scores = read_tsv(output_root / "barneo_ng_algorithm_scores.tsv")
    ng_t1 = read_tsv(output_root / "barneo_ng_t1_candidates.tsv")
    ng_summary = read_tsv(output_root / "barneo_ng_decision_summary.tsv")
    ga_rl_bridge = read_tsv(output_root / "ga_rl_barneo_candidate_scores.tsv")
    ga_rl_t1 = read_tsv(output_root / "ga_rl_barneo_t1_confirmed.tsv")
    ga_rl_t1_unique = read_tsv(output_root / "ga_rl_barneo_t1_unique_candidates.tsv")
    ga_rl_high_blocked = read_tsv(output_root / "ga_rl_barneo_high_ga_blocked_audit.tsv")
    ga_rl_failure = read_tsv(output_root / "ga_rl_barneo_failure_analysis.tsv")
    ga_rl_patient_gate = read_tsv(output_root / "ga_rl_patient_gate_priority.tsv")
    ga_rl_patient_gate_summary = read_tsv(output_root / "ga_rl_patient_gate_priority_summary.tsv")
    ga_rl_heatmap = read_tsv(output_root / "ga_rl_t1_failure_heatmap.tsv")
    ga_rl_neighbors = read_tsv(output_root / "ga_rl_t1_nearest_neighbors.tsv")
    ga_rl_distribution = read_tsv(output_root / "ga_rl_t1_distribution_summary.tsv")
    ga_rl_atlas = read_tsv(output_root / "ga_rl_failure_atlas.tsv")
    ga_rl_atlas_summary = read_tsv(output_root / "ga_rl_failure_atlas_summary.tsv")
    ga_rl_atlas_reason = read_tsv(output_root / "ga_rl_failure_atlas_reason_summary.tsv")
    ga_rl_rescue = read_tsv(output_root / "ga_rl_rescue_quadrant.tsv")
    ga_rl_rescue_summary = read_tsv(output_root / "ga_rl_rescue_quadrant_summary.tsv")
    ga_rl_rescue_reason = read_tsv(output_root / "ga_rl_rescue_quadrant_reason_summary.tsv")
    ga_rl_rescue_evidence = read_tsv(output_root / "ga_rl_rescue_evidence.tsv")
    ga_rl_rescue_evidence_summary = read_tsv(output_root / "ga_rl_rescue_evidence_summary.tsv")
    mhc_ii_sep = read_tsv(output_root / "mhc_ii_separate_benchmark.tsv")
    mhc_ii_sep_summary = read_tsv(output_root / "mhc_ii_separate_benchmark_summary.tsv")
    ga_rl_summary = read_tsv(output_root / "ga_rl_barneo_decision_summary.tsv")
    patient_req = read_tsv(output_root / "patient_gated_clean_neo_metadata_requirements.tsv")
    public_audit = read_tsv(output_root / "clean_neobench_public_tool_overlap_audit.tsv")
    public_intake = read_tsv(output_root / "clean_neobench_public_overlap_audit_intake.tsv")

    clean_winloss = winloss[
        winloss.get("method_role", pd.Series(dtype=str)).isin(["anchor", "internal_candidate"])
    ].copy() if not winloss.empty else winloss
    qk_winloss = winloss[
        winloss.get("method_role", pd.Series(dtype=str)).eq("bounded_fallback")
    ].copy() if not winloss.empty else winloss
    public_winloss = winloss[
        winloss.get("method_role", pd.Series(dtype=str)).eq("caveated_public_comparator")
    ].copy() if not winloss.empty else winloss

    blocked = patient_req[
        patient_req.get("metadata_status", pd.Series(dtype=str)).isin(["absent_from_master", "present_but_empty"])
    ].copy() if not patient_req.empty else patient_req
    public_unresolved = public_audit[
        public_audit.get("training_overlap_audit_status", pd.Series(dtype=str)).astype(str).str.contains("unresolved|missing", case=False, na=False)
    ].copy() if not public_audit.empty else public_audit

    text = f"""# CLEAN-NeoBench / BAR-Neo Current Status KR

## 한 줄 결론

현재 위치는 `new SOTA predictor`가 아니라 **leakage-aware AI neoantigen benchmark + benchmark-adaptive reliability/abstention controller**다. 실용적으로는 agent가 여러 expert/ensemble을 불러오고, BAR-Neo-BMA가 benchmark 성능과 실패 패턴으로 가중치를 조절한 뒤, claim-safe layer가 위험한 후보를 abstain/manual review로 보낸다.

## 현재 run 숫자

| 항목 | 값 |
|---|---:|
| candidates | {metric(summary, "n_candidates")} |
| methods | {metric(summary, "n_methods")} |
| BMA weighted methods | {metric(summary, "n_bma_methods_weighted")} |
| split metric rows | {metric(summary, "n_metric_rows")} |
| BAR-Neo-BMA candidates | {metric(summary, "n_bma_candidates")} |
| BAR-Neo-BMA abstain | {metric(summary, "n_bma_abstain")} |
| BAR-Neo-BMA non-abstain | {metric(summary, "n_bma_nonabstain")} |
| contextual clean claims allowed | {metric(summary, "n_contextual_bma_clean_claim_allowed")} |
| BAR-Neo-X priority review candidates | {metric(summary, "n_barneo_x_priority_review_candidates")} |
| BAR-Neo-X top claim-safe score | {metric(summary, "barneo_x_top_claim_safe_score")} |
| challenge rows | {metric(summary, "n_challenge_pack_rows")} |
| visual dashboard figures | {metric(summary, "n_visual_dashboard_figures")} |
| public clean comparators allowed | {metric(summary, "n_public_clean_comparators_allowed")} |
| stress-guarded candidates | {metric(summary, "n_stress_guarded_candidates")} |
| stress-guarded review queue | {metric(summary, "n_stress_guarded_review_queue")} |
| stress-guarded claim/manual-audit | {metric(summary, "n_stress_guarded_claim_safe_after_manual_audit")} |
| stress-guarded top score | {metric(summary, "stress_guarded_top_score")} |
| high-impact leads | {metric(summary, "n_high_impact_leads")} |
| high-impact claim-safe benchmark positives | {metric(summary, "n_high_impact_claim_safe_benchmark_positives")} |
| reviewer-kill audit rows | {metric(summary, "n_reviewer_kill_audit_rows")} |
| reviewer-kill pass | {metric(summary, "n_reviewer_kill_pass")} |
| reviewer-kill manual audit | {metric(summary, "n_reviewer_kill_manual_audit")} |
| reviewer-kill blocked | {metric(summary, "n_reviewer_kill_blocked")} |
| top-pass reviewer evidence candidates | {metric(summary, "n_top_pass_reviewer_evidence_candidates")} |
| high-impact lab handoff rows | {metric(summary, "n_high_impact_lab_handoff_rows")} |
| T1 translational readiness candidates | {metric(summary, "n_t1_translational_readiness_candidates")} |
| T1 pHLA assay-design ready | {metric(summary, "n_t1_phla_assay_design_ready")} |
| T1 translational review ready | {metric(summary, "n_t1_translational_review_ready")} |
| T1 metadata intake rows | {metric(summary, "n_t1_metadata_intake_rows")} |
| T1 assay design rows | {metric(summary, "n_t1_assay_design_rows")} |
| Nature claim ladder rows | {metric(summary, "n_nature_claim_ladder_rows")} |
| Nature validation gaps | {metric(summary, "n_nature_validation_gap_rows")} |
| BAR-Neo-NG scores | {metric(summary, "n_barneo_ng_scores")} |
| BAR-Neo-NG T1 candidates | {metric(summary, "n_barneo_ng_t1_candidates")} |
| BAR-Neo-NG decision rows | {metric(summary, "n_barneo_ng_decision_rows")} |
| BAR-Neo-NG top score | {metric(summary, "barneo_ng_top_score")} |
| GA/RL bridge rows | {metric(summary, "n_ga_rl_bridge_rows")} |
| GA/RL bridge T1 rows | {metric(summary, "n_ga_rl_bridge_t1_confirmed")} |
| GA/RL bridge unique candidates | {metric(summary, "n_ga_rl_bridge_unique_candidates")} |
| GA/RL bridge unique T1 candidates | {metric(summary, "n_ga_rl_bridge_t1_unique_candidates")} |
| GA/RL high-GA blocked audit | {metric(summary, "n_ga_rl_bridge_high_ga_blocked_audit")} |
| GA/RL failure analysis rows | {metric(summary, "n_ga_rl_bridge_failure_rows")} |
| GA/RL patient-gate rows | {metric(summary, "n_ga_rl_patient_gate_priority_rows")} |
| GA/RL patient-gate candidates | {metric(summary, "n_ga_rl_patient_gate_priority_candidates")} |
| GA/RL bridge decision rows | {metric(summary, "n_ga_rl_bridge_decision_rows")} |
| GA/RL bridge top priority score | {metric(summary, "ga_rl_bridge_top_priority_score")} |

## 우리 알고리즘 위치

- `Structure_LR`: honest local anchor.
- `W7B_stacked`, `W7A_full`, source-balanced/PU RF: 현재 clean internal 후보군.
- `BAR-Neo`: 후보 intrinsic + method rank + reliability + uncertainty + leakage/source/HLA context로 candidate reliability score와 abstention reason을 낸다.
- `BAR-Neo-BMA`: 실용 agent ensemble controller. expert utility, top-k, calibration, source collapse, role prior, public caveat penalty, QK bounded fallback penalty를 가중치에 반영한다.
- `BAR-Neo-X`: leakage/claim-safe reranking과 설명 layer. 지금 clean claim은 막고 priority review만 만든다.
- `BAR-Neo-NG`: 새 알고리즘 본체. stress-guarded score, reviewer kill-audit, pHLA assay readiness, public/fallback penalty, leakage/near-neighbor cap, translational metadata gate를 하나의 claim-gated 점수로 합친다.
- QK 계열: 성능 좋은 slice가 있어도 **bounded fallback/fusion component**다. headline claim이 아니다.

## BAR-Neo-NG: 우리 새 알고리즘 적용 결과

{dataframe_to_markdown(ng_summary if not ng_summary.empty else ng_summary, max_rows=20)}

{dataframe_to_markdown(ng_scores[["candidate_id", "barneo_ng_rank_global", "barneo_ng_score", "barneo_ng_decision", "source_name", "hla_allele_4digit", "peptide", "label", "leakage_risk_level", "reviewer_kill_disposition", "impact_readiness_tier", "barneo_ng_reason"]].head(20) if not ng_scores.empty else ng_scores, max_rows=20)}

{dataframe_to_markdown(ng_t1[["candidate_id", "barneo_ng_rank_global", "barneo_ng_score", "barneo_ng_decision", "source_name", "hla_allele_4digit", "peptide", "label", "barneo_ng_reason"]].head(10) if not ng_t1.empty else ng_t1, max_rows=10)}

판단: BAR-Neo-NG는 점수가 높은 후보를 그냥 올리지 않고, reviewer-kill pass + pHLA assay-ready + 낮은 public/fallback 의존 + leakage gate 통과 후보만 T1로 올린다. 나머지는 T2/manual-audit/blocked/abstain으로 분리한다.

## GA/RL로 찾은 후보: BAR-Neo-NG gate 통과 여부

{dataframe_to_markdown(ga_rl_summary if not ga_rl_summary.empty else ga_rl_summary, max_rows=20)}

{dataframe_to_markdown(ga_rl_t1_unique[["ga_rl_barneo_unique_rank", "candidate_id", "peptide", "hla_allele_4digit", "source_name", "max_ga_rl_score_any_track", "barneo_ng_score", "ga_rl_barneo_priority_score", "supporting_ga_rl_tracks", "handoff_tier", "assay_next_step"]].head(10) if not ga_rl_t1_unique.empty else ga_rl_t1_unique, max_rows=10)}

{dataframe_to_markdown(ga_rl_t1[["ga_rl_barneo_rank_global", "candidate_id", "ga_rl_track", "peptide", "hla_allele_4digit", "ga_rl_score", "barneo_ng_score", "ga_rl_barneo_priority_score", "ga_rl_barneo_decision", "ga_rl_barneo_reason"]].head(12) if not ga_rl_t1.empty else ga_rl_t1, max_rows=12)}

{dataframe_to_markdown(ga_rl_bridge[["ga_rl_barneo_rank_global", "candidate_id", "ga_rl_track", "peptide", "hla_allele_4digit", "ga_rl_score", "barneo_ng_score", "ga_rl_barneo_priority_score", "ga_rl_barneo_decision", "leakage_risk_level", "ga_rl_native_blockers"]].head(20) if not ga_rl_bridge.empty else ga_rl_bridge, max_rows=20)}

{dataframe_to_markdown(ga_rl_high_blocked[["candidate_id", "ga_rl_track", "peptide", "hla_allele_4digit", "ga_rl_score", "barneo_ng_score", "ga_rl_barneo_decision", "leakage_risk_level", "audit_priority", "audit_question"]].head(20) if not ga_rl_high_blocked.empty else ga_rl_high_blocked, max_rows=20)}

{dataframe_to_markdown(ga_rl_failure[["analysis_group", "candidate_id", "ga_rl_track", "peptide", "hla_allele_4digit", "source_name", "ga_rl_score", "barneo_ng_score", "ga_rl_barneo_decision", "analysis_reason", "analysis_overlap_status"]].head(20) if not ga_rl_failure.empty else ga_rl_failure, max_rows=20)}

## GA/RL patient gate

{dataframe_to_markdown(subset_columns(ga_rl_patient_gate_summary, ["patient_gate_unique_rank", "candidate_id", "peptide", "hla_allele_4digit", "source_name", "patient_gate_best_score", "patient_gate_best_weighted_priority", "patient_gate_best_scenario_name", "patient_gate_best_disease", "clinical_use"]).head(10) if not ga_rl_patient_gate_summary.empty else ga_rl_patient_gate_summary, max_rows=10)}

{dataframe_to_markdown(ga_rl_patient_gate[["patient_gate_priority_rank", "candidate_id", "scenario_id", "scenario_name", "disease", "research_priority", "demo_patient_gated_score", "patient_gate_weighted_priority", "patient_gated_rank_within_scenario", "patient_gated_confidence_bin", "hard_gate_status", "abstain", "abstention_reason_primary"]].head(20) if not ga_rl_patient_gate.empty else ga_rl_patient_gate, max_rows=20)}

## GA/RL failure heatmap

{dataframe_to_markdown(ga_rl_distribution[["candidate_id", "peptide", "hla_allele_4digit", "source_name", "label", "ga_rl_barneo_priority_score", "barneo_ng_score", "leakage_risk_level", "split_source_heldout", "split_hla_heldout", "source_positive_prevalence", "source_support_count", "hla_positive_prevalence", "hla_support_count", "top_nn_candidate_id", "top_nn_peptide", "top_nn_similarity", "top_nn_label"]].head(10) if not ga_rl_distribution.empty else ga_rl_distribution, max_rows=10)}

{dataframe_to_markdown(ga_rl_heatmap[["candidate_id", "axis", "status", "value", "detail"]].head(30) if not ga_rl_heatmap.empty else ga_rl_heatmap, max_rows=30)}

{dataframe_to_markdown(ga_rl_neighbors[["query_candidate_id", "candidate_id", "peptide", "source_name", "hla_allele_4digit", "label", "similarity", "same_source", "same_hla", "same_supertype"]].head(18) if not ga_rl_neighbors.empty else ga_rl_neighbors, max_rows=18)}

## GA/RL failure atlas

{dataframe_to_markdown(ga_rl_atlas_summary[["atlas_group", "n_candidates", "n_pos", "mean_priority_score", "mean_barneo_ng_score", "n_high_leakage", "mean_source_prev", "mean_hla_prev", "mean_top_nn_similarity"]].head(10) if not ga_rl_atlas_summary.empty else ga_rl_atlas_summary, max_rows=10)}

{dataframe_to_markdown(ga_rl_atlas_reason[["atlas_group", "atlas_decision", "audit_priority", "n_candidates", "mean_priority_score", "mean_barneo_ng_score", "mean_source_prev", "mean_hla_prev"]].head(20) if not ga_rl_atlas_reason.empty else ga_rl_atlas_reason, max_rows=20)}

{dataframe_to_markdown(ga_rl_atlas[["atlas_order", "atlas_group", "candidate_id", "candidate_tier", "atlas_decision", "peptide", "hla_allele_4digit", "source_name", "ga_rl_barneo_priority_score", "barneo_ng_score", "leakage_risk_level", "audit_priority", "audit_question", "top_nn_candidate_id", "top_nn_similarity"]].head(25) if not ga_rl_atlas.empty else ga_rl_atlas, max_rows=25)}

## GA/RL rescue quadrant

{dataframe_to_markdown(ga_rl_rescue_summary[["rescue_category", "n_candidates", "n_pos", "mean_priority_score", "mean_rescue_priority_score", "mean_top_nn_similarity"]].head(10) if not ga_rl_rescue_summary.empty else ga_rl_rescue_summary, max_rows=10)}

{dataframe_to_markdown(ga_rl_rescue_reason[["rescue_category", "rescue_reason_primary", "n_candidates", "mean_priority_score", "mean_top_nn_similarity"]].head(20) if not ga_rl_rescue_reason.empty else ga_rl_rescue_reason, max_rows=20)}

{dataframe_to_markdown(ga_rl_rescue[["rescue_order", "rescue_category", "candidate_id", "peptide", "hla_allele_4digit", "source_name", "label", "ga_rl_barneo_priority_score", "ga_rl_barneo_decision", "rescue_action", "rescue_reason_primary", "top_nn_candidate_id", "top_nn_similarity"]].head(25) if not ga_rl_rescue.empty else ga_rl_rescue, max_rows=25)}

## GA/RL rescue evidence

{dataframe_to_markdown(ga_rl_rescue_evidence_summary[["evidence_grade", "n_candidates", "mean_missing", "mean_priority", "mean_source_prev", "mean_hla_prev"]].head(10) if not ga_rl_rescue_evidence_summary.empty else ga_rl_rescue_evidence_summary, max_rows=10)}

{dataframe_to_markdown(ga_rl_rescue_evidence[["candidate_id", "peptide", "hla_allele_4digit", "source_name", "missing_metadata_count", "present_metadata_count", "evidence_grade", "missing_metadata", "present_metadata", "evidence_action"]].head(10) if not ga_rl_rescue_evidence.empty else ga_rl_rescue_evidence, max_rows=10)}

## MHC-II separate benchmark

{dataframe_to_markdown(mhc_ii_sep_summary[["group_by", "n_candidates", "mean_ga_rl_score", "mean_priority_score"]].head(20) if not mhc_ii_sep_summary.empty else mhc_ii_sep_summary, max_rows=20)}

{dataframe_to_markdown(mhc_ii_sep[["separate_benchmark_rank", "candidate_id", "ga_rl_track", "source_name", "hla_allele_4digit", "ga_rl_score", "ga_rl_barneo_priority_score", "ga_rl_barneo_decision", "locked_public", "do_not_train", "manual_QA_required"]].head(12) if not mhc_ii_sep.empty else mhc_ii_sep, max_rows=12)}

판단: `GA/RL_hit_confirmed_T1_by_BAR_Neo_NG`는 track 중복 포함 row 수이고, 실제 unique 후보는 현재 `KLMNIQQKL`, `LLVDLAEEL`, `MLGEQLFPL` 3개다. GA/RL raw champion/상위권 중 다수는 high leakage 또는 clean-claim blocker 때문에 audit/failure case로 보내는 게 맞다.

## 대박으로 바꾼 새 layer: BAR-Neo Stress-Guarded

{dataframe_to_markdown(stress_guarded[["candidate_id", "stress_guarded_final_review_score", "stress_guarded_claim_safe_score", "stress_guarded_discovery_score", "stress_guarded_confidence", "stress_guarded_action", "stress_guarded_reason_primary", "source_name", "hla_allele_4digit", "peptide", "label"]].head(20) if not stress_guarded.empty else stress_guarded, max_rows=20)}

## Stress-guarded method weights

{dataframe_to_markdown(stress_guarded_weights[["method_name", "method_role", "mean_AUPRC", "reviewer_safe_score", "source_heldout_median_delta_AUPRC", "hla_stress_median_delta_AUPRC", "stress_guarded_base_weight", "recommended_use"]].head(15) if not stress_guarded_weights.empty else stress_guarded_weights, max_rows=15)}

## High-impact lead evidence pack

{dataframe_to_markdown(high_impact_topk.head(12) if not high_impact_topk.empty else high_impact_topk, max_rows=12)}

{dataframe_to_markdown(high_impact_leads[["candidate_id", "manual_audit_priority", "stress_guarded_rank_global", "stress_guarded_final_review_score", "benchmark_label_status", "source_name", "hla_allele_4digit", "peptide", "public_qk_independence_note"]].head(12) if not high_impact_leads.empty else high_impact_leads, max_rows=12)}

## Reviewer kill audit: 진짜 밀 후보 / 차단 후보

{dataframe_to_markdown(reviewer_kill[["candidate_id", "manual_audit_priority", "stress_guarded_rank_global", "stress_guarded_final_review_score", "source_name", "hla_allele_4digit", "peptide", "leakage_risk_level", "exact_peptide_hla_train_overlap", "near_peptide_train_overlap", "nearest_neighbor_similarity", "source_negatives_at_or_above", "public_weight_fraction", "fallback_weight_fraction", "reviewer_kill_disposition", "reviewer_kill_reason"]].head(20) if not reviewer_kill.empty else reviewer_kill, max_rows=20)}

판단: 현재 automated kill-audit 통과 후보는 `passes_current_reviewer_kill_audit`만이다. `manual_audit_required`는 near-neighbor caveat 때문에 추가 overlap/source 확인이 필요하고, `blocked_from_clean_claim`은 clean benchmark claim에서 제외한다.

## Top-pass reviewer evidence cards

{dataframe_to_markdown(top_pass[["candidate_id", "stress_guarded_rank_global", "stress_guarded_final_review_score", "source_name", "hla_allele_4digit", "peptide", "exact_overlap_clean", "near_neighbor_clean", "decoy_pressure_clean", "public_dependency_clean", "fallback_dependency_clean", "top_method", "recommended_next_action"]].head(10) if not top_pass.empty else top_pass, max_rows=10)}

## High-impact lab handoff queue

{dataframe_to_markdown(lab_handoff[["candidate_id", "handoff_tier", "stress_guarded_rank_global", "source_name", "hla_allele_4digit", "peptide", "reviewer_kill_disposition", "handoff_action", "must_verify_before_claim"]].head(20) if not lab_handoff.empty else lab_handoff, max_rows=20)}

## T1 translational readiness: 임팩트 올리는 지점

{dataframe_to_markdown(t1_readiness[["candidate_id", "stress_guarded_rank_global", "stress_guarded_final_review_score", "peptide", "hla_allele_4digit", "benchmark_evidence_score", "peptide_hla_assay_design_ready", "translational_metadata_score", "impact_readiness_tier", "recommended_next_action"]].head(10) if not t1_readiness.empty else t1_readiness, max_rows=10)}

{dataframe_to_markdown(t1_assay[["candidate_id", "assay_handoff_tier", "peptide", "hla_allele_4digit", "anchor_proxy", "minimum_next_data_packet"]].head(10) if not t1_assay.empty else t1_assay, max_rows=10)}

판단: T1 3개는 pHLA assay-design 후보로는 올라갔다. patient/translational lead claim은 WT peptide, gene/mutation, expression/clonality, patient context, HLA-LOH/B2M/safety metadata가 붙기 전까지 막는다.

## T1 impact upgrade: intake / assay / go-no-go

{dataframe_to_markdown(t1_intake[["candidate_id", "peptide", "hla_allele_4digit", "current_tier", "wt_peptide", "gene", "mutation_id", "expression_tpm", "vaf", "patient_id", "validation_rule"]].head(10) if not t1_intake.empty else t1_intake, max_rows=10)}

{dataframe_to_markdown(t1_impact_assay[["candidate_id", "assay_priority", "assay_id", "assay_or_review", "required_input", "current_status"]].head(24) if not t1_impact_assay.empty else t1_impact_assay, max_rows=24)}

## Nature-grade evolution: claim ladder / validation gap

{dataframe_to_markdown(nature_claims[["claim_level", "claim", "current_status", "evidence_now", "required_to_upgrade", "forbidden_overclaim"]].head(10) if not nature_claims.empty else nature_claims, max_rows=10)}

{dataframe_to_markdown(nature_gaps[["gap_id", "blocks_claim", "current_status", "minimum_fix", "impact_if_fixed"]].head(10) if not nature_gaps.empty else nature_gaps, max_rows=10)}

{dataframe_to_markdown(nature_sprint[["priority", "workstream", "task", "upgrade_unlocked"]].head(10) if not nature_sprint.empty else nature_sprint, max_rows=10)}

## Structure_LR 대비 승패 핵심표

{dataframe_to_markdown(winloss[["method_name", "method_role", "method_disposition", "n_matched_anchor_splits", "n_wins_vs_anchor", "n_losses_vs_anchor", "n_ties_vs_anchor", "win_rate_vs_anchor", "median_delta_AUPRC_vs_anchor", "source_heldout_median_delta_AUPRC", "hla_heldout_median_delta_AUPRC"]].head(20) if not winloss.empty else winloss, max_rows=20)}

## Clean internal 후보만 보면

{dataframe_to_markdown(clean_winloss[["method_name", "n_wins_vs_anchor", "n_losses_vs_anchor", "win_rate_vs_anchor", "median_delta_AUPRC_vs_anchor", "source_heldout_median_delta_AUPRC", "hla_heldout_median_delta_AUPRC"]].head(15) if not clean_winloss.empty else clean_winloss, max_rows=15)}

## QK / fallback 계열 판단

{dataframe_to_markdown(qk_winloss[["method_name", "n_wins_vs_anchor", "n_losses_vs_anchor", "win_rate_vs_anchor", "median_delta_AUPRC_vs_anchor", "source_heldout_median_delta_AUPRC", "hla_heldout_median_delta_AUPRC"]].head(12) if not qk_winloss.empty else qk_winloss, max_rows=12)}

판단: `W7A_QK_only`는 여러 matched split에서 Structure_LR을 이기지만, 역할은 fallback/fusion이다. 보고서에서는 이 결과를 ranking headline이나 quantum claim으로 쓰지 않는다.

## Public pretrained tool 판단

{dataframe_to_markdown(public_winloss[["method_name", "n_wins_vs_anchor", "n_losses_vs_anchor", "win_rate_vs_anchor", "median_delta_AUPRC_vs_anchor", "method_disposition"]].head(12) if not public_winloss.empty else public_winloss, max_rows=12)}

{dataframe_to_markdown(public_unresolved[["method_name", "training_overlap_audit_status", "clean_comparator_allowed_after_audit", "reviewer_disposition"]].head(12) if not public_unresolved.empty else public_unresolved, max_rows=12)}

판단: public tool이 일부 split에서 좋아도 row-level training-corpus overlap audit 전까지는 caveated comparator다.

## Public overlap audit intake

{dataframe_to_markdown(public_intake[["public_tool", "expected_drop_file", "required_key", "current_clean_allowed", "current_disposition"]].head(12) if not public_intake.empty else public_intake, max_rows=12)}

## 어디서 틀리는지

- source shift: 일부 RF/PU/source-balanced 계열은 overall AUPRC가 높지만 source-heldout median delta가 크게 음수다.
- low prevalence: TESLA-like 저 prevalence slice에서 false-positive pressure가 커진다.
- rare / underrepresented HLA: HLA-heldout 또는 Korean-HLA focus에서 method별 collapse가 갈린다.
- public/internal disagreement: public pretrained predictor가 높고 clean internal support가 약하면 manual review 또는 abstain.
- patient gate: PAAD/THCA patient-level fields가 부족해서 현재 patient-gated demo는 clinical-use가 아니라 research triage only다.

## Source/HLA stress audit

{dataframe_to_markdown(stress[["method_name", "method_role", "mean_AUPRC", "source_heldout_median_delta_AUPRC", "hla_stress_median_delta_AUPRC", "korean_hla_median_delta_AUPRC", "stress_flags", "recommended_use"]].head(20) if not stress.empty else stress, max_rows=20)}

## Challenge / distribution audit

{dataframe_to_markdown(challenge[["challenge_axis", "n_unique_candidates", "positive_prevalence", "recommended_split_contract", "success_metric"]].head(20) if not challenge.empty else challenge, max_rows=20)}

## Patient-gated blocker

{dataframe_to_markdown(blocked[["field", "gate_group", "required_for", "metadata_status", "reviewer_safe_default_if_missing"]].head(25) if not blocked.empty else blocked, max_rows=25)}

## 바로 볼 파일

- HTML dossier: `project/papers_hub_2026_05_04/clean_neobench_barneo_dossier_2026_05_09.html`
- Visual dashboard: `project/papers_hub_2026_05_04/clean_neobench_visual_dashboard_2026_05_10.html`
- Win/loss report: `project/results/clean_neobench_barneo_2026_05_09/CLEAN_NEOBENCH_WINLOSS_REPORT_KR.md`
- Win/loss table: `project/results/clean_neobench_barneo_2026_05_09/clean_neobench_winloss_method_summary.tsv`
- Split winner board: `project/results/clean_neobench_barneo_2026_05_09/clean_neobench_split_winner_board.tsv`
- Stress-guarded candidate scores: `project/results/clean_neobench_barneo_2026_05_09/barneo_stress_guarded_candidate_scores.tsv`
- Stress-guarded report: `project/results/clean_neobench_barneo_2026_05_09/BAR_NEO_STRESS_GUARDED_REPORT_KR.md`
- High-impact candidate pack: `project/results/clean_neobench_barneo_2026_05_09/BAR_NEO_HIGH_IMPACT_CANDIDATE_PACK_KR.md`
- High-impact HTML: `project/papers_hub_2026_05_04/barneo_high_impact_candidate_pack_2026_05_10.html`
- Reviewer kill audit: `project/results/clean_neobench_barneo_2026_05_09/BAR_NEO_REVIEWER_KILL_AUDIT_KR.md`
- Reviewer kill audit table: `project/results/clean_neobench_barneo_2026_05_09/barneo_high_impact_reviewer_kill_audit.tsv`
- Top-pass reviewer evidence: `project/results/clean_neobench_barneo_2026_05_09/BAR_NEO_TOP_PASS_REVIEWER_EVIDENCE_PACK_KR.md`
- Top-pass reviewer evidence HTML: `project/papers_hub_2026_05_04/barneo_top_pass_reviewer_evidence_2026_05_10.html`
- High-impact lab handoff queue: `project/results/clean_neobench_barneo_2026_05_09/barneo_high_impact_lab_handoff_queue.tsv`
- T1 translational readiness pack: `project/results/clean_neobench_barneo_2026_05_09/BAR_NEO_T1_TRANSLATIONAL_READINESS_PACK_KR.md`
- T1 translational readiness HTML: `project/papers_hub_2026_05_04/barneo_t1_translational_readiness_2026_05_10.html`
- T1 impact upgrade plan: `project/results/clean_neobench_barneo_2026_05_09/BAR_NEO_T1_IMPACT_UPGRADE_PLAN_KR.md`
- T1 metadata intake template: `project/results/clean_neobench_barneo_2026_05_09/barneo_t1_metadata_intake_template.tsv`
- T1 assay design matrix: `project/results/clean_neobench_barneo_2026_05_09/barneo_t1_assay_design_matrix.tsv`
- BAR-Neo-NG algorithm card: `project/results/clean_neobench_barneo_2026_05_09/BAR_NEO_NG_ALGORITHM_CARD.md`
- BAR-Neo-NG Korean result: `project/results/clean_neobench_barneo_2026_05_09/BAR_NEO_NG_RESULTS_KR.md`
- BAR-Neo-NG scores: `project/results/clean_neobench_barneo_2026_05_09/barneo_ng_algorithm_scores.tsv`
- BAR-Neo-NG T1 candidates: `project/results/clean_neobench_barneo_2026_05_09/barneo_ng_t1_candidates.tsv`
- BAR-Neo-NG HTML: `project/papers_hub_2026_05_04/barneo_ng_algorithm_results_2026_05_10.html`
- GA/RL bridge result: `project/results/clean_neobench_barneo_2026_05_09/GA_RL_BARNEO_BRIDGE_RESULTS_KR.md`
- GA/RL bridge scores: `project/results/clean_neobench_barneo_2026_05_09/ga_rl_barneo_candidate_scores.tsv`
- GA/RL bridge T1: `project/results/clean_neobench_barneo_2026_05_09/ga_rl_barneo_t1_confirmed.tsv`
- GA/RL bridge unique T1 handoff: `project/results/clean_neobench_barneo_2026_05_09/ga_rl_barneo_t1_unique_candidates.tsv`
- GA/RL high-GA blocked audit: `project/results/clean_neobench_barneo_2026_05_09/ga_rl_barneo_high_ga_blocked_audit.tsv`
- GA/RL failure analysis: `project/results/clean_neobench_barneo_2026_05_09/ga_rl_barneo_failure_analysis.tsv`
- GA/RL failure analysis KR: `project/results/clean_neobench_barneo_2026_05_09/GA_RL_BARNEO_FAILURE_ANALYSIS_KR.md`
- GA/RL patient gate: `project/results/clean_neobench_barneo_2026_05_09/ga_rl_patient_gate_priority.tsv`
- GA/RL patient gate summary: `project/results/clean_neobench_barneo_2026_05_09/ga_rl_patient_gate_priority_summary.tsv`
- GA/RL patient gate KR: `project/results/clean_neobench_barneo_2026_05_09/GA_RL_BARNEO_PATIENT_GATE_KR.md`
- GA/RL unique T1 handoff report: `project/results/clean_neobench_barneo_2026_05_09/GA_RL_BARNEO_UNIQUE_T1_HANDOFF_KR.md`
- GA/RL bridge HTML: `project/papers_hub_2026_05_04/ga_rl_barneo_bridge_results_2026_05_11.html`
- Nature-grade evolution plan: `project/results/clean_neobench_barneo_2026_05_09/BAR_NEO_NATURE_GRADE_EVOLUTION_PLAN_KR.md`
- Nature-grade evolution HTML: `project/papers_hub_2026_05_04/barneo_nature_grade_evolution_2026_05_10.html`
- Source/HLA stress audit: `project/results/clean_neobench_barneo_2026_05_09/CLEAN_NEOBENCH_SOURCE_HLA_STRESS_AUDIT_KR.md`
- Korean-HLA method board: `project/results/clean_neobench_barneo_2026_05_09/clean_neobench_korean_hla_method_board.tsv`
- Public overlap audit intake: `project/results/clean_neobench_barneo_2026_05_09/CLEAN_NEOBENCH_PUBLIC_OVERLAP_AUDIT_INTAKE_KR.md`
- Distribution error audit: `project/results/clean_neobench_barneo_2026_05_09/CLEAN_NEOBENCH_DISTRIBUTION_ERROR_AUDIT_KR.md`
- Challenge pack: `project/results/clean_neobench_barneo_2026_05_09/CLEAN_NEOBENCH_CHALLENGE_PACK_KR.md`

## 다음 우선순위

1. Public tool training-corpus row-level overlap audit.
2. Source-heldout collapse가 큰 RF/PU/source-balanced 모델의 source reweighting 또는 abstention 강화.
3. Rare/Korean-HLA focus split에서 W7B/W7A/source-balanced RF를 다시 stress-test.
4. PAAD/THCA patient-gated demo에 실제 disease timing, expression/clonality, immune context, safety metadata 연결.
5. MHC-II는 별도 benchmark로 분리.

## Claim boundary

Allowed: leakage-aware AI neoantigen predictor benchmarking framework, benchmark-adaptive reliability ranking, calibrated candidate prioritization with abstention, reviewer-safe comparison, research triage framework.

Forbidden: clinical vaccine selection, new SOTA predictor, external validation proven, QK headline claim, public tools as clean baselines without overlap audit, Class I and Class II unified predictor.
"""

    path = output_root / OUTPUT
    path.write_text(text.strip() + "\n")
    update_manifest(
        output_root,
        "clean_neobench_current_status_kr",
        {
            "outputs": [OUTPUT],
            "warnings": ["Current status report is a reviewer-safe project snapshot, not a SOTA or clinical claim."],
        },
    )
    print(f"[clean-neobench-status-kr] wrote {path}")


if __name__ == "__main__":
    main()
