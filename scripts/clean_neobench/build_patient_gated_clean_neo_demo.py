#!/usr/bin/env python3
"""Build a reviewer-safe PAAD/THCA patient-gated CLEAN-Neo demo layer.

The current benchmark candidates usually do not carry disease timing,
presentation integrity, immune-context, or safety metadata. This script
therefore writes a transparent demo and requirements audit rather than
pretending to make patient-level clinical decisions.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from common import dataframe_to_markdown, ensure_dir, normalize_empty, update_manifest, write_tsv


REQUIREMENTS = [
    ("patient_id", "identity", "both", "links candidates to patient-level gates", "unknown patient context"),
    ("cancer_type", "disease_context", "both", "selects PAAD vs THCA gate rules", "no disease-specific confidence"),
    ("disease_context", "disease_context", "both", "captures MRD, recurrence, metastatic, ATC, PDTC, RR-DTC", "disease gate unknown"),
    ("tumor_burden_context", "disease_context", "PAAD", "separates resected/MRD low-burden from bulky metastatic disease", "PAAD vaccine-first confidence capped"),
    ("disease_status", "disease_context", "both", "captures resected, MRD-positive, progressive, recurrent, metastatic", "disease gate unknown"),
    ("tumor_stage", "disease_context", "both", "supports high-risk recurrence/progression triage", "stage-based gate unknown"),
    ("histology", "disease_context", "THCA", "separates PTC, PDTC, ATC, RR-DTC", "routine low-risk PTC cannot be separated"),
    ("driver_mutation", "antigen", "both", "supports KRAS/BRAF/RET/NTRK/RAS/TERT context", "driver-aware strategy unavailable"),
    ("kras_mutation", "antigen", "PAAD", "enables separate KRAS public-vaccine logic", "KRAS-public branch unavailable"),
    ("braf_v600e", "disease_context", "THCA", "enables ATC BRAF/MEK backbone logic", "BRAF-specific THCA strategy unavailable"),
    ("ret_ntrk_fusion", "disease_context", "THCA", "supports targeted backbone selection", "fusion-targeted backbone unknown"),
    ("tert_status", "disease_context", "THCA", "supports high-risk recurrence context", "risk context incomplete"),
    ("hla", "presentation", "both", "minimum HLA typing", "presentation gate weak"),
    ("hla_allele_4digit", "presentation", "both", "allele-specific calibration and Korean-HLA focus", "allele support unknown"),
    ("hla_loh", "presentation", "both", "detects antigen-presentation escape", "presentation escape unknown"),
    ("b2m_status", "presentation", "both", "detects MHC-I presentation failure", "presentation escape unknown"),
    ("hla_expression", "presentation", "both", "supports intact HLA expression", "presentation gate unknown"),
    ("antigen_processing_status", "presentation", "both", "TAP/proteasome/presentation integrity", "processing gate unknown"),
    ("expression_tpm", "antigen", "both", "supports expressed antigen evidence", "antigen evidence weak"),
    ("mutant_expression", "antigen", "both", "supports expressed mutant allele", "antigen evidence weak"),
    ("vaf", "antigen", "both", "supports clonality/allele fraction", "clonality unknown"),
    ("clonality", "antigen", "both", "prioritizes clonal antigens", "antigen durability unknown"),
    ("immune_context_score", "immune_context", "both", "immune-competent tumor context", "immune gate unknown"),
    ("tls_score", "immune_context", "THCA", "supports HT/TLS immune-active thyroid context", "THCA immune-active branch weak"),
    ("ifng_score", "immune_context", "both", "supports inflamed/IFNG context", "immune gate unknown"),
    ("cytolytic_score", "immune_context", "both", "supports CD8/cytolytic activity", "immune gate unknown"),
    ("autoimmune_history", "safety", "both", "immune-related toxicity and thyroiditis risk", "safety gate unknown"),
    ("thyroiditis_status", "safety", "THCA", "Hashimoto/GD/thyroid autoimmunity safety context", "THCA safety gate unknown"),
    ("steroid_use", "immune_context", "both", "immune competence and vaccine response risk", "immune gate unknown"),
    ("alc", "immune_context", "both", "lymphocyte competence", "immune gate unknown"),
    ("ecog", "safety", "both", "fitness for combination therapy", "safety gate unknown"),
    ("ctdna_status", "disease_context", "PAAD", "MRD monitoring and response window", "MRD branch unavailable"),
    ("ca19_9", "disease_context", "PAAD", "PAAD disease monitoring", "PAAD monitoring unavailable"),
    ("rai_refractory", "disease_context", "THCA", "progressive RR-DTC branch", "RR-DTC branch unavailable"),
    ("manufacturing_feasible", "manufacturing", "both", "tissue/RNA/QC/turnaround feasibility", "manufacturing gate unknown"),
    ("tumor_tissue_qc", "manufacturing", "both", "input material adequacy", "manufacturing gate unknown"),
    ("wes_rna_qc", "manufacturing", "both", "private neoantigen calling feasibility", "manufacturing gate unknown"),
    ("turnaround_window", "manufacturing", "both", "whether vaccine manufacture fits disease tempo", "manufacturing gate unknown"),
]


SCENARIOS = [
    {
        "scenario_id": "PAAD_RESECTED_PERSONALIZED",
        "disease": "PAAD",
        "disease_segment": "PAAD resected/personalized",
        "scenario_name": "PAAD resected/MRD low-burden personalized vaccine window",
        "disease_context_gate": 1.00,
        "presentation_gate": 0.70,
        "antigen_gate": 0.80,
        "immune_context_gate": 0.70,
        "safety_gate": 0.80,
        "research_priority": "high_if_metadata_complete",
        "required_evidence": "resected or MRD-positive low-burden PDAC; sufficient private neoantigens; intact HLA/presentation; fit for standard therapy",
    },
    {
        "scenario_id": "PAAD_MRD_KRAS_PUBLIC",
        "disease": "PAAD",
        "disease_segment": "PAAD MRD+ KRAS-public",
        "scenario_name": "PAAD MRD-positive KRAS public-vaccine branch",
        "disease_context_gate": 1.00,
        "presentation_gate": 0.70,
        "antigen_gate": 0.65,
        "immune_context_gate": 0.70,
        "safety_gate": 0.80,
        "research_priority": "high_if_kras_hla_match",
        "required_evidence": "MRD-positive PDAC; covered KRAS/NRAS mutation; HLA context; ctDNA/CA19-9 monitoring",
    },
    {
        "scenario_id": "PAAD_RESECTED_KRAS_BROAD",
        "disease": "PAAD",
        "disease_segment": "PAAD resected KRAS broad",
        "scenario_name": "PAAD resected KRAS broad long-peptide branch",
        "disease_context_gate": 0.90,
        "presentation_gate": 0.65,
        "antigen_gate": 0.60,
        "immune_context_gate": 0.65,
        "safety_gate": 0.65,
        "research_priority": "conditional",
        "required_evidence": "resected KRAS-mutant PDAC after standard therapy; immune-toxicity fitness; driver mutation coverage",
    },
    {
        "scenario_id": "THCA_BRAF_V600E_ATC",
        "disease": "THCA",
        "disease_segment": "THCA BRAF V600E ATC",
        "scenario_name": "THCA BRAF V600E ATC with targeted backbone first",
        "disease_context_gate": 0.80,
        "presentation_gate": 0.60,
        "antigen_gate": 0.55,
        "immune_context_gate": 0.70,
        "safety_gate": 0.65,
        "research_priority": "conditional_add_on",
        "required_evidence": "ATC/PDTC context; BRAF V600E; HLA/presentation intact; vaccine is add-on/trial logic, not replacement for backbone",
    },
    {
        "scenario_id": "THCA_RR_DTC_PDTC_HIGH_RISK",
        "disease": "THCA",
        "disease_segment": "THCA RR-DTC/PDTC",
        "scenario_name": "THCA progressive RR-DTC/PDTC/high-risk recurrence",
        "disease_context_gate": 0.80,
        "presentation_gate": 0.65,
        "antigen_gate": 0.60,
        "immune_context_gate": 0.65,
        "safety_gate": 0.70,
        "research_priority": "conditional",
        "required_evidence": "progressive radioiodine-refractory DTC, PDTC, or high-risk recurrence; intact presentation; measurable immune context",
    },
    {
        "scenario_id": "THCA_HT_TLS_IMMUNE_ACTIVE_PTC",
        "disease": "THCA",
        "disease_segment": "THCA HT/TLS immune-active PTC",
        "scenario_name": "THCA immune-active recurrence context, not routine low-risk PTC",
        "disease_context_gate": 0.55,
        "presentation_gate": 0.65,
        "antigen_gate": 0.55,
        "immune_context_gate": 0.90,
        "safety_gate": 0.45,
        "research_priority": "narrow_trial_only",
        "required_evidence": "high-risk or recurrent context; HT/TLS/IFNG evidence; autoimmune/thyroiditis safety review",
    },
    {
        "scenario_id": "PRESENTATION_FAILURE_NO_GO",
        "disease": "Both",
        "disease_segment": "Both HLA/presentation failure",
        "scenario_name": "No-go presentation failure branch",
        "disease_context_gate": 0.00,
        "presentation_gate": 0.00,
        "antigen_gate": 0.00,
        "immune_context_gate": 0.00,
        "safety_gate": 0.00,
        "research_priority": "no_go",
        "required_evidence": "HLA LOH, B2M loss, or major presentation failure means vaccine candidate ranking should be deferred",
    },
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", default=".", help="Repository root")
    parser.add_argument("--output-root", required=True, help="Output directory")
    parser.add_argument("--queue-limit", type=int, default=250, help="Top manual-review rows per scenario")
    return parser.parse_args()


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t") if path.exists() else pd.DataFrame()


def present_status(master: pd.DataFrame, field: str) -> dict[str, Any]:
    if master.empty:
        return {"field_present": False, "missing_or_empty": 0, "non_empty": 0, "metadata_status": "master_unavailable"}
    if field not in master.columns:
        return {"field_present": False, "missing_or_empty": len(master), "non_empty": 0, "metadata_status": "absent_from_master"}
    s = master[field]
    empty = s.isna() | s.astype(str).map(lambda x: normalize_empty(x) == "")
    non_empty = int((~empty).sum())
    missing = int(empty.sum())
    if non_empty == 0:
        status = "present_but_empty"
    elif missing > 0:
        status = "present_partial"
    else:
        status = "present_complete"
    return {"field_present": True, "missing_or_empty": missing, "non_empty": non_empty, "metadata_status": status}


def build_metadata_requirements(master: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for field, group, required_for, how_used, default in REQUIREMENTS:
        row = {
            "field": field,
            "gate_group": group,
            "required_for": required_for,
            "how_used": how_used,
            "reviewer_safe_default_if_missing": default,
        }
        row.update(present_status(master, field))
        rows.append(row)
    return pd.DataFrame(rows)


def build_scenarios(combinations: pd.DataFrame) -> pd.DataFrame:
    combo_map: dict[str, dict[str, Any]] = {}
    if not combinations.empty and "disease_segment" in combinations.columns:
        for _, row in combinations.iterrows():
            combo_map[str(row.get("disease_segment", ""))] = row.to_dict()

    rows = []
    for scenario in SCENARIOS:
        row = dict(scenario)
        multiplier = (
            float(row["disease_context_gate"])
            * float(row["presentation_gate"])
            * float(row["antigen_gate"])
            * float(row["immune_context_gate"])
            * float(row["safety_gate"])
        )
        row["scenario_gate_multiplier"] = multiplier
        combo = combo_map.get(row["disease_segment"], {})
        row["combination_strategy"] = combo.get("combination", "")
        row["evidence_grade"] = combo.get("evidence_grade", "")
        row["main_caveat"] = combo.get("main_caveat", "")
        row["patient_selection_from_strategy_matrix"] = combo.get("patient_selection", "")
        rows.append(row)
    return pd.DataFrame(rows)


def choose_queue(output_root: Path) -> pd.DataFrame:
    for name in [
        "barneo_manual_review_queue.tsv",
        "barneo_failure_aware_candidate_scores.tsv",
        "barneo_bma_candidate_scores.tsv",
        "barneo_candidate_scores.tsv",
    ]:
        path = output_root / name
        df = read_tsv(path)
        if not df.empty:
            df = df.copy()
            df["queue_source_file"] = name
            return df
    return pd.DataFrame()


def numeric_series(df: pd.DataFrame, candidates: list[str], default: float = 0.0) -> pd.Series:
    for col in candidates:
        if col in df.columns:
            return pd.to_numeric(df[col], errors="coerce").fillna(default)
    return pd.Series(default, index=df.index, dtype=float)


def build_candidate_queue(queue: pd.DataFrame, scenarios: pd.DataFrame, limit: int) -> pd.DataFrame:
    if queue.empty or scenarios.empty:
        return pd.DataFrame()

    q = queue.copy()
    priority = numeric_series(
        q,
        [
            "manual_review_priority_score",
            "failure_aware_score",
            "patient_gated_bma_score",
            "barneo_bma_score",
            "barneo_score",
        ],
        0.0,
    )
    confidence = numeric_series(
        q,
        ["failure_aware_confidence", "bma_confidence_score", "confidence_score"],
        0.25,
    ).clip(0, 1)
    q["patient_gate_base_score"] = priority.clip(0, 1)
    q["model_confidence_gate"] = confidence.clip(0.05, 1.0)
    if "manual_review_priority_score" in q.columns:
        q = q.sort_values("manual_review_priority_score", ascending=False)
    elif "failure_aware_score" in q.columns:
        q = q.sort_values("failure_aware_score", ascending=False)
    q = q.head(limit).reset_index(drop=True)

    rows = []
    metadata_completion_cap = 0.25
    for _, srow in scenarios.iterrows():
        for _, crow in q.iterrows():
            base = float(crow.get("patient_gate_base_score", 0.0))
            model_gate = float(crow.get("model_confidence_gate", 0.25))
            scenario_mult = float(srow.get("scenario_gate_multiplier", 0.0))
            score = base * scenario_mult * model_gate * metadata_completion_cap
            reasons = [
                "Missing patient/disease metadata prevents real patient-gated confidence",
                "Research triage only; not clinical vaccine selection",
            ]
            inherited = normalize_empty(crow.get("failure_aware_reason_primary", "")) or normalize_empty(crow.get("bma_abstention_reason_primary", ""))
            if inherited:
                reasons.append(f"Benchmark reliability inherited: {inherited}")
            if srow.get("scenario_id") == "PRESENTATION_FAILURE_NO_GO":
                reasons.append("Presentation failure scenario forces vaccine deferral")
            if srow.get("disease") == "THCA":
                reasons.append("THCA routine low-risk PTC remains low priority unless high-risk/progressive context is documented")
            if srow.get("disease") == "PAAD":
                reasons.append("PAAD confident use requires resected/MRD/low-burden timing or explicit KRAS-public branch evidence")
            row = {
                "scenario_id": srow.get("scenario_id"),
                "scenario_name": srow.get("scenario_name"),
                "disease": srow.get("disease"),
                "research_priority": srow.get("research_priority"),
                "candidate_id": crow.get("candidate_id", ""),
                "peptide": crow.get("peptide", ""),
                "hla_allele_4digit": crow.get("hla_allele_4digit", ""),
                "source_name": crow.get("source_name", ""),
                "label": crow.get("label", np.nan),
                "base_score_source": crow.get("queue_source_file", ""),
                "patient_gate_base_score": base,
                "model_confidence_gate": model_gate,
                "scenario_gate_multiplier": scenario_mult,
                "metadata_completion_cap": metadata_completion_cap,
                "demo_patient_gated_score": score,
                "patient_gated_rank_within_scenario": np.nan,
                "patient_metadata_status": "missing_required_patient_context",
                "patient_gated_confidence_bin": "low",
                "hard_gate_status": "metadata_missing",
                "abstain": True,
                "abstention_reason_primary": "Missing patient/disease metadata prevents real patient-gated confidence",
                "abstention_reason_all": "; ".join(dict.fromkeys(reasons)),
                "research_triage_only": True,
                "clinical_use": False,
                "required_evidence": srow.get("required_evidence", ""),
                "combination_strategy": srow.get("combination_strategy", ""),
                "main_caveat": srow.get("main_caveat", ""),
                "manual_review_priority_bin": crow.get("manual_review_priority_bin", ""),
                "failure_aware_review_tier": crow.get("failure_aware_review_tier", ""),
                "review_action_required": crow.get("review_action_required", ""),
                "leakage_risk_level": crow.get("leakage_risk_level", ""),
                "split_low_prevalence": crow.get("split_low_prevalence", ""),
                "selected_methods": crow.get("selected_methods", ""),
            }
            rows.append(row)

    out = pd.DataFrame(rows)
    if not out.empty:
        out["patient_gated_rank_within_scenario"] = (
            out.groupby("scenario_id")["demo_patient_gated_score"].rank(method="first", ascending=False).astype(int)
        )
        out = out.sort_values(["scenario_id", "patient_gated_rank_within_scenario"]).reset_index(drop=True)
    return out


def write_report(
    output_root: Path,
    requirements: pd.DataFrame,
    scenarios: pd.DataFrame,
    queue: pd.DataFrame,
    triage_matrix: pd.DataFrame,
) -> None:
    missing_critical = requirements[requirements["metadata_status"].isin(["absent_from_master", "present_but_empty"])]
    top_queue = queue.sort_values(["scenario_id", "patient_gated_rank_within_scenario"]).head(30) if not queue.empty else queue
    report = f"""
# Patient-Gated CLEAN-Neo PAAD/THCA Demo

## Status

This module is implemented as a reviewer-safe patient-gated research triage layer. The current CLEAN-NeoBench candidate rows do not contain enough PAAD/THCA patient context for real patient-level confidence. Therefore all current candidate-level patient-gated rows are marked `research_triage_only=true`, `clinical_use=false`, `abstain=true`, and `patient_metadata_status=missing_required_patient_context`.

## Practical Algorithm Position

The practical stack is:

1. Agentic ensemble calls benchmarked internal/public/caveated methods.
2. BAR-Neo-BMA uses benchmark-derived Bayesian-style method weights.
3. The QUBO-style selector is only a bounded expert-selection/fusion heuristic; it is not a quantum-advantage claim.
4. Patient gates multiply the BAR-Neo/BMA score only after disease timing, presentation, antigen, immune-context, safety, and model-confidence evidence are available.

## Scenario Gates

{dataframe_to_markdown(scenarios[["scenario_id", "disease", "research_priority", "disease_context_gate", "presentation_gate", "antigen_gate", "immune_context_gate", "safety_gate", "scenario_gate_multiplier", "combination_strategy", "main_caveat"]], max_rows=20)}

## Missing Metadata Blocking Real Patient Gates

{dataframe_to_markdown(missing_critical[["field", "gate_group", "required_for", "metadata_status", "reviewer_safe_default_if_missing"]], max_rows=80)}

## Current Demo Queue

These rows demonstrate how the gate would cap BAR-Neo/BMA candidates under PAAD/THCA contexts. They are not patient decisions.

{dataframe_to_markdown(top_queue[["scenario_id", "candidate_id", "peptide", "hla_allele_4digit", "patient_gate_base_score", "model_confidence_gate", "scenario_gate_multiplier", "metadata_completion_cap", "demo_patient_gated_score", "patient_gated_confidence_bin", "abstention_reason_primary"]], max_rows=30) if not top_queue.empty else "No queue rows available."}

## Strategy Matrix Source

{dataframe_to_markdown(triage_matrix.head(20), max_rows=20) if not triage_matrix.empty else "No patient triage matrix was found."}

## Claim Boundary

Allowed: research triage framework, calibrated prioritization with abstention, PAAD/THCA patient-gate requirements, benchmark-adaptive reliability ranking.

Forbidden: clinical vaccine selection, new SOTA predictor, external validation proven, quantum advantage, public pretrained tools as clean baselines without row-level overlap audit, and pooled Class I/Class II predictor claims.
"""
    (output_root / "PATIENT_GATED_CLEAN_NEO_PAAD_THCA_DEMO.md").write_text(report.strip() + "\n")


def main() -> None:
    args = parse_args()
    repo_root = Path(args.repo_root)
    output_root = Path(args.output_root)
    ensure_dir(output_root)

    master = read_tsv(output_root / "clean_neobench_master.tsv")
    combinations = read_tsv(repo_root / "project/results/p_cancer_vaccine_strategy_2026_05_09/combination_priority_matrix.tsv")
    triage_matrix = read_tsv(repo_root / "project/results/p_cancer_vaccine_strategy_2026_05_09/patient_triage_matrix.tsv")

    requirements = build_metadata_requirements(master)
    scenarios = build_scenarios(combinations)
    queue_source = choose_queue(output_root)
    queue = build_candidate_queue(queue_source, scenarios, args.queue_limit)

    write_tsv(requirements, output_root / "patient_gated_clean_neo_metadata_requirements.tsv")
    write_tsv(scenarios, output_root / "patient_gated_clean_neo_demo_scenarios.tsv")
    write_tsv(queue, output_root / "patient_gated_clean_neo_candidate_queue.tsv")
    write_report(output_root, requirements, scenarios, queue, triage_matrix)

    critical = requirements[requirements["gate_group"].isin(["disease_context", "presentation", "antigen", "immune_context", "safety"])]
    blocked = critical[critical["metadata_status"].isin(["absent_from_master", "present_but_empty"])]
    update_manifest(
        output_root,
        "patient_gated_clean_neo_demo",
        {
            "outputs": [
                "patient_gated_clean_neo_metadata_requirements.tsv",
                "patient_gated_clean_neo_demo_scenarios.tsv",
                "patient_gated_clean_neo_candidate_queue.tsv",
                "PATIENT_GATED_CLEAN_NEO_PAAD_THCA_DEMO.md",
            ],
            "n_patient_gate_scenarios": int(len(scenarios)),
            "n_patient_gated_candidate_rows": int(len(queue)),
            "n_patient_gate_required_fields": int(len(requirements)),
            "n_patient_gate_blocked_critical_fields": int(len(blocked)),
            "warnings": [
                "Patient-gated PAAD/THCA layer is a demo because benchmark candidates lack required patient/disease metadata.",
                "All current patient-gated candidate rows are research triage only and clinical_use=false.",
            ],
        },
    )
    print(
        "[patient-gated-clean-neo] wrote demo rows="
        f"{len(queue)} scenarios={len(scenarios)} blocked_fields={len(blocked)} to {output_root}"
    )


if __name__ == "__main__":
    main()
