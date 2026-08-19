#!/usr/bin/env python3
"""Build P0 CROSS-Neo experiments from the Kaggle-winner playbook merge.

This is intentionally a deterministic integration layer, not a new model sprint.
It applies validation/weighting/gating tricks that transfer cleanly from Kaggle
winner patterns into the existing BAR-Neo, TCR, MD, and patient-gate outputs.
"""

from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
PLAYBOOK_DIR = ROOT / "project/results/cross_neo_kaggle_winner_playbook_2026_05_10"
OUT_DIR = PLAYBOOK_DIR / "merged_p0_experiments"

BAR_DIR = ROOT / "project/results/clean_neobench_barneo_2026_05_09"
MD_DIR = ROOT / "project/results/cross_neo_md_audit_2026_05_10"
HIGH_IMPACT_DIR = MD_DIR / "high_impact_decision_package"
TCR_DIR = ROOT / "project/results/cross_neo_v2_sota_sprint_2026_05_09/tcr_extension"


PATHS = {
    "barneo_candidates": BAR_DIR / "barneo_stress_guarded_candidate_scores.tsv",
    "barneo_weights": BAR_DIR / "barneo_stress_guarded_method_weights.tsv",
    "split_metrics": BAR_DIR / "clean_neobench_split_metrics.tsv",
    "overlap_flags": BAR_DIR / "clean_neobench_overlap_flags.tsv",
    "patient_queue": BAR_DIR / "patient_gated_clean_neo_candidate_queue.tsv",
    "highimpact_top13": HIGH_IMPACT_DIR / "no_false_positive_top13_candidates.tsv",
    "highimpact_optimized": HIGH_IMPACT_DIR / "optimized_preset_candidate_lists.tsv",
    "wetlab_plate_old": HIGH_IMPACT_DIR / "wetlab_validation_plate_plan.tsv",
    "md_evidence": MD_DIR / "md_evidence_scores.tsv",
    "tcr_unique": TCR_DIR / "tcr_wetlab_candidate_prioritization_unique_pmhc.tsv",
}


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
    return (
        filled.astype(str)
        .str.strip()
        .str.lower()
        .isin({"1", "true", "t", "yes", "y"})
    )


def clip01(value: pd.Series | np.ndarray | float) -> pd.Series | np.ndarray | float:
    return np.clip(value, 0.0, 1.0)


def minmax(series: pd.Series, default: float = 0.0) -> pd.Series:
    s = pd.to_numeric(series, errors="coerce")
    if s.notna().sum() == 0:
        return pd.Series(default, index=series.index, dtype="float64")
    lo = float(s.min())
    hi = float(s.max())
    if abs(hi - lo) < 1e-12:
        return pd.Series(default, index=series.index, dtype="float64")
    return ((s - lo) / (hi - lo)).fillna(default)


def first_numeric(df: pd.DataFrame, cols: Iterable[str], default: float = 0.0) -> pd.Series:
    out = pd.Series(np.nan, index=df.index, dtype="float64")
    for col in cols:
        if col in df.columns:
            out = out.fillna(pd.to_numeric(df[col], errors="coerce"))
    return out.fillna(default)


def first_text(df: pd.DataFrame, cols: Iterable[str], default: str = "") -> pd.Series:
    out = pd.Series(pd.NA, index=df.index, dtype="object")
    for col in cols:
        if col in df.columns:
            values = df[col].where(df[col].notna(), pd.NA)
            out = out.fillna(values)
    return out.fillna(default).astype(str)


def safe_cols(df: pd.DataFrame, cols: list[str]) -> list[str]:
    return [col for col in cols if col in df.columns]


def make_method_weights_v2(weights: pd.DataFrame) -> pd.DataFrame:
    if weights.empty:
        return weights

    w = weights.copy()
    family_count = w.groupby("method_family", dropna=False)["method_name"].transform("count")
    role_count = w.groupby("method_role", dropna=False)["method_name"].transform("count")

    stress_floor_parts = pd.concat(
        [
            numeric(w, "source_heldout_min_AUPRC", np.nan),
            numeric(w, "hla_stress_min_AUPRC", np.nan),
            numeric(w, "low_prevalence_AUPRC", np.nan),
            numeric(w, "overall_mean_AUPRC", np.nan),
        ],
        axis=1,
    )
    w["stress_min_floor"] = stress_floor_parts.min(axis=1).fillna(0.0)

    delta_parts = pd.concat(
        [
            numeric(w, "source_heldout_median_delta_AUPRC", np.nan),
            numeric(w, "hla_stress_median_delta_AUPRC", np.nan),
            numeric(w, "korean_hla_median_delta_AUPRC", np.nan),
            numeric(w, "low_prevalence_median_delta_AUPRC", np.nan),
        ],
        axis=1,
    )
    w["stress_delta_floor"] = delta_parts.min(axis=1).fillna(0.0)

    w["redundancy_penalty_proxy"] = 1.0 / (1.0 + 0.08 * (family_count - 1) + 0.04 * (role_count - 1))

    public = truthy(w, "uses_public_pretraining")
    clean_allowed = truthy(w, "clean_comparator_allowed", True)
    public_penalty = np.select(
        [public & ~clean_allowed, public, ~clean_allowed],
        [0.35, 0.70, 0.60],
        default=1.00,
    )
    w["public_or_clean_comparator_penalty"] = public_penalty

    role = text(w, "method_role").str.lower()
    family = text(w, "method_family").str.lower()
    recommended = text(w, "recommended_use").str.lower()
    w["fallback_role_cap"] = np.select(
        [
            role.str.contains("caveated|public", regex=True) | family.str.contains("public", regex=False),
            role.str.contains("fallback", regex=False) | recommended.str.contains("fallback", regex=False),
            role.str.contains("anchor", regex=False),
            role.str.contains("internal", regex=False),
        ],
        [0.25, 0.45, 0.85, 1.00],
        default=0.75,
    )

    top10_norm = minmax(numeric(w, "mean_top10_precision", 0.0))
    calibration_penalty = 1.0 - numeric(w, "mean_ECE", 0.0).clip(0.0, 0.60) * 0.55
    base = first_numeric(w, ["stress_guarded_base_weight", "stress_safe_score", "reviewer_safe_score"], 0.0)
    w["bma_v2_weight_raw"] = (
        base
        * (0.52 + 0.48 * w["stress_min_floor"].clip(0.0, 1.0))
        * (0.65 + 0.35 * top10_norm)
        * w["redundancy_penalty_proxy"]
        * w["public_or_clean_comparator_penalty"]
        * w["fallback_role_cap"]
        * calibration_penalty
    )
    denom = w["bma_v2_weight_raw"].sum()
    w["bma_v2_weight"] = w["bma_v2_weight_raw"] / denom if denom > 0 else 0.0
    max_raw = w["bma_v2_weight_raw"].max()
    w["bma_v2_utility_norm"] = w["bma_v2_weight_raw"] / max_raw if max_raw > 0 else 0.0

    w["bma_v2_use_class"] = np.select(
        [
            (w["bma_v2_utility_norm"] >= 0.60) & clean_allowed & ~public,
            (w["bma_v2_utility_norm"] >= 0.35) & clean_allowed,
            public | ~clean_allowed,
        ],
        [
            "primary_internal_expert",
            "secondary_diversity_expert",
            "caveated_public_or_overlap_expert",
        ],
        default="bounded_fallback_only",
    )

    order_cols = [
        "method_name",
        "method_role",
        "method_family",
        "bma_v2_weight",
        "bma_v2_utility_norm",
        "bma_v2_use_class",
        "stress_min_floor",
        "stress_delta_floor",
        "redundancy_penalty_proxy",
        "public_or_clean_comparator_penalty",
        "fallback_role_cap",
    ]
    rest = [c for c in w.columns if c not in order_cols]
    return w[order_cols + rest].sort_values(["bma_v2_weight", "stress_min_floor"], ascending=False)


def prepare_highimpact() -> pd.DataFrame:
    frames = []
    for key in ("highimpact_top13", "highimpact_optimized"):
        df = read_tsv(PATHS[key])
        if not df.empty:
            df = df.copy()
            df["source_table_for_merge"] = key
            frames.append(df)
    if not frames:
        return pd.DataFrame()

    hi = pd.concat(frames, ignore_index=True)
    if "row_id" not in hi.columns:
        return pd.DataFrame()
    hi["priority_order_numeric"] = numeric(hi, "priority_order", 9999)
    hi["preset_rank_bonus"] = np.where(text(hi, "source_table_for_merge").eq("highimpact_top13"), 0, 1)
    hi = hi.sort_values(["preset_rank_bonus", "priority_order_numeric", "ultra_priority_score"], ascending=[True, True, False])
    hi = hi.drop_duplicates("row_id", keep="first")
    hi = hi.rename(columns={"row_id": "candidate_id"})

    selected = [
        "candidate_id",
        "pmhc_score_mean",
        "tcr_augmented_score_mean",
        "score_delta_mean",
        "best_tcr_augmented_score",
        "wetlab_priority_score",
        "source_dataset",
        "study_id",
        "assay_type",
        "cancer_type",
        "peptide",
        "wildtype_peptide",
        "hla_4digit",
        "hla_supertype",
        "source_protein",
        "tcr_evidence_sources",
        "tcr_evidence_count",
        "paired_tcr_evidence_count",
        "beta_only_tcr_evidence_count",
        "peptide_hla_tcr_label_count",
        "cancer_context_evidence_count",
        "tcr_binding_positive_count",
        "tcr_binding_negative_count",
        "structure_evidence_count",
        "tcr_link_evidence_basis",
        "tcr_link_claim_status",
        "tcr_link_warning",
        "tcr_evidence_quality_adjustment",
        "wetlab_priority_score_evidence_adjusted",
        "wetlab_priority_tier",
        "n_unique_models",
        "n_model_families",
        "ensemble_mean",
        "ensemble_median",
        "ensemble_std",
        "ensemble_q05",
        "ensemble_q95",
        "bayes_mean",
        "bayes_width_90",
        "perturb_mean",
        "perturb_width_90",
        "dropout_sensitivity",
        "recommendation_tier",
        "ultra_priority_score",
        "md_label",
        "md_score",
        "md_structural_score",
        "control_readiness_score",
        "live_completion_fraction",
        "wt_status",
        "anchor_preserved_decoy",
        "wt_status_control",
        "anchor_preserved_decoy_control",
        "peptide_length",
        "class_i_like_hla",
        "supported_length_main_dl",
        "main_dl_score",
        "dl_uncertainty_width",
        "dropout_unstable",
        "has_tcr_evidence",
        "has_paired_tcr_evidence",
        "source_compatible_tcr",
        "paired_source_compatible_tcr",
        "tcr_rescue_signal",
        "tcr_branch_support",
        "tcr_branch_discordant_low",
        "wt_or_decoy_ready",
        "dl_first_decision",
        "next_action",
        "dl_first_reason",
        "dl_first_priority_score",
        "baker_structural_score",
        "preset_name",
        "call_rule",
        "prediction_outcome",
        "priority_order",
        "source_table_for_merge",
    ]
    selected = safe_cols(hi, selected)
    rename = {
        c: f"{c}_highimpact"
        for c in selected
        if c
        not in {
            "candidate_id",
            "pmhc_score_mean",
            "tcr_augmented_score_mean",
            "score_delta_mean",
            "best_tcr_augmented_score",
            "wetlab_priority_score",
            "wetlab_priority_score_evidence_adjusted",
            "wetlab_priority_tier",
            "tcr_evidence_count",
            "paired_tcr_evidence_count",
            "main_dl_score",
            "md_score",
            "md_label",
            "control_readiness_score",
            "baker_structural_score",
            "has_tcr_evidence",
            "has_paired_tcr_evidence",
            "tcr_branch_discordant_low",
            "wt_or_decoy_ready",
            "prediction_outcome",
            "preset_name",
            "priority_order",
        }
    }
    return hi[selected].rename(columns=rename)


def prepare_patient_aggregate(patient: pd.DataFrame) -> pd.DataFrame:
    if patient.empty or "candidate_id" not in patient.columns:
        return pd.DataFrame()

    p = patient.copy()
    p["patient_missing_required_context"] = text(p, "patient_metadata_status").str.contains("missing", case=False, na=False)
    p["patient_any_abstain"] = truthy(p, "abstain")

    agg = (
        p.groupby("candidate_id", dropna=False)
        .agg(
            patient_scenario_count=("scenario_id", "nunique"),
            patient_best_demo_score=("demo_patient_gated_score", "max"),
            patient_median_demo_score=("demo_patient_gated_score", "median"),
            patient_min_metadata_cap=("metadata_completion_cap", "min"),
            patient_max_model_confidence_gate=("model_confidence_gate", "max"),
            patient_any_missing_required_context=("patient_missing_required_context", "max"),
            patient_any_abstain=("patient_any_abstain", "max"),
            patient_scenarios=("scenario_id", lambda x: "|".join(sorted(set(map(str, x.dropna()))))[:500]),
            patient_required_evidence=("required_evidence", lambda x: "|".join(sorted(set(map(str, x.dropna()))))[:500]),
        )
        .reset_index()
    )
    return agg


def build_candidate_scores_v2(
    candidates: pd.DataFrame,
    overlap: pd.DataFrame,
    highimpact: pd.DataFrame,
    patient_agg: pd.DataFrame,
) -> pd.DataFrame:
    if candidates.empty:
        return pd.DataFrame()

    c = candidates.copy()

    if not overlap.empty and "candidate_id" in overlap.columns:
        overlap_cols = [
            "candidate_id",
            "exact_peptide_train_overlap",
            "exact_peptide_hla_train_overlap",
            "near_peptide_train_overlap",
            "source_protein_window_train_overlap",
            "study_train_overlap",
            "patient_train_overlap",
            "public_tool_training_overlap_any",
            "public_tool_training_overlap_detail",
            "leakage_risk_level",
        ]
        overlap_small = overlap[safe_cols(overlap, overlap_cols)].drop_duplicates("candidate_id")
        c = c.merge(overlap_small, on="candidate_id", how="left", suffixes=("", "_overlap_table"))
        if "leakage_risk_level_overlap_table" in c.columns:
            c["leakage_risk_level"] = c["leakage_risk_level"].fillna(c["leakage_risk_level_overlap_table"])

    if not highimpact.empty:
        c = c.merge(highimpact, on="candidate_id", how="left")

    if not patient_agg.empty:
        c = c.merge(patient_agg, on="candidate_id", how="left")

    clean_score = first_numeric(c, ["stress_guarded_clean_score", "clean_contextual_bma_score", "contextual_bma_score"], 0.0)
    claim_score = first_numeric(c, ["stress_guarded_final_review_score", "stress_guarded_claim_safe_score", "barneo_x_claim_safe_score"], 0.0)
    discovery_score = first_numeric(c, ["stress_guarded_discovery_score", "barneo_x_discovery_score"], 0.0)
    main_dl = first_numeric(c, ["main_dl_score"], np.nan)
    tcr_aug = first_numeric(c, ["tcr_augmented_score_mean", "best_tcr_augmented_score"], np.nan)
    md_score = first_numeric(c, ["md_score"], np.nan)
    baker_score = first_numeric(c, ["baker_structural_score"], np.nan)

    main_dl_fill = main_dl.fillna(clean_score)
    tcr_aug_fill = tcr_aug.fillna(0.50)
    md_fill = md_score.fillna(0.65)
    baker_fill = baker_score.fillna(0.50)

    c["antigen_model_gate"] = clip01(0.50 * clean_score + 0.25 * main_dl_fill + 0.15 * first_numeric(c, ["ensemble_median", "ensemble_mean"], clean_score) + 0.10 * baker_fill)

    hla_support = numeric(c, "hla_allele_support_count", np.nan)
    presentation_gate = 0.35 + 0.65 * clean_score
    presentation_gate = presentation_gate - np.where(truthy(c, "is_underrepresented_hla"), 0.12, 0.0)
    presentation_gate = presentation_gate - np.where(truthy(c, "is_mhc_ii"), 0.10, 0.0)
    presentation_gate = presentation_gate - np.where(hla_support.fillna(9999) < 25, 0.10, 0.0)
    c["presentation_hla_gate"] = clip01(presentation_gate)

    paired_tcr = truthy(c, "has_paired_tcr_evidence") | (numeric(c, "paired_tcr_evidence_count", 0.0) > 0)
    any_tcr = truthy(c, "has_tcr_evidence") | (numeric(c, "tcr_evidence_count", 0.0) > 0)
    tcr_count_bonus = np.log1p(numeric(c, "tcr_evidence_count", 0.0).clip(0, 200)) / np.log1p(200)
    tcr_gate = np.select(
        [paired_tcr, any_tcr],
        [
            0.48 + 0.45 * tcr_aug_fill + 0.07 * tcr_count_bonus,
            0.42 + 0.35 * tcr_aug_fill + 0.05 * tcr_count_bonus,
        ],
        default=0.62,
    )
    tcr_gate = np.minimum(tcr_gate, np.where(truthy(c, "tcr_branch_discordant_low"), 0.45, 1.0))
    c["tcr_or_unknown_gate"] = clip01(tcr_gate)

    md_label = first_text(c, ["md_label"], "").str.upper()
    md_gate = np.select(
        [
            md_label.str.contains("STRONG", regex=False),
            md_label.str.contains("MODERATE", regex=False),
            md_label.str.contains("LOW|WEAK|FAIL", regex=True),
        ],
        [0.92, 0.78, 0.40],
        default=md_fill,
    )
    c["md_structure_gate"] = clip01(md_gate)

    patient_cap = numeric(c, "patient_min_metadata_cap", np.nan)
    patient_missing = truthy(c, "patient_any_missing_required_context") | truthy(c, "patient_any_abstain")
    patient_gate = patient_cap.fillna(0.42)
    patient_gate = np.where(patient_missing, np.minimum(patient_gate, 0.30), patient_gate)
    c["patient_context_gate"] = clip01(patient_gate)

    exact_overlap = truthy(c, "exact_peptide_train_overlap") | truthy(c, "exact_peptide_hla_train_overlap")
    near_overlap = truthy(c, "near_peptide_train_overlap") | truthy(c, "source_protein_window_train_overlap")
    public_overlap = truthy(c, "public_tool_training_overlap_any")
    leakage = text(c, "leakage_risk_level").str.lower()
    c["overlap_clean_cap"] = np.select(
        [
            exact_overlap | leakage.eq("high"),
            near_overlap | leakage.eq("medium"),
            public_overlap,
        ],
        [0.42, 0.55, 0.60],
        default=1.00,
    )

    source_hla_cap = np.ones(len(c), dtype="float64")
    source_hla_cap = np.minimum(source_hla_cap, np.where(truthy(c, "is_low_prevalence_context"), 0.72, 1.00))
    source_hla_cap = np.minimum(source_hla_cap, np.where(truthy(c, "is_underrepresented_hla"), 0.70, 1.00))
    source_hla_cap = np.minimum(source_hla_cap, np.where(leakage.eq("medium"), 0.80, 1.00))
    c["source_hla_generalization_cap"] = clip01(source_hla_cap)

    gate_cols = [
        "antigen_model_gate",
        "presentation_hla_gate",
        "tcr_or_unknown_gate",
        "md_structure_gate",
        "patient_context_gate",
        "overlap_clean_cap",
        "source_hla_generalization_cap",
    ]
    gate_matrix = c[gate_cols].astype(float)
    c["validity_dag_cap"] = gate_matrix.min(axis=1)
    c["primary_claim_blocker"] = gate_matrix.idxmin(axis=1).str.replace("_", " ", regex=False)

    raw_claim_blend = clip01(0.42 * claim_score + 0.20 * clean_score + 0.16 * main_dl_fill + 0.12 * md_fill + 0.10 * tcr_aug_fill)
    c["bma_v2_claim_safe_score"] = np.minimum(raw_claim_blend, c["validity_dag_cap"])
    c["bma_v2_discovery_score"] = clip01(0.38 * discovery_score + 0.22 * main_dl_fill + 0.16 * tcr_aug_fill + 0.12 * md_fill + 0.07 * baker_fill + 0.05 * clean_score)
    c["discovery_minus_claim_gap"] = c["bma_v2_discovery_score"] - c["bma_v2_claim_safe_score"]

    confidence = first_numeric(c, ["stress_guarded_confidence", "contextual_confidence_score"], 0.0)
    c["bma_v2_action"] = np.select(
        [
            c["overlap_clean_cap"] <= 0.42,
            c["patient_context_gate"] <= 0.30,
            c["tcr_or_unknown_gate"] < 0.50,
            (c["bma_v2_claim_safe_score"] >= 0.55) & (c["validity_dag_cap"] >= 0.55) & (confidence >= 0.65),
            c["bma_v2_discovery_score"] >= 0.60,
        ],
        [
            "audit_blocked_overlap_or_leakage",
            "research_triage_only_patient_context_missing",
            "tcr_structure_disagreement_review",
            "claim_safe_candidate_after_manual_audit",
            "wetlab_priority_with_claim_boundary",
        ],
        default="watchlist_or_abstain",
    )
    c["bma_v2_claim_boundary"] = np.select(
        [
            c["bma_v2_action"].eq("claim_safe_candidate_after_manual_audit"),
            c["bma_v2_action"].eq("wetlab_priority_with_claim_boundary"),
            c["bma_v2_action"].eq("research_triage_only_patient_context_missing"),
            c["bma_v2_action"].eq("audit_blocked_overlap_or_leakage"),
        ],
        [
            "benchmark-level claim only after manual overlap and metadata audit",
            "experiment-priority claim; no clinical vaccine selection claim",
            "research triage only until patient metadata is complete",
            "discovery/audit only because overlap or leakage cap is active",
        ],
        default="supporting evidence only",
    )

    c = c.sort_values(["bma_v2_claim_safe_score", "bma_v2_discovery_score"], ascending=False)
    c["bma_v2_rank_claim_safe"] = np.arange(1, len(c) + 1)
    c = c.sort_values(["bma_v2_discovery_score", "bma_v2_claim_safe_score"], ascending=False)
    c["bma_v2_rank_discovery"] = np.arange(1, len(c) + 1)
    c = c.sort_values(["bma_v2_claim_safe_score", "bma_v2_discovery_score"], ascending=False)

    front = [
        "candidate_id",
        "peptide",
        "hla_allele_4digit",
        "label",
        "source_name",
        "source_dataset",
        "leakage_risk_level",
        "bma_v2_claim_safe_score",
        "bma_v2_discovery_score",
        "discovery_minus_claim_gap",
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
        "bma_v2_claim_boundary",
        "bma_v2_rank_claim_safe",
        "bma_v2_rank_discovery",
    ]
    front = safe_cols(c, front)
    return c[front + [col for col in c.columns if col not in front]]


def build_wetlab_plate_v2(highimpact_raw: pd.DataFrame, candidate_v2: pd.DataFrame) -> pd.DataFrame:
    if highimpact_raw.empty:
        base = candidate_v2.head(40).copy()
        if base.empty:
            return base
        base = base.rename(columns={"candidate_id": "row_id", "hla_allele_4digit": "hla_4digit"})
    else:
        base = highimpact_raw.copy()
        if "row_id" not in base.columns and "candidate_id" in base.columns:
            base = base.rename(columns={"candidate_id": "row_id"})
        base["priority_order_numeric"] = numeric(base, "priority_order", 9999)
        base = base.sort_values(["priority_order_numeric", "ultra_priority_score"], ascending=[True, False])
        base = base.drop_duplicates("row_id", keep="first")

    if not candidate_v2.empty and "candidate_id" in candidate_v2.columns:
        v2_cols = [
            "candidate_id",
            "bma_v2_claim_safe_score",
            "bma_v2_discovery_score",
            "bma_v2_action",
            "primary_claim_blocker",
            "validity_dag_cap",
            "patient_context_gate",
            "overlap_clean_cap",
            "bma_v2_claim_boundary",
        ]
        base = base.merge(
            candidate_v2[safe_cols(candidate_v2, v2_cols)].drop_duplicates("candidate_id"),
            left_on="row_id",
            right_on="candidate_id",
            how="left",
        )

    paired = truthy(base, "has_paired_tcr_evidence") | (numeric(base, "paired_tcr_evidence_count", 0.0) > 0)
    any_tcr = truthy(base, "has_tcr_evidence") | (numeric(base, "tcr_evidence_count", 0.0) > 0)
    control_ready = truthy(base, "wt_or_decoy_ready") | (numeric(base, "control_readiness_score", 0.0) >= 0.50)
    md = first_numeric(base, ["md_score", "md_structural_score"], 0.0)
    tcr = first_numeric(base, ["tcr_augmented_score_mean", "best_tcr_augmented_score"], 0.0)
    claim = numeric(base, "bma_v2_claim_safe_score", 0.0)
    discovery = numeric(base, "bma_v2_discovery_score", numeric(base, "ultra_priority_score", 0.0))

    base["plate_v2_order_score"] = clip01(
        0.23 * claim
        + 0.25 * discovery
        + 0.16 * tcr
        + 0.12 * md
        + 0.10 * numeric(base, "main_dl_score", 0.0)
        + 0.07 * numeric(base, "control_readiness_score", 0.0)
        + 0.05 * numeric(base, "baker_structural_score", 0.0)
        + 0.02 * np.where(text(base, "prediction_outcome").eq("TP"), 1.0, 0.0)
    )
    base["plate_v2_tier"] = np.select(
        [
            paired & (md >= 0.55) & control_ready,
            paired & (tcr >= 0.80),
            (numeric(base, "main_dl_score", 0.0) >= 0.65) & (md >= 0.50),
            claim >= 0.50,
            any_tcr,
        ],
        [
            "TIER_1_TCR_MD_CONTROL_READY",
            "TIER_1_TCR_PRIORITY_CONTROL_PENDING",
            "TIER_2_PMHCI_DL_MD_PRIORITY",
            "TIER_2_BENCHMARK_CLAIM_SAFE",
            "TIER_3_TCR_DIAGNOSTIC_ONLY",
        ],
        default="TIER_4_AUDIT_OR_RESERVE",
    )
    base["plate_v2_assay_bundle"] = np.select(
        [
            paired & control_ready,
            paired,
            control_ready,
        ],
        [
            "mutant_vs_wt_decoy_pMHC_binding+TCR_tetramer+T_cell_activation",
            "mutant_pMHC_binding+TCR_tetramer; add wt/decoy before claim",
            "mutant_vs_wt_decoy_pMHC_binding+activation_screen",
        ],
        default="low-cost pMHC binding screen before TCR/MD expansion",
    )
    base["plate_v2_decision_use"] = np.where(
        base["plate_v2_tier"].str.startswith("TIER_1"),
        "P0 wetlab first plate",
        np.where(base["plate_v2_tier"].str.startswith("TIER_2"), "P1 backup or replicate plate", "audit reserve"),
    )
    base["plate_v2_claim_boundary"] = first_text(
        base,
        ["bma_v2_claim_boundary"],
        "experiment-priority only; not a clinical vaccine selection call",
    )

    front = [
        "row_id",
        "peptide",
        "hla_4digit",
        "label_binary",
        "source_dataset",
        "plate_v2_order_score",
        "plate_v2_tier",
        "plate_v2_decision_use",
        "plate_v2_assay_bundle",
        "bma_v2_claim_safe_score",
        "bma_v2_discovery_score",
        "primary_claim_blocker",
        "validity_dag_cap",
        "tcr_augmented_score_mean",
        "paired_tcr_evidence_count",
        "tcr_evidence_count",
        "md_label",
        "md_score",
        "control_readiness_score",
        "wt_status",
        "anchor_preserved_decoy",
        "wt_status_control",
        "anchor_preserved_decoy_control",
        "main_dl_score",
        "baker_structural_score",
        "prediction_outcome",
        "preset_name",
        "next_action",
        "plate_v2_claim_boundary",
    ]
    front = safe_cols(base, front)
    base = base.sort_values(["plate_v2_order_score", "bma_v2_discovery_score"], ascending=False)
    return base[front + [col for col in base.columns if col not in front]]


def classify_split_axis(contract: str, group: str) -> str:
    raw = f"{contract} {group}".lower()
    if "overall" in raw:
        return "overall"
    if "source" in raw:
        return "source_heldout"
    if "study" in raw:
        return "study_heldout"
    if "hla" in raw:
        return "hla_stress"
    if "low_prevalence" in raw or "low prevalence" in raw:
        return "low_prevalence"
    if "korean" in raw:
        return "korean_hla"
    if "overlap" in raw or "exact" in raw or "near" in raw:
        return "overlap_audit"
    return "other"


def build_generalization_contract_v2(split_metrics: pd.DataFrame, method_v2: pd.DataFrame) -> pd.DataFrame:
    if split_metrics.empty:
        base = method_v2.copy()
        if base.empty:
            return base
        base["external_claim_contract_v2"] = "missing_split_metrics_contract_needed"
        return base

    sm = split_metrics.copy()
    sm["contract_axis"] = [
        classify_split_axis(c, g) for c, g in zip(text(sm, "split_contract"), text(sm, "split_group"))
    ]

    rows = []
    for method, mdf in sm.groupby("method_name", dropna=False):
        row = {"method_name": method}
        row["method_family"] = first_text(mdf, ["method_family"], "").iloc[0]
        row["method_role"] = first_text(mdf, ["method_role"], "").iloc[0]
        row["uses_public_pretraining"] = bool(truthy(mdf, "uses_public_pretraining").max())
        row["clean_comparator_allowed"] = bool(truthy(mdf, "clean_comparator_allowed", True).min())
        for axis, adf in mdf.groupby("contract_axis", dropna=False):
            prefix = axis
            row[f"{prefix}_n_splits"] = int(len(adf))
            row[f"{prefix}_median_AUPRC"] = float(numeric(adf, "AUPRC", np.nan).median())
            row[f"{prefix}_min_AUPRC"] = float(numeric(adf, "AUPRC", np.nan).min())
            row[f"{prefix}_median_top10_precision"] = float(numeric(adf, "top10_precision", np.nan).median())
            row[f"{prefix}_median_ECE"] = float(numeric(adf, "calibration_ece", np.nan).median())
        rows.append(row)
    contract = pd.DataFrame(rows)

    if not method_v2.empty:
        merge_cols = [
            "method_name",
            "bma_v2_weight",
            "bma_v2_utility_norm",
            "bma_v2_use_class",
            "stress_min_floor",
            "stress_delta_floor",
            "public_or_clean_comparator_penalty",
            "fallback_role_cap",
        ]
        contract = contract.merge(method_v2[safe_cols(method_v2, merge_cols)], on="method_name", how="left")

    source_ok = numeric(contract, "source_heldout_min_AUPRC", 0.0) >= 0.35
    hla_ok = numeric(contract, "hla_stress_min_AUPRC", 0.0) >= 0.20
    lowprev_ok = numeric(contract, "low_prevalence_median_AUPRC", 0.0) >= 0.35
    top10_ok = numeric(contract, "overall_median_top10_precision", 0.0) >= 0.40
    clean_ok = truthy(contract, "clean_comparator_allowed", True) & ~truthy(contract, "uses_public_pretraining")
    ece_ok = numeric(contract, "overall_median_ECE", 1.0) <= 0.35

    contract["contract_source_heldout_pass"] = source_ok
    contract["contract_hla_stress_pass"] = hla_ok
    contract["contract_low_prevalence_pass"] = lowprev_ok
    contract["contract_topk_pass"] = top10_ok
    contract["contract_calibration_pass"] = ece_ok
    contract["contract_clean_internal_pass"] = clean_ok
    contract["external_claim_contract_v2"] = np.select(
        [
            source_ok & hla_ok & lowprev_ok & top10_ok & ece_ok & clean_ok,
            source_ok & hla_ok & clean_ok,
            clean_ok,
        ],
        [
            "external_claim_allowed_after_case_audit",
            "benchmark_claim_with_low_prevalence_caveat",
            "internal_support_only_until_stress_axes_pass",
        ],
        default="diagnostic_or_ablation_only",
    )
    contract["required_next_evidence"] = np.select(
        [
            ~clean_ok,
            ~source_ok,
            ~hla_ok,
            ~lowprev_ok,
            ~ece_ok,
        ],
        [
            "remove public/overlap dependence or demote to ablation",
            "source-heldout rescue or source-stratified cap",
            "HLA stress rescue or allele support cap",
            "low-prevalence stress rescue",
            "calibration repair before claim use",
        ],
        default="manual overlap and candidate-level metadata audit",
    )
    return contract.sort_values(["external_claim_contract_v2", "bma_v2_weight"], ascending=[True, False])


def prepare_tcr_base(tcr_unique: pd.DataFrame, highimpact_raw: pd.DataFrame) -> pd.DataFrame:
    frames = []
    if not tcr_unique.empty:
        frames.append(tcr_unique.copy())
    if not highimpact_raw.empty:
        frames.append(highimpact_raw.copy())
    if not frames:
        return pd.DataFrame()
    base = pd.concat(frames, ignore_index=True, sort=False)
    if "row_id" not in base.columns and "candidate_id" in base.columns:
        base = base.rename(columns={"candidate_id": "row_id"})
    if "row_id" in base.columns:
        score = first_numeric(base, ["wetlab_priority_score_evidence_adjusted", "ultra_priority_score"], 0.0)
        base["_dedupe_score"] = score
        base = base.sort_values("_dedupe_score", ascending=False).drop_duplicates("row_id", keep="first")
        base = base.drop(columns=["_dedupe_score"])
    return base


def build_tcr_structure_disagreement_queue(
    tcr_unique: pd.DataFrame, highimpact_raw: pd.DataFrame, candidate_v2: pd.DataFrame
) -> pd.DataFrame:
    base = prepare_tcr_base(tcr_unique, highimpact_raw)
    if base.empty:
        return base

    if not candidate_v2.empty:
        v2_cols = [
            "candidate_id",
            "bma_v2_claim_safe_score",
            "bma_v2_discovery_score",
            "primary_claim_blocker",
            "validity_dag_cap",
            "tcr_or_unknown_gate",
            "md_structure_gate",
            "overlap_clean_cap",
            "patient_context_gate",
            "bma_v2_action",
        ]
        base = base.merge(
            candidate_v2[safe_cols(candidate_v2, v2_cols)].drop_duplicates("candidate_id"),
            left_on="row_id",
            right_on="candidate_id",
            how="left",
        )

    pmhc = numeric(base, "pmhc_score_mean", 0.0)
    tcr = first_numeric(base, ["tcr_augmented_score_mean", "best_tcr_augmented_score"], 0.0)
    md = first_numeric(base, ["md_score", "md_structural_score", "md_structure_gate"], np.nan)
    main_dl = numeric(base, "main_dl_score", pmhc)
    paired = truthy(base, "has_paired_tcr_evidence") | (numeric(base, "paired_tcr_evidence_count", 0.0) > 0)

    base["pmhc_low_tcr_high"] = (pmhc < 0.55) & (tcr >= 0.80)
    base["pmhc_high_tcr_low"] = (pmhc >= 0.70) & (tcr < 0.55)
    base["tcr_high_md_low_or_missing"] = (tcr >= 0.80) & (md.fillna(0.0) < 0.50)
    base["model_high_md_low"] = (main_dl >= 0.70) & (md.fillna(0.0) < 0.50)
    base["md_high_model_low"] = (md.fillna(0.0) >= 0.70) & (main_dl < 0.55)
    base["tcr_structure_disagreement_score"] = clip01(
        0.36 * (tcr - pmhc).abs()
        + 0.24 * (main_dl - md.fillna(main_dl)).abs()
        + 0.18 * np.where(base["pmhc_low_tcr_high"] | base["pmhc_high_tcr_low"], 1.0, 0.0)
        + 0.12 * np.where(truthy(base, "tcr_branch_discordant_low"), 1.0, 0.0)
        + 0.10 * np.where(paired, 1.0, 0.4)
    )
    base["tcr_structure_action"] = np.select(
        [
            base["pmhc_low_tcr_high"] & (md.fillna(0.0) >= 0.55) & paired,
            base["pmhc_low_tcr_high"],
            base["pmhc_high_tcr_low"],
            base["tcr_high_md_low_or_missing"],
            base["md_high_model_low"],
            truthy(base, "tcr_branch_discordant_low"),
        ],
        [
            "rescue_with_TCR_MD_but_keep_claim_boundary",
            "manual_TCR_specificity_review_plus_MD",
            "pMHC_only_no_recognition_claim",
            "run_or_repeat_MD_before_plate_claim",
            "MD_supported_model_miss_case_for_error_analysis",
            "deprioritize_or_repeat_assay_due_to_TCR_discordance",
        ],
        default="no_major_TCR_structure_disagreement",
    )

    front = [
        "row_id",
        "peptide",
        "hla_4digit",
        "label_binary",
        "source_dataset",
        "tcr_structure_disagreement_score",
        "tcr_structure_action",
        "pmhc_low_tcr_high",
        "pmhc_high_tcr_low",
        "tcr_high_md_low_or_missing",
        "model_high_md_low",
        "md_high_model_low",
        "pmhc_score_mean",
        "tcr_augmented_score_mean",
        "main_dl_score",
        "md_label",
        "md_score",
        "paired_tcr_evidence_count",
        "tcr_evidence_count",
        "tcr_link_claim_status",
        "tcr_link_warning",
        "bma_v2_discovery_score",
        "bma_v2_claim_safe_score",
        "primary_claim_blocker",
        "bma_v2_action",
    ]
    front = safe_cols(base, front)
    return base.sort_values("tcr_structure_disagreement_score", ascending=False)[front + [c for c in base.columns if c not in front]]


def build_patient_context_gate_v2(patient: pd.DataFrame) -> pd.DataFrame:
    if patient.empty:
        return patient
    p = patient.copy()
    p["patient_missing_required_context"] = text(p, "patient_metadata_status").str.contains("missing", case=False, na=False)
    p["patient_context_cap_v2"] = clip01(
        np.minimum.reduce(
            [
                numeric(p, "metadata_completion_cap", 0.25).to_numpy(),
                numeric(p, "model_confidence_gate", 0.25).to_numpy(),
                np.maximum(numeric(p, "scenario_gate_multiplier", 0.25).to_numpy(), 0.25),
            ]
        )
    )
    p["patient_context_cap_v2"] = np.where(
        p["patient_missing_required_context"] | truthy(p, "abstain"),
        np.minimum(p["patient_context_cap_v2"], 0.25),
        p["patient_context_cap_v2"],
    )
    p["patient_gate_v2_action"] = np.select(
        [
            truthy(p, "clinical_use"),
            p["patient_context_cap_v2"] <= 0.25,
            p["patient_context_cap_v2"] <= 0.50,
        ],
        [
            "manual_override_required_clinical_use_not_allowed_by_default",
            "research_triage_only_missing_patient_context",
            "scenario_demo_only_complete_patient_metadata_needed",
        ],
        default="scenario_supported_research_prioritization_only",
    )
    p["patient_gate_v2_required_unlock"] = np.select(
        [
            p["patient_missing_required_context"],
            numeric(p, "model_confidence_gate", 0.0) < 0.50,
            numeric(p, "scenario_gate_multiplier", 0.0) < 0.50,
        ],
        [
            "fill patient/disease/treatment metadata before patient-level ranking",
            "raise model confidence or keep candidate in abstain bucket",
            "match disease/scenario evidence before prioritization",
        ],
        default="manual overlap audit plus wetlab evidence",
    )
    front = [
        "scenario_id",
        "scenario_name",
        "disease",
        "candidate_id",
        "peptide",
        "hla_allele_4digit",
        "label",
        "patient_context_cap_v2",
        "patient_gate_v2_action",
        "patient_gate_v2_required_unlock",
        "demo_patient_gated_score",
        "patient_gate_base_score",
        "model_confidence_gate",
        "scenario_gate_multiplier",
        "metadata_completion_cap",
        "patient_metadata_status",
        "hard_gate_status",
        "abstain",
        "clinical_use",
        "research_triage_only",
        "required_evidence",
        "combination_strategy",
        "main_caveat",
        "review_action_required",
        "leakage_risk_level",
    ]
    front = safe_cols(p, front)
    return p.sort_values(["patient_context_cap_v2", "demo_patient_gated_score"], ascending=False)[front + [c for c in p.columns if c not in front]]


def build_label_noise_audit_queue(candidate_v2: pd.DataFrame) -> pd.DataFrame:
    if candidate_v2.empty:
        return candidate_v2

    c = candidate_v2.copy()
    label = numeric(c, "label", np.nan)
    discovery = numeric(c, "bma_v2_discovery_score", 0.0)
    claim = numeric(c, "bma_v2_claim_safe_score", 0.0)
    disagreement = numeric(c, "stress_guarded_expert_disagreement", 0.0)
    gap = numeric(c, "discovery_minus_claim_gap", 0.0).clip(lower=0.0)
    leakage_high = text(c, "leakage_risk_level").str.lower().eq("high") | (numeric(c, "overlap_clean_cap", 1.0) <= 0.42)

    model_label_conflict = np.select(
        [
            label.eq(0) & (discovery >= 0.62),
            label.eq(1) & (discovery <= 0.35),
            label.eq(1) & leakage_high & (claim >= 0.35),
            disagreement >= 0.25,
            gap >= 0.30,
        ],
        [1.00, 0.90, 0.80, 0.60, 0.55],
        default=0.0,
    )
    c["label_noise_priority_score"] = clip01(
        0.40 * model_label_conflict
        + 0.20 * disagreement.clip(0.0, 1.0)
        + 0.18 * np.where(leakage_high, 1.0, 0.0)
        + 0.14 * gap.clip(0.0, 1.0)
        + 0.08 * np.where(truthy(c, "is_low_prevalence_context") | truthy(c, "is_underrepresented_hla"), 1.0, 0.0)
    )
    c["label_noise_bucket"] = np.select(
        [
            label.eq(0) & (discovery >= 0.62),
            label.eq(1) & (discovery <= 0.35),
            label.eq(1) & leakage_high,
            disagreement >= 0.25,
            gap >= 0.30,
        ],
        [
            "possible_false_negative_or_assay_miss",
            "possible_false_positive_or_context_mismatch",
            "positive_label_overlap_or_leakage_audit",
            "expert_disagreement_audit",
            "high_discovery_low_claim_boundary_case",
        ],
        default="low_priority_label_noise",
    )
    c["label_noise_next_action"] = np.select(
        [
            c["label_noise_bucket"].eq("possible_false_negative_or_assay_miss"),
            c["label_noise_bucket"].eq("positive_label_overlap_or_leakage_audit"),
            c["label_noise_bucket"].eq("expert_disagreement_audit"),
            c["label_noise_bucket"].eq("high_discovery_low_claim_boundary_case"),
        ],
        [
            "manual assay/source review; consider rescue in wetlab if controls exist",
            "remove from clean claim set until overlap provenance is resolved",
            "inspect method-family contributions and source split behavior",
            "keep as discovery-only and identify the active claim blocker",
        ],
        default="no immediate action",
    )
    front = [
        "candidate_id",
        "peptide",
        "hla_allele_4digit",
        "label",
        "source_name",
        "source_dataset",
        "label_noise_priority_score",
        "label_noise_bucket",
        "label_noise_next_action",
        "bma_v2_discovery_score",
        "bma_v2_claim_safe_score",
        "discovery_minus_claim_gap",
        "stress_guarded_expert_disagreement",
        "leakage_risk_level",
        "overlap_clean_cap",
        "primary_claim_blocker",
        "bma_v2_action",
    ]
    front = safe_cols(c, front)
    return c.sort_values(["label_noise_priority_score", "bma_v2_discovery_score"], ascending=False)[front + [col for col in c.columns if col not in front]]


def make_report(
    method_v2: pd.DataFrame,
    candidate_v2: pd.DataFrame,
    wetlab_v2: pd.DataFrame,
    generalization: pd.DataFrame,
    disagreement: pd.DataFrame,
    patient_v2: pd.DataFrame,
    label_noise: pd.DataFrame,
    summary: dict,
) -> str:
    top_claim = candidate_v2.head(8)[
        safe_cols(
            candidate_v2,
            ["candidate_id", "peptide", "hla_allele_4digit", "bma_v2_claim_safe_score", "bma_v2_discovery_score", "bma_v2_action", "primary_claim_blocker"],
        )
    ]
    top_plate = wetlab_v2.head(8)[
        safe_cols(
            wetlab_v2,
            ["row_id", "peptide", "hla_4digit", "plate_v2_order_score", "plate_v2_tier", "plate_v2_assay_bundle"],
        )
    ]
    top_noise = label_noise.head(8)[
        safe_cols(
            label_noise,
            ["candidate_id", "peptide", "hla_allele_4digit", "label_noise_priority_score", "label_noise_bucket", "label_noise_next_action"],
        )
    ]

    def md_table(df: pd.DataFrame) -> str:
        if df.empty:
            return "_no rows_"
        return df.to_markdown(index=False, floatfmt=".3f")

    action_counts = candidate_v2["bma_v2_action"].value_counts().to_dict() if "bma_v2_action" in candidate_v2 else {}
    plate_counts = wetlab_v2["plate_v2_tier"].value_counts().to_dict() if "plate_v2_tier" in wetlab_v2 else {}
    contract_counts = generalization["external_claim_contract_v2"].value_counts().to_dict() if "external_claim_contract_v2" in generalization else {}
    patient_counts = patient_v2["patient_gate_v2_action"].value_counts().to_dict() if "patient_gate_v2_action" in patient_v2 else {}

    return f"""# CROSS-Neo Kaggle-merged P0 experiment report

Generated: {summary["generated_at"]}

## One-line result

Kaggle winner playbook에서 바로 옮길 수 있는 부분은 새 거대 모델이 아니라 **decision controller**였다. BAR-Neo 점수는 discovery와 claim-safe를 분리했고, TCR/MD/wetlab plate는 올리되 patient metadata와 overlap/leakage는 hard cap으로 잠갔다.

## Activated tricks

- OpenVaccine/Ribonanza style: heterogeneous ensemble은 유지하되 validation stress axis에서 무너지는 expert는 cap.
- BELKA/cheminformatics style: high-throughput 후보는 plate-first queue로 만들고, claim은 gate 통과 후보만 허용.
- Open Problems perturbation style: source/HLA/patient distribution shift를 별도 contract로 잠금.
- MoA/CAFA style: noisy label과 method-family disagreement를 별도 audit queue로 뽑음.
- PANDA pathology style: top-line score보다 reproducible operating point와 audit trail을 우선.

## Output counts

- candidate rows: {summary["candidate_rows"]}
- method rows: {summary["method_rows"]}
- wetlab plate rows: {summary["wetlab_plate_rows"]}
- TCR/structure disagreement rows: {summary["tcr_structure_rows"]}
- label-noise audit rows: {summary["label_noise_rows"]}

## Candidate action counts

```json
{json.dumps(action_counts, ensure_ascii=False, indent=2)}
```

## Plate tier counts

```json
{json.dumps(plate_counts, ensure_ascii=False, indent=2)}
```

## Generalization contract counts

```json
{json.dumps(contract_counts, ensure_ascii=False, indent=2)}
```

## Patient gate counts

```json
{json.dumps(patient_counts, ensure_ascii=False, indent=2)}
```

## Top claim-safe/controller candidates

{md_table(top_claim)}

## P0 wetlab plate v2

{md_table(top_plate)}

## Highest label-noise / source-conflict audits

{md_table(top_noise)}

## Decision

1. `claim_safe_decision_controller`는 작동한다. 현재 데이터에서는 patient metadata와 overlap cap이 가장 강한 제한자라서, 높은 discovery 후보도 clinical/patient claim으로 자동 승격하지 않는다.
2. Wetlab은 `plate_v2_tier` 기준으로 진행한다. Tier 1은 TCR/MD/control-readiness가 같이 올라온 후보, Tier 2는 pMHC/DL 또는 benchmark claim 후보, Tier 3-4는 audit reserve다.
3. 다음 paper-blocking 실험은 label-noise/source-conflict 상위 후보와 TCR-structure disagreement 상위 후보의 수동 provenance audit이다. 새 모델 학습은 이 audit 후에만 의미 있다.

## Files

- `barneo_bma_v2_diversity_validity_scores.tsv`
- `barneo_bma_v2_method_weights.tsv`
- `cross_neo_wetlab_plate_v2.tsv`
- `cross_neo_generalization_contract_v2.tsv`
- `tcr_structure_disagreement_queue.tsv`
- `patient_context_gate_v2.tsv`
- `cross_neo_label_noise_audit_queue.tsv`
- `merged_p0_experiment_summary.json`
"""


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    candidates = read_tsv(PATHS["barneo_candidates"], required=True)
    weights = read_tsv(PATHS["barneo_weights"], required=True)
    split_metrics = read_tsv(PATHS["split_metrics"])
    overlap = read_tsv(PATHS["overlap_flags"])
    patient = read_tsv(PATHS["patient_queue"])
    highimpact_raw = read_tsv(PATHS["highimpact_top13"])
    md_evidence = read_tsv(PATHS["md_evidence"])
    tcr_unique = read_tsv(PATHS["tcr_unique"])

    method_v2 = make_method_weights_v2(weights)
    highimpact = prepare_highimpact()
    patient_agg = prepare_patient_aggregate(patient)
    candidate_v2 = build_candidate_scores_v2(candidates, overlap, highimpact, patient_agg)
    wetlab_v2 = build_wetlab_plate_v2(highimpact_raw, candidate_v2)
    generalization_v2 = build_generalization_contract_v2(split_metrics, method_v2)
    disagreement_v2 = build_tcr_structure_disagreement_queue(tcr_unique, highimpact_raw, candidate_v2)
    patient_v2 = build_patient_context_gate_v2(patient)
    label_noise_v2 = build_label_noise_audit_queue(candidate_v2)

    outputs = {
        "barneo_bma_v2_method_weights.tsv": method_v2,
        "barneo_bma_v2_diversity_validity_scores.tsv": candidate_v2,
        "cross_neo_wetlab_plate_v2.tsv": wetlab_v2,
        "cross_neo_generalization_contract_v2.tsv": generalization_v2,
        "tcr_structure_disagreement_queue.tsv": disagreement_v2,
        "patient_context_gate_v2.tsv": patient_v2,
        "cross_neo_label_noise_audit_queue.tsv": label_noise_v2,
    }
    output_paths = {name: str(write_tsv(df, name).relative_to(ROOT)) for name, df in outputs.items()}

    summary = {
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "output_dir": str(OUT_DIR.relative_to(ROOT)),
        "input_paths": {key: str(path.relative_to(ROOT)) for key, path in PATHS.items() if path.exists()},
        "output_paths": output_paths,
        "candidate_rows": int(len(candidate_v2)),
        "method_rows": int(len(method_v2)),
        "wetlab_plate_rows": int(len(wetlab_v2)),
        "generalization_rows": int(len(generalization_v2)),
        "tcr_structure_rows": int(len(disagreement_v2)),
        "patient_gate_rows": int(len(patient_v2)),
        "label_noise_rows": int(len(label_noise_v2)),
        "md_evidence_rows_read": int(len(md_evidence)),
        "top_candidate_actions": candidate_v2["bma_v2_action"].value_counts().to_dict()
        if "bma_v2_action" in candidate_v2
        else {},
        "top_plate_tiers": wetlab_v2["plate_v2_tier"].value_counts().to_dict()
        if "plate_v2_tier" in wetlab_v2
        else {},
    }

    summary_path = OUT_DIR / "merged_p0_experiment_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n")

    report = make_report(
        method_v2=method_v2,
        candidate_v2=candidate_v2,
        wetlab_v2=wetlab_v2,
        generalization=generalization_v2,
        disagreement=disagreement_v2,
        patient_v2=patient_v2,
        label_noise=label_noise_v2,
        summary=summary,
    )
    (OUT_DIR / "MERGED_P0_EXPERIMENT_REPORT_KR.md").write_text(report)

    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
