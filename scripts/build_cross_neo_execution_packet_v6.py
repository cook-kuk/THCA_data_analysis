#!/usr/bin/env python3
"""Build the CROSS-Neo v6 assay execution packet.

v6 converts the preregistered v5 plan into lab-facing manifests, blank result
entry sheets, a statistical analysis packet, an interpreter dry run, and a live
dashboard. It does not change the v5 thresholds or promote clinical claims.
"""

from __future__ import annotations

import hashlib
import json
import shutil
import subprocess
import zipfile
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "project/results/cross_neo_kaggle_winner_playbook_2026_05_10"
V4_DIR = BASE / "translational_impact_v4_2026_05_10"
V5_DIR = BASE / "preregistered_impact_v5_2026_05_10"
OUT_DIR = BASE / "execution_packet_v6_2026_05_10"
FIG_DIR = OUT_DIR / "figures"

HUB_DIR = ROOT / "project/papers_hub_2026_05_04"
HUB_ASSET_DIR = HUB_DIR / "assets/cross_neo_execution_v6"
LIVE_HUB_DIR = Path("/var/www/papers/papers_hub_2026_05_04")
LIVE_ASSET_DIR = LIVE_HUB_DIR / "assets/cross_neo_execution_v6"

PLATE_V4 = V4_DIR / "wetlab_plate_v4_value_of_information.tsv"
ENDPOINT_V5 = V5_DIR / "preregistered_endpoint_plan_v5.tsv"
WELLMAP_V5 = V5_DIR / "assay_96well_map_preregistered_v5.tsv"
RESULT_TEMPLATE_V5 = V5_DIR / "assay_result_template_v5.tsv"
LADDER_V5 = V5_DIR / "claim_unlock_ladder_v5.tsv"
SAFEGUARDS_V5 = V5_DIR / "reviewer_safeguards_v5.tsv"
INTERPRETER = ROOT / "scripts/interpret_cross_neo_assay_results_v6.py"


ARM_BUNDLE = {
    "A_clean_discovery": "mutant pMHC binding + WT/decoy specificity + technical replicate",
    "B_mechanism_TCR_MD": "mutant pMHC binding + TCR readout + WT/decoy context",
    "C_label_rescue": "mutant pMHC binding + provenance resolution",
    "D_specificity_moat": "hard-negative weak-binding check + specificity controls",
    "E_positive_QC_control": "positive-control dynamic range and QC replicate",
    "F_model_boundary": "TCR/structure disagreement resolution",
}

ARM_SUCCESS = {
    "A_clean_discovery": "PASS if mutant pMHC binding is positive and WT/decoy specificity is acceptable",
    "B_mechanism_TCR_MD": "PASS if mutant pMHC binding is positive and TCR readout supports recognition",
    "C_label_rescue": "PASS if mutant pMHC binding is positive and provenance audit resolves the label-noise hypothesis",
    "D_specificity_moat": "PASS if hard-negative candidate remains weak/non-binding under the same assay context",
    "E_positive_QC_control": "PASS if positive control shows expected dynamic range and QC replicate passes",
    "F_model_boundary": "PASS if the TCR/structure disagreement is resolved by orthogonal readout or review",
}

ARM_OWNER = {
    "A_clean_discovery": "binding_specificity_lane",
    "B_mechanism_TCR_MD": "mechanism_tcr_md_lane",
    "C_label_rescue": "provenance_rescue_lane",
    "D_specificity_moat": "specificity_moat_lane",
    "E_positive_QC_control": "assay_qc_lane",
    "F_model_boundary": "model_boundary_lane",
}


def read_tsv(path: Path, required: bool = True) -> pd.DataFrame:
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


def safe_cols(df: pd.DataFrame, cols: list[str]) -> list[str]:
    return [col for col in cols if col in df.columns]


def fmt(value: object, digits: int = 3) -> str:
    if pd.isna(value):
        return ""
    if isinstance(value, (float, np.floating)):
        return f"{float(value):.{digits}f}"
    return str(value)


def stable_id(*parts: object, prefix: str = "XNV") -> str:
    payload = "|".join(str(p) for p in parts)
    digest = hashlib.sha1(payload.encode("utf-8")).hexdigest()[:8].upper()
    return f"{prefix}-{digest}"


def clean_sequence(value: object, fallback: str) -> str:
    text = "" if pd.isna(value) else str(value).strip()
    if text and text.lower() not in {"nan", "none", "na"}:
        return text
    return fallback


def build_candidate_manifest(plate: pd.DataFrame, endpoint: pd.DataFrame) -> pd.DataFrame:
    merged = plate.merge(
        endpoint[
            [
                "plate_v4_arm",
                "endpoint",
                "threshold_successes",
                "confirmatory_threshold_successes",
                "preregistered_unlock_rule",
                "confirmatory_unlock_rule",
                "main_text_claim_tier",
                "claim_layer_after_unlock",
                "false_unlock_risk_under_null",
                "confirmatory_false_unlock_risk_under_null",
                "endpoint_priority_score",
            ]
        ],
        on="plate_v4_arm",
        how="left",
    )
    merged["execution_blind_id"] = merged.apply(lambda r: stable_id(r["plate_v4_arm"], r["candidate_id"], prefix="CAND"), axis=1)
    merged["assay_bundle_required"] = merged["plate_v4_arm"].map(ARM_BUNDLE).fillna("candidate review")
    merged["candidate_success_definition"] = merged["plate_v4_arm"].map(ARM_SUCCESS).fillna("PASS/FAIL manual review")
    merged["execution_lane_owner"] = merged["plate_v4_arm"].map(ARM_OWNER).fillna("review_lane")
    merged["order_priority_score"] = (
        0.35 * pd.to_numeric(merged.get("v4_claim_unlock_score", 0), errors="coerce").fillna(0)
        + 0.25 * pd.to_numeric(merged.get("v4_expected_information_gain", 0), errors="coerce").fillna(0)
        + 0.20 * pd.to_numeric(merged.get("v4_assay_feasibility_score", 0), errors="coerce").fillna(0)
        + 0.20 * pd.to_numeric(merged.get("endpoint_priority_score", 0), errors="coerce").fillna(0)
    )
    merged["order_batch"] = np.select(
        [
            merged["plate_v4_arm"].isin(["B_mechanism_TCR_MD", "A_clean_discovery"]),
            merged["plate_v4_arm"].isin(["C_label_rescue", "D_specificity_moat"]),
            merged["plate_v4_arm"].eq("E_positive_QC_control"),
        ],
        ["batch_1_claim_unlock", "batch_2_risk_resolution", "batch_0_assay_qc"],
        default="batch_3_boundary_review",
    )
    merged["result_entry_row_id"] = merged.apply(lambda r: stable_id(r["plate_v4_slot"], r["candidate_id"], prefix="RES"), axis=1)
    merged["claim_boundary_execution"] = np.where(
        merged["plate_v4_arm"].eq("E_positive_QC_control"),
        "assay QC only; cannot support novelty",
        merged.get("claim_boundary_v3", "computational/diagnostic boundary remains"),
    )
    cols = [
        "plate_v4_slot",
        "plate_v4_arm",
        "execution_blind_id",
        "result_entry_row_id",
        "candidate_id",
        "peptide",
        "hla_allele_4digit",
        "impact_lane",
        "label",
        "source_name",
        "leakage_risk_level",
        "assay_bundle_required",
        "candidate_success_definition",
        "endpoint",
        "preregistered_unlock_rule",
        "confirmatory_unlock_rule",
        "main_text_claim_tier",
        "claim_layer_after_unlock",
        "false_unlock_risk_under_null",
        "confirmatory_false_unlock_risk_under_null",
        "v4_claim_unlock_score",
        "v4_expected_information_gain",
        "v4_lane_specific_prior",
        "order_priority_score",
        "order_batch",
        "execution_lane_owner",
        "wildtype_peptide",
        "wt_status",
        "anchor_preserved_decoy_control",
        "claim_boundary_execution",
    ]
    return merged[safe_cols(merged, cols)].sort_values(["order_batch", "order_priority_score"], ascending=[True, False])


def build_execution_wellmap(wellmap: pd.DataFrame, manifest: pd.DataFrame) -> pd.DataFrame:
    mcols = [
        "candidate_id",
        "execution_blind_id",
        "result_entry_row_id",
        "assay_bundle_required",
        "candidate_success_definition",
        "order_batch",
        "execution_lane_owner",
    ]
    out = wellmap.merge(manifest[safe_cols(manifest, mcols)], on="candidate_id", how="left")
    out["well_blind_id"] = out.apply(lambda r: stable_id(r["well"], r["candidate_id"], r["assay_condition"], prefix="WELL"), axis=1)
    out["entry_field"] = np.select(
        [
            out["assay_condition"].eq("MUT_pMHC"),
            out["assay_condition"].str.contains("WT|DECOY", regex=True, na=False),
            out["assay_condition"].str.contains("TCR", regex=True, na=False),
            out["assay_condition"].str.contains("PROVENANCE", regex=True, na=False),
            out["assay_condition"].str.contains("QC|TECH", regex=True, na=False),
            out["assay_condition"].str.contains("NEG", regex=True, na=False),
        ],
        [
            "mut_pmhc_call",
            "specificity_pass_call",
            "tcr_or_provenance_pass_call",
            "tcr_or_provenance_pass_call",
            "technical_qc_call",
            "hard_negative_weak_call",
        ],
        default="well_call",
    )
    out["allowed_call_values"] = "PASS|FAIL|PENDING|EXCLUDE"
    out["candidate_success_dependency"] = out["candidate_success_definition"].fillna("manual PASS/FAIL")
    cols = [
        "well",
        "plate_row",
        "plate_col",
        "well_blind_id",
        "plate_v4_slot",
        "plate_v4_arm",
        "execution_blind_id",
        "candidate_id",
        "peptide",
        "hla_allele_4digit",
        "impact_lane",
        "assay_condition",
        "readout",
        "expected_interpretation",
        "entry_field",
        "allowed_call_values",
        "candidate_success_dependency",
        "endpoint",
        "preregistered_unlock_rule",
        "order_batch",
        "execution_lane_owner",
        "claim_boundary_v3",
    ]
    return out[safe_cols(out, cols)].sort_values(["plate_col", "plate_row"])


def build_order_manifest(plate: pd.DataFrame, wellmap: pd.DataFrame, manifest: pd.DataFrame) -> pd.DataFrame:
    condition_lookup = wellmap.groupby("candidate_id")["assay_condition"].apply(lambda s: set(s.dropna().astype(str))).to_dict()
    manifest_lookup = manifest.set_index("candidate_id").to_dict("index")
    rows: list[dict[str, object]] = []

    for _, row in plate.iterrows():
        cid = row["candidate_id"]
        conditions = condition_lookup.get(cid, set())
        mrow = manifest_lookup.get(cid, {})
        base = {
            "candidate_id": cid,
            "plate_v4_arm": row.get("plate_v4_arm", ""),
            "peptide": row.get("peptide", ""),
            "hla_allele_4digit": row.get("hla_allele_4digit", ""),
            "impact_lane": row.get("impact_lane", ""),
            "order_batch": mrow.get("order_batch", ""),
            "execution_blind_id": mrow.get("execution_blind_id", ""),
            "claim_boundary_execution": mrow.get("claim_boundary_execution", ""),
        }

        def add(reagent_type: str, sequence_or_item: str, purpose: str, required: str = "yes") -> None:
            line = dict(base)
            line.update(
                {
                    "order_line_id": stable_id(cid, reagent_type, sequence_or_item, prefix="ORD"),
                    "reagent_type": reagent_type,
                    "sequence_or_item": sequence_or_item,
                    "purpose": purpose,
                    "required_for_endpoint": required,
                }
            )
            rows.append(line)

        add("MUTANT_PEPTIDE", clean_sequence(row.get("peptide", ""), "manual_mutant_sequence_review_needed"), "primary mutant pMHC evidence")
        add("HLA_CONTEXT", clean_sequence(row.get("hla_allele_4digit", ""), "manual_hla_review_needed"), "restriction context")

        if "WT_CONTROL" in conditions:
            add(
                "WT_OR_NATIVE_CONTROL_PEPTIDE",
                clean_sequence(row.get("wildtype_peptide", ""), f"manual_wt_control_needed_for_{cid}"),
                "specificity control",
            )
        if "DECOY_CONTROL" in conditions:
            add(
                "ANCHOR_PRESERVED_DECOY_PEPTIDE",
                clean_sequence(row.get("anchor_preserved_decoy_control", ""), f"manual_decoy_control_needed_for_{cid}"),
                "specificity decoy control",
            )
        if any("TCR" in condition for condition in conditions):
            add(
                "TCR_READOUT_CONTEXT",
                f"paired_TCR_or_template_review; evidence_count={fmt(row.get('tcr_evidence_count', ''))}; md_label={fmt(row.get('md_label', ''))}",
                "orthogonal TCR/structure resolution",
            )
        if any("PROVENANCE" in condition for condition in conditions):
            add("PROVENANCE_AUDIT_ITEM", clean_sequence(row.get("source_name", ""), "manual_source_review_needed"), "label-noise rescue provenance")
        if row.get("plate_v4_arm") == "E_positive_QC_control":
            add("POSITIVE_QC_CONTEXT", "assay dynamic-range positive control", "assay QC only")
        if row.get("plate_v4_arm") == "D_specificity_moat":
            add("HARD_NEGATIVE_CONTEXT", "hard negative weak-binding check", "specificity moat")

    out = pd.DataFrame(rows)
    return out.sort_values(["order_batch", "plate_v4_arm", "candidate_id", "reagent_type"])


def build_candidate_result_entry(manifest: pd.DataFrame) -> pd.DataFrame:
    out = manifest.copy()
    for col in [
        "candidate_success_call",
        "mut_pmhc_call",
        "specificity_pass_call",
        "tcr_or_provenance_pass_call",
        "hard_negative_weak_call",
        "technical_qc_call",
        "exclusion_reason",
        "operator_notes",
    ]:
        out[col] = ""
    out["allowed_call_values"] = "PASS|FAIL|PENDING|EXCLUDE"
    out["interpreter_priority"] = "candidate_success_call overrides derived fields if populated"
    cols = [
        "result_entry_row_id",
        "plate_v4_slot",
        "plate_v4_arm",
        "execution_blind_id",
        "candidate_id",
        "peptide",
        "hla_allele_4digit",
        "impact_lane",
        "endpoint",
        "candidate_success_definition",
        "candidate_success_call",
        "mut_pmhc_call",
        "specificity_pass_call",
        "tcr_or_provenance_pass_call",
        "hard_negative_weak_call",
        "technical_qc_call",
        "allowed_call_values",
        "interpreter_priority",
        "exclusion_reason",
        "operator_notes",
        "preregistered_unlock_rule",
        "confirmatory_unlock_rule",
        "claim_boundary_execution",
    ]
    return out[safe_cols(out, cols)].sort_values("plate_v4_slot")


def build_well_result_entry(wellmap: pd.DataFrame) -> pd.DataFrame:
    out = wellmap.copy()
    for col in ["raw_signal", "normalized_signal", "well_call", "qc_flag", "exclude_reason", "operator_notes"]:
        out[col] = ""
    cols = [
        "well",
        "plate_row",
        "plate_col",
        "well_blind_id",
        "plate_v4_arm",
        "execution_blind_id",
        "candidate_id",
        "assay_condition",
        "readout",
        "entry_field",
        "allowed_call_values",
        "raw_signal",
        "normalized_signal",
        "well_call",
        "qc_flag",
        "exclude_reason",
        "operator_notes",
    ]
    return out[safe_cols(out, cols)].sort_values(["plate_col", "plate_row"])


def build_decision_worksheet(endpoint: pd.DataFrame) -> pd.DataFrame:
    out = endpoint.copy()
    for col in [
        "observed_success_count",
        "completed_candidate_count",
        "exploratory_unlock_call",
        "confirmatory_unlock_call",
        "decision_notes",
    ]:
        out[col] = ""
    cols = [
        "plate_v4_arm",
        "endpoint",
        "n_candidates",
        "threshold_successes",
        "confirmatory_threshold_successes",
        "preregistered_unlock_rule",
        "confirmatory_unlock_rule",
        "main_text_claim_tier",
        "false_unlock_risk_under_null",
        "confirmatory_false_unlock_risk_under_null",
        "power_at_mean_prior",
        "confirmatory_power_at_mean_prior",
        "claim_layer_after_unlock",
        "observed_success_count",
        "completed_candidate_count",
        "exploratory_unlock_call",
        "confirmatory_unlock_call",
        "decision_notes",
    ]
    return out[safe_cols(out, cols)].sort_values("plate_v4_arm")


def build_audit_trail(endpoint: pd.DataFrame, safeguards: pd.DataFrame) -> pd.DataFrame:
    rows = [
        {
            "audit_item": "freeze_candidate_list",
            "source_file": "execution_candidate_manifest_v6.tsv",
            "required_before_unlock": "candidate list unchanged or deviations logged",
            "failure_action": "mark deviation and keep claim exploratory",
        },
        {
            "audit_item": "use_predeclared_thresholds",
            "source_file": "preregistered_endpoint_plan_v5.tsv",
            "required_before_unlock": "exploratory and confirmatory thresholds copied without post-hoc tuning",
            "failure_action": "do not promote endpoint",
        },
        {
            "audit_item": "candidate_success_call_trace",
            "source_file": "candidate_result_entry_v6.tsv + candidate_calls_v6.tsv",
            "required_before_unlock": "PASS/FAIL source documented for every counted candidate",
            "failure_action": "treat missing or ambiguous candidates as pending",
        },
        {
            "audit_item": "positive_controls_not_novelty",
            "source_file": "claim_boundary_execution column",
            "required_before_unlock": "E arm used for QC only",
            "failure_action": "remove from novelty/discovery claims",
        },
        {
            "audit_item": "overlap_blocked_candidates_labeled",
            "source_file": "leakage_risk_level + claim_boundary_execution",
            "required_before_unlock": "source-overlap status retained in final table",
            "failure_action": "mechanistic/support only, no clean benchmark novelty",
        },
    ]
    for _, row in safeguards.iterrows():
        rows.append(
            {
                "audit_item": f"reviewer_risk::{row.get('reviewer_risk', '')}",
                "source_file": row.get("evidence_file", ""),
                "required_before_unlock": row.get("v5_safeguard", ""),
                "failure_action": row.get("residual_risk", ""),
            }
        )
    out = pd.DataFrame(rows)
    out["endpoint_plan_checksum"] = stable_id(
        "|".join(endpoint["plate_v4_arm"].astype(str)),
        "|".join(endpoint["confirmatory_unlock_rule"].astype(str)),
        prefix="PLAN",
    )
    return out


def make_figures(
    manifest: pd.DataFrame,
    order_manifest: pd.DataFrame,
    wellmap: pd.DataFrame,
    endpoint: pd.DataFrame,
) -> list[Path]:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    figures: list[Path] = []

    stage_counts = pd.Series(
        {
            "candidates": len(manifest),
            "order lines": len(order_manifest),
            "plate wells": len(wellmap),
            "endpoints": len(endpoint),
        }
    )
    fig, ax = plt.subplots(figsize=(8.8, 4.8))
    colors = ["#2563eb", "#16a34a", "#f59e0b", "#dc2626"]
    ax.bar(stage_counts.index, stage_counts.values, color=colors)
    for i, value in enumerate(stage_counts.values):
        ax.text(i, value + max(stage_counts.values) * 0.02, str(int(value)), ha="center", va="bottom", fontweight="bold")
    ax.set_title("Execution packet scale")
    ax.set_ylabel("Count")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    path = FIG_DIR / "fig1_execution_packet_scale.png"
    fig.savefig(path, dpi=200)
    plt.close(fig)
    figures.append(path)

    reagent_counts = pd.crosstab(order_manifest["plate_v4_arm"], order_manifest["reagent_type"])
    fig, ax = plt.subplots(figsize=(10.6, 5.6))
    bottom = np.zeros(len(reagent_counts))
    palette = ["#2563eb", "#16a34a", "#f59e0b", "#dc2626", "#7c3aed", "#0f766e", "#64748b", "#db2777"]
    for i, col in enumerate(reagent_counts.columns):
        ax.bar(reagent_counts.index.str.replace("_", " ", regex=False), reagent_counts[col], bottom=bottom, label=col, color=palette[i % len(palette)])
        bottom += reagent_counts[col].to_numpy()
    ax.set_ylabel("Order/audit lines")
    ax.set_title("Reagent and audit manifest by arm")
    ax.tick_params(axis="x", rotation=25)
    ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), fontsize=8, frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    path = FIG_DIR / "fig2_reagent_manifest_by_arm.png"
    fig.savefig(path, dpi=200)
    plt.close(fig)
    figures.append(path)

    field_counts = pd.crosstab(wellmap["plate_v4_arm"], wellmap["entry_field"])
    fig, ax = plt.subplots(figsize=(10.4, 5.6))
    bottom = np.zeros(len(field_counts))
    for i, col in enumerate(field_counts.columns):
        ax.bar(field_counts.index.str.replace("_", " ", regex=False), field_counts[col], bottom=bottom, label=col, color=palette[(i + 2) % len(palette)])
        bottom += field_counts[col].to_numpy()
    ax.set_ylabel("Wells")
    ax.set_title("Result-entry fields embedded in 96-well map")
    ax.tick_params(axis="x", rotation=25)
    ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), fontsize=8, frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    path = FIG_DIR / "fig3_result_entry_fields.png"
    fig.savefig(path, dpi=200)
    plt.close(fig)
    figures.append(path)

    e = endpoint.sort_values("plate_v4_arm")
    y = np.arange(len(e))
    fig, ax = plt.subplots(figsize=(9.6, 5.2))
    ax.barh(y - 0.18, e["threshold_successes"], height=0.32, color="#22c55e", label="exploratory")
    ax.barh(y + 0.18, e["confirmatory_threshold_successes"], height=0.32, color="#f97316", label="confirmatory")
    ax.set_yticks(y)
    ax.set_yticklabels(e["plate_v4_arm"].str.replace("_", " ", regex=False), fontsize=8)
    ax.set_xlabel("Required candidate successes")
    ax.set_title("Frozen unlock thresholds")
    ax.legend(frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    path = FIG_DIR / "fig4_frozen_unlock_thresholds.png"
    fig.savefig(path, dpi=200)
    plt.close(fig)
    figures.append(path)
    return figures


def html_table(df: pd.DataFrame, cols: list[str], n: int = 20) -> str:
    view = df[safe_cols(df, cols)].head(n).copy()
    for col in view.columns:
        if pd.api.types.is_numeric_dtype(view[col]):
            view[col] = view[col].map(lambda x: f"{x:.3f}" if pd.notna(x) else "")
    return view.to_html(index=False, classes="data-table", escape=False)


def make_sap(summary: dict, endpoint: pd.DataFrame, audit: pd.DataFrame) -> str:
    return f"""# CROSS-Neo v6 statistical analysis plan

Generated: {summary["generated_at"]}

## Analysis set

The analysis set is the 24-candidate v4/v5 plate. Candidate additions, exclusions, or replacements must be recorded in `reviewer_audit_trail_v6.tsv` before endpoint interpretation.

## Primary endpoint interpretation

Candidate-level PASS/FAIL calls are entered in `candidate_result_entry_v6.tsv`. The interpreter uses `candidate_success_call` when populated; otherwise it derives a candidate call from the predeclared arm-specific fields.

## Frozen endpoint thresholds

{endpoint.to_markdown(index=False, floatfmt=".3f")}

## Multiple-endpoint boundary

Exploratory thresholds support screening conclusions only. Main-text-ready claims use the confirmatory thresholds. The confirmatory family null-risk sum remains {summary["confirmatory_null_risk_sum"]:.3f}; no endpoint in this packet unlocks a clinical vaccine-selection claim.

## Audit checklist

{audit.to_markdown(index=False)}
"""


def make_report(summary: dict, manifest: pd.DataFrame, endpoint: pd.DataFrame, order_manifest: pd.DataFrame) -> str:
    top = manifest.sort_values("order_priority_score", ascending=False).head(10)
    return f"""# CROSS-Neo execution packet v6

Generated: {summary["generated_at"]}

## What changed

v6 turns the v5 preregistration package into a concrete execution handoff: order manifest, 96-well execution map, blank candidate/well result sheets, frozen decision worksheet, and an automatic endpoint interpreter.

## Headline metrics

- Candidates: {summary["candidate_count"]}
- 96-well entries: {summary["well_count"]}
- Reagent/audit order lines: {summary["order_line_count"]}
- Endpoint arms: {summary["endpoint_count"]}
- Confirmatory family null-risk sum: {summary["confirmatory_null_risk_sum"]:.3f}
- Dry-run pending endpoints: {summary["dry_run_pending_endpoints"]}

## Top execution priorities

{top[safe_cols(top, ["plate_v4_arm", "candidate_id", "peptide", "hla_allele_4digit", "order_batch", "order_priority_score", "claim_layer_after_unlock"])].to_markdown(index=False, floatfmt=".3f")}

## Frozen endpoints

{endpoint[safe_cols(endpoint, ["plate_v4_arm", "n_candidates", "threshold_successes", "confirmatory_threshold_successes", "main_text_claim_tier", "confirmatory_false_unlock_risk_under_null"])].to_markdown(index=False, floatfmt=".3f")}

## Interpreter command

```bash
python scripts/interpret_cross_neo_assay_results_v6.py \\
  --results {OUT_DIR / "candidate_result_entry_v6.tsv"} \\
  --endpoint-plan {V5_DIR / "preregistered_endpoint_plan_v5.tsv"} \\
  --out-dir {OUT_DIR / "interpreted_results_v6"}
```

## Claim boundary

The packet is an assay execution and interpretation scaffold. Positive controls remain QC only, source-overlap candidates remain mechanistic/support only, and patient-specific clinical selection remains locked.
"""


def make_cover_sheet(summary: dict) -> str:
    return f"""# CROSS-Neo v6 lab handoff cover sheet

Generated: {summary["generated_at"]}

This packet contains the files needed to run the v5 assay plan as a traceable execution package.

Core files:

- `peptide_hla_reagent_order_manifest_v6.tsv`: peptide/HLA/control/provenance order and audit lines.
- `plate_96_execution_map_v6.tsv`: 96-well map with blinded well IDs and result-entry fields.
- `candidate_result_entry_v6.tsv`: candidate-level PASS/FAIL sheet consumed by the interpreter.
- `well_result_entry_v6.tsv`: well-level raw/normalized signal entry sheet.
- `decision_worksheet_v6.tsv`: endpoint-level manual decision worksheet.
- `STATISTICAL_ANALYSIS_PLAN_V6.md`: frozen analysis rules.

Do not edit endpoint thresholds after assay results are visible. If a candidate is excluded, keep the row and mark `EXCLUDE` with an explicit reason.
"""


def make_html(summary: dict, manifest: pd.DataFrame, endpoint: pd.DataFrame, audit: pd.DataFrame, figures: list[Path]) -> Path:
    HUB_DIR.mkdir(parents=True, exist_ok=True)
    HUB_ASSET_DIR.mkdir(parents=True, exist_ok=True)
    LIVE_ASSET_DIR.mkdir(parents=True, exist_ok=True)
    for fig in figures:
        shutil.copy2(fig, HUB_ASSET_DIR / fig.name)
        shutil.copy2(fig, LIVE_ASSET_DIR / fig.name)
    fig_cards = "\n".join(
        f'<figure><img src="assets/cross_neo_execution_v6/{fig.name}" alt="{fig.stem}"><figcaption>{fig.stem}</figcaption></figure>'
        for fig in figures
    )
    html = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>CROSS-Neo execution packet v6</title>
<style>
:root {{ color-scheme: dark; --bg:#0b1020; --panel:#111827; --ink:#e5e7eb; --muted:#9ca3af; --line:#334155; --gold:#f59e0b; --green:#22c55e; --blue:#38bdf8; }}
* {{ box-sizing:border-box; }}
body {{ margin:0; background:var(--bg); color:var(--ink); font-family:Inter, Arial, sans-serif; line-height:1.5; }}
header {{ padding:54px clamp(22px,5vw,72px) 30px; border-bottom:1px solid var(--line); background:#0f172a; }}
.kicker {{ color:var(--gold); text-transform:uppercase; letter-spacing:.08em; font-size:12px; font-weight:800; }}
h1 {{ font-family:Georgia, serif; font-size:clamp(36px,6vw,74px); line-height:1.0; margin:10px 0 14px; letter-spacing:0; }}
.lead {{ max-width:980px; color:#cbd5e1; font-size:18px; }}
.stats {{ display:grid; grid-template-columns:repeat(6,minmax(120px,1fr)); gap:12px; margin-top:26px; }}
.stat {{ border:1px solid var(--line); padding:14px; background:#0b1222; border-radius:8px; }}
.stat b {{ display:block; font-size:24px; color:white; }}
.stat span {{ color:var(--muted); font-size:12px; }}
main {{ display:grid; grid-template-columns:260px 1fr; gap:28px; padding:28px clamp(18px,4vw,56px) 56px; }}
nav {{ position:sticky; top:16px; align-self:start; border:1px solid var(--line); border-radius:8px; padding:16px; background:#0f172a; }}
nav a {{ display:block; color:#cbd5e1; text-decoration:none; padding:8px 0; border-bottom:1px solid #1f2937; }}
section {{ margin-bottom:30px; }}
h2 {{ font-size:24px; margin:0 0 12px; }}
.num {{ color:var(--gold); font-weight:800; margin-right:8px; }}
.data-table {{ width:100%; border-collapse:collapse; font-size:13px; }}
.data-table th,.data-table td {{ border-bottom:1px solid #243044; padding:8px 10px; text-align:left; vertical-align:top; }}
.data-table th {{ color:#93c5fd; background:#111827; position:sticky; top:0; }}
.table-wrap {{ overflow-x:auto; border:1px solid var(--line); border-radius:8px; }}
.figgrid {{ display:grid; grid-template-columns:repeat(2,minmax(260px,1fr)); gap:16px; }}
figure {{ margin:0; border:1px solid var(--line); border-radius:8px; background:#0f172a; padding:10px; }}
img {{ width:100%; height:auto; display:block; }}
figcaption {{ color:var(--muted); font-size:12px; margin-top:6px; }}
.warn {{ border-left:4px solid var(--gold); padding:12px 14px; background:#1f2937; color:#e5e7eb; border-radius:6px; }}
@media (max-width:900px) {{ main {{ grid-template-columns:1fr; }} nav {{ position:relative; top:0; }} .stats {{ grid-template-columns:repeat(2,1fr); }} .figgrid {{ grid-template-columns:1fr; }} }}
</style>
</head>
<body>
<header>
  <div class="kicker">CROSS-Neo v6 execution packet</div>
  <h1>Assay handoff with frozen unlock rules</h1>
  <p class="lead">The v5 preregistered plan is now converted into order lines, a 96-well execution map, candidate/well result sheets, audit trail, and an automatic endpoint interpreter.</p>
  <div class="stats">
    <div class="stat"><b>{summary["candidate_count"]}</b><span>candidates</span></div>
    <div class="stat"><b>{summary["well_count"]}</b><span>wells</span></div>
    <div class="stat"><b>{summary["order_line_count"]}</b><span>order/audit lines</span></div>
    <div class="stat"><b>{summary["endpoint_count"]}</b><span>endpoints</span></div>
    <div class="stat"><b>{summary["confirmatory_null_risk_sum"]:.3f}</b><span>confirmatory null-risk sum</span></div>
    <div class="stat"><b>{summary["dry_run_pending_endpoints"]}</b><span>dry-run pending endpoints</span></div>
  </div>
</header>
<main>
<nav>
  <a href="#tldr">TL;DR</a>
  <a href="#endpoint">Frozen endpoints</a>
  <a href="#manifest">Execution manifest</a>
  <a href="#audit">Reviewer safeguards</a>
  <a href="#figures">Figures</a>
  <a href="#files">Files</a>
</nav>
<div>
<section id="tldr">
  <h2><span class="num">01</span>TL;DR</h2>
  <p class="warn">No endpoint threshold was tuned in v6. The packet freezes v5 rules and makes the handoff auditable from reagent/order line to candidate PASS/FAIL to endpoint unlock.</p>
</section>
<section id="endpoint">
  <h2><span class="num">02</span>Frozen endpoints</h2>
  <div class="table-wrap">{html_table(endpoint, ["plate_v4_arm","n_candidates","threshold_successes","confirmatory_threshold_successes","main_text_claim_tier","confirmatory_false_unlock_risk_under_null","confirmatory_power_at_mean_prior"], 10)}</div>
</section>
<section id="manifest">
  <h2><span class="num">03</span>Top execution manifest</h2>
  <div class="table-wrap">{html_table(manifest.sort_values("order_priority_score", ascending=False), ["plate_v4_arm","candidate_id","peptide","hla_allele_4digit","order_batch","order_priority_score","candidate_success_definition"], 14)}</div>
</section>
<section id="audit">
  <h2><span class="num">04</span>Reviewer safeguards</h2>
  <div class="table-wrap">{html_table(audit, ["audit_item","source_file","required_before_unlock","failure_action"], 20)}</div>
</section>
<section id="figures">
  <h2><span class="num">05</span>Figures</h2>
  <div class="figgrid">{fig_cards}</div>
</section>
<section id="files">
  <h2><span class="num">06</span>Sources + paths</h2>
  <p>Output directory: <code>{OUT_DIR}</code></p>
  <p>Primary zip: <code>{OUT_DIR / "cross_neo_execution_packet_v6_2026_05_10.zip"}</code></p>
  <p>Interpreter: <code>{INTERPRETER}</code></p>
  <p>Live page copied to: <code>{LIVE_HUB_DIR / "cross_neo_execution_packet_v6.html"}</code></p>
</section>
</div>
</main>
</body>
</html>
"""
    html_path = HUB_DIR / "cross_neo_execution_packet_v6.html"
    html_path.write_text(html, encoding="utf-8")
    LIVE_HUB_DIR.mkdir(parents=True, exist_ok=True)
    shutil.copy2(html_path, LIVE_HUB_DIR / html_path.name)
    return html_path


def run_interpreter_dry_run(candidate_results: Path) -> tuple[Path, dict]:
    dry_dir = OUT_DIR / "interpreter_dry_run_empty_results_v6"
    cmd = [
        "python",
        str(INTERPRETER),
        "--results",
        str(candidate_results),
        "--endpoint-plan",
        str(ENDPOINT_V5),
        "--out-dir",
        str(dry_dir),
    ]
    completed = subprocess.run(cmd, cwd=ROOT, check=True, capture_output=True, text=True)
    summary = json.loads(completed.stdout)
    return dry_dir, summary


def run_interpreter_confirmatory_smoke(candidate_results: Path, endpoint: pd.DataFrame) -> tuple[Path, Path, dict]:
    smoke_results = OUT_DIR / "interpreter_smoke_confirmatory_threshold_results_v6.tsv"
    df = read_tsv(candidate_results)
    df["candidate_success_call"] = "FAIL"
    for _, row in endpoint.iterrows():
        arm = row["plate_v4_arm"]
        threshold = int(row["confirmatory_threshold_successes"])
        idx = df.index[df["plate_v4_arm"].eq(arm)].tolist()[:threshold]
        df.loc[idx, "candidate_success_call"] = "PASS"
    df.to_csv(smoke_results, sep="\t", index=False)

    smoke_dir = OUT_DIR / "interpreter_smoke_confirmatory_threshold_calls_v6"
    cmd = [
        "python",
        str(INTERPRETER),
        "--results",
        str(smoke_results),
        "--endpoint-plan",
        str(ENDPOINT_V5),
        "--out-dir",
        str(smoke_dir),
    ]
    completed = subprocess.run(cmd, cwd=ROOT, check=True, capture_output=True, text=True)
    summary = json.loads(completed.stdout)
    return smoke_results, smoke_dir, summary


def make_zip(paths: list[Path], figures: list[Path]) -> Path:
    zip_path = OUT_DIR / "cross_neo_execution_packet_v6_2026_05_10.zip"
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for path in paths + figures:
            if path.exists() and path.is_file():
                zf.write(path, arcname=path.relative_to(OUT_DIR))
        for directory in OUT_DIR.glob("interpreter_*_v6"):
            if not directory.is_dir():
                continue
            for path in directory.glob("*"):
                if path.is_file():
                    zf.write(path, arcname=path.relative_to(OUT_DIR))
    return zip_path


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    plate = read_tsv(PLATE_V4)
    endpoint = read_tsv(ENDPOINT_V5)
    wellmap = read_tsv(WELLMAP_V5)
    read_tsv(RESULT_TEMPLATE_V5)
    ladder = read_tsv(LADDER_V5)
    safeguards = read_tsv(SAFEGUARDS_V5)

    manifest = build_candidate_manifest(plate, endpoint)
    execution_wellmap = build_execution_wellmap(wellmap, manifest)
    order_manifest = build_order_manifest(plate, execution_wellmap, manifest)
    candidate_results = build_candidate_result_entry(manifest)
    well_results = build_well_result_entry(execution_wellmap)
    decision = build_decision_worksheet(endpoint)
    audit = build_audit_trail(endpoint, safeguards)

    paths = [
        write_tsv(manifest, "execution_candidate_manifest_v6.tsv"),
        write_tsv(order_manifest, "peptide_hla_reagent_order_manifest_v6.tsv"),
        write_tsv(execution_wellmap, "plate_96_execution_map_v6.tsv"),
        write_tsv(candidate_results, "candidate_result_entry_v6.tsv"),
        write_tsv(well_results, "well_result_entry_v6.tsv"),
        write_tsv(decision, "decision_worksheet_v6.tsv"),
        write_tsv(audit, "reviewer_audit_trail_v6.tsv"),
        write_tsv(ladder, "claim_unlock_ladder_v6_inherited_from_v5.tsv"),
    ]

    dry_dir, dry_summary = run_interpreter_dry_run(paths[3])
    smoke_results, smoke_dir, smoke_summary = run_interpreter_confirmatory_smoke(paths[3], endpoint)
    paths.append(smoke_results)

    summary = {
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "input_plate_v4": str(PLATE_V4),
        "input_endpoint_v5": str(ENDPOINT_V5),
        "input_wellmap_v5": str(WELLMAP_V5),
        "candidate_count": int(len(manifest)),
        "well_count": int(len(execution_wellmap)),
        "order_line_count": int(len(order_manifest)),
        "endpoint_count": int(len(endpoint)),
        "confirmatory_null_risk_sum": float(endpoint["confirmatory_false_unlock_risk_under_null"].sum()),
        "exploratory_null_risk_sum": float(endpoint["false_unlock_risk_under_null"].sum()),
        "dry_run_pending_endpoints": int(dry_summary["pending_endpoints"]),
        "dry_run_dir": str(dry_dir),
        "confirmatory_smoke_dir": str(smoke_dir),
        "confirmatory_smoke_exploratory_unlocked": int(smoke_summary["exploratory_unlocked"]),
        "confirmatory_smoke_confirmatory_unlocked": int(smoke_summary["confirmatory_unlocked"]),
    }

    sap_path = OUT_DIR / "STATISTICAL_ANALYSIS_PLAN_V6.md"
    report_path = OUT_DIR / "EXECUTION_PACKET_V6_REPORT_KR.md"
    cover_path = OUT_DIR / "LAB_HANDOFF_COVER_SHEET_V6.md"
    command_path = OUT_DIR / "interpreter_command_v6.txt"
    summary_path = OUT_DIR / "execution_packet_v6_summary.json"
    sap_path.write_text(make_sap(summary, endpoint, audit), encoding="utf-8")
    report_path.write_text(make_report(summary, manifest, endpoint, order_manifest), encoding="utf-8")
    cover_path.write_text(make_cover_sheet(summary), encoding="utf-8")
    command_path.write_text(
        "\n".join(
            [
                "python scripts/interpret_cross_neo_assay_results_v6.py \\",
                f"  --results {OUT_DIR / 'candidate_result_entry_v6.tsv'} \\",
                f"  --endpoint-plan {ENDPOINT_V5} \\",
                f"  --out-dir {OUT_DIR / 'interpreted_results_v6'}",
                "",
            ]
        ),
        encoding="utf-8",
    )
    paths.extend([sap_path, report_path, cover_path, command_path, summary_path])

    figures = make_figures(manifest, order_manifest, execution_wellmap, endpoint)
    html_path = make_html(summary, manifest, endpoint, audit, figures)
    expected_zip_path = OUT_DIR / "cross_neo_execution_packet_v6_2026_05_10.zip"
    summary["html_path"] = str(html_path)
    summary["live_html_path"] = str(LIVE_HUB_DIR / html_path.name)
    summary["zip_path"] = str(expected_zip_path)
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    zip_path = make_zip(paths, figures)

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
