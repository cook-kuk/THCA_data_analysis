#!/usr/bin/env python3
"""DL-first candidate funnel for CROSS-Neo neoantigen prioritization.

This is intentionally ordered from cheap to expensive:

1. Main CROSS-Neo / pMHC model ensemble screen.
2. Bayesian shrinkage and dropout-like perturbation uncertainty screen.
3. TCR-aware rescue/review only where TCR evidence exists.
4. Structure/MD escalation list only after the cheap gates.

MD evidence is reported for context but is not used to rescue candidates that
fail the cheap DL gates. This keeps expensive simulation as a late audit layer.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


REPO = Path(__file__).resolve().parents[3]
CROSS = REPO / "project/results/cross_neo_v2_sota_sprint_2026_05_09"
MD_OUT = REPO / "project/results/cross_neo_md_audit_2026_05_10"
UNC = MD_OUT / "small_dataset_uncertainty"
OUT = MD_OUT / "dl_first_funnel"
FIG = MD_OUT / "figures"


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t") if path.exists() else pd.DataFrame()


def norm_hla(x: object) -> str:
    if pd.isna(x):
        return ""
    s = str(x).strip().upper()
    return s if s.startswith("HLA-") else f"HLA-{s}"


def add_key(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["peptide_norm"] = out["peptide"].astype(str).str.upper().str.strip()
    out["hla_norm"] = out["hla_4digit"].map(norm_hla)
    out["pmhc_key"] = out["peptide_norm"] + "|" + out["hla_norm"]
    return out


def num(df: pd.DataFrame, col: str, default: float = 0.0) -> pd.Series:
    if col not in df.columns:
        return pd.Series(default, index=df.index, dtype=float)
    return pd.to_numeric(df[col], errors="coerce").fillna(default)


def txt(df: pd.DataFrame, col: str) -> pd.Series:
    if col not in df.columns:
        return pd.Series("", index=df.index, dtype=object)
    return df[col].fillna("").astype(str)


def load_inputs() -> pd.DataFrame:
    base = read_tsv(CROSS / "tcr_extension/tcr_wetlab_candidate_prioritization_unique_pmhc.tsv")
    if base.empty:
        raise SystemExit("Missing tcr_wetlab_candidate_prioritization_unique_pmhc.tsv")
    base = add_key(base)

    pieces = [
        read_tsv(UNC / "model_score_aggregation.tsv"),
        read_tsv(UNC / "bayesian_posterior_scores.tsv"),
        read_tsv(UNC / "dropout_perturbation_summary.tsv"),
    ]
    out = base
    for piece in pieces:
        if not piece.empty and "row_id" in piece.columns:
            out = out.merge(piece, on="row_id", how="left")

    ultra = read_tsv(MD_OUT / "ultra_priority/ultra_wetlab_priority_candidates.tsv")
    if not ultra.empty:
        ultra = add_key(ultra)
        keep = [
            "pmhc_key",
            "recommendation_tier",
            "ultra_priority_score",
            "md_label",
            "md_score",
            "md_structural_score",
            "control_readiness_score",
            "live_completion_fraction",
            "wt_status",
            "anchor_preserved_decoy",
            "why",
        ]
        out = out.merge(ultra[[c for c in keep if c in ultra.columns]].drop_duplicates("pmhc_key"), on="pmhc_key", how="left")

    controls = read_tsv(MD_OUT / "counterfactual_design/candidate_control_sequences.tsv")
    if not controls.empty:
        controls = add_key(controls)
        keep = ["pmhc_key", "tentative_wt_sequence", "wt_status", "anchor_preserved_decoy"]
        present = [c for c in keep if c in controls.columns]
        out = out.merge(controls[present].drop_duplicates("pmhc_key"), on="pmhc_key", how="left", suffixes=("", "_control"))

    return out


def score_and_gate(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()

    out["peptide_length"] = out["peptide_norm"].str.len()
    out["class_i_like_hla"] = out["hla_norm"].str.match(r"^HLA-[ABC]\*").fillna(False)
    out["supported_length_main_dl"] = np.where(
        out["class_i_like_hla"],
        out["peptide_length"].between(8, 11),
        out["peptide_length"].between(9, 25),
    )

    rank_support = (1.0 - num(out, "rank_pct_q25", 1.0)).clip(0, 1)
    out["main_dl_score"] = (
        0.45 * num(out, "ensemble_median")
        + 0.20 * num(out, "ensemble_q75")
        + 0.15 * num(out, "vote_gt_050")
        + 0.10 * rank_support
        + 0.10 * num(out, "pmhc_score_mean")
    ).clip(0, 1)

    out["dl_upper_bound"] = num(out, "ensemble_q95").combine(num(out, "bayes_q95"), max)
    out["dl_lower_bound"] = num(out, "bayes_q05")
    out["dl_uncertainty_width"] = num(out, "bayes_q95") - num(out, "bayes_q05")
    out["dropout_unstable"] = (num(out, "dropout_sensitivity") >= 0.12) | (num(out, "perturb_width_90") >= 0.45)

    out["stage1_main_dl_broad_pass"] = out["supported_length_main_dl"] & (
        (out["main_dl_score"] >= 0.45) | (num(out, "ensemble_q95") >= 0.75) | (num(out, "pmhc_score_mean") >= 0.65)
    )
    out["stage2_uncertainty_pass"] = out["stage1_main_dl_broad_pass"] & (
        (
            (num(out, "bayes_mean") >= 0.50)
            & (num(out, "perturb_prob_gt_050") >= 0.50)
            & (~out["dropout_unstable"])
        )
        | (num(out, "bayes_q95") >= 0.70)
        | (out["main_dl_score"] >= 0.70)
    )
    out["stage3_robust_dl_pass"] = out["stage2_uncertainty_pass"] & (~out["dropout_unstable"]) & (
        (num(out, "bayes_mean") >= 0.55)
        | (num(out, "perturb_prob_gt_050") >= 0.70)
        | (out["main_dl_score"] >= 0.75)
    )

    out["has_tcr_evidence"] = num(out, "tcr_evidence_count") > 0
    out["has_paired_tcr_evidence"] = num(out, "paired_tcr_evidence_count") > 0
    out["source_compatible_tcr"] = (
        out["has_tcr_evidence"]
        & (num(out, "cancer_context_evidence_count") > 0)
        & (num(out, "pathogen_context_evidence_count") == 0)
    )
    out["paired_source_compatible_tcr"] = out["source_compatible_tcr"] & out["has_paired_tcr_evidence"]
    out["tcr_rescue_signal"] = (
        out["supported_length_main_dl"]
        & out["paired_source_compatible_tcr"]
        & (num(out, "tcr_augmented_score_mean") >= 0.85)
        & (num(out, "best_tcr_augmented_score") >= 0.85)
        & (num(out, "bayes_q95") >= 0.50)
    )
    out["tcr_branch_support"] = (num(out, "tcr_augmented_score_mean") >= 0.65) | (num(out, "best_tcr_augmented_score") >= 0.85)
    out["tcr_branch_discordant_low"] = out["source_compatible_tcr"] & (~out["tcr_branch_support"])

    out["wt_or_decoy_ready"] = txt(out, "wt_status").str.len().gt(0) | txt(out, "anchor_preserved_decoy").str.len().gt(0)
    out["md_context_only"] = txt(out, "md_label")

    decisions: list[str] = []
    next_actions: list[str] = []
    reasons: list[str] = []
    for _, r in out.iterrows():
        if not bool(r["supported_length_main_dl"]):
            decisions.append("DL_CULL_UNSUPPORTED_PEPTIDE_LENGTH")
            next_actions.append("do_not_send_to_structure_or_md")
            reasons.append("peptide length is outside the cheap model support range for this HLA class")
        elif not bool(r["stage1_main_dl_broad_pass"]):
            decisions.append("DL_CULL_LOW_MAIN_MODEL_SCORE")
            next_actions.append("hold_unless_new_patient_or_tcr_data")
            reasons.append("main DL ensemble does not clear the broad cheap screen")
        elif not bool(r["stage2_uncertainty_pass"]):
            if bool(r["tcr_rescue_signal"]):
                decisions.append("TCR_RESCUE_REVIEW_AFTER_WEAK_MAIN_DL")
                next_actions.append("curate_tcr_wt_then_structure_before_md")
                reasons.append("main DL is uncertain, but paired source-compatible TCR model signal is strong")
            else:
                decisions.append("DL_CULL_UNCERTAIN_OR_LOW_POSTERIOR")
                next_actions.append("hold_after_bayesian_dropout_screen")
                reasons.append("Bayesian/dropout screen does not support escalation")
        elif bool(r["dropout_unstable"]) and not bool(r["paired_source_compatible_tcr"]):
            decisions.append("DL_CULL_DROPOUT_UNSTABLE")
            next_actions.append("hold_until_more_data_or_model_agreement")
            reasons.append("score moves too much under dropout/weight perturbation without paired TCR support")
        elif bool(r["stage3_robust_dl_pass"]) and bool(r["paired_source_compatible_tcr"]) and not bool(r["tcr_branch_support"]):
            decisions.append("DL_TCR_DISCORDANT_AUDIT_NOT_POSITIVE")
            next_actions.append("use_as_false_positive_audit_or_negative_control_not_md_positive")
            reasons.append("main DL is high, but source-compatible TCR-aware branch is low")
        elif bool(r["stage3_robust_dl_pass"]) and bool(r["source_compatible_tcr"]) and not bool(r["tcr_branch_support"]):
            decisions.append("DL_TCR_DISCORDANT_AUDIT_NOT_POSITIVE")
            next_actions.append("use_as_false_positive_audit_or_negative_control_not_md_positive")
            reasons.append("main DL is high, but available TCR-aware evidence is discordant")
        elif bool(r["stage3_robust_dl_pass"]) and bool(r["paired_source_compatible_tcr"]) and bool(r["tcr_branch_support"]):
            decisions.append("DL_PRIMARY_SURVIVOR_WITH_PAIRED_TCR")
            next_actions.append("prepare_wt_decoy_structure_then_md_if_controls_ready")
            reasons.append("cheap main DL and TCR-aware branch both support escalation with paired source-compatible TCR evidence")
        elif bool(r["stage3_robust_dl_pass"]) and bool(r["source_compatible_tcr"]) and bool(r["tcr_branch_support"]):
            decisions.append("DL_PRIMARY_SURVIVOR_TCR_CONTEXT_ONLY")
            next_actions.append("curate_paired_tcr_or_template_before_md")
            reasons.append("cheap main DL and TCR-aware branch both support review, but paired TCR is not available")
        elif bool(r["stage3_robust_dl_pass"]):
            decisions.append("DL_PRIMARY_SURVIVOR_MODEL_ONLY")
            next_actions.append("curate_tcr_wt_or_use_as_pmhc_only_candidate")
            reasons.append("cheap DL signal survives, but TCR/WT evidence is not yet enough for expensive MD")
        elif bool(r["tcr_rescue_signal"]):
            decisions.append("TCR_RESCUE_REVIEW_AFTER_WEAK_MAIN_DL")
            next_actions.append("curate_tcr_wt_then_structure_before_md")
            reasons.append("TCR-aware cheap branch rescues a main-model-uncertain candidate")
        else:
            decisions.append("DL_HOLD_MODEL_SIGNAL_NOT_ROBUST")
            next_actions.append("hold_for_more_model_or_registry_evidence")
            reasons.append("some model signal exists, but it is not robust enough for the next expensive layer")

    out["dl_first_decision"] = decisions
    out["next_action"] = next_actions
    out["dl_first_reason"] = reasons

    out["dl_first_priority_score"] = (
        0.38 * out["main_dl_score"]
        + 0.22 * num(out, "bayes_mean")
        + 0.16 * num(out, "perturb_prob_gt_050")
        + 0.10 * out["paired_source_compatible_tcr"].astype(float)
        + 0.06 * out["source_compatible_tcr"].astype(float)
        + 0.08 * out["tcr_branch_support"].astype(float)
        + 0.05 * out["wt_or_decoy_ready"].astype(float)
        - 0.10 * out["dropout_unstable"].astype(float)
        - 0.12 * out["tcr_branch_discordant_low"].astype(float)
        - 0.08 * (num(out, "pathogen_context_evidence_count") > num(out, "cancer_context_evidence_count")).astype(float)
    ).clip(0, 1)

    order = {
        "DL_PRIMARY_SURVIVOR_WITH_PAIRED_TCR": 0,
        "TCR_RESCUE_REVIEW_AFTER_WEAK_MAIN_DL": 1,
        "DL_PRIMARY_SURVIVOR_TCR_CONTEXT_ONLY": 2,
        "DL_PRIMARY_SURVIVOR_MODEL_ONLY": 3,
        "DL_TCR_DISCORDANT_AUDIT_NOT_POSITIVE": 4,
        "DL_HOLD_MODEL_SIGNAL_NOT_ROBUST": 5,
        "DL_CULL_UNCERTAIN_OR_LOW_POSTERIOR": 6,
        "DL_CULL_DROPOUT_UNSTABLE": 7,
        "DL_CULL_LOW_MAIN_MODEL_SCORE": 8,
        "DL_CULL_UNSUPPORTED_PEPTIDE_LENGTH": 9,
    }
    out["dl_first_order"] = out["dl_first_decision"].map(order).fillna(99).astype(int)
    return out.sort_values(["dl_first_order", "dl_first_priority_score"], ascending=[True, False])


def write_tables(df: pd.DataFrame) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    cols = [
        "row_id",
        "peptide",
        "hla_4digit",
        "peptide_length",
        "label_binary",
        "dl_first_decision",
        "next_action",
        "dl_first_priority_score",
        "main_dl_score",
        "ensemble_median",
        "ensemble_q95",
        "vote_gt_050",
        "bayes_mean",
        "bayes_q05",
        "bayes_q95",
        "perturb_prob_gt_050",
        "dropout_sensitivity",
        "tcr_augmented_score_mean",
        "best_tcr_augmented_score",
        "tcr_branch_support",
        "tcr_branch_discordant_low",
        "tcr_evidence_count",
        "paired_tcr_evidence_count",
        "cancer_context_evidence_count",
        "pathogen_context_evidence_count",
        "wt_status",
        "anchor_preserved_decoy",
        "md_label",
        "live_completion_fraction",
        "dl_first_reason",
    ]
    present = [c for c in cols if c in df.columns]
    df[present].to_csv(OUT / "dl_first_candidate_funnel.tsv", sep="\t", index=False)

    survivor_mask = df["dl_first_decision"].isin(
        [
            "DL_PRIMARY_SURVIVOR_WITH_PAIRED_TCR",
            "TCR_RESCUE_REVIEW_AFTER_WEAK_MAIN_DL",
            "DL_PRIMARY_SURVIVOR_TCR_CONTEXT_ONLY",
            "DL_PRIMARY_SURVIVOR_MODEL_ONLY",
        ]
    )
    df.loc[survivor_mask, present].to_csv(OUT / "dl_first_survivors_for_tcr_structure_review.tsv", sep="\t", index=False)

    md_mask = df["dl_first_decision"].isin(["DL_PRIMARY_SURVIVOR_WITH_PAIRED_TCR", "TCR_RESCUE_REVIEW_AFTER_WEAK_MAIN_DL"])
    df.loc[md_mask, present].to_csv(OUT / "dl_first_md_escalation_candidates.tsv", sep="\t", index=False)

    wetlab_mask = df["dl_first_decision"].eq("DL_PRIMARY_SURVIVOR_WITH_PAIRED_TCR") & df["wt_or_decoy_ready"]
    df.loc[wetlab_mask, present].to_csv(OUT / "dl_first_wetlab_shortlist.tsv", sep="\t", index=False)

    df.loc[df["dl_first_decision"].str.startswith("DL_CULL"), present].to_csv(OUT / "dl_first_culled_candidates.tsv", sep="\t", index=False)

    counts = df["dl_first_decision"].value_counts().rename_axis("decision").reset_index(name="n")
    counts.to_csv(OUT / "dl_first_decision_counts.tsv", sep="\t", index=False)

    stage_counts = pd.DataFrame(
        [
            {"stage": "00_all_candidates", "n": int(len(df))},
            {"stage": "01_supported_peptide_length", "n": int(df["supported_length_main_dl"].sum())},
            {"stage": "02_main_dl_broad_pass", "n": int(df["stage1_main_dl_broad_pass"].sum())},
            {"stage": "03_bayesian_dropout_pass", "n": int(df["stage2_uncertainty_pass"].sum())},
            {"stage": "04_robust_main_dl_pass", "n": int(df["stage3_robust_dl_pass"].sum())},
            {"stage": "05_source_compatible_tcr_after_dl", "n": int((df["stage3_robust_dl_pass"] & df["source_compatible_tcr"]).sum())},
            {"stage": "06_paired_tcr_after_dl", "n": int((df["stage3_robust_dl_pass"] & df["paired_source_compatible_tcr"]).sum())},
            {"stage": "07_tcr_branch_supported_after_dl", "n": int((df["stage3_robust_dl_pass"] & df["source_compatible_tcr"] & df["tcr_branch_support"]).sum())},
            {"stage": "08_tcr_rescue_reviews", "n": int((df["dl_first_decision"] == "TCR_RESCUE_REVIEW_AFTER_WEAK_MAIN_DL").sum())},
            {"stage": "09_md_escalation_after_cheap_gates", "n": int(md_mask.sum())},
            {"stage": "10_wetlab_shortlist_after_dl", "n": int(wetlab_mask.sum())},
        ]
    )
    stage_counts.to_csv(OUT / "dl_first_stage_counts.tsv", sep="\t", index=False)


def plot(df: pd.DataFrame) -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    stage_counts = read_tsv(OUT / "dl_first_stage_counts.tsv")
    fig, ax = plt.subplots(figsize=(10.8, 5.6))
    ax.plot(stage_counts["stage"], stage_counts["n"], marker="o", color="#2f80ed", lw=2.2)
    ax.set_ylabel("Candidate count")
    ax.set_title("DL-first neoantigen funnel before structure/MD")
    ax.tick_params(axis="x", rotation=35, labelsize=8)
    ax.grid(axis="y", alpha=0.25)
    fig.savefig(FIG / "fig_md18_dl_first_candidate_funnel.png", dpi=220, bbox_inches="tight")
    fig.savefig(FIG / "fig_md18_dl_first_candidate_funnel.pdf", bbox_inches="tight")
    plt.close(fig)

    top = df.head(30).copy()
    label = top["peptide"].astype(str) + "\n" + top["hla_4digit"].astype(str)
    fig, ax = plt.subplots(figsize=(11.5, 7.0))
    ax.scatter(top["main_dl_score"], label, s=55, color="#2f80ed", label="main DL score")
    ax.errorbar(
        top["bayes_mean"],
        label,
        xerr=[top["bayes_mean"] - top["bayes_q05"], top["bayes_q95"] - top["bayes_mean"]],
        fmt="none",
        ecolor="#7cb7ff",
        alpha=0.8,
        label="Bayesian 90%",
    )
    ax.scatter(top["tcr_augmented_score_mean"], label, s=38, marker="x", color="#e76f51", label="TCR-aware cheap branch")
    ax.set_xlim(0, 1.02)
    ax.set_xlabel("Score")
    ax.set_title("Cheap model evidence first; TCR branch is rescue/review")
    ax.legend(frameon=False, loc="lower right")
    ax.grid(axis="x", alpha=0.25)
    fig.savefig(FIG / "fig_md19_dl_tcr_rescue_uncertainty.png", dpi=220, bbox_inches="tight")
    fig.savefig(FIG / "fig_md19_dl_tcr_rescue_uncertainty.pdf", bbox_inches="tight")
    plt.close(fig)


def write_report(df: pd.DataFrame) -> None:
    counts = read_tsv(OUT / "dl_first_decision_counts.tsv")
    stages = read_tsv(OUT / "dl_first_stage_counts.tsv")
    top_cols = [
        "row_id",
        "peptide",
        "hla_4digit",
        "dl_first_decision",
        "next_action",
        "dl_first_priority_score",
        "main_dl_score",
        "bayes_mean",
        "bayes_q05",
        "bayes_q95",
        "perturb_prob_gt_050",
        "tcr_augmented_score_mean",
        "paired_tcr_evidence_count",
        "md_label",
    ]
    top = df[[c for c in top_cols if c in df.columns]].head(25)
    lines = [
        "# DL-First Candidate Funnel",
        "",
        "## Intent",
        "",
        "This funnel runs cheap CROSS-Neo model evidence before expensive structure or MD. MD is a late audit/escalation layer, not an early rescue signal.",
        "",
        "## Stage Counts",
        "",
        stages.to_markdown(index=False),
        "",
        "## Decision Counts",
        "",
        counts.to_markdown(index=False),
        "",
        "## Highest Priority After Cheap Gates",
        "",
        top.to_markdown(index=False),
        "",
        "## Operating Rule",
        "",
        "1. Use the main pMHC/CROSS-Neo ensemble to remove low-scoring or unsupported-length peptides.",
        "2. Use Bayesian shrinkage and dropout-like perturbation to remove candidates whose probability changes too much.",
        "3. Use TCR-aware scores only as an optional rescue/review branch for TCR-available rows.",
        "4. Send only DL survivors or TCR-rescue reviews to WT/decoy structure preparation.",
        "5. Run explicit-solvent MD only on the tiny subset with model support, paired/source-compatible TCR evidence, and controls.",
        "",
        "## Claim Boundary",
        "",
        "These are prioritization scores, not calibrated immunogenicity probabilities. A candidate that survives this funnel is a better wetlab/MD candidate, not a proven positive.",
    ]
    (OUT / "dl_first_funnel_report.md").write_text("\n".join(lines) + "\n")
    (OUT / "dl_first_summary.json").write_text(
        json.dumps(
            {
                "n_candidates": int(len(df)),
                "stage_counts": stages.set_index("stage")["n"].to_dict(),
                "decision_counts": counts.set_index("decision")["n"].to_dict(),
                "top_candidates": top.head(10).to_dict("records"),
            },
            indent=2,
        )
    )


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    df = score_and_gate(load_inputs())
    write_tables(df)
    plot(df)
    write_report(df)
    print("[dl-first-funnel]", (OUT / "dl_first_summary.json").read_text())


if __name__ == "__main__":
    main()
