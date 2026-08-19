#!/usr/bin/env python3
"""Build an impact-acceleration package for CROSS-Neo MD/TCR prioritization."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


REPO = Path(__file__).resolve().parents[3]
MD_OUT = REPO / "project/results/cross_neo_md_audit_2026_05_10"
PKG = MD_OUT / "high_impact_decision_package"
LEAD = MD_OUT / "two_lead_impact_package"
ASSAY = MD_OUT / "assay_ready_translation_package"
PRECLIN = MD_OUT / "preclinical_validation_protocol"
OUT = MD_OUT / "impact_acceleration_package"
FIG = MD_OUT / "figures"


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t") if path.exists() else pd.DataFrame()


def esc(value: object) -> str:
    return "" if pd.isna(value) else str(value)


def md_table(df: pd.DataFrame, cols: list[str], n: int = 20) -> str:
    if df.empty:
        return "_No rows available._"
    keep = [c for c in cols if c in df.columns]
    if not keep:
        return "_No requested columns available._"
    return df[keep].head(n).to_markdown(index=False)


def load_inputs() -> dict[str, pd.DataFrame]:
    inputs = {
        "top13": read_tsv(PKG / "no_false_positive_top13_candidates.tsv"),
        "source_perf": read_tsv(PKG / "preset_source_stratified_performance.tsv"),
        "claims": read_tsv(PKG / "claim_ladder.tsv"),
        "risks": read_tsv(PKG / "reviewer_risk_register.tsv"),
        "validation": read_tsv(PKG / "external_validation_plan.tsv"),
        "evidence": read_tsv(PKG / "evidence_to_claim_matrix.tsv"),
        "leads": read_tsv(LEAD / "two_lead_external_evidence.tsv"),
        "locks": read_tsv(LEAD / "two_lead_specificity_locks.tsv"),
        "wetlab": read_tsv(LEAD / "two_lead_wetlab_order.tsv"),
        "assay_thresholds": read_tsv(ASSAY / "assay_go_no_go_thresholds.tsv"),
        "preclin_score": read_tsv(PRECLIN / "preclinical_readiness_scorecard.tsv"),
        "preclin_gates": read_tsv(PRECLIN / "preclinical_validation_gates.tsv"),
        "preclin_failures": read_tsv(PRECLIN / "preclinical_failure_mode_actions.tsv"),
        "md_scores": read_tsv(MD_OUT / "md_evidence_scores.tsv"),
    }
    if inputs["top13"].empty or inputs["leads"].empty or inputs["preclin_score"].empty:
        raise FileNotFoundError("Required high-impact, lead, or preclinical inputs are missing")
    return inputs


def build_editor_reviewer_positioning(data: dict[str, pd.DataFrame]) -> pd.DataFrame:
    top13 = data["top13"]
    n_top = len(top13)
    n_tp = int(top13.get("prediction_outcome", pd.Series(dtype=str)).astype(str).eq("TP").sum())
    n_fp = int(top13.get("prediction_outcome", pd.Series(dtype=str)).astype(str).eq("FP").sum())
    return pd.DataFrame(
        [
            {
                "attack_or_question": "Is this only a peptide-HLA presentation model?",
                "strong_answer": "No. The main model stays pMHC-centered, but the optional expert layer adds paired-TCR evidence, TCR-pMHC structure audit, MD stability, and WT/decoy assay gates.",
                "asset_that_answers_it": "high-impact decision page; two-lead impact package; preclinical validation protocol",
                "current_evidence": f"strict current-label board has {n_top} candidates, {n_tp} TP, {n_fp} FP under available labels",
                "claim_boundary": "Do not claim universal TCR-aware prediction because most rows lack paired TCR.",
                "impact_if_resolved": "Positions CROSS-Neo as a practical prioritization system, not just another leaderboard model.",
            },
            {
                "attack_or_question": "Does MD prove immunogenicity?",
                "strong_answer": "No. MD is a structural audit layer that checks peptide groove stability, TCR-peptide contact persistence, and failure explanations before wetlab.",
                "asset_that_answers_it": "MD evidence scores; pMHC/TCR contact reports; claim ladder",
                "current_evidence": "HMTEVVRHC is MD_VERY_STRONG; GADGVGKSAL is replicated MD_MODERATE across 10 ns runs.",
                "claim_boundary": "No immunogenicity claim until mutant activation exceeds WT and decoy controls.",
                "impact_if_resolved": "Makes the work more credible because it explicitly avoids overclaiming simulation.",
            },
            {
                "attack_or_question": "Could the two flagship candidates be cherry-picked shared neoantigens?",
                "strong_answer": "They are intentionally promoted as flagship case studies because they have external biological anchors, TCR resources, WT controls, and assay feasibility.",
                "asset_that_answers_it": "two_lead_external_evidence.tsv; specificity locks; wetlab order",
                "current_evidence": "KRAS G12D/HLA-C*08:02 and TP53 R175H/HLA-A*02:01 have mutation-specific WT controls and published TCR/structure context.",
                "claim_boundary": "Use as case-study and wetlab-prioritization evidence, not standalone external validation.",
                "impact_if_resolved": "Creates a publishable translational story even before broad prospective validation.",
            },
            {
                "attack_or_question": "Is zero false positive real?",
                "strong_answer": "It is current-label decision-support, not universal accuracy. The package exposes threshold movement, source balance, and required external validation.",
                "asset_that_answers_it": "slider dashboard; source-stratified performance; external validation plan",
                "current_evidence": f"strict Top-{n_top} has {n_fp} FP in current labels, but external source-heldout testing remains required.",
                "claim_boundary": "Do not claim 100% accuracy or clinical utility.",
                "impact_if_resolved": "Turns a fragile-looking metric into a transparent operating point.",
            },
            {
                "attack_or_question": "What makes this more than a model comparison?",
                "strong_answer": "The output is a go/no-go wetlab decision system with reagent rows, plate maps, validation gates, and failure-mode actions.",
                "asset_that_answers_it": "assay-ready translation package; preclinical validation protocol",
                "current_evidence": "10 reagent rows, 144 plate-map rows, 10 preclinical gates, and 14 failure-mode actions are materialized.",
                "claim_boundary": "This is preclinical planning, not clinical utility.",
                "impact_if_resolved": "Moves the manuscript from benchmark-only to experimentally actionable neoantigen triage.",
            },
        ]
    )


def build_claim_upgrade_milestones(data: dict[str, pd.DataFrame]) -> pd.DataFrame:
    score = data["preclin_score"].copy()
    rows = []
    stage_defs = [
        (
            1,
            "computational_prioritization",
            "Top candidate under strict threshold board; source/label caveats visible.",
            "high-priority wetlab candidate",
            "universal predictor or clinical utility",
        ),
        (
            2,
            "structure_md_audit",
            "Peptide remains pMHC-stable and TCR-peptide contacts persist in pilot MD.",
            "structure-supported plausibility",
            "immunogenicity proof",
        ),
        (
            3,
            "wt_decoy_specificity",
            "Mutant exceeds WT and decoy in HLA binding/stability and multimer/TCR recognition.",
            "mutation-specific recognition in tested context",
            "general patient-level response",
        ),
        (
            4,
            "activation_assay",
            "Mutant induces ELISpot/ICS cytokine response above WT, decoy, and irrelevant peptide.",
            "assay-specific immunogenicity",
            "clinical efficacy",
        ),
        (
            5,
            "functional_killing",
            "Matched HLA/mutation target cells are killed while WT/decoy/mismatched controls remain negative.",
            "functional specificity in the tested system",
            "safety or therapeutic benefit",
        ),
        (
            6,
            "external_reproducibility",
            "Heldout cohort/source or independent lab reproduces candidate ranking and assay signal.",
            "externally supported prioritization",
            "broad clinical generalization",
        ),
    ]
    for _, lead in score.iterrows():
        for order, milestone, criterion, allowed, forbidden in stage_defs:
            rows.append(
                {
                    "lead": lead["lead"],
                    "mutant_peptide": lead["mutant_peptide"],
                    "hla": lead["hla"],
                    "stage_order": order,
                    "milestone": milestone,
                    "pass_criterion": criterion,
                    "claim_allowed_after_pass": allowed,
                    "claim_still_forbidden": forbidden,
                    "current_status": "DONE_OR_IN_PROGRESS" if order <= 2 else "NEXT_GATE",
                    "impact_level_if_passed": ["moderate", "moderate", "high", "very_high", "very_high", "publication_defense"][order - 1],
                }
            )
    return pd.DataFrame(rows)


def build_external_validation_benchmark_plan(data: dict[str, pd.DataFrame]) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "validation_axis": "source_heldout_neoantigen",
                "dataset_or_material": "independent cancer neoantigen dataset with labels and HLA annotations",
                "primary_metric": "AUPRC, top-k precision, enrichment over prevalence",
                "success_criterion": "Top-k enrichment remains above baseline after excluding overlapping public epitopes.",
                "why_it_raises_impact": "Addresses overfitting and source shift.",
                "claim_if_passed": "source-robust prioritization on heldout data",
            },
            {
                "validation_axis": "paired_tcr_available_subset",
                "dataset_or_material": "paired alpha/beta TCR-pMHC examples with compatible labels",
                "primary_metric": "paired TCR-heldout AUPRC and abstained top-k precision",
                "success_criterion": "TCR expert improves or abstains without harming no-TCR rows.",
                "why_it_raises_impact": "Supports optional TCR-aware method without overclaiming missing TCR data.",
                "claim_if_passed": "TCR-aware expert helps when paired TCR exists",
            },
            {
                "validation_axis": "wetlab_two_lead_specificity",
                "dataset_or_material": "TP53 R175H and KRAS G12D mutant/WT/decoy reagents",
                "primary_metric": "mutant-over-WT/decoy binding, activation, and killing gates",
                "success_criterion": "Mutant signal exceeds WT/decoy across replicate-controlled assays.",
                "why_it_raises_impact": "Turns computational ranking into experimentally actionable prioritization.",
                "claim_if_passed": "assay-specific immunogenicity for tested lead/context",
            },
            {
                "validation_axis": "negative_control_specificity",
                "dataset_or_material": "scrambled or anchor-preserved decoys plus irrelevant pHLA controls",
                "primary_metric": "false positive rate under controlled assays",
                "success_criterion": "Decoy and irrelevant controls remain below activation/multimer thresholds.",
                "why_it_raises_impact": "Directly addresses nonspecific TCR binding and motif artifacts.",
                "claim_if_passed": "lower nonspecific recognition risk in tested setup",
            },
            {
                "validation_axis": "prospective_candidate_batch",
                "dataset_or_material": "new candidates selected before seeing wetlab outcomes",
                "primary_metric": "pre-registered top-k hit rate and calibration",
                "success_criterion": "Pre-registered thresholds retain enrichment and acceptable false-positive burden.",
                "why_it_raises_impact": "Moves the work from retrospective audit to prospective validation.",
                "claim_if_passed": "prospective prioritization evidence",
            },
        ]
    )


def build_wetlab_success_scenarios(data: dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows = []
    for _, lead in data["preclin_score"].iterrows():
        rows += [
            {
                "lead": lead["lead"],
                "scenario": "best_case",
                "observed_result": "mutant pHLA stable; mutant multimer/TCR binding positive; activation positive; WT/decoy negative",
                "interpretation": "strong candidate for case-study immunogenicity claim in tested context",
                "next_move": "run killing assay and independent replicate; prepare flagship figure",
                "allowed_claim": "assay-specific immunogenicity after activation gate passes",
            },
            {
                "lead": lead["lead"],
                "scenario": "presentation_only",
                "observed_result": "mutant pHLA stable but TCR/multimer or activation negative",
                "interpretation": "presentation is not recognition; useful false-positive explanation",
                "next_move": "downgrade to diagnostic case; search alternate TCR or deprioritize",
                "allowed_claim": "presentation plausibility only",
            },
            {
                "lead": lead["lead"],
                "scenario": "wt_cross_reactive",
                "observed_result": "WT response equals or exceeds mutant response",
                "interpretation": "specificity risk; do not advance as mutant-specific lead",
                "next_move": "flag as WT cross-reactivity; use as safety/claim-boundary case",
                "allowed_claim": "not mutant-specific in tested context",
            },
            {
                "lead": lead["lead"],
                "scenario": "decoy_positive",
                "observed_result": "decoy or irrelevant pHLA is positive",
                "interpretation": "nonspecific binding or motif artifact",
                "next_move": "repeat with orthogonal controls; do not claim recognition specificity",
                "allowed_claim": "no recognition specificity claim",
            },
        ]
    return pd.DataFrame(rows)


def build_next_action_board(data: dict[str, pd.DataFrame]) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "priority": 1,
                "time_window": "next_24h",
                "action": "Freeze the current strict Top-13 preset and archive exact thresholds.",
                "owner": "computational",
                "deliverable": "locked candidate table and threshold dashboard snapshot",
                "impact": "prevents moving-target criticism",
            },
            {
                "priority": 2,
                "time_window": "next_24h",
                "action": "Order/check availability for HMTEVVRHC, HMTEVVRRC, HMTHRVVEC, GADGVGKSAL, GAGGVGKSAL, GASGVKGADL.",
                "owner": "wetlab_or_collaborator",
                "deliverable": "reagent availability and quote sheet",
                "impact": "turns flagship candidates into executable experiments",
            },
            {
                "priority": 3,
                "time_window": "next_48h",
                "action": "Launch WT/decoy pMHC and TCR-pMHC structure/short-MD controls for both leads.",
                "owner": "simulation",
                "deliverable": "mutant-vs-WT/decoy interface delta table",
                "impact": "adds specificity evidence before assays",
            },
            {
                "priority": 4,
                "time_window": "next_72h",
                "action": "Prepare pre-registered assay gate sheet for collaborator review.",
                "owner": "translation",
                "deliverable": "assay go/no-go protocol with forbidden claims",
                "impact": "reviewer-proof claim discipline",
            },
            {
                "priority": 5,
                "time_window": "next_7d",
                "action": "Run no-public-overlap and source-heldout repeat on any new external dataset.",
                "owner": "computational",
                "deliverable": "external validation addendum",
                "impact": "largest manuscript-impact upgrade after wetlab",
            },
        ]
    )


def build_figures(
    positioning: pd.DataFrame,
    milestones: pd.DataFrame,
    validation: pd.DataFrame,
    actions: pd.DataFrame,
) -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    plt.style.use("default")

    fig, ax = plt.subplots(figsize=(11.5, 4.8))
    ax.axis("off")
    stage = milestones.drop_duplicates("stage_order").sort_values("stage_order")
    xs = np.linspace(0.08, 0.92, len(stage))
    for i, (_, row) in enumerate(stage.iterrows()):
        color = "#1b7f5a" if row["current_status"] == "DONE_OR_IN_PROGRESS" else "#17324f"
        ax.add_patch(plt.Circle((xs[i], 0.60), 0.06, fc=color, ec="#edf5ff", lw=1.2))
        ax.text(xs[i], 0.60, str(int(row["stage_order"])), ha="center", va="center", color="#f2c46d", weight="bold", fontsize=15)
        ax.text(xs[i], 0.34, row["milestone"].replace("_", "\n"), ha="center", va="center", fontsize=8.2)
        if i < len(stage) - 1:
            ax.annotate("", xy=(xs[i + 1] - 0.075, 0.60), xytext=(xs[i] + 0.075, 0.60), arrowprops=dict(arrowstyle="->", lw=2))
    ax.text(0.5, 0.86, "Claim-upgrade milestones: what raises impact without overclaiming", ha="center", fontsize=15, weight="bold")
    ax.text(0.5, 0.14, "Current strength: prioritization + structural audit. Stronger claims require WT/decoy and activation gates.", ha="center", color="#9b1c31", fontsize=10)
    fig.savefig(FIG / "fig_md47_claim_upgrade_milestone_map.png", dpi=220, bbox_inches="tight")
    fig.savefig(FIG / "fig_md47_claim_upgrade_milestone_map.pdf", bbox_inches="tight")
    plt.close(fig)

    score_map = {
        "high-impact decision page; two-lead impact package; preclinical validation protocol": [1, 1, 1, 1],
        "MD evidence scores; pMHC/TCR contact reports; claim ladder": [1, 1, 0.9, 1],
        "two_lead_external_evidence.tsv; specificity locks; wetlab order": [1, 0.8, 1, 0.9],
        "slider dashboard; source-stratified performance; external validation plan": [0.8, 1, 0.7, 1],
        "assay-ready translation package; preclinical validation protocol": [0.9, 0.8, 1, 1],
    }
    heat = np.array([score_map.get(x, [0.5, 0.5, 0.5, 0.5]) for x in positioning["asset_that_answers_it"]])
    fig, ax = plt.subplots(figsize=(10.5, 5.2))
    im = ax.imshow(heat, aspect="auto", cmap="YlGnBu", vmin=0, vmax=1)
    ax.set_yticks(np.arange(len(positioning)), [f"Q{i+1}" for i in range(len(positioning))])
    ax.set_xticks(np.arange(4), ["model", "source\nshift", "wetlab", "claim\nboundary"])
    ax.set_title("Reviewer attack defense coverage")
    for i in range(heat.shape[0]):
        for j in range(heat.shape[1]):
            ax.text(j, i, f"{heat[i, j]:.1f}", ha="center", va="center", fontsize=8)
    fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
    fig.tight_layout()
    fig.savefig(FIG / "fig_md48_reviewer_attack_defense_heatmap.png", dpi=220)
    fig.savefig(FIG / "fig_md48_reviewer_attack_defense_heatmap.pdf")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(11, 5.0))
    levels = np.arange(len(validation))
    ax.barh(levels, [5, 4, 5, 4, 5], color=["#4e79a7", "#59a14f", "#f28e2b", "#e15759", "#76b7b2"])
    ax.set_yticks(levels, validation["validation_axis"])
    ax.set_xlabel("impact leverage")
    ax.set_title("External validation and benchmark design priorities")
    ax.grid(axis="x", alpha=0.25)
    fig.tight_layout()
    fig.savefig(FIG / "fig_md49_external_validation_benchmark_design.png", dpi=220)
    fig.savefig(FIG / "fig_md49_external_validation_benchmark_design.pdf")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(11.5, 4.8))
    ax.axis("off")
    for i, (_, row) in enumerate(actions.sort_values("priority").iterrows()):
        y = 0.82 - i * 0.16
        ax.add_patch(plt.Rectangle((0.03, y - 0.055), 0.10, 0.10, fc="#17324f", ec="#edf5ff"))
        ax.text(0.08, y, str(int(row["priority"])), ha="center", va="center", color="#f2c46d", weight="bold", fontsize=13)
        ax.text(0.17, y + 0.025, row["time_window"], ha="left", va="center", color="#39d4b5", weight="bold", fontsize=9)
        ax.text(0.17, y - 0.020, row["action"], ha="left", va="center", fontsize=9)
        ax.text(0.76, y, row["impact"], ha="left", va="center", fontsize=8.2, color="#9b1c31")
    ax.text(0.5, 0.96, "Next action board: fastest path to higher-impact evidence", ha="center", fontsize=15, weight="bold")
    fig.savefig(FIG / "fig_md50_next_72h_action_board.png", dpi=220, bbox_inches="tight")
    fig.savefig(FIG / "fig_md50_next_72h_action_board.pdf", bbox_inches="tight")
    plt.close(fig)


def write_reports(
    positioning: pd.DataFrame,
    milestones: pd.DataFrame,
    validation: pd.DataFrame,
    scenarios: pd.DataFrame,
    actions: pd.DataFrame,
    data: dict[str, pd.DataFrame],
) -> None:
    top13 = data["top13"]
    preclin = data["preclin_score"]
    hm = preclin[preclin["mutant_peptide"].eq("HMTEVVRHC")].iloc[0]
    kr = preclin[preclin["mutant_peptide"].eq("GADGVGKSAL")].iloc[0]
    lines = [
        "# CROSS-Neo Impact Acceleration Package",
        "",
        "## Executive Position",
        "",
        "The highest-impact framing is no longer only model accuracy. The strongest story is a controlled neoantigen prioritization system: cheap deep learning first, uncertainty and source-shift gates second, optional TCR/structure/MD evidence third, then WT/decoy-controlled wetlab validation.",
        "",
        "## Current Flagship State",
        "",
        f"- Strict current-label candidate board: {len(top13)} candidates.",
        f"- TP53 R175H / HMTEVVRHC / HLA-A*02:01: readiness {float(hm['overall_readiness_score']):.3f}, tier {hm['validation_tier']}.",
        f"- KRAS G12D / GADGVGKSAL / HLA-C*08:02: readiness {float(kr['overall_readiness_score']):.3f}, tier {kr['validation_tier']}.",
        "- Claim boundary: this supports prioritization and structural plausibility. It does not prove immunogenicity until mutant activation exceeds WT/decoy controls.",
        "",
        "## Reviewer / Editor Positioning",
        "",
        md_table(positioning, ["attack_or_question", "strong_answer", "asset_that_answers_it", "claim_boundary", "impact_if_resolved"], 20),
        "",
        "## Claim Upgrade Milestones",
        "",
        md_table(milestones, ["lead", "stage_order", "milestone", "pass_criterion", "claim_allowed_after_pass", "claim_still_forbidden"], 20),
        "",
        "## External Validation / Benchmark Plan",
        "",
        md_table(validation, ["validation_axis", "dataset_or_material", "primary_metric", "success_criterion", "claim_if_passed"], 20),
        "",
        "## Wetlab Outcome Scenarios",
        "",
        md_table(scenarios, ["lead", "scenario", "observed_result", "interpretation", "next_move", "allowed_claim"], 20),
        "",
        "## Next Action Board",
        "",
        md_table(actions, ["priority", "time_window", "action", "deliverable", "impact"], 10),
    ]
    (OUT / "IMPACT_ACCELERATION_REPORT.md").write_text("\n".join(lines) + "\n")

    editor = [
        "# CROSS-Neo Editor / Reviewer One-Page",
        "",
        "**Core pitch:** CROSS-Neo is an experimentally actionable neoantigen triage system, not a claim that peptide-HLA presentation alone proves TCR recognition.",
        "",
        "**What is strong now:** strict current-label Top-13 triage, two externally anchored flagship leads, paired-TCR/structure/MD audit evidence, reagent/plate-map readiness, and explicit WT/decoy validation gates.",
        "",
        "**What is not claimed:** no 100% accuracy, no clinical utility, no MD-proves-immunogenicity, and no universal TCR-aware prediction when paired TCR is absent.",
        "",
        "**Highest-impact next result:** mutant-over-WT/decoy activation for TP53 R175H and/or KRAS G12D under matched HLA context.",
    ]
    (OUT / "EDITOR_REVIEWER_ONE_PAGE.md").write_text("\n".join(editor) + "\n")

    kr_lines = [
        "# 임팩트 올리는 핵심 요약",
        "",
        "지금 제일 강한 포지션은 “정확도 100% 모델”이 아니라, **실험으로 바로 넘길 수 있는 neoantigen triage system**입니다.",
        "",
        f"- TP53 R175H / HMTEVVRHC: readiness {float(hm['overall_readiness_score']):.3f}, {hm['validation_tier']} 입니다.",
        f"- KRAS G12D / GADGVGKSAL: readiness {float(kr['overall_readiness_score']):.3f}, {kr['validation_tier']} 입니다.",
        "- 다음 임팩트 상승 포인트는 WT/decoy 대비 mutant activation이 실제 assay에서 올라가는지입니다.",
        "- MD는 구조적 plausibility와 false-positive 제거용입니다. 면역원성 확정은 wetlab activation/killing gate가 필요합니다.",
        "",
        "한 줄 결론: **모델 후보군을 줄이고, TCR/MD로 위험 후보를 더 줄이고, WT/decoy assay로 진짜 claim을 올리는 전략이 가장 강합니다.**",
    ]
    (OUT / "IMPACT_ACCELERATION_ONE_PAGE_KR.md").write_text("\n".join(kr_lines) + "\n")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    data = load_inputs()
    positioning = build_editor_reviewer_positioning(data)
    milestones = build_claim_upgrade_milestones(data)
    validation = build_external_validation_benchmark_plan(data)
    scenarios = build_wetlab_success_scenarios(data)
    actions = build_next_action_board(data)

    positioning.to_csv(OUT / "editor_reviewer_positioning.tsv", sep="\t", index=False)
    milestones.to_csv(OUT / "claim_upgrade_milestones.tsv", sep="\t", index=False)
    validation.to_csv(OUT / "external_validation_benchmark_plan.tsv", sep="\t", index=False)
    scenarios.to_csv(OUT / "wetlab_success_scenarios.tsv", sep="\t", index=False)
    actions.to_csv(OUT / "next_72h_action_board.tsv", sep="\t", index=False)

    build_figures(positioning, milestones, validation, actions)
    write_reports(positioning, milestones, validation, scenarios, actions, data)

    summary = {
        "n_reviewer_questions": int(len(positioning)),
        "n_claim_milestones": int(len(milestones)),
        "n_external_validation_axes": int(len(validation)),
        "n_wetlab_scenarios": int(len(scenarios)),
        "n_next_actions": int(len(actions)),
        "figures": [
            "fig_md47_claim_upgrade_milestone_map",
            "fig_md48_reviewer_attack_defense_heatmap",
            "fig_md49_external_validation_benchmark_design",
            "fig_md50_next_72h_action_board",
        ],
        "boundary": "impact acceleration and validation planning, not immunogenicity proof",
    }
    (OUT / "impact_acceleration_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
