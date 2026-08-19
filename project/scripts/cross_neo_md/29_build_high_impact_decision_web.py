#!/usr/bin/env python3
"""Build a high-impact decision page for CROSS-Neo wetlab triage."""

from __future__ import annotations

import html
import json
import shutil
from pathlib import Path

import pandas as pd


REPO = Path(__file__).resolve().parents[3]
OUT = REPO / "project/results/cross_neo_md_audit_2026_05_10"
PKG = OUT / "high_impact_decision_package"
SIM = OUT / "immunogenicity_simulation_escalation"
LEAD = OUT / "two_lead_impact_package"
ASSAY = OUT / "assay_ready_translation_package"
PRECLIN = OUT / "preclinical_validation_protocol"
IMPACT = OUT / "impact_acceleration_package"
IMMUNO = OUT / "immunogenicity_prediction_upgrade"
FEEDBACK = OUT / "assay_feedback_learner"
SCENARIO = OUT / "assay_feedback_scenarios"
PHYSICS = OUT / "physics_gate_package"
FIG = OUT / "figures"
HUB = REPO / "project/papers_hub_2026_05_04"
WEB = Path("/var/www/papers/papers_hub_2026_05_04")
ASSET = HUB / "assets/cross_neo_md_audit"
WEB_ASSET = WEB / "assets/cross_neo_md_audit"
PAGE = HUB / "cross_neo_high_impact_decision.html"
WEB_PAGE = WEB / "cross_neo_high_impact_decision.html"


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t") if path.exists() else pd.DataFrame()


def esc(value: object) -> str:
    return html.escape("" if pd.isna(value) else str(value))


def fmt(value: object, digits: int = 3) -> str:
    try:
        if pd.isna(value):
            return "NA"
        return f"{float(value):.{digits}f}"
    except Exception:
        return esc(value)


def table(df: pd.DataFrame, cols: list[str], max_rows: int = 20) -> str:
    if df.empty:
        return "<p class='muted'>No rows.</p>"
    cols = [c for c in cols if c in df.columns]
    if not cols:
        return "<p class='muted'>Requested columns are absent.</p>"
    rows = ["<table><thead><tr>"]
    rows += [f"<th>{esc(c.replace('_', ' '))}</th>" for c in cols]
    rows.append("</tr></thead><tbody>")
    for _, record in df.head(max_rows).iterrows():
        rows.append("<tr>")
        for c in cols:
            v = record.get(c, "")
            rows.append(f"<td>{fmt(v) if isinstance(v, (int, float)) else esc(v)}</td>")
        rows.append("</tr>")
    rows.append("</tbody></table>")
    return "\n".join(rows)


def metric(df: pd.DataFrame, column: str, default: float = 0.0) -> float:
    if df.empty or column not in df.columns:
        return default
    return float(pd.to_numeric(df[column], errors="coerce").fillna(0).sum())


def safe_copy2(src: Path, dst: Path) -> None:
    try:
        shutil.copy2(src, dst)
    except PermissionError:
        print(f"[high-impact-web] skip permission-denied asset {dst}")


def copy_assets() -> None:
    ASSET.mkdir(parents=True, exist_ok=True)
    WEB_ASSET.mkdir(parents=True, exist_ok=True)
    for src in PKG.glob("*"):
        if src.is_file():
            safe_copy2(src, ASSET / src.name)
            safe_copy2(src, WEB_ASSET / src.name)
    if SIM.exists():
        for src in SIM.glob("*"):
            if src.is_file():
                safe_copy2(src, ASSET / src.name)
                safe_copy2(src, WEB_ASSET / src.name)
    if LEAD.exists():
        for src in LEAD.glob("*"):
            if src.is_file():
                safe_copy2(src, ASSET / src.name)
                safe_copy2(src, WEB_ASSET / src.name)
    if ASSAY.exists():
        for src in ASSAY.glob("*"):
            if src.is_file():
                safe_copy2(src, ASSET / src.name)
                safe_copy2(src, WEB_ASSET / src.name)
    if PRECLIN.exists():
        for src in PRECLIN.glob("*"):
            if src.is_file():
                safe_copy2(src, ASSET / src.name)
                safe_copy2(src, WEB_ASSET / src.name)
    if IMPACT.exists():
        for src in IMPACT.glob("*"):
            if src.is_file():
                safe_copy2(src, ASSET / src.name)
                safe_copy2(src, WEB_ASSET / src.name)
    if IMMUNO.exists():
        for src in IMMUNO.glob("*"):
            if src.is_file():
                safe_copy2(src, ASSET / src.name)
                safe_copy2(src, WEB_ASSET / src.name)
    if FEEDBACK.exists():
        for src in FEEDBACK.glob("*"):
            if src.is_file():
                safe_copy2(src, ASSET / src.name)
                safe_copy2(src, WEB_ASSET / src.name)
    if SCENARIO.exists():
        for src in SCENARIO.glob("*"):
            if src.is_file():
                safe_copy2(src, ASSET / src.name)
                safe_copy2(src, WEB_ASSET / src.name)
    if PHYSICS.exists():
        for src in PHYSICS.glob("*"):
            if src.is_file():
                safe_copy2(src, ASSET / src.name)
                safe_copy2(src, WEB_ASSET / src.name)
    for fig in [
        "fig_md27_no_fp_top13_candidate_evidence",
        "fig_md28_source_stratified_preset_performance",
        "fig_md29_claim_ladder",
        "fig_md30_wetlab_validation_workflow",
        "fig_md31_immunogenicity_simulation_ladder",
        "fig_md32_candidate_simulation_budget",
        "fig_md33_assay_to_simulation_bridge",
        "fig_md34_immunogenicity_simulation_readiness",
        "fig_md35_two_lead_evidence_moat",
        "fig_md36_specificity_lock_matrix",
        "fig_md37_impact_upgrade_flow",
        "fig_md38_reagent_order_sheet",
        "fig_md39_control_coverage",
        "fig_md40_assay_go_no_go_ladder",
        "fig_md41_reviewer_defense_moat",
        "fig_md46_evidence_to_claim_matrix",
        "fig_md42_preclinical_readiness_scorecard",
        "fig_md43_preclinical_validation_gates",
        "fig_md44_failure_mode_action_map",
        "fig_md45_preclinical_translation_flow",
        "fig_md47_claim_upgrade_milestone_map",
        "fig_md48_reviewer_attack_defense_heatmap",
        "fig_md49_external_validation_benchmark_design",
        "fig_md50_next_72h_action_board",
        "fig_md51_cross_neo_i_architecture",
        "fig_md52_label_feature_contract",
        "fig_md53_cross_neo_i_readiness_ranking",
        "fig_md54_immunogenicity_benchmark_ladder",
        "fig_md55_immunogenicity_risk_control_board",
        "fig_md56_assay_feedback_closed_loop",
        "fig_md57_prior_posterior_update",
        "fig_md58_active_learning_queue",
        "fig_md59_claim_state_distribution",
        "fig_md60_assay_scenario_posterior_heatmap",
        "fig_md61_assay_scenario_claim_transitions",
        "fig_md62_scenario_posterior_delta",
        "fig_md63_assay_scenario_simulator_flow",
        "fig_md64_physics_priority_rank",
        "fig_md65_physics_gate_ladder",
        "fig_md66_physics_budget",
        "fig_md67_physics_gate_rules",
        "fig_md22_threshold_precision_recall_frontier",
        "fig_md23_threshold_preset_confusion",
        "fig_md24_story_overview_flow",
        "fig_md26_representative_candidate_board",
    ]:
        for ext in (".png", ".pdf"):
            src = FIG / f"{fig}{ext}"
            if src.exists():
                safe_copy2(src, ASSET / src.name)
                safe_copy2(src, WEB_ASSET / src.name)


def build_guides(top13: pd.DataFrame, source_perf: pd.DataFrame, wetlab: pd.DataFrame) -> None:
    guide = pd.DataFrame(
        [
            {
                "id": "fig_md27",
                "title": "Top-13 no-FP candidate evidence",
                "main_message": "Strict current-label preset selects 13 positives and no false positives in the available labels.",
                "use": "Main result / wetlab triage opening figure",
                "data_used": "no_false_positive_top13_candidates.tsv",
                "claim_boundary": "Current-label enrichment only; needs external source-heldout validation.",
            },
            {
                "id": "fig_md28",
                "title": "Source-stratified preset performance",
                "main_message": "The strict preset has 0 FP within each current source, with different recall by source.",
                "use": "Reviewer robustness / source-shift figure",
                "data_used": "preset_source_stratified_performance.tsv",
                "claim_boundary": "Source-stratified current labels are not independent external validation.",
            },
            {
                "id": "fig_md29",
                "title": "Claim ladder",
                "main_message": "Separate allowed triage claims from wetlab-only and forbidden claims.",
                "use": "Discussion / reviewer-defense slide",
                "data_used": "claim_ladder.tsv",
                "claim_boundary": "Prevents overclaiming MD, structure, or TCR evidence.",
            },
            {
                "id": "fig_md30",
                "title": "Wetlab validation workflow",
                "main_message": "Candidates advance only through HLA binding, T-cell activation, and WT/decoy controls.",
                "use": "Methods / validation plan",
                "data_used": "wetlab_validation_plate_plan.tsv",
                "claim_boundary": "Wetlab assay context determines immunogenicity; model/MD alone cannot.",
            },
            {
                "id": "fig_md46",
                "title": "Evidence-to-claim matrix",
                "main_message": "Each impact claim is tied to actual data and a blocked overclaim.",
                "use": "Reviewer-defense / high-impact summary figure",
                "data_used": "evidence_to_claim_matrix.tsv",
                "claim_boundary": "Strong current-label prioritization is not external validation or immunogenicity proof.",
            },
        ]
    )
    guide_path = PKG / "high_impact_figure_table_caption_guide.tsv"
    guide.to_csv(guide_path, sep="\t", index=False)
    safe_copy2(guide_path, ASSET / guide_path.name)
    safe_copy2(guide_path, WEB_ASSET / guide_path.name)

    lines = [
        "# CROSS-Neo high-impact one-page summary",
        "",
        "## 핵심 결론",
        "",
        f"- 현재 label 기준 strict preset은 Top-{len(top13)} 후보를 남기고, 이 후보들은 모두 positive label입니다.",
        "- 이 결과는 wetlab 후보 선별용 decision-support evidence이며, 외부검증/임상효용/면역원성 증명은 아닙니다.",
        "- 가장 바로 실험할 후보는 GADGVGKSAL / HLA-C*08:02, HMTEVVRHC / HLA-A*02:01, KLILWRGLK / HLA-A*03:01, ILDKVLVHL / HLA-A*02:01입니다.",
        "",
        "## 왜 하이임팩트인가",
        "",
        "- cheap DL + Bayesian/dropout uncertainty로 먼저 줄입니다.",
        "- TCR/구조/MD는 비싼 최종 audit layer로만 씁니다.",
        "- 후보마다 required control과 go/no-go 기준이 붙어 wetlab plate로 바로 전환됩니다.",
        "",
        "## Source별 현재 성능",
        "",
    ]
    for _, row in source_perf[source_perf.get("preset_name", "").eq("NO_FALSE_POSITIVE_MAX_TP")].iterrows():
        lines.append(
            f"- {row['source_dataset']}: called={int(row['called_positive'])}, TP={int(row['TP'])}, FP={int(row['FP'])}, recall={float(row['recall']):.3f}"
        )
    lines += [
        "",
        "## 실험 우선순위",
        "",
    ]
    for _, row in wetlab.head(8).iterrows():
        lines.append(f"- {int(row['plate_order'])}. {row['peptide']} / {row['hla_4digit']}: {row['tier']} - {row['why_this_candidate']}")
    lines += [
        "",
        "## 금지 claim",
        "",
        "- MD가 immunogenicity를 증명한다고 말하지 않습니다.",
        "- 현재 label 최적화 결과를 외부검증이라고 말하지 않습니다.",
        "- TCR/structure가 없는 후보에서 TCR recognition을 증명했다고 말하지 않습니다.",
    ]
    summary_path = PKG / "HIGH_IMPACT_ONE_PAGE_KR.md"
    summary_path.write_text("\n".join(lines) + "\n")
    safe_copy2(summary_path, ASSET / summary_path.name)
    safe_copy2(summary_path, WEB_ASSET / summary_path.name)


def build_html(
    top13: pd.DataFrame,
    source_perf: pd.DataFrame,
    wetlab: pd.DataFrame,
    claims: pd.DataFrame,
    validation: pd.DataFrame,
    risks: pd.DataFrame,
    sim_matrix: pd.DataFrame,
    sim_bridge: pd.DataFrame,
    sim_rules: pd.DataFrame,
    sim_budget: pd.DataFrame,
    sim_summary: dict,
    lead_external: pd.DataFrame,
    lead_locks: pd.DataFrame,
    lead_wetlab: pd.DataFrame,
    evidence_matrix: pd.DataFrame,
    figure_plan: pd.DataFrame,
    lead_summary: dict,
    assay_reagents: pd.DataFrame,
    assay_plate: pd.DataFrame,
    assay_thresholds: pd.DataFrame,
    assay_moat: pd.DataFrame,
    assay_summary: dict,
    preclin_score: pd.DataFrame,
    preclin_gates: pd.DataFrame,
    preclin_failures: pd.DataFrame,
    preclin_materials: pd.DataFrame,
    preclin_summary: dict,
    impact_positioning: pd.DataFrame,
    impact_milestones: pd.DataFrame,
    impact_validation: pd.DataFrame,
    impact_scenarios: pd.DataFrame,
    impact_actions: pd.DataFrame,
    impact_summary: dict,
    immuno_layers: pd.DataFrame,
    immuno_labels: pd.DataFrame,
    immuno_readiness: pd.DataFrame,
    immuno_features: pd.DataFrame,
    immuno_benchmark: pd.DataFrame,
    immuno_summary: dict,
    physics_plan: pd.DataFrame,
    physics_budget: pd.DataFrame,
    physics_modules: pd.DataFrame,
    physics_rules: pd.DataFrame,
    physics_summary: dict,
    summary: dict,
) -> str:
    tp = int((top13.get("prediction_outcome", pd.Series(dtype=str)).astype(str).eq("TP")).sum())
    fp = int((top13.get("prediction_outcome", pd.Series(dtype=str)).astype(str).eq("FP")).sum())
    tier1 = int(wetlab.get("tier", pd.Series(dtype=str)).astype(str).str.contains("TIER_1", na=False).sum())
    paired = int((pd.to_numeric(top13.get("paired_tcr_evidence_count", pd.Series(dtype=float)), errors="coerce").fillna(0) > 0).sum())
    md_supported = int(top13.get("md_label", pd.Series(dtype=str)).astype(str).str.contains("MODERATE|STRONG|VERY_STRONG", regex=True, na=False).sum())
    sim_modules = int(sim_summary.get("planned_modules", len(sim_matrix)))
    sim_ns = float(sim_summary.get("total_requested_short_md_ns", 0.0))
    lead_count = int(len(lead_external))
    assay_rows = int(assay_summary.get("n_plate_rows", len(assay_plate)))
    preclin_gates_n = int(preclin_summary.get("n_gates", len(preclin_gates)))
    impact_axes = int(impact_summary.get("n_external_validation_axes", len(impact_validation)))
    immuno_layers_n = int(immuno_summary.get("n_model_layers", len(immuno_layers)))
    physics_modules_n = int(physics_summary.get("n_physics_modules", len(physics_modules)))
    physics_flagships_n = int(physics_summary.get("n_flagships", 0))
    fig_cards = "\n".join(
        f"""
        <article>
          <img src="assets/cross_neo_md_audit/{name}.png" alt="{esc(title)}">
          <h3>{esc(title)}</h3>
          <p>{esc(desc)}</p>
          <a href="assets/cross_neo_md_audit/{name}.png">PNG</a>
          <a href="assets/cross_neo_md_audit/{name}.pdf">PDF</a>
        </article>
        """
        for name, title, desc in [
            ("fig_md27_no_fp_top13_candidate_evidence", "Top-13 No-FP Evidence", "Evidence board for the strict current-label Top-13 candidate list."),
            ("fig_md28_source_stratified_preset_performance", "Source-Stratified Performance", "Per-source TP/FP/FN behavior for optimized decision presets."),
            ("fig_md29_claim_ladder", "Claim Ladder", "Allowed, future, and forbidden claim levels for model/TCR/MD evidence."),
            ("fig_md30_wetlab_validation_workflow", "Wetlab Workflow", "Go/no-go validation sequence from peptide-HLA binding to WT/decoy-controlled T-cell assays."),
            ("fig_md31_immunogenicity_simulation_ladder", "Immunogenicity Simulation Ladder", "Which cheap, short-MD, specificity-control, and expensive modules are added before wetlab testing."),
            ("fig_md32_candidate_simulation_budget", "Simulation Budget", "Requested short-MD trajectory length by candidate before wetlab escalation."),
            ("fig_md33_assay_to_simulation_bridge", "Assay Bridge", "How HLA binding, multimer, ELISpot/ICS, and killing assays map to simulation evidence."),
            ("fig_md34_immunogenicity_simulation_readiness", "Readiness Status", "Which simulation modules are ready, queued, held, or blocked by WT/decoy/TCR curation."),
            ("fig_md35_two_lead_evidence_moat", "Two-Lead Evidence Moat", "External literature, paired TCR, MD, WT controls, and reagent actionability for the two flagship candidates."),
            ("fig_md36_specificity_lock_matrix", "Specificity Locks", "Mutation, HLA, structure, MD, WT/decoy, and assay status before any immunogenicity claim."),
            ("fig_md37_impact_upgrade_flow", "Impact Upgrade Flow", "How candidate ranking becomes a controlled immunogenicity test plan."),
            ("fig_md38_reagent_order_sheet", "Reagent Order Sheet", "Mutant, WT, decoy, and pHLA reagent rows needed to start controlled assays."),
            ("fig_md39_control_coverage", "Control Coverage", "Which controls are explicitly represented for each flagship lead."),
            ("fig_md40_assay_go_no_go_ladder", "Assay Go/No-Go Ladder", "Presentation, multimer, activation, and killing gates before any strong claim."),
            ("fig_md41_reviewer_defense_moat", "Reviewer Defense Moat", "Likely reviewer objections and the asset that answers each one."),
            ("fig_md42_preclinical_readiness_scorecard", "Preclinical Readiness", "Readiness scorecard for moving ranked candidates into controlled validation."),
            ("fig_md43_preclinical_validation_gates", "Validation Gates", "Identity, presentation, recognition, activation, and killing gates before stronger claims."),
            ("fig_md44_failure_mode_action_map", "Failure Actions", "How to react if presentation, recognition, activation, or killing gates fail."),
            ("fig_md45_preclinical_translation_flow", "Preclinical Translation Flow", "End-to-end flow from CROSS-Neo ranking to validation-ready experiments."),
            ("fig_md46_evidence_to_claim_matrix", "Evidence-To-Claim Matrix", "Actual evidence, supported claim, and blocked overclaim in one reviewer-facing board."),
            ("fig_md47_claim_upgrade_milestone_map", "Claim Upgrade Milestones", "What must pass before each stronger claim becomes defensible."),
            ("fig_md48_reviewer_attack_defense_heatmap", "Reviewer Defense Coverage", "Where model, source-shift, wetlab, and claim-boundary objections are answered."),
            ("fig_md49_external_validation_benchmark_design", "External Validation Design", "Benchmark axes that would turn the story from retrospective audit into stronger validation."),
            ("fig_md50_next_72h_action_board", "Next 72h Action Board", "Fastest actions for turning the current package into higher-impact evidence."),
            ("fig_md51_cross_neo_i_architecture", "CROSS-Neo-I Architecture", "Immunogenicity prediction stack from presentation to calibration."),
            ("fig_md52_label_feature_contract", "Label/Feature Contract", "Which labels train which model layers and what cannot be used as a positive."),
            ("fig_md53_cross_neo_i_readiness_ranking", "CROSS-Neo-I Readiness", "Top candidates ranked by readiness for immunogenicity-label generation."),
            ("fig_md54_immunogenicity_benchmark_ladder", "Immunogenicity Benchmark Ladder", "Heldout and prospective benchmarks required before stronger prediction claims."),
            ("fig_md55_immunogenicity_risk_control_board", "Risk-Control Board", "Leakage, source-shift, and claim-risk controls for the immunogenicity model."),
            ("fig_md56_assay_feedback_closed_loop", "Assay Feedback Loop", "Closed-loop route from CROSS-Neo prior to assay posterior and next experiment."),
            ("fig_md57_prior_posterior_update", "Prior-To-Posterior Update", "How real assay labels will update candidate readiness."),
            ("fig_md58_active_learning_queue", "Active-Learning Queue", "Next experiment utility ranking for the candidate set."),
            ("fig_md59_claim_state_distribution", "Claim-State Distribution", "Allowed claim states before and after assay feedback."),
            ("fig_md60_assay_scenario_posterior_heatmap", "Assay Scenario Heatmap", "What-if posterior movement for flagship leads under different assay outcomes."),
            ("fig_md61_assay_scenario_claim_transitions", "Scenario Claim Transitions", "How claim states change under simulated assay outcomes."),
            ("fig_md62_scenario_posterior_delta", "Scenario Posterior Delta", "Which assay outcomes raise or lower readiness most."),
            ("fig_md63_assay_scenario_simulator_flow", "Scenario Simulator Flow", "Planning-only route from what-if assay outcome to next decision."),
            ("fig_md64_physics_priority_rank", "Physics Priority Rank", "Which candidates justify endpoint energy, PMF, or FEP escalation."),
            ("fig_md65_physics_gate_ladder", "Physics Gate Ladder", "Escalation ladder from cheap geometry checks to high-cost physics."),
            ("fig_md66_physics_budget", "Physics Budget", "Estimated GPU-hour budget by module and candidate."),
            ("fig_md67_physics_gate_rules", "Physics Gate Rules", "When to launch, hold, or stop each physics module."),
        ]
    )
    best = top13.head(4)
    best_cards = "\n".join(
        f"""
        <article class="candidate">
          <h3>{esc(row.get('peptide'))} <span>{esc(row.get('hla_4digit'))}</span></h3>
          <p>{esc(row.get('source_dataset'))} · label {esc(row.get('actual_label'))} · {esc(row.get('recommendation_tier'))}</p>
          <dl>
            <dt>Main DL</dt><dd>{fmt(row.get('main_dl_score'))}</dd>
            <dt>TCR branch</dt><dd>{fmt(row.get('tcr_augmented_score_mean'))}</dd>
            <dt>Paired TCR evidence</dt><dd>{fmt(row.get('paired_tcr_evidence_count'), 0)}</dd>
            <dt>Structure/MD</dt><dd>{fmt(row.get('baker_structural_score'))} / {esc(row.get('md_label'))}</dd>
          </dl>
        </article>
        """
        for _, row in best.iterrows()
    )
    return f"""<!doctype html>
<html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>CROSS-Neo High-Impact Decision Package</title>
<style>
body{{margin:0;background:#07111f;color:#edf5ff;font-family:Inter,Arial,sans-serif;line-height:1.55}}
a{{color:#39d4b5;text-decoration:none}} a:hover{{text-decoration:underline}}
header{{padding:42px 32px;background:#0f2035;border-bottom:1px solid #26364d}}
.wrap{{max-width:1360px;margin:0 auto;padding:24px 30px 72px}}
h1{{font-family:Georgia,serif;font-size:46px;margin:0 0 10px;color:#fff8e8}}
h2{{font-family:Georgia,serif;color:#fff2d0;font-size:28px;margin:0 0 14px}}
h3{{margin:4px 0 8px;color:#fff8e8}} h3 span{{color:#f2c46d;font-size:15px}}
.lead{{max-width:1050px;color:#cbd7e7;font-size:17px}}
.stats{{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:12px;margin:20px 0 6px}}
.stat{{background:#101d30;border:1px solid #293b55;border-radius:8px;padding:13px}}
.stat b{{display:block;color:#f2c46d;font-size:26px}}
section{{border-bottom:1px solid #26364d;padding:24px 0}}
table{{border-collapse:collapse;width:100%;font-size:12px}} th,td{{border:1px solid #293b55;padding:7px;vertical-align:top}} th{{background:#13243a;color:#f2c46d}}
.links{{display:flex;flex-wrap:wrap;gap:10px}} .pill{{border:1px solid #395170;border-radius:999px;padding:7px 10px;background:#101d30}}
.grid2{{display:grid;grid-template-columns:repeat(2,1fr);gap:16px}} .grid4{{display:grid;grid-template-columns:repeat(4,1fr);gap:14px}}
article{{background:#101d30;border:1px solid #293b55;border-radius:8px;padding:13px}}
img{{width:100%;background:white;border-radius:6px}}
.warn{{border-left:4px solid #ff7b72;background:#111d2e;padding:13px 15px;margin:14px 0}}
.ok{{border-left:4px solid #39d4b5;background:#111d2e;padding:13px 15px;margin:14px 0}}
.muted{{color:#9fb0c7}} dl{{display:grid;grid-template-columns:1fr 1fr;gap:4px 10px;margin:10px 0 0}} dt{{color:#9fb0c7}} dd{{margin:0;color:#edf5ff;font-weight:700}}
@media(max-width:1000px){{.stats,.grid2,.grid4{{grid-template-columns:1fr}} h1{{font-size:34px}}}}
</style></head><body>
<header><div class="wrap">
<h1>CROSS-Neo High-Impact Decision Package</h1>
<p class="lead">A compact decision layer for choosing wetlab neoantigen candidates: cheap DL first, uncertainty gates second, TCR/structure/MD only for high-value escalation.</p>
<div class="stats">
<div class="stat"><b>{len(top13)}</b><span>strict Top candidates</span></div>
<div class="stat"><b>{tp}/{fp}</b><span>TP / FP in current labels</span></div>
<div class="stat"><b>{tier1}</b><span>Tier-1 immediate assays</span></div>
<div class="stat"><b>{paired}</b><span>Top rows with paired TCR evidence</span></div>
<div class="stat"><b>{md_supported}</b><span>MD-supported Top rows</span></div>
<div class="stat"><b>{assay_rows}</b><span>assay plate-map rows</span></div>
<div class="stat"><b>{preclin_gates_n}</b><span>preclinical validation gates</span></div>
<div class="stat"><b>{impact_axes}</b><span>external validation axes</span></div>
<div class="stat"><b>{immuno_layers_n}</b><span>CROSS-Neo-I model layers</span></div>
<div class="stat"><b>{physics_modules_n}/{physics_flagships_n}</b><span>physics modules / flagships</span></div>
</div>
</div></header>
<main class="wrap">
<section><h2>Decision</h2>
<div class="ok">Use the strict Top-13 list as the immediate wetlab triage board, with GADGVGKSAL and HMTEVVRHC as the highest-value structure/TCR case studies. Keep the claim as prioritization until WT/decoy assays and external validation are done.</div>
<div class="warn">This is not 100% accuracy, clinical utility, or immunogenicity proof. The zero-FP result is measured on the currently available labels and must be stress-tested externally.</div>
<div class="links">
<a class="pill" href="cross_neo_dl_first_slider_dashboard.html">Interactive threshold dashboard</a>
<a class="pill" href="cross_neo_decision_storyboard.html">Decision storyboard</a>
<a class="pill" href="cross_neo_md_audit_dossier.html">Full MD/TCR dossier</a>
<a class="pill" href="assets/cross_neo_md_audit/HIGH_IMPACT_DECISION_REPORT.md">Full report</a>
<a class="pill" href="assets/cross_neo_md_audit/HIGH_IMPACT_ONE_PAGE_KR.md">One-page Korean summary</a>
<a class="pill" href="assets/cross_neo_md_audit/KAKAO_ONE_SHOT_HIGH_IMPACT_KR.md">Kakao one-shot KR</a>
<a class="pill" href="assets/cross_neo_md_audit/high_impact_figure_table_caption_guide.tsv">Figure/table guide</a>
<a class="pill" href="assets/cross_neo_md_audit/evidence_to_claim_matrix.tsv">Evidence-to-claim matrix</a>
<a class="pill" href="assets/cross_neo_md_audit/manuscript_figure_table_plan.tsv">Manuscript figure/table plan</a>
<a class="pill" href="assets/cross_neo_md_audit/IMMUNOGENICITY_SIMULATION_ESCALATION_REPORT.md">Immunogenicity simulation report</a>
<a class="pill" href="assets/cross_neo_md_audit/immunogenicity_simulation_run_commands.md">Simulation commands</a>
<a class="pill" href="assets/cross_neo_md_audit/RUNPOD_IMMUNOGENICITY_SIM_STATUS.md">RunPod sim status</a>
<a class="pill" href="assets/cross_neo_md_audit/TWO_LEAD_IMPACT_UPGRADE_REPORT.md">Two-lead impact report</a>
<a class="pill" href="assets/cross_neo_md_audit/TWO_LEAD_MANUSCRIPT_INSERT.md">Two-lead manuscript insert</a>
<a class="pill" href="assets/cross_neo_md_audit/TWO_LEAD_ONE_PAGE_KR.md">Two-lead Korean brief</a>
<a class="pill" href="assets/cross_neo_md_audit/ASSAY_READY_TRANSLATION_REPORT.md">Assay-ready report</a>
<a class="pill" href="assets/cross_neo_md_audit/ASSAY_READY_ONE_PAGE_KR.md">Assay-ready Korean brief</a>
<a class="pill" href="assets/cross_neo_md_audit/PRECLINICAL_VALIDATION_PROTOCOL.md">Preclinical validation protocol</a>
<a class="pill" href="assets/cross_neo_md_audit/COLLABORATION_BRIEF.md">Collaboration brief</a>
<a class="pill" href="assets/cross_neo_md_audit/PRECLINICAL_VALIDATION_ONE_PAGE_KR.md">Preclinical Korean brief</a>
<a class="pill" href="assets/cross_neo_md_audit/IMPACT_ACCELERATION_REPORT.md">Impact acceleration report</a>
<a class="pill" href="assets/cross_neo_md_audit/EDITOR_REVIEWER_ONE_PAGE.md">Editor/reviewer one-page</a>
<a class="pill" href="assets/cross_neo_md_audit/IMPACT_ACCELERATION_ONE_PAGE_KR.md">Impact Korean brief</a>
<a class="pill" href="assets/cross_neo_md_audit/CROSS_NEO_I_IMMUNOGENICITY_UPGRADE_REPORT.md">CROSS-Neo-I immunogenicity report</a>
<a class="pill" href="assets/cross_neo_md_audit/CROSS_NEO_I_ONE_PAGE_KR.md">CROSS-Neo-I Korean brief</a>
<a class="pill" href="assets/cross_neo_md_audit/KAKAO_CROSS_NEO_I_IMMUNOGENICITY_KR.md">CROSS-Neo-I Kakao brief</a>
<a class="pill" href="cross_neo_assay_feedback_loop.html">Assay feedback loop</a>
<a class="pill" href="cross_neo_assay_feedback_scenarios.html">Assay scenario simulator</a>
<a class="pill" href="cross_neo_physics_gate.html">Physics gate</a>
<a class="pill" href="cross_neo_physics_launch_sheet.html">Physics launch sheet</a>
<a class="pill" href="cross_neo_endpoint_energy_gate.html">Endpoint energy gate</a>
<a class="pill" href="cross_neo_physics_case_study_board.html">Physics case study board</a>
<a class="pill" href="cross_neo_integrated_decision_matrix.html">Integrated decision matrix</a>
<a class="pill" href="cross_neo_decision_atlas.html">Decision atlas</a>
<a class="pill" href="cross_neo_executive_impact_console.html">Executive impact console</a>
<a class="pill" href="cross_neo_flagship_execution_packet.html">Flagship execution packet</a>
<a class="pill" href="cross_neo_flagship_command_center.html">Flagship command center</a>
<a class="pill" href="cross_neo_flagship_decision_tower.html">Flagship decision tower</a>
<a class="pill" href="assets/cross_neo_md_audit/ASSAY_FEEDBACK_LEARNER_REPORT.md">Assay feedback report</a>
<a class="pill" href="assets/cross_neo_md_audit/ASSAY_FEEDBACK_SCENARIO_SIMULATOR_REPORT.md">Assay scenario report</a>
<a class="pill" href="assets/cross_neo_md_audit/PHYSICS_GATE_REPORT.md">Physics gate report</a>
<a class="pill" href="assets/cross_neo_md_audit/PHYSICS_GATE_ONE_PAGE_KR.md">Physics gate Korean brief</a>
<a class="pill" href="assets/cross_neo_md_audit/assay_feedback_template.tsv">Assay feedback TSV template</a>
</div></section>
<section><h2>Top Case Cards</h2><div class="grid4">{best_cards}</div></section>
<section><h2>Representative Figures</h2><div class="grid2">{fig_cards}</div></section>
<section><h2>Two-Lead Impact Upgrade</h2>
<p class="muted">The strongest upgrade is no longer only a ranked list. It is a controlled two-lead case-study plan: KRAS G12D / HLA-C*08:02 and TP53 R175H / HLA-A*02:01, each with WT/decoy specificity requirements before immunogenicity claims.</p>
{table(lead_external, ['lead','mutation','mutant_peptide','wildtype_peptide','anchor_preserved_decoy','hla','external_support_level','crossneo_md_label','crossneo_md_score','paired_tcr_evidence_count','claim_upgrade'], 4)}
{table(lead_locks, ['lead','specificity_lock','status','detail'], 20)}
{table(lead_wetlab, ['order','lead','assay','matched_controls','why','expected_output','stop_or_hold_rule'], 12)}
</section>
<section><h2>Assay-Ready Translation Package</h2>
<p class="muted">This converts the two-lead case study into a practical wetlab package: reagent order rows, plate-map rows, go/no-go gates, and reviewer-defense mapping.</p>
{table(assay_reagents, ['lead','item_type','item_name','sequence','hla_context','minimum_spec','purpose','priority','ordering_note'], 20)}
{table(assay_thresholds, ['lead','gate','go_threshold','hold_or_no_go','claim_allowed_if_go'], 10)}
{table(assay_plate, ['plate','well','lead','assay','antigen_type','sequence','hla','concentration_level','replicate','readout'], 24)}
{table(assay_moat, ['reviewer_attack','defense_asset','where_addressed'], 8)}
</section>
<section><h2>Preclinical Validation Protocol</h2>
<p class="muted">This layer turns the top two candidates into a go/no-go validation program. It keeps the allowed claim narrow until WT, decoy, multimer/activation, and killing controls pass.</p>
{table(preclin_score, ['lead','mutant_peptide','hla','wildtype_peptide','decoy_peptide','overall_readiness_score','validation_tier','recommended_next_action','blocked_claim'], 8)}
{table(preclin_gates, ['lead','gate_order','gate','go_criterion','hold_or_no_go_action','claim_allowed_if_pass'], 12)}
{table(preclin_materials, ['lead','material','requirement','why_needed','priority'], 12)}
{table(preclin_failures, ['lead','failure_mode','interpretation','next_analysis'], 14)}
</section>
<section><h2>Impact Acceleration / Reviewer Defense</h2>
<p class="muted">This section translates the current evidence into a stronger manuscript strategy: what can be claimed now, what must pass next, which reviewer attacks are already covered, and what raises the story most quickly.</p>
{table(impact_positioning, ['attack_or_question','strong_answer','asset_that_answers_it','current_evidence','claim_boundary','impact_if_resolved'], 8)}
{table(impact_milestones, ['lead','stage_order','milestone','pass_criterion','claim_allowed_after_pass','claim_still_forbidden','impact_level_if_passed'], 14)}
{table(impact_validation, ['validation_axis','dataset_or_material','primary_metric','success_criterion','claim_if_passed'], 8)}
{table(impact_actions, ['priority','time_window','action','deliverable','impact'], 8)}
{table(impact_scenarios, ['lead','scenario','observed_result','interpretation','next_move','allowed_claim'], 12)}
</section>
<section><h2>CROSS-Neo-I Immunogenicity Prediction Upgrade</h2>
<p class="muted">This layer raises the impact from candidate prioritization to an explicit immunogenicity-prediction model contract. The endpoint is mutant-specific T-cell activation above WT/decoy controls; current values are readiness and benchmark design, not final immunogenicity probabilities.</p>
{table(immuno_layers, ['layer','model_or_feature_family','needed_inputs','target_signal','what_it_adds','claim_boundary'], 12)}
{table(immuno_labels, ['label_name','positive_definition','negative_definition','not_allowed_as_positive','used_for'], 8)}
{table(immuno_readiness, ['peptide','hla_4digit','presentation_proxy','paired_tcr_proxy','md_proxy','wt_decoy_specificity_ready','cross_neo_i_readiness_score','cross_neo_i_tier','next_data_needed'], 13)}
{table(immuno_benchmark, ['benchmark','split','metric','success_criterion','claim_if_passed'], 8)}
</section>
<section><h2>Evidence To Claim Matrix</h2>
{table(evidence_matrix, ['evidence_layer','actual_data','supports_claim','allowed_strength_now','blocked_overclaim','figure_or_table'], 12)}
</section>
<section><h2>Manuscript Figure / Table Plan</h2>
{table(figure_plan, ['slot','asset','purpose','data_used','claim_boundary'], 12)}
</section>
<section><h2>Immunogenicity Simulation Escalation</h2>
<p class="muted">Simulation is added as a pre-wetlab culling and interpretation layer. The current plan contains {sim_modules} modules, including {sim_ns:.1f} ns of short explicit-MD work across the strict Top-13.</p>
{table(sim_matrix, ['priority_order','peptide','hla_4digit','assay_readout_goal','simulation_module','escalation_tier','complex_scope','control_type','replicates','ns_per_rep','total_requested_ns','launch_status','blocking_dependency','go_no_go_signal'], 28)}
</section>
<section><h2>Physics Gate</h2>
<p class="muted">This is the expensive skepticism layer. It ranks which candidates justify endpoint energies, PMF, or FEP, and which should be held back. Physics helps reject or rank; it does not prove immunogenicity.</p>
{table(physics_modules, ['module','physics_question','when_to_use','input_requirement','output','claim_boundary','cost_class'], 8)}
{table(physics_plan, ['row_id','peptide','hla_4digit','physics_priority_score','physics_tier','recommended_physics_module','recommended_next_step','cost_class','expected_runtime_class'], 13)}
{table(physics_budget, ['peptide','hla_4digit','module','estimated_gpu_hours_low','estimated_gpu_hours_high','decision_use'], 12)}
{table(physics_rules, ['observed_pattern','physics_gate_response'], 12)}
</section>
<section><h2>Simulation Budget And Go/No-Go</h2>
{table(sim_budget, ['peptide','hla_4digit','planned_jobs','requested_ns','gpu_hours_low','gpu_hours_high','blocked_jobs'], 13)}
{table(sim_rules, ['stage','go_rule','no_go_or_hold'], 8)}
</section>
<section><h2>Wetlab Assay To Simulation Bridge</h2>
{table(sim_bridge, ['wetlab_assay','simulation_support','decision_use','not_supported_claim'], 8)}
</section>
<section><h2>Top-13 No-FP Candidates</h2>
{table(top13, ['priority_order','row_id','peptide','hla_4digit','source_dataset','actual_label','prediction_outcome','main_dl_score','bayes_mean','tcr_augmented_score_mean','paired_tcr_evidence_count','baker_structural_score','md_label','recommendation_tier','why'], 13)}
</section>
<section><h2>Wetlab Plate Plan</h2>
{table(wetlab, ['plate_order','tier','row_id','peptide','hla_4digit','source_dataset','recommended_assays','required_controls','sample_requirement','go_no_go','why_this_candidate'], 13)}
</section>
<section><h2>Source-Stratified Performance</h2>
{table(source_perf, ['preset_name','source_dataset','n','positive','negative','called_positive','TP','TN','FP','FN','precision','recall','FPR'], 18)}
</section>
<section><h2>Claim Ladder</h2>
{table(claims, ['claim_level','claim','supporting_artifacts','required_before_stronger_claim','risk'], 8)}
</section>
<section><h2>External Validation Plan</h2>
{table(validation, ['validation_axis','dataset_needed','metric','success_criterion','owner'], 8)}
</section>
<section><h2>Reviewer Risk Register</h2>
{table(risks, ['risk','why_reviewer_will_ask','defense','next_action'], 8)}
</section>
<section><h2>Data Paths</h2>
<p class="muted">{esc(PKG)}</p>
<div class="links">
<a class="pill" href="assets/cross_neo_md_audit/no_false_positive_top13_candidates.tsv">no_false_positive_top13_candidates.tsv</a>
<a class="pill" href="assets/cross_neo_md_audit/wetlab_validation_plate_plan.tsv">wetlab_validation_plate_plan.tsv</a>
<a class="pill" href="assets/cross_neo_md_audit/preset_source_stratified_performance.tsv">preset_source_stratified_performance.tsv</a>
<a class="pill" href="assets/cross_neo_md_audit/claim_ladder.tsv">claim_ladder.tsv</a>
<a class="pill" href="assets/cross_neo_md_audit/external_validation_plan.tsv">external_validation_plan.tsv</a>
<a class="pill" href="assets/cross_neo_md_audit/reviewer_risk_register.tsv">reviewer_risk_register.tsv</a>
<a class="pill" href="assets/cross_neo_md_audit/evidence_to_claim_matrix.tsv">evidence_to_claim_matrix.tsv</a>
<a class="pill" href="assets/cross_neo_md_audit/manuscript_figure_table_plan.tsv">manuscript_figure_table_plan.tsv</a>
<a class="pill" href="assets/cross_neo_md_audit/immunogenicity_simulation_matrix.tsv">immunogenicity_simulation_matrix.tsv</a>
<a class="pill" href="assets/cross_neo_md_audit/candidate_simulation_budget.tsv">candidate_simulation_budget.tsv</a>
<a class="pill" href="assets/cross_neo_md_audit/wetlab_assay_simulation_bridge.tsv">wetlab_assay_simulation_bridge.tsv</a>
<a class="pill" href="assets/cross_neo_md_audit/two_lead_external_evidence.tsv">two_lead_external_evidence.tsv</a>
<a class="pill" href="assets/cross_neo_md_audit/two_lead_specificity_locks.tsv">two_lead_specificity_locks.tsv</a>
<a class="pill" href="assets/cross_neo_md_audit/two_lead_wetlab_order.tsv">two_lead_wetlab_order.tsv</a>
<a class="pill" href="assets/cross_neo_md_audit/assay_reagent_order_sheet.tsv">assay_reagent_order_sheet.tsv</a>
<a class="pill" href="assets/cross_neo_md_audit/assay_plate_map.tsv">assay_plate_map.tsv</a>
<a class="pill" href="assets/cross_neo_md_audit/assay_go_no_go_thresholds.tsv">assay_go_no_go_thresholds.tsv</a>
<a class="pill" href="assets/cross_neo_md_audit/preclinical_readiness_scorecard.tsv">preclinical_readiness_scorecard.tsv</a>
<a class="pill" href="assets/cross_neo_md_audit/preclinical_validation_gates.tsv">preclinical_validation_gates.tsv</a>
<a class="pill" href="assets/cross_neo_md_audit/sample_material_request.tsv">sample_material_request.tsv</a>
<a class="pill" href="assets/cross_neo_md_audit/preclinical_failure_mode_actions.tsv">preclinical_failure_mode_actions.tsv</a>
<a class="pill" href="assets/cross_neo_md_audit/editor_reviewer_positioning.tsv">editor_reviewer_positioning.tsv</a>
<a class="pill" href="assets/cross_neo_md_audit/claim_upgrade_milestones.tsv">claim_upgrade_milestones.tsv</a>
<a class="pill" href="assets/cross_neo_md_audit/external_validation_benchmark_plan.tsv">external_validation_benchmark_plan.tsv</a>
<a class="pill" href="assets/cross_neo_md_audit/wetlab_success_scenarios.tsv">wetlab_success_scenarios.tsv</a>
<a class="pill" href="assets/cross_neo_md_audit/next_72h_action_board.tsv">next_72h_action_board.tsv</a>
<a class="pill" href="assets/cross_neo_md_audit/immunogenicity_model_layers.tsv">immunogenicity_model_layers.tsv</a>
<a class="pill" href="assets/cross_neo_md_audit/immunogenicity_label_contract.tsv">immunogenicity_label_contract.tsv</a>
<a class="pill" href="assets/cross_neo_md_audit/cross_neo_i_candidate_readiness.tsv">cross_neo_i_candidate_readiness.tsv</a>
<a class="pill" href="assets/cross_neo_md_audit/immunogenicity_benchmark_contract.tsv">immunogenicity_benchmark_contract.tsv</a>
<a class="pill" href="assets/cross_neo_md_audit/physics_modules.tsv">physics_modules.tsv</a>
<a class="pill" href="assets/cross_neo_md_audit/physics_gate_candidate_plan.tsv">physics_gate_candidate_plan.tsv</a>
<a class="pill" href="assets/cross_neo_md_audit/physics_gate_budget.tsv">physics_gate_budget.tsv</a>
<a class="pill" href="assets/cross_neo_md_audit/physics_gate_rules.tsv">physics_gate_rules.tsv</a>
<a class="pill" href="assets/cross_neo_md_audit/assay_feedback_schema.tsv">assay_feedback_schema.tsv</a>
<a class="pill" href="assets/cross_neo_md_audit/assay_feedback_template.tsv">assay_feedback_template.tsv</a>
<a class="pill" href="assets/cross_neo_md_audit/candidate_posterior_updates.tsv">candidate_posterior_updates.tsv</a>
<a class="pill" href="assets/cross_neo_md_audit/active_learning_next_experiments.tsv">active_learning_next_experiments.tsv</a>
<a class="pill" href="assets/cross_neo_md_audit/claim_state_updates.tsv">claim_state_updates.tsv</a>
</div></section>
</main></body></html>"""


def main() -> None:
    copy_assets()
    top13 = read_tsv(PKG / "no_false_positive_top13_candidates.tsv")
    source_perf = read_tsv(PKG / "preset_source_stratified_performance.tsv")
    wetlab = read_tsv(PKG / "wetlab_validation_plate_plan.tsv")
    claims = read_tsv(PKG / "claim_ladder.tsv")
    validation = read_tsv(PKG / "external_validation_plan.tsv")
    risks = read_tsv(PKG / "reviewer_risk_register.tsv")
    sim_matrix = read_tsv(SIM / "immunogenicity_simulation_matrix.tsv")
    sim_bridge = read_tsv(SIM / "wetlab_assay_simulation_bridge.tsv")
    sim_rules = read_tsv(SIM / "simulation_go_no_go_rules.tsv")
    sim_budget = read_tsv(SIM / "candidate_simulation_budget.tsv")
    sim_summary_path = SIM / "immunogenicity_simulation_summary.json"
    sim_summary = json.loads(sim_summary_path.read_text()) if sim_summary_path.exists() else {}
    lead_external = read_tsv(LEAD / "two_lead_external_evidence.tsv")
    lead_locks = read_tsv(LEAD / "two_lead_specificity_locks.tsv")
    lead_wetlab = read_tsv(LEAD / "two_lead_wetlab_order.tsv")
    evidence_matrix = read_tsv(PKG / "evidence_to_claim_matrix.tsv")
    figure_plan = read_tsv(PKG / "manuscript_figure_table_plan.tsv")
    lead_summary_path = LEAD / "two_lead_impact_summary.json"
    lead_summary = json.loads(lead_summary_path.read_text()) if lead_summary_path.exists() else {}
    assay_reagents = read_tsv(ASSAY / "assay_reagent_order_sheet.tsv")
    assay_plate = read_tsv(ASSAY / "assay_plate_map.tsv")
    assay_thresholds = read_tsv(ASSAY / "assay_go_no_go_thresholds.tsv")
    assay_moat = read_tsv(ASSAY / "assay_reviewer_defense_moat.tsv")
    assay_summary_path = ASSAY / "assay_ready_summary.json"
    assay_summary = json.loads(assay_summary_path.read_text()) if assay_summary_path.exists() else {}
    preclin_score = read_tsv(PRECLIN / "preclinical_readiness_scorecard.tsv")
    preclin_gates = read_tsv(PRECLIN / "preclinical_validation_gates.tsv")
    preclin_failures = read_tsv(PRECLIN / "preclinical_failure_mode_actions.tsv")
    preclin_materials = read_tsv(PRECLIN / "sample_material_request.tsv")
    preclin_summary_path = PRECLIN / "preclinical_validation_summary.json"
    preclin_summary = json.loads(preclin_summary_path.read_text()) if preclin_summary_path.exists() else {}
    impact_positioning = read_tsv(IMPACT / "editor_reviewer_positioning.tsv")
    impact_milestones = read_tsv(IMPACT / "claim_upgrade_milestones.tsv")
    impact_validation = read_tsv(IMPACT / "external_validation_benchmark_plan.tsv")
    impact_scenarios = read_tsv(IMPACT / "wetlab_success_scenarios.tsv")
    impact_actions = read_tsv(IMPACT / "next_72h_action_board.tsv")
    impact_summary_path = IMPACT / "impact_acceleration_summary.json"
    impact_summary = json.loads(impact_summary_path.read_text()) if impact_summary_path.exists() else {}
    immuno_layers = read_tsv(IMMUNO / "immunogenicity_model_layers.tsv")
    immuno_labels = read_tsv(IMMUNO / "immunogenicity_label_contract.tsv")
    immuno_readiness = read_tsv(IMMUNO / "cross_neo_i_candidate_readiness.tsv")
    immuno_features = read_tsv(IMMUNO / "immunogenicity_feature_contract.tsv")
    immuno_benchmark = read_tsv(IMMUNO / "immunogenicity_benchmark_contract.tsv")
    immuno_summary_path = IMMUNO / "immunogenicity_prediction_upgrade_summary.json"
    immuno_summary = json.loads(immuno_summary_path.read_text()) if immuno_summary_path.exists() else {}
    physics_plan = read_tsv(PHYSICS / "physics_gate_candidate_plan.tsv")
    physics_budget = read_tsv(PHYSICS / "physics_gate_budget.tsv")
    physics_modules = read_tsv(PHYSICS / "physics_modules.tsv")
    physics_rules = read_tsv(PHYSICS / "physics_gate_rules.tsv")
    physics_summary_path = PHYSICS / "physics_gate_summary.json"
    physics_summary = json.loads(physics_summary_path.read_text()) if physics_summary_path.exists() else {}
    summary_path = PKG / "high_impact_summary.json"
    summary = json.loads(summary_path.read_text()) if summary_path.exists() else {}
    build_guides(top13, source_perf, wetlab)
    html_text = build_html(
        top13,
        source_perf,
        wetlab,
        claims,
        validation,
        risks,
        sim_matrix,
        sim_bridge,
        sim_rules,
        sim_budget,
        sim_summary,
        lead_external,
        lead_locks,
        lead_wetlab,
        evidence_matrix,
        figure_plan,
        lead_summary,
        assay_reagents,
        assay_plate,
        assay_thresholds,
        assay_moat,
        assay_summary,
        preclin_score,
        preclin_gates,
        preclin_failures,
        preclin_materials,
        preclin_summary,
        impact_positioning,
        impact_milestones,
        impact_validation,
        impact_scenarios,
        impact_actions,
        impact_summary,
        immuno_layers,
        immuno_labels,
        immuno_readiness,
        immuno_features,
        immuno_benchmark,
        immuno_summary,
        physics_plan,
        physics_budget,
        physics_modules,
        physics_rules,
        physics_summary,
        summary,
    )
    PAGE.write_text(html_text)
    WEB.mkdir(parents=True, exist_ok=True)
    safe_copy2(PAGE, WEB_PAGE)
    print(f"[high-impact-web] wrote {PAGE}")
    print(f"[high-impact-web] deploy target {WEB_PAGE}")


if __name__ == "__main__":
    main()
