#!/usr/bin/env python3
"""Build a standalone CROSS-Neo decision storyboard web page and guide files."""

from __future__ import annotations

import html
import shutil
from pathlib import Path

import pandas as pd


REPO = Path(__file__).resolve().parents[3]
MD_OUT = REPO / "project/results/cross_neo_md_audit_2026_05_10"
STORY = MD_OUT / "decision_storyboard"
FIG = MD_OUT / "figures"
HUB = REPO / "project/papers_hub_2026_05_04"
WEB = Path("/var/www/papers/papers_hub_2026_05_04")
ASSET = HUB / "assets/cross_neo_md_audit"
WEB_ASSET = WEB / "assets/cross_neo_md_audit"
PAGE = HUB / "cross_neo_decision_storyboard.html"
WEB_PAGE = WEB / "cross_neo_decision_storyboard.html"


def safe_copy2(src: Path, dst: Path) -> None:
    try:
        shutil.copy2(src, dst)
    except PermissionError:
        print(f"[storyboard-web] skip permission-denied asset {dst}")


FIGURE_CATALOG = [
    ("fig_md1_pipeline", "MD audit pipeline", "CROSS-Neo 후보에서 구조/MD audit까지 가는 전체 절차.", "MD run inventory; workflow rules", "Overview / methods"),
    ("fig_md2_peptide_rmsd", "Peptide RMSD over time", "peptide backbone RMSD로 trajectory 안정성을 확인.", "qc/peptide_rmsd.tsv", "MD QC"),
    ("fig_md3_pmhc_contact_heatmap", "pMHC contact occupancy", "peptide-HLA 접촉 occupancy heatmap.", "contacts/pmhc_contact_occupancy.tsv", "pMHC stability evidence"),
    ("fig_md4_tcr_peptide_contact_heatmap", "TCR-peptide contact occupancy", "TCR-peptide 접촉이 어느 peptide position에서 유지되는지.", "contacts/tcr_peptide_contact_occupancy.tsv", "TCR recognition plausibility"),
    ("fig_md5_mutant_wt_interface_delta", "WT/decoy boundary", "WT/decoy trajectory가 아직 필요한 부분과 claim boundary.", "counterfactual/counterfactual_todo_manifest.tsv", "Limit / control plan"),
    ("fig_md6_replicate_consistency", "Replicate consistency", "replicate 수와 RMSD consistency 요약.", "replicates/replicate_consistency_summary.tsv", "Robustness"),
    ("fig_md7_model_md_agreement", "Model-MD agreement", "DL/TCR 후보와 MD evidence alignment.", "md_evidence_scores.tsv; integrated/cross_neo_md_joined.tsv", "Model + MD bridge"),
    ("fig_md8_claim_boundary", "Claim boundary", "무엇을 말할 수 있고 무엇은 말하면 안 되는지.", "decision reports", "Reviewer defense"),
    ("fig_md9_live_progress_temperature", "Live progress and temperature", "OpenMM run progress, temperature, speed.", "md_status_summary.tsv; md_state_traces.tsv", "Run QC"),
    ("fig_md10_6uon_10ns_rmsd_contact_combo", "GADGVGKSAL 10 ns combo", "GADGVGKSAL RMSD/contact combined evidence.", "qc/peptide_rmsd.tsv; contacts/*.tsv", "Candidate case study"),
    ("fig_md11_anchor_tcr_contact_summary", "Anchor and TCR contact summary", "anchor/pMHC/TCR contact summary across runs.", "contacts/pmhc_contact_summary.tsv; contacts/tcr_contact_summary.tsv", "Structural audit"),
    ("fig_md12_evidence_score_components", "MD evidence score components", "MD score 구성요소별 기여도.", "md_evidence_scores.tsv", "Evidence scoring"),
    ("fig_md13_visual_summary_dashboard", "MD visual dashboard", "MD status, RMSD, contacts, evidence label dashboard.", "md_status_summary.tsv; md_evidence_scores.tsv", "Summary dashboard"),
    ("fig_md14_ultra_priority_ranking", "Ultra priority ranking", "DL/TCR/MD 통합 우선순위 ranking.", "ultra_priority/ultra_wetlab_priority_candidates.tsv", "Wetlab prioritization"),
    ("fig_md15_ultra_priority_evidence_blocks", "Ultra priority evidence blocks", "대표 후보별 evidence block.", "ultra_priority/*.tsv", "Candidate evidence"),
    ("fig_md16_small_dataset_uncertainty_funnel", "Small-dataset uncertainty funnel", "Bayesian/dropout uncertainty gate 결과.", "small_dataset_uncertainty/*.tsv", "Uncertainty gate"),
    ("fig_md17_bayesian_dropout_uncertainty", "Bayesian/dropout uncertainty", "posterior interval과 perturbation sensitivity.", "small_dataset_uncertainty/*.tsv", "Uncertainty diagnostics"),
    ("fig_md18_dl_first_candidate_funnel", "DL-first candidate funnel", "cheap model gate부터 TCR/MD escalation까지 후보 감소.", "dl_first_funnel/*.tsv", "Main funnel"),
    ("fig_md19_dl_tcr_rescue_uncertainty", "DL/TCR rescue uncertainty", "main DL 약하지만 TCR evidence가 살리는 후보를 설명.", "dl_first_funnel/*.tsv", "TCR rescue"),
    ("fig_md20_baker_rosetta_fallback_scores", "Baker/Rosetta fallback scores", "Rosetta 미설치 상황에서 static contact fallback filter.", "baker_rosetta_filter/*.tsv", "Cheap structural filter"),
    ("fig_md21_baker_contact_triage_map", "Baker contact triage map", "pMHC/TCR contact 기반 cheap triage map.", "baker_rosetta_filter/*.tsv", "Structure triage"),
    ("fig_md22_threshold_precision_recall_frontier", "Threshold precision-recall frontier", "slider threshold preset의 precision/recall tradeoff.", "dl_threshold_optimization/*.tsv", "Operating point"),
    ("fig_md23_threshold_preset_confusion", "Threshold preset confusion", "preset별 TP/TN/FP/FN/called count.", "dl_threshold_optimization/*.tsv", "Operating point"),
    ("fig_md24_story_overview_flow", "Story overview flow", "전체 후보가 어떤 gate를 거쳐 대표 후보로 좁혀지는지.", "decision_storyboard/*.tsv", "Opening overview"),
    ("fig_md25_data_usage_matrix", "Data usage matrix", "어떤 데이터가 DL/uncertainty/TCR/structure/MD/UI에 쓰이는지.", "decision_storyboard/data_usage_manifest.tsv", "Data provenance"),
    ("fig_md26_representative_candidate_board", "Representative candidate board", "GADGVGKSAL/HMTEVVRHC 대표 후보별 evidence panel.", "decision_storyboard/representative_candidate_evidence_table.tsv", "Candidate board"),
    ("fig_md27_no_fp_top13_candidate_evidence", "Top-13 no-FP candidate evidence", "strict current-label preset으로 남은 Top-13 후보 evidence heatmap.", "high_impact_decision_package/no_false_positive_top13_candidates.tsv", "High-impact wetlab triage"),
    ("fig_md28_source_stratified_preset_performance", "Source-stratified preset performance", "추천 preset이 source별로 어떤 TP/FP/FN behavior를 보이는지.", "high_impact_decision_package/preset_source_stratified_performance.tsv", "Reviewer robustness"),
    ("fig_md29_claim_ladder", "Claim ladder", "지금 가능한 claim, wetlab 후 가능한 claim, 금지 claim을 분리.", "high_impact_decision_package/claim_ladder.tsv", "Reviewer defense"),
    ("fig_md30_wetlab_validation_workflow", "Wetlab validation workflow", "peptide-HLA binding부터 WT/decoy-controlled assay까지의 go/no-go 흐름.", "high_impact_decision_package/wetlab_validation_plate_plan.tsv", "Validation plan"),
    ("fig_md31_immunogenicity_simulation_ladder", "Immunogenicity simulation ladder", "immunogenicity assay 전 어떤 simulation layer를 추가할지 단계화.", "immunogenicity_simulation_escalation/immunogenicity_simulation_matrix.tsv", "Simulation escalation"),
    ("fig_md32_candidate_simulation_budget", "Candidate simulation budget", "Top 후보별 요청 short-MD ns와 GPU-hour budget.", "immunogenicity_simulation_escalation/candidate_simulation_budget.tsv", "Resource planning"),
    ("fig_md33_assay_to_simulation_bridge", "Assay-to-simulation bridge", "wetlab readout과 simulation evidence가 각각 무엇을 support하는지 연결.", "immunogenicity_simulation_escalation/wetlab_assay_simulation_bridge.tsv", "Assay interpretation"),
    ("fig_md34_immunogenicity_simulation_readiness", "Immunogenicity simulation readiness", "ready/queued/blocked/hold module 상태.", "immunogenicity_simulation_escalation/immunogenicity_simulation_matrix.tsv", "Execution dashboard"),
    ("fig_md35_two_lead_evidence_moat", "Two-lead evidence moat", "GADGVGKSAL/HMTEVVRHC 두 flagship 후보의 external/TCR/MD/control/actionability moat.", "two_lead_impact_package/two_lead_external_evidence.tsv", "Flagship lead case"),
    ("fig_md36_specificity_lock_matrix", "Specificity lock matrix", "mutation/HLA/structure/MD/WT/decoy/assay lock 상태.", "two_lead_impact_package/two_lead_specificity_locks.tsv", "Specificity defense"),
    ("fig_md37_impact_upgrade_flow", "Impact upgrade flow", "ranking에서 controlled immunogenicity test plan으로 넘어가는 흐름.", "two_lead_impact_package/two_lead_wetlab_order.tsv", "Impact framing"),
    ("fig_md38_reagent_order_sheet", "Reagent order sheet", "mutant/WT/decoy/pHLA reagent order row를 lead별로 정리.", "assay_ready_translation_package/assay_reagent_order_sheet.tsv", "Assay-ready execution"),
    ("fig_md39_control_coverage", "Control coverage", "flagship lead마다 어떤 control이 명시되어 있는지.", "assay_ready_translation_package/assay_reviewer_defense_moat.tsv", "Assay control defense"),
    ("fig_md40_assay_go_no_go_ladder", "Assay go/no-go ladder", "presentation, multimer, activation, killing gate 순서와 stop rule.", "assay_ready_translation_package/assay_go_no_go_thresholds.tsv", "Wetlab decision ladder"),
    ("fig_md41_reviewer_defense_moat", "Reviewer defense moat", "reviewer objection별 대응 asset을 연결.", "assay_ready_translation_package/assay_reviewer_defense_moat.tsv", "Reviewer defense"),
    ("fig_md42_preclinical_readiness_scorecard", "Preclinical readiness scorecard", "lead 후보별 preclinical readiness를 gate별로 점수화.", "preclinical_validation_protocol/preclinical_readiness_scorecard.tsv", "Preclinical validation"),
    ("fig_md43_preclinical_validation_gates", "Preclinical validation gates", "preclinical validation gate와 go/no-go 조건.", "preclinical_validation_protocol/preclinical_validation_gates.tsv", "Preclinical validation"),
    ("fig_md44_failure_mode_action_map", "Failure-mode action map", "실패 모드별 다음 action을 정리.", "preclinical_validation_protocol/preclinical_failure_mode_actions.tsv", "Risk/action map"),
    ("fig_md45_preclinical_translation_flow", "Preclinical translation flow", "assay-ready 패키지에서 preclinical validation으로 넘어가는 흐름.", "preclinical_validation_protocol/preclinical_validation_plan.tsv", "Translation flow"),
    ("fig_md46_evidence_to_claim_matrix", "Evidence-to-claim matrix", "각 evidence layer가 어떤 claim을 support하고 어떤 overclaim을 막는지.", "high_impact_decision_package/evidence_to_claim_matrix.tsv", "Reviewer claim boundary"),
    ("fig_md47_claim_upgrade_milestone_map", "Claim upgrade milestones", "강한 claim으로 올라가기 위해 통과해야 하는 milestone.", "impact_acceleration_package/claim_upgrade_milestones.tsv", "Impact acceleration"),
    ("fig_md48_reviewer_attack_defense_heatmap", "Reviewer attack-defense heatmap", "reviewer 공격 포인트별 방어 asset coverage.", "impact_acceleration_package/editor_reviewer_positioning.tsv", "Reviewer defense"),
    ("fig_md49_external_validation_benchmark_design", "External validation benchmark design", "retrospective audit에서 external validation으로 가는 benchmark 축.", "impact_acceleration_package/external_validation_benchmark_plan.tsv", "External validation"),
    ("fig_md50_next_72h_action_board", "Next 72h action board", "가장 빠르게 임팩트를 올릴 action board.", "impact_acceleration_package/next_72h_action_board.tsv", "Action board"),
    ("fig_md51_cross_neo_i_architecture", "CROSS-Neo-I architecture", "presentation/foreignness/TCR/MD/context/calibration으로 이루어진 면역원성 예측 구조.", "immunogenicity_prediction_upgrade/immunogenicity_model_layers.tsv", "Immunogenicity prediction"),
    ("fig_md52_label_feature_contract", "Label-feature contract", "어떤 label이 어떤 model layer를 학습시키는지.", "immunogenicity_prediction_upgrade/immunogenicity_label_contract.tsv", "Model contract"),
    ("fig_md53_cross_neo_i_readiness_ranking", "CROSS-Neo-I readiness ranking", "Top 후보별 면역원성 label-generation readiness.", "immunogenicity_prediction_upgrade/cross_neo_i_candidate_readiness.tsv", "Candidate readiness"),
    ("fig_md54_immunogenicity_benchmark_ladder", "Immunogenicity benchmark ladder", "면역원성 예측 claim을 위해 필요한 heldout/prospective benchmark 단계.", "immunogenicity_prediction_upgrade/immunogenicity_benchmark_contract.tsv", "Benchmark contract"),
    ("fig_md55_immunogenicity_risk_control_board", "Immunogenicity risk-control board", "leakage/source-shift/claim-risk control board.", "immunogenicity_prediction_upgrade/immunogenicity_benchmark_contract.tsv", "Risk control"),
]


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t") if path.exists() else pd.DataFrame()


def esc(x: object) -> str:
    return html.escape("" if pd.isna(x) else str(x))


def fmt(x: object, digits: int = 3) -> str:
    try:
        if pd.isna(x):
            return "NA"
        return f"{float(x):.{digits}f}"
    except Exception:
        return esc(x)


def table(df: pd.DataFrame, cols: list[str], max_rows: int = 20) -> str:
    if df.empty:
        return "<p class='muted'>No rows.</p>"
    cols = [c for c in cols if c in df.columns]
    rows = ["<table><thead><tr>"]
    rows += [f"<th>{esc(c.replace('_', ' '))}</th>" for c in cols]
    rows.append("</tr></thead><tbody>")
    for _, r in df.head(max_rows).iterrows():
        rows.append("<tr>")
        for c in cols:
            v = r.get(c, "")
            rows.append(f"<td>{fmt(v) if isinstance(v, (int, float)) else esc(v)}</td>")
        rows.append("</tr>")
    rows.append("</tbody></table>")
    return "\n".join(rows)


def md_summary_for(reps: pd.DataFrame, peptide: str, hla: str) -> tuple[str, str]:
    if reps.empty:
        return "NA", "NA"
    mask = (reps["peptide"].astype(str).str.upper() == peptide.upper()) & (
        reps["hla_4digit"].astype(str).str.upper() == hla.upper()
    )
    if not mask.any():
        return "NA", "NA"
    row = reps[mask].iloc[0]
    label = row.get("current_md_evidence_label") or row.get("md_label") or "NA"
    try:
        score = f"{float(row.get('current_md_evidence_score')):.3f}"
    except Exception:
        score = "NA"
    return str(label), score


def build_slide_plan(stage: pd.DataFrame, presets: pd.DataFrame, reps: pd.DataFrame) -> pd.DataFrame:
    nofp = presets[presets["preset_name"].eq("NO_FALSE_POSITIVE_MAX_TP")] if not presets.empty and "preset_name" in presets.columns else pd.DataFrame()
    nofp_msg = "13 TP / 0 FP preset" if nofp.empty else f"{int(nofp.iloc[0]['TP'])} TP / {int(nofp.iloc[0]['FP'])} FP preset"
    hm_label, hm_score = md_summary_for(reps, "HMTEVVRHC", "HLA-A*02:01")
    gad_label, gad_score = md_summary_for(reps, "GADGVGKSAL", "HLA-C*08:02")
    rows = [
        {
            "slide": 1,
            "title": "Problem: presentation is not recognition",
            "figure_or_table": "fig_md24_story_overview_flow",
            "what_to_say": "CROSS-Neo starts with cheap prediction and only later escalates to TCR/structure/MD. This prevents spending MD/wetlab resources on broad noisy candidates.",
            "data_used": "TCR wetlab candidate table; DL-first stage counts",
            "decision": "Use as opening overview.",
        },
        {
            "slide": 2,
            "title": "Candidate universe and labels",
            "figure_or_table": "Table 1 stage_label_composition_for_slides",
            "what_to_say": "The starting set has 649 candidates, 161 positive labels and 488 negative labels. The question is not positive-only discovery; it is robust positive triage plus negative/false-positive removal.",
            "data_used": "tcr_wetlab_candidate_prioritization_unique_pmhc.tsv",
            "decision": "Show label reality before model results.",
        },
        {
            "slide": 3,
            "title": "DL-first funnel",
            "figure_or_table": "fig_md18_dl_first_candidate_funnel; slider dashboard",
            "what_to_say": "Cheap DL and uncertainty reduce candidates before any expensive structure/MD step.",
            "data_used": "all_predictions.tsv; Bayesian/dropout uncertainty tables",
            "decision": "Main model operating flow.",
        },
        {
            "slide": 4,
            "title": "Threshold optimizer",
            "figure_or_table": "fig_md22_threshold_precision_recall_frontier; Table 2 presets",
            "what_to_say": f"Instead of hand-tuning, labels are used to find operating presets. The strongest current no-FP setting gives {nofp_msg}.",
            "data_used": "threshold_random_search_results.tsv; threshold_recommended_presets.tsv",
            "decision": "Use this to argue the UI is decision-support, not arbitrary sliders.",
        },
        {
            "slide": 5,
            "title": "Data lineage",
            "figure_or_table": "fig_md25_data_usage_matrix; Table 3 data usage",
            "what_to_say": "Every score and figure is traceable to a specific table, with a claim boundary.",
            "data_used": "data_usage_manifest.tsv",
            "decision": "Reproducibility slide.",
        },
        {
            "slide": 6,
            "title": "Representative candidates",
            "figure_or_table": "fig_md26_representative_candidate_board",
            "what_to_say": f"HMTEVVRHC now has completed 10 ns {hm_label} structural evidence (score {hm_score}) despite weaker main DL; GADGVGKSAL remains the balanced DL+TCR+MD candidate ({gad_label}, score {gad_score}) with a completed 10 ns replicate/control MD run.",
            "data_used": "dl_first_md_escalation_candidates.tsv; Baker fallback scores; md_evidence_scores.tsv; MD status",
            "decision": "Use as candidate prioritization and mechanism-audit slide.",
        },
        {
            "slide": 7,
            "title": "Claim boundary",
            "figure_or_table": "fig_md8_claim_boundary",
            "what_to_say": "The pipeline prioritizes wetlab candidates. It does not prove immunogenicity or clinical utility.",
            "data_used": "decision reports; MD reports",
            "decision": "Use this as reviewer defense.",
        },
    ]
    return pd.DataFrame(rows)


def build_caption_guides(figtab: pd.DataFrame, data_usage: pd.DataFrame, slide_plan: pd.DataFrame) -> tuple[str, str]:
    fig_lines = [
        "# CROSS-Neo Figure 설명 가이드",
        "",
        "## 전체 메시지",
        "",
        "CROSS-Neo는 cheap DL/uncertainty로 후보를 먼저 줄이고, TCR/구조/MD는 마지막 검증/해석 계층으로만 쓴다.",
        "",
    ]
    for _, r in figtab.iterrows():
        fig_lines += [
            f"## {r['id']} — {r['title']}",
            "",
            f"- 보여주는 것: {r['what_it_shows']}",
            f"- 핵심 메시지: {r['main_message']}",
            f"- 쓰는 위치: {r['use_in_slide']}",
            f"- 주의 문장: {r['claim_boundary']}",
            "",
        ]

    table_lines = [
        "# CROSS-Neo Table / Data 설명 가이드",
        "",
        "## 데이터 사용 원칙",
        "",
        "각 데이터는 특정 계층에서만 쓰며, label/score/MD/structure evidence를 서로 ground truth처럼 섞지 않는다.",
        "",
        "## Slide Talk Track",
        "",
        slide_plan.to_markdown(index=False),
        "",
        "## Data Usage",
        "",
        data_usage[["artifact", "rows", "used_for", "pipeline_stage", "claim_boundary"]].to_markdown(index=False),
        "",
    ]
    return "\n".join(fig_lines) + "\n", "\n".join(table_lines) + "\n"


def build_all_figure_catalog(figtab: pd.DataFrame, data_usage: pd.DataFrame, reps: pd.DataFrame) -> tuple[pd.DataFrame, str]:
    rows = []
    available = {p.stem for p in FIG.glob("fig_md*.png")}
    for fig_id, title, explanation, source_files, use_case in FIGURE_CATALOG:
        rows.append(
            {
                "figure_id": fig_id,
                "png_path": f"figures/{fig_id}.png",
                "exists": fig_id in available,
                "title": title,
                "what_it_shows_kr": explanation,
                "source_data": source_files,
                "where_to_use": use_case,
                "claim_boundary": "prioritization/structural audit evidence; not immunogenicity proof",
            }
        )
    catalog = pd.DataFrame(rows)
    hm_label, hm_score = md_summary_for(reps, "HMTEVVRHC", "HLA-A*02:01")
    gad_label, gad_score = md_summary_for(reps, "GADGVGKSAL", "HLA-C*08:02")
    lines = [
        "# CROSS-Neo 전체 Figure / Table 설명 패키지",
        "",
        "## 전체 흐름",
        "",
        "1. 후보 universe와 label 현실을 먼저 보여준다.",
        "2. cheap DL/uncertainty gate로 후보를 줄인다.",
        "3. TCR branch와 structure/MD layer로 recognition plausibility를 본다.",
        "4. 대표 후보 보드에서 GADGVGKSAL/HMTEVVRHC를 비교한다.",
        "5. data usage matrix로 어떤 데이터가 어디에 쓰였는지 방어한다.",
        "",
        "## 현재 핵심 업데이트",
        "",
        f"- HMTEVVRHC/HLA-A*02:01: 10 ns explicit-solvent MD 완료, `{hm_label}` (score {hm_score}).",
        f"- GADGVGKSAL/HLA-C*08:02: 10 ns explicit-solvent MD 완료, `{gad_label}` (score {gad_score}); 추가 10 ns replicate/control MD도 완료.",
        "",
        "## 전체 Figure Catalog",
        "",
        catalog[["figure_id", "exists", "title", "what_it_shows_kr", "source_data", "where_to_use", "claim_boundary"]].to_markdown(index=False),
        "",
        "## 대표 Figure/Table Master",
        "",
        figtab[["id", "title", "source_files", "main_message", "use_in_slide", "claim_boundary"]].to_markdown(index=False),
        "",
        "## Data Usage Manifest",
        "",
        data_usage[["artifact", "rows", "path", "used_for", "pipeline_stage", "claim_boundary"]].to_markdown(index=False),
    ]
    return catalog, "\n".join(lines) + "\n"


def build_html_page(stage: pd.DataFrame, presets: pd.DataFrame, data_usage: pd.DataFrame, figtab: pd.DataFrame, reps: pd.DataFrame, slide_plan: pd.DataFrame) -> str:
    fig_cards = [
        ("fig_md24_story_overview_flow.png", "Overview Flow", "전체 후보가 어떤 gate를 거쳐 줄어드는지 보여주는 대표 그림."),
        ("fig_md18_dl_first_candidate_funnel.png", "DL-First Funnel", "cheap DL/uncertainty stage별 후보 감소."),
        ("fig_md22_threshold_precision_recall_frontier.png", "Threshold Frontier", "slider threshold operating point 탐색 결과."),
        ("fig_md23_threshold_preset_confusion.png", "Preset Confusion", "추천 preset별 TP/FP/FN/called count."),
        ("fig_md25_data_usage_matrix.png", "Data Usage Matrix", "어떤 데이터가 어느 계층에 쓰였는지."),
        ("fig_md26_representative_candidate_board.png", "Candidate Board", "대표 후보별 DL/Bayes/TCR/structure/MD evidence."),
        ("fig_md20_baker_rosetta_fallback_scores.png", "Baker Fallback", "Rosetta 전 static contact 구조 sanity check."),
        ("fig_md13_visual_summary_dashboard.png", "MD Visual Summary", "MD evidence dashboard."),
    ]
    cards_html = "\n".join(
        f"<article><img src='assets/cross_neo_md_audit/{esc(name)}'><h3>{esc(title)}</h3><p>{esc(desc)}</p><a href='assets/cross_neo_md_audit/{esc(name)}'>PNG</a></article>"
        for name, title, desc in fig_cards
    )
    return f"""<!doctype html>
<html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>CROSS-Neo Decision Storyboard</title>
<style>
body{{margin:0;background:#07111f;color:#edf5ff;font-family:Inter,Arial,sans-serif;line-height:1.55}}
a{{color:#39d4b5}} header{{padding:38px 32px;background:#0f2035;border-bottom:1px solid #26364d}}
.wrap{{max-width:1320px;margin:0 auto;padding:24px 30px 70px}} h1{{font-family:Georgia,serif;font-size:44px;margin:0 0 10px}}
h2{{font-family:Georgia,serif;color:#fff2d0;font-size:28px}} .lead{{color:#cbd7e7;max-width:1040px}}
.stats{{display:grid;grid-template-columns:repeat(5,1fr);gap:12px;margin-top:18px}} .stat{{background:#101d30;border:1px solid #293b55;border-radius:8px;padding:13px}} .stat b{{display:block;color:#f2c46d;font-size:25px}}
section{{border-bottom:1px solid #26364d;padding:22px 0}} table{{border-collapse:collapse;width:100%;font-size:12px}} th,td{{border:1px solid #293b55;padding:7px;vertical-align:top}} th{{background:#13243a;color:#f2c46d}}
.figs{{display:grid;grid-template-columns:repeat(2,1fr);gap:16px}} article{{background:#101d30;border:1px solid #293b55;border-radius:8px;padding:12px}} img{{width:100%;background:white;border-radius:6px}}
.links{{display:flex;flex-wrap:wrap;gap:10px}} .pill{{border:1px solid #395170;border-radius:999px;padding:7px 10px;background:#101d30}} .muted{{color:#9fb0c7}}
@media(max-width:900px){{.stats,.figs{{grid-template-columns:1fr}}}}
</style></head><body>
<header><div class="wrap">
<h1>CROSS-Neo Decision Storyboard</h1>
<p class="lead">A presentation-ready explanation package for the DL-first, TCR-aware, structure-audited neoantigen prioritization workflow.</p>
<div class="stats">
<div class="stat"><b>649</b><span>candidate rows</span></div>
<div class="stat"><b>161 / 488</b><span>positive / negative labels</span></div>
<div class="stat"><b>13 / 0</b><span>optimizer TP / FP preset</span></div>
<div class="stat"><b>2</b><span>MD escalation rows</span></div>
<div class="stat"><b>1</b><span>current wetlab shortlist</span></div>
</div>
</div></header>
<main class="wrap">
<section><h2>Files</h2><div class="links">
<a class="pill" href="cross_neo_dl_first_slider_dashboard.html">Interactive slider dashboard</a>
<a class="pill" href="cross_neo_md_audit_dossier.html">Full MD/TCR dossier</a>
<a class="pill" href="assets/cross_neo_md_audit/CROSS_Neo_decision_storyboard_overview_KR.md">Korean overview</a>
<a class="pill" href="assets/cross_neo_md_audit/figure_caption_guide_KR.md">Figure guide KR</a>
<a class="pill" href="assets/cross_neo_md_audit/table_data_explanation_guide_KR.md">Table/data guide KR</a>
<a class="pill" href="assets/cross_neo_md_audit/figure_table_explanation_master.tsv">Figure/Table master TSV</a>
<a class="pill" href="assets/cross_neo_md_audit/data_usage_manifest.tsv">Data usage TSV</a>
</div></section>
<section><h2>Recommended Slide Flow</h2>{table(slide_plan, ['slide','title','figure_or_table','what_to_say','data_used','decision'], 20)}</section>
<section><h2>Stage Label Composition</h2>{table(stage, ['stage','n','positive','negative','positive_fraction'], 20)}</section>
<section><h2>Threshold Presets</h2>{table(presets, ['preset_name','call_rule','called_positive','TP','TN','FP','FN','precision','recall','F1'], 10)}</section>
<section><h2>Representative Candidate Evidence</h2>{table(reps, ['row_id','peptide','hla_4digit','label_binary','dl_first_decision','main_dl_score','bayes_mean','tcr_augmented_score_mean','baker_structural_score','current_md_evidence_label','current_md_evidence_score','current_md_run_id','current_md_peptide_rmsd_final_nm','current_md_tcr_peptide_contacts_tail','completion_fraction'], 8)}</section>
<section><h2>Figure / Table Explanation Master</h2>{table(figtab, ['id','title','main_message','use_in_slide','claim_boundary'], 12)}</section>
<section><h2>Data Usage Manifest</h2>{table(data_usage, ['artifact','rows','used_for','pipeline_stage','claim_boundary'], 12)}</section>
<section><h2>Figure Gallery</h2><div class="figs">{cards_html}</div></section>
<section><h2>Claim Boundary</h2><p class="muted">This is a decision-support and prioritization package. It does not prove immunogenicity, clinical utility, or universal TCR-aware SOTA.</p></section>
</main></body></html>"""


def main() -> None:
    STORY.mkdir(parents=True, exist_ok=True)
    ASSET.mkdir(parents=True, exist_ok=True)
    WEB_ASSET.mkdir(parents=True, exist_ok=True)
    stage = read_tsv(STORY / "stage_label_composition_for_slides.tsv")
    presets = read_tsv(STORY / "threshold_preset_summary_for_slides.tsv")
    data_usage = read_tsv(STORY / "data_usage_manifest.tsv")
    figtab = read_tsv(STORY / "figure_table_explanation_master.tsv")
    reps = read_tsv(STORY / "representative_candidate_evidence_table.tsv")
    slide_plan = build_slide_plan(stage, presets, reps)
    slide_plan.to_csv(STORY / "slide_order_and_talk_track.tsv", sep="\t", index=False)
    fig_guide, table_guide = build_caption_guides(figtab, data_usage, slide_plan)
    all_catalog, all_overview = build_all_figure_catalog(figtab, data_usage, reps)
    all_catalog.to_csv(STORY / "all_md_figure_catalog.tsv", sep="\t", index=False)
    (STORY / "figure_caption_guide_KR.md").write_text(fig_guide)
    (STORY / "table_data_explanation_guide_KR.md").write_text(table_guide)
    (STORY / "ALL_FIGURE_TABLE_OVERVIEW_KR.md").write_text(all_overview)
    html_text = build_html_page(stage, presets, data_usage, figtab, reps, slide_plan)
    PAGE.write_text(html_text)
    WEB.mkdir(parents=True, exist_ok=True)
    safe_copy2(PAGE, WEB_PAGE)
    for p in [
        STORY / "slide_order_and_talk_track.tsv",
        STORY / "figure_caption_guide_KR.md",
        STORY / "table_data_explanation_guide_KR.md",
        STORY / "all_md_figure_catalog.tsv",
        STORY / "ALL_FIGURE_TABLE_OVERVIEW_KR.md",
    ]:
        safe_copy2(p, ASSET / p.name)
        safe_copy2(p, WEB_ASSET / p.name)
    print(f"[storyboard-web] wrote {PAGE}")
    print(f"[storyboard-web] deploy target {WEB_PAGE}")


if __name__ == "__main__":
    main()
