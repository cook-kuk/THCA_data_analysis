#!/usr/bin/env python3
from __future__ import annotations

import html
import shutil
from pathlib import Path

import numpy as np
import pandas as pd

import build_paper1_style_visual_atlas as atlas


HTML_NAME = "kthyro_ctc_emt_representative_overview.html"
REPORT_NAME = "REPRESENTATIVE_OVERVIEW_GUIDE_KR.md"
FIG_TABLE_NAME = "representative_figure_overview.tsv"
TAB_TABLE_NAME = "representative_table_overview.tsv"
DATA_MAP_NAME = "data_to_figure_claim_map.tsv"
OVERVIEW_FIG = "PF15_grand_overview_data_to_claim_map.png"
ASSET = f"assets/{atlas.ASSET}"


def esc(x) -> str:
    return html.escape("" if x is None else str(x), quote=True)


FIGURE_ROWS = [
    {
        "slot": "Grand Overview",
        "representative_asset": OVERVIEW_FIG,
        "figure_title": "Data-to-Claim Grand Overview",
        "role": "전체 흐름 그림. clinical perturbation, public-data priors, representative figures, control tables, current/future claims를 한 장에 연결한다.",
        "data_used": "All local evidence tables; public-pilot output summaries; Yu 2024 published anchor; future Samsung/Yu cohort design.",
        "where_used": "교수님 첫 설명, proposal opening schematic, manuscript overview figure.",
        "claim": "This is a staged Science case: public data builds priors; prospective serial CTC/cfDNA/tissue NGS tests the thyroidectomy perturbation question.",
        "claim_boundary": "Schematic is not a new statistical result and does not claim postoperative CTC transition has already been discovered.",
    },
    {
        "slot": "Main Fig. 1",
        "representative_asset": "PF01_manuscript_flow.png;PF02_data_provenance_map.png",
        "figure_title": "Study Logic And Evidence Provenance",
        "role": "Yu 2024, public-data priors, and future hospital cohort를 구분한다.",
        "data_used": "Yu 2024 published CTC feasibility; TCGA/spatial/scRNA/proteomics public pilot outputs; proposal design tables.",
        "where_used": "Introduction end, proposal Slide 3-4, reviewer response to 'what data do you have?'",
        "claim": "Evidence layers have distinct roles: public priors, published feasibility, future prospective test.",
        "claim_boundary": "Yu 2024 is not local raw data; public data is not used as CTC proof.",
    },
    {
        "slot": "Main Fig. 2",
        "representative_asset": "F03_tcga_braf_state_split.png;F04_tcga_driver_variance_boundary.png",
        "figure_title": "Driver Mutation Alone Is Not Enough",
        "role": "BRAF-mutant PTC도 하나의 biological state로 접히지 않음을 보여준다.",
        "data_used": "TCGA-THCA bulk RNA/mutation/clinical public-pilot tables; BRAF-only vulnerability label distribution; axis variance models.",
        "where_used": "Results 1, Samsung Slide 5, rationale for phenotype + genotype integration.",
        "claim": "Public PTC tumor states split beyond mutation group; genotype alone is insufficient as a state descriptor.",
        "claim_boundary": "Does not prove CTC shedding, postoperative transition, or recurrence prediction.",
    },
    {
        "slot": "Main Fig. 3",
        "representative_asset": "F05_spatial_tissue_state_organization.png;spatial_coherence_summary.png",
        "figure_title": "Spatial Tissue-State Organization",
        "role": "조직 state가 random이 아니라 spatially organized되어 marker/ROI logic을 지지한다.",
        "data_used": "GSE250521 spatial transcriptomics; 16 slides; 57,144 spots; coherence/adjacency/hotspot/interface summaries.",
        "where_used": "Results 2, tissue-context rationale, reviewer response to marker-origin logic.",
        "claim": "Spatial data support organized tissue-state context for marker selection.",
        "claim_boundary": "Spatial organization is not direct proof of CTC origin or shedding source.",
    },
    {
        "slot": "Main Fig. 4",
        "representative_asset": "F06_scrna_marker_context.png;F07_proteomics_dediff_direction.png",
        "figure_title": "Marker Module Cross-Layer Support",
        "role": "Epithelial, EM/M, thyroid-lineage, survival/stress marker modules를 tissue/cell/protein context에서 정리한다. PF10은 별도 supplementary planning matrix.",
        "data_used": "PTC scRNA marker-context summary; bulk proteomics proxy. (PF10 supplementary planning matrix은 별도 슬롯.)",
        "where_used": "Results 3, CTC marker-panel design, Inocras/vendor assay discussion.",
        "claim": "Candidate CTC-EMT marker modules are biologically interpretable across public layers.",
        "claim_boundary": "Marker support is hypothesis-generating; not a validated clinical CTC assay.",
    },
    {
        "slot": "Supplementary Planning Matrix",
        "representative_asset": "PF10_marker_module_support_heatmap.png",
        "figure_title": "PF10 - Marker Module Planning Support Scores (hand-coded; not data-derived)",
        "role": "Marker module별 prior 구성 유용성을 0/1/2로 손코딩한 manuscript-planning matrix. 분석 산출물이 아닌 기획 점수표다.",
        "data_used": "marker_module_cross_layer_support_matrix.tsv (hand-coded). TCGA-THCA tissue (no CTC); GSE250521 spatial; scRNA module summary; bulk proteomics proxy; Yu 2024 published anchor summary (no unpublished/raw Yu CTC data).",
        "where_used": "Supplement only. Reviewer-visible planning artifact.",
        "claim": "수기로 매긴 manuscript-planning judgement이다.",
        "claim_boundary": "측정값/validation metric/cross-layer 분석 결과가 아니다. TCGA는 CTC 데이터를 포함하지 않으며 Yu 2024 컬럼은 published summary만 인용한 것이다.",
    },
    {
        "slot": "Main Fig. 5",
        "representative_asset": "PF03_marker_to_ngs_map.png;F08_ctc_emt_ngs_prior_panel.png",
        "figure_title": "Phenotype-To-Genotype NGS Bridge",
        "role": "CTC phenotype만으로 부족하므로 matched tumor-normal NGS/cfDNA anchoring이 필요함을 보여준다.",
        "data_used": "CTC-EMT marker prior table; marker-to-NGS map; Samsung proposal NGS logic.",
        "where_used": "Results 4, proposal Aim 2, Inocras questionnaire, assay workflow figure.",
        "claim": "CTC-EMT phenotype should be interpreted with genetic anchoring.",
        "claim_boundary": "Does not assume cfDNA/CTC-enriched NGS will work in every early PTC patient.",
    },
    {
        "slot": "Main Fig. 6",
        "representative_asset": "PF04_dataset_gap_bridge.png;PF12_future_data_unlock_map.png;PF05_claim_boundary_matrix.png",
        "figure_title": "Missing Dataset And Future Data Unlock",
        "role": "왜 병원/Yu serial CTC/cfDNA/tissue NGS 데이터가 필요한지 figure-driven으로 설명한다.",
        "data_used": "serial CTC-NGS dataset gap map; future data unlock table; claim boundary matrix; hospital CRF.",
        "where_used": "Discussion/future prospective study, Samsung final ask, Prof. Yu data request.",
        "claim": "No public dataset currently tests the full serial postoperative CTC-EMT genetic-map question.",
        "claim_boundary": "Dataset absence justifies the prospective study; it is not itself a biological discovery.",
    },
    {
        "slot": "Defense/Supp",
        "representative_asset": "PF07_reviewer_to_figure_map.png;PF14_analysis_control_board.png;F09_publishability_decision.png",
        "figure_title": "Reviewer Defense And Execution Control",
        "role": "예상 공격, kill criteria, safe claim sentence를 한 번에 방어한다.",
        "data_used": "reviewer risk heatmap; statistical analysis plan; stage-gate kill criteria; safe sentence bank.",
        "where_used": "Supplement, internal review, Samsung Q&A rehearsal.",
        "claim": "The proposal is ambitious but bounded and stage-gated.",
        "claim_boundary": "Defense artifacts are not new biological results.",
    },
]


TABLE_ROWS = [
    ("data_dictionary_master.tsv", "Data Definition Master", "각 데이터셋이 무엇이고 무엇에 쓰였는지 정의", "Methods/Data section; reviewer data provenance", "forbidden_claim 컬럼을 반드시 같이 본다."),
    ("figure_caption_and_results_scaffold.tsv", "Caption And Results Scaffold", "각 figure의 caption, results paragraph, methods note, reviewer sentence", "Figure legend drafting; manuscript Results scaffold", "문장 scaffold이지 final claim expansion이 아니다."),
    ("representative_figure_overview.tsv", "Representative Figure Overview", "대표 figure spine과 data/claim/boundary", "교수님 briefing; manuscript figure plan", "이 파일 자체는 guide다."),
    ("data_to_figure_claim_map.tsv", "Data-To-Figure-Claim Map", "데이터별 어디에 썼고 어떤 claim을 허용하는지", "Methods, rebuttal, grant review defense", "future data rows는 current result가 아니다."),
    ("claim_boundary_matrix.tsv", "Claim Boundary Matrix", "Allowed/forbidden statement를 명시", "Reviewer defense and final editing", "forbidden 문장은 본문에 쓰지 않는다."),
    ("hospital_data_request_crf.tsv", "Hospital Data Request CRF", "Yu/hospital에서 받아야 할 필드", "7-day execution, IRB/CRF", "요청 필드이지 이미 확보한 데이터가 아니다."),
    ("statistical_analysis_plan.tsv", "Statistical Analysis Plan", "현재 public-data analysis와 future prospective test 분리", "Methods and Samsung review", "future SAP는 현재 결과가 아니다."),
    ("stage_gate_kill_criteria.tsv", "Stage-Gate Kill Criteria", "assay/NGS/follow-up 실패 시 downgrade rule", "Budget/risk management", "kill criteria를 숨기지 않는다."),
    ("safe_sentence_bank.tsv", "Safe Sentence Bank", "본문에 바로 쓸 수 있는 안전 문장", "Writing and reviewer answer", "safe sentence 외 확장 시 claim boundary 재확인."),
]


DATA_ROWS = [
    {
        "dataset_or_source": "TCGA-THCA public pilot",
        "data_type": "bulk tumor RNA/mutation/clinical public cohort",
        "used_in_figures": "Main Fig. 2; Evidence Figures; Slide 5",
        "used_in_tables": "data_dictionary_master; key_number_ledger; data_to_figure_claim_map",
        "allowed_use": "PTC tissue-state heterogeneity beyond driver mutation; BRAF-mutant state split; genotype-alone insufficiency rationale.",
        "not_allowed": "CTC shedding, postoperative CTC transition, recurrence prediction, clinical diagnostic claim.",
    },
    {
        "dataset_or_source": "GSE250521 spatial transcriptomics",
        "data_type": "public thyroid spatial transcriptomics",
        "used_in_figures": "Main Fig. 3; Evidence Figures; spatial coherence panels",
        "used_in_tables": "data_dictionary_master; representative_figure_overview",
        "allowed_use": "Spatial tissue-state organization; marker/ROI context rationale.",
        "not_allowed": "CTC origin proof or direct tumor-cell shedding source.",
    },
    {
        "dataset_or_source": "PTC scRNA marker context",
        "data_type": "single-cell marker module summary",
        "used_in_figures": "Main Fig. 4 (F06); Supplementary Planning Matrix (PF10 - hand-coded planning score, not data-derived)",
        "used_in_tables": "marker_module_cross_layer_support_matrix; data_dictionary_master",
        "allowed_use": "Cell/tissue context for epithelial, mesenchymal, thyroid-lineage, stress/survival modules.",
        "not_allowed": "Definition of circulating CTC state without blood data.",
    },
    {
        "dataset_or_source": "Bulk proteomics proxy",
        "data_type": "protein-direction support layer",
        "used_in_figures": "Main Fig. 4; F07; proteomics supplementary figures",
        "used_in_tables": "marker_module_cross_layer_support_matrix",
        "allowed_use": "Protein-direction support for dedifferentiation/marker-prior logic.",
        "not_allowed": "Direct CTC protein phenotype or clinical assay validation.",
    },
    {
        "dataset_or_source": "Yu 2024 published CTC study",
        "data_type": "published prospective PTC serial CTC feasibility anchor",
        "used_in_figures": "Main Fig. 1; Breakthrough Case; Slide 3; proposal text",
        "used_in_tables": "data_dictionary_master; claim_boundary_matrix",
        "allowed_use": "Feasibility anchor: serial PTC CTC-EMT measurement around thyroidectomy has been reported.",
        "not_allowed": "Local raw reanalysis; genetic-map evidence; our cohort result.",
    },
    {
        "dataset_or_source": "Future Yu/hospital serial CTC dataset",
        "data_type": "prospective blood timepoints T0/T1/T2/T3",
        "used_in_figures": "Future Main Fig. serial transition; Samsung Aim 1; Future Data Unlock Map",
        "used_in_tables": "hospital_data_request_crf; statistical_analysis_plan; stage_gate_kill_criteria",
        "allowed_use": "Required to test postoperative CTC-EMT state transition.",
        "not_allowed": "Do not write as existing result until data are obtained and QC passes.",
    },
    {
        "dataset_or_source": "Matched tumor-normal NGS",
        "data_type": "future tissue/blood genetic anchor",
        "used_in_figures": "Main Fig. 5; Samsung Aim 2; NGS genetic map",
        "used_in_tables": "Inocras questionnaire; CRF; SAP; kill criteria",
        "allowed_use": "Required to anchor CTC/cfDNA signals to resected tumor clone.",
        "not_allowed": "Do not assume all low-input CTC NGS succeeds.",
    },
    {
        "dataset_or_source": "Postoperative follow-up",
        "data_type": "Tg/US/pathology/recurrence-risk clinical fields",
        "used_in_figures": "Future Aim 3; data unlock map; reviewer defense",
        "used_in_tables": "hospital_data_request_crf; statistical_analysis_plan",
        "allowed_use": "Secondary association with residual-risk features after adequate follow-up.",
        "not_allowed": "Immediate recurrence prediction claim before validation.",
    },
]


def draw_box(img, x1, y1, x2, y2, color, text_lines, header=None):
    import cv2

    cv2.rectangle(img, (x1, y1), (x2, y2), color, -1)
    cv2.rectangle(img, (x1, y1), (x2, y2), (25, 35, 48), 3)
    y = y1 + 42
    if header:
        cv2.putText(img, header, (x1 + 22, y), cv2.FONT_HERSHEY_SIMPLEX, 0.85, (143, 45, 37), 2, cv2.LINE_AA)
        y += 38
    for line in text_lines:
        cv2.putText(img, line, (x1 + 22, y), cv2.FONT_HERSHEY_SIMPLEX, 0.62, (16, 32, 51), 2, cv2.LINE_AA)
        y += 32


def draw_arrow(img, p1, p2, color=(143, 45, 37), thickness=4):
    import cv2

    cv2.arrowedLine(img, p1, p2, color, thickness, tipLength=0.04, line_type=cv2.LINE_AA)


def build_overview_png(path: Path) -> None:
    import cv2

    w, h = 2600, 1600
    img = np.full((h, w, 3), (232, 242, 248), dtype=np.uint8)
    # BGR palette
    cream = (241, 250, 255)
    good = (232, 243, 233)
    warn = (212, 241, 255)
    blue = (235, 226, 214)
    redlite = (219, 222, 247)
    coal = (47, 33, 23)
    red = (37, 45, 143)
    green = (80, 107, 66)

    cv2.putText(img, "PTC CTC-EMT-NGS Grand Overview", (70, 95), cv2.FONT_HERSHEY_SIMPLEX, 1.65, coal, 4, cv2.LINE_AA)
    cv2.putText(img, "Public priors -> Representative figures -> Control tables -> Prospective thyroidectomy perturbation test", (75, 145), cv2.FONT_HERSHEY_SIMPLEX, 0.72, red, 2, cv2.LINE_AA)

    draw_box(img, 80, 230, 520, 500, cream, ["PTC diagnosis", "Thyroidectomy", "Post-op surveillance", "Residual-risk uncertainty"], "Clinical problem")
    draw_box(img, 80, 620, 520, 930, warn, ["T0 pre-op blood", "T1 2 weeks", "T2 3 months", "T3 6-12 months", "CTC/cfDNA + tissue NGS"], "Future cohort")

    draw_box(img, 660, 220, 1110, 430, good, ["TCGA-THCA", "BRAF state split", "driver variance boundary"], "Public prior 1")
    draw_box(img, 660, 480, 1110, 700, good, ["GSE250521 spatial", "16 slides / 57,144 spots", "tissue organization"], "Public prior 2")
    draw_box(img, 660, 750, 1110, 980, good, ["scRNA marker context", "bulk proteomics proxy", "marker module support"], "Public prior 3")
    draw_box(img, 660, 1030, 1110, 1250, redlite, ["Yu 2024: feasibility", "62 PTC patients", "serial CTC-EMT around surgery"], "Published anchor")

    draw_box(img, 1260, 210, 1720, 410, cream, ["Fig 1 provenance", "Fig 2 mutation not enough", "Fig 3 spatial context"], "Representative figs")
    draw_box(img, 1260, 470, 1720, 700, cream, ["Fig 4 marker modules", "Fig 5 phenotype-to-NGS", "Fig 6 missing dataset"], "Representative figs")
    draw_box(img, 1260, 760, 1720, 990, blue, ["Data dictionary", "Caption scaffold", "Claim boundary", "CRF / SAP / kill criteria"], "Control tables")

    draw_box(img, 1880, 250, 2420, 520, good, ["Safe current paper:", "public-data prior/framework", "tissue-state + marker + NGS", "interpretation scaffold"], "Write now")
    draw_box(img, 1880, 650, 2420, 940, warn, ["Prospective breakthrough:", "Does thyroidectomy perturb", "CTC E -> E/M -> M state?", "Can persistent signal be", "genetically anchored?"], "Test next")
    draw_box(img, 1880, 1080, 2420, 1340, redlite, ["Do not claim yet:", "recurrence prediction", "ready clinical diagnostic", "public pilot proves CTC", "vendor tech is novelty"], "Boundary")

    for y in [350, 760]:
        draw_arrow(img, (520, y), (650, y))
    draw_arrow(img, (1110, 320), (1260, 310))
    draw_arrow(img, (1110, 590), (1260, 575))
    draw_arrow(img, (1110, 870), (1260, 575))
    draw_arrow(img, (1110, 1140), (1260, 310))
    draw_arrow(img, (1720, 310), (1875, 380))
    draw_arrow(img, (1720, 600), (1875, 790))
    draw_arrow(img, (1720, 890), (1875, 1190))
    draw_arrow(img, (520, 780), (1875, 790), color=green, thickness=5)

    cv2.putText(img, "Core rule: current public data builds the map; Yu/hospital serial blood tests the state transition.", (130, 1500), cv2.FONT_HERSHEY_SIMPLEX, 0.78, coal, 2, cv2.LINE_AA)
    cv2.putText(img, "Category: Science, not Technology/ICT. AI/ICT remains analysis infrastructure.", (130, 1540), cv2.FONT_HERSHEY_SIMPLEX, 0.72, red, 2, cv2.LINE_AA)
    cv2.imwrite(str(path), img)


def write_tables() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    fig_df = pd.DataFrame(FIGURE_ROWS)
    table_df = pd.DataFrame(TABLE_ROWS, columns=["table_file", "table_title", "role", "where_used", "claim_boundary"])
    data_df = pd.DataFrame(DATA_ROWS)
    for name, df in [(FIG_TABLE_NAME, fig_df), (TAB_TABLE_NAME, table_df), (DATA_MAP_NAME, data_df)]:
        out = atlas.TAB / name
        hub = atlas.HUB / name
        live = atlas.LIVE / name
        df.to_csv(out, sep="\t", index=False)
        df.to_csv(hub, sep="\t", index=False)
        shutil.copy2(hub, live)
    return fig_df, table_df, data_df


def fig_img(name: str) -> str:
    return f"{ASSET}/{esc(name)}"


def first_asset(asset_text: str) -> str:
    return asset_text.split(";")[0].strip()


def figure_card(row: dict) -> str:
    img = first_asset(row["representative_asset"])
    all_assets = [x.strip() for x in row["representative_asset"].split(";") if x.strip()]
    links = " ".join(f"<a href='{fig_img(x)}'>{esc(x)}</a>" for x in all_assets)
    return f"""
    <article class="fig-card">
      <a class="thumb" href="{fig_img(img)}"><img src="{fig_img(img)}" alt="{esc(row['figure_title'])}" loading="lazy" /></a>
      <div class="fig-body">
        <span>{esc(row['slot'])}</span>
        <h3>{esc(row['figure_title'])}</h3>
        <p>{esc(row['role'])}</p>
        <dl>
          <dt>Data used</dt><dd>{esc(row['data_used'])}</dd>
          <dt>Where used</dt><dd>{esc(row['where_used'])}</dd>
          <dt>Claim</dt><dd>{esc(row['claim'])}</dd>
          <dt>Boundary</dt><dd>{esc(row['claim_boundary'])}</dd>
          <dt>Assets</dt><dd class="asset-links">{links}</dd>
        </dl>
      </div>
    </article>
    """


def table_html(df: pd.DataFrame) -> str:
    cols = list(df.columns)
    rows = ["<div class='table-wrap'><table><thead><tr>"]
    rows += [f"<th>{esc(c)}</th>" for c in cols]
    rows.append("</tr></thead><tbody>")
    for _, r in df.iterrows():
        rows.append("<tr>")
        for c in cols:
            rows.append(f"<td>{esc(r[c])}</td>")
        rows.append("</tr>")
    rows.append("</tbody></table></div>")
    return "".join(rows)


def html_page(fig_df: pd.DataFrame, table_df: pd.DataFrame, data_df: pd.DataFrame) -> str:
    return f"""<!DOCTYPE html>
<html lang="ko">
<head>
<meta charset="UTF-8" />
<meta name="viewport" content="width=device-width, initial-scale=1.0" />
<title>K-Thyro CTC-EMT Representative Overview</title>
<link rel="preconnect" href="https://fonts.googleapis.com" />
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin />
<link href="https://fonts.googleapis.com/css2?family=Cormorant+Garamond:wght@500;600;700&family=Newsreader:ital,wght@0,400;0,500;0,600;1,400&family=JetBrains+Mono:wght@400;500;700&family=Noto+Sans+KR:wght@400;500;700;800&display=swap" rel="stylesheet" />
<style>
  :root{{--ink:#102033;--muted:#52606f;--paper:#f6f0e6;--card:#fffaf1;--line:#d4c6b3;--line2:#eadfcc;--red:#8f2d25;--blue:#244e73;--green:#426b50;--gold:#b58534;--coal:#17212f;--good:#e9f3e8;--bad:#f7dedb}}
  *{{box-sizing:border-box}} html{{scroll-behavior:smooth}}
  body{{margin:0;background:linear-gradient(135deg,#f8f2e8 0%,#efe3d1 58%,#f7efe0 100%);color:var(--ink);font-family:"Newsreader","Noto Sans KR",serif;line-height:1.6;font-size:17px}}
  a{{color:var(--red);text-decoration:none;border-bottom:1px solid rgba(143,45,37,.25)}} .wrap{{max-width:1250px;margin:0 auto;padding:32px 26px 84px}}
  .top{{display:flex;gap:10px;flex-wrap:wrap;margin-bottom:22px;font-family:"JetBrains Mono",monospace;font-size:12px}} .top a{{border:1px solid var(--line);background:rgba(255,255,255,.72);border-radius:999px;padding:8px 11px;color:var(--ink)}}
  .hero{{border:2px solid var(--coal);background:linear-gradient(145deg,#fffaf1 0%,#f0dfc4 100%);padding:44px 46px;box-shadow:10px 10px 0 rgba(23,33,47,.16)}}
  .kicker{{font-family:"JetBrains Mono",monospace;letter-spacing:.17em;text-transform:uppercase;font-size:12px;color:var(--red);font-weight:900;margin-bottom:14px}}
  h1{{font-family:"Cormorant Garamond",serif;font-size:60px;line-height:1.0;margin:0 0 14px;color:var(--coal);max-width:1040px}} .subtitle{{font-size:23px;color:var(--muted);max-width:1020px;font-style:italic}}
  .layout{{display:grid;grid-template-columns:250px 1fr;gap:26px;margin-top:30px}} .toc{{position:sticky;top:18px;align-self:start;background:rgba(255,250,241,.92);border:1px solid var(--line);border-radius:12px;padding:16px}} .toc h3{{font-family:"Cormorant Garamond",serif;color:var(--red);font-size:24px;margin:0 0 10px}} .toc a{{display:block;color:var(--ink);border-bottom:1px dotted var(--line);padding:8px 0;font-size:14px}}
  section{{background:rgba(255,250,241,.9);border:1px solid var(--line);border-radius:14px;padding:28px;margin-bottom:24px}} h2{{font-family:"Cormorant Garamond",serif;font-size:42px;line-height:1.05;margin:0 0 12px;color:var(--coal)}} .section-note{{color:var(--muted);font-style:italic;margin:0 0 18px}}
  .overview-fig{{display:block;background:#fff;border:1px solid var(--line);border-radius:14px;padding:14px}} .overview-fig img{{display:block;width:100%;height:auto;border-radius:8px;background:#f4eadb}}
  .rule-grid{{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));gap:10px}} .rule-grid div{{background:#fff;border:1px solid var(--line2);border-radius:12px;padding:14px}} .rule-grid b{{font-family:"JetBrains Mono",monospace;font-size:10px;color:var(--red);text-transform:uppercase;letter-spacing:.08em}} .rule-grid p{{margin:8px 0 0}}
  .fig-card{{display:grid;grid-template-columns:340px 1fr;gap:16px;background:#fff;border:1px solid var(--line2);border-radius:14px;padding:14px;margin:16px 0}} .thumb{{border:0}} .thumb img{{display:block;width:100%;height:auto;border-radius:8px;background:#f4eadb}} .fig-body span{{font-family:"JetBrains Mono",monospace;font-size:10px;background:var(--coal);color:#fff;border-radius:999px;padding:5px 8px}} .fig-body h3{{font-family:"Cormorant Garamond",serif;font-size:31px;line-height:1.05;margin:13px 0 8px;color:var(--coal)}} .fig-body p{{margin:0 0 10px;color:var(--muted)}} dl{{display:grid;grid-template-columns:120px 1fr;gap:6px 10px;margin:0}} dt{{font-family:"JetBrains Mono",monospace;font-size:10px;text-transform:uppercase;letter-spacing:.06em;color:var(--red);padding-top:2px}} dd{{margin:0;color:var(--ink)}} .asset-links a{{display:inline-block;margin:0 6px 4px 0;font-size:12px}}
  .table-wrap{{overflow-x:auto}} table{{width:100%;border-collapse:collapse;background:#fff;border:1px solid var(--line);font-size:14px}} th{{background:var(--coal);color:#fff;text-align:left;font-family:"JetBrains Mono",monospace;font-size:11px;letter-spacing:.06em;text-transform:uppercase;padding:10px}} td{{border-top:1px solid var(--line2);padding:10px;vertical-align:top}}
  .link-grid{{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));gap:12px}} .link-card{{display:block;background:#fff;border:1px solid var(--line2);border-radius:12px;padding:16px}} .link-card b{{display:block;font-family:"Cormorant Garamond",serif;font-size:27px;color:var(--red)}} .link-card span{{display:block;color:var(--muted)}}
  .footer{{font-family:"JetBrains Mono",monospace;font-size:12px;color:var(--muted);text-align:center;margin-top:34px;padding-top:18px;border-top:1px solid var(--line)}}
  @media(max-width:1000px){{.layout{{grid-template-columns:1fr}}.toc{{position:relative;top:auto}}.fig-card{{grid-template-columns:1fr}}.rule-grid,.link-grid{{grid-template-columns:1fr}}h1{{font-size:42px}}.hero{{padding:30px 22px}}}}
</style>
</head>
<body>
<div class="wrap">
  <nav class="top">
    <a href="kthyro_ctc_emt_all_pages_launchpad.html">All Pages</a>
    <a href="kthyro_ctc_emt_breakthrough_case.html">Breakthrough Case</a>
    <a href="kthyro_ctc_emt_connected_reader.html">Connected Reader</a>
    <a href="kthyro_ctc_emt_manuscript_storyboard_paper1_style.html">Storyboard</a>
    <a href="kthyro_ctc_emt_visual_atlas_evidence.html">Evidence</a>
  </nav>
  <header class="hero">
    <div class="kicker">K-Thyro CTC-EMT · representative figure/table/data overview</div>
    <h1>Representative Overview</h1>
    <p class="subtitle">대표 figure, 대표 table, 전체 overview 그림, figure 흐름 설명, 각 데이터가 어디에 쓰였는지와 어디까지 주장 가능한지를 한 페이지에 고정한다.</p>
  </header>
  <div class="layout">
    <aside class="toc">
      <h3>Read Order</h3>
      <a href="#overview">Grand Overview</a>
      <a href="#rules">Core Rules</a>
      <a href="#figures">Representative Figures</a>
      <a href="#tables">Representative Tables</a>
      <a href="#data-map">Data Used / Where Used</a>
      <a href="#downloads">Downloads</a>
    </aside>
    <main>
      <section id="overview">
        <h2>Grand Overview Figure</h2>
        <p class="section-note">이 한 장이 전체 흐름이다: public priors가 representative figures로 들어가고, control tables가 claim boundary를 잡고, hospital/Yu data가 prospective breakthrough를 unlock한다.</p>
        <a class="overview-fig" href="{fig_img(OVERVIEW_FIG)}"><img src="{fig_img(OVERVIEW_FIG)}" alt="Grand overview" /></a>
      </section>
      <section id="rules">
        <h2>Core Rules</h2>
        <div class="rule-grid">
          <div><b>Write now</b><p>public-data prior/framework: mutation-state heterogeneity, spatial context, marker module support, NGS anchoring logic.</p></div>
          <div><b>Test next</b><p>thyroidectomy-induced CTC-EMT state transition with serial hospital/Yu blood and matched NGS.</p></div>
          <div><b>Do not claim</b><p>recurrence prediction, clinical diagnostic readiness, public pilot proves CTC biology.</p></div>
          <div><b>Category</b><p>Science-track perturbation biology. AI/ICT is analysis infrastructure only.</p></div>
        </div>
      </section>
      <section id="figures">
        <h2>Representative Figure Spine</h2>
        <p class="section-note">각 figure가 어떤 데이터에서 왔고, 어디에 쓰이고, 무엇을 말하며, 어디서 멈춰야 하는지.</p>
        {''.join(figure_card(r.to_dict()) for _, r in fig_df.iterrows())}
      </section>
      <section id="tables">
        <h2>Representative Table Spine</h2>
        <p class="section-note">글쓰기와 심사 방어에서 어떤 표를 어디에 쓰는지.</p>
        {table_html(table_df)}
      </section>
      <section id="data-map">
        <h2>Data Used / Where Used / Claim Boundary</h2>
        <p class="section-note">리뷰어가 “이 데이터 어디에 썼냐”라고 물을 때 바로 답하는 표.</p>
        {table_html(data_df)}
      </section>
      <section id="downloads">
        <h2>Downloads</h2>
        <div class="link-grid">
          <a class="link-card" href="{FIG_TABLE_NAME}"><b>Figure TSV</b><span>대표 figure overview</span></a>
          <a class="link-card" href="{TAB_TABLE_NAME}"><b>Table TSV</b><span>대표 table overview</span></a>
          <a class="link-card" href="{DATA_MAP_NAME}"><b>Data Map TSV</b><span>data-to-figure-claim map</span></a>
          <a class="link-card" href="{REPORT_NAME}"><b>Guide MD</b><span>보고서 버전</span></a>
        </div>
      </section>
    </main>
  </div>
  <div class="footer">Representative overview: ambitious science story, bounded claims, explicit data provenance.</div>
</div>
</body>
</html>"""


def report_md(fig_df: pd.DataFrame, table_df: pd.DataFrame, data_df: pd.DataFrame) -> str:
    lines = [
        "# Representative Figure/Table/Data Overview",
        "",
        "Purpose: 대표 figure, 대표 table, 전체 overview 그림, figure 흐름 설명, 데이터 사용 위치, claim boundary를 한 문서에 고정한다.",
        "",
        "Core rule: public data는 prior/framework이고, thyroidectomy-induced CTC-EMT transition은 Yu/hospital serial blood + matched NGS로 prospective test한다.",
        "",
        "## Representative Figures",
    ]
    for _, r in fig_df.iterrows():
        lines += [
            f"### {r['slot']}. {r['figure_title']}",
            f"- Asset: {r['representative_asset']}",
            f"- Role: {r['role']}",
            f"- Data used: {r['data_used']}",
            f"- Where used: {r['where_used']}",
            f"- Claim: {r['claim']}",
            f"- Boundary: {r['claim_boundary']}",
            "",
        ]
    lines += ["## Representative Tables"]
    for _, r in table_df.iterrows():
        lines.append(f"- {r['table_title']} ({r['table_file']}): {r['role']} / Where used: {r['where_used']} / Boundary: {r['claim_boundary']}")
    lines += ["", "## Data Used / Where Used"]
    for _, r in data_df.iterrows():
        lines.append(f"- {r['dataset_or_source']}: used in {r['used_in_figures']}; allowed: {r['allowed_use']}; not allowed: {r['not_allowed']}")
    return "\n".join(lines)


def update_index_handoff_launchpad() -> None:
    idx = atlas.HUB / "index.html"
    if idx.exists():
        text = idx.read_text(encoding="utf-8")
        if HTML_NAME not in text:
            needle = '    <a href="kthyro_ctc_emt_breakthrough_case.html"'
            insert = f'    <a href="{HTML_NAME}" style="background:linear-gradient(135deg,#fffaf1 0%,#e9f3e8 100%);border:3px solid #426b50;color:#102033;font-weight:1000;font-size:18px;padding:16px 20px">★★ CTC-EMT REPRESENTATIVE OVERVIEW · figures tables data map</a>\n'
            pos = text.find(needle)
            text = text[:pos] + insert + text[pos:] if pos >= 0 else text + "\n" + insert
            idx.write_text(text, encoding="utf-8")
        shutil.copy2(idx, atlas.LIVE / "index.html")

    for path in [atlas.REP / "VISUAL_WEB_HANDOFF_KR.md", atlas.HUB / "VISUAL_WEB_HANDOFF_KR.md"]:
        if not path.exists():
            continue
        text = path.read_text(encoding="utf-8")
        line = "- Representative Overview: http://40.82.129.113:8012/papers_hub_2026_05_04/kthyro_ctc_emt_representative_overview.html\n"
        if HTML_NAME not in text:
            text = text.replace("## 최우선 링크\n\n", "## 최우선 링크\n\n" + line)
            path.write_text(text, encoding="utf-8")
    handoff = atlas.HUB / "VISUAL_WEB_HANDOFF_KR.md"
    if handoff.exists():
        shutil.copy2(handoff, atlas.LIVE / "VISUAL_WEB_HANDOFF_KR.md")

    launch = atlas.HUB / "kthyro_ctc_emt_all_pages_launchpad.html"
    if launch.exists() and HTML_NAME not in launch.read_text(encoding="utf-8"):
        text = launch.read_text(encoding="utf-8")
        insert = f'<a href="{HTML_NAME}"><b>대표 Overview</b><span>Figure/Table/Data Map</span></a>'
        text = text.replace('<div class="route-grid">', '<div class="route-grid">' + insert, 1)
        launch.write_text(text, encoding="utf-8")
        shutil.copy2(launch, atlas.LIVE / "kthyro_ctc_emt_all_pages_launchpad.html")


def main() -> None:
    atlas.FIG = atlas.OUT / "figures" if not hasattr(atlas, "FIG") else atlas.FIG
    fig_dir = atlas.OUT / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    asset_dir = atlas.HUB / "assets" / atlas.ASSET
    live_asset_dir = atlas.LIVE / "assets" / atlas.ASSET
    overview_path = fig_dir / OVERVIEW_FIG
    build_overview_png(overview_path)
    shutil.copy2(overview_path, asset_dir / OVERVIEW_FIG)
    shutil.copy2(overview_path, live_asset_dir / OVERVIEW_FIG)

    fig_df, table_df, data_df = write_tables()
    html_text = html_page(fig_df, table_df, data_df)
    src = atlas.HUB / HTML_NAME
    live = atlas.LIVE / HTML_NAME
    local = atlas.REP / HTML_NAME
    src.write_text(html_text, encoding="utf-8")
    local.write_text(html_text, encoding="utf-8")
    shutil.copy2(src, live)

    report_text = report_md(fig_df, table_df, data_df)
    report = atlas.REP / REPORT_NAME
    hub_report = atlas.HUB / REPORT_NAME
    report.write_text(report_text, encoding="utf-8")
    hub_report.write_text(report_text, encoding="utf-8")
    shutil.copy2(hub_report, atlas.LIVE / REPORT_NAME)

    update_index_handoff_launchpad()
    print(f"HTML={src}")
    print(f"LIVE={live}")
    print(f"FIGURE={overview_path}")
    print(f"REPORT={report}")
    print(f"TABLES={atlas.TAB / FIG_TABLE_NAME}, {atlas.TAB / TAB_TABLE_NAME}, {atlas.TAB / DATA_MAP_NAME}")


if __name__ == "__main__":
    main()
