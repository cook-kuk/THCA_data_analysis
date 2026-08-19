#!/usr/bin/env python3
from __future__ import annotations

import html
import shutil
from pathlib import Path

import numpy as np
import pandas as pd

import build_paper1_style_visual_atlas as atlas
import build_representative_overview_pack as overview


HTML_NAME = "kthyro_ctc_emt_paper_ready_figure_table_pack.html"
REPORT_NAME = "PAPER_READY_FIGURE_TABLE_PACK_KR.md"
FIG_LEGENDS_TSV = "paper_ready_main_figure_legends.tsv"
TABLE_LEGENDS_TSV = "paper_ready_table_legends.tsv"
SECTION_MAP_TSV = "manuscript_section_to_asset_map.tsv"
FLOW_FIG = "PF16_paper_ready_assembly_flow.png"
ASSET = f"assets/{atlas.ASSET}"


def esc(x) -> str:
    return html.escape("" if x is None else str(x), quote=True)


SECTION_MAP = [
    {
        "manuscript_section": "Title / Abstract",
        "main_message": "A public-data prior and prospective Science proposal for serial PTC CTC-EMT-NGS mapping after thyroidectomy.",
        "primary_assets": "PF15_grand_overview_data_to_claim_map.png;PF11_manuscript_blueprint_board.png",
        "tables": "data_dictionary_master.tsv;claim_boundary_matrix.tsv",
        "write_now_sentence": "We define a staged framework for testing whether thyroidectomy perturbs CTC-EMT state and whether persistent postoperative genetic signals can be anchored to the resected tumor clone.",
        "do_not_write": "Do not say the postoperative CTC-EMT transition has already been discovered in our cohort.",
    },
    {
        "manuscript_section": "Introduction",
        "main_message": "PTC is usually curable, but postoperative residual-risk biology is not directly read in real time.",
        "primary_assets": "PF01_manuscript_flow.png;slide01_problem.png;slide02_hypothesis.png",
        "tables": "manuscript_flow_outline.tsv;safe_sentence_bank.tsv",
        "write_now_sentence": "Thyroidectomy provides a biologically interpretable perturbation around which serial blood and matched tissue genotyping can be organized.",
        "do_not_write": "Do not frame this as a broad AI/spatial/drug platform.",
    },
    {
        "manuscript_section": "Results 1",
        "main_message": "Driver mutation alone is insufficient to define PTC biological state.",
        "primary_assets": "F03_tcga_braf_state_split.png;F04_tcga_driver_variance_boundary.png",
        "tables": "data_to_figure_claim_map.tsv;key_number_ledger.tsv",
        "write_now_sentence": "BRAF-mutant PTC tumors distributed across multiple state labels, supporting the need for phenotype-state interpretation beyond driver genotype.",
        "do_not_write": "Do not infer CTC shedding or recurrence from TCGA tissue-state heterogeneity.",
    },
    {
        "manuscript_section": "Results 2",
        "main_message": "Public spatial data support organized tissue-state context.",
        "primary_assets": "F05_spatial_tissue_state_organization.png;spatial_coherence_summary.png",
        "tables": "representative_figure_overview.tsv;data_dictionary_master.tsv",
        "write_now_sentence": "Spatial tissue-state organization supports marker and ROI logic for prospective CTC-EMT studies.",
        "do_not_write": "Do not claim spatial hotspots are direct CTC shedding sources.",
    },
    {
        "manuscript_section": "Results 3",
        "main_message": "Marker modules have cross-layer tissue/cell/protein support.",
        "primary_assets": "F06_scrna_marker_context.png;F07_proteomics_dediff_direction.png",
        "tables": "marker_module_cross_layer_support_matrix.tsv;representative_table_overview.tsv",
        "write_now_sentence": "Epithelial, hybrid/mesenchymal, thyroid-lineage, stress/survival, and immune/stress modules are interpretable across public evidence layers. A separate supplementary planning matrix (PF10) records reviewer-visible manuscript-planning judgement on each layer and is not used as data here.",
        "do_not_write": "Do not present the marker panel as a validated clinical CTC assay. Do not cite the PF10 supplementary planning matrix as a cross-layer measurement; it is hand-coded judgement only.",
    },
    {
        "manuscript_section": "Results 4",
        "main_message": "Phenotype must be genetically anchored by matched tumor-normal NGS and cfDNA/CTC-enriched sequencing.",
        "primary_assets": "PF03_marker_to_ngs_map.png;F08_ctc_emt_ngs_prior_panel.png",
        "tables": "ctc_emt_marker_prior_table.tsv;statistical_analysis_plan.tsv",
        "write_now_sentence": "The proposed CTC-EMT readout is only interpretable as residual-risk biology when anchored to tissue-derived genetic context.",
        "do_not_write": "Do not assume low-input CTC sequencing succeeds in every early PTC case.",
    },
    {
        "manuscript_section": "Results 5 / Discussion",
        "main_message": "The exact public dataset needed for serial postoperative CTC-EMT genetic mapping is missing.",
        "primary_assets": "PF04_dataset_gap_bridge.png;PF12_future_data_unlock_map.png;PF05_claim_boundary_matrix.png",
        "tables": "hospital_data_request_crf.tsv;stage_gate_kill_criteria.tsv",
        "write_now_sentence": "The missing joint dataset defines the prospective Samsung/Yu cohort requirement and its stage-gated endpoints.",
        "do_not_write": "Do not use dataset absence as a substitute for prospective biological evidence.",
    },
    {
        "manuscript_section": "Supplement / Reviewer Defense",
        "main_message": "Reviewer attack, claim boundary, and kill criteria are explicitly pre-specified.",
        "primary_assets": "PF07_reviewer_to_figure_map.png;PF14_analysis_control_board.png;F09_publishability_decision.png",
        "tables": "reviewer_risk_heatmap.tsv;safe_sentence_bank.tsv;stage_gate_kill_criteria.tsv",
        "write_now_sentence": "The package is designed to be ambitious but falsifiable, with downgrade criteria for assay or data failure.",
        "do_not_write": "Do not hide feasibility risks or overstate current public-data claims.",
    },
    {
        "manuscript_section": "Supplementary Planning Matrix",
        "main_message": "Supplementary, reviewer-visible hand-coded planning judgement on each public/published evidence layer for the marker modules (not a data analysis).",
        "primary_assets": "PF10_marker_module_support_heatmap.png",
        "tables": "marker_module_cross_layer_support_matrix.tsv",
        "write_now_sentence": "We hand-coded a 0/1/2 manuscript-planning judgement of each evidence layer's usefulness for prior construction. The TCGA column records judgement about TCGA-THCA tumor tissue only (TCGA has no CTC, blood, or liquid biopsy data). The Yu 2024 column scores the published Yu 2024 paper's reported feasibility summary; no unpublished or raw Yu CTC data was used.",
        "do_not_write": "Do not cite this matrix as a cross-layer measurement, an assay-performance metric, or evidence of CTC biology. It is a manuscript-planning judgement only.",
    },
]


def build_figure_legends() -> pd.DataFrame:
    rows = []
    for r in overview.FIGURE_ROWS:
        slot = r["slot"]
        title = r["figure_title"]
        if slot == "Grand Overview":
            legend = (
                "Grand overview of the K-Thyro CTC-EMT-NGS proposal logic. "
                "Public tissue, spatial, single-cell, proteomic, and published CTC evidence are separated from the future prospective serial blood experiment. "
                "The figure maps how public priors support representative figures and control tables, while Yu/hospital serial CTC, cfDNA, and matched tissue-normal NGS data are required to test thyroidectomy-induced CTC-EMT state transition."
            )
            result = "The overview fixes the current manuscript as a public-data prior/framework and the Samsung proposal as a prospective perturbation test."
        elif slot == "Main Fig. 1":
            legend = (
                "Study logic and evidence provenance. Yu 2024 is used as a published feasibility anchor for serial CTC-EMT measurement in PTC, whereas TCGA, spatial, scRNA, and proteomic resources are used as public marker and tissue-state priors."
            )
            result = "Evidence sources were assigned distinct functions so that public omics are not treated as substitutes for serial CTC data."
        elif slot == "Main Fig. 2":
            legend = (
                "Driver mutation alone does not define a single biological state in PTC. "
                "TCGA-THCA public-pilot outputs show that BRAF-mutant tumors distribute across multiple state/vulnerability labels, with the largest BRAF label accounting for 33.2% of the BRAF-mutant subset."
            )
            result = "This supports phenotype-state plus genotype interpretation rather than driver-only postoperative blood-signal interpretation."
        elif slot == "Main Fig. 3":
            legend = (
                "Spatial tissue-state organization in public thyroid spatial transcriptomics. "
                "GSE250521-derived coherence and spatial summary outputs support non-random organization of tissue states across 16 slides and 57,144 spots."
            )
            result = "Spatial organization supports marker/ROI logic for prospective CTC-EMT studies, not direct CTC-origin proof."
        elif slot == "Main Fig. 4":
            legend = (
                "Cross-layer support for candidate CTC-EMT marker modules. "
                "Single-cell marker context, proteomic direction, and marker-module matrices support epithelial, hybrid/mesenchymal, thyroid-lineage, stress/survival, and immune/stress modules."
            )
            result = "The marker panel is biologically organized and suitable for prospective testing but remains hypothesis-generating."
        elif slot == "Main Fig. 5":
            legend = (
                "Phenotype-to-genotype bridge for CTC-EMT-NGS interpretation. "
                "CTC epithelial, hybrid E/M, and mesenchymal states are mapped to matched tissue-normal NGS, cfDNA, and CTC-enriched sequencing logic."
            )
            result = "The proposed blood readout requires genetic anchoring because phenotype alone can be nonspecific in early PTC."
        elif slot == "Main Fig. 6":
            legend = (
                "Missing dataset and prospective data unlock map. "
                "No public dataset jointly provides serial pre/post-thyroidectomy CTC-EMT phenotype, matched tumor-normal NGS, cfDNA or CTC-enriched NGS, and postoperative follow-up."
            )
            result = "The data gap justifies the prospective Samsung/Yu cohort and defines required CRF and assay fields."
        else:
            legend = (
                "Reviewer defense and execution-control figure set. "
                "Expected reviewer attacks, analysis-control tables, publishability boundaries, and stage-gated kill criteria are mapped to supporting figures and tables."
            )
            result = "The defense layer keeps the proposal ambitious but bounded and falsifiable."
        rows.append(
            {
                "slot": slot,
                "figure_title": title,
                "assets": r["representative_asset"],
                "paper_ready_legend": legend,
                "results_scaffold": result,
                "data_used": r["data_used"],
                "where_used": r["where_used"],
                "claim_boundary": r["claim_boundary"],
            }
        )
    return pd.DataFrame(rows)


def build_table_legends() -> pd.DataFrame:
    rows = []
    for table_file, table_title, role, where_used, boundary in overview.TABLE_ROWS:
        rows.append(
            {
                "table_file": table_file,
                "table_title": table_title,
                "paper_ready_caption": f"{table_title}. {role}. This table is used for {where_used}.",
                "where_used": where_used,
                "claim_boundary": boundary,
            }
        )
    return pd.DataFrame(rows)


def draw_flow_png(path: Path) -> None:
    import cv2

    w, h = 2600, 1500
    img = np.full((h, w, 3), (238, 244, 248), dtype=np.uint8)
    coal = (47, 33, 23)
    red = (37, 45, 143)
    green = (80, 107, 66)
    cream = (241, 250, 255)
    good = (232, 243, 233)
    warn = (212, 241, 255)
    bad = (219, 222, 247)

    cv2.putText(img, "Paper-Ready Assembly Flow", (70, 95), cv2.FONT_HERSHEY_SIMPLEX, 1.75, coal, 4, cv2.LINE_AA)
    cv2.putText(img, "Representative overview -> figure legends -> table captions -> manuscript sections -> reviewer defense", (75, 145), cv2.FONT_HERSHEY_SIMPLEX, 0.74, red, 2, cv2.LINE_AA)

    def box(x1, y1, x2, y2, color, head, lines):
        cv2.rectangle(img, (x1, y1), (x2, y2), color, -1)
        cv2.rectangle(img, (x1, y1), (x2, y2), coal, 3)
        cv2.putText(img, head, (x1 + 24, y1 + 48), cv2.FONT_HERSHEY_SIMPLEX, 0.82, red, 2, cv2.LINE_AA)
        y = y1 + 92
        for line in lines:
            cv2.putText(img, line, (x1 + 24, y), cv2.FONT_HERSHEY_SIMPLEX, 0.62, coal, 2, cv2.LINE_AA)
            y += 34

    def arrow(x1, y1, x2, y2, color=red):
        cv2.arrowedLine(img, (x1, y1), (x2, y2), color, 4, tipLength=0.04, line_type=cv2.LINE_AA)

    box(80, 230, 560, 480, cream, "1. Grand overview", ["PF15 data-to-claim map", "one-slide storyline", "Science category lock"])
    box(700, 230, 1180, 480, good, "2. Figure spine", ["Main Fig 1-6", "Defense/Supp figures", "legend + results scaffold"])
    box(1320, 230, 1800, 480, warn, "3. Table spine", ["data dictionary", "claim boundary", "CRF / SAP / kill criteria"])
    box(1940, 230, 2420, 480, good, "4. Manuscript map", ["section-by-section use", "where each asset goes", "safe sentence per section"])

    box(360, 660, 860, 930, cream, "Current paper", ["public-data prior/framework", "mutation-state heterogeneity", "spatial/marker/NGS logic"])
    box(1040, 660, 1540, 930, good, "Samsung proposal", ["thyroidectomy perturbation", "serial CTC-EMT readout", "matched NGS genetic anchor"])
    box(1720, 660, 2220, 930, bad, "Do not claim", ["recurrence prediction now", "ready diagnostic", "public data proves CTC"])

    box(520, 1080, 2080, 1280, warn, "Final handoff", ["Open Representative Overview first, then Breakthrough Case, then Storyboard.", "Every figure/table has data used, where used, and claim boundary."])

    arrow(560, 355, 700, 355)
    arrow(1180, 355, 1320, 355)
    arrow(1800, 355, 1940, 355)
    arrow(560, 480, 610, 660)
    arrow(1180, 480, 1290, 660)
    arrow(1800, 480, 1970, 660)
    arrow(860, 930, 1000, 1080, green)
    arrow(1540, 930, 1440, 1080, green)
    arrow(2220, 930, 1900, 1080)

    cv2.putText(img, "Rule: A high-impact story is acceptable only when each claim remains inside its data layer.", (105, 1410), cv2.FONT_HERSHEY_SIMPLEX, 0.76, coal, 2, cv2.LINE_AA)
    cv2.imwrite(str(path), img)


def table_html(df: pd.DataFrame, cls: str = "") -> str:
    out = [f"<div class='table-wrap {cls}'><table><thead><tr>"]
    for c in df.columns:
        out.append(f"<th>{esc(c)}</th>")
    out.append("</tr></thead><tbody>")
    for _, r in df.iterrows():
        out.append("<tr>")
        for c in df.columns:
            out.append(f"<td>{esc(r[c])}</td>")
        out.append("</tr>")
    out.append("</tbody></table></div>")
    return "".join(out)


def asset_link(name: str) -> str:
    return f"{ASSET}/{esc(name)}"


def figure_blocks(fig_df: pd.DataFrame) -> str:
    blocks = []
    for _, r in fig_df.iterrows():
        img = str(r["assets"]).split(";")[0].strip()
        blocks.append(
            f"""
            <article class="legend-card">
              <a class="legend-img" href="{asset_link(img)}"><img src="{asset_link(img)}" alt="{esc(r['figure_title'])}" loading="lazy" /></a>
              <div class="legend-body">
                <span>{esc(r['slot'])}</span>
                <h3>{esc(r['figure_title'])}</h3>
                <h4>Paper-ready legend</h4>
                <p>{esc(r['paper_ready_legend'])}</p>
                <h4>Results scaffold</h4>
                <p>{esc(r['results_scaffold'])}</p>
                <div class="mini-grid">
                  <div><b>Data used</b>{esc(r['data_used'])}</div>
                  <div><b>Where used</b>{esc(r['where_used'])}</div>
                  <div class="limit"><b>Claim boundary</b>{esc(r['claim_boundary'])}</div>
                </div>
              </div>
            </article>
            """
        )
    return "\n".join(blocks)


def html_page(fig_df: pd.DataFrame, table_df: pd.DataFrame, section_df: pd.DataFrame) -> str:
    return f"""<!DOCTYPE html>
<html lang="ko">
<head>
<meta charset="UTF-8" />
<meta name="viewport" content="width=device-width, initial-scale=1.0" />
<title>K-Thyro CTC-EMT Paper-Ready Figure/Table Pack</title>
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
  .legend-card{{display:grid;grid-template-columns:330px 1fr;gap:16px;background:#fff;border:1px solid var(--line2);border-radius:14px;padding:14px;margin:16px 0}} .legend-img{{border:0}} .legend-img img{{display:block;width:100%;height:auto;border-radius:8px;background:#f4eadb}} .legend-body span{{font-family:"JetBrains Mono",monospace;font-size:10px;background:var(--coal);color:#fff;border-radius:999px;padding:5px 8px}} .legend-body h3{{font-family:"Cormorant Garamond",serif;font-size:31px;line-height:1.05;margin:13px 0 8px;color:var(--coal)}} .legend-body h4{{font-family:"JetBrains Mono",monospace;font-size:11px;letter-spacing:.08em;text-transform:uppercase;color:var(--red);margin:12px 0 4px}} .legend-body p{{margin:0 0 8px;color:var(--ink)}}
  .mini-grid{{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:8px;margin-top:10px}} .mini-grid div{{background:#fff8e8;border:1px solid var(--line2);border-radius:8px;padding:10px;font-size:14px}} .mini-grid b{{font-family:"JetBrains Mono",monospace;font-size:10px;color:var(--red);display:block;text-transform:uppercase;letter-spacing:.06em;margin-bottom:4px}} .mini-grid .limit{{background:#f8ead7}}
  .table-wrap{{overflow-x:auto}} table{{width:100%;border-collapse:collapse;background:#fff;border:1px solid var(--line);font-size:14px}} th{{background:var(--coal);color:#fff;text-align:left;font-family:"JetBrains Mono",monospace;font-size:11px;letter-spacing:.06em;text-transform:uppercase;padding:10px}} td{{border-top:1px solid var(--line2);padding:10px;vertical-align:top}}
  .link-grid{{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));gap:12px}} .link-card{{display:block;background:#fff;border:1px solid var(--line2);border-radius:12px;padding:16px}} .link-card b{{display:block;font-family:"Cormorant Garamond",serif;font-size:27px;color:var(--red)}} .link-card span{{display:block;color:var(--muted)}}
  .footer{{font-family:"JetBrains Mono",monospace;font-size:12px;color:var(--muted);text-align:center;margin-top:34px;padding-top:18px;border-top:1px solid var(--line)}}
  @media(max-width:1000px){{.layout{{grid-template-columns:1fr}}.toc{{position:relative;top:auto}}.legend-card{{grid-template-columns:1fr}}.rule-grid,.link-grid,.mini-grid{{grid-template-columns:1fr}}h1{{font-size:42px}}.hero{{padding:30px 22px}}}}
</style>
</head>
<body>
<div class="wrap">
  <nav class="top">
    <a href="kthyro_ctc_emt_all_pages_launchpad.html">All Pages</a>
    <a href="kthyro_ctc_emt_representative_overview.html">Representative Overview</a>
    <a href="kthyro_ctc_emt_breakthrough_case.html">Breakthrough Case</a>
    <a href="kthyro_ctc_emt_connected_reader.html">Connected Reader</a>
    <a href="kthyro_ctc_emt_manuscript_storyboard_paper1_style.html">Storyboard</a>
  </nav>
  <header class="hero">
    <div class="kicker">K-Thyro CTC-EMT · paper-ready figure/table/section pack</div>
    <h1>Paper-Ready Figure/Table Pack</h1>
    <p class="subtitle">대표 figure와 table을 실제 manuscript/proposal에 바로 붙일 수 있게 legend, results scaffold, section map, data used, where used, do-not-write boundary까지 완성한 패키지.</p>
  </header>
  <div class="layout">
    <aside class="toc">
      <h3>Read Order</h3>
      <a href="#flow">Assembly Flow</a>
      <a href="#rules">Use Rules</a>
      <a href="#figures">Figure Legends</a>
      <a href="#tables">Table Legends</a>
      <a href="#sections">Section Map</a>
      <a href="#downloads">Downloads</a>
    </aside>
    <main>
      <section id="flow">
        <h2>Assembly Flow Figure</h2>
        <p class="section-note">Representative Overview에서 원고/제안서로 들어가는 조립 흐름.</p>
        <a class="overview-fig" href="{asset_link(FLOW_FIG)}"><img src="{asset_link(FLOW_FIG)}" alt="Paper-ready assembly flow" /></a>
      </section>
      <section id="rules">
        <h2>Use Rules</h2>
        <div class="rule-grid">
          <div><b>First</b><p>Representative Overview로 전체 그림을 보여준다.</p></div>
          <div><b>Then</b><p>이 페이지의 figure legend와 section map을 원고/슬라이드에 배치한다.</p></div>
          <div><b>Always</b><p>Data used와 claim boundary를 같이 붙인다.</p></div>
          <div><b>Never</b><p>Yu/hospital data 없이 CTC transition discovery로 쓰지 않는다.</p></div>
        </div>
      </section>
      <section id="figures">
        <h2>Paper-Ready Figure Legends</h2>
        <p class="section-note">각 대표 figure의 legend, Results scaffold, 데이터, 사용 위치, boundary.</p>
        {figure_blocks(fig_df)}
      </section>
      <section id="tables">
        <h2>Paper-Ready Table Captions</h2>
        <p class="section-note">각 표가 manuscript/proposal 어디에 들어가고 어디까지 주장 가능한지.</p>
        {table_html(table_df)}
      </section>
      <section id="sections">
        <h2>Manuscript Section To Asset Map</h2>
        <p class="section-note">원고 섹션별 main message, asset, table, write-now sentence, do-not-write sentence.</p>
        {table_html(section_df)}
      </section>
      <section id="downloads">
        <h2>Downloads</h2>
        <div class="link-grid">
          <a class="link-card" href="{FIG_LEGENDS_TSV}"><b>Figure Legends TSV</b><span>paper-ready figure legends</span></a>
          <a class="link-card" href="{TABLE_LEGENDS_TSV}"><b>Table Captions TSV</b><span>paper-ready table captions</span></a>
          <a class="link-card" href="{SECTION_MAP_TSV}"><b>Section Map TSV</b><span>section-to-asset map</span></a>
          <a class="link-card" href="{REPORT_NAME}"><b>Guide MD</b><span>markdown handoff</span></a>
        </div>
      </section>
    </main>
  </div>
  <div class="footer">Paper-ready pack: stronger story, explicit data provenance, no unsupported CTC overclaim.</div>
</div>
</body>
</html>"""


def report_md(fig_df: pd.DataFrame, table_df: pd.DataFrame, section_df: pd.DataFrame) -> str:
    lines = [
        "# Paper-Ready Figure/Table Pack",
        "",
        "Purpose: 대표 figure/table을 실제 manuscript/proposal에 바로 배치할 수 있게 legend, results scaffold, section map, data used, where used, claim boundary를 완성한다.",
        "",
        "Core boundary: 현재 public-data output은 tissue-state/marker/NGS-anchor prior framework다. Yu/hospital serial CTC/cfDNA/tissue NGS data 없이 postoperative CTC-EMT transition discovery라고 쓰지 않는다.",
        "",
        "## Figure Legends",
    ]
    for _, r in fig_df.iterrows():
        lines += [
            f"### {r['slot']}. {r['figure_title']}",
            f"- Assets: {r['assets']}",
            f"- Legend: {r['paper_ready_legend']}",
            f"- Results scaffold: {r['results_scaffold']}",
            f"- Data used: {r['data_used']}",
            f"- Where used: {r['where_used']}",
            f"- Boundary: {r['claim_boundary']}",
            "",
        ]
    lines += ["## Table Captions"]
    for _, r in table_df.iterrows():
        lines.append(f"- {r['table_title']} ({r['table_file']}): {r['paper_ready_caption']} Boundary: {r['claim_boundary']}")
    lines += ["", "## Section Map"]
    for _, r in section_df.iterrows():
        lines.append(f"- {r['manuscript_section']}: {r['main_message']} Assets: {r['primary_assets']} Write now: {r['write_now_sentence']} Do not write: {r['do_not_write']}")
    return "\n".join(lines)


def copy_table(df: pd.DataFrame, name: str) -> None:
    out = atlas.TAB / name
    hub = atlas.HUB / name
    live = atlas.LIVE / name
    df.to_csv(out, sep="\t", index=False)
    df.to_csv(hub, sep="\t", index=False)
    shutil.copy2(hub, live)


def update_index_handoff_launchpad() -> None:
    idx = atlas.HUB / "index.html"
    if idx.exists():
        text = idx.read_text(encoding="utf-8")
        if HTML_NAME not in text:
            needle = '    <a href="kthyro_ctc_emt_representative_overview.html"'
            insert = f'    <a href="{HTML_NAME}" style="background:linear-gradient(135deg,#fffaf1 0%,#f0dfc4 100%);border:3px solid #b58534;color:#102033;font-weight:1000;font-size:18px;padding:16px 20px">★★ CTC-EMT PAPER-READY FIGURE/TABLE PACK · legends sections data</a>\n'
            pos = text.find(needle)
            text = text[:pos] + insert + text[pos:] if pos >= 0 else text + "\n" + insert
            idx.write_text(text, encoding="utf-8")
        shutil.copy2(idx, atlas.LIVE / "index.html")

    for path in [atlas.REP / "VISUAL_WEB_HANDOFF_KR.md", atlas.HUB / "VISUAL_WEB_HANDOFF_KR.md"]:
        if not path.exists():
            continue
        text = path.read_text(encoding="utf-8")
        line = "- Paper-Ready Figure/Table Pack: http://40.82.129.113:8012/papers_hub_2026_05_04/kthyro_ctc_emt_paper_ready_figure_table_pack.html\n"
        if HTML_NAME not in text:
            text = text.replace("## 최우선 링크\n\n", "## 최우선 링크\n\n" + line)
            path.write_text(text, encoding="utf-8")
    handoff = atlas.HUB / "VISUAL_WEB_HANDOFF_KR.md"
    if handoff.exists():
        shutil.copy2(handoff, atlas.LIVE / "VISUAL_WEB_HANDOFF_KR.md")

    launch = atlas.HUB / "kthyro_ctc_emt_all_pages_launchpad.html"
    if launch.exists() and HTML_NAME not in launch.read_text(encoding="utf-8"):
        text = launch.read_text(encoding="utf-8")
        insert = f'<a href="{HTML_NAME}"><b>원고 조립</b><span>Figure/Table Pack</span></a>'
        text = text.replace('<div class="route-grid">', '<div class="route-grid">' + insert, 1)
        launch.write_text(text, encoding="utf-8")
        shutil.copy2(launch, atlas.LIVE / "kthyro_ctc_emt_all_pages_launchpad.html")


def main() -> None:
    fig_dir = atlas.OUT / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    asset_dir = atlas.HUB / "assets" / atlas.ASSET
    live_asset_dir = atlas.LIVE / "assets" / atlas.ASSET
    flow_path = fig_dir / FLOW_FIG
    draw_flow_png(flow_path)
    shutil.copy2(flow_path, asset_dir / FLOW_FIG)
    shutil.copy2(flow_path, live_asset_dir / FLOW_FIG)

    fig_df = build_figure_legends()
    table_df = build_table_legends()
    section_df = pd.DataFrame(SECTION_MAP)
    copy_table(fig_df, FIG_LEGENDS_TSV)
    copy_table(table_df, TABLE_LEGENDS_TSV)
    copy_table(section_df, SECTION_MAP_TSV)

    html_text = html_page(fig_df, table_df, section_df)
    src = atlas.HUB / HTML_NAME
    live = atlas.LIVE / HTML_NAME
    local = atlas.REP / HTML_NAME
    src.write_text(html_text, encoding="utf-8")
    local.write_text(html_text, encoding="utf-8")
    shutil.copy2(src, live)

    md = report_md(fig_df, table_df, section_df)
    report = atlas.REP / REPORT_NAME
    hub_report = atlas.HUB / REPORT_NAME
    report.write_text(md, encoding="utf-8")
    hub_report.write_text(md, encoding="utf-8")
    shutil.copy2(hub_report, atlas.LIVE / REPORT_NAME)

    update_index_handoff_launchpad()
    print(f"HTML={src}")
    print(f"LIVE={live}")
    print(f"FIGURE={flow_path}")
    print(f"REPORT={report}")
    print(f"TABLES={atlas.TAB / FIG_LEGENDS_TSV}, {atlas.TAB / TABLE_LEGENDS_TSV}, {atlas.TAB / SECTION_MAP_TSV}")


if __name__ == "__main__":
    main()
