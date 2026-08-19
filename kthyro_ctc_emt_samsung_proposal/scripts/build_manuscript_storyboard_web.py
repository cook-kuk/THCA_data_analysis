#!/usr/bin/env python3
from __future__ import annotations

import html
import shutil
from pathlib import Path

import pandas as pd

import build_paper1_style_visual_atlas as atlas
import build_paper1_style_visual_atlas_split as split


HTML_NAME = "kthyro_ctc_emt_manuscript_storyboard_paper1_style.html"
TABLE_NAME = "manuscript_storyboard_web.tsv"
REPORT_NAME = "MANUSCRIPT_STORYBOARD_WEB_KR.md"


def esc(x) -> str:
    return html.escape("" if x is None else str(x), quote=True)


STORY_ROWS = [
    {
        "slot": "Main Fig. 1",
        "title": "Study Logic And Data Provenance",
        "assets": "PF01_manuscript_flow.png;PF02_data_provenance_map.png",
        "results_heading": "A staged public-data prior map for a future serial CTC-EMT perturbation study",
        "write_now": "We first separated the evidence layers that can be analyzed now from the prospective serial blood experiment that remains to be performed. Public TCGA, spatial, single-cell, and proteomic resources support tissue-state and marker-prior logic; Yu 2024 supports feasibility of serial CTC-EMT measurement; hospital serial CTC/cfDNA/tissue NGS data are required for the actual transition test.",
        "data_used": "TCGA-THCA public pilot; GSE250521 spatial summary; scRNA/proteomic marker summaries; Yu 2024 published CTC feasibility; local proposal reports.",
        "claim_boundary": "This is a framework/provenance figure. It does not claim that public data measured postoperative CTC transition.",
        "reviewer_answer": "The study is honest because it explicitly labels which evidence is public prior, which is published feasibility, and which is future prospective validation.",
    },
    {
        "slot": "Main Fig. 2",
        "title": "Driver Mutation Is Not Sufficient To Define State",
        "assets": "F03_tcga_braf_state_split.png;F04_tcga_driver_variance_boundary.png",
        "results_heading": "BRAF-mutant PTC remains state-heterogeneous in public tumor data",
        "write_now": "In the TCGA-THCA public pilot context, BRAF-mutant tumors distributed across multiple state/vulnerability labels, with the largest BRAF label accounting for only 33.2% of the BRAF-mutant subset. This supports the need to interpret postoperative blood signals with phenotype-state information rather than driver genotype alone.",
        "data_used": "TCGA-THCA bulk tumor RNA/mutation/clinical outputs; BRAF-only vulnerability label distribution; driver/state variance models.",
        "claim_boundary": "This justifies genotype-plus-state interpretation. It is not evidence that CTC EMT states predict recurrence.",
        "reviewer_answer": "The result is narrow: mutation alone is insufficient as a state descriptor in tumor tissue; no blood biology is inferred directly.",
    },
    {
        "slot": "Main Fig. 3",
        "title": "Spatial Tissue-State Organization",
        "assets": "F05_spatial_tissue_state_organization.png;spatial_coherence_summary.png;spatial_niche_adjacency_enrichment_heatmap.png",
        "results_heading": "Public spatial data support non-random tissue-state organization",
        "write_now": "Spatial thyroid tissue data showed organized tissue-state structure rather than random mixing. In the current package, the spatial layer is used to support ROI and marker-selection logic for future CTC-EMT studies, not to prove the origin of circulating cells.",
        "data_used": "GSE250521 spatial transcriptomics summary; 16 slides; 57,144 spots; spatial coherence, adjacency, hotspot, and interface outputs.",
        "claim_boundary": "Spatial organization is a tissue-context prior only. It does not identify the shedding source of postoperative CTCs.",
        "reviewer_answer": "The spatial analysis is deliberately positioned as marker/ROI rationale, not as direct CTC evidence.",
    },
    {
        "slot": "Main Fig. 4",
        "title": "Marker Module Support Across Data Layers",
        "assets": "F06_scrna_marker_context.png;F07_proteomics_dediff_direction.png",
        "results_heading": "Candidate CTC-EMT marker modules have interpretable tissue and molecular context",
        "write_now": "Single-cell and proteomic summaries provide supporting context for epithelial, mesenchymal, thyroid-lineage, stress/survival, and immune-evasion marker modules. These data support rational panel construction for future CTC-EMT phenotyping and NGS anchoring. A separate supplementary planning matrix (PF10) records our hand-coded judgement of how useful each public/published evidence layer is for prior construction; that matrix is a planning score, not a data analysis result.",
        "data_used": "PTC scRNA marker-context summary; bulk proteomics proxy. The supplementary planning matrix (PF10) is shown separately and is not part of this figure.",
        "claim_boundary": "This is marker prioritization, not a validated clinical CTC assay and not a completed blood-based discovery. PF10 (supplementary) is a manuscript-planning judgement, not a cross-layer measurement.",
        "reviewer_answer": "The marker panel is hypothesis-generating but biologically organized across independent public layers; the supplementary PF10 planning matrix records reviewer-visible judgement and does not introduce new data.",
    },
    {
        "slot": "Supplementary Planning Matrix",
        "title": "PF10 - Marker Module Planning Support Scores (not data-derived)",
        "assets": "PF10_marker_module_support_heatmap.png",
        "results_heading": "Supplementary, reviewer-visible record of manuscript-planning judgement on marker-module priors",
        "write_now": "As a supplementary planning matrix, we hand-coded (0/1/2) how useful each public/published evidence layer was for constructing each marker-module prior. These integers are manuscript-planning judgements, not measurements or validated statistics. The TCGA column records judgement about TCGA-THCA tumor tissue context only (TCGA contains no CTC, blood, or liquid biopsy data). The Yu 2024 column scores the published Yu 2024 paper's reported feasibility summary (62 PTC, 87% CTC detection); no unpublished or raw Yu CTC data was used.",
        "data_used": "Hand-coded planning matrix marker_module_cross_layer_support_matrix.tsv; reference layers cited as context only (TCGA-THCA tissue [no CTC], GSE250521 spatial, scRNA module summary, bulk proteomics proxy, Yu 2024 published anchor summary).",
        "claim_boundary": "These are planning scores, not validation metrics, not assay-performance values, and not a cross-layer measurement. TCGA does not contain CTC data and was never analyzed as such here.",
        "reviewer_answer": "We disclose this matrix as a supplementary planning artifact so reviewers can see the prior-construction judgement explicitly, and so it cannot be misread as a data-derived cross-layer analysis.",
    },
    {
        "slot": "Main Fig. 5",
        "title": "Phenotype-To-Genotype CTC-EMT-NGS Framework",
        "assets": "PF03_marker_to_ngs_map.png;F08_ctc_emt_ngs_prior_panel.png",
        "results_heading": "A matched genetic anchor is required to interpret postoperative blood signals",
        "write_now": "We define a research framework in which epithelial, hybrid E/M, and mesenchymal CTC phenotypes are interpreted together with matched tumor-normal NGS and serial cfDNA. The genetic anchor is essential because CTC phenotype alone may be nonspecific in early PTC.",
        "data_used": "Marker-to-NGS map; CTC-EMT-NGS prior panel; claim-boundary and Samsung proposal logic.",
        "claim_boundary": "The framework does not assume cfDNA or CTC-enriched NGS will work in all early PTC patients.",
        "reviewer_answer": "This figure lowers overclaim risk by making genotype anchoring a requirement rather than an afterthought.",
    },
    {
        "slot": "Main Fig. 6",
        "title": "Missing Dataset And Prospective Unlock",
        "assets": "PF04_dataset_gap_bridge.png;PF12_future_data_unlock_map.png;PF05_claim_boundary_matrix.png",
        "results_heading": "No public dataset currently tests the serial postoperative CTC-EMT genetic-map question",
        "write_now": "We found no public dataset that jointly provides serial pre/post-thyroidectomy blood, CTC EMT phenotype transition, matched tissue-normal NGS, cfDNA or CTC-enriched NGS, and postoperative follow-up. This absence defines the scientific need for the prospective Samsung/Yu cohort.",
        "data_used": "Dataset-gap bridge; future data unlock map; claim boundary matrix; hospital CRF table.",
        "claim_boundary": "The missing-dataset result justifies the prospective study; it is not a substitute for the prospective study.",
        "reviewer_answer": "The proposal is needed because the exact joint dataset required for the central biological question does not exist publicly.",
    },
    {
        "slot": "Supplement/Defense",
        "title": "Reviewer Defense And Execution Control",
        "assets": "PF07_reviewer_to_figure_map.png;F09_publishability_decision.png;PF13_execution_board.png;PF14_analysis_control_board.png",
        "results_heading": "Defense materials define the safe manuscript lane",
        "write_now": "The supplementary defense layer maps expected reviewer attacks to figures, sets publishability boundaries, and preserves execution control through key-number, statistical-plan, and kill-criteria tables.",
        "data_used": "Reviewer risk heatmap; publishability decision table; execution board; analysis control package.",
        "claim_boundary": "These are control/defense artifacts, not additional biological discoveries.",
        "reviewer_answer": "The defense package is included to prevent claim creep and to make the paper auditable.",
    },
]


def asset_img(name: str) -> str:
    return f"assets/{atlas.ASSET}/{esc(name)}"


def panel_images(assets: str) -> str:
    names = [x.strip() for x in assets.split(";") if x.strip()]
    imgs = []
    for name in names:
        imgs.append(
            f"""
            <a class="panel-img" href="{asset_img(name)}">
              <img src="{asset_img(name)}" alt="{esc(name)}" loading="lazy" />
              <span>{esc(name)}</span>
            </a>
            """
        )
    return "<div class='panel-grid'>" + "".join(imgs) + "</div>"


def story_block(row: dict) -> str:
    return f"""
    <section id="{esc(row['slot'].lower().replace(' ', '-').replace('.', '').replace('/', '-'))}">
      <div class="slot-row">
        <span>{esc(row['slot'])}</span>
        <h2>{esc(row['title'])}</h2>
      </div>
      {panel_images(row['assets'])}
      <div class="story-grid">
        <div><b>Results heading</b><p>{esc(row['results_heading'])}</p></div>
        <div><b>Write-now paragraph</b><p>{esc(row['write_now'])}</p></div>
        <div><b>Data used</b><p>{esc(row['data_used'])}</p></div>
        <div class="limit"><b>Claim boundary</b><p>{esc(row['claim_boundary'])}</p></div>
        <div class="wide-note"><b>Reviewer answer</b><p>{esc(row['reviewer_answer'])}</p></div>
      </div>
    </section>
    """


def storyboard_table() -> pd.DataFrame:
    return pd.DataFrame(STORY_ROWS)


def markdown_report(df: pd.DataFrame) -> str:
    lines = [
        "# Manuscript Storyboard Web",
        "",
        "Purpose: 논문 본문 흐름에 맞춰 figure, data used, results scaffold, claim boundary를 한 줄로 고정한다.",
        "",
        "Core boundary: 현재 public-data output은 CTC transition discovery가 아니라 tissue-state/marker/genetic-anchor prior framework다.",
        "",
    ]
    for _, r in df.iterrows():
        lines.extend(
            [
                f"## {r['slot']}. {r['title']}",
                f"- Results heading: {r['results_heading']}",
                f"- Data used: {r['data_used']}",
                f"- Write-now paragraph: {r['write_now']}",
                f"- Claim boundary: {r['claim_boundary']}",
                f"- Reviewer answer: {r['reviewer_answer']}",
                "",
            ]
        )
    return "\n".join(lines)


def html_page(df: pd.DataFrame) -> str:
    toc = "".join(
        f"<a href='#{esc(r['slot'].lower().replace(' ', '-').replace('.', '').replace('/', '-'))}'>{esc(r['slot'])}</a>"
        for _, r in df.iterrows()
    )
    cards = """
    <section id="identity">
      <h2>Paper Identity Lock</h2>
      <p class="section-note">This page is intentionally stricter than the deck. It decides exactly what can be written as a paper now.</p>
      <div class="grid g3">
        <div class="card"><h4>Write now</h4><div class="big">Prior map</div><p class="small">Public-data in silico tissue-state / marker / NGS-anchor framework.</p></div>
        <div class="card"><h4>Do not write now</h4><div class="big">CTC discovery</div><p class="small">No postoperative serial CTC transition or recurrence model without Yu/hospital data.</p></div>
        <div class="card"><h4>Unlock</h4><div class="big">Cohort</div><p class="small">Serial CTC phenotype + cfDNA + matched tissue-normal NGS + follow-up.</p></div>
      </div>
      <div class="callout good"><b>Practical rule.</b> Every Results paragraph must end inside its data layer. TCGA stays tissue. Spatial stays tissue context. Marker modules stay prioritization. CTC transition remains a prospective biological test.</div>
    </section>
    """
    body = cards + "\n".join(story_block(r.to_dict()) for _, r in df.iterrows())
    body += f"""
    <section id="downloads">
      <h2>Downloads And Adjacent Pages</h2>
      <div class="grid g3">
        <a class="card" href="{TABLE_NAME}"><h4>Storyboard TSV</h4><p>Machine-readable figure/story table.</p></a>
        <a class="card" href="{REPORT_NAME}"><h4>Storyboard MD</h4><p>Korean writing scaffold.</p></a>
        <a class="card" href="kthyro_ctc_emt_visual_atlas_split_hub.html"><h4>Split atlas</h4><p>Large figure-by-figure review pages.</p></a>
      </div>
    </section>
    """
    return split.page(
        "Manuscript Storyboard",
        "논문 Figure 1-6 흐름에 맞춰 data used, results scaffold, reviewer answer, claim boundary를 한 화면에서 고정한다.",
        "hub",
        body,
        "<a href='#identity'>Identity lock</a>" + toc + "<a href='#downloads'>Downloads</a>",
    ).replace(
        "</style>",
        """
        .slot-row{display:flex;align-items:flex-start;gap:14px;border-bottom:1px solid var(--line2);margin-bottom:16px;padding-bottom:10px}
        .slot-row span{font-family:"JetBrains Mono",monospace;background:var(--coal);color:#fff;border-radius:999px;padding:7px 10px;font-size:11px;white-space:nowrap}
        .slot-row h2{margin:0}
        .panel-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(270px,1fr));gap:12px;margin:16px 0}
        .panel-img{display:block;background:#fff;border:1px solid var(--line2);border-radius:10px;padding:10px}
        .panel-img img{display:block;width:100%;height:auto;border-radius:7px;background:#f4eadb}
        .panel-img span{display:block;font-family:"JetBrains Mono",monospace;font-size:10px;color:var(--muted);margin-top:6px;overflow-wrap:anywhere}
        .story-grid{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:10px;margin-top:12px}
        .story-grid div{background:#fff8e8;border:1px solid var(--line2);border-radius:8px;padding:12px}
        .story-grid b{font-family:"JetBrains Mono",monospace;font-size:10px;text-transform:uppercase;letter-spacing:.08em;color:var(--red);display:block;margin-bottom:4px}
        .story-grid p{margin:0}
        .story-grid .limit{background:#f8ead7}
        .story-grid .wide-note{grid-column:1/-1;background:var(--good)}
        @media(max-width:900px){.story-grid{grid-template-columns:1fr}.slot-row{display:block}.slot-row span{display:inline-block;margin-bottom:8px}}
        </style>""",
    )


def update_index() -> None:
    idx = atlas.HUB / "index.html"
    if not idx.exists():
        return
    text = idx.read_text(encoding="utf-8")
    if HTML_NAME not in text:
        needle = '    <a href="kthyro_ctc_emt_visual_atlas_split_hub.html"'
        insert = f'    <a href="{HTML_NAME}" style="background:linear-gradient(135deg,#fffaf1 0%,#f0dfc4 100%);border:2px solid #426b50;color:#102033;font-weight:900;font-size:16px;padding:14px 18px">★ CTC-EMT MANUSCRIPT STORYBOARD · figure-to-results flow</a>\n'
        pos = text.find(needle)
        text = text[:pos] + insert + text[pos:] if pos >= 0 else text + "\n" + insert
        idx.write_text(text, encoding="utf-8")
    shutil.copy2(idx, atlas.LIVE / "index.html")


def main() -> None:
    df = storyboard_table()
    table_path = atlas.HUB / TABLE_NAME
    live_table = atlas.LIVE / TABLE_NAME
    report_path = atlas.HUB / REPORT_NAME
    live_report = atlas.LIVE / REPORT_NAME
    out_table = atlas.TAB / TABLE_NAME
    out_report = atlas.REP / REPORT_NAME

    df.to_csv(table_path, sep="\t", index=False)
    df.to_csv(out_table, sep="\t", index=False)
    report = markdown_report(df)
    report_path.write_text(report, encoding="utf-8")
    out_report.write_text(report, encoding="utf-8")
    shutil.copy2(table_path, live_table)
    shutil.copy2(report_path, live_report)

    html_text = html_page(df)
    src = atlas.HUB / HTML_NAME
    local = atlas.REP / HTML_NAME
    live = atlas.LIVE / HTML_NAME
    src.write_text(html_text, encoding="utf-8")
    local.write_text(html_text, encoding="utf-8")
    shutil.copy2(src, live)
    update_index()

    print(f"HTML={src}")
    print(f"LIVE={live}")
    print(f"TABLE={out_table}")
    print(f"REPORT={out_report}")


if __name__ == "__main__":
    main()
