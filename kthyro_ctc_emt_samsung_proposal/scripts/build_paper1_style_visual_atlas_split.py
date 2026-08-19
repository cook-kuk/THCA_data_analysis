#!/usr/bin/env python3
from __future__ import annotations

import shutil
from pathlib import Path

import build_paper1_style_visual_atlas as base


PAGES = {
    "reader": "kthyro_ctc_emt_connected_reader.html",
    "hub": "kthyro_ctc_emt_visual_atlas_split_hub.html",
    "manuscript": "kthyro_ctc_emt_visual_atlas_manuscript.html",
    "evidence": "kthyro_ctc_emt_visual_atlas_evidence.html",
    "defense": "kthyro_ctc_emt_visual_atlas_defense.html",
    "deck": "kthyro_ctc_emt_visual_atlas_deck.html",
    "tables": "kthyro_ctc_emt_visual_atlas_tables.html",
}


FRONT_FIGS = [
    ("PF11_manuscript_blueprint_board.png", "Manuscript Blueprint Board"),
    ("PF09_manuscript_claim_scorecard.png", "Claim-Evidence Scorecard"),
    ("PF14_analysis_control_board.png", "Analysis Control Board"),
    ("PF13_execution_board.png", "Execution Board"),
    ("PF12_future_data_unlock_map.png", "Future Data Unlock Map"),
]
EVIDENCE_FIGS = [
    ("F03_tcga_braf_state_split.png", "BRAF-mutant PTC is not one state"),
    ("F04_tcga_driver_variance_boundary.png", "Driver variance boundary"),
    ("F05_spatial_tissue_state_organization.png", "Spatial tissue states are organized"),
    ("F06_scrna_marker_context.png", "Single-cell marker context"),
    ("F07_proteomics_dediff_direction.png", "Proteomics dedifferentiation direction"),
    ("PF03_marker_to_ngs_map.png", "Marker-to-NGS map"),
    ("F08_ctc_emt_ngs_prior_panel.png", "CTC-EMT-NGS prior panel"),
]
DEFENSE_FIGS = [
    ("PF02_data_provenance_map.png", "Data provenance map"),
    ("PF04_dataset_gap_bridge.png", "Missing dataset bridge"),
    ("PF05_claim_boundary_matrix.png", "Claim boundary matrix"),
    ("PF07_reviewer_to_figure_map.png", "Reviewer attack map"),
    ("F09_publishability_decision.png", "Publishability decision"),
    ("PF10_marker_module_support_heatmap.png", "Supplementary planning matrix - marker module planning scores (hand-coded; TCGA has no CTC; Yu 2024 = published anchor only)"),
]
DECK_FIGS = [
    ("slide01_problem.png", "Slide 1. The Narrow Problem"),
    ("slide02_hypothesis.png", "Slide 2. Core Hypothesis"),
    ("slide03_literature_gap.png", "Slide 3. Yu 2024 anchor and gap"),
    ("slide04_missing_dataset.png", "Slide 4. Missing serial dataset"),
    ("slide05_driver_not_enough.png", "Slide 5. Driver mutation is not enough"),
    ("slide06_spatial_context.png", "Slide 6. Tissue context"),
    ("slide07_ctc_ngs_panel.png", "Slide 7. CTC-EMT-NGS panel"),
    ("slide08_sampling_timeline.png", "Slide 8. Sampling timeline"),
    ("slide09_genetic_map.png", "Slide 9. NGS genetic map"),
    ("slide10_claim_boundary.png", "Slide 10. Claim boundaries"),
    ("slide11_budget_roadmap.png", "Slide 11. Budget roadmap"),
    ("slide12_final_ask.png", "Slide 12. Final ask"),
]


def css() -> str:
    return """
    :root{--ink:#102033;--muted:#52606f;--paper:#f6f0e6;--card:#fffaf1;--line:#d4c6b3;--line2:#eadfcc;--red:#8f2d25;--blue:#244e73;--green:#426b50;--gold:#b58534;--coal:#17212f;--good:#e9f3e8;--bad:#f7dedb}
    *{box-sizing:border-box} html{scroll-behavior:smooth}
    body{margin:0;background:linear-gradient(135deg,#f8f2e8 0%,#efe3d1 56%,#f7efe0 100%);color:var(--ink);font-family:"Newsreader","Noto Sans KR",serif;line-height:1.66;font-size:17px}
    a{color:var(--red);text-decoration:none;border-bottom:1px solid rgba(143,45,37,.32)}
    .wrap{max-width:1180px;margin:0 auto;padding:34px 28px 80px}
    .topnav{display:flex;gap:10px;align-items:center;flex-wrap:wrap;margin-bottom:24px;font-family:"JetBrains Mono",monospace;font-size:12px}
    .topnav a{border:1px solid var(--line);background:rgba(255,255,255,.62);padding:8px 11px;border-radius:999px;color:var(--ink)}
    .topnav a.active{background:var(--coal);color:#fff;border-color:var(--coal)}
    .hero{border:2px solid var(--coal);background:linear-gradient(145deg,#fffaf1 0%,#f0dfc4 100%);padding:42px 46px;box-shadow:10px 10px 0 rgba(23,33,47,.16)}
    .kicker,.mono{font-family:"JetBrains Mono",monospace}.kicker{letter-spacing:.17em;text-transform:uppercase;font-size:12px;color:var(--red);font-weight:700;margin-bottom:16px}
    h1{font-family:"Cormorant Garamond",serif;font-size:56px;line-height:1.02;margin:0 0 12px;color:var(--coal)}
    .subtitle{font-size:22px;color:var(--muted);max-width:960px;font-style:italic}
    .chips{display:flex;gap:8px;flex-wrap:wrap;margin-top:18px}.chip{font-family:"JetBrains Mono",monospace;font-size:11px;border:1px solid var(--line);background:#fff;padding:7px 9px;border-radius:6px}.chip.hot{background:var(--bad);color:var(--red);font-weight:800}.chip.good{background:var(--good);color:var(--green);font-weight:800}
    .layout{display:grid;grid-template-columns:250px 1fr;gap:26px;margin-top:30px}.toc{position:sticky;top:18px;align-self:start;background:rgba(255,250,241,.9);border:1px solid var(--line);border-radius:12px;padding:16px}.toc h3{font-family:"Cormorant Garamond",serif;margin:0 0 10px;font-size:24px;color:var(--red)}.toc a{display:block;color:var(--ink);padding:8px 0;border-bottom:1px dotted var(--line);font-size:14px}
    section{background:rgba(255,250,241,.88);border:1px solid var(--line);padding:30px;margin-bottom:24px;border-radius:14px}h2{font-family:"Cormorant Garamond",serif;font-size:38px;line-height:1.08;margin:0 0 10px;color:var(--coal)}h3{font-family:"Cormorant Garamond",serif;font-size:27px;margin:24px 0 10px;color:var(--red)}.section-note{color:var(--muted);font-style:italic;margin-bottom:18px}
    .grid{display:grid;gap:16px}.g2{grid-template-columns:repeat(2,minmax(0,1fr))}.g3{grid-template-columns:repeat(3,minmax(0,1fr))}.g4{grid-template-columns:repeat(4,minmax(0,1fr))}
    .card{display:block;background:#fff;border:1px solid var(--line2);padding:18px;border-radius:12px}.card h4{margin:0 0 8px;font-family:"JetBrains Mono",monospace;font-size:12px;letter-spacing:.08em;text-transform:uppercase;color:var(--red)}.big{font-family:"Cormorant Garamond",serif;font-size:40px;font-weight:700;color:var(--red);line-height:1}.small{font-size:14px;color:var(--muted)}.callout{border-left:5px solid var(--red);background:#fff;padding:18px 20px;margin:18px 0;border-radius:0 10px 10px 0}.callout.good{border-color:var(--green);background:var(--good)}
    figure.figure-block{margin:26px 0;background:#fff;border:1px solid var(--line);border-radius:14px;padding:16px;box-shadow:0 6px 20px rgba(23,33,47,.06)}.fig-title-row{display:flex;gap:12px;align-items:flex-start;border-bottom:1px solid var(--line2);margin-bottom:14px;padding-bottom:8px}.fig-title-row h3{margin:0;line-height:1.1}.fig-id{font-family:"JetBrains Mono",monospace;font-size:11px;background:var(--coal);color:#fff;border-radius:999px;padding:6px 9px;white-space:nowrap}
    figure img{display:block;width:100%;height:auto;border-radius:8px;background:#f4eadb;cursor:zoom-in}figcaption{font-size:14px;color:var(--muted);padding:12px 6px 4px;line-height:1.55}figcaption strong{color:var(--ink)}
    .fig-notes{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:8px;margin-top:10px}.fig-notes div{background:#fff8e8;border:1px solid var(--line2);border-radius:8px;padding:10px;font-size:14px;line-height:1.45}.fig-notes b{font-family:"JetBrains Mono",monospace;font-size:10px;text-transform:uppercase;letter-spacing:.08em;color:var(--red);display:block;margin-bottom:4px}.fig-notes .limit{background:#f8ead7}.fig-notes .wide-note{grid-column:1/-1;background:var(--good)}
    .table-wrap{overflow-x:auto}table{width:100%;border-collapse:collapse;margin:18px 0;background:#fff;font-size:14px;border:1px solid var(--line)}th{background:var(--coal);color:#fff;text-align:left;font-family:"JetBrains Mono",monospace;font-size:11px;letter-spacing:.06em;text-transform:uppercase;padding:10px}td{border-top:1px solid var(--line2);padding:10px;vertical-align:top}.footer{font-family:"JetBrains Mono",monospace;font-size:12px;color:var(--muted);text-align:center;margin-top:34px;padding-top:18px;border-top:1px solid var(--line)}
    .lightbox-overlay{position:fixed;inset:0;background:rgba(10,16,24,.88);z-index:9999;display:none;align-items:center;justify-content:center;padding:34px;cursor:zoom-out}.lightbox-overlay.show{display:flex}.lightbox-overlay img{max-width:96vw;max-height:88vh;object-fit:contain;background:#fff;border-radius:10px;box-shadow:0 20px 80px rgba(0,0,0,.48)}.lightbox-close{position:fixed;top:18px;right:24px;color:#fff;font-size:42px;line-height:1;font-family:Arial,sans-serif}.lightbox-caption{position:fixed;left:50%;bottom:18px;transform:translateX(-50%);max-width:min(900px,88vw);background:rgba(255,250,241,.96);color:var(--ink);border-radius:10px;padding:10px 14px;font-size:13px}
    @media(max-width:900px){.layout{grid-template-columns:1fr}.toc{position:relative;top:auto}.g2,.g3,.g4,.fig-notes{grid-template-columns:1fr}h1{font-size:39px}.hero{padding:28px 22px}.wrap{padding:18px 14px}.lightbox-caption{display:none}}
    """


def nav(active: str) -> str:
    links = [
        ("reader", "Connected Reader"),
        ("hub", "Atlas Hub"),
        ("manuscript", "Manuscript"),
        ("evidence", "Evidence"),
        ("defense", "Defense"),
        ("deck", "Deck"),
        ("tables", "Tables"),
    ]
    items = ['<a href="index.html">8-Papers Hub</a>', '<a href="paper1.html">paper1 style ref</a>']
    for key, label in links:
        cls = " active" if key == active else ""
        items.append(f'<a class="{cls.strip()}" href="{PAGES[key]}">{label}</a>')
    return "<nav class='topnav'>" + "\n".join(items) + "</nav>"


def page(title: str, subtitle: str, active: str, body: str, toc: str = "") -> str:
    toc_html = toc or "<a href='#top'>Top</a>"
    return f"""<!DOCTYPE html>
<html lang="ko">
<head>
<meta charset="UTF-8" />
<meta name="viewport" content="width=device-width, initial-scale=1.0" />
<title>{base.esc(title)} | K-Thyro CTC-EMT</title>
<link rel="preconnect" href="https://fonts.googleapis.com" />
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin />
<link href="https://fonts.googleapis.com/css2?family=Cormorant+Garamond:wght@500;600;700&family=Newsreader:ital,wght@0,400;0,500;0,600;1,400&family=JetBrains+Mono:wght@400;500;700&family=Noto+Sans+KR:wght@400;500;700;800&display=swap" rel="stylesheet" />
<style>{css()}</style>
</head>
<body>
<div class="wrap" id="top">
  {nav(active)}
  <header class="hero">
    <div class="kicker">K-Thyro CTC-EMT · split visual atlas · paper1-style</div>
    <h1>{base.esc(title)}</h1>
    <p class="subtitle">{base.esc(subtitle)}</p>
    <div class="chips">
      <span class="chip hot">Do not overclaim CTC biology from public data</span>
      <span class="chip good">Use as figure-by-figure manuscript scaffold</span>
      <span class="chip">Caption And Results Scaffold per figure</span>
    </div>
  </header>
  <div class="layout">
    <aside class="toc"><h3>Contents</h3>{toc_html}</aside>
    <main>{body}</main>
  </div>
  <div class="footer">Generated from kthyro_ctc_emt_paper_flow assets. Current claim remains a public-data prior/framework paper.</div>
</div>
<div class="lightbox-overlay" id="lightbox"><span class="lightbox-close">&times;</span><img alt="expanded figure"/><div class="lightbox-caption"></div></div>
<script>
const lb=document.getElementById('lightbox');const lbImg=lb.querySelector('img');const lbCap=lb.querySelector('.lightbox-caption');
document.querySelectorAll('figure img').forEach(img=>img.addEventListener('click',e=>{{e.preventDefault();lbImg.src=img.src;const fig=img.closest('figure');lbCap.textContent=fig?fig.querySelector('h3').textContent:img.alt;lb.classList.add('show');}}));
lb.addEventListener('click',()=>lb.classList.remove('show'));
</script>
</body>
</html>"""


def load_context() -> tuple[dict[str, dict], dict[str, dict]]:
    explain = base.row_by_image(base.read_tsv("figure_explanation_master.tsv"))
    captions = base.row_by_image(base.read_tsv("figure_caption_and_results_scaffold.tsv"))
    return explain, captions


def fig_list(items: list[tuple[str, str]], explain: dict[str, dict], captions: dict[str, dict]) -> str:
    return "\n".join(base.fig_block(img, title, explain, captions) for img, title in items)


def section(id_: str, title: str, note: str, body: str) -> str:
    return base.section(id_, title, note, body)


def hub_body() -> str:
    cards = [
        ("manuscript", "Manuscript Boards", "논문 identity, claim scorecard, execution board를 먼저 본다.", len(FRONT_FIGS)),
        ("evidence", "Evidence Figures", "현재 public-data로 쓸 수 있는 실제 in silico 결과.", len(EVIDENCE_FIGS)),
        ("defense", "Defense Figures", "심사위원 공격, 데이터 gap, claim boundary 방어.", len(DEFENSE_FIGS)),
        ("deck", "Samsung Deck Atlas", "PPT 슬라이드 그림을 큰 화면으로 분리 검토.", len(DECK_FIGS)),
        ("tables", "Control Tables", "숫자/CRF/SAP/kill criteria/safe sentence 관리.", 8),
    ]
    card_html = []
    for key, title, desc, n in cards:
        card_html.append(
            f"""<a class="card" href="{PAGES[key]}">
              <h4>{base.esc(title)}</h4>
              <div class="big">{n}</div>
              <p class="small">{base.esc(desc)}</p>
            </a>"""
        )
    return section(
        "hub",
        "Split Atlas Hub",
        "CV2 contact sheet 대신, 검토 목적별로 분리한 paper1-style visual atlas.",
        f"""
        <div class="callout good"><b>추천 순서.</b> Manuscript Boards → Evidence Figures → Defense Figures → Tables → Deck. 한 페이지에서 다 보려면 기존 one-page atlas를 열면 된다.</div>
        <div class="grid g3">{''.join(card_html)}</div>
        <h3>Direct Links</h3>
        <div class="grid g3">
          <a class="card" href="{PAGES['reader']}"><h4>Connected reader</h4><p>5분/15분/회의용 경로로 전체 자료를 연결해서 보기.</p></a>
          <a class="card" href="kthyro_ctc_emt_manuscript_storyboard_paper1_style.html"><h4>Manuscript storyboard</h4><p>Figure 1-6 → Results paragraph → claim boundary 흐름.</p></a>
          <a class="card" href="kthyro_ctc_emt_visual_atlas_paper1_style.html"><h4>One-page atlas</h4><p>모든 그림을 한 페이지에서 보는 긴 버전.</p></a>
          <a class="card" href="kthyro_ctc_emt_paper_flow_master_dossier.html"><h4>Master dossier</h4><p>표/claim/논문 흐름 전체 통제판.</p></a>
          <a class="card" href="assets/{base.ASSET}/PF08_cv2_visual_atlas.png"><h4>Old CV2 sheet</h4><p>원래 한 장짜리 contact sheet.</p></a>
        </div>
        """,
    )


def tables_body() -> str:
    tables = [
        ("key", "Key Number Ledger", "본문 숫자 drift 방지.", base.read_tsv("key_number_ledger.tsv")),
        ("layout", "Main/Supp Figure Layout", "논문 figure 배치.", base.read_tsv("main_supp_figure_layout.tsv")),
        ("sap", "Statistical Analysis Plan", "현재 분석과 미래 병원-data 분석 분리.", base.read_tsv("statistical_analysis_plan.tsv")),
        ("gates", "Stage-Gate Kill Criteria", "실패 시 downgrade/kill 기준.", base.read_tsv("stage_gate_kill_criteria.tsv")),
        ("crf", "Hospital Data Request CRF", "유형원 교수님/병원에 요청할 필드.", base.read_tsv("hospital_data_request_crf.tsv")),
        ("risk", "Reviewer Risk Heatmap", "예상 공격의 severity와 대응.", base.read_tsv("reviewer_risk_heatmap.tsv")),
        ("safe", "Safe Sentence Bank", "논문에 바로 넣을 수 있는 안전 문장.", base.read_tsv("safe_sentence_bank.tsv")),
        ("scenario", "Scenario Decision Tree", "Yu/hospital data 유무별 manuscript 경로.", base.read_tsv("scenario_decision_tree.tsv")),
    ]
    chunks = []
    for anchor, title, note, df in tables:
        chunks.append(section(anchor, title, note, base.table_html(df, 40)))
    return "\n".join(chunks)


def write_page(name: str, html: str) -> None:
    src = base.HUB / name
    live = base.LIVE / name
    local = base.REP / name
    src.write_text(html, encoding="utf-8")
    local.write_text(html, encoding="utf-8")
    shutil.copy2(src, live)


def update_index() -> None:
    idx = base.HUB / "index.html"
    if not idx.exists():
        return
    text = idx.read_text(encoding="utf-8")
    if PAGES["hub"] not in text:
        needle = '    <a href="kthyro_ctc_emt_visual_atlas_paper1_style.html"'
        insert = f'    <a href="{PAGES["hub"]}" style="background:linear-gradient(135deg,#fffaf1 0%,#f0dfc4 100%);border:2px solid #244e73;color:#102033;font-weight:900;font-size:16px;padding:14px 18px">★ CTC-EMT SPLIT VISUAL ATLAS · readable section pages</a>\n'
        pos = text.find(needle)
        text = text[:pos] + insert + text[pos:] if pos >= 0 else text + "\n" + insert
        idx.write_text(text, encoding="utf-8")
    shutil.copy2(idx, base.LIVE / "index.html")


def main() -> None:
    explain, captions = load_context()
    write_page(
        PAGES["hub"],
        page("CTC-EMT-NGS Split Visual Atlas", "보기 어려운 CV2 wall을 목적별 하위 페이지로 쪼갔다. 각 하위 페이지는 큰 그림, 결과 scaffold, 데이터 출처, claim boundary를 같이 보여준다.", "hub", hub_body()),
    )
    write_page(
        PAGES["manuscript"],
        page("Manuscript Boards", "논문이 무엇인지, 무엇이 아닌지, 어떤 순서로 써야 하는지를 먼저 고정하는 페이지.", "manuscript", section("figures", "Manuscript Boards", "Identity, claim strength, execution, future data unlock.", fig_list(FRONT_FIGS, explain, captions)), "<a href='#figures'>Figures</a>"),
    )
    write_page(
        PAGES["evidence"],
        page("Evidence Figures", "현재 공개데이터로 쓸 수 있는 in silico 결과만 모은 페이지. CTC transition discovery로 과장하지 않는다.", "evidence", section("figures", "Evidence Figures", "Mutation/state heterogeneity, spatial organization, marker context, prior-panel logic.", fig_list(EVIDENCE_FIGS, explain, captions)), "<a href='#figures'>Figures</a>"),
    )
    write_page(
        PAGES["defense"],
        page("Reviewer Defense Figures", "삼성 Science 리뷰어가 공격할 지점을 먼저 열어두고, 어디까지 주장 가능한지 분리한다.", "defense", section("figures", "Reviewer Defense Figures", "Provenance, dataset gap, claim boundary, reviewer attack map.", fig_list(DEFENSE_FIGS, explain, captions)), "<a href='#figures'>Figures</a>"),
    )
    write_page(
        PAGES["deck"],
        page("Samsung Deck Atlas", "12장 PPT 슬라이드 그림을 크게 분리해 검토하는 페이지.", "deck", section("figures", "Samsung Deck Figures", "Science-track proposal flow: problem, perturbation, Yu anchor, missing dataset, design, NGS map, claim boundaries.", fig_list(DECK_FIGS, explain, captions)), "<a href='#figures'>Figures</a>"),
    )
    write_page(
        PAGES["tables"],
        page("Control Tables", "숫자, 분석계획, CRF, kill criteria, safe sentence를 표로 통제하는 페이지.", "tables", tables_body(), "".join(f"<a href='#{a}'>{base.esc(t)}</a>" for a, t in [
            ("key", "Key numbers"), ("layout", "Figure layout"), ("sap", "SAP"), ("gates", "Gates"), ("crf", "CRF"), ("risk", "Risk"), ("safe", "Safe sentences"), ("scenario", "Scenario")
        ])),
    )
    update_index()
    print(f"HUB={base.HUB / PAGES['hub']}")
    for key, name in PAGES.items():
        print(f"{key}={base.LIVE / name}")


if __name__ == "__main__":
    main()
