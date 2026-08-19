#!/usr/bin/env python3
from __future__ import annotations

import html
import shutil
from pathlib import Path

import pandas as pd


ROOT = Path("/home/seungho/personal/THCA_data_analysis")
PROP = ROOT / "kthyro_ctc_emt_samsung_proposal"
OUT = PROP / "outputs/paper_flow_master_dossier"
TAB = OUT / "tables"
REP = OUT / "reports"
HUB = ROOT / "project/papers_hub_2026_05_04"
LIVE = Path("/var/www/papers/papers_hub_2026_05_04")
ASSET = "kthyro_ctc_emt_paper_flow"
HTML = HUB / "kthyro_ctc_emt_visual_atlas_paper1_style.html"
LIVE_HTML = LIVE / "kthyro_ctc_emt_visual_atlas_paper1_style.html"
LOCAL_HTML = REP / "kthyro_ctc_emt_visual_atlas_paper1_style.html"


def esc(x) -> str:
    return html.escape("" if x is None else str(x), quote=True)


def read_tsv(name: str) -> pd.DataFrame:
    p = TAB / name
    if not p.exists():
        return pd.DataFrame()
    return pd.read_csv(p, sep="\t")


def row_by_image(df: pd.DataFrame) -> dict[str, dict]:
    if df.empty or "image" not in df.columns:
        return {}
    return {str(r["image"]): r.to_dict() for _, r in df.iterrows()}


def table_html(df: pd.DataFrame, max_rows: int = 50) -> str:
    if df.empty:
        return "<p class='small'>No table available.</p>"
    df = df.head(max_rows)
    out = ["<div class='table-wrap'><table><thead><tr>"]
    for c in df.columns:
        out.append(f"<th>{esc(c)}</th>")
    out.append("</tr></thead><tbody>")
    for _, row in df.iterrows():
        out.append("<tr>")
        for c in df.columns:
            out.append(f"<td>{esc(row[c])}</td>")
        out.append("</tr>")
    out.append("</tbody></table></div>")
    return "".join(out)


def fig_block(img: str, title: str, explain: dict[str, dict], captions: dict[str, dict], slot: str = "") -> str:
    e = explain.get(img, {})
    c = captions.get(img, {})
    figure_id = e.get("figure", slot or img.replace(".png", ""))
    role = e.get("role", "")
    data_used = e.get("data_used", "")
    how = e.get("how_to_read", "")
    paper_use = e.get("paper_use", "")
    narration = e.get("narration", "")
    boundary = e.get("boundary", "")
    long_caption = c.get("long_caption", role)
    results = c.get("results_scaffold", narration)
    methods = c.get("methods_note", data_used)
    reviewer = c.get("reviewer_sentence", boundary)
    return f"""
    <figure class="figure-block" id="{esc(figure_id)}">
      <div class="fig-title-row">
        <span class="fig-id">{esc(figure_id)}</span>
        <h3>{esc(title)}</h3>
      </div>
      <a href="assets/{ASSET}/{esc(img)}"><img src="assets/{ASSET}/{esc(img)}" alt="{esc(title)}" loading="lazy" /></a>
      <figcaption><strong>Caption And Results Scaffold.</strong> {esc(long_caption)}</figcaption>
      <div class="fig-notes">
        <div><b>Data used</b>{esc(data_used or methods)}</div>
        <div><b>How to read</b>{esc(how)}</div>
        <div><b>Results scaffold</b>{esc(results)}</div>
        <div><b>Methods / provenance</b>{esc(methods)}</div>
        <div><b>Paper use</b>{esc(paper_use)}</div>
        <div class="limit"><b>Claim boundary</b>{esc(boundary)}</div>
        <div class="wide-note"><b>Reviewer-safe sentence</b>{esc(reviewer)}</div>
      </div>
    </figure>
    """


def metric_card(label: str, value: str, note: str, cls: str = "") -> str:
    return f"<div class='card {esc(cls)}'><h4>{esc(label)}</h4><div class='big'>{esc(value)}</div><p class='small'>{esc(note)}</p></div>"


def section(id_: str, title: str, note: str, body: str) -> str:
    return f"""
    <section id="{esc(id_)}">
      <h2>{esc(title)}</h2>
      <p class="section-note">{esc(note)}</p>
      {body}
    </section>
    """


def build_html() -> str:
    explain_df = read_tsv("figure_explanation_master.tsv")
    captions_df = read_tsv("figure_caption_and_results_scaffold.tsv")
    key_numbers = read_tsv("key_number_ledger.tsv")
    main_figs = read_tsv("main_supp_figure_layout.tsv")
    sap = read_tsv("statistical_analysis_plan.tsv")
    gates = read_tsv("stage_gate_kill_criteria.tsv")
    crf = read_tsv("hospital_data_request_crf.tsv")
    risks = read_tsv("reviewer_risk_heatmap.tsv")
    safe = read_tsv("safe_sentence_bank.tsv")
    scenario = read_tsv("scenario_decision_tree.tsv")

    explain = row_by_image(explain_df)
    captions = row_by_image(captions_df)

    front_figs = [
        ("PF11_manuscript_blueprint_board.png", "Manuscript Blueprint Board"),
        ("PF09_manuscript_claim_scorecard.png", "Claim-Evidence Scorecard"),
        ("PF14_analysis_control_board.png", "Analysis Control Board"),
        ("PF13_execution_board.png", "Execution Board"),
        ("PF12_future_data_unlock_map.png", "Future Data Unlock Map"),
    ]
    evidence_figs = [
        ("F03_tcga_braf_state_split.png", "BRAF-mutant PTC is not one state"),
        ("F04_tcga_driver_variance_boundary.png", "Driver variance boundary"),
        ("F05_spatial_tissue_state_organization.png", "Spatial tissue states are organized"),
        ("F06_scrna_marker_context.png", "Single-cell marker context"),
        ("F07_proteomics_dediff_direction.png", "Proteomics dedifferentiation direction"),
        ("PF03_marker_to_ngs_map.png", "Marker-to-NGS map"),
        ("F08_ctc_emt_ngs_prior_panel.png", "CTC-EMT-NGS prior panel"),
    ]
    defense_figs = [
        ("PF02_data_provenance_map.png", "Data provenance map"),
        ("PF04_dataset_gap_bridge.png", "Missing dataset bridge"),
        ("PF05_claim_boundary_matrix.png", "Claim boundary matrix"),
        ("PF07_reviewer_to_figure_map.png", "Reviewer attack map"),
        ("F09_publishability_decision.png", "Publishability decision"),
        ("PF10_marker_module_support_heatmap.png", "Supplementary planning matrix - marker module planning scores (hand-coded; TCGA has no CTC; Yu 2024 = published anchor only)"),
    ]
    deck_figs = [(f"slide{i:02d}_{name}.png", title) for i, name, title in [
        (1, "problem", "Slide 1. The Narrow Problem"),
        (2, "hypothesis", "Slide 2. Core Hypothesis"),
        (3, "literature_gap", "Slide 3. Yu 2024 anchor and gap"),
        (4, "missing_dataset", "Slide 4. Missing serial dataset"),
        (5, "driver_not_enough", "Slide 5. Driver mutation is not enough"),
        (6, "spatial_context", "Slide 6. Tissue context"),
        (7, "ctc_ngs_panel", "Slide 7. CTC-EMT-NGS panel"),
        (8, "sampling_timeline", "Slide 8. Sampling timeline"),
        (9, "genetic_map", "Slide 9. NGS genetic map"),
        (10, "claim_boundary", "Slide 10. Claim boundaries"),
        (11, "budget_roadmap", "Slide 11. Budget roadmap"),
        (12, "final_ask", "Slide 12. Final ask"),
    ]]

    def fig_list(items: list[tuple[str, str]]) -> str:
        return "\n".join(fig_block(img, title, explain, captions) for img, title in items)

    hero_metrics = "".join(
        [
            metric_card("Page purpose", "Readable", "CV2 atlas split into paper1-style sections.", "reading-card"),
            metric_card("Current paper", "Prior map", "Public-data framework, not CTC discovery.", "reading-card"),
            metric_card("Figures", str(len(front_figs) + len(evidence_figs) + len(defense_figs) + len(deck_figs)), "Large separated figures with captions."),
            metric_card("Claim rule", "Bounded", "Every figure carries what not to claim.", "reading-card"),
        ]
    )

    body = f"""
    <div class="wrap">
      <nav class="topnav">
        <a href="index.html">8-Papers Hub</a>
        <a href="kthyro_ctc_emt_paper_flow_master_dossier.html">Master Dossier</a>
        <a href="assets/{ASSET}/PF08_cv2_visual_atlas.png">Old CV2 Atlas</a>
      </nav>
      <header class="hero">
        <div class="kicker">K-Thyro CTC-EMT · visual atlas · paper1-style readable page</div>
        <h1>CTC-EMT-NGS Visual Atlas</h1>
        <p class="subtitle">한 장짜리 CV2 wall을 버리고, 논문 검토용으로 큰 그림과 설명을 분리한 페이지. 각 그림은 caption, results paragraph, methods/provenance, reviewer-safe sentence, claim boundary를 바로 옆에 갖는다.</p>
        <div class="meta">
          <span class="chip hot">Current claim: public-data prior/framework paper</span>
          <span class="chip good">Future claim requires Yu/hospital serial CTC data</span>
          <span class="chip">Science-track, not ICT identity</span>
          <span class="chip">PTC · thyroidectomy · CTC-EMT · matched NGS</span>
        </div>
        <div class="grid g4" style="margin-top:24px">{hero_metrics}</div>
      </header>

      <div class="layout">
        <aside class="toc">
          <h3>Contents</h3>
          <a href="#how">How to read</a>
          <a href="#front">Manuscript boards</a>
          <a href="#evidence">Evidence figures</a>
          <a href="#defense">Defense figures</a>
          <a href="#deck">Samsung deck atlas</a>
          <a href="#tables">Control tables</a>
          <a href="#downloads">Downloads</a>
        </aside>
        <main>
          {section("how", "How To Read", "This page replaces the dense CV2 wall with readable figure-by-figure sections.", f'''
            <div class="callout good"><b>읽는 순서.</b> 먼저 Manuscript Blueprint와 Analysis Control Board로 논문 범위를 고정하고, Evidence figures에서 실제 public-data 근거를 본 뒤, Defense figures에서 심사위원 방어 문장을 확인한다.</div>
            <div class="logic-strip">
              <div><b>Now</b><p>public-data prior/framework paper</p></div>
              <div><b>Not now</b><p>postoperative CTC transition discovery</p></div>
              <div><b>Unlock</b><p>Yu/hospital serial CTC + matched NGS</p></div>
              <div><b>Rule</b><p>every result has a boundary</p></div>
            </div>
          ''')}
          {section("front", "Manuscript Boards", "Start here. These boards decide what the paper is, what it is not, and what to do next.", fig_list(front_figs))}
          {section("evidence", "Evidence Figures", "These are the public-data in silico results that can support the current paper.", fig_list(evidence_figs))}
          {section("defense", "Reviewer Defense Figures", "These figures keep the proposal/paper honest under Samsung-style review.", fig_list(defense_figs))}
          {section("deck", "Samsung Deck Atlas", "Large separated deck figures for transfer and review. These are not cramped into one contact sheet.", fig_list(deck_figs))}
          {section("tables", "Control Tables", "The working tables that prevent numerical drift and overclaiming.", f'''
            <h3>Key Number Ledger</h3>{table_html(key_numbers, 20)}
            <h3>Main/Supp Figure Layout</h3>{table_html(main_figs, 20)}
            <h3>Statistical Analysis Plan</h3>{table_html(sap, 20)}
            <h3>Stage-Gate Kill Criteria</h3>{table_html(gates, 20)}
            <h3>Hospital Data Request CRF</h3>{table_html(crf, 20)}
            <h3>Reviewer Risk Heatmap</h3>{table_html(risks, 20)}
            <h3>Safe Sentence Bank</h3>{table_html(safe, 20)}
            <h3>Scenario Decision Tree</h3>{table_html(scenario, 20)}
          ''')}
          {section("downloads", "Downloads", "Direct links for handoff.", f'''
            <div class="grid g3">
              <a class="card" href="assets/{ASSET}/kthyro_ctc_emt_paper_flow_packet_2026_05_09.zip"><h4>Packet ZIP</h4><p>All figures/tables/reports.</p></a>
              <a class="card" href="kthyro_ctc_emt_paper_flow_master_dossier.html"><h4>Master dossier</h4><p>Dense control page.</p></a>
              <a class="card" href="assets/{ASSET}/kthyro_ctc_emt_samsung_science_pitch_v2_reference_style.pptx"><h4>Samsung PPTX</h4><p>Science-track deck.</p></a>
            </div>
          ''')}
        </main>
      </div>
      <div class="footer">Generated from kthyro_ctc_emt_paper_flow assets. Public data supports prior-map logic only, not CTC transition discovery.</div>
    </div>
    <div class="lightbox-overlay" id="lightbox"><span class="lightbox-close">&times;</span><img alt="expanded figure"/><div class="lightbox-caption"></div></div>
    <script>
      const lb = document.getElementById('lightbox');
      const lbImg = lb.querySelector('img');
      const lbCap = lb.querySelector('.lightbox-caption');
      document.querySelectorAll('figure img').forEach(img => {{
        img.addEventListener('click', e => {{
          e.preventDefault();
          lbImg.src = img.src;
          const fig = img.closest('figure');
          lbCap.textContent = fig ? fig.querySelector('h3').textContent : img.alt;
          lb.classList.add('show');
        }});
      }});
      lb.addEventListener('click', () => lb.classList.remove('show'));
    </script>
    """

    return f"""<!DOCTYPE html>
<html lang="ko">
<head>
<meta charset="UTF-8" />
<meta name="viewport" content="width=device-width, initial-scale=1.0" />
<title>K-Thyro CTC-EMT Visual Atlas | Paper1 Style</title>
<link rel="preconnect" href="https://fonts.googleapis.com" />
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin />
<link href="https://fonts.googleapis.com/css2?family=Cormorant+Garamond:wght@500;600;700&family=Newsreader:ital,wght@0,400;0,500;0,600;1,400&family=JetBrains+Mono:wght@400;500;700&family=Noto+Sans+KR:wght@400;500;700;800&display=swap" rel="stylesheet" />
<style>
:root{{
  --ink:#102033; --muted:#52606f; --paper:#f6f0e6; --card:#fffaf1;
  --line:#d4c6b3; --line2:#eadfcc; --red:#8f2d25; --red2:#b64b3f;
  --blue:#244e73; --green:#426b50; --gold:#b58534; --coal:#17212f;
  --cream:#fff7e4; --good:#e9f3e8; --warn:#fff1d4; --bad:#f7dedb;
}}
*{{box-sizing:border-box}}
html{{scroll-behavior:smooth}}
body{{margin:0;background:radial-gradient(circle at 8% 0%, rgba(181,133,52,.18), transparent 28rem),linear-gradient(135deg,#f8f2e8 0%,#efe3d1 55%,#f7efe0 100%);color:var(--ink);font-family:"Newsreader","Noto Sans KR",serif;line-height:1.68;font-size:17px}}
a{{color:var(--red);text-decoration:none;border-bottom:1px solid rgba(143,45,37,.35)}}
.wrap{{max-width:1180px;margin:0 auto;padding:34px 28px 80px}}
.topnav{{display:flex;justify-content:space-between;gap:18px;align-items:center;font-family:"JetBrains Mono",monospace;font-size:12px;margin-bottom:28px;color:var(--muted);flex-wrap:wrap}}
.topnav a{{border:1px solid var(--line);padding:8px 12px;background:rgba(255,255,255,.55);border-radius:999px}}
.hero{{border:2px solid var(--coal);background:linear-gradient(145deg,#fffaf1 0%,#f0dfc4 100%);padding:44px 46px;box-shadow:10px 10px 0 rgba(23,33,47,.16);position:relative;overflow:hidden}}
.hero:after{{content:"";position:absolute;right:-80px;top:-80px;width:260px;height:260px;border:38px solid rgba(143,45,37,.12);border-radius:50%}}
.kicker,.mono{{font-family:"JetBrains Mono",monospace}}
.kicker{{letter-spacing:.18em;text-transform:uppercase;font-size:12px;color:var(--red);font-weight:700;margin-bottom:18px}}
h1{{font-family:"Cormorant Garamond",serif;font-size:58px;line-height:1.02;margin:0 0 14px;color:var(--coal);max-width:920px}}
.subtitle{{font-size:22px;color:var(--muted);max-width:940px;font-style:italic}}
.meta{{display:flex;flex-wrap:wrap;gap:10px;margin-top:26px}}
.chip{{font-family:"JetBrains Mono",monospace;font-size:11px;border:1px solid var(--line);background:#fff;padding:7px 10px;border-radius:6px}}
.chip.hot{{background:var(--bad);border-color:#d9aaa3;color:var(--red);font-weight:700}}
.chip.good{{background:var(--good);border-color:#b8d0b5;color:var(--green);font-weight:700}}
.layout{{display:grid;grid-template-columns:260px 1fr;gap:28px;margin-top:34px}}
.toc{{position:sticky;top:18px;align-self:start;background:rgba(255,250,241,.9);border:1px solid var(--line);padding:18px;border-radius:12px}}
.toc h3{{font-family:"Cormorant Garamond",serif;margin:0 0 12px;font-size:23px;color:var(--red)}}
.toc a{{display:block;border-bottom:1px dotted var(--line);padding:8px 0;color:var(--ink);font-size:14px}}
section{{background:rgba(255,250,241,.86);border:1px solid var(--line);padding:32px;margin-bottom:24px;border-radius:14px}}
h2{{font-family:"Cormorant Garamond",serif;font-size:38px;line-height:1.1;margin:0 0 12px;color:var(--coal)}}
h3{{font-family:"Cormorant Garamond",serif;font-size:27px;margin:26px 0 10px;color:var(--red)}}
.section-note{{color:var(--muted);font-style:italic;margin-bottom:20px}}
.grid{{display:grid;gap:16px}}
.g3{{grid-template-columns:repeat(3,minmax(0,1fr))}}
.g4{{grid-template-columns:repeat(4,minmax(0,1fr))}}
.card{{background:#fff;border:1px solid var(--line2);padding:18px;border-radius:12px}}
.card h4{{margin:0 0 8px;font-family:"JetBrains Mono",monospace;font-size:12px;letter-spacing:.08em;text-transform:uppercase;color:var(--red)}}
.big{{font-family:"Cormorant Garamond",serif;font-size:42px;font-weight:700;color:var(--red);line-height:1}}
.small{{font-size:14px;color:var(--muted)}}
.reading-card{{border-top:4px solid var(--blue)}}
.callout{{border-left:5px solid var(--red);background:#fff;padding:20px 22px;margin:20px 0;border-radius:0 10px 10px 0}}
.callout.good{{border-color:var(--green);background:var(--good)}}
.logic-strip{{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));gap:10px;margin:16px 0}}
.logic-strip div{{background:#fff;border:1px solid var(--line2);border-radius:10px;padding:12px}}
.logic-strip b{{font-family:"JetBrains Mono",monospace;font-size:10px;color:var(--red);text-transform:uppercase;letter-spacing:.08em}}
figure{{margin:26px 0;background:#fff;border:1px solid var(--line);border-radius:14px;padding:16px;box-shadow:0 6px 20px rgba(23,33,47,.06)}}
.fig-title-row{{display:flex;gap:12px;align-items:flex-start;border-bottom:1px solid var(--line2);margin-bottom:14px;padding-bottom:8px}}
.fig-title-row h3{{margin:0;line-height:1.1}}
.fig-id{{font-family:"JetBrains Mono",monospace;font-size:11px;background:var(--coal);color:#fff;border-radius:999px;padding:6px 9px;white-space:nowrap}}
figure img{{display:block;width:100%;height:auto;border-radius:8px;background:#f4eadb;cursor:zoom-in}}
figcaption{{font-size:14px;color:var(--muted);padding:12px 6px 4px;line-height:1.55}}
figcaption strong{{color:var(--ink)}}
.fig-notes{{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:8px;margin-top:10px}}
.fig-notes div{{background:#fff8e8;border:1px solid var(--line2);border-radius:8px;padding:10px;color:var(--ink);font-size:14px;line-height:1.45}}
.fig-notes b{{font-family:"JetBrains Mono",monospace;font-size:10px;text-transform:uppercase;letter-spacing:.08em;color:var(--red);display:block;margin-bottom:4px}}
.fig-notes .limit{{background:#f8ead7}}
.fig-notes .wide-note{{grid-column:1/-1;background:var(--good)}}
.table-wrap{{overflow-x:auto}}
table{{width:100%;border-collapse:collapse;margin:18px 0;background:#fff;font-size:14px;border:1px solid var(--line)}}
th{{background:var(--coal);color:#fff;text-align:left;font-family:"JetBrains Mono",monospace;font-size:11px;letter-spacing:.06em;text-transform:uppercase;padding:10px}}
td{{border-top:1px solid var(--line2);padding:10px;vertical-align:top}}
.footer{{font-family:"JetBrains Mono",monospace;font-size:12px;color:var(--muted);text-align:center;margin-top:34px;padding-top:18px;border-top:1px solid var(--line)}}
.lightbox-overlay{{position:fixed;inset:0;background:rgba(10,16,24,.88);z-index:9999;display:none;align-items:center;justify-content:center;padding:34px;cursor:zoom-out}}
.lightbox-overlay.show{{display:flex}}
.lightbox-overlay img{{max-width:96vw;max-height:88vh;object-fit:contain;background:#fff;border-radius:10px;box-shadow:0 20px 80px rgba(0,0,0,.48);cursor:zoom-out}}
.lightbox-close{{position:fixed;top:18px;right:24px;color:#fff;font-size:42px;line-height:1;font-family:Arial,sans-serif;cursor:pointer;z-index:10000}}
.lightbox-caption{{position:fixed;left:50%;bottom:18px;transform:translateX(-50%);max-width:min(900px,88vw);background:rgba(255,250,241,.96);color:var(--ink);border-radius:10px;padding:10px 14px;font-size:13px;line-height:1.45}}
@media(max-width:900px){{.layout{{grid-template-columns:1fr}}.toc{{position:relative;top:auto}}.g3,.g4,.fig-notes,.logic-strip{{grid-template-columns:1fr}}h1{{font-size:40px}}.hero{{padding:30px 24px}}.wrap{{padding:20px 14px}}.lightbox-overlay{{padding:14px}}.lightbox-caption{{display:none!important}}}}
@media print{{body{{background:#fff}}.toc,.topnav{{display:none}}.layout{{display:block}}.hero,section,figure,.card{{break-inside:avoid;box-shadow:none}}}}
</style>
</head>
<body>
{body}
</body>
</html>
"""


def update_index() -> None:
    idx = HUB / "index.html"
    if not idx.exists():
        return
    text = idx.read_text(encoding="utf-8")
    if "kthyro_ctc_emt_visual_atlas_paper1_style.html" in text:
        return
    needle = '    <a href="kthyro_ctc_emt_paper_flow_master_dossier.html"'
    insert = '    <a href="kthyro_ctc_emt_visual_atlas_paper1_style.html" style="background:linear-gradient(135deg,#fffaf1 0%,#f0dfc4 100%);border:2px solid #8f2d25;color:#102033;font-weight:900;font-size:16px;padding:14px 18px">★ CTC-EMT VISUAL ATLAS · paper1-style readable figures</a>\n'
    pos = text.find(needle)
    if pos >= 0:
        text = text[:pos] + insert + text[pos:]
    else:
        text += "\n" + insert
    idx.write_text(text, encoding="utf-8")
    shutil.copy2(idx, LIVE / "index.html")


def main() -> None:
    html_text = build_html()
    HTML.write_text(html_text, encoding="utf-8")
    LOCAL_HTML.write_text(html_text, encoding="utf-8")
    shutil.copy2(HTML, LIVE_HTML)
    update_index()
    print(f"HTML={HTML}")
    print(f"LIVE={LIVE_HTML}")
    print(f"LOCAL={LOCAL_HTML}")


if __name__ == "__main__":
    main()
