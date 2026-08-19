#!/usr/bin/env python3
from __future__ import annotations

import html
import shutil
from pathlib import Path

import numpy as np
import pandas as pd

import build_paper1_style_visual_atlas as atlas


HTML_NAME = "kthyro_ctc_emt_impact_upgrade.html"
REPORT_NAME = "IMPACT_UPGRADE_PACK_KR.md"
IMPACT_TSV = "impact_upgrade_core_arguments.tsv"
SCENARIO_TSV = "impact_upgrade_outcome_scenarios.tsv"
ONE_LINERS_TSV = "impact_upgrade_reviewer_one_liners.tsv"
FIG_IMPACT = "PF17_samsung_impact_canvas.png"
FIG_SCENARIOS = "PF18_three_outcome_scenario_matrix.png"
ASSET = f"assets/{atlas.ASSET}"


CORE_ARGUMENTS = [
    {
        "impact_axis": "Not another biomarker study",
        "strong_line": "이 proposal의 단위는 biomarker가 아니라 thyroidectomy라는 human perturbation이다.",
        "why_it_matters": "PTC는 예후가 좋아 보이지만, 수술 후 residual-risk biology를 실시간 cellular/genetic state로 읽는 틀이 부족하다.",
        "support": "Breakthrough Case; Representative Overview; PF15; slide02_hypothesis",
        "boundary": "현재 data로 transition을 발견했다고 쓰지 않고, prospective primary endpoint로 둔다.",
    },
    {
        "impact_axis": "Genotype plus state, not mutation alone",
        "strong_line": "BRAF 하나로 환자 상태를 읽는 시대는 끝났다; 혈액 신호도 phenotype-state와 genetic anchor가 같이 필요하다.",
        "why_it_matters": "TCGA public pilot에서 BRAF-mutant PTC도 여러 tissue-state로 갈라진다. 수술 후 blood signal을 mutation-only로 해석하면 biology를 놓친다.",
        "support": "F03; F04; PF03; F08",
        "boundary": "TCGA tissue result를 CTC/recurrence result로 넘기지 않는다.",
    },
    {
        "impact_axis": "Every outcome is informative",
        "strong_line": "성공 시 discovery, 실패 시에도 assay/data boundary를 명확히 남기는 stage-gated science다.",
        "why_it_matters": "CTC collapse, persistence without genetic anchor, persistence with matched genetic signal의 세 시나리오 모두 다음 연구 의사결정을 바꾼다.",
        "support": "PF14; PF18; stage_gate_kill_criteria.tsv",
        "boundary": "실패 가능성을 숨기지 않는다. clinical claim은 validation 전 금지한다.",
    },
    {
        "impact_axis": "Samsung Science fit",
        "strong_line": "넓은 플랫폼이 아니라 좁고 선명한 biological perturbation 질문이다.",
        "why_it_matters": "one disease, one perturbation, one transition, one genetic anchor로 좁혔기 때문에 Science-track reviewer가 공격하기 어렵다.",
        "support": "PF11; slide10_claim_boundary; reviewer defense table",
        "boundary": "AI/ICT는 해석 인프라이지 category identity가 아니다.",
    },
]


SCENARIOS = [
    {
        "scenario": "A. CTC collapses after thyroidectomy",
        "interpretation": "Tumor shedding is surgery-sensitive; postoperative residual cellular signal is low.",
        "publishable_value": "Defines expected perturbation response and negative-control kinetics for PTC CTC-EMT monitoring.",
        "next_step": "Focus on high-risk subset, cfDNA sensitivity, and assay detection limits.",
        "claim_boundary": "Does not prove zero residual disease.",
    },
    {
        "scenario": "B. EM/M CTC persists but no genetic anchor",
        "interpretation": "Phenotype signal may include nonspecific EMT-like cells, assay background, or tumor-unmatched biology.",
        "publishable_value": "Important falsification of phenotype-only CTC interpretation; supports need for matched NGS.",
        "next_step": "Tighten enrichment/QC, add orthogonal marker validation, restrict claim to phenotype feasibility.",
        "claim_boundary": "Do not call this residual tumor biology.",
    },
    {
        "scenario": "C. EM/M CTC persists with tissue/cfDNA genetic match",
        "interpretation": "Strongest residual-risk biology signal: persistent EMT-state cellular phenotype genetically anchored to tumor clone.",
        "publishable_value": "High-impact discovery candidate for postoperative biological monitoring in PTC.",
        "next_step": "Scale cohort, test association with LVI/LNM/ETE/Tg/US/follow-up, prepare validation study.",
        "claim_boundary": "Still not a ready clinical diagnostic until externally validated.",
    },
]


ONE_LINERS = [
    ("One-sentence pitch", "We turn thyroidectomy into a controlled human perturbation experiment to read whether PTC CTC-EMT states collapse, persist, or genetically match residual-risk biology."),
    ("Why this is new", "The novelty is not CTC counting; it is serial CTC-EMT state transition plus matched tumor-normal/cfDNA genetic anchoring around surgery."),
    ("Why public data matters", "Public omics do not prove CTC biology; they prove mutation alone is too weak and provide the marker/genetic map needed for the prospective test."),
    ("Why Science, not ICT", "AI models help interpret longitudinal clone-state transitions, but the category identity is a biological perturbation question."),
    ("Why fundable", "The project is narrow, falsifiable, stage-gated, and capable of producing interpretable biology even when the strongest scenario fails."),
]


def esc(x) -> str:
    return html.escape("" if x is None else str(x), quote=True)


def draw_impact_canvas(path: Path) -> None:
    import cv2

    w, h = 2600, 1500
    img = np.full((h, w, 3), (238, 244, 248), dtype=np.uint8)
    coal = (47, 33, 23)
    red = (37, 45, 143)
    green = (80, 107, 66)
    cream = (241, 250, 255)
    good = (232, 243, 233)
    warn = (212, 241, 255)
    rose = (219, 222, 247)

    def text(line, x, y, scale=0.62, color=coal, thick=2):
        cv2.putText(img, line, (x, y), cv2.FONT_HERSHEY_SIMPLEX, scale, color, thick, cv2.LINE_AA)

    def box(x1, y1, x2, y2, color, header, lines, header_color=red):
        cv2.rectangle(img, (x1, y1), (x2, y2), color, -1)
        cv2.rectangle(img, (x1, y1), (x2, y2), coal, 3)
        text(header, x1 + 24, y1 + 52, 0.82, header_color, 2)
        y = y1 + 96
        for line in lines:
            text(line, x1 + 24, y, 0.60, coal, 2)
            y += 34

    def arrow(p1, p2, color=red, thick=4):
        cv2.arrowedLine(img, p1, p2, color, thick, tipLength=0.04, line_type=cv2.LINE_AA)

    text("Samsung Science Impact Canvas", 70, 95, 1.65, coal, 4)
    text("High impact without overclaim: one perturbation, one transition, one genetic anchor", 75, 145, 0.75, red, 2)

    box(80, 250, 720, 560, cream, "1. Biological perturbation", [
        "Thyroidectomy abruptly reduces tumor burden",
        "Serial blood becomes a time-resolved readout",
        "Question: collapse, persistence, or shift?"
    ])
    box(80, 770, 720, 1080, good, "2. Public-data prior", [
        "BRAF-mutant PTC remains state-heterogeneous",
        "Spatial tissue states are organized",
        "Marker modules have tissue/cell context"
    ])
    box(940, 250, 1660, 1080, warn, "Core Science Case", [
        "PTC CTC-EMT state transition",
        "before and after thyroidectomy",
        "+ matched tumor-normal NGS",
        "+ cfDNA / CTC-enriched NGS subset",
        "+ postoperative residual-risk features",
        "",
        "Primary endpoint:",
        "postoperative CTC E / E-M / M transition"
    ], header_color=green)
    box(1880, 250, 2520, 560, good, "3. Genetic anchor", [
        "Phenotype alone is not enough",
        "Tumor-normal NGS defines clone map",
        "cfDNA/CTC-enriched NGS tests persistence"
    ])
    box(1880, 770, 2520, 1080, rose, "4. Claim discipline", [
        "Not a ready diagnostic",
        "Not public-data CTC proof",
        "Not AI/ICT category identity",
        "Stage-gated and falsifiable"
    ])
    box(430, 1240, 2170, 1390, cream, "Reviewer take-home", [
        "This is not a broad platform. It is a narrow human perturbation experiment that can expose residual-risk biology in PTC."
    ])

    arrow((720, 405), (940, 420))
    arrow((720, 925), (940, 920))
    arrow((1660, 420), (1880, 405))
    arrow((1660, 920), (1880, 925))
    arrow((1300, 1080), (1300, 1235), green, 5)

    text("Best final title: Serial Genetic Mapping of CTC-EMT State Transitions After Thyroidectomy in Papillary Thyroid Cancer", 90, 1460, 0.62, coal, 2)
    cv2.imwrite(str(path), img)


def draw_scenario_matrix(path: Path) -> None:
    import cv2

    w, h = 2600, 1550
    img = np.full((h, w, 3), (238, 244, 248), dtype=np.uint8)
    coal = (47, 33, 23)
    red = (37, 45, 143)
    green = (80, 107, 66)
    cream = (241, 250, 255)
    good = (232, 243, 233)
    warn = (212, 241, 255)
    rose = (219, 222, 247)

    def text(line, x, y, scale=0.58, color=coal, thick=2):
        cv2.putText(img, line, (x, y), cv2.FONT_HERSHEY_SIMPLEX, scale, color, thick, cv2.LINE_AA)

    text("Three Outcome Scenario Matrix", 70, 95, 1.65, coal, 4)
    text("A strong proposal stays valuable under multiple honest biological outcomes", 75, 145, 0.75, red, 2)

    xs = [90, 915, 1740]
    colors = [cream, warn, good]
    headers = ["A. CTC collapses", "B. EM/M persists, no anchor", "C. EM/M persists + genetic match"]
    bodies = [
        ["Expected perturbation response", "Surgery-sensitive tumor shedding", "Defines normal postoperative kinetics", "Boundary: not zero residual disease"],
        ["Phenotype-only signal is unsafe", "Could be nonspecific/background", "Validates need for matched NGS", "Boundary: not residual tumor biology"],
        ["Strongest discovery scenario", "Persistent EMT-state blood signal", "Genetically anchored to tumor clone", "Boundary: needs external validation"],
    ]
    for x, c, head, lines in zip(xs, colors, headers, bodies):
        cv2.rectangle(img, (x, 240), (x + 760, 1180), c, -1)
        cv2.rectangle(img, (x, 240), (x + 760, 1180), coal, 3)
        text(head, x + 28, 305, 0.78, red if head[0] != "C" else green, 2)
        y = 390
        for line in lines:
            text(line, x + 34, y, 0.62, coal, 2)
            y += 56
        cv2.rectangle(img, (x + 40, 720), (x + 720, 1120), (255, 255, 255), -1)
        cv2.rectangle(img, (x + 40, 720), (x + 720, 1120), coal, 2)
        text("Publishable value", x + 65, 780, 0.62, red, 2)
        if head[0] == "A":
            vals = ["negative-control kinetics", "assay sensitivity boundary", "high-risk subset design"]
        elif head[0] == "B":
            vals = ["falsifies phenotype-only claims", "tightens enrichment/QC", "supports genetic anchor necessity"]
        else:
            vals = ["residual-risk biology candidate", "scale-up cohort rationale", "validation-ready hypothesis"]
        yy = 850
        for v in vals:
            text("- " + v, x + 70, yy, 0.54, coal, 2)
            yy += 54

    cv2.rectangle(img, (260, 1280), (2340, 1435), rose, -1)
    cv2.rectangle(img, (260, 1280), (2340, 1435), coal, 3)
    text("Reviewer-safe impact statement", 300, 1340, 0.72, red, 2)
    text("The experiment is high-impact because all three outcomes change how PTC postoperative liquid biopsy should be interpreted.", 300, 1400, 0.62, coal, 2)
    cv2.imwrite(str(path), img)


def table_html(df: pd.DataFrame) -> str:
    out = ["<div class='table-wrap'><table><thead><tr>"]
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


def html_page(args_df: pd.DataFrame, scenario_df: pd.DataFrame, line_df: pd.DataFrame) -> str:
    return f"""<!DOCTYPE html>
<html lang="ko">
<head>
<meta charset="UTF-8" />
<meta name="viewport" content="width=device-width, initial-scale=1.0" />
<title>K-Thyro CTC-EMT Impact Upgrade</title>
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
  h1{{font-family:"Cormorant Garamond",serif;font-size:62px;line-height:1.0;margin:0 0 14px;color:var(--coal);max-width:1040px}} .subtitle{{font-size:23px;color:var(--muted);max-width:1020px;font-style:italic}}
  .thesis{{margin-top:24px;background:#fff;border:1px solid var(--line);border-left:7px solid var(--red);border-radius:0 12px 12px 0;padding:22px;font-size:22px}} .thesis b{{color:var(--red)}}
  .layout{{display:grid;grid-template-columns:250px 1fr;gap:26px;margin-top:30px}} .toc{{position:sticky;top:18px;align-self:start;background:rgba(255,250,241,.92);border:1px solid var(--line);border-radius:12px;padding:16px}} .toc h3{{font-family:"Cormorant Garamond",serif;color:var(--red);font-size:24px;margin:0 0 10px}} .toc a{{display:block;color:var(--ink);border-bottom:1px dotted var(--line);padding:8px 0;font-size:14px}}
  section{{background:rgba(255,250,241,.9);border:1px solid var(--line);border-radius:14px;padding:28px;margin-bottom:24px}} h2{{font-family:"Cormorant Garamond",serif;font-size:42px;line-height:1.05;margin:0 0 12px;color:var(--coal)}} .section-note{{color:var(--muted);font-style:italic;margin:0 0 18px}}
  .big-fig{{display:block;background:#fff;border:1px solid var(--line);border-radius:14px;padding:14px;margin:14px 0 20px}} .big-fig img{{display:block;width:100%;height:auto;border-radius:8px;background:#f4eadb}}
  .impact-grid{{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:14px}} .impact-card{{background:#fff;border:1px solid var(--line2);border-radius:12px;padding:16px}} .impact-card b{{display:block;font-family:"JetBrains Mono",monospace;font-size:10px;color:#fff;background:var(--red);width:max-content;border-radius:999px;padding:5px 8px;text-transform:uppercase;letter-spacing:.08em}} .impact-card h3{{font-family:"Cormorant Garamond",serif;font-size:30px;line-height:1.05;margin:13px 0 8px;color:var(--coal)}} .impact-card p{{margin:0 0 8px}} .boundary{{background:#f8ead7;border:1px solid var(--line2);border-radius:8px;padding:10px;font-size:14px;color:var(--ink)}}
  .table-wrap{{overflow-x:auto}} table{{width:100%;border-collapse:collapse;background:#fff;border:1px solid var(--line);font-size:14px}} th{{background:var(--coal);color:#fff;text-align:left;font-family:"JetBrains Mono",monospace;font-size:11px;letter-spacing:.06em;text-transform:uppercase;padding:10px}} td{{border-top:1px solid var(--line2);padding:10px;vertical-align:top}}
  .link-grid{{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));gap:12px}} .link-card{{display:block;background:#fff;border:1px solid var(--line2);border-radius:12px;padding:16px}} .link-card b{{display:block;font-family:"Cormorant Garamond",serif;font-size:27px;color:var(--red)}} .link-card span{{display:block;color:var(--muted)}}
  .footer{{font-family:"JetBrains Mono",monospace;font-size:12px;color:var(--muted);text-align:center;margin-top:34px;padding-top:18px;border-top:1px solid var(--line)}}
  @media(max-width:1000px){{.layout{{grid-template-columns:1fr}}.toc{{position:relative;top:auto}}.impact-grid,.link-grid{{grid-template-columns:1fr}}h1{{font-size:42px}}.hero{{padding:30px 22px}}}}
</style>
</head>
<body>
<div class="wrap">
  <nav class="top">
    <a href="kthyro_ctc_emt_all_pages_launchpad.html">All Pages</a>
    <a href="kthyro_ctc_emt_paper_ready_figure_table_pack.html">Paper-Ready Pack</a>
    <a href="kthyro_ctc_emt_representative_overview.html">Representative Overview</a>
    <a href="kthyro_ctc_emt_breakthrough_case.html">Breakthrough Case</a>
  </nav>
  <header class="hero">
    <div class="kicker">K-Thyro CTC-EMT · impact upgrade · high impact, no overclaim</div>
    <h1>Impact Upgrade Pack</h1>
    <p class="subtitle">가능성 문구를 줄이고, 삼성 Science reviewer가 바로 보는 impact canvas, 세 가지 outcome scenario, reviewer one-liner로 제안서의 설득력을 올린다.</p>
    <div class="thesis"><b>Impact thesis.</b> The project is fundable because it converts thyroidectomy into a narrow, falsifiable, serial liquid-biopsy perturbation experiment with a required genetic anchor.</div>
  </header>
  <div class="layout">
    <aside class="toc">
      <h3>Read Order</h3>
      <a href="#canvas">Impact Canvas</a>
      <a href="#arguments">Core Arguments</a>
      <a href="#scenarios">Outcome Scenarios</a>
      <a href="#liners">Reviewer One-Liners</a>
      <a href="#downloads">Downloads</a>
    </aside>
    <main>
      <section id="canvas">
        <h2>Impact Canvas</h2>
        <p class="section-note">한 장으로 reviewer에게 보여줄 high-impact logic.</p>
        <a class="big-fig" href="{ASSET}/{FIG_IMPACT}"><img src="{ASSET}/{FIG_IMPACT}" alt="Samsung impact canvas" /></a>
      </section>
      <section id="arguments">
        <h2>Core Arguments</h2>
        <p class="section-note">대박처럼 보이게 하는 문장과 그 boundary를 같이 둔다.</p>
        <div class="impact-grid">
        {''.join(f'<article class="impact-card"><b>{esc(r["impact_axis"])}</b><h3>{esc(r["strong_line"])}</h3><p>{esc(r["why_it_matters"])}</p><p><b style="background:var(--coal)">Support</b> {esc(r["support"])}</p><div class="boundary">{esc(r["boundary"])}</div></article>' for _, r in args_df.iterrows())}
        </div>
      </section>
      <section id="scenarios">
        <h2>Three Outcome Scenario Matrix</h2>
        <p class="section-note">가장 강한 점: 세 가지 결과 모두 해석 가능하고 다음 의사결정을 바꾼다.</p>
        <a class="big-fig" href="{ASSET}/{FIG_SCENARIOS}"><img src="{ASSET}/{FIG_SCENARIOS}" alt="Three outcome scenario matrix" /></a>
        {table_html(scenario_df)}
      </section>
      <section id="liners">
        <h2>Reviewer One-Liners</h2>
        <p class="section-note">회의/카톡/슬라이드에 바로 넣을 수 있는 문장.</p>
        {table_html(line_df)}
      </section>
      <section id="downloads">
        <h2>Downloads</h2>
        <div class="link-grid">
          <a class="link-card" href="{IMPACT_TSV}"><b>Arguments TSV</b><span>impact argument table</span></a>
          <a class="link-card" href="{SCENARIO_TSV}"><b>Scenario TSV</b><span>outcome matrix</span></a>
          <a class="link-card" href="{ONE_LINERS_TSV}"><b>One-liners TSV</b><span>reviewer lines</span></a>
          <a class="link-card" href="{REPORT_NAME}"><b>Guide MD</b><span>impact handoff</span></a>
        </div>
      </section>
    </main>
  </div>
  <div class="footer">Impact upgraded by sharpening the Science question and the falsifiable outcome logic, not by overclaiming current data.</div>
</div>
</body>
</html>"""


def report_md(args_df: pd.DataFrame, scenario_df: pd.DataFrame, line_df: pd.DataFrame) -> str:
    lines = [
        "# Impact Upgrade Pack",
        "",
        "Purpose: 가능성 문구를 줄이고, 삼성 Science reviewer가 바로 이해할 수 있는 impact canvas, outcome scenarios, one-liners로 설득력을 높인다.",
        "",
        "Impact thesis: thyroidectomy is a narrow human perturbation; serial CTC-EMT plus matched NGS makes it a falsifiable Science question.",
        "",
        "Boundary: high impact does not mean current public data proves postoperative CTC transition.",
        "",
        "## Core arguments",
    ]
    for _, r in args_df.iterrows():
        lines.append(f"- {r['impact_axis']}: {r['strong_line']} Support: {r['support']} Boundary: {r['boundary']}")
    lines += ["", "## Outcome scenarios"]
    for _, r in scenario_df.iterrows():
        lines.append(f"- {r['scenario']}: {r['interpretation']} Value: {r['publishable_value']} Boundary: {r['claim_boundary']}")
    lines += ["", "## Reviewer one-liners"]
    for _, r in line_df.iterrows():
        lines.append(f"- {r['context']}: {r['one_liner']}")
    return "\n".join(lines)


def copy_table(df: pd.DataFrame, name: str) -> None:
    for p in [atlas.TAB / name, atlas.HUB / name]:
        df.to_csv(p, sep="\t", index=False)
    shutil.copy2(atlas.HUB / name, atlas.LIVE / name)


def update_index_handoff_launchpad() -> None:
    idx = atlas.HUB / "index.html"
    if idx.exists():
        text = idx.read_text(encoding="utf-8")
        if HTML_NAME not in text:
            needle = '    <a href="kthyro_ctc_emt_paper_ready_figure_table_pack.html"'
            insert = f'    <a href="{HTML_NAME}" style="background:linear-gradient(135deg,#8f2d25 0%,#fffaf1 100%);border:3px solid #102033;color:#102033;font-weight:1000;font-size:18px;padding:16px 20px">★★ CTC-EMT IMPACT UPGRADE · Samsung reviewer canvas</a>\n'
            pos = text.find(needle)
            text = text[:pos] + insert + text[pos:] if pos >= 0 else text + "\n" + insert
            idx.write_text(text, encoding="utf-8")
        shutil.copy2(idx, atlas.LIVE / "index.html")

    for path in [atlas.REP / "VISUAL_WEB_HANDOFF_KR.md", atlas.HUB / "VISUAL_WEB_HANDOFF_KR.md"]:
        if not path.exists():
            continue
        text = path.read_text(encoding="utf-8")
        line = "- Impact Upgrade Pack: http://40.82.129.113:8012/papers_hub_2026_05_04/kthyro_ctc_emt_impact_upgrade.html\n"
        if HTML_NAME not in text:
            text = text.replace("## 최우선 링크\n\n", "## 최우선 링크\n\n" + line)
            path.write_text(text, encoding="utf-8")
    handoff = atlas.HUB / "VISUAL_WEB_HANDOFF_KR.md"
    if handoff.exists():
        shutil.copy2(handoff, atlas.LIVE / "VISUAL_WEB_HANDOFF_KR.md")

    launch = atlas.HUB / "kthyro_ctc_emt_all_pages_launchpad.html"
    if launch.exists() and HTML_NAME not in launch.read_text(encoding="utf-8"):
        text = launch.read_text(encoding="utf-8")
        insert = f'<a href="{HTML_NAME}"><b>임팩트 업</b><span>Reviewer Impact Canvas</span></a>'
        text = text.replace('<div class="route-grid">', '<div class="route-grid">' + insert, 1)
        launch.write_text(text, encoding="utf-8")
        shutil.copy2(launch, atlas.LIVE / "kthyro_ctc_emt_all_pages_launchpad.html")


def main() -> None:
    fig_dir = atlas.OUT / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    asset_dir = atlas.HUB / "assets" / atlas.ASSET
    live_asset_dir = atlas.LIVE / "assets" / atlas.ASSET

    impact_path = fig_dir / FIG_IMPACT
    scenario_path = fig_dir / FIG_SCENARIOS
    draw_impact_canvas(impact_path)
    draw_scenario_matrix(scenario_path)
    for p in [impact_path, scenario_path]:
        shutil.copy2(p, asset_dir / p.name)
        shutil.copy2(p, live_asset_dir / p.name)

    args_df = pd.DataFrame(CORE_ARGUMENTS)
    scenario_df = pd.DataFrame(SCENARIOS)
    line_df = pd.DataFrame(ONE_LINERS, columns=["context", "one_liner"])
    copy_table(args_df, IMPACT_TSV)
    copy_table(scenario_df, SCENARIO_TSV)
    copy_table(line_df, ONE_LINERS_TSV)

    html_text = html_page(args_df, scenario_df, line_df)
    src = atlas.HUB / HTML_NAME
    live = atlas.LIVE / HTML_NAME
    local = atlas.REP / HTML_NAME
    src.write_text(html_text, encoding="utf-8")
    local.write_text(html_text, encoding="utf-8")
    shutil.copy2(src, live)

    md = report_md(args_df, scenario_df, line_df)
    report = atlas.REP / REPORT_NAME
    hub_report = atlas.HUB / REPORT_NAME
    report.write_text(md, encoding="utf-8")
    hub_report.write_text(md, encoding="utf-8")
    shutil.copy2(hub_report, atlas.LIVE / REPORT_NAME)

    update_index_handoff_launchpad()
    print(f"HTML={src}")
    print(f"LIVE={live}")
    print(f"FIGURES={impact_path}, {scenario_path}")
    print(f"REPORT={report}")
    print(f"TABLES={atlas.TAB / IMPACT_TSV}, {atlas.TAB / SCENARIO_TSV}, {atlas.TAB / ONE_LINERS_TSV}")


if __name__ == "__main__":
    main()
