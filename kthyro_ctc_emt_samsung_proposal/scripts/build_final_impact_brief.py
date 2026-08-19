#!/usr/bin/env python3
from __future__ import annotations

import html
import shutil
from pathlib import Path

import numpy as np
import pandas as pd

import build_paper1_style_visual_atlas as atlas


HTML_NAME = "kthyro_ctc_emt_final_impact_brief.html"
REPORT_NAME = "FINAL_IMPACT_BRIEF_KR.md"
OPENING_TSV = "final_impact_opening_lines.tsv"
DECISION_TSV = "final_impact_reviewer_decision_scorecard.tsv"
FATAL_TSV = "final_impact_fatal_flaw_closure.tsv"
FIG_PYRAMID = "PF19_final_funding_logic_pyramid.png"
FIG_SCORECARD = "PF20_reviewer_decision_scorecard.png"
ASSET = f"assets/{atlas.ASSET}"


OPENING_LINES = [
    ("Opening thesis", "이 과제는 CTC를 한 번 세는 연구가 아니라, 갑상선절제술이라는 인간 내 perturbation을 이용해 PTC의 circulating EMT state가 어떻게 무너지고, 남고, 유전적으로 연결되는지를 읽는 Science 과제입니다."),
    ("Why Samsung Science", "질문이 좁습니다: one disease, one perturbation, one state transition, one genetic anchor입니다."),
    ("Why now", "Yu 2024가 serial PTC CTC-EMT 측정 가능성을 열었지만, matched tissue-normal NGS/cfDNA/CTC-enriched NGS와 연결된 postoperative genetic map은 비어 있습니다."),
    ("Why public pilot matters", "우리 public pilot은 CTC biology를 증명하지 않습니다. 대신 BRAF/mutation alone이 부족하고 phenotype-state와 genotype을 함께 읽어야 한다는 사전 지도를 제공합니다."),
    ("Most fundable sentence", "수술 후 EM/M CTC가 사라지는지, 남지만 유전적으로 맞지 않는지, 혹은 tumor clone과 맞아 residual-risk biology 후보가 되는지를 stage-gated로 검정하겠습니다."),
    ("Boundary sentence", "현재는 clinical diagnostic이나 recurrence predictor가 아니라, postoperative biological state transition을 검정하는 research-grade perturbation study입니다."),
]


DECISION_ROWS = [
    {
        "review_axis": "Science category fit",
        "score": "High",
        "why_reviewer_cares": "Narrow biological perturbation, not broad ICT/platform.",
        "support": "PF19; PF11; slide10",
        "remaining_risk": "Reviewer may still see many assays unless primary endpoint is repeated clearly.",
        "tight_answer": "Primary endpoint is one thing: postoperative CTC E/EM/M transition.",
    },
    {
        "review_axis": "Novelty",
        "score": "High if matched NGS is secured",
        "why_reviewer_cares": "CTC count alone is not enough; serial EMT state plus genetic map is the novelty.",
        "support": "PF03; F08; PF20",
        "remaining_risk": "CTC-enriched NGS may fail in early PTC.",
        "tight_answer": "Stage-gated subset analysis; matched tumor-normal and cfDNA remain genetic anchors.",
    },
    {
        "review_axis": "Feasibility",
        "score": "Moderate-high",
        "why_reviewer_cares": "Yu 2024 supports serial CTC-EMT feasibility; hospital data access is decisive.",
        "support": "Slide 3; CRF; future data unlock map",
        "remaining_risk": "No local raw Yu dataset yet.",
        "tight_answer": "Proposal asks exactly for the missing serial phenotype/genetic/follow-up fields.",
    },
    {
        "review_axis": "Impact",
        "score": "High",
        "why_reviewer_cares": "A positive genetic-match scenario opens residual-risk biology; negative scenarios still constrain liquid-biopsy interpretation.",
        "support": "PF18; impact scenario table",
        "remaining_risk": "Clinical utility requires later validation.",
        "tight_answer": "This phase is biological monitoring discovery, not clinical deployment.",
    },
    {
        "review_axis": "Overclaim control",
        "score": "Strong",
        "why_reviewer_cares": "Reviewer trust increases when forbidden claims are explicit.",
        "support": "PF05; PF14; safe sentence bank",
        "remaining_risk": "Slides may drift into diagnostic language.",
        "tight_answer": "Every figure carries claim boundary and kill criteria.",
    },
]


FATAL_FLAWS = [
    ("No local CTC data yet", "Do not pretend we have it", "Use Yu 2024 as feasibility; make hospital serial data the explicit first unlock", "CRF, data unlock map"),
    ("CTC phenotype can be nonspecific", "Do not rely on phenotype alone", "Require matched tumor-normal NGS/cfDNA/CTC-enriched subset anchoring", "PF03, F08"),
    ("Early PTC may have weak cfDNA signal", "Do not guarantee success", "Stage-gate cfDNA to high-risk/high-signal subset and predefine QC failure downgrade", "PF14, kill criteria"),
    ("PTC recurrence is slow", "Do not promise early recurrence prediction", "Use biological state transition as primary endpoint; clinical association as secondary/exploratory", "SAP, claim boundary"),
    ("Too broad / platform-like", "Remove spatial/drug/AI identity", "Keep one disease, one surgery, one transition, one genetic map", "PF19, slide10"),
]


def esc(x) -> str:
    return html.escape("" if x is None else str(x), quote=True)


def draw_pyramid(path: Path) -> None:
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

    def text(s, x, y, scale=0.6, color=coal, thick=2):
        cv2.putText(img, s, (x, y), cv2.FONT_HERSHEY_SIMPLEX, scale, color, thick, cv2.LINE_AA)

    text("Final Funding Logic Pyramid", 70, 95, 1.75, coal, 4)
    text("How this stops being a possibility and becomes a fundable Samsung Science question", 75, 145, 0.75, red, 2)

    levels = [
        (410, 1200, 2190, 1370, cream, "Clinical gap", ["Postoperative residual-risk biology is not directly read in real time"]),
        (520, 1010, 2080, 1170, good, "Perturbation", ["Thyroidectomy abruptly changes tumor burden and shedding context"]),
        (650, 820, 1950, 980, warn, "Primary readout", ["Serial CTC epithelial / hybrid E-M / mesenchymal state transition"]),
        (790, 630, 1810, 790, good, "Genetic anchor", ["Matched tumor-normal NGS + cfDNA + CTC-enriched subset where feasible"]),
        (950, 440, 1650, 600, rose, "Breakthrough test", ["Collapse vs persistence vs genetic match after surgery"]),
        (1100, 250, 1500, 400, cream, "Fundable Science", ["Residual-risk biology in PTC"]),
    ]
    for x1, y1, x2, y2, color, head, lines in levels:
        cv2.rectangle(img, (x1, y1), (x2, y2), color, -1)
        cv2.rectangle(img, (x1, y1), (x2, y2), coal, 3)
        text(head, x1 + 24, y1 + 48, 0.72, red if head != "Fundable Science" else green, 2)
        y = y1 + 96
        for line in lines:
            text(line, x1 + 24, y, 0.58, coal, 2)
            y += 34

    cv2.rectangle(img, (85, 360), (670, 700), (255, 255, 255), -1)
    cv2.rectangle(img, (85, 360), (670, 700), coal, 3)
    text("Why reviewer says yes", 115, 420, 0.72, red, 2)
    for i, line in enumerate(["narrow question", "clear perturbation", "published feasibility", "genetic anchoring", "stage-gated risk"]):
        text("- " + line, 125, 490 + i * 38, 0.56, coal, 2)

    cv2.rectangle(img, (1930, 360), (2515, 700), (255, 255, 255), -1)
    cv2.rectangle(img, (1930, 360), (2515, 700), coal, 3)
    text("Why reviewer trusts it", 1960, 420, 0.72, red, 2)
    for i, line in enumerate(["no recurrence overclaim", "no public-data CTC claim", "no AI category drift", "kill criteria explicit", "clinical validation later"]):
        text("- " + line, 1970, 490 + i * 38, 0.56, coal, 2)

    text("One-line ask: Fund the first serial genetic map of PTC CTC-EMT state transition after thyroidectomy.", 135, 1490, 0.72, coal, 2)
    cv2.imwrite(str(path), img)


def draw_scorecard(path: Path) -> None:
    import cv2

    w, h = 2600, 1500
    img = np.full((h, w, 3), (238, 244, 248), dtype=np.uint8)
    coal = (47, 33, 23)
    red = (37, 45, 143)
    green = (80, 107, 66)
    cream = (241, 250, 255)
    good = (232, 243, 233)
    warn = (212, 241, 255)

    def text(s, x, y, scale=0.58, color=coal, thick=2):
        cv2.putText(img, s, (x, y), cv2.FONT_HERSHEY_SIMPLEX, scale, color, thick, cv2.LINE_AA)

    text("Samsung Reviewer Decision Scorecard", 70, 95, 1.65, coal, 4)
    text("Fundable when impact is high and claim boundaries are explicit", 75, 145, 0.74, red, 2)

    x0, y0 = 80, 230
    widths = [400, 300, 680, 520, 520]
    headers = ["Review axis", "Score", "Why it matters", "Remaining risk", "Tight answer"]
    x = x0
    for head, ww in zip(headers, widths):
        cv2.rectangle(img, (x, y0), (x + ww, y0 + 70), coal, -1)
        text(head, x + 18, y0 + 45, 0.55, (255, 255, 255), 2)
        x += ww

    y = y0 + 70
    for idx, row in enumerate(DECISION_ROWS):
        color = good if idx % 2 == 0 else cream
        x = x0
        vals = [
            row["review_axis"],
            row["score"],
            row["why_reviewer_cares"],
            row["remaining_risk"],
            row["tight_answer"],
        ]
        for val, ww in zip(vals, widths):
            cv2.rectangle(img, (x, y), (x + ww, y + 190), color, -1)
            cv2.rectangle(img, (x, y), (x + ww, y + 190), coal, 2)
            words = str(val).split()
            line = ""
            yy = y + 42
            max_chars = max(14, int(ww / 15))
            for word in words:
                if len(line) + len(word) + 1 > max_chars:
                    text(line, x + 14, yy, 0.46, red if ww == widths[1] else coal, 1)
                    yy += 30
                    line = word
                else:
                    line = (line + " " + word).strip()
            if line:
                text(line, x + 14, yy, 0.46, red if ww == widths[1] else coal, 1)
            x += ww
        y += 190

    cv2.rectangle(img, (260, 1285), (2340, 1425), warn, -1)
    cv2.rectangle(img, (260, 1285), (2340, 1425), coal, 3)
    text("Decision line", 300, 1340, 0.72, red, 2)
    text("Fund if the proposal stays narrow: thyroidectomy perturbation + serial CTC-EMT + matched genetic map.", 300, 1400, 0.62, coal, 2)
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


def html_page(opening_df: pd.DataFrame, decision_df: pd.DataFrame, fatal_df: pd.DataFrame) -> str:
    return f"""<!DOCTYPE html>
<html lang="ko">
<head>
<meta charset="UTF-8" />
<meta name="viewport" content="width=device-width, initial-scale=1.0" />
<title>K-Thyro CTC-EMT Final Impact Brief</title>
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
  .ask{{margin-top:24px;background:#fff;border:1px solid var(--line);border-left:7px solid var(--green);border-radius:0 12px 12px 0;padding:22px;font-size:24px}} .ask b{{color:var(--green)}}
  .layout{{display:grid;grid-template-columns:250px 1fr;gap:26px;margin-top:30px}} .toc{{position:sticky;top:18px;align-self:start;background:rgba(255,250,241,.92);border:1px solid var(--line);border-radius:12px;padding:16px}} .toc h3{{font-family:"Cormorant Garamond",serif;color:var(--red);font-size:24px;margin:0 0 10px}} .toc a{{display:block;color:var(--ink);border-bottom:1px dotted var(--line);padding:8px 0;font-size:14px}}
  section{{background:rgba(255,250,241,.9);border:1px solid var(--line);border-radius:14px;padding:28px;margin-bottom:24px}} h2{{font-family:"Cormorant Garamond",serif;font-size:42px;line-height:1.05;margin:0 0 12px;color:var(--coal)}} .section-note{{color:var(--muted);font-style:italic;margin:0 0 18px}}
  .big-fig{{display:block;background:#fff;border:1px solid var(--line);border-radius:14px;padding:14px;margin:14px 0 20px}} .big-fig img{{display:block;width:100%;height:auto;border-radius:8px;background:#f4eadb}}
  .line-grid{{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:12px}} .line-card{{background:#fff;border:1px solid var(--line2);border-radius:12px;padding:16px}} .line-card b{{display:block;font-family:"JetBrains Mono",monospace;font-size:10px;background:var(--coal);color:#fff;width:max-content;border-radius:999px;padding:5px 8px;text-transform:uppercase;letter-spacing:.08em}} .line-card p{{font-size:18px;margin:12px 0 0}}
  .table-wrap{{overflow-x:auto}} table{{width:100%;border-collapse:collapse;background:#fff;border:1px solid var(--line);font-size:14px}} th{{background:var(--coal);color:#fff;text-align:left;font-family:"JetBrains Mono",monospace;font-size:11px;letter-spacing:.06em;text-transform:uppercase;padding:10px}} td{{border-top:1px solid var(--line2);padding:10px;vertical-align:top}}
  .link-grid{{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));gap:12px}} .link-card{{display:block;background:#fff;border:1px solid var(--line2);border-radius:12px;padding:16px}} .link-card b{{display:block;font-family:"Cormorant Garamond",serif;font-size:27px;color:var(--red)}} .link-card span{{display:block;color:var(--muted)}}
  .footer{{font-family:"JetBrains Mono",monospace;font-size:12px;color:var(--muted);text-align:center;margin-top:34px;padding-top:18px;border-top:1px solid var(--line)}}
  @media(max-width:1000px){{.layout{{grid-template-columns:1fr}}.toc{{position:relative;top:auto}}.line-grid,.link-grid{{grid-template-columns:1fr}}h1{{font-size:42px}}.hero{{padding:30px 22px}}}}
</style>
</head>
<body>
<div class="wrap">
  <nav class="top">
    <a href="kthyro_ctc_emt_all_pages_launchpad.html">All Pages</a>
    <a href="kthyro_ctc_emt_impact_upgrade.html">Impact Upgrade</a>
    <a href="kthyro_ctc_emt_paper_ready_figure_table_pack.html">Paper-Ready Pack</a>
    <a href="kthyro_ctc_emt_representative_overview.html">Representative Overview</a>
  </nav>
  <header class="hero">
    <div class="kicker">K-Thyro CTC-EMT · final impact brief · first 30 seconds</div>
    <h1>Final Impact Brief</h1>
    <p class="subtitle">교수님이나 삼성 reviewer가 처음 열었을 때 바로 fundable Science case로 읽히게 하는 최종 impact opening. 강하게 말하되, 현재 데이터의 경계를 더 명확히 한다.</p>
    <div class="ask"><b>One-line ask.</b> Fund the first serial genetic map of PTC CTC-EMT state transition after thyroidectomy.</div>
  </header>
  <div class="layout">
    <aside class="toc">
      <h3>Read Order</h3>
      <a href="#pyramid">Funding Pyramid</a>
      <a href="#scorecard">Reviewer Scorecard</a>
      <a href="#opening">Opening Lines</a>
      <a href="#fatal">Fatal Flaw Closure</a>
      <a href="#downloads">Downloads</a>
    </aside>
    <main>
      <section id="pyramid">
        <h2>Funding Logic Pyramid</h2>
        <p class="section-note">왜 이 과제가 가능성 단계가 아니라 fundable Science 질문인지 한 장으로 압축.</p>
        <a class="big-fig" href="{ASSET}/{FIG_PYRAMID}"><img src="{ASSET}/{FIG_PYRAMID}" alt="Final funding logic pyramid" /></a>
      </section>
      <section id="scorecard">
        <h2>Reviewer Decision Scorecard</h2>
        <p class="section-note">리뷰어가 보는 축별로 강점, 남은 리스크, 한 줄 답변을 정리.</p>
        <a class="big-fig" href="{ASSET}/{FIG_SCORECARD}"><img src="{ASSET}/{FIG_SCORECARD}" alt="Reviewer decision scorecard" /></a>
        {table_html(decision_df)}
      </section>
      <section id="opening">
        <h2>Opening Lines</h2>
        <p class="section-note">제안서 첫 문단, 발표 첫 30초, 교수님 카톡에 바로 쓰는 문장.</p>
        <div class="line-grid">
          {''.join(f'<article class="line-card"><b>{esc(r["context"])}</b><p>{esc(r["line"])}</p></article>' for _, r in opening_df.iterrows())}
        </div>
      </section>
      <section id="fatal">
        <h2>Fatal Flaw Closure</h2>
        <p class="section-note">심사위원의 치명타 질문을 숨기지 않고 먼저 닫는다.</p>
        {table_html(fatal_df)}
      </section>
      <section id="downloads">
        <h2>Downloads</h2>
        <div class="link-grid">
          <a class="link-card" href="{OPENING_TSV}"><b>Opening TSV</b><span>first 30-sec lines</span></a>
          <a class="link-card" href="{DECISION_TSV}"><b>Scorecard TSV</b><span>reviewer decision scorecard</span></a>
          <a class="link-card" href="{FATAL_TSV}"><b>Fatal Flaw TSV</b><span>attack closure table</span></a>
          <a class="link-card" href="{REPORT_NAME}"><b>Guide MD</b><span>final impact memo</span></a>
        </div>
      </section>
    </main>
  </div>
  <div class="footer">Final impact brief: stronger first impression by narrowing the science and exposing the risk controls.</div>
</div>
</body>
</html>"""


def report_md(opening_df: pd.DataFrame, decision_df: pd.DataFrame, fatal_df: pd.DataFrame) -> str:
    lines = [
        "# Final Impact Brief",
        "",
        "Purpose: 첫 30초에 fundable Science case로 읽히게 하는 최종 impact opening.",
        "",
        "One-line ask: Fund the first serial genetic map of PTC CTC-EMT state transition after thyroidectomy.",
        "",
        "Boundary: current public data builds the prior map; prospective Yu/hospital data are required for the actual postoperative CTC-EMT transition test.",
        "",
        "## Opening lines",
    ]
    for _, r in opening_df.iterrows():
        lines.append(f"- {r['context']}: {r['line']}")
    lines += ["", "## Reviewer scorecard"]
    for _, r in decision_df.iterrows():
        lines.append(f"- {r['review_axis']}: {r['score']}. Risk: {r['remaining_risk']} Answer: {r['tight_answer']}")
    lines += ["", "## Fatal flaw closure"]
    for _, r in fatal_df.iterrows():
        lines.append(f"- {r['fatal_flaw']}: {r['closure_strategy']} Support: {r['support']}")
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
            needle = '    <a href="kthyro_ctc_emt_impact_upgrade.html"'
            insert = f'    <a href="{HTML_NAME}" style="background:linear-gradient(135deg,#102033 0%,#f0dfc4 100%);border:3px solid #8f2d25;color:#102033;font-weight:1000;font-size:18px;padding:16px 20px">★★ CTC-EMT FINAL IMPACT BRIEF · first 30 seconds</a>\n'
            pos = text.find(needle)
            text = text[:pos] + insert + text[pos:] if pos >= 0 else text + "\n" + insert
            idx.write_text(text, encoding="utf-8")
        shutil.copy2(idx, atlas.LIVE / "index.html")

    for path in [atlas.REP / "VISUAL_WEB_HANDOFF_KR.md", atlas.HUB / "VISUAL_WEB_HANDOFF_KR.md"]:
        if not path.exists():
            continue
        text = path.read_text(encoding="utf-8")
        line = "- Final Impact Brief: http://40.82.129.113:8012/papers_hub_2026_05_04/kthyro_ctc_emt_final_impact_brief.html\n"
        if HTML_NAME not in text:
            text = text.replace("## 최우선 링크\n\n", "## 최우선 링크\n\n" + line)
            path.write_text(text, encoding="utf-8")
    handoff = atlas.HUB / "VISUAL_WEB_HANDOFF_KR.md"
    if handoff.exists():
        shutil.copy2(handoff, atlas.LIVE / "VISUAL_WEB_HANDOFF_KR.md")

    launch = atlas.HUB / "kthyro_ctc_emt_all_pages_launchpad.html"
    if launch.exists() and HTML_NAME not in launch.read_text(encoding="utf-8"):
        text = launch.read_text(encoding="utf-8")
        insert = f'<a href="{HTML_NAME}"><b>최종 임팩트</b><span>Final 30-sec Brief</span></a>'
        text = text.replace('<div class="route-grid">', '<div class="route-grid">' + insert, 1)
        launch.write_text(text, encoding="utf-8")
        shutil.copy2(launch, atlas.LIVE / "kthyro_ctc_emt_all_pages_launchpad.html")


def main() -> None:
    fig_dir = atlas.OUT / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    asset_dir = atlas.HUB / "assets" / atlas.ASSET
    live_asset_dir = atlas.LIVE / "assets" / atlas.ASSET

    pyramid = fig_dir / FIG_PYRAMID
    scorecard = fig_dir / FIG_SCORECARD
    draw_pyramid(pyramid)
    draw_scorecard(scorecard)
    for p in [pyramid, scorecard]:
        shutil.copy2(p, asset_dir / p.name)
        shutil.copy2(p, live_asset_dir / p.name)

    opening_df = pd.DataFrame(OPENING_LINES, columns=["context", "line"])
    decision_df = pd.DataFrame(DECISION_ROWS)
    fatal_df = pd.DataFrame(FATAL_FLAWS, columns=["fatal_flaw", "wrong_move", "closure_strategy", "support"])
    copy_table(opening_df, OPENING_TSV)
    copy_table(decision_df, DECISION_TSV)
    copy_table(fatal_df, FATAL_TSV)

    html_text = html_page(opening_df, decision_df, fatal_df)
    src = atlas.HUB / HTML_NAME
    live = atlas.LIVE / HTML_NAME
    local = atlas.REP / HTML_NAME
    src.write_text(html_text, encoding="utf-8")
    local.write_text(html_text, encoding="utf-8")
    shutil.copy2(src, live)

    md = report_md(opening_df, decision_df, fatal_df)
    report = atlas.REP / REPORT_NAME
    hub_report = atlas.HUB / REPORT_NAME
    report.write_text(md, encoding="utf-8")
    hub_report.write_text(md, encoding="utf-8")
    shutil.copy2(hub_report, atlas.LIVE / REPORT_NAME)

    update_index_handoff_launchpad()
    print(f"HTML={src}")
    print(f"LIVE={live}")
    print(f"FIGURES={pyramid}, {scorecard}")
    print(f"REPORT={report}")
    print(f"TABLES={atlas.TAB / OPENING_TSV}, {atlas.TAB / DECISION_TSV}, {atlas.TAB / FATAL_TSV}")


if __name__ == "__main__":
    main()
