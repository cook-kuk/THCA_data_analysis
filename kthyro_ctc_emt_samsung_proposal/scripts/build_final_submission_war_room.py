#!/usr/bin/env python3
from __future__ import annotations

import html
import shutil
from pathlib import Path

import numpy as np
import pandas as pd

import build_paper1_style_visual_atlas as atlas


HTML_NAME = "kthyro_ctc_emt_submission_war_room.html"
REPORT_NAME = "FINAL_SUBMISSION_WAR_ROOM_KR.md"
SEND_ORDER_TSV = "final_submission_send_order.tsv"
ACTION_TSV = "final_prof_yu_inocras_action_table.tsv"
REHEARSAL_TSV = "final_samsung_reviewer_rehearsal_qa.tsv"
FIG_WAR = "PF21_final_submission_war_room.png"
FIG_7DAY = "PF22_7day_execution_lock.png"
ASSET = f"assets/{atlas.ASSET}"


SEND_ORDER = [
    ("1", "Final Impact Brief", "첫 30초 funding logic", "kthyro_ctc_emt_final_impact_brief.html", "교수님/삼성 reviewer에게 제일 먼저"),
    ("2", "Impact Upgrade", "impact canvas + outcome scenarios", "kthyro_ctc_emt_impact_upgrade.html", "대박 논리와 실패해도 남는 value"),
    ("3", "Paper-Ready Pack", "figure legend/table caption/section map", "kthyro_ctc_emt_paper_ready_figure_table_pack.html", "원고/제안서 조립용"),
    ("4", "Representative Overview", "figure/table/data map", "kthyro_ctc_emt_representative_overview.html", "데이터 어디 썼는지 설명"),
    ("5", "Connected Reader", "쉬운 guided path", "kthyro_ctc_emt_connected_reader.html", "처음 보는 사람용"),
    ("6", "Transfer ZIP", "전체 산출물", f"assets/{atlas.ASSET}/kthyro_ctc_emt_paper_flow_packet_2026_05_09.zip", "파일로 넘길 때"),
]


ACTION_ROWS = [
    {
        "owner": "Prof. Yu / Hospital",
        "ask": "Yu 2024 serial CTC raw phenotype table availability",
        "needed_fields": "patient ID key, timepoint T0/T1/T2, total CTC count, E/EM/M subtype, marker intensity, QC flags",
        "why_needed": "Primary endpoint: postoperative CTC-EMT state transition.",
        "deadline": "Day 1-2",
        "risk_if_missing": "Proposal remains strong but discovery manuscript cannot claim transition.",
    },
    {
        "owner": "Prof. Yu / Hospital",
        "ask": "Matched pathology and follow-up CRF",
        "needed_fields": "LVI, LNM, ETE, tumor size, multifocality, RAI, Tg/anti-Tg, US, persistence/recurrence, follow-up months",
        "why_needed": "Aim 3 secondary association with residual-risk features.",
        "deadline": "Day 2-3",
        "risk_if_missing": "Only biological transition feasibility can be discussed; risk association delayed.",
    },
    {
        "owner": "Inocras",
        "ask": "Tumor-normal NGS workflow and minimum input/QC",
        "needed_fields": "FFPE input, normal source, panel/exome choice, TERT/BRAF/fusion support, VAF threshold, failure rate",
        "why_needed": "Genetic anchor for cfDNA/CTC signal interpretation.",
        "deadline": "Day 3",
        "risk_if_missing": "Genetic map becomes under-specified and reviewer attack strengthens.",
    },
    {
        "owner": "Inocras",
        "ask": "cfDNA and CTC-enriched/low-input feasibility",
        "needed_fields": "blood volume, extraction method, UMI support, limit of detection, duplicate error model, report format",
        "why_needed": "Aim 2 feasibility and stage-gated NGS subset.",
        "deadline": "Day 3-4",
        "risk_if_missing": "Must downgrade CTC/cfDNA NGS to exploratory option.",
    },
    {
        "owner": "Internal",
        "ask": "12-slide deck compression",
        "needed_fields": "Final Impact Brief -> Aim flow -> budget -> claim boundaries",
        "why_needed": "Samsung Science reviewer must see one narrow question.",
        "deadline": "Day 5",
        "risk_if_missing": "Proposal may look like broad platform again.",
    },
    {
        "owner": "Internal",
        "ask": "Reviewer rehearsal",
        "needed_fields": "fatal flaw closure table, expected attack answers, kill criteria",
        "why_needed": "Preempt Science reviewer objections.",
        "deadline": "Day 6-7",
        "risk_if_missing": "Overclaim/feasibility attack remains unhandled.",
    },
]


REHEARSAL_QA = [
    ("Why should Samsung fund this if PTC prognosis is usually good?", "Because recurrence/persistence is a subset problem and current surveillance does not directly read postoperative residual-risk cellular/genetic state. The primary endpoint is biological transition, not immediate clinical deployment.", "PF19, slide01, PF20"),
    ("Is this just CTC counting?", "No. The novelty is serial CTC-EMT state transition after thyroidectomy plus matched tumor-normal/cfDNA genetic anchoring.", "PF17, PF03, F08"),
    ("Do you already have postoperative CTC transition data?", "No, and we should say that clearly. Yu 2024 supports feasibility; the proposal funds the missing serial genetic map.", "PF04, PF12, Final Impact Brief"),
    ("Why not Technology/ICT?", "No proprietary device or AI primitive is the core. The biological perturbation is thyroidectomy; ICT is interpretation infrastructure.", "PF11, slide10"),
    ("What if CTC/cfDNA NGS fails?", "The design is stage-gated. Failure becomes a defined feasibility boundary and does not invalidate the primary CTC-EMT transition endpoint.", "PF14, PF18"),
    ("What is the strongest positive outcome?", "Persistent postoperative EM/M CTC with tissue/cfDNA genetic match. That becomes a residual-risk biology candidate requiring scale-up validation.", "PF18, PF20"),
]


def esc(x) -> str:
    return html.escape("" if x is None else str(x), quote=True)


def draw_war_room(path: Path) -> None:
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

    def text(s, x, y, scale=0.58, color=coal, thick=2):
        cv2.putText(img, s, (x, y), cv2.FONT_HERSHEY_SIMPLEX, scale, color, thick, cv2.LINE_AA)

    def box(x1, y1, x2, y2, color, head, lines):
        cv2.rectangle(img, (x1, y1), (x2, y2), color, -1)
        cv2.rectangle(img, (x1, y1), (x2, y2), coal, 3)
        text(head, x1 + 24, y1 + 50, 0.78, red, 2)
        yy = y1 + 95
        for line in lines:
            text(line, x1 + 24, yy, 0.55, coal, 2)
            yy += 34

    def arrow(p1, p2, color=red):
        cv2.arrowedLine(img, p1, p2, color, 4, tipLength=0.04, line_type=cv2.LINE_AA)

    text("Final Submission War Room", 70, 95, 1.7, coal, 4)
    text("One page to send, ask, rehearse, and lock the Samsung Science proposal", 75, 145, 0.74, red, 2)

    box(90, 250, 680, 570, cream, "Send first", [
        "Final Impact Brief",
        "Impact Upgrade",
        "Paper-Ready Pack",
        "Representative Overview"
    ])
    box(90, 760, 680, 1080, good, "Ask Prof. Yu", [
        "serial CTC phenotype",
        "marker intensity / QC flags",
        "pathology + follow-up CRF",
        "sample/timepoint definitions"
    ])
    box(1010, 250, 1590, 570, warn, "Ask Inocras", [
        "tumor-normal NGS QC",
        "cfDNA limit of detection",
        "CTC-enriched low-input NGS",
        "report fields / turnaround"
    ])
    box(1010, 760, 1590, 1080, good, "Lock proposal", [
        "one primary endpoint",
        "stage-gated NGS",
        "claim boundary table",
        "30억/3년 execution"
    ])
    box(1920, 250, 2500, 570, rose, "Rehearse attacks", [
        "No local CTC data?",
        "CTC phenotype nonspecific?",
        "Why Science not ICT?",
        "What if NGS fails?"
    ])
    box(1920, 760, 2500, 1080, cream, "Final output", [
        "Samsung Science deck",
        "3-page preliminary results",
        "IRB/data checklist",
        "vendor questionnaire"
    ])

    arrow((680, 410), (1010, 410))
    arrow((1590, 410), (1920, 410))
    arrow((680, 920), (1010, 920), green)
    arrow((1590, 920), (1920, 920), green)
    arrow((1295, 570), (1295, 760), green)

    cv2.rectangle(img, (260, 1240), (2340, 1405), warn, -1)
    cv2.rectangle(img, (260, 1240), (2340, 1405), coal, 3)
    text("Final rule", 300, 1300, 0.72, red, 2)
    text("Everything must point back to one question: does thyroidectomy perturb PTC CTC-EMT state, and can persistent signal be genetically anchored?", 300, 1360, 0.58, coal, 2)
    cv2.imwrite(str(path), img)


def draw_7day(path: Path) -> None:
    import cv2

    w, h = 2600, 1500
    img = np.full((h, w, 3), (238, 244, 248), dtype=np.uint8)
    coal = (47, 33, 23)
    red = (37, 45, 143)
    green = (80, 107, 66)
    cream = (241, 250, 255)
    good = (232, 243, 233)
    warn = (212, 241, 255)

    def text(s, x, y, scale=0.55, color=coal, thick=2):
        cv2.putText(img, s, (x, y), cv2.FONT_HERSHEY_SIMPLEX, scale, color, thick, cv2.LINE_AA)

    text("7-Day Execution Lock", 70, 95, 1.65, coal, 4)
    text("Convert impact package into submission-ready Samsung Science proposal", 75, 145, 0.74, red, 2)

    days = [
        ("D1", "Prof. Yu", "serial CTC raw table"),
        ("D2", "CRF", "pathology/follow-up lock"),
        ("D3", "Inocras", "NGS/cfDNA feasibility"),
        ("D4", "Assay", "CTC marker panel lock"),
        ("D5", "Deck", "12-slide final flow"),
        ("D6", "Budget/IRB", "30억 + checklist"),
        ("D7", "Mock review", "attack answers locked"),
    ]
    x = 90
    y = 330
    for i, (d, head, body) in enumerate(days):
        color = good if i in [0, 1, 2] else cream if i in [3, 4] else warn
        cv2.rectangle(img, (x, y), (x + 330, y + 580), color, -1)
        cv2.rectangle(img, (x, y), (x + 330, y + 580), coal, 3)
        text(d, x + 28, y + 70, 1.1, red, 3)
        text(head, x + 28, y + 135, 0.70, coal, 2)
        yy = y + 215
        for line in body.split("/"):
            text(line, x + 28, yy, 0.52, coal, 2)
            yy += 38
        if i < len(days) - 1:
            cv2.arrowedLine(img, (x + 330, y + 290), (x + 390, y + 290), red, 4, tipLength=0.12, line_type=cv2.LINE_AA)
        x += 360

    cv2.rectangle(img, (220, 1080), (2380, 1280), (255, 255, 255), -1)
    cv2.rectangle(img, (220, 1080), (2380, 1280), coal, 3)
    text("Go/no-go by Day 7", 260, 1140, 0.78, red, 2)
    text("If Yu serial phenotype and Inocras NGS feasibility are aligned, submit as narrow Samsung Science perturbation proposal.", 260, 1200, 0.62, coal, 2)
    text("If not aligned, keep proposal but downgrade genetic-map claims to stage-gated feasibility.", 260, 1245, 0.58, coal, 2)
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


def html_page(send_df: pd.DataFrame, action_df: pd.DataFrame, qa_df: pd.DataFrame) -> str:
    return f"""<!DOCTYPE html>
<html lang="ko">
<head>
<meta charset="UTF-8" />
<meta name="viewport" content="width=device-width, initial-scale=1.0" />
<title>K-Thyro CTC-EMT Submission War Room</title>
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
  .layout{{display:grid;grid-template-columns:250px 1fr;gap:26px;margin-top:30px}} .toc{{position:sticky;top:18px;align-self:start;background:rgba(255,250,241,.92);border:1px solid var(--line);border-radius:12px;padding:16px}} .toc h3{{font-family:"Cormorant Garamond",serif;color:var(--red);font-size:24px;margin:0 0 10px}} .toc a{{display:block;color:var(--ink);border-bottom:1px dotted var(--line);padding:8px 0;font-size:14px}}
  section{{background:rgba(255,250,241,.9);border:1px solid var(--line);border-radius:14px;padding:28px;margin-bottom:24px}} h2{{font-family:"Cormorant Garamond",serif;font-size:42px;line-height:1.05;margin:0 0 12px;color:var(--coal)}} .section-note{{color:var(--muted);font-style:italic;margin:0 0 18px}}
  .big-fig{{display:block;background:#fff;border:1px solid var(--line);border-radius:14px;padding:14px;margin:14px 0 20px}} .big-fig img{{display:block;width:100%;height:auto;border-radius:8px;background:#f4eadb}}
  .table-wrap{{overflow-x:auto}} table{{width:100%;border-collapse:collapse;background:#fff;border:1px solid var(--line);font-size:14px}} th{{background:var(--coal);color:#fff;text-align:left;font-family:"JetBrains Mono",monospace;font-size:11px;letter-spacing:.06em;text-transform:uppercase;padding:10px}} td{{border-top:1px solid var(--line2);padding:10px;vertical-align:top}}
  .link-grid{{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));gap:12px}} .link-card{{display:block;background:#fff;border:1px solid var(--line2);border-radius:12px;padding:16px}} .link-card b{{display:block;font-family:"Cormorant Garamond",serif;font-size:27px;color:var(--red)}} .link-card span{{display:block;color:var(--muted)}}
  .footer{{font-family:"JetBrains Mono",monospace;font-size:12px;color:var(--muted);text-align:center;margin-top:34px;padding-top:18px;border-top:1px solid var(--line)}}
  @media(max-width:1000px){{.layout{{grid-template-columns:1fr}}.toc{{position:relative;top:auto}}.link-grid{{grid-template-columns:1fr}}h1{{font-size:42px}}.hero{{padding:30px 22px}}}}
</style>
</head>
<body>
<div class="wrap">
  <nav class="top">
    <a href="kthyro_ctc_emt_all_pages_launchpad.html">All Pages</a>
    <a href="kthyro_ctc_emt_final_impact_brief.html">Final Impact</a>
    <a href="kthyro_ctc_emt_impact_upgrade.html">Impact Upgrade</a>
    <a href="kthyro_ctc_emt_paper_ready_figure_table_pack.html">Paper-Ready Pack</a>
  </nav>
  <header class="hero">
    <div class="kicker">K-Thyro CTC-EMT · final submission war room</div>
    <h1>Submission War Room</h1>
    <p class="subtitle">이제 볼 자료가 아니라 보낼 자료다. 교수님께 보낼 순서, 병원/Yu/Inocras 요청, 삼성 reviewer 리허설 Q&A, 7일 실행 lock을 한 페이지로 묶었다.</p>
  </header>
  <div class="layout">
    <aside class="toc">
      <h3>Read Order</h3>
      <a href="#war">War Room Map</a>
      <a href="#send">Send Order</a>
      <a href="#actions">Action Table</a>
      <a href="#qa">Reviewer Q&A</a>
      <a href="#seven">7-Day Lock</a>
      <a href="#downloads">Downloads</a>
    </aside>
    <main>
      <section id="war">
        <h2>War Room Map</h2>
        <a class="big-fig" href="{ASSET}/{FIG_WAR}"><img src="{ASSET}/{FIG_WAR}" alt="Final submission war room" /></a>
      </section>
      <section id="send">
        <h2>Send Order</h2>
        <p class="section-note">교수님께는 이 순서로 링크를 보내면 된다.</p>
        {table_html(send_df)}
      </section>
      <section id="actions">
        <h2>Prof. Yu / Inocras Action Table</h2>
        <p class="section-note">무슨 데이터를 누구에게 왜 요청하는지.</p>
        {table_html(action_df)}
      </section>
      <section id="qa">
        <h2>Samsung Reviewer Rehearsal Q&A</h2>
        <p class="section-note">리뷰어 공격에 대한 짧고 방어 가능한 답변.</p>
        {table_html(qa_df)}
      </section>
      <section id="seven">
        <h2>7-Day Execution Lock</h2>
        <a class="big-fig" href="{ASSET}/{FIG_7DAY}"><img src="{ASSET}/{FIG_7DAY}" alt="7-day execution lock" /></a>
      </section>
      <section id="downloads">
        <h2>Downloads</h2>
        <div class="link-grid">
          <a class="link-card" href="{SEND_ORDER_TSV}"><b>Send Order TSV</b><span>링크 전송 순서</span></a>
          <a class="link-card" href="{ACTION_TSV}"><b>Action TSV</b><span>Yu/Inocras 요청표</span></a>
          <a class="link-card" href="{REHEARSAL_TSV}"><b>Q&A TSV</b><span>reviewer rehearsal</span></a>
          <a class="link-card" href="{REPORT_NAME}"><b>Guide MD</b><span>final handoff memo</span></a>
        </div>
      </section>
    </main>
  </div>
  <div class="footer">Submission War Room: from impressive package to executable submission plan.</div>
</div>
</body>
</html>"""


def report_md(send_df: pd.DataFrame, action_df: pd.DataFrame, qa_df: pd.DataFrame) -> str:
    lines = [
        "# Final Submission War Room",
        "",
        "Purpose: 교수님께 보낼 순서, 병원/Yu/Inocras 요청, 삼성 reviewer Q&A, 7일 실행 lock을 한 문서에 묶는다.",
        "",
        "Core rule: everything must point back to one question: does thyroidectomy perturb PTC CTC-EMT state, and can persistent signal be genetically anchored?",
        "",
        "## Send order",
    ]
    for _, r in send_df.iterrows():
        lines.append(f"- {r['order']}. {r['page']}: {r['url']} — {r['purpose']}")
    lines += ["", "## Actions"]
    for _, r in action_df.iterrows():
        lines.append(f"- {r['owner']}: {r['ask']} / Needed: {r['needed_fields']} / Deadline: {r['deadline']}")
    lines += ["", "## Reviewer Q&A"]
    for _, r in qa_df.iterrows():
        lines.append(f"- Q: {r['attack']} A: {r['answer']} Support: {r['support']}")
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
            needle = '    <a href="kthyro_ctc_emt_final_impact_brief.html"'
            insert = f'    <a href="{HTML_NAME}" style="background:linear-gradient(135deg,#426b50 0%,#fffaf1 100%);border:3px solid #102033;color:#102033;font-weight:1000;font-size:18px;padding:16px 20px">★★ CTC-EMT SUBMISSION WAR ROOM · send ask rehearse execute</a>\n'
            pos = text.find(needle)
            text = text[:pos] + insert + text[pos:] if pos >= 0 else text + "\n" + insert
            idx.write_text(text, encoding="utf-8")
        shutil.copy2(idx, atlas.LIVE / "index.html")

    for path in [atlas.REP / "VISUAL_WEB_HANDOFF_KR.md", atlas.HUB / "VISUAL_WEB_HANDOFF_KR.md"]:
        if not path.exists():
            continue
        text = path.read_text(encoding="utf-8")
        line = "- Submission War Room: http://40.82.129.113:8012/papers_hub_2026_05_04/kthyro_ctc_emt_submission_war_room.html\n"
        if HTML_NAME not in text:
            text = text.replace("## 최우선 링크\n\n", "## 최우선 링크\n\n" + line)
            path.write_text(text, encoding="utf-8")
    handoff = atlas.HUB / "VISUAL_WEB_HANDOFF_KR.md"
    if handoff.exists():
        shutil.copy2(handoff, atlas.LIVE / "VISUAL_WEB_HANDOFF_KR.md")

    launch = atlas.HUB / "kthyro_ctc_emt_all_pages_launchpad.html"
    if launch.exists() and HTML_NAME not in launch.read_text(encoding="utf-8"):
        text = launch.read_text(encoding="utf-8")
        insert = f'<a href="{HTML_NAME}"><b>제출 War Room</b><span>Send/Ask/Rehearse</span></a>'
        text = text.replace('<div class="route-grid">', '<div class="route-grid">' + insert, 1)
        launch.write_text(text, encoding="utf-8")
        shutil.copy2(launch, atlas.LIVE / "kthyro_ctc_emt_all_pages_launchpad.html")


def main() -> None:
    fig_dir = atlas.OUT / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    asset_dir = atlas.HUB / "assets" / atlas.ASSET
    live_asset_dir = atlas.LIVE / "assets" / atlas.ASSET

    war = fig_dir / FIG_WAR
    seven = fig_dir / FIG_7DAY
    draw_war_room(war)
    draw_7day(seven)
    for p in [war, seven]:
        shutil.copy2(p, asset_dir / p.name)
        shutil.copy2(p, live_asset_dir / p.name)

    send_df = pd.DataFrame(SEND_ORDER, columns=["order", "page", "purpose", "url", "when_to_send"])
    action_df = pd.DataFrame(ACTION_ROWS)
    qa_df = pd.DataFrame(REHEARSAL_QA, columns=["attack", "answer", "support"])
    copy_table(send_df, SEND_ORDER_TSV)
    copy_table(action_df, ACTION_TSV)
    copy_table(qa_df, REHEARSAL_TSV)

    html_text = html_page(send_df, action_df, qa_df)
    src = atlas.HUB / HTML_NAME
    live = atlas.LIVE / HTML_NAME
    local = atlas.REP / HTML_NAME
    src.write_text(html_text, encoding="utf-8")
    local.write_text(html_text, encoding="utf-8")
    shutil.copy2(src, live)

    md = report_md(send_df, action_df, qa_df)
    report = atlas.REP / REPORT_NAME
    hub_report = atlas.HUB / REPORT_NAME
    report.write_text(md, encoding="utf-8")
    hub_report.write_text(md, encoding="utf-8")
    shutil.copy2(hub_report, atlas.LIVE / REPORT_NAME)

    update_index_handoff_launchpad()
    print(f"HTML={src}")
    print(f"LIVE={live}")
    print(f"FIGURES={war}, {seven}")
    print(f"REPORT={report}")
    print(f"TABLES={atlas.TAB / SEND_ORDER_TSV}, {atlas.TAB / ACTION_TSV}, {atlas.TAB / REHEARSAL_TSV}")


if __name__ == "__main__":
    main()
