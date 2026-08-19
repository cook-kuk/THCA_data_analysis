#!/usr/bin/env python3
from __future__ import annotations

import shutil

import build_paper1_style_visual_atlas as atlas


HTML_NAME = "kthyro_ctc_emt_all_pages_launchpad.html"
REPORT_NAME = "ALL_PAGES_LAUNCHPAD_GUIDE_KR.md"


GROUPS = [
    (
        "Start Here",
        "처음 보는 사람은 이 줄만 보면 된다.",
        [
            ("Submission War Room", "교수님 전송 순서, Yu/Inocras 요청, 삼성 reviewer Q&A, 7일 실행표.", "kthyro_ctc_emt_submission_war_room.html", "best"),
            ("Final Impact Brief", "첫 30초용 funding logic pyramid와 reviewer scorecard.", "kthyro_ctc_emt_final_impact_brief.html", "best"),
            ("Impact Upgrade Pack", "삼성 reviewer용 impact canvas와 세 가지 outcome scenario.", "kthyro_ctc_emt_impact_upgrade.html", "best"),
            ("Paper-Ready Figure/Table Pack", "원고/제안서에 바로 꽂는 figure legend, table caption, section map.", "kthyro_ctc_emt_paper_ready_figure_table_pack.html", "best"),
            ("Representative Overview", "대표 figure/table/data map. 어떤 데이터가 어디에 쓰였는지 한방 정리.", "kthyro_ctc_emt_representative_overview.html", "best"),
            ("Breakthrough Case", "가능성에서 fundable Science case로 압축한 대박 논리.", "kthyro_ctc_emt_breakthrough_case.html", "best"),
            ("Connected Reader", "가장 쉬운 guided path. 5분/15분/리뷰어 모드.", "kthyro_ctc_emt_connected_reader.html", "best"),
            ("Manuscript Storyboard", "Figure 1-6 → Results → boundary 논문 흐름.", "kthyro_ctc_emt_manuscript_storyboard_paper1_style.html", "best"),
            ("Visual Web Handoff", "교수님께 보낼 링크와 보는 순서.", "VISUAL_WEB_HANDOFF_KR.md", "best"),
        ],
    ),
    (
        "Figure Pages",
        "그림을 크게 나눠 보는 페이지들.",
        [
            ("Split Visual Atlas Hub", "분리형 atlas의 입구.", "kthyro_ctc_emt_visual_atlas_split_hub.html", ""),
            ("Evidence Figures", "public-data in silico 결과만.", "kthyro_ctc_emt_visual_atlas_evidence.html", ""),
            ("Reviewer Defense Figures", "예상 공격, claim boundary, dataset gap.", "kthyro_ctc_emt_visual_atlas_defense.html", ""),
            ("Manuscript Boards", "blueprint, scorecard, execution board.", "kthyro_ctc_emt_visual_atlas_manuscript.html", ""),
            ("Samsung Deck Atlas", "PPT 12장 그림 크게 보기.", "kthyro_ctc_emt_visual_atlas_deck.html", ""),
            ("Control Tables", "CRF, SAP, kill criteria, safe sentence.", "kthyro_ctc_emt_visual_atlas_tables.html", ""),
            ("One-page Visual Atlas", "모든 그림을 한 페이지에 길게.", "kthyro_ctc_emt_visual_atlas_paper1_style.html", ""),
            ("Old CV2 Contact Sheet", "예전 한 장짜리 atlas.", f"assets/{atlas.ASSET}/PF08_cv2_visual_atlas.png", "old"),
        ],
    ),
    (
        "Master / Proposal",
        "전체 통제판과 삼성 proposal 산출물.",
        [
            ("Master Dossier", "표, claim, 그림, 논문 흐름 전체.", "kthyro_ctc_emt_paper_flow_master_dossier.html", ""),
            ("Samsung Science Dossier", "Science-track proposal dossier.", "kthyro_ctc_emt_samsung_science_dossier.html", ""),
            ("No-Yu In Silico Dossier", "유형원 데이터 없이 가능한 분석 범위.", "kthyro_no_yu_in_silico_ctc_emt_dossier.html", ""),
            ("Samsung PPTX", "제안서 발표용 PPTX.", f"assets/{atlas.ASSET}/kthyro_ctc_emt_samsung_science_pitch_v2_reference_style.pptx", "file"),
            ("Transfer ZIP", "전체 figures/tables/reports packet.", f"assets/{atlas.ASSET}/kthyro_ctc_emt_paper_flow_packet_2026_05_09.zip", "file"),
        ],
    ),
    (
        "Reports / Tables",
        "글쓰기와 데이터 정의에 바로 쓰는 파일.",
        [
            ("Connected Reader Guide", "연결 읽기 경로 markdown.", "CONNECTED_READER_GUIDE_KR.md", ""),
            ("Storyboard TSV", "Figure-to-results machine table.", "manuscript_storyboard_web.tsv", ""),
            ("Storyboard MD", "Figure-to-results writing scaffold.", "MANUSCRIPT_STORYBOARD_WEB_KR.md", ""),
            ("All Pages Guide", "이 launchpad 설명.", REPORT_NAME, ""),
            ("Data Definition Master", "쓴 데이터와 claim boundary.", f"assets/{atlas.ASSET}/DATA_DEFINITION_MASTER_KR.md", ""),
            ("Caption/Results Scaffold", "각 그림 caption/results scaffold.", f"assets/{atlas.ASSET}/FIGURE_CAPTION_AND_RESULTS_SCAFFOLD_KR.md", ""),
            ("Manuscript Blueprint", "논문 identity와 Results blocks.", f"assets/{atlas.ASSET}/MANUSCRIPT_BLUEPRINT_KR.md", ""),
            ("Execution Package", "next action과 reviewer risk.", f"assets/{atlas.ASSET}/MANUSCRIPT_EXECUTION_PACKAGE_KR.md", ""),
        ],
    ),
]


def card(title: str, desc: str, href: str, tag: str) -> str:
    label = {
        "best": "START",
        "old": "OLD",
        "file": "FILE",
    }.get(tag, "OPEN")
    cls = f" card {tag}" if tag else "card"
    return f"""
    <a class="{cls}" href="{href}">
      <span class="tag">{label}</span>
      <h3>{title}</h3>
      <p>{desc}</p>
      <code>{href}</code>
    </a>
    """


def group_html(title: str, note: str, items: list[tuple[str, str, str, str]]) -> str:
    return f"""
    <section id="{title.lower().replace(' ', '-').replace('/', '-')}">
      <div class="section-head">
        <h2>{title}</h2>
        <p>{note}</p>
      </div>
      <div class="card-grid">{''.join(card(*x) for x in items)}</div>
    </section>
    """


def html_page() -> str:
    total = sum(len(g[2]) for g in GROUPS)
    toc = "".join(
        f"<a href='#{title.lower().replace(' ', '-').replace('/', '-')}'>{title}</a>" for title, _, _ in GROUPS
    )
    quick = """
    <div class="route-grid">
      <a href="kthyro_ctc_emt_submission_war_room.html"><b>제출 War Room</b><span>Send/Ask/Rehearse</span></a>
      <a href="kthyro_ctc_emt_final_impact_brief.html"><b>최종 임팩트</b><span>Final 30-sec Brief</span></a>
      <a href="kthyro_ctc_emt_impact_upgrade.html"><b>임팩트 업</b><span>Reviewer Impact Canvas</span></a>
      <a href="kthyro_ctc_emt_paper_ready_figure_table_pack.html"><b>원고 조립</b><span>Figure/Table Pack</span></a>
      <a href="kthyro_ctc_emt_representative_overview.html"><b>대표 Overview</b><span>Figure/Table/Data Map</span></a>
      <a href="kthyro_ctc_emt_breakthrough_case.html"><b>대박 논리</b><span>Breakthrough Case</span></a>
      <a href="kthyro_ctc_emt_connected_reader.html"><b>처음 설명</b><span>Connected Reader</span></a>
      <a href="kthyro_ctc_emt_manuscript_storyboard_paper1_style.html"><b>논문 흐름</b><span>Storyboard</span></a>
      <a href="kthyro_ctc_emt_visual_atlas_evidence.html"><b>그림 근거</b><span>Evidence</span></a>
      <a href="kthyro_ctc_emt_visual_atlas_defense.html"><b>심사 방어</b><span>Defense</span></a>
      <a href="kthyro_ctc_emt_visual_atlas_tables.html"><b>표/CRF</b><span>Tables</span></a>
      <a href="assets/kthyro_ctc_emt_paper_flow/kthyro_ctc_emt_paper_flow_packet_2026_05_09.zip"><b>전체 파일</b><span>ZIP</span></a>
    </div>
    """
    return f"""<!DOCTYPE html>
<html lang="ko">
<head>
<meta charset="UTF-8" />
<meta name="viewport" content="width=device-width, initial-scale=1.0" />
<title>K-Thyro CTC-EMT All Pages Launchpad</title>
<link rel="preconnect" href="https://fonts.googleapis.com" />
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin />
<link href="https://fonts.googleapis.com/css2?family=Cormorant+Garamond:wght@500;600;700&family=Newsreader:ital,wght@0,400;0,500;0,600;1,400&family=JetBrains+Mono:wght@400;500;700&family=Noto+Sans+KR:wght@400;500;700;800&display=swap" rel="stylesheet" />
<style>
  :root{{--ink:#102033;--muted:#52606f;--paper:#f6f0e6;--card:#fffaf1;--line:#d4c6b3;--line2:#eadfcc;--red:#8f2d25;--blue:#244e73;--green:#426b50;--gold:#b58534;--coal:#17212f;--good:#e9f3e8;--bad:#f7dedb}}
  *{{box-sizing:border-box}} html{{scroll-behavior:smooth}}
  body{{margin:0;background:linear-gradient(135deg,#f8f2e8 0%,#efe3d1 58%,#f7efe0 100%);color:var(--ink);font-family:"Newsreader","Noto Sans KR",serif;line-height:1.6;font-size:17px}}
  a{{color:inherit;text-decoration:none}} .wrap{{max-width:1240px;margin:0 auto;padding:32px 26px 80px}}
  .top{{display:flex;gap:10px;align-items:center;flex-wrap:wrap;margin-bottom:22px;font-family:"JetBrains Mono",monospace;font-size:12px}}
  .top a{{border:1px solid var(--line);background:rgba(255,255,255,.68);border-radius:999px;padding:8px 11px}}
  .hero{{border:2px solid var(--coal);background:linear-gradient(145deg,#fffaf1 0%,#f0dfc4 100%);padding:42px 46px;box-shadow:10px 10px 0 rgba(23,33,47,.16)}}
  .kicker{{font-family:"JetBrains Mono",monospace;letter-spacing:.17em;text-transform:uppercase;font-size:12px;color:var(--red);font-weight:800;margin-bottom:14px}}
  h1{{font-family:"Cormorant Garamond",serif;font-size:58px;line-height:1.02;margin:0 0 12px;color:var(--coal)}}
  .subtitle{{font-size:22px;color:var(--muted);max-width:980px;font-style:italic}}
  .stats{{display:flex;flex-wrap:wrap;gap:10px;margin-top:20px}} .stat{{font-family:"JetBrains Mono",monospace;font-size:11px;border:1px solid var(--line);background:#fff;border-radius:6px;padding:7px 10px}} .stat.hot{{background:var(--bad);color:var(--red);font-weight:800}} .stat.good{{background:var(--good);color:var(--green);font-weight:800}}
  .layout{{display:grid;grid-template-columns:250px 1fr;gap:26px;margin-top:30px}} .toc{{position:sticky;top:18px;align-self:start;background:rgba(255,250,241,.92);border:1px solid var(--line);border-radius:12px;padding:16px}} .toc h3{{font-family:"Cormorant Garamond",serif;font-size:24px;color:var(--red);margin:0 0 10px}} .toc a{{display:block;border-bottom:1px dotted var(--line);padding:8px 0;font-size:14px}}
  section{{background:rgba(255,250,241,.9);border:1px solid var(--line);border-radius:14px;padding:28px;margin-bottom:24px}} .section-head{{display:flex;justify-content:space-between;gap:18px;align-items:flex-end;border-bottom:1px solid var(--line2);margin-bottom:16px;padding-bottom:10px}} .section-head h2{{font-family:"Cormorant Garamond",serif;font-size:39px;line-height:1;margin:0;color:var(--coal)}} .section-head p{{margin:0;color:var(--muted);font-style:italic;max-width:560px}}
  .route-grid{{display:grid;grid-template-columns:repeat(6,minmax(0,1fr));gap:10px;margin-top:22px}} .route-grid a{{background:#fff;border:1px solid var(--line2);border-radius:12px;padding:15px;min-height:96px}} .route-grid b{{display:block;color:var(--red);font-family:"JetBrains Mono",monospace;font-size:12px;text-transform:uppercase;letter-spacing:.06em}} .route-grid span{{display:block;font-family:"Cormorant Garamond",serif;font-size:25px;margin-top:8px}}
  .card-grid{{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:14px}} .card{{display:block;background:#fff;border:1px solid var(--line2);border-radius:12px;padding:16px;min-height:184px;position:relative;box-shadow:0 5px 16px rgba(23,33,47,.05)}} .card.best{{border-top:5px solid var(--red)}} .card.file{{border-top:5px solid var(--green)}} .card.old{{opacity:.78}} .tag{{font-family:"JetBrains Mono",monospace;font-size:10px;background:var(--coal);color:#fff;border-radius:999px;padding:5px 8px}} .card.best .tag{{background:var(--red)}} .card.file .tag{{background:var(--green)}} .card h3{{font-family:"Cormorant Garamond",serif;font-size:29px;line-height:1.05;margin:14px 0 8px;color:var(--coal)}} .card p{{margin:0;color:var(--muted)}} .card code{{display:block;margin-top:12px;font-family:"JetBrains Mono",monospace;font-size:10px;color:var(--blue);overflow-wrap:anywhere;background:#f8f1e5;border:1px solid var(--line2);border-radius:6px;padding:7px}}
  .callout{{border-left:5px solid var(--red);background:#fff;padding:18px 20px;margin-top:18px;border-radius:0 10px 10px 0}} .callout b{{color:var(--red)}} .footer{{font-family:"JetBrains Mono",monospace;font-size:12px;color:var(--muted);text-align:center;margin-top:34px;padding-top:18px;border-top:1px solid var(--line)}}
  @media(max-width:1000px){{.layout{{grid-template-columns:1fr}}.toc{{position:relative;top:auto}}.route-grid,.card-grid{{grid-template-columns:repeat(2,minmax(0,1fr))}}h1{{font-size:42px}}.hero{{padding:30px 22px}}.section-head{{display:block}}}}
  @media(max-width:620px){{.route-grid,.card-grid{{grid-template-columns:1fr}}.wrap{{padding:18px 14px}}}}
</style>
</head>
<body>
  <div class="wrap">
    <nav class="top">
      <a href="index.html">8-Papers Hub</a>
      <a href="kthyro_ctc_emt_connected_reader.html">Connected Reader</a>
      <a href="kthyro_ctc_emt_manuscript_storyboard_paper1_style.html">Storyboard</a>
      <a href="kthyro_ctc_emt_visual_atlas_split_hub.html">Split Atlas</a>
    </nav>
    <header class="hero">
      <div class="kicker">K-Thyro CTC-EMT · one-click all-pages launchpad</div>
      <h1>All Pages Launchpad</h1>
      <p class="subtitle">모든 웹페이지, 보고서, 표, PPTX, ZIP을 한 화면에서 바로 들어가게 만든 최상단 포털. 처음 보는 사람은 위에서 아래로 누르면 된다.</p>
      <div class="stats">
        <span class="stat hot">Main entry page</span>
        <span class="stat good">{total} direct links</span>
        <span class="stat">Science-track boundary preserved</span>
        <span class="stat">Public-data prior, not CTC discovery</span>
      </div>
      {quick}
    </header>
    <div class="layout">
      <aside class="toc"><h3>All Links</h3>{toc}</aside>
      <main>
        <section id="rule">
          <div class="section-head"><h2>Use Rule</h2><p>한 문장 boundary를 먼저 보고 이동한다.</p></div>
          <div class="callout"><b>Core boundary.</b> 현재 산출물은 public-data tissue-state / marker / NGS-anchor prior framework다. Yu/hospital serial CTC data 없이 postoperative CTC-EMT transition discovery나 recurrence prediction으로 쓰지 않는다.</div>
        </section>
        {''.join(group_html(*g) for g in GROUPS)}
      </main>
    </div>
    <div class="footer">Generated launchpad for K-Thyro CTC-EMT Samsung Science proposal package.</div>
  </div>
</body>
</html>"""


def guide_md() -> str:
    lines = [
        "# All Pages Launchpad Guide",
        "",
        "Purpose: 모든 K-Thyro CTC-EMT 웹페이지와 파일을 한 화면에서 들어가게 만든 최상단 포털.",
        "",
        "Main URL: http://40.82.129.113:8012/papers_hub_2026_05_04/kthyro_ctc_emt_all_pages_launchpad.html",
        "",
        "Recommended first path: Connected Reader → Manuscript Storyboard → Evidence Figures → Reviewer Defense Figures → Control Tables.",
        "",
        "Core boundary: public data는 prior/framework이고, postoperative CTC-EMT transition discovery는 Yu/hospital serial CTC data가 필요하다.",
        "",
    ]
    for title, note, items in GROUPS:
        lines += [f"## {title}", note, ""]
        for name, desc, href, _ in items:
            lines.append(f"- {name}: {href} — {desc}")
        lines.append("")
    return "\n".join(lines)


def update_index_and_handoff() -> None:
    idx = atlas.HUB / "index.html"
    if idx.exists():
        text = idx.read_text(encoding="utf-8")
        if HTML_NAME not in text:
            needle = '    <a href="kthyro_ctc_emt_connected_reader.html"'
            insert = f'    <a href="{HTML_NAME}" style="background:linear-gradient(135deg,#fffaf1 0%,#f0dfc4 100%);border:3px solid #102033;color:#102033;font-weight:1000;font-size:18px;padding:16px 20px">★★ CTC-EMT ALL PAGES LAUNCHPAD · one-click portal</a>\n'
            pos = text.find(needle)
            text = text[:pos] + insert + text[pos:] if pos >= 0 else text + "\n" + insert
            idx.write_text(text, encoding="utf-8")
        shutil.copy2(idx, atlas.LIVE / "index.html")

    for path in [atlas.REP / "VISUAL_WEB_HANDOFF_KR.md", atlas.HUB / "VISUAL_WEB_HANDOFF_KR.md"]:
        if not path.exists():
            continue
        text = path.read_text(encoding="utf-8")
        line = "- All Pages Launchpad: http://40.82.129.113:8012/papers_hub_2026_05_04/kthyro_ctc_emt_all_pages_launchpad.html\n"
        if "kthyro_ctc_emt_all_pages_launchpad.html" not in text:
            text = text.replace("## 최우선 링크\n\n", "## 최우선 링크\n\n" + line)
            path.write_text(text, encoding="utf-8")
    handoff = atlas.HUB / "VISUAL_WEB_HANDOFF_KR.md"
    if handoff.exists():
        shutil.copy2(handoff, atlas.LIVE / "VISUAL_WEB_HANDOFF_KR.md")


def main() -> None:
    src = atlas.HUB / HTML_NAME
    live = atlas.LIVE / HTML_NAME
    local = atlas.REP / HTML_NAME
    guide = atlas.REP / REPORT_NAME
    hub_guide = atlas.HUB / REPORT_NAME

    html = html_page()
    src.write_text(html, encoding="utf-8")
    local.write_text(html, encoding="utf-8")
    shutil.copy2(src, live)

    guide_text = guide_md()
    guide.write_text(guide_text, encoding="utf-8")
    hub_guide.write_text(guide_text, encoding="utf-8")
    shutil.copy2(hub_guide, atlas.LIVE / REPORT_NAME)
    update_index_and_handoff()

    print(f"HTML={src}")
    print(f"LIVE={live}")
    print(f"GUIDE={guide}")


if __name__ == "__main__":
    main()
