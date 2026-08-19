#!/usr/bin/env python3
from __future__ import annotations

import shutil

import build_paper1_style_visual_atlas as atlas


HTML_NAME = "kthyro_ctc_emt_breakthrough_case.html"
REPORT_NAME = "BREAKTHROUGH_CASE_KR.md"


ASSET = f"assets/{atlas.ASSET}"


BREAKTHROUGH = [
    {
        "label": "Big Idea",
        "title": "수술을 치료가 아니라 perturbation으로 읽는다",
        "body": "갑상선절제술은 종양부하를 갑자기 줄이는 인간 내 controlled perturbation이다. 이 전후로 CTC-EMT 상태가 collapse, persist, shift 중 무엇을 보이는지 읽으면 PTC residual-risk biology를 시간축에서 볼 수 있다.",
        "boundary": "현재 public data로 이 transition을 발견했다고 쓰지 않는다. Samsung/Yu cohort에서 검정할 핵심 biological question이다.",
    },
    {
        "label": "Why Now",
        "title": "Yu 2024가 feasibility를 열었고, genetic map이 비어 있다",
        "body": "Yu 2024는 PTC에서 수술 전, 2주, 3개월 serial CTC-EMT 측정이 가능하고 수술 후 CTC count가 감소한다는 anchor를 제공한다. 그러나 matched tissue-normal NGS, cfDNA, CTC-enriched NGS와 연결된 clone-state map은 없다.",
        "boundary": "Yu 2024를 우리 raw data처럼 쓰지 않는다. published feasibility anchor로만 둔다.",
    },
    {
        "label": "Public Pilot Value",
        "title": "mutation alone is not enough를 public data로 먼저 잠근다",
        "body": "TCGA public pilot은 BRAF-mutant PTC도 여러 tissue-state label로 갈라지고, largest BRAF label이 33.2%에 그친다는 점을 보여준다. 따라서 postoperative blood signal은 driver mutation만이 아니라 phenotype-state와 함께 해석해야 한다.",
        "boundary": "이 결과는 CTC shedding, recurrence, viable metastasis를 증명하지 않는다.",
    },
    {
        "label": "Killer Experiment",
        "title": "serial CTC-EMT phenotype에 matched genetic anchor를 붙인다",
        "body": "T0 pre-op, T1 2 weeks, T2 3 months, T3 6-12 months blood를 CTC phenotype, cfDNA, matched tumor-normal NGS, FFPE/mIHC와 연결한다. 지속되는 EM/M CTC와 tissue-matched genetic signal이 residual-risk biology 후보가 된다.",
        "boundary": "clinical diagnostic이나 recurrence predictor로 바로 주장하지 않는다. primary endpoint는 biological state transition이다.",
    },
]


KILLER_FIGS = [
    ("PF11_manuscript_blueprint_board.png", "Fig. A", "Project identity lock", "Science-track biological perturbation으로 고정한다."),
    ("slide02_hypothesis.png", "Fig. B", "Thyroidectomy perturbation", "수술 후 CTC state collapse/persistence/shift를 핵심 readout으로 만든다."),
    ("F03_tcga_braf_state_split.png", "Fig. C", "Mutation is not enough", "BRAF subset도 state-heterogeneous하므로 phenotype-genotype integration이 필요하다."),
    ("PF03_marker_to_ngs_map.png", "Fig. D", "Phenotype-to-genotype bridge", "CTC-EMT marker module과 matched NGS를 한 지도에 묶는다."),
    ("PF04_dataset_gap_bridge.png", "Fig. E", "Missing dataset", "공개 데이터에 없는 joint serial map이 proposal 필요성을 만든다."),
    ("PF12_future_data_unlock_map.png", "Fig. F", "Data unlock plan", "병원/Yu 데이터가 들어오면 어떤 figure와 claim이 unlock되는지 정한다."),
]


REVIEWER_ATTACKS = [
    ("이건 CTC 결과가 아니라 public tissue omics 아닌가?", "맞다. 그래서 current claim은 prior/framework다. CTC transition은 prospective endpoint로 둔다.", "PF02, PF04, PF05"),
    ("AI/ICT 아닌가?", "아니다. AI는 clone-state/longitudinal modeling 인프라일 뿐이고, category identity는 thyroidectomy perturbation biology다.", "PF11, slide10"),
    ("PTC는 예후 좋은데 왜 필요한가?", "대부분 curable이지만 residual/persistent risk subset에서 postoperative biology를 직접 읽는 도구가 부족하다.", "slide01, slide08"),
    ("CTC phenotype은 nonspecific할 수 있다", "그래서 matched tumor-normal NGS, cfDNA, and CTC-enriched NGS subset anchoring을 핵심으로 둔다.", "PF03, F08"),
    ("early PTC에서 cfDNA/CTC NGS 안 될 수 있다", "stage-gated feasibility와 kill criteria를 미리 둔다. 실패하면 clinical claim이 아니라 assay feasibility result로 downgrade한다.", "PF14, tables"),
]


SEVEN_DAY = [
    ("Day 1", "Prof. Yu meeting", "serial CTC raw count/state table, marker intensity, timepoint definitions 요청"),
    ("Day 2", "Hospital CRF lock", "LVI/LNM/ETE, Tg/US, recurrence/persistence, RAI, pathology fields 고정"),
    ("Day 3", "Inocras check", "tumor-normal NGS, cfDNA, low-input CTC-enriched NGS feasibility와 QC thresholds 확인"),
    ("Day 4", "Assay design", "E/E-M/M + thyroid-lineage + survival/stress panel version lock"),
    ("Day 5", "Stat plan", "primary endpoint: postoperative CTC-EMT state transition; secondary: genetic concordance"),
    ("Day 6", "Deck compression", "12 slides를 Science-track biological perturbation flow로 압축"),
    ("Day 7", "Samsung reviewer rehearsal", "expected attack table로 30분 mock review"),
]


def fig_card(name: str, slot: str, title: str, role: str) -> str:
    return f"""
    <article class="fig-card">
      <a href="{ASSET}/{name}"><img src="{ASSET}/{name}" alt="{title}" loading="lazy" /></a>
      <div class="fig-meta"><span>{slot}</span><h3>{title}</h3><p>{role}</p></div>
    </article>
    """


def impact_card(item: dict) -> str:
    return f"""
    <article class="impact-card">
      <span>{item['label']}</span>
      <h3>{item['title']}</h3>
      <p>{item['body']}</p>
      <div class="boundary"><b>Boundary</b>{item['boundary']}</div>
    </article>
    """


def attack_rows() -> str:
    rows = []
    for attack, answer, support in REVIEWER_ATTACKS:
        rows.append(f"<tr><td>{attack}</td><td>{answer}</td><td>{support}</td></tr>")
    return "\n".join(rows)


def day_rows() -> str:
    return "".join(f"<tr><td>{d}</td><td>{focus}</td><td>{task}</td></tr>" for d, focus, task in SEVEN_DAY)


def html_page() -> str:
    return f"""<!DOCTYPE html>
<html lang="ko">
<head>
<meta charset="UTF-8" />
<meta name="viewport" content="width=device-width, initial-scale=1.0" />
<title>K-Thyro CTC-EMT Breakthrough Case</title>
<link rel="preconnect" href="https://fonts.googleapis.com" />
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin />
<link href="https://fonts.googleapis.com/css2?family=Cormorant+Garamond:wght@500;600;700&family=Newsreader:ital,wght@0,400;0,500;0,600;1,400&family=JetBrains+Mono:wght@400;500;700&family=Noto+Sans+KR:wght@400;500;700;800&display=swap" rel="stylesheet" />
<style>
  :root{{--ink:#102033;--muted:#52606f;--paper:#f6f0e6;--card:#fffaf1;--line:#d4c6b3;--line2:#eadfcc;--red:#8f2d25;--blue:#244e73;--green:#426b50;--gold:#b58534;--coal:#17212f;--good:#e9f3e8;--bad:#f7dedb}}
  *{{box-sizing:border-box}} html{{scroll-behavior:smooth}}
  body{{margin:0;background:linear-gradient(135deg,#f8f2e8 0%,#efe3d1 58%,#f7efe0 100%);color:var(--ink);font-family:"Newsreader","Noto Sans KR",serif;line-height:1.6;font-size:17px}}
  a{{color:inherit;text-decoration:none}} .wrap{{max-width:1240px;margin:0 auto;padding:32px 26px 84px}}
  .top{{display:flex;gap:10px;flex-wrap:wrap;margin-bottom:22px;font-family:"JetBrains Mono",monospace;font-size:12px}} .top a{{border:1px solid var(--line);background:rgba(255,255,255,.72);border-radius:999px;padding:8px 11px}}
  .hero{{border:2px solid var(--coal);background:linear-gradient(145deg,#fffaf1 0%,#f0dfc4 100%);padding:46px;box-shadow:10px 10px 0 rgba(23,33,47,.16)}}
  .kicker{{font-family:"JetBrains Mono",monospace;letter-spacing:.17em;text-transform:uppercase;font-size:12px;color:var(--red);font-weight:900;margin-bottom:14px}}
  h1{{font-family:"Cormorant Garamond",serif;font-size:62px;line-height:1.0;margin:0 0 14px;color:var(--coal);max-width:1000px}}
  .subtitle{{font-size:23px;color:var(--muted);max-width:990px;font-style:italic}}
  .thesis{{margin-top:24px;background:#fff;border:1px solid var(--line);border-left:7px solid var(--red);border-radius:0 12px 12px 0;padding:22px;font-size:22px}} .thesis b{{color:var(--red)}}
  .layout{{display:grid;grid-template-columns:250px 1fr;gap:26px;margin-top:30px}} .toc{{position:sticky;top:18px;align-self:start;background:rgba(255,250,241,.92);border:1px solid var(--line);border-radius:12px;padding:16px}} .toc h3{{font-family:"Cormorant Garamond",serif;color:var(--red);font-size:24px;margin:0 0 10px}} .toc a{{display:block;border-bottom:1px dotted var(--line);padding:8px 0;font-size:14px}}
  section{{background:rgba(255,250,241,.9);border:1px solid var(--line);border-radius:14px;padding:28px;margin-bottom:24px}} h2{{font-family:"Cormorant Garamond",serif;font-size:42px;line-height:1.05;margin:0 0 12px;color:var(--coal)}} .section-note{{color:var(--muted);font-style:italic;margin:0 0 18px}}
  .impact-grid{{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:14px}} .impact-card{{background:#fff;border:1px solid var(--line2);border-radius:12px;padding:18px}} .impact-card>span{{font-family:"JetBrains Mono",monospace;font-size:10px;background:var(--red);color:#fff;border-radius:999px;padding:5px 8px}} .impact-card h3{{font-family:"Cormorant Garamond",serif;font-size:30px;line-height:1.05;margin:14px 0 8px;color:var(--coal)}} .impact-card p{{margin:0;color:var(--ink)}} .boundary{{margin-top:14px;background:#f8ead7;border:1px solid var(--line2);border-radius:8px;padding:10px;font-size:14px}} .boundary b{{font-family:"JetBrains Mono",monospace;font-size:10px;color:var(--red);letter-spacing:.08em;text-transform:uppercase;display:block;margin-bottom:4px}}
  .fig-grid{{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:16px}} .fig-card{{background:#fff;border:1px solid var(--line2);border-radius:12px;overflow:hidden}} .fig-card img{{display:block;width:100%;height:auto;background:#f4eadb}} .fig-meta{{padding:14px}} .fig-meta span{{font-family:"JetBrains Mono",monospace;font-size:10px;background:var(--coal);color:#fff;border-radius:999px;padding:5px 8px}} .fig-meta h3{{font-family:"Cormorant Garamond",serif;font-size:28px;margin:12px 0 6px;color:var(--coal)}} .fig-meta p{{margin:0;color:var(--muted)}}
  .claim-ladder{{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));gap:10px}} .claim-ladder div{{background:#fff;border:1px solid var(--line2);border-radius:12px;padding:15px;position:relative}} .claim-ladder div:not(:last-child):after{{content:"→";position:absolute;right:-13px;top:40%;font-family:"JetBrains Mono",monospace;color:var(--red);font-weight:900}} .claim-ladder b{{font-family:"JetBrains Mono",monospace;font-size:10px;color:var(--red);text-transform:uppercase;letter-spacing:.08em}} .claim-ladder p{{margin:8px 0 0}}
  table{{width:100%;border-collapse:collapse;background:#fff;border:1px solid var(--line);font-size:14px}} th{{background:var(--coal);color:#fff;text-align:left;font-family:"JetBrains Mono",monospace;font-size:11px;letter-spacing:.06em;text-transform:uppercase;padding:10px}} td{{border-top:1px solid var(--line2);padding:10px;vertical-align:top}}
  .link-grid{{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:12px}} .link-card{{display:block;background:#fff;border:1px solid var(--line2);border-radius:12px;padding:16px}} .link-card b{{display:block;font-family:"Cormorant Garamond",serif;font-size:27px;color:var(--red)}} .link-card span{{display:block;color:var(--muted)}}
  .footer{{font-family:"JetBrains Mono",monospace;font-size:12px;color:var(--muted);text-align:center;margin-top:34px;padding-top:18px;border-top:1px solid var(--line)}}
  @media(max-width:980px){{.layout{{grid-template-columns:1fr}}.toc{{position:relative;top:auto}}.impact-grid,.fig-grid,.claim-ladder,.link-grid{{grid-template-columns:1fr}}.claim-ladder div:after{{display:none}}h1{{font-size:42px}}.hero{{padding:30px 22px}}}}
</style>
</head>
<body>
<div class="wrap">
  <nav class="top">
    <a href="kthyro_ctc_emt_all_pages_launchpad.html">All Pages</a>
    <a href="kthyro_ctc_emt_connected_reader.html">Connected Reader</a>
    <a href="kthyro_ctc_emt_manuscript_storyboard_paper1_style.html">Storyboard</a>
    <a href="kthyro_ctc_emt_visual_atlas_evidence.html">Evidence</a>
    <a href="kthyro_ctc_emt_visual_atlas_defense.html">Defense</a>
  </nav>
  <header class="hero">
    <div class="kicker">K-Thyro CTC-EMT · breakthrough case · no overclaim</div>
    <h1>From Possibility To A Fundable Science Case</h1>
    <p class="subtitle">대박은 “이미 다 발견했다”가 아니다. 대박은 리뷰어가 보자마자 biological perturbation, serial liquid biopsy, matched genetic anchor가 한 질문으로 잠겼다고 느끼게 만드는 것이다.</p>
    <div class="thesis"><b>Final thesis.</b> Thyroidectomy can be used as a controlled human perturbation to test whether PTC CTC-EMT states collapse, persist, or genetically match residual-risk biology after surgery.</div>
  </header>
  <div class="layout">
    <aside class="toc">
      <h3>Read Order</h3>
      <a href="#impact">Impact Case</a>
      <a href="#ladder">Claim Ladder</a>
      <a href="#figures">Killer Figures</a>
      <a href="#attacks">Reviewer Attacks</a>
      <a href="#seven">7-Day Push</a>
      <a href="#links">Next Pages</a>
    </aside>
    <main>
      <section id="impact">
        <h2>Impact Case</h2>
        <p class="section-note">가능성에서 멈추지 않게, 왜 큰 질문인지와 어디서 멈춰야 하는지를 동시에 고정한다.</p>
        <div class="impact-grid">{''.join(impact_card(x) for x in BREAKTHROUGH)}</div>
      </section>
      <section id="ladder">
        <h2>Claim Ladder</h2>
        <p class="section-note">리뷰어가 허용할 수 있는 claim만 계단식으로 올린다.</p>
        <div class="claim-ladder">
          <div><b>Public result</b><p>PTC tissue state heterogeneity exists beyond driver mutation.</p></div>
          <div><b>Design inference</b><p>Blood signals need phenotype plus genetic anchoring.</p></div>
          <div><b>Science question</b><p>Does thyroidectomy perturb CTC-EMT state?</p></div>
          <div><b>Prospective test</b><p>Serial CTC/cfDNA/tissue NGS cohort tests residual-risk biology.</p></div>
        </div>
      </section>
      <section id="figures">
        <h2>Killer Figure Set</h2>
        <p class="section-note">이 6개만 연결되면 proposal spine이 선다.</p>
        <div class="fig-grid">{''.join(fig_card(*x) for x in KILLER_FIGS)}</div>
      </section>
      <section id="attacks">
        <h2>Reviewer Attacks</h2>
        <p class="section-note">삼성 Science reviewer가 때릴 지점과 답변.</p>
        <table><thead><tr><th>Attack</th><th>Answer</th><th>Support</th></tr></thead><tbody>{attack_rows()}</tbody></table>
      </section>
      <section id="seven">
        <h2>7-Day Push</h2>
        <p class="section-note">대박은 슬로건이 아니라 7일 안에 데이터/assay/claim을 잠그는 것이다.</p>
        <table><thead><tr><th>Day</th><th>Focus</th><th>Concrete output</th></tr></thead><tbody>{day_rows()}</tbody></table>
      </section>
      <section id="links">
        <h2>Next Pages</h2>
        <div class="link-grid">
          <a class="link-card" href="kthyro_ctc_emt_all_pages_launchpad.html"><b>All Pages</b><span>전체 포털</span></a>
          <a class="link-card" href="kthyro_ctc_emt_connected_reader.html"><b>Connected Reader</b><span>쉽게 따라가는 읽기 경로</span></a>
          <a class="link-card" href="kthyro_ctc_emt_manuscript_storyboard_paper1_style.html"><b>Storyboard</b><span>Figure 1-6 논문 흐름</span></a>
        </div>
      </section>
    </main>
  </div>
  <div class="footer">Breakthrough page keeps the claim ambitious but bounded: Science-track perturbation proposal, not unsupported CTC discovery.</div>
</div>
</body>
</html>"""


def report_md() -> str:
    lines = [
        "# Breakthrough Case",
        "",
        "Purpose: 가능성에서 멈추지 않고, 삼성 Science reviewer가 납득할 수 있는 fundable biological perturbation case로 압축한다.",
        "",
        "Final thesis: Thyroidectomy can be used as a controlled human perturbation to test whether PTC CTC-EMT states collapse, persist, or genetically match residual-risk biology after surgery.",
        "",
        "Core boundary: current public data supports prior/framework logic only. Prospective Yu/hospital serial CTC/cfDNA/tissue NGS data are required for the actual transition discovery.",
        "",
        "## Impact pillars",
    ]
    for item in BREAKTHROUGH:
        lines.extend([f"- {item['label']}: {item['title']} — {item['body']} Boundary: {item['boundary']}"])
    lines.extend(["", "## Killer figures"])
    for fig, slot, title, role in KILLER_FIGS:
        lines.append(f"- {slot}: {title} ({fig}) — {role}")
    lines.extend(["", "## Seven-day push"])
    for d, focus, task in SEVEN_DAY:
        lines.append(f"- {d}: {focus} — {task}")
    return "\n".join(lines)


def update_index_handoff_launchpad() -> None:
    idx = atlas.HUB / "index.html"
    if idx.exists():
        text = idx.read_text(encoding="utf-8")
        if HTML_NAME not in text:
            needle = '    <a href="kthyro_ctc_emt_all_pages_launchpad.html"'
            insert = f'    <a href="{HTML_NAME}" style="background:linear-gradient(135deg,#8f2d25 0%,#f0dfc4 100%);border:3px solid #102033;color:#102033;font-weight:1000;font-size:18px;padding:16px 20px">★★ CTC-EMT BREAKTHROUGH CASE · from possibility to fundable science</a>\n'
            pos = text.find(needle)
            text = text[:pos] + insert + text[pos:] if pos >= 0 else text + "\n" + insert
            idx.write_text(text, encoding="utf-8")
        shutil.copy2(idx, atlas.LIVE / "index.html")

    for path in [atlas.REP / "VISUAL_WEB_HANDOFF_KR.md", atlas.HUB / "VISUAL_WEB_HANDOFF_KR.md"]:
        if not path.exists():
            continue
        text = path.read_text(encoding="utf-8")
        line = "- Breakthrough Case: http://40.82.129.113:8012/papers_hub_2026_05_04/kthyro_ctc_emt_breakthrough_case.html\n"
        if HTML_NAME not in text:
            text = text.replace("## 최우선 링크\n\n", "## 최우선 링크\n\n" + line)
            path.write_text(text, encoding="utf-8")
    handoff = atlas.HUB / "VISUAL_WEB_HANDOFF_KR.md"
    if handoff.exists():
        shutil.copy2(handoff, atlas.LIVE / "VISUAL_WEB_HANDOFF_KR.md")

    launch = atlas.HUB / "kthyro_ctc_emt_all_pages_launchpad.html"
    if launch.exists() and HTML_NAME not in launch.read_text(encoding="utf-8"):
        text = launch.read_text(encoding="utf-8")
        insert = f'<a href="{HTML_NAME}"><b>대박 논리</b><span>Breakthrough Case</span></a>'
        text = text.replace('<div class="route-grid">', '<div class="route-grid">' + insert, 1)
        launch.write_text(text, encoding="utf-8")
        shutil.copy2(launch, atlas.LIVE / "kthyro_ctc_emt_all_pages_launchpad.html")


def main() -> None:
    src = atlas.HUB / HTML_NAME
    live = atlas.LIVE / HTML_NAME
    local = atlas.REP / HTML_NAME
    report = atlas.REP / REPORT_NAME
    hub_report = atlas.HUB / REPORT_NAME

    text = html_page()
    src.write_text(text, encoding="utf-8")
    local.write_text(text, encoding="utf-8")
    shutil.copy2(src, live)

    md = report_md()
    report.write_text(md, encoding="utf-8")
    hub_report.write_text(md, encoding="utf-8")
    shutil.copy2(hub_report, atlas.LIVE / REPORT_NAME)

    update_index_handoff_launchpad()
    print(f"HTML={src}")
    print(f"LIVE={live}")
    print(f"REPORT={report}")


if __name__ == "__main__":
    main()
