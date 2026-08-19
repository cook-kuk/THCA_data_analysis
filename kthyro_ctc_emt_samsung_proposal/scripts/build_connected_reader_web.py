#!/usr/bin/env python3
from __future__ import annotations

import shutil
from pathlib import Path

import build_paper1_style_visual_atlas as atlas
import build_paper1_style_visual_atlas_split as split


HTML_NAME = "kthyro_ctc_emt_connected_reader.html"
REPORT_NAME = "CONNECTED_READER_GUIDE_KR.md"


def asset(name: str) -> str:
    return f"assets/{atlas.ASSET}/{name}"


READING_STEPS = [
    {
        "n": "01",
        "title": "한 문장으로 먼저 고정",
        "question": "이 프로젝트가 무엇인가?",
        "answer": "갑상선절제술을 controlled human perturbation으로 보고, 수술 전후 CTC-EMT 상태와 matched genetic map을 읽는 Science-track 제안이다.",
        "open": "kthyro_ctc_emt_manuscript_storyboard_paper1_style.html",
        "open_label": "Manuscript Storyboard",
        "fig": "PF11_manuscript_blueprint_board.png",
        "say": "현재 public-data output은 tissue-state / marker / NGS-anchor prior framework다.",
        "dont": "AI platform, drug platform, broad spatial platform으로 말하지 않는다.",
    },
    {
        "n": "02",
        "title": "왜 mutation만으로 부족한가",
        "question": "왜 CTC phenotype과 NGS를 같이 봐야 하나?",
        "answer": "TCGA public pilot에서 BRAF-mutant PTC도 하나의 상태로 접히지 않는다. driver mutation은 필요하지만 residual-risk state를 단독으로 설명하지 못한다.",
        "open": "kthyro_ctc_emt_visual_atlas_evidence.html#figures",
        "open_label": "Evidence Figures",
        "fig": "F03_tcga_braf_state_split.png",
        "say": "BRAF-mutant tumors split across multiple state labels; largest BRAF label is 33.2%.",
        "dont": "이 결과를 CTC biology나 recurrence prediction으로 넘기지 않는다.",
    },
    {
        "n": "03",
        "title": "조직 상태가 실제로 조직화되어 있는가",
        "question": "EMT/lineage/stress marker를 왜 tissue context에서 골랐나?",
        "answer": "공개 spatial data는 tissue-state organization을 보여준다. 이는 CTC origin proof가 아니라 marker/ROI logic을 강화하는 근거다.",
        "open": "kthyro_ctc_emt_visual_atlas_evidence.html#F05",
        "open_label": "Spatial Evidence",
        "fig": "F05_spatial_tissue_state_organization.png",
        "say": "Spatial organization supports marker and ROI selection.",
        "dont": "Spatial hotspot이 곧 CTC shedding source라고 말하지 않는다.",
    },
    {
        "n": "04",
        "title": "CTC-EMT-NGS panel로 연결",
        "question": "어떤 marker와 NGS 축으로 미래 cohort를 읽을 것인가?",
        "answer": "epithelial, hybrid/mesenchymal, thyroid-lineage, survival/stemness, immune/stress module을 matched tumor-normal NGS와 연결한다.",
        "open": "kthyro_ctc_emt_visual_atlas_evidence.html#F08",
        "open_label": "Panel Figure",
        "fig": "F08_ctc_emt_ngs_prior_panel.png",
        "say": "Phenotype alone is nonspecific; genetic anchoring is required.",
        "dont": "CTC/cfDNA NGS가 모든 early PTC에서 성공한다고 보장하지 않는다.",
    },
    {
        "n": "05",
        "title": "왜 병원/Yu 데이터가 필요한가",
        "question": "public data로 왜 끝낼 수 없는가?",
        "answer": "serial pre/post thyroidectomy blood, CTC EMT phenotype, matched tissue-normal NGS, cfDNA/CTC-enriched NGS, follow-up을 동시에 가진 public dataset이 없다.",
        "open": "kthyro_ctc_emt_visual_atlas_defense.html#PF04",
        "open_label": "Dataset Gap",
        "fig": "PF04_dataset_gap_bridge.png",
        "say": "The missing joint dataset is the proposal justification.",
        "dont": "missing dataset 분석을 prospective discovery 결과처럼 쓰지 않는다.",
    },
    {
        "n": "06",
        "title": "교수님/병원에 무엇을 요청할지",
        "question": "다음 7일에 무엇을 받아야 하나?",
        "answer": "serial CTC raw phenotype, marker intensity, matched tissue-normal NGS fields, cfDNA availability, pathology/follow-up CRF를 요청한다.",
        "open": "kthyro_ctc_emt_visual_atlas_tables.html#crf",
        "open_label": "CRF Table",
        "fig": "PF12_future_data_unlock_map.png",
        "say": "Data request is figure-driven: each requested field unlocks a specific future figure.",
        "dont": "불필요한 broad spatial/drug/AI platform data request로 넓히지 않는다.",
    },
]


def step_card(step: dict) -> str:
    return f"""
    <section id="step-{step['n']}">
      <div class="step-head">
        <span>{step['n']}</span>
        <div>
          <h2>{step['title']}</h2>
          <p class="section-note">{step['question']}</p>
        </div>
      </div>
      <div class="connected-row">
        <a class="hero-fig" href="{asset(step['fig'])}">
          <img src="{asset(step['fig'])}" alt="{step['title']}" loading="lazy" />
        </a>
        <div class="explain-panel">
          <h3>쉽게 말하면</h3>
          <p>{step['answer']}</p>
          <div class="mini-grid">
            <div class="ok"><b>Say</b><p>{step['say']}</p></div>
            <div class="no"><b>Do not say</b><p>{step['dont']}</p></div>
          </div>
          <a class="next-button" href="{step['open']}">Open {step['open_label']}</a>
        </div>
      </div>
    </section>
    """


def quick_paths() -> str:
    return """
    <section id="paths">
      <h2>Three Reading Modes</h2>
      <p class="section-note">상황별로 어디부터 볼지 정한다. 교수님께는 5분 경로가 가장 안전하다.</p>
      <div class="grid g3">
        <div class="card"><h4>5 min</h4><div class="big">1→5→6</div><p class="small">정체성, missing dataset, 병원 요청만 빠르게 확인.</p></div>
        <div class="card"><h4>15 min</h4><div class="big">1→2→3→4→5</div><p class="small">논문 흐름과 public-data 근거를 순서대로 확인.</p></div>
        <div class="card"><h4>Review</h4><div class="big">2→5→Defense</div><p class="small">심사위원이 공격할 overclaim boundary부터 점검.</p></div>
      </div>
      <div class="callout good"><b>연결 규칙.</b> Public data는 “왜 이 질문이 필요한지”를 만든다. Yu/hospital serial CTC data가 들어와야 “수술 후 CTC-EMT transition”을 검정할 수 있다.</div>
    </section>
    """


def page_body() -> str:
    ladder = "".join(
        f"<a href='#step-{s['n']}'><span>{s['n']}</span>{s['title']}</a>" for s in READING_STEPS
    )
    body = f"""
    <section id="map">
      <h2>One Connected Map</h2>
      <p class="section-note">흩어진 그림을 하나의 논문/제안서 논리로 연결한다.</p>
      <div class="path-ladder">{ladder}</div>
    </section>
    {quick_paths()}
    {''.join(step_card(s) for s in READING_STEPS)}
    <section id="handoff">
      <h2>Send This First</h2>
      <div class="grid g3">
        <a class="card" href="kthyro_ctc_emt_manuscript_storyboard_paper1_style.html"><h4>Manuscript Storyboard</h4><p>논문 흐름과 Figure 1-6.</p></a>
        <a class="card" href="kthyro_ctc_emt_visual_atlas_split_hub.html"><h4>Split Visual Atlas</h4><p>그림을 섹션별로 크게 보기.</p></a>
        <a class="card" href="VISUAL_WEB_HANDOFF_KR.md"><h4>Handoff Memo</h4><p>교수님께 보낼 링크 모음.</p></a>
      </div>
    </section>
    """
    return body


def html_page() -> str:
    toc = (
        "<a href='#map'>Connected map</a>"
        "<a href='#paths'>Reading modes</a>"
        + "".join(f"<a href='#step-{s['n']}'>{s['n']}. {s['title']}</a>" for s in READING_STEPS)
        + "<a href='#handoff'>Send first</a>"
    )
    page = split.page(
        "Connected Reader",
        "흩어진 그림과 표를 하나의 읽기 경로로 묶었다. 처음 보는 사람도 5분 안에 왜 Science-track이고, 무엇을 주장할 수 있고, 병원 데이터가 왜 필요한지 따라가게 만든다.",
        "reader",
        page_body(),
        toc,
    )
    return page.replace(
        "</style>",
        """
        .path-ladder{display:grid;grid-template-columns:repeat(6,minmax(0,1fr));gap:10px}
        .path-ladder a{background:#fff;border:1px solid var(--line2);border-radius:12px;padding:14px;min-height:118px;color:var(--ink);position:relative}
        .path-ladder a:not(:last-child):after{content:"→";position:absolute;right:-13px;top:42%;font-family:"JetBrains Mono",monospace;color:var(--red);font-weight:800;z-index:2}
        .path-ladder span{display:block;font-family:"JetBrains Mono",monospace;font-size:11px;background:var(--coal);color:#fff;border-radius:999px;width:max-content;padding:5px 8px;margin-bottom:8px}
        .step-head{display:flex;gap:14px;align-items:flex-start;border-bottom:1px solid var(--line2);padding-bottom:10px;margin-bottom:16px}
        .step-head>span{font-family:"JetBrains Mono",monospace;background:var(--red);color:#fff;border-radius:999px;padding:8px 11px;font-size:13px;font-weight:800}
        .step-head h2{margin:0}
        .connected-row{display:grid;grid-template-columns:minmax(0,1.2fr) minmax(320px,.8fr);gap:18px;align-items:start}
        .hero-fig{display:block;background:#fff;border:1px solid var(--line2);border-radius:12px;padding:12px}
        .hero-fig img{display:block;width:100%;height:auto;border-radius:8px;background:#f4eadb}
        .explain-panel{background:#fff;border:1px solid var(--line2);border-radius:12px;padding:18px;position:sticky;top:18px}
        .explain-panel h3{margin-top:0}
        .mini-grid{display:grid;grid-template-columns:1fr;gap:10px;margin:16px 0}
        .mini-grid div{border:1px solid var(--line2);border-radius:9px;padding:11px}
        .mini-grid b{font-family:"JetBrains Mono",monospace;font-size:10px;letter-spacing:.08em;text-transform:uppercase;display:block;margin-bottom:4px}
        .mini-grid p{margin:0}.ok{background:var(--good)}.no{background:#f8ead7}
        .next-button{display:inline-block;border:0;background:var(--coal);color:#fff;border-radius:999px;padding:10px 14px;font-family:"JetBrains Mono",monospace;font-size:12px;font-weight:700}
        @media(max-width:1000px){.path-ladder{grid-template-columns:repeat(2,minmax(0,1fr))}.path-ladder a:after{display:none}.connected-row{grid-template-columns:1fr}.explain-panel{position:relative;top:auto}}
        </style>""",
    )


def markdown_report() -> str:
    lines = [
        "# Connected Reader Guide",
        "",
        "Purpose: 처음 보는 사람이 자료를 쉽게 따라가도록 Storyboard, Evidence, Defense, Tables, Deck을 하나의 읽기 경로로 연결한다.",
        "",
        "Core claim boundary: public data는 prior/framework이고, postoperative CTC-EMT transition discovery는 Yu/hospital serial data가 필요하다.",
        "",
    ]
    for s in READING_STEPS:
        lines.extend(
            [
                f"## {s['n']}. {s['title']}",
                f"- Question: {s['question']}",
                f"- Easy answer: {s['answer']}",
                f"- Say: {s['say']}",
                f"- Do not say: {s['dont']}",
                f"- Link: {s['open']}",
                "",
            ]
        )
    return "\n".join(lines)


def update_index() -> None:
    idx = atlas.HUB / "index.html"
    if not idx.exists():
        return
    text = idx.read_text(encoding="utf-8")
    if HTML_NAME not in text:
        needle = '    <a href="kthyro_ctc_emt_manuscript_storyboard_paper1_style.html"'
        insert = f'    <a href="{HTML_NAME}" style="background:linear-gradient(135deg,#fffaf1 0%,#f0dfc4 100%);border:2px solid #8f2d25;color:#102033;font-weight:900;font-size:16px;padding:14px 18px">★ CTC-EMT CONNECTED READER · easiest guided path</a>\n'
        pos = text.find(needle)
        text = text[:pos] + insert + text[pos:] if pos >= 0 else text + "\n" + insert
        idx.write_text(text, encoding="utf-8")
    shutil.copy2(idx, atlas.LIVE / "index.html")


def update_handoff() -> None:
    for path in [atlas.REP / "VISUAL_WEB_HANDOFF_KR.md", atlas.HUB / "VISUAL_WEB_HANDOFF_KR.md"]:
        if not path.exists():
            continue
        text = path.read_text(encoding="utf-8")
        line = "- Connected Reader: http://40.82.129.113:8012/papers_hub_2026_05_04/kthyro_ctc_emt_connected_reader.html\n"
        if "kthyro_ctc_emt_connected_reader.html" not in text:
            text = text.replace("## 최우선 링크\n\n", "## 최우선 링크\n\n" + line)
            path.write_text(text, encoding="utf-8")
    src = atlas.HUB / "VISUAL_WEB_HANDOFF_KR.md"
    if src.exists():
        shutil.copy2(src, atlas.LIVE / "VISUAL_WEB_HANDOFF_KR.md")


def main() -> None:
    src = atlas.HUB / HTML_NAME
    live = atlas.LIVE / HTML_NAME
    local = atlas.REP / HTML_NAME
    report = atlas.REP / REPORT_NAME
    hub_report = atlas.HUB / REPORT_NAME
    src.write_text(html_page(), encoding="utf-8")
    local.write_text(src.read_text(encoding="utf-8"), encoding="utf-8")
    shutil.copy2(src, live)
    report.write_text(markdown_report(), encoding="utf-8")
    hub_report.write_text(report.read_text(encoding="utf-8"), encoding="utf-8")
    shutil.copy2(hub_report, atlas.LIVE / REPORT_NAME)
    update_index()
    update_handoff()
    print(f"HTML={src}")
    print(f"LIVE={live}")
    print(f"REPORT={report}")


if __name__ == "__main__":
    main()
