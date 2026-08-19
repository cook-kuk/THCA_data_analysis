from __future__ import annotations

import csv
import html
import shutil
from pathlib import Path


ROOT = Path("/home/seungho/personal/THCA_data_analysis")
RESULT = ROOT / "project/results/cell_therapy_caption_atlas_2026_05_11"
HUB = ROOT / "project/papers_hub_2026_05_04"
LIVE_HUB = Path("/var/www/papers/papers_hub_2026_05_04")
ASSET_DIR = HUB / "assets/cell_therapy_caption_atlas_2026_05_11"
LIVE_ASSET_DIR = LIVE_HUB / "assets/cell_therapy_caption_atlas_2026_05_11"


SECTIONS = [
    {
        "section": "CAR-T",
        "anchor": "car_t",
        "theme": "Response durability, signaling, manufacturing, lineage state, and solid-tumor engineering",
        "accessions": ["GSE290722", "GSE313971", "GSE308384", "GSE236468", "GSE304795", "GSE312384", "GSE292859", "GSE284026"],
        "figures": [
            ("Fig 1", "CAR-T response durability landscape", "Longitudinal CAR-T sampling across leukapheresis, early post-infusion, and later follow-up time points; contrasts durable responders against early relapsers and highlights the native cytotoxic T-cell repertoire as a companion state variable."),
            ("Fig 2", "CAR-T tonic signaling and binding-domain redesign", "Construct-level comparison of tonic signaling, antigen sensitivity, and fitness across CD22-CAR designs; useful for showing why binding-domain architecture can rescue performance without changing the target antigen."),
            ("Fig 3", "CAR-T manufacturing and expansion QC", "Process panel linking nanoparticle-assisted expansion, cell-yield, phenotype preservation, and anti-lymphoma potency; ideal for demonstrating how manufacturing changes propagate into a deployable product profile."),
            ("Fig 4", "CAR-T translation and persistence in solid tumors", "Translation-control / persistence panel that frames pTα or related translational rewiring as a route to solid-tumor activity rather than a pure in vitro biology finding."),
            ("Fig 5", "CAR-T lineage state and exhaustion rescue", "Chromatin-state / lineage-state view centered on BACH2 and exhaustion-related transitions; good for arguing that transcriptional brakes and memory skew explain potency differences."),
            ("Fig 6", "Armored CAR-T payload localization", "Tumor-restricted payload delivery schematic showing how armored CAR design can raise local cytokine exposure while keeping systemic toxicity lower than blanket secretion."),
        ],
        "tables": [
            ("Table 1", "CAR-T candidate map", "Accession, modality, score, best task, and practical caveat for every CAR-T-related GEO record in the bundle."),
            ("Table 2", "CAR-T validation corridor", "Pairs each CAR-T anchor with the best external or orthogonal validation cohort, so the reader can see which claims are strong versus still mechanistic."),
            ("Table 3", "CAR-T risk and boundary table", "Summarizes what must not be overclaimed: mouse-only findings, process-only signals, construct-specific effects, and antigen-low model dependence."),
        ],
    },
    {
        "section": "CAR-NKT",
        "anchor": "car_nkt",
        "theme": "Off-the-shelf design, responder phenotype, and bridge-to-clinic logic",
        "accessions": ["GSE270430", "GSE223071"],
        "figures": [
            ("Fig 1", "Allogeneic CAR-NKT design", "Architecture showing how HSPC-engineered CD33-targeted CAR-NKT cells are positioned as an off-the-shelf alternative to bespoke CAR-T manufacturing."),
            ("Fig 2", "CAR-NKT plus hypomethylating-agent synergy", "Response panel that combines target engagement, HMA co-treatment, and anti-malignancy activity to show the translation hook rather than only the antigen construct."),
            ("Fig 3", "Responder phenotype from CAR-NKT clinical series", "Clinical readout from the neuroblastoma CAR-NKT series highlighting which responder states are most useful for validation and why the platform is still thin."),
            ("Fig 4", "CAR-NKT translational corridor", "Cross-disease bridge figure showing why AML/MDS and neuroblastoma should be treated as adjacent validation lanes rather than one shared proof."),
        ],
        "tables": [
            ("Table 1", "CAR-NKT validation sources", "Accession-linked list of the primary CAR-NKT dataset and the adjacent clinical series used to support the platform hypothesis."),
            ("Table 2", "CAR-NKT caveat table", "Notes where the evidence is strong enough for a platform pitch and where it still depends on one disease setting."),
        ],
    },
    {
        "section": "TIL",
        "anchor": "til",
        "theme": "Exhaustion rescue, product potency, organoid/xenograft support",
        "accessions": ["GSE303442"],
        "figures": [
            ("Fig 1", "TIL exhaustion rescue via RELB", "Mechanistic panel showing how RELB reprograms exhausted TILs toward a more therapy-competent state."),
            ("Fig 2", "TIL product potency and expansion", "Product-quality panel that links rescue to expansion, cytotoxicity, and the kind of potency metrics used for adoptive cell therapy release decisions."),
            ("Fig 3", "Organoid and xenograft support", "Orthogonal validation panel that moves the claim from expression change toward a more convincing functional readout."),
            ("Fig 4", "TIL manufacturing and persistence logic", "A practical QC figure that separates product expansion, persistence, and re-infusion readiness from direct clinical efficacy."),
        ],
        "tables": [
            ("Table 1", "TIL validation anchors", "Lists the primary TIL GEO study and the orthogonal support classes that strengthen the productization story."),
            ("Table 2", "TIL boundary table", "Keeps the reader from overcalling a manufacturing/potency result as a proven outcome effect."),
        ],
    },
    {
        "section": "Blinatumomab",
        "anchor": "blinatumomab",
        "theme": "Dose scheduling, T-cell dysfunction, memory shift, and cytokine-risk management",
        "accessions": ["GSE294241", "GSE292621"],
        "figures": [
            ("Fig 1", "Blinatumomab response durability", "Longitudinal patient sampling across therapy windows showing how T-cell dysfunction evolves during blinatumomab exposure."),
            ("Fig 2", "Re-dosing drives TCF7-high central memory state", "Single-cell multi-omics panel showing how re-dosing shifts CD8 T cells toward a more central-memory phenotype with reduced cytokine output."),
            ("Fig 3", "Cytokine-risk and efficacy tradeoff", "Side-by-side panel that links improved memory state to lower cytokine production, making the scheduling problem concrete rather than purely descriptive."),
            ("Fig 4", "Dose window and sample timing", "A practical clinical figure that lays out where pre/post samples should land if the goal is response prediction rather than a generic biology readout."),
        ],
        "tables": [
            ("Table 1", "Blinatumomab biomarkers and scheduling cues", "Compiles the response-state markers most relevant to who should be re-dosed and when."),
            ("Table 2", "Blinatumomab safety boundary table", "Highlights why cytokine-risk biology must be separated from efficacy claims."),
        ],
    },
]


def esc(s: object) -> str:
    return html.escape(str(s), quote=True)


def ensure_dirs() -> None:
    RESULT.mkdir(parents=True, exist_ok=True)
    ASSET_DIR.mkdir(parents=True, exist_ok=True)
    LIVE_ASSET_DIR.mkdir(parents=True, exist_ok=True)
    HUB.mkdir(parents=True, exist_ok=True)
    LIVE_HUB.mkdir(parents=True, exist_ok=True)


def write_tsv(path: Path) -> None:
    headers = ["section", "type", "label", "title", "description", "accessions"]
    rows = []
    for sec in SECTIONS:
        for label, title, desc in sec["figures"]:
            rows.append({"section": sec["section"], "type": "figure", "label": label, "title": title, "description": desc, "accessions": ", ".join(sec["accessions"])})
        for label, title, desc in sec["tables"]:
            rows.append({"section": sec["section"], "type": "table", "label": label, "title": title, "description": desc, "accessions": ", ".join(sec["accessions"])})
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=headers, delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def write_summary(path: Path) -> None:
    lines = [
        "# Cell Therapy Caption Atlas",
        "",
        "Date: 2026-05-11",
        "",
        "This atlas consolidates CAR-T, CAR-NKT, TIL, and blinatumomab into a sendable caption bank for figures and tables.",
        "",
        "Snapshot:",
        f"- Sections: {len(SECTIONS)}",
        f"- Figures: {sum(len(s['figures']) for s in SECTIONS)}",
        f"- Tables: {sum(len(s['tables']) for s in SECTIONS)}",
        "",
        "Use:",
        "- Copy these captions into a report, email, or manuscript planning note.",
        "- The figures are explanation templates, not claims of new wet-lab results.",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def modality_cards() -> str:
    return "".join(
        f"<a class='card' href='#{esc(sec['anchor'])}'><h3>{esc(sec['section'])}</h3><p>{esc(sec['theme'])}</p><p class='small'>{len(sec['figures'])} figure captions · {len(sec['tables'])} table captions</p></a>"
        for sec in SECTIONS
    )


def section_html(sec: dict) -> str:
    fig_cards = "".join(
        f"<div class='caption-card'><div class='pill figure'>{esc(label)}</div><h3>{esc(title)}</h3><p>{esc(desc)}</p></div>"
        for label, title, desc in sec["figures"]
    )
    table_cards = "".join(
        f"<div class='caption-card'><div class='pill table'>{esc(label)}</div><h3>{esc(title)}</h3><p>{esc(desc)}</p></div>"
        for label, title, desc in sec["tables"]
    )
    acc_list = "".join(f"<li>{esc(a)}</li>" for a in sec["accessions"])
    return f"""
<section id="{esc(sec['anchor'])}">
  <h2><span class="num">{esc(sec['section'])}</span> {esc(sec['section'])}</h2>
  <div class="sub">{esc(sec['theme'])}</div>
  <div class="box good"><b>Anchor accessions:</b> {', '.join(esc(a) for a in sec['accessions'])}</div>
  <div class="grid">{fig_cards}</div>
  <h3 class="subhead">Table captions</h3>
  <div class="grid">{table_cards}</div>
  <div class="box warn"><b>Boundary.</b> {esc('Keep the claims in this section aligned with the anchor accessions and avoid generalizing beyond the validation lane.')}</div>
  <div class="smallwrap"><ul>{acc_list}</ul></div>
</section>
"""


def build_html() -> str:
    cards = modality_cards()
    sections = "".join(section_html(sec) for sec in SECTIONS)
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8" />
<meta name="viewport" content="width=device-width,initial-scale=1.0" />
<title>Cell Therapy Caption Atlas · CAR-T / CAR-NKT / TIL / Blinatumomab</title>
<style>
:root{{--bg:#061019;--panel:#0f1726;--line:#22324b;--ink:#edf3fb;--muted:#9aabc0;--gold:#ffd28a;--teal:#35d39d;--blue:#7eb6ff;--red:#ff8a6b;--violet:#c79cff}}
*{{box-sizing:border-box}}
html,body{{margin:0;background:var(--bg);color:var(--ink);font:14px/1.55 Inter,"Noto Sans KR",sans-serif}}
a{{color:inherit;text-decoration:none}}
h1,h2,h3{{font-family:"Cormorant Garamond","Newsreader",serif;letter-spacing:-.02em;color:#fff7dc}}
.hero{{padding:64px 34px 38px;background:radial-gradient(circle at top left,#172946 0%,#0b1220 60%,#050810 100%);border-bottom:1px solid var(--line)}}
.hero-inner{{max-width:1500px;margin:0 auto}}
.kicker{{font:700 11px "JetBrains Mono",monospace;color:var(--gold);letter-spacing:.24em;text-transform:uppercase}}
h1{{font-size:64px;line-height:.96;margin:12px 0}}
.lead{{max-width:1100px;color:#cad5e5;font:18px/1.65 "Newsreader",serif}}
.stats{{display:grid;grid-template-columns:repeat(5,1fr);gap:10px;margin-top:24px}}
.stat{{background:rgba(255,255,255,.04);border:1px solid rgba(255,210,138,.18);padding:13px 14px;border-radius:10px}}
.stat b{{display:block;color:var(--gold);font:700 28px "Cormorant Garamond",serif;line-height:1}}
.stat span{{display:block;color:var(--muted);font:700 10px "JetBrains Mono",monospace;text-transform:uppercase;letter-spacing:.08em;margin-top:4px}}
.wrap{{max-width:1500px;margin:0 auto;padding:0 34px 80px}}
section{{padding:30px 0;border-bottom:1px solid var(--line)}}
h2{{font-size:36px;margin:0 0 6px}}
h2 .num{{font:700 13px "JetBrains Mono",monospace;color:var(--gold);letter-spacing:.18em;margin-right:10px}}
.sub{{font:700 10px "JetBrains Mono",monospace;color:var(--muted);letter-spacing:.12em;text-transform:uppercase;margin-bottom:16px}}
.subhead{{font:700 12px "JetBrains Mono",monospace;color:var(--gold);letter-spacing:.12em;margin:26px 0 12px;text-transform:uppercase}}
.grid{{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:14px}}
.card{{display:block;background:var(--panel);border:1px solid var(--line);border-radius:12px;padding:16px 18px;min-height:120px}}
.card:hover{{border-color:var(--gold)}}
.card p{{margin:8px 0 0}}
.caption-card{{background:var(--panel);border:1px solid var(--line);border-radius:12px;padding:14px 16px}}
.caption-card h3{{margin:10px 0 6px}}
.pill{{display:inline-block;padding:3px 9px;border-radius:999px;font:700 10px "JetBrains Mono",monospace;text-transform:uppercase;letter-spacing:.08em}}
.pill.figure{{background:rgba(53,211,157,.12);color:var(--teal);border:1px solid rgba(53,211,157,.22)}}
.pill.table{{background:rgba(126,182,255,.12);color:var(--blue);border:1px solid rgba(126,182,255,.22)}}
.box{{background:var(--panel);border:1px solid var(--line);border-left:4px solid var(--gold);border-radius:0 8px 8px 0;padding:14px 16px;margin:14px 0}}
.box.good{{border-left-color:var(--teal);background:#0a1a16}}
.box.warn{{border-left-color:var(--red);background:#1d1310}}
.small{{color:var(--muted);font-size:12px}}
.smallwrap ul{{margin:10px 0 0 18px}}
.smallwrap li{{margin:4px 0}}
@media(max-width:1100px){{.grid{{grid-template-columns:1fr 1fr}}}}
@media(max-width:980px){{.stats{{grid-template-columns:repeat(2,1fr)}}.grid{{grid-template-columns:1fr}}h1{{font-size:40px}}.wrap{{padding:0 20px 60px}}}}
</style>
</head>
<body>
<header class="hero">
  <div class="hero-inner">
    <div class="kicker">Cell therapy caption atlas · 2026-05-11</div>
    <h1>CAR-T / CAR-NKT / TIL /<br/><em>Blinatumomab caption bank</em></h1>
    <p class="lead">
      This is the sendable figure-and-table explanation package.
      It gives you caption-ready language for CAR-T, CAR-NKT, TIL, and blinatumomab in one place.
    </p>
    <div style="margin-top:18px;color:#9aabc0;font:12px 'JetBrains Mono',monospace;">
      <a href="index.html" style="color:#35d39d">Hub index</a> ·
      <a href="car_t_bundle_2026_05_11.html" style="color:#35d39d">CAR-T bundle</a> ·
      <a href="car_nkt_platform_2026_05_11.html" style="color:#7eb6ff">CAR-NKT</a> ·
      <a href="therapy_spectrum_topic_report_2026_05_11.html" style="color:#c79cff">topic report</a>
    </div>
    <div class="stats">
      <div class="stat"><b>{len(SECTIONS)}</b><span>Modalities</span></div>
      <div class="stat"><b>{sum(len(s['figures']) for s in SECTIONS)}</b><span>Figure captions</span></div>
      <div class="stat"><b>{sum(len(s['tables']) for s in SECTIONS)}</b><span>Table captions</span></div>
      <div class="stat"><b>Sendable</b><span>One-page bank</span></div>
      <div class="stat"><b>GEO</b><span>Public anchors</span></div>
    </div>
  </div>
</header>
<div class="wrap">
<section>
  <h2><span class="num">01</span>Overview</h2>
  <div class="sub">What is included</div>
  <div class="grid">{cards}</div>
  <div class="box good"><b>Use case:</b> if you need a single page to send or paste into a brief, this is it.</div>
</section>
{sections}
</div>
</body>
</html>"""


def main() -> None:
    ensure_dirs()
    tsv = RESULT / "cell_therapy_caption_atlas.tsv"
    summary = RESULT / "SUMMARY.md"
    html_path = HUB / "cell_therapy_caption_atlas_2026_05_11.html"

    write_tsv(tsv)
    write_summary(summary)
    html_path.write_text(build_html(), encoding="utf-8")

    shutil.copy2(tsv, ASSET_DIR / tsv.name)
    shutil.copy2(summary, ASSET_DIR / summary.name)
    shutil.copy2(tsv, LIVE_ASSET_DIR / tsv.name)
    shutil.copy2(summary, LIVE_ASSET_DIR / summary.name)
    shutil.copy2(html_path, LIVE_HUB / html_path.name)


if __name__ == "__main__":
    main()
