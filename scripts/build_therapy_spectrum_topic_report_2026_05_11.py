from __future__ import annotations

import csv
import html
import shutil
from dataclasses import dataclass
from pathlib import Path


ROOT = Path("/home/seungho/personal/THCA_data_analysis")
SRC = ROOT / "project/results/therapy_spectrum_promising_master_2026_05_11/therapy_spectrum_promising_master.tsv"
RESULT = ROOT / "project/results/therapy_spectrum_topic_report_2026_05_11"
HUB = ROOT / "project/papers_hub_2026_05_04"
LIVE_HUB = Path("/var/www/papers/papers_hub_2026_05_04")
ASSET_DIR = HUB / "assets/therapy_spectrum_topic_report_2026_05_11"
LIVE_ASSET_DIR = LIVE_HUB / "assets/therapy_spectrum_topic_report_2026_05_11"


@dataclass(frozen=True)
class Topic:
    slug: str
    title: str
    kicker: str
    hypothesis: str
    anchor: str
    verdict: str
    summary: str
    validation: list[str]
    next_steps: list[str]
    candidates: list[str]
    caveat: str


TOPICS = [
    Topic(
        slug="car_t_durability",
        title="CAR-T Durability",
        kicker="Response longevity / native repertoire / persistence",
        hypothesis="Durable CAR-T response is driven not just by infused CAR persistence, but by the native cytotoxic T-cell repertoire and baseline immune state.",
        anchor="GSE290722",
        verdict="strong",
        summary="The strongest candidate in the set: longitudinal LBCL data, direct durable-vs-early-relapse comparison, and an obvious validation corridor from apheresis to post-infusion cohorts.",
        validation=["GSE273170", "GSE243973", "GSE246342", "GSE145007", "GSE162975", "GSE223655", "GSE158676"],
        next_steps=[
            "Build a shared response-score model across baseline and post-infusion cohorts.",
            "Test whether monocyte burden, TSCM, and native clonal expansion explain durability better than CAR peak alone.",
            "Package as a response-prediction / monitoring module rather than a CAR invention claim.",
        ],
        candidates=["GSE290722", "GSE273170", "GSE243973", "GSE246342", "GSE145007", "GSE162975", "GSE223655"],
        caveat="LBCL-centric; treat as a CAR-T response framework, not a universal cell-therapy law.",
    ),
    Topic(
        slug="mrna_lnp_platform",
        title="mRNA-LNP Platform",
        kicker="Formulation / innate response / antigen production",
        hypothesis="mRNA-LNP efficacy is jointly controlled by nucleoside chemistry, lipid composition, and local innate sensing at the injection site.",
        anchor="GSE303153",
        verdict="strong",
        summary="This is the cleanest platform-level validation lane in the set, with human and non-human primate data plus independent injection-site / RNA-stability corroboration.",
        validation=["GSE239574", "GSE98211", "GSE233059", "PMID 39191748"],
        next_steps=[
            "Rank formulation combinations by antigen output versus reactogenicity.",
            "Use injection-site IFN-β / fibroblast activation as a companion readout.",
            "Keep the report explicit that this is platform chemistry, not cancer-specific evidence.",
        ],
        candidates=["GSE303153", "GSE239574", "GSE98211", "GSE233059"],
        caveat="Platform-level evidence is strong, but the oncology translation still needs a disease-specific bridge.",
    ),
    Topic(
        slug="cdk4_6_resistance",
        title="CDK4/6 Resistance and Biomarkers",
        kicker="Late progression / endocrine resistance / immune suppression",
        hypothesis="Late progression on CDK4/6 inhibitors is governed by a reproducible resistance program that can be read out from tumor and immune compartments.",
        anchor="GSE274141",
        verdict="strong",
        summary="This is the best validation-ready oncology biomarker lane in the set: the GEO record already names two independent validation cohorts, and the surrounding literature supports the same resistance axes.",
        validation=["GSE150997", "GSE117742", "GSE223700", "GSE74391", "GSE282705"],
        next_steps=[
            "Build a cross-cohort biomarker ranker for late progression.",
            "Add PLK1 / IFN / immune suppression as candidate axes.",
            "Translate into a therapy-switch or combination decision aid.",
        ],
        candidates=["GSE274141", "GSE150997", "GSE117742", "GSE223700", "GSE74391", "GSE282705"],
        caveat="Breast-cancer specific; do not overclaim across unrelated tumor types.",
    ),
    Topic(
        slug="car_nkt_platform",
        title="CAR-NKT Platform",
        kicker="Off-the-shelf cell therapy / responder signals",
        hypothesis="Allogeneic CAR-NKT can solve manufacturing and selection pain points if the target biology and responder phenotype are matched tightly enough.",
        anchor="GSE270430",
        verdict="moderate",
        summary="Good translational hook with real responder signals, but the field remains thin, so the report should emphasize platform potential rather than broad proof.",
        validation=["GSE223071", "GSE235923"],
        next_steps=[
            "Compare responder phenotypes across AML/MDS and neuroblastoma CAR-NKT series.",
            "Extract common markers for persistence and hyporesponsiveness.",
            "Keep the BD angle tied to manufacturing and target selection, not broad efficacy claims.",
        ],
        candidates=["GSE270430", "GSE223071", "GSE235923"],
        caveat="Independent disease settings are still needed before a strong generalization claim.",
    ),
    Topic(
        slug="pdac_sdc1",
        title="SDC1 / PDAC Antibody",
        kicker="Macropinocytosis / KRAS biology / combination therapy",
        hypothesis="Anti-SDC1 blocks nutrient salvage in KRAS-driven PDAC and can be positioned as a combination-enabling therapeutic antibody.",
        anchor="GSE311267",
        verdict="moderate",
        summary="Mechanistically strong, with the right target biology behind it, but validation is more orthogonal than same-lane GEO-based.",
        validation=["PMID 30918400", "PMID 26293675", "GSE138437"],
        next_steps=[
            "Use SDC1-high / KRAS-driven PDAC as the priority indication.",
            "Treat macropinocytosis as the mechanistic readout and ADCC as the effector readout.",
            "Require a second dataset before calling it a general antibody platform.",
        ],
        candidates=["GSE311267", "GSE138437"],
        caveat="PDAC/KRAS-centric; not a broad pan-cancer antibody claim.",
    ),
    Topic(
        slug="car_engineering_potency",
        title="CAR Engineering and Potency",
        kicker="Tonic signaling / manufacturing / exhaustion rescue / product QC",
        hypothesis="CAR potency is tunable through signaling architecture, manufacturing process, and exhaustion-state engineering, and these axes can be converted into QC-like decision rules.",
        anchor="GSE313971",
        verdict="moderate",
        summary="This is a composite engineering theme that covers construct design, expansion chemistry, TIL rescue, and persistence engineering. Very strong for methods and product design, but the biological claims need careful narrowing.",
        validation=["GSE145007", "GSE158676", "GSE162975", "GSE223655", "GSE245427", "GSE103632", "GSE303442", "GSE308384", "GSE304795", "GSE312384", "GSE292859"],
        next_steps=[
            "Split the theme into construct design, manufacturing, and exhaustion rescue subpages if needed.",
            "Treat GSE313971 as the construct-design anchor and GSE308384 / GSE303442 as process/product anchors.",
            "Keep pTα, BACH2, and armored CAR-T as submodules, not one over-broad claim.",
        ],
        candidates=["GSE313971", "GSE308384", "GSE303442", "GSE304795", "GSE312384", "GSE292859", "GSE236468", "GSE291443"],
        caveat="This is the most heterogeneous topic, so the page must keep subclaims explicitly separated.",
    ),
]


def esc(s: object) -> str:
    return html.escape(str(s), quote=True)


def ensure_dirs() -> None:
    RESULT.mkdir(parents=True, exist_ok=True)
    ASSET_DIR.mkdir(parents=True, exist_ok=True)
    LIVE_ASSET_DIR.mkdir(parents=True, exist_ok=True)
    HUB.mkdir(parents=True, exist_ok=True)
    LIVE_HUB.mkdir(parents=True, exist_ok=True)


def read_promising_rows() -> dict[str, dict[str, str]]:
    rows: dict[str, dict[str, str]] = {}
    with SRC.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            rows[row["accession"]] = row
    return rows


def write_overview(rows: dict[str, dict[str, str]]) -> None:
    cards = []
    for t in TOPICS:
        cand = rows.get(t.anchor, {})
        cards.append(
            f"<a class='topic-card' href='{esc(t.slug)}.html'>"
            f"<div class='topic-kicker'>{esc(t.kicker)}</div>"
            f"<h3>{esc(t.title)}</h3>"
            f"<p class='topic-summary'>{esc(t.summary)}</p>"
            f"<p class='topic-meta'><b>Anchor:</b> {esc(t.anchor)} | <b>Verdict:</b> {esc(t.verdict)}</p>"
            f"<p class='topic-meta small'>{esc(t.caveat)}</p>"
            f"</a>"
        )
    overview = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8" />
<meta name="viewport" content="width=device-width,initial-scale=1.0" />
<title>Therapy Spectrum Topic Report · Overview</title>
<style>
:root{{--bg:#061019;--panel:#0f1726;--line:#22324b;--ink:#edf3fb;--muted:#9aabc0;--gold:#ffd28a;--teal:#35d39d;--blue:#7eb6ff;--red:#ff8a6b;--violet:#c79cff}}
*{{box-sizing:border-box}}
html,body{{margin:0;background:var(--bg);color:var(--ink);font:14px/1.55 Inter,"Noto Sans KR",sans-serif}}
a{{color:inherit;text-decoration:none}}
h1,h2,h3{{font-family:"Cormorant Garamond","Newsreader",serif;letter-spacing:-.02em;color:#fff7dc}}
.hero{{padding:64px 34px 38px;background:radial-gradient(circle at top left,#172946 0%,#0b1220 60%,#050810 100%);border-bottom:1px solid var(--line)}}
.hero-inner{{max-width:1440px;margin:0 auto}}
.kicker{{font:700 11px "JetBrains Mono",monospace;color:var(--gold);letter-spacing:.24em;text-transform:uppercase}}
h1{{font-size:64px;line-height:.96;margin:12px 0}}
.lead{{max-width:1040px;color:#cad5e5;font:18px/1.65 "Newsreader",serif}}
.stats{{display:grid;grid-template-columns:repeat(5,1fr);gap:10px;margin-top:24px}}
.stat{{background:rgba(255,255,255,.04);border:1px solid rgba(255,210,138,.18);padding:13px 14px;border-radius:10px}}
.stat b{{display:block;color:var(--gold);font:700 28px "Cormorant Garamond",serif;line-height:1}}
.stat span{{display:block;color:var(--muted);font:700 10px "JetBrains Mono",monospace;text-transform:uppercase;letter-spacing:.08em;margin-top:4px}}
.wrap{{max-width:1440px;margin:0 auto;padding:0 34px 80px}}
section{{padding:30px 0;border-bottom:1px solid var(--line)}}
h2{{font-size:36px;margin:0 0 6px}}
h2 .num{{font:700 13px "JetBrains Mono",monospace;color:var(--gold);letter-spacing:.18em;margin-right:10px}}
.grid{{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:14px}}
.topic-card{{display:block;background:var(--panel);border:1px solid var(--line);border-radius:12px;padding:16px 18px;transition:transform .12s ease,border-color .12s ease}}
.topic-card:hover{{transform:translateY(-2px);border-color:var(--gold)}}
.topic-kicker{{font:700 11px "JetBrains Mono",monospace;color:var(--gold);letter-spacing:.16em;text-transform:uppercase}}
.topic-summary{{color:#cad5e5;margin:10px 0 8px}}
.topic-meta{{color:var(--muted);margin:5px 0}}
.small{{font-size:12px}}
@media(max-width:980px){{.stats{{grid-template-columns:repeat(2,1fr)}}.grid{{grid-template-columns:1fr}}h1{{font-size:40px}}.wrap{{padding:0 20px 60px}}}}
</style>
</head>
<body>
<header class="hero">
  <div class="hero-inner">
    <div class="kicker">Topic report · 2026-05-11</div>
    <h1>Therapy Spectrum<br/><em>Topic Pages</em></h1>
    <p class="lead">One page per hypothesis family. Use this as the top-level report for validation, paper drafting, and BD triage.</p>
    <div style="margin-top:18px;color:#9aabc0;font:12px 'JetBrains Mono',monospace;">
      <a href="index.html" style="color:#35d39d">Hub index</a> ·
      <a href="therapy_spectrum_promising_master_2026_05_11.html" style="color:#c79cff">promising master</a> ·
      <a href="therapy_spectrum_validation_matrix_2026_05_11.html" style="color:#35d39d">validation matrix</a>
    </div>
    <div class="stats">
      <div class="stat"><b>{len(TOPICS)}</b><span>Topic pages</span></div>
      <div class="stat"><b>3</b><span>Strong</span></div>
      <div class="stat"><b>3</b><span>Moderate</span></div>
      <div class="stat"><b>GEO + PubMed</b><span>Evidence base</span></div>
      <div class="stat"><b>Hub</b><span>Auto-deployed</span></div>
    </div>
  </div>
</header>
<div class="wrap">
<section>
  <h2><span class="num">01</span>Topic Index</h2>
  <div class="grid">{''.join(cards)}</div>
</section>
</div>
</body>
</html>"""
    overview_path = HUB / "therapy_spectrum_topic_report_2026_05_11.html"
    overview_path.write_text(overview, encoding="utf-8")
    shutil.copy2(overview_path, LIVE_HUB / overview_path.name)
    shutil.copy2(overview_path, ASSET_DIR / overview_path.name)
    shutil.copy2(overview_path, LIVE_ASSET_DIR / overview_path.name)


def topic_page(topic: Topic) -> str:
    validation_cards = "".join(f"<li>{esc(v)}</li>" for v in topic.validation)
    candidate_cards = "".join(f"<li>{esc(c)}</li>" for c in topic.candidates)
    next_cards = "".join(f"<li>{esc(n)}</li>" for n in topic.next_steps)
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8" />
<meta name="viewport" content="width=device-width,initial-scale=1.0" />
<title>Therapy Spectrum · {esc(topic.title)}</title>
<style>
:root{{--bg:#061019;--panel:#0f1726;--line:#22324b;--ink:#edf3fb;--muted:#9aabc0;--gold:#ffd28a;--teal:#35d39d;--blue:#7eb6ff;--red:#ff8a6b}}
*{{box-sizing:border-box}}
html,body{{margin:0;background:var(--bg);color:var(--ink);font:14px/1.55 Inter,"Noto Sans KR",sans-serif}}
a{{color:var(--teal);text-decoration:none}} a:hover{{color:var(--gold);text-decoration:underline}}
h1,h2,h3{{font-family:"Cormorant Garamond","Newsreader",serif;letter-spacing:-.02em;color:#fff7dc}}
.hero{{padding:64px 34px 34px;background:radial-gradient(circle at top left,#172946 0%,#0b1220 60%,#050810 100%);border-bottom:1px solid var(--line)}}
.hero-inner{{max-width:1320px;margin:0 auto}}
.kicker{{font:700 11px "JetBrains Mono",monospace;color:var(--gold);letter-spacing:.24em;text-transform:uppercase}}
h1{{font-size:60px;line-height:.98;margin:12px 0}}
.lead{{max-width:1040px;color:#cad5e5;font:18px/1.65 "Newsreader",serif}}
.meta{{display:grid;grid-template-columns:repeat(4,1fr);gap:10px;margin-top:24px}}
.meta div{{background:rgba(255,255,255,.04);border:1px solid rgba(255,210,138,.18);padding:13px 14px;border-radius:10px}}
.meta b{{display:block;color:var(--gold);font:700 28px "Cormorant Garamond",serif;line-height:1}}
.meta span{{display:block;color:var(--muted);font:700 10px "JetBrains Mono",monospace;text-transform:uppercase;letter-spacing:.08em;margin-top:4px}}
.wrap{{max-width:1320px;margin:0 auto;padding:0 34px 80px}}
section{{padding:28px 0;border-bottom:1px solid var(--line)}}
h2{{font-size:34px;margin:0 0 6px}}
.sub{{font:700 10px "JetBrains Mono",monospace;color:var(--muted);letter-spacing:.12em;text-transform:uppercase;margin-bottom:14px}}
.box{{background:var(--panel);border:1px solid var(--line);border-left:4px solid var(--gold);border-radius:0 8px 8px 0;padding:14px 16px;margin:14px 0}}
.box.good{{border-left-color:var(--teal);background:#0a1a16}}
.box.warn{{border-left-color:var(--red);background:#1d1310}}
.grid{{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:14px}}
.card{{background:var(--panel);border:1px solid var(--line);border-radius:10px;padding:14px 16px}}
ul{{margin:10px 0 0 18px}}
li{{margin:4px 0}}
@media(max-width:980px){{.meta{{grid-template-columns:repeat(2,1fr)}}.grid{{grid-template-columns:1fr}}h1{{font-size:40px}}.wrap{{padding:0 20px 60px}}}}
</style>
</head>
<body>
<header class="hero">
  <div class="hero-inner">
    <div class="kicker">Topic page · 2026-05-11</div>
    <h1>{esc(topic.title)}<br/><em>{esc(topic.kicker)}</em></h1>
    <p class="lead">{esc(topic.hypothesis)}</p>
    <div style="margin-top:18px;color:#9aabc0;font:12px 'JetBrains Mono',monospace;">
      <a href="therapy_spectrum_topic_report_2026_05_11.html" style="color:#35d39d">Topic overview</a> ·
      <a href="therapy_spectrum_promising_master_2026_05_11.html" style="color:#c79cff">promising master</a> ·
      <a href="therapy_spectrum_validation_matrix_2026_05_11.html" style="color:#35d39d">validation matrix</a>
    </div>
    <div class="meta">
      <div><b>{esc(topic.anchor)}</b><span>Anchor GEO</span></div>
      <div><b>{esc(topic.verdict)}</b><span>Verdict</span></div>
      <div><b>{len(topic.validation)}</b><span>Validation links</span></div>
      <div><b>{len(topic.candidates)}</b><span>Candidate set</span></div>
    </div>
  </div>
</header>
<div class="wrap">
<section>
  <h2>01 Overview</h2>
  <div class="sub">Why this topic matters</div>
  <div class="box good">{esc(topic.summary)}</div>
  <div class="box warn"><b>Caveat.</b> {esc(topic.caveat)}</div>
</section>
<section>
  <h2>02 Validation</h2>
  <div class="sub">Independent / orthogonal data that strengthen the case</div>
  <div class="grid">
    <div class="card"><h3>Validation data</h3><ul>{validation_cards}</ul></div>
    <div class="card"><h3>Topic candidates</h3><ul>{candidate_cards}</ul></div>
  </div>
</section>
<section>
  <h2>03 Next Steps</h2>
  <div class="sub">What a real validation package should do next</div>
  <ul>{next_cards}</ul>
</section>
</div>
</body>
</html>"""


def main() -> None:
    ensure_dirs()
    rows = read_promising_rows()
    overview_path = HUB / "therapy_spectrum_topic_report_2026_05_11.html"
    write_overview(rows)
    for t in TOPICS:
        path = HUB / f"{t.slug}_2026_05_11.html"
        path.write_text(topic_page(t), encoding="utf-8")
        shutil.copy2(path, LIVE_HUB / path.name)
        shutil.copy2(path, ASSET_DIR / path.name)
        shutil.copy2(path, LIVE_ASSET_DIR / path.name)
    # copy overview from helper
    shutil.copy2(overview_path, LIVE_HUB / overview_path.name)


if __name__ == "__main__":
    main()
