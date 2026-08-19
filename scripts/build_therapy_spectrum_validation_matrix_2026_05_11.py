from __future__ import annotations

import csv
import html
import shutil
from collections import Counter
from pathlib import Path


ROOT = Path("/home/seungho/personal/THCA_data_analysis")
RESULT = ROOT / "project/results/therapy_spectrum_validation_matrix_2026_05_11"
HUB = ROOT / "project/papers_hub_2026_05_04"
LIVE_HUB = Path("/var/www/papers/papers_hub_2026_05_04")
ASSET_DIR = HUB / "assets/therapy_spectrum_validation_matrix_2026_05_11"
LIVE_ASSET_DIR = LIVE_HUB / "assets/therapy_spectrum_validation_matrix_2026_05_11"


ROWS = [
    {
        "accession": "GSE290722",
        "candidate": "CAR-T response durability",
        "verdict": "strong",
        "primary_evidence": "GEO index study with 32 LBCL patients; longitudinal scRNA/scTCR across leukapheresis, 4w, 6m, 12m; durable responders vs early relapse contrasted directly",
        "external_validation": "GSE273170, GSE243973, GSE246342, GSE145007, GSE162975, GSE223655, GSE158676",
        "why": "There is now a clear validation corridor from baseline immune-state cohorts to post-infusion durability cohorts, plus costimulation and checkpoint-mechanism studies.",
        "caveat": "LBCL-only; use as CAR-T response framework, not universal biology.",
        "source": "turn3search0; turn3search1; turn3search2; turn3search3",
    },
    {
        "accession": "GSE303153",
        "candidate": "mRNA-LNP vaccine platform",
        "verdict": "strong",
        "primary_evidence": "Human + cynomolgus macaque RNA-seq in a formulation/mechanism study; reduced innate immune recognition and synergies with lipid composition.",
        "external_validation": "GSE239574, GSE98211, GSE233059, PMID 39191748",
        "why": "Platform behavior is corroborated by injection-site, macaque, and in vivo RNA stability / re-adenylation studies.",
        "caveat": "Not cancer-specific; validation is platform-level rather than oncology-specific.",
        "source": "turn0search5; turn0search6; turn1search0; turn1search4; turn2search0",
    },
    {
        "accession": "GSE270430",
        "candidate": "CAR-NKT / off-the-shelf cell therapy",
        "verdict": "moderate",
        "primary_evidence": "Human AML/MDS characterization leading to allogeneic CD33 CAR-NKT plus HMA synergy.",
        "external_validation": "GSE223071, GSE235923",
        "why": "Clinical CAR-NKT responder behavior is visible in the neuroblastoma phase-1 series; the platform also sits in a wider AML/NKT transplant and relapse biology context.",
        "caveat": "Still a thin field; independent disease settings matter a lot.",
        "source": "turn1search5; turn1search3; turn2search2",
    },
    {
        "accession": "GSE274141",
        "candidate": "CDK4/6 biomarker / endocrine resistance",
        "verdict": "strong",
        "primary_evidence": "HR+/HER2- metastatic breast cancer scRNA-seq with explicit early vs late progression and marker validation.",
        "external_validation": "MDACC n=89, Korean n=61, GSE150997, GSE117742, GSE223700, GSE74391, GSE282705",
        "why": "Built-in validation cohorts are already stated in the GEO record, and multiple independent CDK4/6 resistance datasets support the same resistance/IFN/PLK1 axis.",
        "caveat": "Breast-cancer specific; do not extrapolate beyond CDK4/6/endocrine context without new data.",
        "source": "turn0search7; turn1search2; turn2search5; turn2search7; turn2search9",
    },
    {
        "accession": "GSE311267",
        "candidate": "SDC1 therapeutic antibody / PDAC",
        "verdict": "moderate",
        "primary_evidence": "Human-specific anti-SDC1 mAb RNA-seq on PDAC cells; blocks macropinocytosis, drives ADCC, and synergizes with chemotherapy/KRAS/inmunotherapy.",
        "external_validation": "PMID 30918400; PMID 26293675; PDAC macropinocytosis biology; GSE138437",
        "why": "Target biology is strongly supported, but independent GEO-scale validation is still sparse; still, the KRAS/macropinocytosis axis is well established in PDAC.",
        "caveat": "PDAC and KRAS-driven biology are central; this is not a broad pan-cancer antibody platform yet.",
        "source": "turn1search1; turn2search1; turn2search8; turn2search10; turn2search11",
    },
    {
        "accession": "GSE303442",
        "candidate": "TIL adoptive therapy / RELB",
        "verdict": "moderate",
        "primary_evidence": "Human TIL product expansion/rescue with organoid and xenograft support.",
        "external_validation": "TIL literature + clinical TIL clonotype dynamics (not GEO-specific)",
        "why": "Excellent mechanistic/productization signal, but the validation path is mostly orthogonal rather than another same-lane GEO study.",
        "caveat": "Manufacturing/potency layer, not direct clinical endpoint proof.",
        "source": "turn0search1",
    },
    {
        "accession": "GSE313971",
        "candidate": "CD22-CAR tonic signaling",
        "verdict": "moderate",
        "primary_evidence": "Human CD22-CAR design optimization with tonic signaling/fitness rationale.",
        "external_validation": "GSE145007, GSE158676, GSE162975, GSE223655, GSE245427, GSE103632",
        "why": "Good mechanistic anchor, and the broader CAR signaling / persistence / memory literature gives a strong validation corridor even if the exact CD22 construct is specific.",
        "caveat": "Construct-specific and preclinical.",
        "source": "turn3search0; turn3search4",
    },
    {
        "accession": "GSE308384",
        "candidate": "CAR-T manufacturing / T-cell expansion NPs",
        "verdict": "moderate",
        "primary_evidence": "Process paper with direct expansion/persistence readouts and in vivo anti-lymphoma activity.",
        "external_validation": "CAR manufacturing literature / prior bead-based expansion benchmarks / CD62L targeting and less-differentiated T-cell delivery concepts",
        "why": "Strong operational value, but validation is process-centric rather than disease-centric; still, it fits a broader manufacturing optimization field.",
        "caveat": "Not proof of superior clinical efficacy on its own.",
        "source": "turn0search6",
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


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    headers = [
        "accession",
        "candidate",
        "verdict",
        "primary_evidence",
        "external_validation",
        "why",
        "caveat",
        "source",
    ]
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=headers, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_summary(path: Path, rows: list[dict[str, object]]) -> None:
    counts = Counter(r["verdict"] for r in rows)
    lines = [
        "# Therapy Spectrum Validation Matrix",
        "",
        "Date: 2026-05-11",
        "",
        "This matrix summarizes how much external data exist to validate the strongest therapy-spectrum candidates.",
        "",
        "Snapshot:",
        f"- Rows: {len(rows)}",
        f"- Verdict counts: {', '.join(f'{k}={v}' for k, v in counts.items())}",
        "",
        "Interpretation:",
        "- strong = direct GEO plus a clear validation corridor",
        "- moderate = promising GEO + orthogonal/adjacent validation, but still thinner",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def table(rows: list[dict[str, object]]) -> str:
    return "\n".join(
        f"<tr>"
        f"<td><a href='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={esc(r['accession'])}'>{esc(r['accession'])}</a></td>"
        f"<td>{esc(r['candidate'])}</td>"
        f"<td><span class='tag {esc(r['verdict'])}'>{esc(r['verdict'])}</span></td>"
        f"<td>{esc(r['primary_evidence'])}</td>"
        f"<td>{esc(r['external_validation'])}</td>"
        f"<td>{esc(r['why'])}</td>"
        f"<td>{esc(r['caveat'])}</td>"
        f"<td><code>{esc(r['source'])}</code></td>"
        f"</tr>"
        for r in rows
    )


def build_html(rows: list[dict[str, object]]) -> str:
    counts = Counter(r["verdict"] for r in rows)
    strong = [r for r in rows if r["verdict"] == "strong"]
    moderate = [r for r in rows if r["verdict"] == "moderate"]
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8" />
<meta name="viewport" content="width=device-width,initial-scale=1.0" />
<title>Therapy Spectrum Validation Matrix</title>
<style>
:root{{--bg:#071019;--panel:#0f1726;--line:#22324b;--ink:#edf3fb;--muted:#9aabc0;--gold:#ffd28a;--teal:#35d39d;--blue:#7eb6ff;--red:#ff8a6b}}
*{{box-sizing:border-box}}
html,body{{margin:0;background:var(--bg);color:var(--ink);font:14px/1.55 Inter,"Noto Sans KR",sans-serif}}
a{{color:var(--teal);text-decoration:none}} a:hover{{color:var(--gold);text-decoration:underline}}
code{{font-family:"JetBrains Mono",monospace;font-size:12px}}
h1,h2,h3{{font-family:"Cormorant Garamond","Newsreader",serif;letter-spacing:-.02em;color:#fff7dc}}
.hero{{padding:64px 34px 38px;background:radial-gradient(circle at top left,#172946 0%,#0b1220 60%,#050810 100%);border-bottom:1px solid var(--line)}}
.hero-inner{{max-width:1440px;margin:0 auto}}
.kicker{{font:700 11px "JetBrains Mono",monospace;color:var(--gold);letter-spacing:.24em;text-transform:uppercase}}
h1{{font-size:62px;line-height:.96;margin:12px 0}}
.lead{{max-width:1060px;color:#cad5e5;font:18px/1.65 "Newsreader",serif}}
.stats{{display:grid;grid-template-columns:repeat(4,1fr);gap:10px;margin-top:24px}}
.stat{{background:rgba(255,255,255,.04);border:1px solid rgba(255,210,138,.18);padding:13px 14px;border-radius:10px}}
.stat b{{display:block;color:var(--gold);font:700 28px "Cormorant Garamond",serif;line-height:1}}
.stat span{{display:block;color:var(--muted);font:700 10px "JetBrains Mono",monospace;text-transform:uppercase;letter-spacing:.08em;margin-top:4px}}
.wrap{{max-width:1440px;margin:0 auto;padding:0 34px 80px}}
section{{padding:30px 0;border-bottom:1px solid var(--line)}}
h2{{font-size:36px;margin:0 0 6px}}
h2 .num{{font:700 13px "JetBrains Mono",monospace;color:var(--gold);letter-spacing:.18em;margin-right:10px}}
.sub{{font:700 10px "JetBrains Mono",monospace;color:var(--muted);letter-spacing:.12em;text-transform:uppercase;margin-bottom:16px}}
.t{{width:100%;border-collapse:collapse;margin:12px 0 16px;font-size:12px}}
.t th,.t td{{border:1px solid var(--line);padding:8px 9px;vertical-align:top}}
.t th{{background:#16243a;color:var(--gold);font:700 10px "JetBrains Mono",monospace;letter-spacing:.05em;text-transform:uppercase;text-align:left}}
.t tr:nth-child(even) td{{background:rgba(255,255,255,.02)}}
.tag{{display:inline-block;padding:2px 8px;border-radius:99px;font:700 10px "JetBrains Mono",monospace;letter-spacing:.06em;text-transform:uppercase;border:1px solid transparent}}
.tag.strong{{background:rgba(53,211,157,.12);border-color:rgba(53,211,157,.28);color:var(--teal)}}
.tag.moderate{{background:rgba(126,182,255,.12);border-color:rgba(126,182,255,.28);color:var(--blue)}}
.box{{background:var(--panel);border:1px solid var(--line);border-left:4px solid var(--gold);border-radius:0 8px 8px 0;padding:14px 16px;margin:14px 0}}
.box.good{{border-left-color:var(--teal);background:#0a1a16}}
.box.warn{{border-left-color:var(--red);background:#1d1310}}
.grid{{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:14px}}
.card{{background:var(--panel);border:1px solid var(--line);border-radius:10px;padding:14px 16px}}
@media(max-width:980px){{.stats{{grid-template-columns:repeat(2,1fr)}}.grid{{grid-template-columns:1fr}}h1{{font-size:40px}}.wrap{{padding:0 20px 60px}}}}
</style>
</head>
<body>
<header class="hero">
  <div class="hero-inner">
    <div class="kicker">Validation matrix · 2026-05-11</div>
    <h1>Therapy Spectrum<br/><em>Validation Matrix</em></h1>
    <p class="lead">
      This is the cross-check layer: which promising candidates already have independent or orthogonal data that make them believable, and which ones still need the next dataset.
    </p>
    <div style="margin-top:18px;color:#9aabc0;font:12px 'JetBrains Mono',monospace;">
      <a href="index.html" style="color:#35d39d">Hub index</a> ·
      <a href="therapy_spectrum_promising_master_2026_05_11.html" style="color:#c79cff">promising master</a> ·
      <a href="therapy_spectrum_ready_split_2026_05_11.html" style="color:#f2c46d">ready split</a>
    </div>
    <div class="stats">
      <div class="stat"><b>{len(rows)}</b><span>Rows</span></div>
      <div class="stat"><b>{counts.get('strong', 0)}</b><span>Strong</span></div>
      <div class="stat"><b>{counts.get('moderate', 0)}</b><span>Moderate</span></div>
      <div class="stat"><b>GEO + PubMed</b><span>Source basis</span></div>
    </div>
  </div>
</header>
<div class="wrap">
<section>
  <h2><span class="num">01</span>TL;DR</h2>
  <div class="sub">Which candidates are actually defensible now</div>
  <div class="grid">
    <div class="card"><h3>Strong</h3><p>{len(strong)} candidates with direct GEO plus clear validation path.</p></div>
    <div class="card"><h3>Moderate</h3><p>{len(moderate)} candidates with good primary data but thinner or more orthogonal validation.</p></div>
  </div>
  <div class="box good"><b>Best validation anchors:</b> GSE274141, GSE303153, and GSE290722.</div>
</section>

<section>
  <h2><span class="num">02</span>Full Matrix</h2>
  <div class="sub">Candidate × evidence × verdict</div>
  <table class="t">
    <thead><tr><th>Accession</th><th>Candidate</th><th>Verdict</th><th>Primary evidence</th><th>External validation</th><th>Why</th><th>Caveat</th><th>Source</th></tr></thead>
    <tbody>{table(rows)}</tbody>
  </table>
  <div class="box warn"><b>Boundary.</b> “Strong” means defensible for a paper or validation package, not clinically proven. “Moderate” still needs a cleaner second cohort or orthogonal assay.</div>
</section>
</div>
</body>
</html>"""


def main() -> None:
    ensure_dirs()
    tsv = RESULT / "therapy_spectrum_validation_matrix.tsv"
    summary = RESULT / "SUMMARY.md"
    html_path = HUB / "therapy_spectrum_validation_matrix_2026_05_11.html"

    write_tsv(tsv, ROWS)
    write_summary(summary, ROWS)
    html_path.write_text(build_html(ROWS), encoding="utf-8")

    shutil.copy2(tsv, ASSET_DIR / tsv.name)
    shutil.copy2(summary, ASSET_DIR / summary.name)
    shutil.copy2(tsv, LIVE_ASSET_DIR / tsv.name)
    shutil.copy2(summary, LIVE_ASSET_DIR / summary.name)
    shutil.copy2(html_path, LIVE_HUB / html_path.name)


if __name__ == "__main__":
    main()
