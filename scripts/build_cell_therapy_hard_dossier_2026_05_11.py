from __future__ import annotations

import csv
import html
import shutil
from pathlib import Path


ROOT = Path("/home/seungho/personal/THCA_data_analysis")
PROMISING = ROOT / "project/results/therapy_spectrum_promising_master_2026_05_11/therapy_spectrum_promising_master.tsv"
VALMAT = ROOT / "project/results/therapy_spectrum_validation_matrix_2026_05_11/therapy_spectrum_validation_matrix.tsv"
RESULT = ROOT / "project/results/cell_therapy_hard_dossier_2026_05_11"
HUB = ROOT / "project/papers_hub_2026_05_04"
LIVE_HUB = Path("/var/www/papers/papers_hub_2026_05_04")
ASSET_DIR = HUB / "assets/cell_therapy_hard_dossier_2026_05_11"
LIVE_ASSET_DIR = LIVE_HUB / "assets/cell_therapy_hard_dossier_2026_05_11"


FIGURE_BANK = {
    "CAR-T": [
        ("Fig 1", "Durable response versus early relapse", "Directly contrast durable responders and early relapsers using longitudinal CAR-T sampling, then show how native cytotoxic repertoire expansion tracks the outcome split."),
        ("Fig 2", "Construct-level signaling tradeoff", "Lay out how CD22-CAR linker length, tonic signaling, and antigen sensitivity trade off against each other; the point is to show that sequence-level design can rescue function without overclaiming universal effects."),
        ("Fig 3", "Manufacturing QC and dose efficiency", "Show the expansion process, phenotype preservation, and anti-lymphoma potency of the nanoparticle-expansion workflow; include the dose-sparing readout as a manufacturing argument rather than a clinical cure claim."),
        ("Fig 4", "Persistence and translational rewiring", "Place translational control, lineage state, and solid-tumor persistence in one panel so the reader sees that CAR-T engineering is not one axis but several separable levers."),
        ("Fig 5", "Exhaustion rescue and memory skew", "Use the BACH2 / lineage-state result to show how stem-like or memory-like programs can be reweighted away from dysfunctional exhaustion."),
        ("Fig 6", "Armored payload localization", "Show how local payload delivery improves tumor-restricted activity while the boundary box reminds the reader that this remains preclinical design support."),
    ],
    "CAR-NKT": [
        ("Fig 1", "Off-the-shelf CAR-NKT architecture", "Explain why HSPC-engineered, IL-15-enhanced CAR-NKT is a manufacturing and selection solution rather than just a new antigen receptor."),
        ("Fig 2", "Myeloid-malignancy synergy", "Show the CAR-NKT plus hypomethylating-agent interaction as the actual translation hook, not just the CAR construct."),
        ("Fig 3", "Responder phenotype from the clinical series", "Use the neuroblastoma phase-1 readout as the patient-side anchor that keeps the platform from floating away from real outcomes."),
        ("Fig 4", "Cross-disease validation corridor", "Compare AML/MDS and neuroblastoma as adjacent validation lanes, making the claim boundary explicit."),
    ],
    "TIL": [
        ("Fig 1", "RELB rescue of exhausted TILs", "Show how RELB shifts exhausted TILs toward a proliferative, memory/costimulatory-like state."),
        ("Fig 2", "TCR diversity preservation", "Use the TCR-seq result to show that expansion does not have to mean diversity collapse."),
        ("Fig 3", "Organoid killing and xenograft control", "Move from expression to function by pairing organoid polyfunctionality with tumor control in mouse xenografts."),
        ("Fig 4", "Potency and release logic", "Make the manufacturing implication explicit: which states should be in a TIL release assay and which should be treated as liability."),
    ],
    "Blinatumomab": [
        ("Fig 1", "T-cell dysfunction during continuous dosing", "Show that persistent blinatumomab exposure can erode T-cell function and therefore make duration a biological variable, not just a protocol detail."),
        ("Fig 2", "Re-dosing and TCF7-high memory state", "Show how re-dosing shifts the T-cell compartment toward central-memory-like cells with lower cytokine output."),
        ("Fig 3", "Cytokine-risk versus response window", "Contrast efficacy gain against cytokine control and make clear why these data are useful for schedule optimization."),
        ("Fig 4", "Clinical sampling map", "Lay out pre-, on-treatment, and re-dose timing so the reader can see where the biomarker has to be collected to be useful."),
    ],
}


TABLE_BANK = {
    "CAR-T": [
        ("Table 1", "Anchor datasets and tiers", "List every CAR-T dataset with score, tier, and the best task to keep the suite from being a pile of near-duplicate CAR-T papers."),
        ("Table 2", "Validation corridor", "Match each anchor to the strongest validation cohort or orthogonal support and state whether the evidence is strong, moderate, or still fragile."),
        ("Table 3", "Claim boundary table", "Explicitly mark what must not be said: mouse-only, process-only, construct-specific, or antigen-low model dependent."),
    ],
    "CAR-NKT": [
        ("Table 1", "Clinical and preclinical anchor map", "Put the AML/MDS CAR-NKT data and the neuroblastoma clinical series side by side so the reader sees the bridge and the gap."),
        ("Table 2", "Platform caveats", "State where the allogeneic/off-the-shelf story is directly supported and where it is still mostly a hypothesis."),
    ],
    "TIL": [
        ("Table 1", "Rescue versus potency metrics", "Pair proliferation, TCR diversity, organoid killing, and xenograft control to show the adoption path from biology to product."),
        ("Table 2", "Manufacturing risk table", "Flag where heterogeneity or exhaustion could break translation and what assay should catch it."),
    ],
    "Blinatumomab": [
        ("Table 1", "Biomarker and scheduling axes", "Organize the readouts by what they tell you about who should be re-dosed and when."),
        ("Table 2", "Safety boundary table", "Separate cytokine-risk biology from efficacy claims so the report stays claim-safe."),
    ],
}


MODALITY_SUMMARY = {
    "CAR-T": {
        "verdict": "strong",
        "anchor": "GSE290722",
        "why": "Best-supported axis in the bundle; direct response durability plus multiple validation corridors.",
        "boundary": "Keep the claim focused on response prediction, persistence, and engineering, not universal CAR-T biology.",
        "validation": ["GSE273170", "GSE243973", "GSE246342", "GSE145007", "GSE162975", "GSE223655", "GSE158676"],
    },
    "CAR-NKT": {
        "verdict": "moderate",
        "anchor": "GSE270430",
        "why": "Real translational hook with a clinical CAR-NKT series, but still a thin field.",
        "boundary": "Do not overgeneralize across diseases; keep it as an off-the-shelf cell-therapy platform pitch.",
        "validation": ["GSE223071", "GSE235923"],
    },
    "TIL": {
        "verdict": "moderate",
        "anchor": "GSE303442",
        "why": "Excellent rescue/potency signal and orthogonal functional support, but not a mature validation corridor.",
        "boundary": "Frame as product potency and manufacturing rescue rather than direct outcome proof.",
        "validation": ["Organoids", "Xenografts", "TCR-seq diversity retention", "Memory/costimulatory shift"],
    },
    "Blinatumomab": {
        "verdict": "moderate",
        "anchor": "GSE294241",
        "why": "Best for dosing-window / dysfunction / cytokine-risk optimization, not for a broad immunotherapy claim.",
        "boundary": "Keep it disease- and schedule-specific to B-ALL; it is a decision-support lane.",
        "validation": ["GSE196463", "PMCID 7923532", "B-ALL longitudinal response studies"],
    },
}


def esc(s: object) -> str:
    return html.escape(str(s), quote=True)


def ensure_dirs() -> None:
    RESULT.mkdir(parents=True, exist_ok=True)
    ASSET_DIR.mkdir(parents=True, exist_ok=True)
    LIVE_ASSET_DIR.mkdir(parents=True, exist_ok=True)
    HUB.mkdir(parents=True, exist_ok=True)
    LIVE_HUB.mkdir(parents=True, exist_ok=True)


def read_scores() -> dict[str, dict[str, str]]:
    rows: dict[str, dict[str, str]] = {}
    for src in [PROMISING, VALMAT]:
        with src.open("r", encoding="utf-8", newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                acc = row["accession"]
                prev = rows.get(acc, {})
                merged = dict(prev)
                for key, value in row.items():
                    if value not in ("", None):
                        merged[key] = value
                rows[acc] = merged
    return rows


def build_caption_cards(section: str, type_name: str, items: list[tuple[str, str, str]]) -> str:
    cards = "".join(
        f"<div class='caption-card'><div class='pill {type_name}'>{esc(label)}</div><h3>{esc(title)}</h3><p>{esc(desc)}</p></div>"
        for label, title, desc in items
    )
    return cards


def build_modal_section(modality: str, rows: dict[str, dict[str, str]]) -> str:
    meta = MODALITY_SUMMARY[modality]
    row = rows.get(meta["anchor"], {})
    figs = build_caption_cards(modality, "figure", FIGURE_BANK[modality])
    tabs = build_caption_cards(modality, "table", TABLE_BANK[modality])
    val = "".join(f"<li>{esc(v)}</li>" for v in meta["validation"])
    evidence_rows = [
        ("Anchor accession", meta["anchor"]),
        ("Anchor title", row.get("title", "n/a")),
        ("Score total", row.get("score_total", "n/a")),
        ("Tier", row.get("tier", "n/a")),
        ("Verdict", meta["verdict"]),
        ("Best task", row.get("best_task", "n/a")),
        ("Validation depth", f"{len(meta['validation'])} links"),
        ("Claim boundary", meta["boundary"]),
    ]
    evidence_html = "".join(
        f"<tr><th>{esc(k)}</th><td>{esc(v)}</td></tr>" for k, v in evidence_rows
    )
    return f"""
<section id="{modality.lower().replace('-','').replace('/','').replace(' ','_')}">
  <h2><span class="num">{esc(modality)}</span> {esc(modality)}</h2>
  <div class="sub">{esc(meta['anchor'])} | {esc(meta['verdict'])}</div>
  <div class="grid2">
    <div class="box good"><b>Why it matters.</b> {esc(meta['why'])}</div>
    <div class="box warn"><b>Claim boundary.</b> {esc(meta['boundary'])}</div>
  </div>
  <div class="subhead">Evidence matrix</div>
  <div class="box">
    <table class="evidence">
      <tbody>{evidence_html}</tbody>
    </table>
  </div>
  <div class="box"><b>Validation ladder.</b><ul>{val}</ul></div>
  <div class="subhead">Figure captions</div>
  <div class="grid">{figs}</div>
  <div class="subhead">Table captions</div>
  <div class="grid">{tabs}</div>
</section>
"""


def build_html(rows: dict[str, dict[str, str]]) -> str:
    cards = "".join(
        f"<a class='topic-card' href='#{mod.lower().replace('-','').replace('/','').replace(' ','_')}'><div class='topic-kicker'>{esc(mod)}</div><h3>{esc(meta['anchor'])}</h3><p class='topic-summary'>{esc(meta['why'])}</p><p class='topic-meta'><b>Verdict:</b> {esc(meta['verdict'])}</p><p class='topic-meta small'>{esc(meta['boundary'])}</p></a>"
        for mod, meta in MODALITY_SUMMARY.items()
    )
    sections = "".join(build_modal_section(mod, rows) for mod in ["CAR-T", "CAR-NKT", "TIL", "Blinatumomab"])
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8" />
<meta name="viewport" content="width=device-width,initial-scale=1.0" />
<title>Cell Therapy Hard Dossier · CAR-T / CAR-NKT / TIL / Blinatumomab</title>
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
.lead{{max-width:1120px;color:#cad5e5;font:18px/1.65 "Newsreader",serif}}
.stats{{display:grid;grid-template-columns:repeat(5,1fr);gap:10px;margin-top:24px}}
.stat{{background:rgba(255,255,255,.04);border:1px solid rgba(255,210,138,.18);padding:13px 14px;border-radius:10px}}
.stat b{{display:block;color:var(--gold);font:700 28px "Cormorant Garamond",serif;line-height:1}}
.stat span{{display:block;color:var(--muted);font:700 10px "JetBrains Mono",monospace;text-transform:uppercase;letter-spacing:.08em;margin-top:4px}}
.wrap{{max-width:1500px;margin:0 auto;padding:0 34px 80px}}
section{{padding:30px 0;border-bottom:1px solid var(--line)}}
h2{{font-size:36px;margin:0 0 6px}}
h2 .num{{font:700 13px "JetBrains Mono",monospace;color:var(--gold);letter-spacing:.18em;margin-right:10px}}
.sub{{font:700 10px "JetBrains Mono",monospace;color:var(--muted);letter-spacing:.12em;text-transform:uppercase;margin-bottom:16px}}
.subhead{{font:700 12px "JetBrains Mono",monospace;color:var(--gold);letter-spacing:.12em;margin:24px 0 12px;text-transform:uppercase}}
.grid{{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:14px}}
.grid2{{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:14px}}
.topic-card{{display:block;background:var(--panel);border:1px solid var(--line);border-radius:12px;padding:16px 18px;min-height:120px}}
.topic-summary{{color:#cad5e5;margin:10px 0 8px}}
.topic-meta{{color:var(--muted);margin:5px 0}}
.small{{font-size:12px}}
.caption-card{{background:var(--panel);border:1px solid var(--line);border-radius:12px;padding:14px 16px}}
.caption-card h3{{margin:10px 0 6px}}
.pill{{display:inline-block;padding:3px 9px;border-radius:999px;font:700 10px "JetBrains Mono",monospace;text-transform:uppercase;letter-spacing:.08em}}
.pill.figure{{background:rgba(53,211,157,.12);color:var(--teal);border:1px solid rgba(53,211,157,.22)}}
.pill.table{{background:rgba(126,182,255,.12);color:var(--blue);border:1px solid rgba(126,182,255,.22)}}
.box{{background:var(--panel);border:1px solid var(--line);border-left:4px solid var(--gold);border-radius:0 8px 8px 0;padding:14px 16px;margin:14px 0}}
.box.good{{border-left-color:var(--teal);background:#0a1a16}}
.box.warn{{border-left-color:var(--red);background:#1d1310}}
.box ul{{margin:8px 0 0 18px}}
.box li{{margin:4px 0}}
.evidence{{width:100%;border-collapse:collapse}}
.evidence th,.evidence td{{text-align:left;vertical-align:top;padding:7px 10px;border-bottom:1px solid rgba(255,255,255,.08)}}
.evidence th{{width:210px;color:var(--gold);font:700 11px "JetBrains Mono",monospace;text-transform:uppercase;letter-spacing:.08em}}
.evidence td{{color:#dde7f5}}
@media(max-width:1100px){{.grid,.grid2{{grid-template-columns:1fr 1fr}}}}
@media(max-width:980px){{.stats{{grid-template-columns:repeat(2,1fr)}}.grid,.grid2{{grid-template-columns:1fr}}h1{{font-size:40px}}.wrap{{padding:0 20px 60px}}}}
</style>
</head>
<body>
<header class="hero">
  <div class="hero-inner">
    <div class="kicker">Hard dossier · 2026-05-11</div>
    <h1>Cell Therapy Hard Dossier<br/><em>CAR-T / CAR-NKT / TIL / Blinatumomab</em></h1>
    <p class="lead">
      This is the harder version: each modality gets a claim boundary, validation ladder, figure-caption bank, and table-caption bank.
      It is built for sending, but it is also explicit about what can and cannot be claimed from the data.
    </p>
    <div style="margin-top:18px;color:#9aabc0;font:12px 'JetBrains Mono',monospace;">
      <a href="index.html" style="color:#35d39d">Hub index</a> ·
      <a href="cell_therapy_caption_atlas_2026_05_11.html" style="color:#c79cff">caption atlas</a> ·
      <a href="car_t_bundle_2026_05_11.html" style="color:#35d39d">CAR-T bundle</a>
    </div>
    <div class="stats">
      <div class="stat"><b>4</b><span>Modalities</span></div>
      <div class="stat"><b>18</b><span>Figure captions</span></div>
      <div class="stat"><b>9</b><span>Table captions</span></div>
      <div class="stat"><b>3</b><span>Strong / 1 moderate</span></div>
      <div class="stat"><b>Sendable</b><span>Hard pack</span></div>
    </div>
  </div>
</header>
<div class="wrap">
<section>
  <h2><span class="num">01</span>Index</h2>
  <div class="sub">What is in the hard pack</div>
  <div class="grid">{cards}</div>
</section>
<section>
  <h2><span class="num">02</span>Decision matrix</h2>
  <div class="sub">What survives a hard read</div>
  <div class="box">
    <table class="evidence">
      <tbody>
        <tr><th>Pattern</th><td>CAR-T is the only strong lane; the other three are useful but stay below universal-claim level.</td></tr>
        <tr><th>Send now</th><td>GSE290722 for response durability; GSE303442 for TIL potency; GSE270430 for off-the-shelf CAR-NKT; GSE294241 for blinatumomab schedule/dysfunction.</td></tr>
        <tr><th>Do not claim</th><td>Universal CAR-T biology, disease-agnostic CAR-NKT effect, direct TIL clinical efficacy, or broad blinatumomab immunotherapy generalization.</td></tr>
        <tr><th>Validation floor</th><td>CAR-T has a multi-cohort corridor; mRNA-LNP and CDK4/6 lanes remain stronger outside this cell-therapy pack but are not the sendable focus here.</td></tr>
      </tbody>
    </table>
  </div>
</section>
{sections}
</div>
</body>
</html>"""


def main() -> None:
    ensure_dirs()
    rows = read_scores()
    tsv = RESULT / "cell_therapy_hard_dossier.tsv"
    summary = RESULT / "SUMMARY.md"
    html_path = HUB / "cell_therapy_hard_dossier_2026_05_11.html"

    # simple TSV combining the selected accessions and strong/moderate verdicts
    with tsv.open("w", encoding="utf-8", newline="") as f:
        fieldnames = ["accession", "modality", "tier", "score_total", "verdict", "claim_boundary", "validation_links"]
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        for modality, meta in MODALITY_SUMMARY.items():
            anchor = meta["anchor"]
            row = rows.get(anchor, {})
            w.writerow(
                {
                    "accession": anchor,
                    "modality": modality,
                    "tier": row.get("tier", ""),
                    "score_total": row.get("score_total", ""),
                    "verdict": meta["verdict"],
                    "claim_boundary": meta["boundary"],
                    "validation_links": ", ".join(meta["validation"]),
                }
            )

    summary.write_text(
        "\n".join(
                [
                "# Cell Therapy Hard Dossier",
                "",
                "Date: 2026-05-11",
                "",
                "Harder report for CAR-T, CAR-NKT, TIL, and blinatumomab with explicit claim boundaries and validation ladders.",
                "",
                "Snapshot:",
                f"- Modalities: {len(MODALITY_SUMMARY)}",
                f"- Figure captions: {sum(len(v) for v in FIGURE_BANK.values())}",
                f"- Table captions: {sum(len(v) for v in TABLE_BANK.values())}",
                f"- Strong: {sum(1 for m in MODALITY_SUMMARY.values() if m['verdict']=='strong')}",
                f"- Moderate: {sum(1 for m in MODALITY_SUMMARY.values() if m['verdict']=='moderate')}",
                "",
                "Use:",
                "- Send this when you need a rigorous, not fluffy, cell-therapy summary.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    html_path.write_text(build_html(rows), encoding="utf-8")

    for src in [tsv, summary]:
        shutil.copy2(src, ASSET_DIR / src.name)
        shutil.copy2(src, LIVE_ASSET_DIR / src.name)
    shutil.copy2(html_path, LIVE_HUB / html_path.name)
    shutil.copy2(html_path, ASSET_DIR / html_path.name)
    shutil.copy2(html_path, LIVE_ASSET_DIR / html_path.name)


if __name__ == "__main__":
    main()
