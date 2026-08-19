from __future__ import annotations

import csv
import html
import shutil
from pathlib import Path


ROOT = Path("/home/seungho/personal/THCA_data_analysis")
RESULT = ROOT / "project/results/therapy_spectrum_triage_2026_05_11"
HUB = ROOT / "project/papers_hub_2026_05_04"
LIVE_HUB = Path("/var/www/papers/papers_hub_2026_05_04")
ASSET_DIR = HUB / "assets/therapy_spectrum_triage_2026_05_11"
LIVE_ASSET_DIR = LIVE_HUB / "assets/therapy_spectrum_triage_2026_05_11"


MODALITIES = [
    {
        "modality": "neoantigen / mRNA vaccine / ICI combo",
        "group": "immune-vaccine",
        "paperability": 5,
        "validation_readiness": 5,
        "bd_readiness": 5,
        "default_status": "GO",
        "best_in_silico_task": "response ranking, T cell state, clonality, combo prioritization",
        "example_geo_accessions": "GSE222011; GSE249674; GSE303153; GSE260817",
        "search_keywords": "neoantigen, mRNA vaccine, ICI, combo, TCR, scRNA",
        "why_it_works": "response-linked paired samples, TCR/scRNA hooks, and a direct companion biomarker path",
        "caveat": "avoid claiming generalized vaccine efficacy from a single cohort",
        "agent_query": "site:ncbi.nlm.nih.gov/geo neoantigen vaccine ICI mRNA TCR scRNA",
    },
    {
        "modality": "CAR-T / CAR-NK / CAR-NKT / TIL / TCR-T",
        "group": "cell-therapy",
        "paperability": 5,
        "validation_readiness": 4,
        "bd_readiness": 5,
        "default_status": "GO",
        "best_in_silico_task": "potency, exhaustion, manufacturing QC, response stratification",
        "example_geo_accessions": "GSE266618; GSE278506; GSE223071; GSE217387; GSE185205",
        "search_keywords": "CAR-T, CAR-NK, CAR-NKT, TIL, TCR-T, exhaustion, potency",
        "why_it_works": "strong transcriptomic signal, clear ex vivo readouts, and direct translational relevance",
        "caveat": "BD should be potency/QC/selection, not de novo cell-product invention",
        "agent_query": "site:ncbi.nlm.nih.gov/geo CAR-T CAR-NK CAR-NKT TIL TCR-T exhaustion potency",
    },
    {
        "modality": "T cell engager / bispecific / checkpoint / agonist antibody",
        "group": "antibody-therapy",
        "paperability": 4,
        "validation_readiness": 3,
        "bd_readiness": 4,
        "default_status": "WATCH",
        "best_in_silico_task": "target selection, TME rewiring, resistance mapping",
        "example_geo_accessions": "GSE189359; GSE260674; GSE292621; GSE260817",
        "search_keywords": "bispecific, checkpoint, CD3 engager, agonist antibody, TME",
        "why_it_works": "actionable mechanism and combo logic, often paired with immune signatures",
        "caveat": "validation often depends on model system, not direct patient response",
        "agent_query": "site:ncbi.nlm.nih.gov/geo bispecific CD3 engager checkpoint agonist antibody",
    },
    {
        "modality": "ADC / radioligand / cell-surface targeting",
        "group": "targeted-delivery",
        "paperability": 4,
        "validation_readiness": 3,
        "bd_readiness": 4,
        "default_status": "WATCH",
        "best_in_silico_task": "target expression, lineage plasticity, payload sensitivity",
        "example_geo_accessions": "GSE274871; GSE289746; GSE236808; GSE264573",
        "search_keywords": "ADC, radioligand, cell-surface target, lineage plasticity",
        "why_it_works": "strong for target prioritization and resistance logic",
        "caveat": "BD strongest when framed as companion diagnostics or target selection",
        "agent_query": "site:ncbi.nlm.nih.gov/geo antibody-drug conjugate radioligand lineage plasticity",
    },
    {
        "modality": "PARP / DNA damage response / radiosensitizer",
        "group": "small-molecule-ddr",
        "paperability": 5,
        "validation_readiness": 5,
        "bd_readiness": 4,
        "default_status": "GO",
        "best_in_silico_task": "resistance, sensitization, combo design, biomarker discovery",
        "example_geo_accessions": "GSE289437; GSE315883; GSE233820; GSE285648",
        "search_keywords": "PARP inhibitor, DDR, radiosensitizer, resistance, sensitivity",
        "why_it_works": "large public footprint, repeatable perturbation logic, and easy cross-cohort tests",
        "caveat": "BD is narrower unless tied to a clear biomarker or drug-sensitivity panel",
        "agent_query": "site:ncbi.nlm.nih.gov/geo PARP inhibitor DNA damage response radiosensitizer resistance",
    },
    {
        "modality": "kinase / MAPK / RTK inhibition",
        "group": "small-molecule-kinase",
        "paperability": 5,
        "validation_readiness": 5,
        "bd_readiness": 4,
        "default_status": "GO",
        "best_in_silico_task": "pathway response, resistance, combo prioritization",
        "example_geo_accessions": "GSE289437; GSE315883; GSE260551",
        "search_keywords": "MEK, ERK, RAF, RTK, kinase inhibitor, MAPK",
        "why_it_works": "clean perturbational logic and high reuse across cancer types",
        "caveat": "separate pathway response from direct therapeutic efficacy claims",
        "agent_query": "site:ncbi.nlm.nih.gov/geo MEK ERK RAF kinase inhibitor MAPK",
    },
    {
        "modality": "hormone / endocrine / AR / ER",
        "group": "endocrine-therapy",
        "paperability": 4,
        "validation_readiness": 4,
        "bd_readiness": 3,
        "default_status": "GO",
        "best_in_silico_task": "lineage-state shifts, resistance, biomarker discovery",
        "example_geo_accessions": "GSE80077; GSE289437",
        "search_keywords": "hormone therapy, AR, ER, endocrine resistance, lineage state",
        "why_it_works": "well-defined treatment states with straightforward responder/non-responder logic",
        "caveat": "BD needs a sharply defined clinical indication, not a pan-cancer claim",
        "agent_query": "site:ncbi.nlm.nih.gov/geo hormone therapy AR ER endocrine resistance",
    },
    {
        "modality": "epigenetic priming / chromatin therapy",
        "group": "epigenetic",
        "paperability": 4,
        "validation_readiness": 4,
        "bd_readiness": 3,
        "default_status": "GO",
        "best_in_silico_task": "priming, reactivation, immune-rewiring signatures",
        "example_geo_accessions": "GSE168225; GSE179674; GSE188249",
        "search_keywords": "DNMT inhibitor, HDAC inhibitor, epigenetic priming, chromatin",
        "why_it_works": "works well as a response-enabling layer with transcriptomic readouts",
        "caveat": "avoid overcalling direct causality from methylation/transcript changes alone",
        "agent_query": "site:ncbi.nlm.nih.gov/geo epigenetic priming DNMT HDAC chromatin",
    },
    {
        "modality": "oncolytic virus",
        "group": "viral-therapy",
        "paperability": 3,
        "validation_readiness": 3,
        "bd_readiness": 3,
        "default_status": "RESERVE",
        "best_in_silico_task": "permissivity, immune remodeling, viral susceptibility",
        "example_geo_accessions": "GSE244034; GSE260516; GSE223395; GSE222002",
        "search_keywords": "oncolytic virus, viral permissivity, immune remodeling",
        "why_it_works": "fascinating biology, but cohort sizes and model diversity vary widely",
        "caveat": "best used as a special-case module, not the platform center",
        "agent_query": "site:ncbi.nlm.nih.gov/geo oncolytic virus permissivity immune remodeling",
    },
    {
        "modality": "cytokine engineering / costimulation / immune agonist",
        "group": "immune-engineering",
        "paperability": 4,
        "validation_readiness": 3,
        "bd_readiness": 4,
        "default_status": "WATCH",
        "best_in_silico_task": "potency augmentation, T cell expansion, exhaustion rescue",
        "example_geo_accessions": "GSE249674; GSE260817; GSE189359",
        "search_keywords": "IL-15, IL-21, 4-1BB, cytokine engineering, costimulation",
        "why_it_works": "good for mechanism and combo ranking, especially in immune contexts",
        "caveat": "validation is often preclinical and needs orthogonal support",
        "agent_query": "site:ncbi.nlm.nih.gov/geo cytokine engineering 4-1BB IL-15 IL-21 costimulation",
    },
    {
        "modality": "innate immune agonist / myeloid reprogramming",
        "group": "innate-immunity",
        "paperability": 4,
        "validation_readiness": 3,
        "bd_readiness": 3,
        "default_status": "WATCH",
        "best_in_silico_task": "TME reprogramming, myeloid suppression, combo selection",
        "example_geo_accessions": "GSE281467; GSE184009; GSE222002",
        "search_keywords": "myeloid, innate immune agonist, CD40, TLR, STING, macrophage",
        "why_it_works": "routes directly into TME remodeling and IO-combo logic",
        "caveat": "interpretation often depends on mouse model choice and dosing",
        "agent_query": "site:ncbi.nlm.nih.gov/geo CD40 TLR STING myeloid reprogramming",
    },
    {
        "modality": "radiotherapy / radiosensitizer / IO combo",
        "group": "combo-radiation",
        "paperability": 4,
        "validation_readiness": 4,
        "bd_readiness": 3,
        "default_status": "GO",
        "best_in_silico_task": "synergy ranking, immune activation, chemokine induction",
        "example_geo_accessions": "GSE233820; GSE260674; GSE189359",
        "search_keywords": "radiotherapy, radiosensitizer, anti-PD-L1, chemokine, synergy",
        "why_it_works": "combination logic is explicit and the transcriptional response is assayable",
        "caveat": "avoid extrapolating from model-level synergy to universal clinical benefit",
        "agent_query": "site:ncbi.nlm.nih.gov/geo radiotherapy radiosensitizer immune combo chemokine",
    },
    {
        "modality": "cell therapy manufacturing / bioreactor / expansion",
        "group": "platform-qc",
        "paperability": 4,
        "validation_readiness": 5,
        "bd_readiness": 5,
        "default_status": "GO",
        "best_in_silico_task": "manufacturing QC, expansion efficiency, process optimization",
        "example_geo_accessions": "GSE261103; GSE295031; GSE266618",
        "search_keywords": "manufacturing, bioreactor, expansion, QC, process optimization",
        "why_it_works": "clear process metric and direct industrial relevance",
        "caveat": "BD should focus on process analytics, not therapy claim escalation",
        "agent_query": "site:ncbi.nlm.nih.gov/geo CAR T manufacturing bioreactor expansion QC",
    },
    {
        "modality": "adoptive cell state / exhaustion / TCR dynamics",
        "group": "state-modeling",
        "paperability": 4,
        "validation_readiness": 4,
        "bd_readiness": 4,
        "default_status": "GO",
        "best_in_silico_task": "state scoring, clone expansion, post-treatment trajectory",
        "example_geo_accessions": "GSE123813; GSE156728; GSE185205; GSE221776",
        "search_keywords": "T cell state, TCR dynamics, exhaustion, clonotype, trajectory",
        "why_it_works": "excellent for trajectory modeling and response-state interpretation",
        "caveat": "clone expansion is not itself proof of better treatment efficacy",
        "agent_query": "site:ncbi.nlm.nih.gov/geo T cell state TCR dynamics exhaustion clonotype",
    },
    {
        "modality": "stromal / TME reprogramming",
        "group": "tme-reprogramming",
        "paperability": 4,
        "validation_readiness": 4,
        "bd_readiness": 3,
        "default_status": "WATCH",
        "best_in_silico_task": "fibroblast / endothelial / stromal response mapping",
        "example_geo_accessions": "GSE196065; GSE168225; GSE260674",
        "search_keywords": "stroma, fibroblast, endothelial, TME, reprogramming",
        "why_it_works": "turns resistance and microenvironment state into a computable axis",
        "caveat": "needs careful separation of stroma-specific and tumor-intrinsic effects",
        "agent_query": "site:ncbi.nlm.nih.gov/geo stroma fibroblast endothelial TME reprogramming",
    },
    {
        "modality": "resistance / combo prioritization / response modeling",
        "group": "meta-platform",
        "paperability": 5,
        "validation_readiness": 5,
        "bd_readiness": 5,
        "default_status": "GO",
        "best_in_silico_task": "therapy triage, responder ranking, decision support across modalities",
        "example_geo_accessions": "GSE130157; GSE249630; GSE217387; GSE185205",
        "search_keywords": "response modeling, resistance, therapy prioritization, combo, triage",
        "why_it_works": "highest generality and best fit for an agentic ranking framework",
        "caveat": "best positioned as decision support, not therapy invention",
        "agent_query": "site:ncbi.nlm.nih.gov/geo response modeling resistance therapy prioritization combo",
    },
]


def compute_total(row: dict[str, object]) -> int:
    return int(row["paperability"]) + int(row["validation_readiness"]) + int(row["bd_readiness"])


def ensure_dirs() -> None:
    RESULT.mkdir(parents=True, exist_ok=True)
    HUB.mkdir(parents=True, exist_ok=True)
    ASSET_DIR.mkdir(parents=True, exist_ok=True)
    LIVE_HUB.mkdir(parents=True, exist_ok=True)
    LIVE_ASSET_DIR.mkdir(parents=True, exist_ok=True)


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    headers = [
        "modality",
        "group",
        "paperability",
        "validation_readiness",
        "bd_readiness",
        "score_total",
        "default_status",
        "best_in_silico_task",
        "example_geo_accessions",
        "search_keywords",
        "why_it_works",
        "caveat",
        "agent_query",
    ]
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=headers, delimiter="\t")
        writer.writeheader()
        for row in rows:
            out = dict(row)
            out["score_total"] = compute_total(out)
            writer.writerow(out)


def write_query_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    headers = ["modality", "default_status", "agent_query", "search_keywords", "example_geo_accessions"]
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=headers, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row[k] for k in headers})


def write_summary(path: Path, rows: list[dict[str, object]]) -> None:
    go = sum(1 for r in rows if r["default_status"] == "GO")
    watch = sum(1 for r in rows if r["default_status"] == "WATCH")
    reserve = sum(1 for r in rows if r["default_status"] == "RESERVE")
    best = sorted(rows, key=lambda r: (-compute_total(r), -int(r["paperability"]), r["modality"]))
    lines = [
        "# Therapy Spectrum Triage",
        "",
        "Date: 2026-05-11",
        "",
        "This folder extends the GEO therapy-spectrum triage into a larger modality map and adds query templates for automation.",
        "",
        "Files:",
        "- `therapy_spectrum_triage.tsv` - expanded modality triage table",
        "- `therapy_spectrum_query_templates.tsv` - GEO query templates for the agent screener",
        "",
        "Snapshot:",
        f"- Rows: {len(rows)}",
        f"- GO / WATCH / RESERVE: {go} / {watch} / {reserve}",
        f"- Top total-score lanes: {best[0]['modality']}; {best[1]['modality']}; {best[2]['modality']}",
        "",
        "Operating rule:",
        "- Frame the platform as therapy triage + response prediction + biomarker selection + potency/QC + combo prioritization.",
        "- Promote only studies with intervention + response readout + deployable output.",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def esc(s: object) -> str:
    return html.escape(str(s), quote=True)


def build_table(rows: list[dict[str, object]]) -> str:
    trs = []
    for r in rows:
        status_cls = {"GO": "go", "WATCH": "watch", "RESERVE": "reserve"}[r["default_status"]]
        trs.append(
            "<tr>"
            f"<td>{esc(r['modality'])}</td>"
            f"<td>{esc(r['group'])}</td>"
            f"<td>{r['paperability']}</td>"
            f"<td>{r['validation_readiness']}</td>"
            f"<td>{r['bd_readiness']}</td>"
            f"<td>{compute_total(r)}</td>"
            f"<td><span class='tag {status_cls}'>{esc(r['default_status'])}</span></td>"
            f"<td>{esc(r['best_in_silico_task'])}</td>"
            f"<td>{esc(r['example_geo_accessions'])}</td>"
            f"<td>{esc(r['search_keywords'])}</td>"
            f"<td>{esc(r['why_it_works'])}</td>"
            f"<td>{esc(r['caveat'])}</td>"
            f"<td><code>{esc(r['agent_query'])}</code></td>"
            "</tr>"
        )
    return "\n".join(trs)


def build_html(rows: list[dict[str, object]]) -> str:
    sorted_rows = sorted(rows, key=lambda r: (-compute_total(r), -int(r["paperability"]), r["modality"]))
    top = sorted_rows[:6]
    html_rows = build_table(sorted_rows)
    top_cards = "".join(
        f"<div class='card'><h3>{i+1}. {esc(r['modality'])}</h3><p><b>Score:</b> {compute_total(r)} | <b>Status:</b> {esc(r['default_status'])}</p><p>{esc(r['best_in_silico_task'])}</p><p class='small'>{esc(r['why_it_works'])}</p></div>"
        for i, r in enumerate(top[:6])
    )
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8" />
<meta name="viewport" content="width=device-width,initial-scale=1.0" />
<title>Therapy Spectrum Triage · Expanded GEO Agent Map</title>
<style>
:root{{--bg:#061019;--panel:#0d1726;--line:#23324a;--ink:#eaf1fb;--muted:#94a5ba;--gold:#ffd28a;--teal:#35d39d;--red:#ff8a6b;--blue:#7eb6ff}}
*{{box-sizing:border-box}}
html,body{{margin:0;background:var(--bg);color:var(--ink);font:14px/1.55 Inter,"Noto Sans KR",-apple-system,sans-serif}}
a{{color:var(--teal);text-decoration:none}} a:hover{{color:var(--gold);text-decoration:underline}}
code,pre,.mono{{font-family:"JetBrains Mono","SF Mono",Menlo,monospace;font-size:12px}}
h1,h2,h3{{font-family:"Cormorant Garamond","Newsreader",serif;letter-spacing:-.02em;color:#fff7dd}}
.hero{{padding:64px 34px 44px;background:radial-gradient(circle at top left,#162540 0%,#0a1220 60%,#050810 100%);border-bottom:1px solid var(--line);position:relative;overflow:hidden}}
.hero:before{{content:"";position:absolute;inset:0;background:url("data:image/svg+xml;utf8,<svg xmlns='http://www.w3.org/2000/svg' width='220' height='220'><g fill='none' stroke='%23ffd28a' stroke-width='0.45' opacity='0.12'><circle cx='55' cy='55' r='38'/><circle cx='165' cy='50' r='24'/><circle cx='68' cy='168' r='30'/><circle cx='164' cy='162' r='48'/><path d='M55 55 L165 162 M165 50 L68 168'/></g></svg>");pointer-events:none}}
.hero-inner{{position:relative;max-width:1380px;margin:0 auto;z-index:2}}
.kicker{{font:700 11px "JetBrains Mono",monospace;color:var(--gold);letter-spacing:.24em;text-transform:uppercase}}
h1{{font-size:64px;line-height:.98;margin:12px 0 12px}}
.lead{{max-width:1040px;color:#cad5e5;font:18px/1.65 "Newsreader",serif}}
.lead b{{color:#fff7dd}}
.stats{{display:grid;grid-template-columns:repeat(6,1fr);gap:10px;margin-top:24px}}
.stat{{background:rgba(255,255,255,.04);border:1px solid rgba(255,210,138,.18);padding:13px 14px;border-radius:10px}}
.stat b{{display:block;color:var(--gold);font:700 28px "Cormorant Garamond",serif;line-height:1}}
.stat span{{display:block;color:var(--muted);font:700 10px "JetBrains Mono",monospace;text-transform:uppercase;letter-spacing:.08em;margin-top:4px}}
.crumbs{{margin-top:18px;color:var(--muted);font:12px "JetBrains Mono",monospace}}
.wrap{{max-width:1380px;margin:0 auto;display:grid;grid-template-columns:260px 1fr;gap:34px;padding:0 34px 80px}}
.toc{{position:sticky;top:0;align-self:start;max-height:100vh;overflow:auto;border-right:1px solid var(--line);padding:28px 18px 28px 0}}
.toc h4{{font:700 10px "JetBrains Mono",monospace;color:var(--gold);letter-spacing:.18em;text-transform:uppercase;margin:0 0 12px}}
.toc a{{display:block;color:#cbd5e5;padding:5px 0;font-size:12px}}
main{{min-width:0;padding-top:24px}}
section{{padding:32px 0;border-bottom:1px solid var(--line)}}
h2{{font-size:36px;margin:0 0 6px}}
h2 .num{{font:700 13px "JetBrains Mono",monospace;color:var(--gold);letter-spacing:.18em;margin-right:10px}}
.sub{{font:700 10px "JetBrains Mono",monospace;color:var(--muted);letter-spacing:.12em;text-transform:uppercase;margin-bottom:16px}}
.t{{width:100%;border-collapse:collapse;margin:12px 0 16px;font-size:12px}}
.t th,.t td{{border:1px solid var(--line);padding:8px 9px;vertical-align:top}}
.t th{{background:#16243a;color:var(--gold);font:700 10px "JetBrains Mono",monospace;letter-spacing:.05em;text-transform:uppercase;text-align:left}}
.t tr:nth-child(even) td{{background:rgba(255,255,255,.02)}}
.tag{{display:inline-block;padding:2px 8px;border-radius:99px;font:700 10px "JetBrains Mono",monospace;letter-spacing:.06em;text-transform:uppercase;margin-right:4px;border:1px solid transparent}}
.tag.go{{background:rgba(53,211,157,.12);border-color:rgba(53,211,157,.28);color:var(--teal)}}
.tag.watch{{background:rgba(255,210,138,.12);border-color:rgba(255,210,138,.28);color:var(--gold)}}
.tag.reserve{{background:rgba(255,138,107,.12);border-color:rgba(255,138,107,.28);color:var(--red)}}
.box{{background:var(--panel);border:1px solid var(--line);border-left:4px solid var(--gold);border-radius:0 8px 8px 0;padding:14px 16px;margin:14px 0}}
.box.good{{border-left-color:var(--teal);background:#0a1a16}}
.box.warn{{border-left-color:var(--red);background:#1d1310}}
.grid{{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:14px}}
.card{{background:var(--panel);border:1px solid var(--line);border-radius:10px;padding:14px 16px}}
.card h3{{margin:0 0 8px}}
pre{{white-space:pre-wrap;background:#08101a;border:1px solid var(--line);padding:14px;border-radius:8px;color:#d8e2f0;overflow:auto}}
.small{{color:var(--muted);font-size:12px}}
@media(max-width:1100px){{.grid{{grid-template-columns:1fr 1fr}}}}
@media(max-width:980px){{.wrap{{display:block;padding:0 20px 60px}}.toc{{display:none}}.stats{{grid-template-columns:repeat(2,1fr)}}h1{{font-size:42px}}.grid{{grid-template-columns:1fr}}}}
</style>
</head>
<body>
<header class="hero">
  <div class="hero-inner">
    <div class="kicker">GEO / therapy-spectrum triage · expanded · 2026-05-11</div>
    <h1>Therapy Spectrum<br/><em>Expanded for Agentic GEO Scanning</em></h1>
    <p class="lead">
      This is the expanded version of the triage pack. The platform should be driven by
      <b>therapy triage, response prediction, biomarker selection, potency/QC, and combo prioritization</b>.
      The strongest lanes are the neoantigen / mRNA vaccine / ICI family, engineered cell therapy, PARP/DDR and kinase inhibition, and the cross-modal response / resistance meta-platform.
    </p>
    <div class="stats">
      <div class="stat"><b>{len(rows)}</b><span>Modalities triaged</span></div>
      <div class="stat"><b>{sum(1 for r in rows if r["default_status"] == "GO")}</b><span>Go lanes</span></div>
      <div class="stat"><b>{sum(1 for r in rows if r["default_status"] == "WATCH")}</b><span>Watch lanes</span></div>
      <div class="stat"><b>{sum(1 for r in rows if r["default_status"] == "RESERVE")}</b><span>Reserve lanes</span></div>
      <div class="stat"><b>GEO</b><span>Daily-updated source</span></div>
      <div class="stat"><b>agent</b><span>Query-driven screening</span></div>
    </div>
    <div class="crumbs"><a href="index.html">← Hub index</a> · <a href="therapy_spectrum_recent_registry_2026_05_11.html">recent GEO registry</a> · <a href="therapy_spectrum_go_shortlist_2026_05_11.html">GO shortlist</a> · <a href="therapy_spectrum_ready_split_2026_05_11.html">paper/validation/BD split</a> · <a href="assets/therapy_spectrum_triage_2026_05_11/therapy_spectrum_triage.tsv">raw TSV</a> · <a href="assets/therapy_spectrum_triage_2026_05_11/therapy_spectrum_query_templates.tsv">query TSV</a> · <a href="assets/therapy_spectrum_triage_2026_05_11/SUMMARY.md">summary</a></div>
  </div>
</header>
<div class="wrap">
<nav class="toc">
  <h4>Contents</h4>
  <a href="#tldr">01 TL;DR</a>
  <a href="#rules">02 Rules</a>
  <a href="#top">03 Top lanes</a>
  <a href="#matrix">04 Modality matrix</a>
  <a href="#pipeline">05 Agent pipeline</a>
  <a href="#queries">06 Query templates</a>
  <a href="#decision">07 Decision summary</a>
</nav>
<main>
<section id="tldr">
  <h2><span class="num">01</span>TL;DR</h2>
  <div class="sub">The right question is not "what is the coolest therapy"; it is "what becomes a reproducible decision system"</div>
  <table class="t">
    <thead><tr><th>Lane</th><th>Total score</th><th>Paperability</th><th>Validation</th><th>BD</th><th>Disposition</th></tr></thead>
    <tbody>
      <tr><td>neoantigen / mRNA vaccine / ICI combo</td><td>15</td><td>5</td><td>5</td><td>5</td><td><span class="tag go">Go</span></td></tr>
      <tr><td>CAR-T / CAR-NK / CAR-NKT / TIL / TCR-T</td><td>14</td><td>5</td><td>4</td><td>5</td><td><span class="tag go">Go</span></td></tr>
      <tr><td>PARP / DNA damage response / radiosensitizer</td><td>14</td><td>5</td><td>5</td><td>4</td><td><span class="tag go">Go</span></td></tr>
      <tr><td>kinase / MAPK / RTK inhibition</td><td>14</td><td>5</td><td>5</td><td>4</td><td><span class="tag go">Go</span></td></tr>
      <tr><td>resistance / combo prioritization / response modeling</td><td>15</td><td>5</td><td>5</td><td>5</td><td><span class="tag go">Go</span></td></tr>
    </tbody>
  </table>
  <div class="box good"><b>Best business frame:</b> sell a decision engine that ranks therapies, selects companion biomarkers, and prioritizes potency/QC or combo follow-up. That is easier to validate and productize than generic GEO mining.</div>
</section>

<section id="rules">
  <h2><span class="num">02</span>Rules</h2>
  <div class="sub">Promotion requires intervention plus response readout plus a deployable output</div>
  <div class="grid">
    <div class="card"><h3>Paperable</h3><p>Coherent mechanism, replication, and a figure/table story.</p><p class="small">Paired pre/post, responder vs non-responder, or orthogonal validation strengthen this lane.</p></div>
    <div class="card"><h3>Validation-ready</h3><p>Cheap orthogonal assay or external cohort is feasible.</p><p class="small">This is the sweet spot for biomarkers, potency/QC, and combo ranking.</p></div>
    <div class="card"><h3>BD-ready</h3><p>Output changes a decision or becomes a product artifact.</p><p class="small">Companion diagnostic, batch release, candidate ranking, or responder selection.</p></div>
  </div>
</section>

<section id="top">
  <h2><span class="num">03</span>Top Lanes</h2>
  <div class="sub">The six lanes that should be searched first every time</div>
  <div class="grid">
    {top_cards}
  </div>
</section>

<section id="matrix">
  <h2><span class="num">04</span>Modality Matrix</h2>
  <div class="sub">Expanded lane map with query hooks</div>
  <table class="t">
    <thead>
      <tr><th>Modality</th><th>Group</th><th>Total</th><th>Status</th><th>Best task</th><th>Example GEO</th><th>Query keywords</th><th>Why it works</th><th>Caveat</th><th>Agent query</th></tr>
    </thead>
    <tbody>
      {html_rows}
    </tbody>
  </table>
</section>

<section id="pipeline">
  <h2><span class="num">05</span>Agent Pipeline</h2>
  <div class="sub">Search, filter, score, route, gate</div>
  <pre>1. Search GEO by modality keyword families
   - vaccine / neoantigen / mRNA / ICI / TCR
   - CAR-T / CAR-NK / TIL / TCR-T / exhaustion / potency
   - bispecific / checkpoint / ADC / radioligand
   - PARP / MEK / ERK / hormone / epigenetic
   - oncolytic / radiotherapy / cytokine / myeloid / stroma

2. Filter by study geometry
   - keep intervention + control
   - keep pre/post or responder/non-responder
   - prefer paired human samples, then humanized models, then informative preclinical systems

3. Score evidence
   - paired design
   - sample size
   - orthogonal validation
   - independent replication
   - clinical actionability

4. Route to task
   - response ranking
   - potency / QC
   - resistance mapping
   - target prioritization
   - combo prioritization

5. Gate output
   - paperable only if mechanism + replication + clear figure exist
   - validation-ready only if a cheap orthogonal test exists
   - BD-ready only if the output guides a decision or workflow</pre>
  <div class="box warn"><b>Important.</b> Atlas-only or descriptive studies stay background unless they convert into a decision artifact. This keeps the agent from generating pretty but non-actionable GEO summaries.</div>
</section>

<section id="queries">
  <h2><span class="num">06</span>Query Templates</h2>
  <div class="sub">Directly reusable search strings</div>
  <table class="t">
    <thead><tr><th>Modality</th><th>Agent query</th><th>Use case</th></tr></thead>
    <tbody>
      {''.join(f"<tr><td>{esc(r['modality'])}</td><td><code>{esc(r['agent_query'])}</code></td><td>{esc(r['best_in_silico_task'])}</td></tr>" for r in sorted_rows[:12])}
    </tbody>
  </table>
</section>

<section id="decision">
  <h2><span class="num">07</span>Decision Summary</h2>
  <div class="sub">Build order</div>
  <table class="t">
    <thead><tr><th>Priority</th><th>Build</th><th>Reason</th></tr></thead>
    <tbody>
      <tr><td>1</td><td>neoantigen / vaccine / ICI triage engine</td><td>Best overlap between paperability, validation, and BD</td></tr>
      <tr><td>2</td><td>cell-therapy potency / QC scorer</td><td>Clear industrial use case and strong translational value</td></tr>
      <tr><td>3</td><td>PARP / kinase / endocrine / epigenetic response ranker</td><td>Best reusable cross-cancer validation layer</td></tr>
      <tr><td>4</td><td>bispecific / ADC / target-selection module</td><td>Useful commercial and mechanism module, but should stay modular</td></tr>
      <tr><td>5</td><td>oncolytic / radioligand reserve module</td><td>Interesting but lower-density and more heterogeneous</td></tr>
    </tbody>
  </table>
  <div class="box good"><b>Bottom line.</b> GEO is a source for building a therapy triage engine, not a magic therapy discovery engine. The winning product is a reproducible decision system with modality-specific routes.</div>
</section>
</main>
</div>
</body>
</html>"""


def copy_outputs(paths: list[Path]) -> None:
    for path in paths:
        shutil.copy2(path, ASSET_DIR / path.name)
        shutil.copy2(path, LIVE_ASSET_DIR / path.name)


def main() -> None:
    ensure_dirs()
    rows = MODALITIES
    tsv = RESULT / "therapy_spectrum_triage.tsv"
    queries = RESULT / "therapy_spectrum_query_templates.tsv"
    summary = RESULT / "SUMMARY.md"
    html_path = HUB / "therapy_spectrum_triage_2026_05_11.html"

    write_tsv(tsv, rows)
    write_query_tsv(queries, rows)
    write_summary(summary, rows)
    html_path.write_text(build_html(rows), encoding="utf-8")

    copy_outputs([tsv, queries, summary])
    shutil.copy2(html_path, LIVE_HUB / html_path.name)


if __name__ == "__main__":
    main()
