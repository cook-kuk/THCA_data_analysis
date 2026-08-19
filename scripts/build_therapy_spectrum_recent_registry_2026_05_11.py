from __future__ import annotations

import csv
import html
import shutil
from pathlib import Path


ROOT = Path("/home/seungho/personal/THCA_data_analysis")
RESULT = ROOT / "project/results/therapy_spectrum_recent_registry_2026_05_11"
HUB = ROOT / "project/papers_hub_2026_05_04"
LIVE_HUB = Path("/var/www/papers/papers_hub_2026_05_04")
ASSET_DIR = HUB / "assets/therapy_spectrum_recent_registry_2026_05_11"
LIVE_ASSET_DIR = LIVE_HUB / "assets/therapy_spectrum_recent_registry_2026_05_11"


ROWS = [
    {
        "accession": "GSE308384",
        "title": "Dextran-based T-cell expansion nanoparticles for manufacturing CAR T cells with augmented efficacy",
        "date": "2025-11-19",
        "modality": "CAR-T manufacturing",
        "status": "GO",
        "species": "Homo sapiens",
        "n": "NA",
        "best_task": "manufacturing QC / expansion / process optimization",
        "task_fit": 5,
        "paperability": 4,
        "validation": 5,
        "bd": 5,
        "why": "clear process metric and direct industrial relevance; expansion nanoparticles are close to a deployable workflow",
        "caveat": "process paper, not proof of superior clinical efficacy",
        "source": "turn4view0 L23-L28",
    },
    {
        "accession": "GSE270430",
        "title": "Allogeneic HSPC-engineered CD33-targeting CAR-NKT cells synergize with hypomethylating agents for effective and safe treatment of myeloid malignancies",
        "date": "2025-01-01",
        "modality": "CAR-NKT",
        "status": "GO",
        "species": "Homo sapiens",
        "n": "NA",
        "best_task": "cell-therapy prioritization / off-the-shelf design",
        "task_fit": 5,
        "paperability": 5,
        "validation": 4,
        "bd": 5,
        "why": "off-the-shelf CAR-NKT plus HMA synergy is a strong translational and BD hook",
        "caveat": "still preclinical; selection and safety need orthogonal validation",
        "source": "turn4view1 L23-L28",
    },
    {
        "accession": "GSE312384",
        "title": "BACH2 regulates T cell lineage states to enhance CAR T cell function [ATAC-seq]",
        "date": "2025-12-08",
        "modality": "CAR-T lineage state",
        "status": "GO",
        "species": "Homo sapiens",
        "n": "15",
        "best_task": "lineage-state / exhaustion rescue / epigenetic regulation",
        "task_fit": 5,
        "paperability": 5,
        "validation": 4,
        "bd": 4,
        "why": "BACH2 is a central lever for CAR function and clinical outcome association appears in the series",
        "caveat": "ATAC-seq mechanism support, not direct patient release assay",
        "source": "turn2view2 L28-L47",
    },
    {
        "accession": "GSE236468",
        "title": "Restoration of LAT signaling through the expression of a novel Adjunctive LAT-Activating CAR (ALA-CART) enhances the antigen-sensitivity and persistence of CAR T cells",
        "date": "2025-02-04",
        "modality": "CAR-T signaling",
        "status": "GO",
        "species": "Homo sapiens",
        "n": "NA",
        "best_task": "antigen-sensitivity / persistence / signal engineering",
        "task_fit": 5,
        "paperability": 5,
        "validation": 4,
        "bd": 4,
        "why": "directly addresses antigen-low relapse, a classic translational bottleneck",
        "caveat": "effects are model-dependent and need safety framing",
        "source": "turn4view3 L23-L33",
    },
    {
        "accession": "GSE292859",
        "title": "CRISPR engineering of armored CAR T cells enables tumor-restricted payload delivery with enhanced safety and efficacy",
        "date": "2025-06-18",
        "modality": "armored CAR-T",
        "status": "GO",
        "species": "Homo sapiens",
        "n": "2 per group",
        "best_task": "payload localization / tumor-restricted delivery / safety",
        "task_fit": 5,
        "paperability": 4,
        "validation": 4,
        "bd": 5,
        "why": "tumor-localized cytokine delivery is highly BD-friendly and process-adjacent",
        "caveat": "n is small; keep it as design support rather than definitive validation",
        "source": "turn4view2 L23-L28",
    },
    {
        "accession": "GSE287125",
        "title": "Immunotherapy-related cognitive impairment after CAR T cell therapy in mice [scRNA-seq]",
        "date": "2025-02-18",
        "modality": "CAR-T safety / toxicity",
        "status": "WATCH",
        "species": "Homo sapiens",
        "n": "7",
        "best_task": "safety signal mapping / adverse-effect biomarker discovery",
        "task_fit": 4,
        "paperability": 4,
        "validation": 3,
        "bd": 3,
        "why": "valuable as a safety companion module for CAR-T programs",
        "caveat": "toxicity biology, not efficacy; do not mix with potency claims",
        "source": "turn4view4 L23-L31",
    },
    {
        "accession": "GSE291443",
        "title": "Engineering an in vivo charging station for CAR-redirected invariant natural killer T cells to enhance cancer therapy",
        "date": "2025-05-07",
        "modality": "CAR-iNKT",
        "status": "GO",
        "species": "Homo sapiens",
        "n": "NA",
        "best_task": "persistence / recruitment / local activation",
        "task_fit": 5,
        "paperability": 4,
        "validation": 4,
        "bd": 4,
        "why": "in vivo activation architecture is a strong cell-therapy platform concept",
        "caveat": "therapy is still a platform design, not a disease-specific clinical proof",
        "source": "turn4view5 L23-L30",
    },
    {
        "accession": "GSE304795",
        "title": "pTα enhances mRNA translation and potentiates CAR T cells for solid tumor eradication [eCLIP-seq]",
        "date": "2026-01-22",
        "modality": "CAR-T translation / solid tumor",
        "status": "GO",
        "species": "Homo sapiens",
        "n": "NA",
        "best_task": "translation control / persistence / solid tumor potency",
        "task_fit": 5,
        "paperability": 5,
        "validation": 4,
        "bd": 4,
        "why": "one of the strongest recent examples for engineering CAR-T into solid tumors",
        "caveat": "published late; keep claim boundary at platform support",
        "source": "turn2view8 L23-L30",
    },
    {
        "accession": "GSE313971",
        "title": "Rational redesign of antigen binding domain improves in vivo efficacy of the CD22-CAR",
        "date": "2025-12-17",
        "modality": "CD22-CAR tonic signaling",
        "status": "GO",
        "species": "Homo sapiens",
        "n": "12",
        "best_task": "tonic signaling / antigen-binding redesign / CAR fitness",
        "task_fit": 5,
        "paperability": 5,
        "validation": 5,
        "bd": 4,
        "why": "direct CAR design optimization with a concrete clinical failure mode",
        "caveat": "human sample size is modest; keep efficacy claims platform-level",
        "source": "turn3search0; turn3search4",
    },
    {
        "accession": "GSE284026",
        "title": "NR2F6 deletion revives CAR-T cell function and induces antigen-agnostic immune memory in solid tumors",
        "date": "2026-01-23",
        "modality": "NR2F6 CAR-T solid tumor",
        "status": "GO",
        "species": "Mus musculus",
        "n": "NA",
        "best_task": "solid-tumor persistence / immune memory / transcriptional brake removal",
        "task_fit": 5,
        "paperability": 5,
        "validation": 5,
        "bd": 4,
        "why": "recent solid-tumor CAR-T engineering with a durable memory readout",
        "caveat": "mouse system; use as mechanism and platform support, not clinical proof",
        "source": "turn3search2; turn3search6",
    },
    {
        "accession": "GSE303442",
        "title": "RELB Reprograms Exhausted Tumor-Infiltrating Lymphocytes for Improved Adoptive Cell Therapy",
        "date": "2025-12-02",
        "modality": "TIL adoptive therapy",
        "status": "GO",
        "species": "Homo sapiens",
        "n": "NA",
        "best_task": "TIL expansion / exhaustion rescue / product potency",
        "task_fit": 5,
        "paperability": 5,
        "validation": 5,
        "bd": 4,
        "why": "clean adoptive-cell-therapy productization signal with organoid and xenograft support",
        "caveat": "still a manufacturing/productivity layer, not direct clinical outcome evidence",
        "source": "turn3search3; turn3search8",
    },
    {
        "accession": "GSE290722",
        "title": "Longitudinal single-cell immunoprofiling links durable CAR T response to sustained activation and clonotypic expansion of the native cytotoxic T cell repertoire",
        "date": "2025-03-04",
        "modality": "CAR-T response durability",
        "status": "GO",
        "species": "Homo sapiens",
        "n": "32",
        "best_task": "response durability / native repertoire expansion / outcome prediction",
        "task_fit": 5,
        "paperability": 5,
        "validation": 5,
        "bd": 5,
        "why": "longitudinal patient data ties native T-cell dynamics to durable CAR-T response",
        "caveat": "LBCL-specific; use as response framework, not universal CAR-T biology",
        "source": "turn1search1; turn1search7; turn1search9",
    },
    {
        "accession": "GSE311267",
        "title": "Syndecan-1-targeted therapeutic antibody inhibits macropinocytosis and induces antitumor immunity in pancreatic cancer",
        "date": "2025-12-02",
        "modality": "therapeutic antibody / PDAC",
        "status": "GO",
        "species": "Homo sapiens",
        "n": "NA",
        "best_task": "antibody target prioritization / KRAS-driven PDAC combination logic",
        "task_fit": 5,
        "paperability": 5,
        "validation": 4,
        "bd": 5,
        "why": "human-specific antibody with combination synergy against KRAS-driven PDAC",
        "caveat": "PDAC-specific biology; adapt carefully outside KRAS/SDC1 contexts",
        "source": "turn1search0; turn1search3",
    },
    {
        "accession": "GSE274141",
        "title": "Single-Cell RNA Sequencing Identifies Molecular Biomarkers Predicting Late Progression to CDK4/6 Inhibition in Patients with HR+/HER2- Metastatic Breast Cancer [FFPE]",
        "date": "2025-02-18",
        "modality": "CDK4/6 biomarker / endocrine resistance",
        "status": "GO",
        "species": "Homo sapiens",
        "n": "455",
        "best_task": "late-progression biomarker / endocrine resistance / PFS prediction",
        "task_fit": 5,
        "paperability": 5,
        "validation": 5,
        "bd": 4,
        "why": "large human FFPE cohort with direct progression biomarker utility",
        "caveat": "breast-cancer specific; keep claims tied to CDK4/6/endocrine context",
        "source": "turn1search2; turn1search10",
    },
    {
        "accession": "GSE298158",
        "title": "Single dose of CD137-Antibody Drug Conjugate protects non-human primate allogeneic HCT recipients against acute GVHD",
        "date": "2025-05-27",
        "modality": "CD137-ADC",
        "status": "GO",
        "species": "Homo sapiens / Macaca mulatta / Mus musculus",
        "n": "NA",
        "best_task": "activated-T-cell depletion / safety / allogeneic transplant control",
        "task_fit": 4,
        "paperability": 4,
        "validation": 5,
        "bd": 4,
        "why": "ADC crosses human, NHP, and mouse and gives a strong safety-design readout",
        "caveat": "transplant lane, not oncology; keep separate from efficacy-only claims",
        "source": "turn3search2",
    },
    {
        "accession": "GSE294241",
        "title": "T-cell dysfunction during blinatumomab therapy in pediatric acute lymphoblastic leukemia",
        "date": "2025-04-17",
        "modality": "blinatumomab / therapy response",
        "status": "GO",
        "species": "Homo sapiens",
        "n": "11 patients",
        "best_task": "response durability / T-cell dysfunction / dosing window",
        "task_fit": 5,
        "paperability": 5,
        "validation": 4,
        "bd": 4,
        "why": "direct therapy-response study with longitudinal patient sampling and a practical dosing question",
        "caveat": "disease-specific and pediatric; do not overgeneralize beyond B-ALL",
        "source": "turn3search1; turn3search9",
    },
    {
        "accession": "GSE303153",
        "title": "Synergistic Effect of Nucleoside Modification and Ionizable Lipid Composition on Translation and Immune Responses to mRNA Vaccines",
        "date": "2025-09-11",
        "modality": "mRNA vaccine platform",
        "status": "GO",
        "species": "Macaca fascicularis / Homo sapiens",
        "n": "292",
        "best_task": "vaccine translation / innate response / formulation ranking",
        "task_fit": 5,
        "paperability": 5,
        "validation": 5,
        "bd": 5,
        "why": "very large, translational mRNA-LNP dataset with direct formulation readouts",
        "caveat": "not cancer-specific; use as platform chemistry support",
        "source": "turn3view0 L23-L29; turn3view3 L32-L32",
    },
    {
        "accession": "GSE283041",
        "title": "Vaccination against oncofetal targets in the tumor vasculature",
        "date": "2025-04-09",
        "modality": "cancer vaccine / vasculature",
        "status": "WATCH",
        "species": "Mus musculus",
        "n": "NA",
        "best_task": "target selection / vascular immunotherapy",
        "task_fit": 4,
        "paperability": 4,
        "validation": 3,
        "bd": 4,
        "why": "useful bridge between vaccine logic and targetable tumor vasculature",
        "caveat": "mouse-only, so keep it in watch lane unless replicated",
        "source": "turn2view6 L23-L29",
    },
    {
        "accession": "GSE244034",
        "title": "Tumor-targeting oncolytic viruses engineered for glioblastoma immunotherapy",
        "date": "2025-07-14",
        "modality": "oncolytic virus",
        "status": "RESERVE",
        "species": "Mus musculus",
        "n": "NA",
        "best_task": "permissivity / immune remodeling / viral susceptibility",
        "task_fit": 3,
        "paperability": 3,
        "validation": 3,
        "bd": 3,
        "why": "interesting but still model-heavy and heterogeneous",
        "caveat": "reserve module only; not a core platform bet",
        "source": "turn4view6 L23-L29",
    },
    {
        "accession": "GSE288846",
        "title": "Engineering a novel HSV-1 strain for oncolytic therapy of solid tumors",
        "date": "2025-04-09",
        "modality": "oncolytic virus",
        "status": "RESERVE",
        "species": "Homo sapiens",
        "n": "NA",
        "best_task": "viral engineering / safety / tumor targeting",
        "task_fit": 3,
        "paperability": 3,
        "validation": 3,
        "bd": 3,
        "why": "engineering-heavy and interesting, but still not the main platform lane",
        "caveat": "keep as reserve unless a specific therapeutic angle is needed",
        "source": "turn4view8 L23-L28",
    },
    {
        "accession": "GSE288907",
        "title": "Viro-immunotherapy using oncolytic viruses represents a promising treatment for PDAC",
        "date": "2025-02-11",
        "modality": "oncolytic virus / PDAC",
        "status": "WATCH",
        "species": "NA",
        "n": "10 PDOs",
        "best_task": "patient-derived sensitivity / virotherapy ranking",
        "task_fit": 4,
        "paperability": 4,
        "validation": 3,
        "bd": 3,
        "why": "PDO-based drug-style readout is more screening-friendly than mouse-only OV papers",
        "caveat": "species/model details must be verified before reuse",
        "source": "turn1search5 snippet",
    },
    {
        "accession": "GSE292621",
        "title": "Single cell multi-omics reveals re-dosing with CD3 bispecific antibody induces a TCF7 high central memory CD8+ T cell population associated with reduced cytokine production",
        "date": "2025-03-28",
        "modality": "bispecific / CD3 engager",
        "status": "GO",
        "species": "Homo sapiens",
        "n": "NA",
        "best_task": "dose scheduling / memory-state shift / cytokine risk",
        "task_fit": 5,
        "paperability": 5,
        "validation": 4,
        "bd": 4,
        "why": "directly informs re-dosing strategy and cytokine-risk management",
        "caveat": "use as IO schedule support, not universal efficacy proof",
        "source": "turn4view9 L23-L28",
    },
    {
        "accession": "GSE289437",
        "title": "PARP inhibitors elicit distinct transcriptional programs in homologous recombination repair competent castrate resistant prostate cancer",
        "date": "2025-07-21",
        "modality": "PARP inhibitor",
        "status": "GO",
        "species": "Homo sapiens",
        "n": "NA",
        "best_task": "resistance / sensitization / biomarker discovery",
        "task_fit": 5,
        "paperability": 5,
        "validation": 5,
        "bd": 4,
        "why": "clean perturbation logic and strong cross-cohort reuse potential",
        "caveat": "commercial relevance depends on a narrow biomarker route",
        "source": "turn4view7 L23-L27",
    },
    {
        "accession": "GSE290830",
        "title": "Harnessing STING Signaling and Natural Killer Cells to overcome PARP Inhibitor Resistance in Homologous Recombination Repair altered Breast Cancer [scRNA-seq]",
        "date": "2025-03-03",
        "modality": "PARP combo / STING / NK",
        "status": "GO",
        "species": "NA",
        "n": "NA",
        "best_task": "combo prioritization / immune rescue / resistance reversal",
        "task_fit": 5,
        "paperability": 5,
        "validation": 4,
        "bd": 4,
        "why": "exactly the sort of resistance-to-combo bridge that can become a decision engine",
        "caveat": "species and experimental system should be checked before a claim",
        "source": "turn1search10 snippet",
    },
]


def score(row: dict[str, object]) -> int:
    return int(row["paperability"]) + int(row["validation"]) + int(row["bd"])


def ensure_dirs() -> None:
    RESULT.mkdir(parents=True, exist_ok=True)
    ASSET_DIR.mkdir(parents=True, exist_ok=True)
    LIVE_ASSET_DIR.mkdir(parents=True, exist_ok=True)
    HUB.mkdir(parents=True, exist_ok=True)
    LIVE_HUB.mkdir(parents=True, exist_ok=True)


def esc(s: object) -> str:
    return html.escape(str(s), quote=True)


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    headers = [
        "accession",
        "title",
        "date",
        "modality",
        "status",
        "species",
        "n",
        "score_total",
        "best_task",
        "task_fit",
        "paperability",
        "validation",
        "bd",
        "why",
        "caveat",
        "source",
    ]
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=headers, delimiter="\t")
        w.writeheader()
        for r in rows:
            out = dict(r)
            out["score_total"] = score(out)
            w.writerow(out)


def write_summary(path: Path, rows: list[dict[str, object]]) -> None:
    go = sum(1 for r in rows if r["status"] == "GO")
    watch = sum(1 for r in rows if r["status"] == "WATCH")
    reserve = sum(1 for r in rows if r["status"] == "RESERVE")
    top = sorted(rows, key=lambda r: (-score(r), -int(r["paperability"]), r["accession"]))[:5]
    lines = [
        "# Therapy Spectrum Recent Registry",
        "",
        "Date: 2026-05-11",
        "",
        "This registry collects recent GEO therapy-related records surfaced from current queries and ranks them by platform usefulness.",
        "",
        "Snapshot:",
        f"- Rows: {len(rows)}",
        f"- GO / WATCH / RESERVE: {go} / {watch} / {reserve}",
        f"- Top five: {top[0]['accession']} ({top[0]['modality']}); {top[1]['accession']} ({top[1]['modality']}); {top[2]['accession']} ({top[2]['modality']}); {top[3]['accession']} ({top[3]['modality']}); {top[4]['accession']} ({top[4]['modality']})",
        "",
        "Rule:",
        "- Prefer intervention + response readout + a deployable decision artifact.",
        "- Use watch/reserve for mouse-only, toxicity-only, or platform-engineering studies that need more validation.",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_table(rows: list[dict[str, object]]) -> str:
    ordered = sorted(rows, key=lambda r: (-score(r), -int(r["paperability"]), r["accession"]))
    trs = []
    for r in ordered:
        cls = {"GO": "go", "WATCH": "watch", "RESERVE": "reserve"}[r["status"]]
        trs.append(
            "<tr>"
            f"<td><a href='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={esc(r['accession'])}'>{esc(r['accession'])}</a></td>"
            f"<td>{esc(r['title'])}</td>"
            f"<td>{esc(r['date'])}</td>"
            f"<td>{esc(r['modality'])}</td>"
            f"<td><span class='tag {cls}'>{esc(r['status'])}</span></td>"
            f"<td>{esc(r['species'])}</td>"
            f"<td>{esc(r['n'])}</td>"
            f"<td>{score(r)}</td>"
            f"<td>{esc(r['best_task'])}</td>"
            f"<td>{r['task_fit']}</td>"
            f"<td>{r['paperability']}</td>"
            f"<td>{r['validation']}</td>"
            f"<td>{r['bd']}</td>"
            f"<td>{esc(r['why'])}</td>"
            f"<td>{esc(r['caveat'])}</td>"
            f"<td><code>{esc(r['source'])}</code></td>"
            "</tr>"
        )
    return "\n".join(trs)


def build_html(rows: list[dict[str, object]]) -> str:
    ordered = sorted(rows, key=lambda r: (-score(r), -int(r["paperability"]), r["accession"]))
    top = ordered[:6]
    cards = "".join(
        f"<div class='card'><h3>{i+1}. {esc(r['accession'])}</h3><p><b>{esc(r['modality'])}</b></p><p>Score {score(r)} | {esc(r['status'])} | {esc(r['date'])}</p><p>{esc(r['best_task'])}</p><p class='small'>{esc(r['why'])}</p></div>"
        for i, r in enumerate(top)
    )
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8" />
<meta name="viewport" content="width=device-width,initial-scale=1.0" />
<title>Therapy Spectrum Recent Registry · GEO Latest Scan</title>
<style>
:root{{--bg:#061019;--panel:#0d1726;--line:#23324a;--ink:#eaf1fb;--muted:#94a5ba;--gold:#ffd28a;--teal:#35d39d;--red:#ff8a6b}}
*{{box-sizing:border-box}}
html,body{{margin:0;background:var(--bg);color:var(--ink);font:14px/1.55 Inter,"Noto Sans KR",-apple-system,sans-serif}}
a{{color:var(--teal);text-decoration:none}} a:hover{{color:var(--gold);text-decoration:underline}}
code,pre,.mono{{font-family:"JetBrains Mono","SF Mono",Menlo,monospace;font-size:12px}}
h1,h2,h3{{font-family:"Cormorant Garamond","Newsreader",serif;letter-spacing:-.02em;color:#fff7dd}}
.hero{{padding:64px 34px 44px;background:radial-gradient(circle at top left,#162540 0%,#0a1220 60%,#050810 100%);border-bottom:1px solid var(--line)}}
.hero-inner{{max-width:1380px;margin:0 auto}}
.kicker{{font:700 11px "JetBrains Mono",monospace;color:var(--gold);letter-spacing:.24em;text-transform:uppercase}}
h1{{font-size:64px;line-height:.98;margin:12px 0 12px}}
.lead{{max-width:1040px;color:#cad5e5;font:18px/1.65 "Newsreader",serif}}
.stats{{display:grid;grid-template-columns:repeat(5,1fr);gap:10px;margin-top:24px}}
.stat{{background:rgba(255,255,255,.04);border:1px solid rgba(255,210,138,.18);padding:13px 14px;border-radius:10px}}
.stat b{{display:block;color:var(--gold);font:700 28px "Cormorant Garamond",serif;line-height:1}}
.stat span{{display:block;color:var(--muted);font:700 10px "JetBrains Mono",monospace;text-transform:uppercase;letter-spacing:.08em;margin-top:4px}}
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
    <div class="kicker">GEO latest scan · 2026-05-11</div>
    <h1>Recent GEO Therapy<br/><em>Candidate Registry</em></h1>
    <p class="lead">
      This registry captures recent GEO records surfaced from the current therapy-spectrum scan and ranks them by how directly they support a reproducible in silico decision system.
      The strongest records are CAR-T manufacturing/engineering, CAR-NKT, mRNA-LNP platform chemistry, and bispecific re-dosing state dynamics.
    </p>
    <div class="crumbs" style="margin-top:18px;color:#94a5ba;font:12px 'JetBrains Mono',monospace;">
      <a href="index.html" style="color:#35d39d">← Hub index</a> ·
      <a href="therapy_spectrum_triage_2026_05_11.html" style="color:#7eb6ff">therapy-spectrum triage</a> ·
      <a href="therapy_spectrum_go_shortlist_2026_05_11.html" style="color:#ffd28a">GO shortlist</a> ·
      <a href="therapy_spectrum_ready_split_2026_05_11.html" style="color:#f2c46d">paper/validation/BD split</a>
    </div>
    <div class="stats">
      <div class="stat"><b>{len(rows)}</b><span>Recent records</span></div>
      <div class="stat"><b>{sum(1 for r in rows if r['status']=='GO')}</b><span>GO</span></div>
      <div class="stat"><b>{sum(1 for r in rows if r['status']=='WATCH')}</b><span>WATCH</span></div>
      <div class="stat"><b>{sum(1 for r in rows if r['status']=='RESERVE')}</b><span>RESERVE</span></div>
      <div class="stat"><b>GEO</b><span>Latest public scan</span></div>
    </div>
  </div>
</header>
<div class="wrap">
<nav class="toc">
  <h4>Contents</h4>
  <a href="#tldr">01 TL;DR</a>
  <a href="#top">02 Top records</a>
  <a href="#table">03 Full table</a>
  <a href="#rules">04 Rules</a>
</nav>
<main>
<section id="tldr">
  <h2><span class="num">01</span>TL;DR</h2>
  <div class="sub">Recent records worth immediate inspection</div>
  <table class="t">
    <thead><tr><th>Accession</th><th>Modality</th><th>Score</th><th>Status</th><th>Date</th><th>Best task</th></tr></thead>
    <tbody>
      {''.join(f"<tr><td><a href='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={esc(r['accession'])}'>{esc(r['accession'])}</a></td><td>{esc(r['modality'])}</td><td>{score(r)}</td><td><span class='tag { {'GO':'go','WATCH':'watch','RESERVE':'reserve'}[r['status']] }'>{esc(r['status'])}</span></td><td>{esc(r['date'])}</td><td>{esc(r['best_task'])}</td></tr>" for r in ordered[:8])}
    </tbody>
  </table>
  <div class="box good"><b>Immediate emphasis:</b> CAR-T manufacturing and engineering papers are the strongest direct productization lane; mRNA-LNP chemistry is the strongest platform-chemistry lane; bispecific and PARP-combo records are the strongest response/ranking lanes.</div>
</section>

<section id="top">
  <h2><span class="num">02</span>Top Records</h2>
  <div class="grid">{cards}</div>
</section>

<section id="table">
  <h2><span class="num">03</span>Full Table</h2>
  <div class="sub">Ranked by score_total = paperability + validation + BD</div>
  <table class="t">
    <thead><tr><th>Accession</th><th>Title</th><th>Date</th><th>Modality</th><th>Status</th><th>Species</th><th>n</th><th>Score</th><th>Best task</th><th>Task fit</th><th>P</th><th>V</th><th>BD</th><th>Why</th><th>Caveat</th><th>Source</th></tr></thead>
    <tbody>{build_table(rows)}</tbody>
  </table>
</section>

<section id="rules">
  <h2><span class="num">04</span>Rules</h2>
  <div class="sub">How to promote candidates</div>
  <pre>GO:
 - human or human-relevant
 - intervention + response readout
 - clear decision artifact

WATCH:
 - mouse-only, toxicity-only, or platform-engineering
 - informative, but needs replication or orthogonal support

RESERVE:
 - biology is interesting, but the lane is too sparse / heterogeneous for a core bet</pre>
  <div class="box warn"><b>Boundary.</b> This registry is a selection aid, not a therapeutic endorsement. Use it to prioritize analysis and companion biomarker work, not to claim clinical efficacy from GEO alone.</div>
</section>
</main>
</div>
</body>
</html>"""


def main() -> None:
    ensure_dirs()
    tsv = RESULT / "therapy_spectrum_recent_registry.tsv"
    summary = RESULT / "SUMMARY.md"
    html_path = HUB / "therapy_spectrum_recent_registry_2026_05_11.html"

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
