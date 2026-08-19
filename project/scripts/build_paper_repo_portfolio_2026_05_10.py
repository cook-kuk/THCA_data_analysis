#!/usr/bin/env python3
"""Build a local repo-style paper portfolio map.

This script is intentionally local-only. It writes under project/results and
does not copy to /var/www.
"""

from __future__ import annotations

import csv
import html
import json
from datetime import datetime
from pathlib import Path
from typing import Iterable


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "project" / "results"
HUB = ROOT / "project" / "papers_hub_2026_0509"
OUT = RESULTS / "paper_repo_portfolio_2026_05_10"


BRANCHES = [
    {
        "branch": "paper1-main",
        "title": "Paper 1 / DM1 molecular dark matter",
        "status": "SUBMIT_NOW",
        "repo_role": "trunk",
        "claim": "Fusion-driven, TERT-linked thyroid dedifferentiation state anchored by the 8-gene DM1 panel.",
        "data": "TCGA-THCA RNA/HM450/clinical/SV, Korean K2, Lee 2024, GSE76039, Lu 2023, Pu 2021, DepMap/PRISM reserve.",
        "software": "Python, pandas, scipy, statsmodels, scikit-learn, NuSVR, matplotlib, static HTML dossier builder.",
        "interpretation": "Treat as correlative multi-cohort mechanism and reviewer-defense trunk; Q9 remains motivates-not-confirms.",
        "boundary": "No direct mechanism-of-silencing proof; no clinical diagnostic claim without prospective validation.",
        "tokens": ["paper1", "dm1", "deconv", "dark_matter", "v17", "k2", "tds", "fig8", "mechanism"],
        "parents": [],
        "children": ["paper1b-braf-axis", "paper9-ic50-perturbation", "paper11-pancancer", "paper2-wsi-image"],
        "priority": 1,
    },
    {
        "branch": "paper1b-braf-axis",
        "title": "Paper 1B / BRAF cPTC DM2 vs DM1xTERT axis",
        "status": "MERGE_OR_STANDALONE_DECISION",
        "repo_role": "feature branch",
        "claim": "BRAF-classic PTC axis can become a companion mechanism story if it adds survival, RAI, MAPK, immune and methylation depth.",
        "data": "TCGA BRAF/RAS/TERT, Lee 2024, HM450, spatial BRAF inference, single-cell support, ICI panels.",
        "software": "Python, survival models, Fisher tests, deconvolution, spatial scoring, HTML synthesis.",
        "interpretation": "Useful as Paper 1 enhancer or standalone only if it avoids competing with the DM1 trunk.",
        "boundary": "Do not let BRAF-axis dilute the main DM1 novelty unless it resolves a reviewer attack.",
        "tokens": ["braf", "dm2", "tert", "p2_braf", "mapk", "hm450", "rai_response"],
        "parents": ["paper1-main"],
        "children": ["paper9-ic50-perturbation"],
        "priority": 3,
    },
    {
        "branch": "paper2-wsi-image",
        "title": "Paper 2 / image-DM1 and Path2Space-style pathology spatial bridge",
        "status": "VIABLE_PILOT",
        "repo_role": "application branch",
        "claim": "H&E and spatial-pathology features can carry a DM1-like / thyroid-state signal, but current evidence is pilot-level.",
        "data": "TCGA WSI embedded slides, GSE230424, GSE250521, GSE248205, pathology tile atlases and spatial overlays.",
        "software": "DINOv2/UNI embeddings, CLAM, scikit-learn, spatial module scoring, RunPod GPU jobs.",
        "interpretation": "Show image/spatial feasibility; keep it as pilot unless slide/site/sample split robustness improves.",
        "boundary": "No clinical deployment, no universal pathology biomarker claim, no leakage-prone slide-level inflation.",
        "tokens": ["wsi", "image", "clam", "path2space", "gse230424", "gse250521", "spatial", "cv2"],
        "parents": ["paper1-main"],
        "children": ["paper2b-ht-immune", "paper10-12-network-atlas"],
        "priority": 4,
    },
    {
        "branch": "paper2b-ht-immune",
        "title": "Paper 2B / Hashimoto thyroiditis immune-HLA context",
        "status": "BIOLOGY_STRONG_ALLELE_WEAK",
        "repo_role": "sidecar branch",
        "claim": "HT/PTC immune context is real in expression/spatial layers, but HLA allele replication is underpowered.",
        "data": "TCGA HT labels, GSE286332, GSE248205, GSE250521, HLA deep-dive tracks.",
        "software": "arcasHLA outputs, Fisher/meta-analysis, spatial burden summaries, module scoring.",
        "interpretation": "Use as immune-context enhancer, not as a genetic association manuscript yet.",
        "boundary": "No Korean/allele-level HLA claim until adult Korean GD/PTC NGS exists.",
        "tokens": ["hashimoto", "ht", "hla", "gse286332", "aitd", "panasian"],
        "parents": ["paper2-wsi-image"],
        "children": ["paper4-gd-hla"],
        "priority": 7,
    },
    {
        "branch": "paper3-ici-readiness",
        "title": "Paper 3 / ICI vulnerability and immune-readiness map",
        "status": "FROZEN_DESIGN",
        "repo_role": "reserve branch",
        "claim": "DM1/ATC-like immune state can be organized against public ICI cohorts, but thyroid-specific response prediction is not yet proven.",
        "data": "Public ICI cohorts, ATC/THCA immune modules, HLA/IFNG tracks, phase2/phase4 ICI summaries.",
        "software": "Cohort pooling, module scoring, response association tests, forest summaries.",
        "interpretation": "Keep as readiness map and external context; do not oversell as predictor.",
        "boundary": "No thyroid ICI response model without matched thyroid ICI outcomes.",
        "tokens": ["ici", "phase_C_ICI", "paper3", "immune", "ifng", "track18"],
        "parents": ["paper1-main"],
        "children": [],
        "priority": 8,
    },
    {
        "branch": "paper4-gd-hla",
        "title": "Paper 4 / Korean GD and Pan-Asian HLA",
        "status": "UPGRADE_IF_DATA_ARRIVES",
        "repo_role": "data-gated branch",
        "claim": "Pan-Asian HLA anchors are plausible, but Korean replication is the gate.",
        "data": "Korean PTC pool, Chen 2018 Han Chinese GD, GSE286332 arcasHLA, HLA deep-dive tracks.",
        "software": "arcasHLA, allele-frequency harmonization, Fisher tests, meta-analysis forest plots.",
        "interpretation": "Useful as grant/data-request map now; manuscript only after Korean adult GD NGS.",
        "boundary": "No Korean GD association claim from proxy PTC/HT data.",
        "tokens": ["gd", "hla", "panasian", "arcas", "finn", "ukbb", "gse286332"],
        "parents": ["paper2b-ht-immune"],
        "children": [],
        "priority": 9,
    },
    {
        "branch": "paper9-ic50-perturbation",
        "title": "Paper 9 / IC50, PRISM, DepMap and perturbation vulnerability",
        "status": "WETLAB_REQUIRED",
        "repo_role": "therapeutic branch",
        "claim": "DM1-high thyroid/pan-cancer states show a MAPK-heavy drug-response proxy and selected metabolic dependencies.",
        "data": "PRISM drug response, DepMap CRISPR, thyroid cell line annotations, spatial drug target overlays.",
        "software": "PRISM differential ranking, DepMap dependency tests, target-state overlap matrices, spatial overlay scoring.",
        "interpretation": "This is an IC50/proxy map for hypothesis prioritization, not a validated therapeutic claim.",
        "boundary": "No synthetic-lethality or patient-treatment claim without wet-lab IC50/viability validation.",
        "tokens": ["ic50", "prism", "depmap", "synthetic", "lethality", "drug", "vulnerability", "perturbation"],
        "parents": ["paper1-main", "paper1b-braf-axis"],
        "children": ["paper11-pancancer"],
        "priority": 5,
    },
    {
        "branch": "paper10-12-network-atlas",
        "title": "Paper 10/12 / network and atlas companion",
        "status": "ABSORB_INTO_P11",
        "repo_role": "atlas branch",
        "claim": "Network, module and atlas views help explain portability and candidate submodules.",
        "data": "paper12 network edges/modules, cancer gene matrix, HLA atlas, source atlas, figure/table atlases.",
        "software": "Correlation networks, module summaries, static atlas pages, visualization tables.",
        "interpretation": "Better as navigation and supplement layer than as separate manuscript today.",
        "boundary": "Do not split into low-novelty standalone paper unless it gains external validation.",
        "tokens": ["paper12", "network", "atlas", "module", "matrix", "cancer_gene"],
        "parents": ["paper2-wsi-image", "paper11-pancancer"],
        "children": [],
        "priority": 10,
    },
    {
        "branch": "paper11-pancancer",
        "title": "Paper 11 / pan-cancer portable DM1 biology",
        "status": "NEXT_MAJOR_AFTER_P1",
        "repo_role": "release branch",
        "claim": "DM1-like dedifferentiation has portable pan-cancer structure with lineage-specific boundaries.",
        "data": "TCGA pan-cancer, DepMap, PRISM, hallmark scores, survival meta, Korean multi-cohort support.",
        "software": "Pan-cancer module scoring, survival meta-analysis, lineage mitigation, PRISM/MAPK enrichment.",
        "interpretation": "Strong next manuscript after Paper 1 if lineage-confounding mitigation remains central.",
        "boundary": "Must present lineage-specificity and residual confounding up front.",
        "tokens": ["paper11", "pancancer", "lineage", "mitigation", "nature_master", "phase_H"],
        "parents": ["paper1-main", "paper9-ic50-perturbation"],
        "children": ["paper10-12-network-atlas"],
        "priority": 2,
    },
    {
        "branch": "neoantigen-clean-neo",
        "title": "Neoantigen / CLEAN-Neo, CROSS-Neo and MD-audited candidate selection",
        "status": "SEPARATE_UNIVERSE",
        "repo_role": "separate repository",
        "claim": "Leakage-aware neoantigen benchmarking, abstention and candidate selection can become an independent methods paper.",
        "data": "TESLA, CEDAR, dbPepNeo2, NeoDB, McPAS, IEDB, MD/TCR/pMHC audit outputs.",
        "software": "CLEAN-Neo pipeline, QK/rule-gated models, ESM2 features, external predictors, PMTnet, OpenMM/MD audit.",
        "interpretation": "Separate from THCA DM1 unless used only as vaccine-method appendix.",
        "boundary": "No THCA vaccine-readiness claim; public predictors are comparator-only when licenses/data leakage are uncertain.",
        "tokens": ["neo", "neobench", "barneo", "cross_neo", "tcr", "pmhc", "md", "vaccine"],
        "parents": [],
        "children": [],
        "priority": 6,
    },
    {
        "branch": "paper13-spatial-drug",
        "title": "Paper 13 / spatial drug vulnerability and target-celltype overlap",
        "status": "VISUAL_STRONG_DATA_NEEDS_LOCK",
        "repo_role": "visual application branch",
        "claim": "Spatial maps can show drug-target/state overlap and candidate hotspot burden across thyroid samples.",
        "data": "GSE250521/GSE230424 spatial overlays, target gene expression, drug-state scorecards, CV2 validation ROI packs.",
        "software": "Spatial module scoring, target-state overlap, ROI gallery generation, static visual dossiers.",
        "interpretation": "Excellent visual application layer; needs locked quantitative controls before manuscript escalation.",
        "boundary": "Spatial co-localization is not drug efficacy; match/mismatch boards must stay explicit.",
        "tokens": ["spatial_drug", "target_celltype", "vulnerability", "cv2", "roi", "selpercatinib", "encorafenib"],
        "parents": ["paper9-ic50-perturbation", "paper2-wsi-image"],
        "children": [],
        "priority": 11,
    },
]


METHODS = [
    {
        "method": "8-gene DM1/module scoring",
        "software": "Python, pandas, scipy, statsmodels",
        "inputs": "RNA expression matrices, cohort metadata, driver annotations",
        "outputs": "DM1 score, cohort rank, contrast tables, survival/driver overlays",
        "interpretation": "Relative molecular-state axis; strongest when replicated across cohorts and split logic.",
        "failure_mode": "Batch, lineage and label leakage can make clean biology look stronger than it is.",
    },
    {
        "method": "Methylation-expression support",
        "software": "pandas, scipy correlation, statsmodels",
        "inputs": "TCGA HM450 probes, gene expression, DM calls",
        "outputs": "gene-level methylation-expression correlation and group contrasts",
        "interpretation": "Supports epigenetic compatibility, not causal silencing.",
        "failure_mode": "Probe mapping and tumor purity can confound direction.",
    },
    {
        "method": "Deconvolution / cell-state attribution",
        "software": "NuSVR-style deconvolution, marker panels, custom Python",
        "inputs": "bulk RNA, marker signatures, spatial/scRNA references",
        "outputs": "cell-state fractions, marker module overlays, MAPK/TDS panel relationships",
        "interpretation": "Explains which compartment might carry a signal.",
        "failure_mode": "Reference mismatch and collinearity; fractions are estimates, not sorted-cell truth.",
    },
    {
        "method": "Pathology foundation modeling",
        "software": "DINOv2/UNI, CLAM, scikit-learn, RunPod GPU",
        "inputs": "H&E whole-slide tiles, slide labels, spatial labels for bridge analyses",
        "outputs": "slide/tile embeddings, AUC, attention maps, pathology-spatial overlays",
        "interpretation": "Image feasibility and localization signal.",
        "failure_mode": "Slide/site leakage, small N, inflated AUC without leave-one-site/sample tests.",
    },
    {
        "method": "Spatial transcriptomics overlay",
        "software": "scanpy-like matrices, custom scoring, matplotlib",
        "inputs": "Visium/spatial spot matrices, histology coordinates, module genes",
        "outputs": "spot maps, hotspot burden, target-state overlap, stage trend tables",
        "interpretation": "Shows tissue localization and co-occurrence.",
        "failure_mode": "Spot resolution mixes cells; spatial co-occurrence is not causal interaction.",
    },
    {
        "method": "PRISM / DepMap IC50 proxy mapping",
        "software": "Python, pandas, scipy, multiple-testing correction",
        "inputs": "PRISM drug response, DepMap CRISPR, model lineage metadata, DM1 scores",
        "outputs": "drug rankings, MAPK enrichment, dependency tables, spatial target overlays",
        "interpretation": "Hypothesis-prioritization map for wet-lab IC50/viability testing.",
        "failure_mode": "Cell-line N and lineage confounding; proxy response is not clinical efficacy.",
    },
    {
        "method": "Neoantigen leakage-aware benchmarking",
        "software": "CLEAN-Neo, rule-gated/QK models, external predictors, ESM2, OpenMM/MD audit",
        "inputs": "TESLA/CEDAR/dbPepNeo2/NeoDB/McPAS/IEDB, TCR/pMHC structures, public predictors",
        "outputs": "locked split metrics, abstention queues, candidate scorecards, MD evidence cards",
        "interpretation": "Methods and candidate-prioritization framework.",
        "failure_mode": "Source overlap, allele leakage, license constraints and incomplete immunogenicity labels.",
    },
    {
        "method": "Static paper hub / dossier generation",
        "software": "Python HTML builders, CSS, TSV/JSON inventories",
        "inputs": "local result files, figure assets, paper metadata",
        "outputs": "paper1-style HTML dossiers, page inventories, figure/table/claim ledgers",
        "interpretation": "Review and navigation layer, not primary analysis.",
        "failure_mode": "Can overstate if not tied to source paths and claim boundaries.",
    },
]


EDGES = [
    ("paper1-main", "paper1b-braf-axis", "shared TCGA/DM axis; BRAF branch tests whether driver subtype explains or complements DM1"),
    ("paper1-main", "paper9-ic50-perturbation", "DM1/MAPK state seeds perturbation and PRISM/DepMap prioritization"),
    ("paper1-main", "paper11-pancancer", "DM1 trunk becomes pan-cancer portability hypothesis"),
    ("paper1-main", "paper2-wsi-image", "molecular DM1 labels seed image/pathology predictors"),
    ("paper2-wsi-image", "paper2b-ht-immune", "spatial/pathology immune context branches into HT/HLA questions"),
    ("paper2b-ht-immune", "paper4-gd-hla", "HT/HLA context defines Korean GD data gate"),
    ("paper9-ic50-perturbation", "paper13-spatial-drug", "drug-response proxies are localized with target-state spatial overlays"),
    ("paper9-ic50-perturbation", "paper11-pancancer", "PRISM/DepMap pan-cancer signal informs Paper 11 therapeutic appendix"),
    ("paper11-pancancer", "paper10-12-network-atlas", "network/atlas layer supports portability interpretation"),
]


TOPIC_QUERIES = {
    "paper1-main": ["p_deconv_2026_05_08", "dm1_robustness_v2026_05_08", "v17", "audit_2026_04_30", "p_external_expression_validation"],
    "paper1b-braf-axis": ["p2_braf_nature_sprint_2026_05_09", "braf"],
    "paper2-wsi-image": ["p2_image_dm1_v2_foundation_clam_2026_05_07", "p_wsi_dm1_2026_05_09", "path2space"],
    "paper2b-ht-immune": ["htptc", "hashimoto", "hla_deepdive", "d4p2"],
    "paper3-ici-readiness": ["paper3_ici", "phase_C_ICI", "track18_ici"],
    "paper4-gd-hla": ["gd", "panasian", "track1_paper4", "d4p1"],
    "paper9-ic50-perturbation": ["p_synthetic_lethality_2026_05_09", "prism", "depmap", "paper9"],
    "paper10-12-network-atlas": ["paper12_network", "atlas", "network"],
    "paper11-pancancer": ["paper11_pancancer", "nature_master"],
    "neoantigen-clean-neo": ["p_neo_bayesian_2026_05_09", "clean_neobench", "cross_neo", "cancer_vaccine"],
    "paper13-spatial-drug": ["spatial_drug", "target_celltype", "vulnerability"],
}


def e(value: object) -> str:
    return html.escape(str(value), quote=True)


def read_inventory() -> list[dict[str, str]]:
    json_path = HUB / "page_inventory_2026_0509.json"
    if not json_path.exists():
        return []
    try:
        data = json.loads(json_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return []
    if isinstance(data, list):
        return [x for x in data if isinstance(x, dict)]
    if isinstance(data, dict):
        pages = data.get("pages", [])
        return [x for x in pages if isinstance(x, dict)]
    return []


def all_candidate_files() -> list[Path]:
    keep_suffix = {".md", ".json", ".tsv", ".csv", ".png", ".pdf", ".html"}
    files: list[Path] = []
    for path in RESULTS.rglob("*"):
        if path.is_file() and path.suffix.lower() in keep_suffix:
            files.append(path)
    return files


def rel(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def matched_files(branch: dict[str, object], files: list[Path], limit: int = 18) -> list[str]:
    queries = [str(x).lower() for x in TOPIC_QUERIES.get(str(branch["branch"]), [])]
    tokens = [str(x).lower() for x in branch.get("tokens", [])]
    hits: list[tuple[int, str]] = []
    for path in files:
        s = rel(path).lower()
        score = 0
        for q in queries:
            if q in s:
                score += 8
        for t in tokens:
            if t and t in s:
                score += 1
        if path.name in {"SUMMARY.md"}:
            score += 4
        if "summary" in path.name.lower():
            score += 2
        if "report" in path.name.lower():
            score += 1
        if score:
            hits.append((-score, rel(path)))
    hits.sort()
    dedup: list[str] = []
    seen = set()
    for _, name in hits:
        if name not in seen:
            dedup.append(name)
            seen.add(name)
        if len(dedup) >= limit:
            break
    return dedup


def page_matches(branch: dict[str, object], pages: Iterable[dict[str, str]], limit: int = 10) -> list[dict[str, str]]:
    tokens = [str(x).lower() for x in branch.get("tokens", [])]
    branch_name = str(branch["branch"]).lower()
    matches: list[tuple[int, dict[str, str]]] = []
    for page in pages:
        name = str(page.get("page") or page.get("file") or page.get("path") or page.get("href") or "").lower()
        title = str(page.get("title") or "").lower()
        blob = f"{name} {title}"
        score = 0
        for token in tokens:
            if token in blob:
                score += 2
        if branch_name.replace("-", "_") in blob:
            score += 4
        if score:
            matches.append((-score, page))
    matches.sort(key=lambda x: (x[0], str(x[1].get("page") or x[1].get("file") or x[1].get("path") or "")))
    return [p for _, p in matches[:limit]]


def write_tsv(path: Path, rows: list[dict[str, object]], columns: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({col: row.get(col, "") for col in columns})


def table(headers: list[str], rows: Iterable[Iterable[object]]) -> str:
    head = "".join(f"<th>{e(h)}</th>" for h in headers)
    body = []
    for row in rows:
        body.append("<tr>" + "".join(f"<td>{e(cell)}</td>" for cell in row) + "</tr>")
    return f"<table><thead><tr>{head}</tr></thead><tbody>{''.join(body)}</tbody></table>"


def css() -> str:
    return """
    :root{--bg:#0b1118;--panel:#111a24;--panel2:#162231;--ink:#f4efe4;--muted:#a8b4c0;--line:#2b3d50;--gold:#e8bf64;--cyan:#73d6ff;--red:#ff8a7a;--green:#9bdf9f}
    *{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--ink);font-family:Inter,Arial,sans-serif;line-height:1.58}
    a{color:var(--cyan);text-decoration:none}a:hover{text-decoration:underline}
    .wrap{max-width:1320px;margin:0 auto;padding:30px 22px 70px}
    .hero{padding:44px 0 24px;border-bottom:1px solid var(--line)}
    .kicker{color:var(--gold);font-size:12px;letter-spacing:.08em;text-transform:uppercase;font-weight:800}
    h1{font-family:Georgia,serif;font-size:clamp(36px,5vw,72px);line-height:1.02;margin:10px 0 14px;letter-spacing:0}
    h2{font-family:Georgia,serif;font-size:30px;margin:34px 0 12px;border-top:1px solid var(--line);padding-top:22px}
    h3{font-size:18px;margin:22px 0 8px;color:var(--gold)}
    .lead{max-width:980px;color:#d9e2ea;font-size:18px}
    .stats{display:grid;grid-template-columns:repeat(6,minmax(0,1fr));gap:10px;margin-top:24px}
    .stat{background:linear-gradient(180deg,#162231,#101923);border:1px solid var(--line);border-radius:8px;padding:13px}
    .stat b{display:block;font-size:24px;color:white}.stat span{color:var(--muted);font-size:12px}
    .grid{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:14px}
    .card{background:var(--panel);border:1px solid var(--line);border-radius:8px;padding:16px}
    .card strong{color:#fff}.card .status{display:inline-block;color:#101923;background:var(--gold);border-radius:4px;padding:2px 6px;font-weight:800;font-size:11px;margin-bottom:8px}
    .muted{color:var(--muted)}.warn{color:var(--red);font-weight:700}.ok{color:var(--green);font-weight:700}
    table{width:100%;border-collapse:collapse;background:var(--panel);border:1px solid var(--line);border-radius:8px;overflow:hidden;margin:12px 0 22px}
    th,td{padding:10px 11px;border-bottom:1px solid var(--line);vertical-align:top;font-size:13px}
    th{background:#192738;color:var(--gold);text-align:left;font-size:12px;text-transform:uppercase}
    tr:last-child td{border-bottom:0}.mono{font-family:JetBrains Mono,Consolas,monospace;font-size:12px}
    .flow{display:flex;flex-wrap:wrap;gap:8px;margin:12px 0 20px}.pill{border:1px solid var(--line);background:var(--panel2);border-radius:6px;padding:8px 10px;font-size:13px}
    .src{max-height:180px;overflow:auto;background:#0f1721;border:1px solid var(--line);border-radius:8px;padding:10px}
    @media(max-width:900px){.stats,.grid{grid-template-columns:1fr}.wrap{padding:22px 14px 50px}h1{font-size:38px}}
    """


def page_shell(title: str, lead: str, body: str, stats: list[tuple[str, str]]) -> str:
    stat_html = "".join(f"<div class='stat'><b>{e(v)}</b><span>{e(k)}</span></div>" for k, v in stats)
    nav = (
        "<div class='flow'>"
        "<a class='pill' href='paper_repo_map_v2_2026_05_10.html'>repo map</a>"
        "<a class='pill' href='software_method_ledger_v2_2026_05_10.html'>software ledger</a>"
        "<a class='pill' href='ic50_prediction_map_v2_2026_05_10.html'>IC50 map</a>"
        "<a class='pill' href='deconvolution_repo_v2_2026_05_10.html'>deconvolution</a>"
        "<a class='pill' href='neoantigen_repo_v2_2026_05_10.html'>neoantigen</a>"
        "</div>"
    )
    return f"""<!doctype html>
<html lang="ko"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>{e(title)}</title><style>{css()}</style></head>
<body><main class="wrap"><section class="hero"><div class="kicker">THCA paper repo portfolio / local build 2026-05-10</div>
<h1>{e(title)}</h1><p class="lead">{e(lead)}</p>{stat_html}</section>{nav}{body}</main></body></html>"""


def branch_card(branch: dict[str, object], files: list[str], pages: list[dict[str, str]]) -> str:
    links = []
    for page in pages[:6]:
        href = page.get("page") or page.get("file") or page.get("path") or page.get("href") or ""
        title = page.get("title") or href
        if href:
            links.append(f"<a href='../../papers_hub_2026_0509/{e(href)}'>{e(title)}</a>")
    page_html = "<br>".join(links) if links else "<span class='muted'>No page inventory match yet</span>"
    src_html = "<br>".join(e(x) for x in files[:8]) if files else "<span class='muted'>No source match</span>"
    return f"""
    <article class="card">
      <div class="status">{e(branch['status'])}</div>
      <h3>{e(branch['branch'])}</h3>
      <p><strong>{e(branch['title'])}</strong></p>
      <p>{e(branch['claim'])}</p>
      <p class="muted"><b>Data:</b> {e(branch['data'])}</p>
      <p class="muted"><b>Software:</b> {e(branch['software'])}</p>
      <p><b>Meaning:</b> {e(branch['interpretation'])}</p>
      <p class="warn">Boundary: {e(branch['boundary'])}</p>
      <div class="src mono">{src_html}</div>
      <p class="mono">{page_html}</p>
    </article>"""


def build_main(branch_rows: list[dict[str, object]], method_rows: list[dict[str, object]], page_count: int, file_count: int) -> str:
    cards = "".join(
        branch_card(row, row.get("source_files", []), row.get("hub_pages", []))
        for row in sorted(branch_rows, key=lambda x: int(x["priority"]))
    )
    edge_table = table(["parent", "child", "why this dependency matters"], EDGES)
    branch_table = table(
        ["priority", "branch", "status", "repo role", "claim boundary"],
        [[b["priority"], b["branch"], b["status"], b["repo_role"], b["boundary"]] for b in sorted(branch_rows, key=lambda x: int(x["priority"]))],
    )
    method_table = table(
        ["method", "software", "inputs", "how to interpret"],
        [[m["method"], m["software"], m["inputs"], m["interpretation"]] for m in method_rows],
    )
    body = f"""
    <h2>One Repo Map</h2>
    <p>이 지도는 paper 후보를 git repo처럼 trunk, feature branch, release branch, data-gated branch로 나눈다. 각 branch는 부모/자식, source files, hub pages, data, software, 해석 경계가 같이 붙어야 한다.</p>
    {branch_table}
    <h2>Branch History / Dependency Graph</h2>
    {edge_table}
    <h2>Paper Branch Cards</h2>
    <div class="grid">{cards}</div>
    <h2>Software Interpretation Ledger</h2>
    {method_table}
    <h2>Operational Rule</h2>
    <p><span class="ok">Go:</span> Paper 1 trunk, Paper 11 next-major, Paper 9 wet-lab-priority map, Neo separate methods repo. <span class="warn">Hold:</span> GD-HLA Korean claim, thyroid ICI predictor, clinical WSI deployment, spatial drug efficacy.</p>
    """
    return page_shell(
        "Paper Repo Map v2",
        "현재 results와 hub를 기준으로 paper 후보를 git-style branch/history/claim-boundary 구조로 재정리한 local-only 지도.",
        body,
        [("branches", str(len(branch_rows))), ("hub pages scanned", str(page_count)), ("source files scanned", str(file_count)), ("software methods", str(len(method_rows))), ("edges", str(len(EDGES))), ("deploy", "no")],
    )


def build_software_page(method_rows: list[dict[str, object]]) -> str:
    body = table(
        ["method", "software", "inputs", "outputs", "meaning", "failure mode"],
        [[m["method"], m["software"], m["inputs"], m["outputs"], m["interpretation"], m["failure_mode"]] for m in method_rows],
    )
    body += "<h2>How To Read This</h2><p>소프트웨어 이름은 주장의 강도가 아니라 계산 계층을 뜻한다. every result must cite input data, transformation, output table/figure and claim boundary.</p>"
    return page_shell(
        "Software / Method Ledger",
        "각 paper branch에서 어떤 데이터가 어떤 소프트웨어를 통과해 어떤 결과와 주장을 만드는지 해석 규칙까지 묶은 장부.",
        body,
        [("methods", str(len(method_rows))), ("primary use", "interpretation"), ("deploy", "no")],
    )


def build_topic_page(key: str, title: str, lead: str, branch_rows: list[dict[str, object]], method_rows: list[dict[str, object]]) -> str:
    selected = [b for b in branch_rows if b["branch"] == key or key in b.get("parents", []) or key in b.get("children", [])]
    cards = "".join(branch_card(b, b.get("source_files", []), b.get("hub_pages", [])) for b in selected)
    relevant_methods = []
    keywords = {
        "paper9-ic50-perturbation": ["PRISM", "DepMap", "Spatial"],
        "paper1-main": ["Deconvolution", "8-gene", "Methylation", "Spatial"],
        "neoantigen-clean-neo": ["Neoantigen"],
    }.get(key, [])
    for method in method_rows:
        blob = json.dumps(method, ensure_ascii=False)
        if any(k.lower() in blob.lower() for k in keywords):
            relevant_methods.append(method)
    body = f"<h2>Branch Scope</h2><div class='grid'>{cards}</div>"
    if relevant_methods:
        body += "<h2>Method Contract</h2>" + table(
            ["method", "software", "inputs", "outputs", "meaning", "failure mode"],
            [[m["method"], m["software"], m["inputs"], m["outputs"], m["interpretation"], m["failure_mode"]] for m in relevant_methods],
        )
    body += "<h2>Decision Boundary</h2>"
    if key == "paper9-ic50-perturbation":
        body += "<p>IC50/PRISM/DepMap branch는 wet-lab 후보 선별 지도다. MAPK-heavy signal은 strong prioritization이고, synthetic lethality/therapy claim은 wet-lab IC50 and viability validation 뒤에만 가능하다.</p>"
    elif key == "paper1-main":
        body += "<p>Deconvolution branch는 DM1 mechanism defense의 설명 계층이다. cell fraction/state estimates는 marker/reference-dependent estimates이므로 co-located caveat와 source tables가 항상 같이 있어야 한다.</p>"
    elif key == "neoantigen-clean-neo":
        body += "<p>Neoantigen branch는 THCA paper가 아니라 별도 methods repo다. Leakage-aware split, abstention, MD/TCR evidence는 candidate prioritization이지 immunogenicity proof가 아니다.</p>"
    else:
        body += "<p>Boundary must stay attached to every claim.</p>"
    return page_shell(title, lead, body, [("branches", str(len(selected))), ("deploy", "no"), ("source mode", "local")])


def build_markdown(branch_rows: list[dict[str, object]], method_rows: list[dict[str, object]], page_count: int, file_count: int) -> str:
    lines = [
        "# Paper Repo Map v2 (2026-05-10)",
        "",
        "Local-only build. No /var/www deployment was performed.",
        "",
        f"- Branches: {len(branch_rows)}",
        f"- Hub pages scanned: {page_count}",
        f"- Results files scanned: {file_count}",
        f"- Software/method contracts: {len(method_rows)}",
        "",
        "## Branches",
        "",
        "| priority | branch | status | role | data | software | boundary |",
        "|---:|---|---|---|---|---|---|",
    ]
    for b in sorted(branch_rows, key=lambda x: int(x["priority"])):
        lines.append(f"| {b['priority']} | {b['branch']} | {b['status']} | {b['repo_role']} | {b['data']} | {b['software']} | {b['boundary']} |")
    lines += [
        "",
        "## History / Dependency Edges",
        "",
        "| parent | child | reason |",
        "|---|---|---|",
    ]
    for parent, child, reason in EDGES:
        lines.append(f"| {parent} | {child} | {reason} |")
    lines += [
        "",
        "## Software Interpretation Ledger",
        "",
        "| method | software | inputs | outputs | interpretation | failure mode |",
        "|---|---|---|---|---|---|",
    ]
    for m in method_rows:
        lines.append(f"| {m['method']} | {m['software']} | {m['inputs']} | {m['outputs']} | {m['interpretation']} | {m['failure_mode']} |")
    lines += [
        "",
        "## Current Paper Count Decision",
        "",
        "Broad count: 11 active/near-active branches plus the Neoantigen separate repository branch.",
        "Immediate manuscript lane: Paper 1 trunk. Next high-value lane: Paper 11. Wet-lab priority lane: Paper 9 IC50/PRISM/DepMap. Visual-pilot lane: Paper 2/13. Data-gated lane: Paper 4 GD-HLA. Separate methods lane: CLEAN-Neo/CROSS-Neo.",
    ]
    return "\n".join(lines) + "\n"


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    pages = read_inventory()
    files = all_candidate_files()
    branch_rows: list[dict[str, object]] = []
    for branch in BRANCHES:
        row = dict(branch)
        row["source_files"] = matched_files(branch, files)
        row["hub_pages"] = page_matches(branch, pages)
        row["source_file_count"] = len(row["source_files"])
        row["hub_page_count"] = len(row["hub_pages"])
        branch_rows.append(row)

    branch_flat = []
    for row in branch_rows:
        flat = dict(row)
        flat["source_files"] = "; ".join(row.get("source_files", []))
        flat["hub_pages"] = "; ".join((p.get("page") or p.get("file") or p.get("path") or p.get("href") or "") for p in row.get("hub_pages", []))
        flat["parents"] = "; ".join(row.get("parents", []))
        flat["children"] = "; ".join(row.get("children", []))
        branch_flat.append(flat)

    write_tsv(
        OUT / "paper_repo_branches_v2.tsv",
        branch_flat,
        ["priority", "branch", "title", "status", "repo_role", "claim", "data", "software", "interpretation", "boundary", "parents", "children", "source_file_count", "source_files", "hub_page_count", "hub_pages"],
    )
    write_tsv(
        OUT / "paper_repo_edges_v2.tsv",
        [{"parent": a, "child": b, "reason": c} for a, b, c in EDGES],
        ["parent", "child", "reason"],
    )
    write_tsv(
        OUT / "software_method_ledger_v2.tsv",
        METHODS,
        ["method", "software", "inputs", "outputs", "interpretation", "failure_mode"],
    )
    write_tsv(
        OUT / "claim_boundary_ledger_v2.tsv",
        [{"branch": b["branch"], "status": b["status"], "claim": b["claim"], "boundary": b["boundary"], "interpretation": b["interpretation"]} for b in branch_rows],
        ["branch", "status", "claim", "boundary", "interpretation"],
    )
    (OUT / "paper_repo_branches_v2.json").write_text(json.dumps(branch_rows, ensure_ascii=False, indent=2), encoding="utf-8")
    (OUT / "PAPER_REPO_MAP_V2_KR.md").write_text(build_markdown(branch_rows, METHODS, len(pages), len(files)), encoding="utf-8")
    (OUT / "paper_repo_map_v2_2026_05_10.html").write_text(build_main(branch_rows, METHODS, len(pages), len(files)), encoding="utf-8")
    (OUT / "software_method_ledger_v2_2026_05_10.html").write_text(build_software_page(METHODS), encoding="utf-8")
    (OUT / "ic50_prediction_map_v2_2026_05_10.html").write_text(
        build_topic_page(
            "paper9-ic50-perturbation",
            "IC50 / PRISM / DepMap Prediction Map",
            "Drug-response and dependency proxy branch for wet-lab prioritization, not treatment proof.",
            branch_rows,
            METHODS,
        ),
        encoding="utf-8",
    )
    (OUT / "deconvolution_repo_v2_2026_05_10.html").write_text(
        build_topic_page(
            "paper1-main",
            "Deconvolution Repo Map",
            "Paper 1 mechanism-defense branch that links marker panels, cell-state estimates, MAPK/TDS and spatial support.",
            branch_rows,
            METHODS,
        ),
        encoding="utf-8",
    )
    (OUT / "neoantigen_repo_v2_2026_05_10.html").write_text(
        build_topic_page(
            "neoantigen-clean-neo",
            "Neoantigen / CLEAN-Neo Repo Map",
            "Separate methods repository for leakage-aware ranking, abstention, TCR/pMHC and MD-audited candidate selection.",
            branch_rows,
            METHODS,
        ),
        encoding="utf-8",
    )
    manifest_rows = [
        {"file": rel(p), "bytes": p.stat().st_size, "built_at": datetime.now().isoformat(timespec="seconds")}
        for p in sorted(OUT.iterdir())
        if p.is_file()
    ]
    write_tsv(OUT / "MANIFEST.tsv", manifest_rows, ["file", "bytes", "built_at"])
    print(json.dumps({"out": rel(OUT), "branches": len(branch_rows), "pages": len(pages), "source_files": len(files), "files_written": len(manifest_rows) + 1}, ensure_ascii=False))


if __name__ == "__main__":
    main()
