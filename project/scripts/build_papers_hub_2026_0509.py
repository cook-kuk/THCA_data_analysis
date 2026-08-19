#!/usr/bin/env python3
"""Build the 2026-05-09 renewed papers hub.

The builder is intentionally conservative:
- it never edits the 2026_05_04 source hub;
- it copies static assets/pages into project/papers_hub_2026_0509;
- it injects a v5 per-page data-definition block at the very top of every
  HTML page;
- it appends a machine-built figure/table/source audit annex so every visual
  has an explicit data/source row even when the original page caption was thin.
"""

from __future__ import annotations

import html
import json
import re
import shutil
import argparse
import csv
from dataclasses import dataclass, field
from datetime import datetime
from html.parser import HTMLParser
from pathlib import Path
from typing import Iterable


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "project" / "papers_hub_2026_05_04"
TARGET = ROOT / "project" / "papers_hub_2026_0509"
DEPLOY = Path("/var/www/papers/papers_hub_2026_0509")
BUILD_DATE = "2026-05-09"


DATA_TERMS = [
    {
        "key": "TCGA-THCA",
        "needles": ["tcga", "thca", "tert", "hm450", "cBioPortal", "braf", "ras"],
        "type": "cohort",
        "definition": "TCGA thyroid carcinoma cohort; RNA-seq, clinical, mutation, methylation and cBioPortal-derived mutation/clinical slices depending on page.",
        "use": "Paper 1 DM1/DM2 derivation, survival, promoter methylation, driver strata, and several downstream reuse pages.",
    },
    {
        "key": "TCGA-THCA HM450",
        "needles": ["hm450", "methylation", "promoter", "epigenetic", "fig8"],
        "type": "cohort",
        "definition": "Illumina 450K methylation subset used for thyroid-lineage promoter beta values.",
        "use": "Paper 1 Fig 8 and mechanism-layer pages.",
    },
    {
        "key": "K2 / PRJEB11591",
        "needles": ["k2", "korean", "prjeb11591", "nbnr"],
        "type": "cohort",
        "definition": "Korean PTC RNA-seq cohort used as an external Korean replication/calibration anchor.",
        "use": "Korean portability, mini-index calibration, and NBNR/DM-like comparisons.",
    },
    {
        "key": "GSE213647 / Lee 2024",
        "needles": ["gse213647", "lee 2024", "lee"],
        "type": "cohort",
        "definition": "Korean PTC RNA-seq cohort used for larger external replication and Hashimoto-like comparison.",
        "use": "Paper 1 cross-cohort generalization and HLA/deconvolution side pages.",
    },
    {
        "key": "GSE286332",
        "needles": ["gse286332", "ptc+ht", "hashimoto", "ht-isolated", "ptcht"],
        "type": "cohort",
        "definition": "Small Korean PTC +/- Hashimoto RNA-seq cohort used for HT-isolated signature discovery.",
        "use": "Paper 2 HT/HLA mechanism pages and TCGA transfer checks.",
    },
    {
        "key": "GSE250521",
        "needles": ["gse250521", "spatial", "visium", "lu 2023"],
        "type": "cohort",
        "definition": "Thyroid spatial transcriptomics / Visium resource spanning thyroid progression contexts.",
        "use": "Spatial validation, spatial failure/rescue, lag pockets, and image-DM1 spatial checks.",
    },
    {
        "key": "GSE76039 / Landa advanced thyroid",
        "needles": ["gse76039", "pdtc", "atc", "advanced", "landa"],
        "type": "cohort",
        "definition": "Advanced thyroid cancer expression resource used as a PDTC/ATC lineage-collapse benchmark.",
        "use": "Mechanism heterogeneity, advanced-disease comparisons, and reverse-causality boundary checks.",
    },
    {
        "key": "GSE184362 / GSE193581 / GSE241184",
        "needles": ["gse184362", "gse193581", "gse241184", "single-cell", "scrna", "thyrocyte"],
        "type": "cohort",
        "definition": "Single-cell thyroid resources used for thyrocyte-intrinsic and cross-author validation checks.",
        "use": "Paper 1 single-cell support, figure sourcing, and cell-state caveats.",
    },
    {
        "key": "Hugo/Riaz melanoma ICI pool",
        "needles": ["hugo", "riaz", "ici", "melanoma", "gse78220", "gse91061"],
        "type": "cohort",
        "definition": "External melanoma immunotherapy RNA/clinical resources used as ICI-response context, not thyroid predictor proof.",
        "use": "Paper 3 vulnerability-only evidence and claim-boundary pages.",
    },
    {
        "key": "AFND / Korean HLA references",
        "needles": ["afnd", "hla", "dqb1", "drb1", "dpb1", "graves", "gd"],
        "type": "reference",
        "definition": "Allele-frequency and literature-derived HLA reference material; some rows are registry/approval-gated rather than executable analysis inputs.",
        "use": "Paper 2 HLA and Paper 4 GD/Pan-Asian backlog pages.",
    },
    {
        "key": "8-gene mini-index",
        "needles": ["8-gene", "8gene", "mini-index", "dm1", "rai"],
        "type": "score",
        "definition": "Curated thyroid-lineage/dedifferentiation panel used to split molecular dark-matter states; higher DM1-like score means more dedifferentiated.",
        "use": "Paper 1 core axis and many downstream transfer/image/vulnerability pages.",
    },
    {
        "key": "RAI_8 / thyroid differentiation score",
        "needles": ["rai_8", "rai", "thyroid_diff", "tds", "differentiation"],
        "type": "score",
        "definition": "Mean z-score over thyroid differentiation genes; direction is differentiated / RAI-avid unless explicitly negated.",
        "use": "Lineage preservation/collapse comparisons and treatment-rationale pages.",
    },
    {
        "key": "MAPK output score",
        "needles": ["mapk", "mek", "braf", "ras", "ret", "ntrk"],
        "type": "score",
        "definition": "MAPK pathway output or driver-class proxy used to compare pathway activity against thyroid-lineage scores.",
        "use": "Paper 1 Fig 8 mechanism extensions v9-v18 and PRISM/DepMap overlays.",
    },
    {
        "key": "HLA-II / IFN-gamma module",
        "needles": ["hla-ii", "ifn", "ifng", "antigen", "autoimmunity"],
        "type": "score",
        "definition": "Inflammation / antigen-presentation module used for HT-overlap and ICI-vulnerability context.",
        "use": "Paper 2/Paper 3/Paper 4 HLA and immune-context pages.",
    },
    {
        "key": "Cohen's d / Spearman rho / Cox HR",
        "needles": ["cohen", "spearman", "rho", "cox", "hr", "forest", "auc", "or"],
        "type": "statistic",
        "definition": "Effect-size and association statistics reported from local TSV/JSON/HTML tables; direction must be read from each table caption.",
        "use": "Most dossier and validation pages.",
    },
    {
        "key": "PRISM / DepMap",
        "needles": ["prism", "depmap", "drug", "dependency", "synthetic"],
        "type": "external screen",
        "definition": "Cancer dependency/drug-screen resources used for vulnerability prioritization and MAPK-inhibitor overlays.",
        "use": "Paper 1/9/10/11 synthetic-lethality and druggability pages.",
    },
    {
        "key": "Neoantigen / pMHC / TCR resources",
        "needles": ["neoantigen", "pmhc", "tcr", "vdjdb", "mcpas", "esm", "kras", "vaccine"],
        "type": "external resource",
        "definition": "Neoantigen, protein-language-model, structural, TCR and literature-corroboration resources used in the vaccine pages.",
        "use": "Cancer vaccine, neoantigen, TCR and pMHC pages.",
    },
]


PAPER_GROUPS = [
    ("Paper 1", ["paper1", "deconv", "fig8", "spatial_lag", "synthetic_lethality"], "DM1 dark-matter / 8-gene / Fig 8 mechanism"),
    ("Paper 2", ["paper2", "image_dm1", "ht", "hla"], "H&E-DM1 and HT/HLA validation"),
    ("Paper 3", ["paper3", "ici"], "ICI vulnerability with predictor boundary"),
    ("Paper 4", ["paper4", "gd", "graves"], "Korean GD / Pan-Asian HLA backlog"),
    ("Paper 5-12", ["paper5", "paper6", "paper7", "paper9", "paper10", "paper11", "paper12"], "Atlas, network, synthetic-lethality and transfer pages"),
    ("Neoantigen/Vaccine", ["neo", "vaccine", "pmhc", "tcr", "hla_atlas"], "Neoantigen, pMHC, TCR and vaccine tooling"),
    ("Operations", ["master", "portfolio", "index", "situation", "plan", "audit", "report"], "Hub, portfolio, advisor and operational pages"),
]


@dataclass
class PageMeta:
    path: Path
    title: str = ""
    headings: list[tuple[str, str]] = field(default_factory=list)
    images: list[dict[str, str]] = field(default_factory=list)
    captions: list[str] = field(default_factory=list)
    table_captions: list[str] = field(default_factory=list)
    table_count: int = 0
    links: list[str] = field(default_factory=list)
    text_sample: str = ""


class MetaParser(HTMLParser):
    def __init__(self) -> None:
        super().__init__(convert_charrefs=True)
        self.title: list[str] = []
        self.headings: list[tuple[str, str]] = []
        self.images: list[dict[str, str]] = []
        self.captions: list[str] = []
        self.table_captions: list[str] = []
        self.links: list[str] = []
        self.table_count = 0
        self._capture: str | None = None
        self._buf: list[str] = []
        self._text: list[str] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        attrs_d = {k: (v or "") for k, v in attrs}
        if tag == "title":
            self._capture = "title"
            self._buf = []
        elif tag in {"h1", "h2", "h3"}:
            self._capture = tag
            self._buf = []
        elif tag in {"figcaption", "caption"}:
            self._capture = tag
            self._buf = []
        elif tag == "p" and "cap" in attrs_d.get("class", ""):
            self._capture = "figcaption"
            self._buf = []
        elif tag == "img":
            self.images.append(
                {
                    "src": attrs_d.get("src", ""),
                    "alt": attrs_d.get("alt", ""),
                    "class": attrs_d.get("class", ""),
                }
            )
        elif tag == "table":
            self.table_count += 1
        elif tag == "a" and attrs_d.get("href"):
            self.links.append(attrs_d["href"])

    def handle_data(self, data: str) -> None:
        if self._capture:
            self._buf.append(data)
        if len(self._text) < 1000:
            s = data.strip()
            if s:
                self._text.append(s)

    def handle_endtag(self, tag: str) -> None:
        if tag != self._capture:
            return
        text = normalize(" ".join(self._buf))
        if text:
            if tag == "title":
                self.title.append(text)
            elif tag in {"h1", "h2", "h3"}:
                self.headings.append((tag, text))
            elif tag == "figcaption":
                self.captions.append(text)
            elif tag == "caption":
                self.table_captions.append(text)
        self._capture = None
        self._buf = []

    @property
    def text_sample(self) -> str:
        return normalize(" ".join(self._text[:400]))


def normalize(value: str) -> str:
    return re.sub(r"\s+", " ", value or "").strip()


def strip_old_blocks(raw: str) -> str:
    raw = re.sub(
        r"<!-- DATA_DEFINITIONS_BLOCK_START.*?DATA_DEFINITIONS_BLOCK_END.*?-->\s*",
        "",
        raw,
        flags=re.S,
    )
    raw = re.sub(
        r"<!-- RENEWAL_0509_BLOCK_START.*?RENEWAL_0509_BLOCK_END.*?-->\s*",
        "",
        raw,
        flags=re.S,
    )
    raw = re.sub(
        r"<!-- RENEWAL_0509_ANNEX_START.*?RENEWAL_0509_ANNEX_END.*?-->\s*",
        "",
        raw,
        flags=re.S,
    )
    return raw


def parse_page(path: Path, raw: str) -> PageMeta:
    parser = MetaParser()
    parser.feed(raw)
    return PageMeta(
        path=path,
        title=parser.title[0] if parser.title else path.stem,
        headings=parser.headings,
        images=parser.images,
        captions=parser.captions,
        table_captions=parser.table_captions,
        table_count=parser.table_count,
        links=parser.links,
        text_sample=parser.text_sample,
    )


def infer_paper_group(meta: PageMeta) -> str:
    """Assign a page family from stable file names first.

    Page text often contains broad terms like HLA, HT, or paper links in the
    global navigation. Filename-first grouping avoids sweeping unrelated pages
    into Paper 2 just because a shared glossary appears in the body.
    """
    name = meta.path.name.lower()
    stem = meta.path.stem.lower()
    title = meta.title.lower()

    if name in {"index.html", "data_definitions_2026_0509.html", "_data_dictionary.html"}:
        return "Operations"
    if stem.endswith("_atlas_2026_0509") or stem == "family_index_2026_0509" or (stem.startswith("family_") and stem.endswith("_dossier_2026_0509")):
        return "Operations"
    if stem.startswith("master") or stem.startswith("portfolio_") or any(x in stem for x in ["situation", "advisor", "send_package", "execution_board", "upgrade"]):
        if "hla" in stem or "yu_send" in stem or "pillar" in stem or "k2_" in stem:
            return "HLA / Paper 2-4"
        return "Operations"
    if stem == "paper1" or stem.startswith("paper1_") or stem.startswith("portfolio_paper1") or stem.startswith("gpl570") or stem.startswith("intro_strategy"):
        return "Paper 1"
    if stem.startswith("paper2") or stem.startswith("portfolio_paper2") or "image_dm1" in stem or "ht_vs_gd" in stem or stem.startswith("pillar1"):
        return "Paper 2"
    if stem.startswith("paper3") or stem.startswith("portfolio_paper3") or "ici" in stem:
        return "Paper 3"
    if stem.startswith("paper4") or stem.startswith("portfolio_paper4") or stem.startswith("paperg") or stem.startswith("portfolio_paperg") or "gd_hla" in stem or "graves" in stem:
        return "Paper 4"
    if stem.startswith("paper5"):
        return "Paper 5"
    if stem.startswith("paper6") or stem.startswith("paperb") or stem.startswith("papern") or stem.startswith("portfolio_paperb") or stem.startswith("portfolio_papern") or "dial" in stem:
        return "Paper 6 / DIAL"
    if stem.startswith("paper7") or "agentic" in stem:
        return "Paper 7"
    if stem.startswith("paper9") or stem.startswith("portfolio_paper9"):
        return "Paper 9"
    if stem.startswith("paper10"):
        return "Paper 10"
    if stem.startswith("paper11") or "pancancer" in stem:
        return "Paper 11"
    if stem.startswith("paper12") or "network" in stem:
        return "Paper 12"
    if any(x in stem for x in ["neo", "vaccine", "tcr", "pmhc", "vdjdb", "alphafold", "lumenix", "paad", "egfr", "tp53", "anchor", "curator", "conflict", "public_epitopes", "rmsd", "esm2", "property", "batch30", "autoimmunity", "driver_gene"]):
        return "Neoantigen / Vaccine"
    if "spatial" in stem or "trajectory" in stem:
        return "Spatial / Trajectory"
    if "paper 1" in title:
        return "Paper 1"
    if "paper 2" in title:
        return "Paper 2"
    if "paper 3" in title:
        return "Paper 3"
    if "paper 4" in title:
        return "Paper 4"
    return "Other"


def infer_terms(meta: PageMeta) -> list[dict[str, str]]:
    hay = f"{meta.path.name} {meta.title} {meta.text_sample} " + " ".join(
        [img.get("src", "") + " " + img.get("alt", "") for img in meta.images[:80]]
    )
    hay_l = hay.lower()
    selected = []
    for term in DATA_TERMS:
        if any(n.lower() in hay_l for n in term["needles"]):
            selected.append(term)
    if not selected:
        selected = [DATA_TERMS[0], DATA_TERMS[-1]]
    return selected


def current_status_for_group(group: str) -> tuple[str, str, str]:
    """Current claim boundary for the 2026-05-09 marathon review surface."""
    rules = {
        "Paper 1": (
            "Active manuscript spine",
            "DM1 / 8-gene / Fig 8 mechanism pages are current only when read with the v13/v14 reviewer-defense layer; v15B-v18 spatial pages are caveat/autopsy support.",
            "Use for Paper 1 review first. Keep mechanism interpretation prose in author-protected manuscript sections.",
        ),
        "Paper 2": (
            "Bounded pilot / separated HLA track",
            "Image-DM1 is pilot-level until external H&E/pathologist QC passes. HT/PTC immune-expression pages are separate from GD HLA genetics.",
            "Use as evidence archive and validation queue, not as a flagship diagnostic claim.",
        ),
        "Paper 3": (
            "ICI vulnerability, not response predictor",
            "Current evidence supports readiness/prioritization framing only; it does not prove thyroid ICI response prediction.",
            "Keep Track A frozen unless raw thyroid ICI response data are explicitly added.",
        ),
        "Paper 4": (
            "Validation-gated GD / Pan-Asian HLA backlog",
            "HLA genetics claims require adult Korean GD germline validation before manuscript-level lock.",
            "Keep GD allele genetics out of Paper 2 cancer-cohort claims.",
        ),
        "Paper 5": (
            "Paused reserve",
            "RAI biomarker package remains useful but is paused to avoid same-data submission conflict with Paper 1/Yu priority.",
            "Reopen only after priority/submission routing is decided.",
        ),
        "Paper 6 / DIAL": (
            "Forensic correction required",
            "v52 LODO leakage finding blocks downstream use of earlier weak performance numbers.",
            "Treat old performance tables as audit history unless rebuilt under v52 rules.",
        ),
        "Paper 7": (
            "Framework / operations layer",
            "Useful for process discipline and agentic framework only; it does not authorize protected prose generation.",
            "Use as operating protocol during the marathon window.",
        ),
        "Paper 9": (
            "Hypothesis / prioritization layer",
            "Synthetic-lethality results are atlas candidates, not validated clinical actionability.",
            "Use to rank follow-up analyses; keep cell-line and cohort-size caveats attached.",
        ),
        "Paper 10": (
            "Hypothesis / prioritization layer",
            "Drug/network transfer pages are prioritization evidence only.",
            "Use as reserve support after Paper 1 decisions.",
        ),
        "Paper 11": (
            "Pan-cancer transfer / atlas layer",
            "Pan-cancer and DepMap/PRISM overlays support target ranking, not thyroid-specific therapeutic proof.",
            "Use for cross-cancer plausibility with explicit lineage-distance caveats.",
        ),
        "Paper 12": (
            "Network atlas layer",
            "Network pages are bulk/co-expression and dependency-prioritization support.",
            "Do not convert to causal mechanism language without validation.",
        ),
        "Neoantigen / Vaccine": (
            "Dossier ready, proof bounded",
            "Neoantigen/pMHC/TCR models support traceable ranking and abstention, not definitive immunogenicity.",
            "Frame PAAD vaccine-first only for resected/MRD low-burden contexts; routine THCA vaccine-first is not defensible.",
        ),
        "HLA / Paper 2-4": (
            "Separated HLA queue",
            "HT/PTC expression ecology and GD germline HLA genetics must remain separated.",
            "Route cancer expression evidence to Paper 2 and GD genetics to Paper 4.",
        ),
        "Spatial / Trajectory": (
            "Context / caveat support",
            "Spatial pages are valuable for anatomy and failure-mode explanation, but current GSE250521 evidence is not the primary Paper 1 lock.",
            "Use to explain heterogeneity and caveats, not as the main mechanism proof.",
        ),
        "Operations": (
            "Hub / review control plane",
            "These pages organize evidence and status; they do not add new biological claims.",
            "Use as navigation, reproducibility, and decision-state entry points.",
        ),
        "Other": (
            "Archive / evidence page",
            "No specific current paper family was inferred from filename/title; read this page as supporting archive unless promoted manually.",
            "Check the figure/table annex and source links before reusing claims.",
        ),
    }
    return rules.get(group, rules["Other"])


def result_inventory() -> dict[str, list[str]]:
    idx: dict[str, list[str]] = {}
    roots = [
        ROOT / "project" / "results" / "p_deconv_2026_05_08",
        ROOT / "project" / "results" / "p2_image_dm1_v2_foundation_clam_2026_05_07",
        ROOT / "project" / "results" / "p_cancer_vaccine_strategy_2026_05_09",
        ROOT / "project" / "results" / "p_neo_bayesian_2026_05_09",
        ROOT / "project" / "results" / "p_synthetic_lethality_2026_05_09",
        ROOT / "project" / "results" / "paper11_pancancer",
        ROOT / "project" / "results" / "hla_two_paper_synthesis_2026_05_09",
    ]
    suffixes = {".png", ".jpg", ".jpeg", ".pdf", ".svg", ".tsv", ".csv", ".json", ".md", ".html"}
    for root in roots:
        if not root.exists():
            continue
        for f in root.rglob("*"):
            if f.is_file() and f.suffix.lower() in suffixes:
                idx.setdefault(f.name.lower(), []).append(str(f.relative_to(ROOT)))
    return idx


def asset_matches(src: str, idx: dict[str, list[str]]) -> list[str]:
    name = Path(src.split("#", 1)[0].split("?", 1)[0]).name.lower()
    if not name:
        return []
    matches = idx.get(name, [])[:6]
    stem = Path(name).stem
    if not matches and stem:
        for key, paths in idx.items():
            if stem and stem in key:
                matches.extend(paths[:3])
                if len(matches) >= 6:
                    break
    return matches[:6]


def guess_visual_data(meta: PageMeta, img: dict[str, str]) -> str:
    text = f"{meta.path.name} {meta.title} {img.get('src','')} {img.get('alt','')}".lower()
    rules = [
        ("paper1" in text or "deconv" in text or "fig_sx" in text, "Paper 1 THCA lineage axis: TCGA/K2/Lee/GSE76039/spatial/DepMap inputs as named by asset and caption."),
        ("paper2_hla" in text or "hla_ncomm" in text, "HLA allele-frequency/literature extraction tables and external GEO validation figures."),
        ("image_dm1" in text or "clam" in text or "pathology" in text, "H&E/image-DM1 audit outputs from TCGA slides, pathology/spatial scout tables, and CLAM/encoder experiments."),
        ("paper3_ici" in text or "ici" in text, "Paper 3 ICI vulnerability resources: thyroid-side RAI/HLA axes plus external ICI cohorts where explicitly named."),
        ("paper4_hla" in text or "graves" in text or "gd" in text, "Paper 4 GD/Pan-Asian HLA registry and readiness matrices."),
        ("nature_master" in text or "paper11" in text or "pancancer" in text, "Paper 11 pan-cancer transfer / DepMap / PRISM / hallmark outputs."),
        ("spatial_full" in text or "visium" in text, "Spatial transcriptomics overlays and spot-level module summaries."),
        ("neo" in text or "vaccine" in text or "pmhc" in text or "tcr" in text, "Neoantigen/pMHC/TCR/vaccine benchmark resources and derived model outputs."),
    ]
    for cond, label in rules:
        if cond:
            return label
    return "Inline page visual; source path and caption are the binding provenance row for this audit."


def visual_reading_notes(meta: PageMeta, img: dict[str, str], caption: str) -> tuple[str, str, str]:
    text = f"{meta.path.name} {meta.title} {img.get('src','')} {img.get('alt','')} {caption}".lower()
    if any(tok in text for tok in ["forest", "hr", "or", "cox", "meta"]):
        return (
            "Effect-size / forest evidence",
            "Read the point estimate, confidence interval, cohort label, and direction together before reusing the claim.",
            "A narrow CI supports precision, but cohort definitions and endpoint choice still set the claim boundary.",
        )
    if any(tok in text for tok in ["heatmap", "matrix", "cluster", "clustering"]):
        return (
            "Matrix / clustering evidence",
            "Read the row/column annotations and color scale first; the visual is about structure, not a single p-value.",
            "Cluster appearance should be paired with the linked table or statistic before manuscript-level use.",
        )
    if any(tok in text for tok in ["spatial", "visium", "spot", "pocket", "lag"]):
        return (
            "Spatial / anatomy evidence",
            "Read the tissue region, spot-level module score, and neighborhood definition together.",
            "Spatial pages are current as context/caveat support unless explicitly promoted by the current situation matrix.",
        )
    if any(tok in text for tok in ["umap", "tsne", "embedding", "trajectory"]):
        return (
            "Embedding / trajectory evidence",
            "Use the embedding to inspect separation, continuity, and outlier structure rather than absolute distance.",
            "Embedding geometry is qualitative; use the companion table for quantitative claims.",
        )
    if any(tok in text for tok in ["survival", "kaplan", "km", "os", "dfs", "pfs"]):
        return (
            "Clinical endpoint evidence",
            "Check endpoint, risk group definition, censoring context, and whether the page reports Cox or log-rank support.",
            "Low-event settings should remain caveated even when the visual separation is strong.",
        )
    if any(tok in text for tok in ["bar", "box", "violin", "scatter", "volcano", "correlation", "rho"]):
        return (
            "Association / distribution evidence",
            "Read group labels, axis direction, effect size, and whether the plot is residualized or unadjusted.",
            "Association plots do not establish mechanism without the linked validation layer.",
        )
    if any(tok in text for tok in ["hla", "allele", "afnd", "dqb", "drb", "dpb"]):
        return (
            "HLA / allele-frequency evidence",
            "Read allele resolution, ancestry panel, registry/literature source, and validation gate before using the result.",
            "Cancer expression ecology and germline GD genetics must remain separated.",
        )
    if any(tok in text for tok in ["neo", "vaccine", "pmhc", "tcr", "esm", "affinity"]):
        return (
            "Neoantigen / pMHC ranking evidence",
            "Read candidate identity, HLA context, model score, and orthogonal support as a ranking stack.",
            "Model ranking is not definitive immunogenicity proof without assay or clinical validation.",
        )
    return (
        "Page visual evidence",
        "Read the asset path, caption, and surrounding section heading as the binding context for the figure.",
        "If no upstream TSV/JSON match is listed, treat the embedded page asset as the immediate provenance row.",
    )


def build_flow_cards(meta: PageMeta, limit: int = 8) -> str:
    headings = [(tag, normalize(text)) for tag, text in meta.headings if normalize(text)]
    if not headings:
        return '<article class="renewal0509-step"><b>01</b><span>No h1-h3 headings detected; use the table and figure audit below as the reading order.</span></article>'
    cards = []
    for i, (tag, text) in enumerate(headings[:limit], 1):
        cards.append(
            '<article class="renewal0509-step">'
            f"<b>{i:02d}</b><em>{html.escape(tag.upper())}</em><span>{html.escape(text[:220])}</span>"
            "</article>"
        )
    if len(headings) > limit:
        cards.append(
            '<article class="renewal0509-step renewal0509-step-more">'
            f"<b>+{len(headings) - limit}</b><span>More section headings are mapped in the bottom annex.</span>"
            "</article>"
        )
    return "\n".join(cards)


def build_figure_cards(meta: PageMeta, idx: dict[str, list[str]] | None = None, limit: int | None = None) -> str:
    images = meta.images if limit is None else meta.images[:limit]
    if not images:
        return '<article class="renewal0509-figcard"><div><b>No figure assets detected</b><p>This page has no image tags. Use the table cards and source links below.</p></div></article>'
    cards = []
    for i, img in enumerate(images, 1):
        caption = meta.captions[i - 1] if i - 1 < len(meta.captions) else img.get("alt", "")
        caption = normalize(caption) or "No caption/alt detected in the source HTML."
        src = img.get("src", "")
        role, how, caveat = visual_reading_notes(meta, img, caption)
        matches = asset_matches(src, idx or {}) if idx is not None else []
        match_text = "; ".join(matches) if matches else "No same-basename result match listed."
        src_attr = html.escape(src, quote=True)
        alt_attr = html.escape(caption[:180], quote=True)
        cards.append(
            '<article class="renewal0509-figcard">'
            f'<div class="renewal0509-figframe"><img src="{src_attr}" alt="{alt_attr}" loading="lazy"></div>'
            '<div class="renewal0509-figbody">'
            f'<div class="renewal0509-mini">Figure {i:02d} · {html.escape(role)}</div>'
            f"<h4>{html.escape(Path(src).name or 'inline visual')}</h4>"
            f"<p>{html.escape(caption[:520])}</p>"
            '<dl>'
            f"<dt>Data entered</dt><dd>{html.escape(guess_visual_data(meta, img))}</dd>"
            f"<dt>How to read</dt><dd>{html.escape(how)}</dd>"
            f"<dt>Boundary</dt><dd>{html.escape(caveat)}</dd>"
            f"<dt>Source match</dt><dd>{html.escape(match_text)}</dd>"
            "</dl></div></article>"
        )
    return "\n".join(cards)


def build_table_cards(meta: PageMeta, limit: int | None = None) -> str:
    count = meta.table_count if limit is None else min(meta.table_count, limit)
    if count == 0:
        return '<article class="renewal0509-tablecard"><b>No HTML tables detected</b><span>The page evidence is figure/text driven; source links are still audited below.</span></article>'
    cards = []
    for i in range(1, count + 1):
        caption = meta.table_captions[i - 1] if i - 1 < len(meta.table_captions) else f"Inline HTML table {i}"
        cards.append(
            '<article class="renewal0509-tablecard">'
            f"<b>Table {i:02d}</b>"
            f"<strong>{html.escape(normalize(caption)[:260])}</strong>"
            "<span>Rows are embedded in this HTML file. Treat the table body as the immediate data source; upstream TSV/JSON is listed when linked or named in the page.</span>"
            "</article>"
        )
    if meta.table_count > count:
        cards.append(
            '<article class="renewal0509-tablecard renewal0509-step-more">'
            f"<b>+{meta.table_count - count}</b><strong>Additional tables</strong><span>Full table inventory appears in the bottom annex.</span>"
            "</article>"
        )
    return "\n".join(cards)


METRIC_PATTERN = re.compile(
    r"\b(?:HR|OR|AUC|AUROC|rho|r|Cohen'?s d|d|p|FDR|q|I2|I²|n)\s*[=:<>]\s*[-+~]?\d+(?:\.\d+)?(?:e[-+]?\d+)?%?"
    r"|[-+]?\d+(?:\.\d+)?\s*%"
    r"|\b\d+\s*/\s*\d+\b"
    r"|\b\d+(?:\.\d+)?\s*\[[^\]]{3,60}\]",
    re.I,
)


GENERATED_PAGE_NAMES = {
    "index.html",
    "data_definitions_2026_0509.html",
    "current_situation_2026_0509.html",
    "figure_atlas_2026_0509.html",
    "table_atlas_2026_0509.html",
    "flow_atlas_2026_0509.html",
    "metric_atlas_2026_0509.html",
    "source_atlas_2026_0509.html",
    "reviewer_tour_2026_0509.html",
    "top22_review_pack_2026_0509.html",
    "hub_search_2026_0509.html",
    "evidence_qc_dashboard_2026_0509.html",
    "claim_boundary_ledger_2026_0509.html",
    "paper1_critical_path_2026_0509.html",
    "figure_data_matrix_2026_0509.html",
    "paper_repo_control_room_2026_05_10.html",
    "high_impact_board_2026_05_10.html",
    "impact_escalation_matrix_2026_05_10.html",
    "killer_figure_storyboard_2026_05_10.html",
    "project_inventory_2026_05_10.html",
    "paper_repo_map_2026_05_10.html",
    "paper_branch_manifest_2026_05_10.html",
    "paper_date_topic_index_2026_05_10.html",
    "software_method_ledger_2026_05_10.html",
    "ic50_prediction_map_2026_05_10.html",
    "deconvolution_repo_2026_05_10.html",
    "neoantigen_repo_2026_05_10.html",
    "family_index_2026_0509.html",
}


def is_generated_hub_page(meta: PageMeta) -> bool:
    name = meta.path.name
    return (
        name in GENERATED_PAGE_NAMES
        or name.startswith("family_") and name.endswith("_dossier_2026_0509.html")
        or name.startswith("branch_dossier_") and name.endswith("_2026_05_10.html")
    )


def audit_pages(pages: list[PageMeta]) -> list[PageMeta]:
    return [p for p in pages if not is_generated_hub_page(p)]


def extract_metric_tokens(meta: PageMeta, limit: int = 14) -> list[str]:
    hay = " ".join(
        [meta.title, meta.text_sample]
        + [caption for caption in meta.captions[:20]]
        + [caption for caption in meta.table_captions[:20]]
    )
    seen: set[str] = set()
    out: list[str] = []
    for match in METRIC_PATTERN.finditer(hay):
        token = normalize(match.group(0))
        key = token.lower()
        if len(token) < 2 or key in seen:
            continue
        seen.add(key)
        out.append(token)
        if len(out) >= limit:
            break
    return out


def build_metric_strip(meta: PageMeta, limit: int = 12) -> str:
    metrics = extract_metric_tokens(meta, limit=limit)
    if not metrics:
        return '<span class="renewal0509-metric muted">No numeric/statistical token extracted from visible page text.</span>'
    return " ".join(f'<span class="renewal0509-metric">{html.escape(m)}</span>' for m in metrics)


def slugify_family(value: str) -> str:
    slug = re.sub(r"[^a-z0-9]+", "_", value.lower()).strip("_")
    return slug or "other"


def review_priority_score(meta: PageMeta) -> int:
    group = infer_paper_group(meta)
    boosts = {
        "Paper 1": 120,
        "Paper 2": 70,
        "Paper 3": 55,
        "Paper 4": 45,
        "Paper 6 / DIAL": 35,
        "Paper 11": 30,
        "Neoantigen / Vaccine": 25,
    }
    name = meta.path.stem.lower()
    name_boost = 0
    for token, boost in [
        ("paper1", 40),
        ("fig8", 35),
        ("reviewer", 35),
        ("defense", 32),
        ("methods", 28),
        ("dossier", 24),
        ("audit", 18),
        ("current", 16),
        ("status", 16),
        ("nature", 14),
        ("spatial", 10),
    ]:
        if token in name:
            name_boost += boost
    return (
        boosts.get(group, 10)
        + name_boost
        + len(meta.images) * 6
        + meta.table_count * 4
        + len(meta.headings) * 2
        + len(extract_metric_tokens(meta, limit=22)) * 5
    )


def source_chips_for_page(meta: PageMeta, idx: dict[str, list[str]], limit: int = 18) -> list[str]:
    out: list[str] = []
    seen: set[str] = set()
    for href in meta.links:
        if href.startswith("#") or href.startswith("javascript:"):
            continue
        if any(tok in href.lower() for tok in ["result", "asset", ".tsv", ".csv", ".json", ".md", ".png", ".pdf", ".html"]):
            if href not in seen:
                seen.add(href)
                out.append(href)
    for img in meta.images:
        for match in asset_matches(img.get("src", ""), idx):
            if match not in seen:
                seen.add(match)
                out.append(match)
            if len(out) >= limit:
                return out
    return out[:limit]


def quality_items(meta: PageMeta, idx: dict[str, list[str]]) -> dict[str, bool]:
    return {
        "figures": len(meta.images) > 0,
        "tables": meta.table_count > 0,
        "flow": len(meta.headings) >= 3,
        "metrics": len(extract_metric_tokens(meta, limit=4)) > 0,
        "sources": len(source_chips_for_page(meta, idx, limit=4)) > 0,
        "definitions": len(infer_terms(meta)) > 0,
    }


def quality_score(meta: PageMeta, idx: dict[str, list[str]]) -> tuple[int, str]:
    items = quality_items(meta, idx)
    score = round(sum(items.values()) / len(items) * 100)
    grade = "A" if score >= 84 else "B" if score >= 67 else "C" if score >= 50 else "D"
    return score, grade


PAPER_BRANCHES = [
    {
        "branch": "paper1-main",
        "paper": "Paper 1",
        "title": "DM1 molecular dark matter / 8-gene thyroid-lineage axis",
        "state": "SUBMIT_NOW",
        "parent": "root",
        "history": "v2/v3/v5 deconvolution -> v13 TDS16/MAPK cherry-pick lock -> v14 cross-cohort forest -> v18 spatial caveat layer",
        "data": "TCGA-THCA RNA/HM450/clinical/SV, K2/PRJEB11591, Lee/GSE213647, GSE76039, Lu 2023/GSE193581, Pu 2021 single-cell, DepMap/PRISM reserve",
        "software": "Python, pandas, scipy, statsmodels/sklearn, lifelines-style survival tables, scipy NNLS, sklearn NuSVR, matplotlib/seaborn, HTML dossier builder",
        "claim": "8-gene DM1/DM2 axis defines fusion-enriched, thyroid-lineage-silenced molecular dark matter with reviewer-defense evidence.",
        "meaning": "This is the trunk branch. Other THCA papers either support it, fork from it, or stay parked until Paper 1 is submitted.",
        "boundary": "Mechanism is correlative and motivates rather than confirms causality; protected manuscript prose remains author-keyboard only.",
        "paths": "project/results/p_deconv_2026_05_08/SUMMARY.md; project/HANDOFF_2026_05_09.md; project/papers_hub_2026_0509/paper1.html",
    },
    {
        "branch": "paper1b-braf-axis",
        "paper": "Paper 1B",
        "title": "BRAF-cPTC DM2 / DM1xTERT escape branch",
        "state": "MERGE_OR_STANDALONE_DECISION",
        "parent": "paper1-main",
        "history": "BRAF sprint -> takeover dossier -> survival/TERT/RAI axis audit",
        "data": "TCGA-THCA BRAF/RAS/TERT/clinical, Landa/GSE76039 advanced thyroid, Paper 1 DM labels, image rescue negative controls",
        "software": "Python survival/effect-size scripts, pandas, scipy/statsmodels, HTML dossier builder",
        "claim": "BRAF-cPTC risk is not generic DM1; DM2 and a DM1xTERT escape state carry the adverse arm.",
        "meaning": "High-value fork. It can strengthen Paper 1 or become a standalone follow-up if it does not cannibalize Paper 1.",
        "boundary": "Do not state that generic DM1 is aggressive in BRAF-cPTC; the aggressive subcase is DM1xTERT/escape framing.",
        "paths": "project/results/p2_braf_nature_sprint_2026_05_09/; project/papers_hub_2026_0509/paper1_paper2_braf_axis_takeover_dossier.html",
    },
    {
        "branch": "paper2-wsi-image",
        "paper": "Paper 2",
        "title": "H&E / WSI projection of DM1",
        "state": "VIABLE_PILOT",
        "parent": "paper1-main",
        "history": "foundation CLAM audit -> RunPod WSI recovery -> 50-slide DINOv2 pilot",
        "data": "TCGA-THCA WSI balanced manifest, 50 embedded slides, DM1/DM2 labels from Paper 1, RunPod split artifacts",
        "software": "OpenSlide-style tiling, DINOv2 ViT-L fallback embeddings, CLAM/foundation-model pipeline, sklearn LOSO models, RunPod A6000 execution",
        "claim": "H&E contains a moderate image-DM1 signal in a small corrected pilot.",
        "meaning": "Useful digital-pathology fork, but it is not the trunk. It is a pilot/methods/audit paper unless external H&E validation arrives.",
        "boundary": "No clinical deployment claim; no robust BRAF/cPTC image-biology claim; small N and center/TSS confounding remain.",
        "paths": "project/results/p_wsi_dm1_2026_05_09/SUMMARY.md; project/papers_hub_2026_0509/wsi_dm1_recovery_dossier_2026_05_09.html",
    },
    {
        "branch": "paper2b-ht-immune",
        "paper": "Paper 2B",
        "title": "HT-overlap PTC immune/HLA expression context",
        "state": "BIOLOGY_STRONG_ALLELE_WEAK",
        "parent": "paper1-main",
        "history": "HT signature transfer -> HLA two-paper synthesis -> spatial/single-cell validation manifests",
        "data": "TCGA HT signature transfer, GSE286332 PTC+HT, Lu 2023 single-cell, GSE250521 spatial, Lee/GSE213647, HLA expression modules",
        "software": "Python module scoring, ssGSEA/singscore-style scoring, Fisher/OR tests, spatial score summaries, HTML execution board",
        "claim": "HT-overlap/immune-active PTC has HLA-II/IFN/TLS/APC biology that can be tracked across expression and spatial resources.",
        "meaning": "Biology branch is real; allele-genetics branch is not ready.",
        "boundary": "Do not claim HLA allele association from cancer cohorts or n=9 vs n=9 GSE286332 allele tests.",
        "paths": "project/results/paper2_ht_biology/; project/results/hla_two_paper_synthesis_2026_05_09/; project/papers_hub_2026_0509/paper2_hla.html",
    },
    {
        "branch": "paper3-ici-readiness",
        "paper": "Paper 3",
        "title": "ICI vulnerability/readiness in advanced thyroid context",
        "state": "FROZEN_DESIGN",
        "parent": "paper1-main",
        "history": "Track A registry -> Track B-lite public-cohort module scoring -> frozen response-prediction boundary",
        "data": "Hugo/Riaz melanoma ICI cohorts, GSE76039/GSE33630/GSE60542/GSE126698 thyroid cohorts, GSE151179 RAI before/after, ATC cell-line bulk",
        "software": "Python module scoring, PCA/hclust, meta-analysis effect-size tables, paired tests, sign-consistency tables",
        "claim": "The current data support immune vulnerability/readiness patterns, not thyroid ICI response prediction.",
        "meaning": "Future branch or Paper 11 reserve. It should not slow Paper 1.",
        "boundary": "No thyroid ICI response predictor claim without thyroid ICI-treated response data.",
        "paths": "project/results/paper3_ici_track_b_lite/phase2_summary.json; project/results/paper3_ici_track_b_lite/phase4_summary.json; project/papers_hub_2026_0509/paper3_ici_status_dossier_2026_05_09.html",
    },
    {
        "branch": "paper4-gd-hla",
        "paper": "Paper 4",
        "title": "Korean GD / Pan-Asian HLA architecture",
        "state": "UPGRADE_IF_DATA_ARRIVES",
        "parent": "paper2b-ht-immune",
        "history": "Pan-Asian literature extraction -> Chu/Shin/Chen anchors -> Korean replication readiness ledger",
        "data": "Chu 2018, Shin 2019, Chen 2011 HLA allele rows, AFND/literature references, candidate Korean replication cohort plan",
        "software": "Python meta-analysis, random-effects screening, BH-FDR tables, source extractability ledger",
        "claim": "Pan-Asian GD HLA architecture is coherent around DPB1*05:01/B*46:01/A*02:07/C*01:02 signals.",
        "meaning": "High-upside branch if Korean adult GD germline NGS arrives.",
        "boundary": "Current result is not a definitive Korean GD discovery; de novo matched case-control data are needed.",
        "paths": "project/results/paper4_gd_hla/; project/results/hla_two_paper_synthesis_2026_05_09/; project/papers_hub_2026_0509/hla_ncomm_upgrade_2026_05_09.html",
    },
    {
        "branch": "paper9-ic50-perturbation",
        "paper": "Paper 9",
        "title": "IC50 / PRISM / synthetic-lethality prioritization map",
        "state": "WETLAB_REQUIRED",
        "parent": "paper1-main",
        "history": "Paper 9 perturbation evidence -> synthetic lethality summary -> PRISM/DepMap actionability reserve",
        "data": "DepMap dependency, PRISM drug-response LFC/AUC proxies, MAPK-inhibitor hits, SLC1A5/GLUD1/GLS metabolic candidates",
        "software": "Python correlation/effect-size scripts, FDR enrichment, Spearman rho, Cohen's d, PRISM/DepMap matrix parsing",
        "claim": "DM1-high/dedifferentiated state prioritizes MAPK-inhibitor sensitivity and glutamine/NAMPT/MYC candidate vulnerabilities.",
        "meaning": "This is the IC50/drug-response map branch. It is useful as actionability/reviewer reserve.",
        "boundary": "PRISM LFC is a drug-response proxy, not direct clinical IC50. Do not claim validated synthetic lethality without perturbation/rescue wet-lab.",
        "paths": "project/results/p_synthetic_lethality_2026_05_09/synthetic_lethality_summary.json; project/results/paper9_perturbation_2026_05_06/paper9_perturbation_evidence_summary.md; project/papers_hub_2026_0509/paper1_synthetic_lethality_dossier.html",
    },
    {
        "branch": "paper11-pancancer",
        "paper": "Paper 11",
        "title": "Pan-cancer portable DM1 axis",
        "state": "NEXT_MAJOR_AFTER_P1",
        "parent": "paper1-main",
        "history": "Pan-cancer summary -> phase A-H modules -> Nature master dossier",
        "data": "TCGA pan-cancer 11k samples across 32-33 lineages, DepMap, PRISM, hallmark/pathway overlays, survival/ICI reserve",
        "software": "Python matrix scoring, lineage-wise correlations, Cox/survival tables, DepMap/PRISM overlays, matplotlib figures",
        "claim": "The thyroid-derived DM1/dedifferentiation architecture ports to multiple cancer lineages and links to dependency/drug-response axes.",
        "meaning": "Best major branch after Paper 1; Paper 10/12 should support this rather than fragment into separate papers.",
        "boundary": "Do not overstate full methylation causality without full HM450/wet-lab layer.",
        "paths": "project/results/paper11_pancancer/summary.json; project/results/paper11_pancancer/nature_master_summary_v2.json; project/papers_hub_2026_0509/paper11_pancancer.html",
    },
    {
        "branch": "paper10-12-network-atlas",
        "paper": "Paper 10/12",
        "title": "Atlas / network companion",
        "state": "ABSORB_INTO_P11",
        "parent": "paper11-pancancer",
        "history": "cell-line atlas -> network architecture -> candidate module support",
        "data": "THCA/cell-line atlas, network edges/modules, 35 candidates, 7 modules, 137 top edges",
        "software": "Python network construction, module detection, edge ranking, atlas HTML generation",
        "claim": "Network and atlas modules explain candidate architecture around LYN/NAMPT/TACSTD2/KCNN4/PAX8.",
        "meaning": "Support branch for Paper 11, not a separate near-term manuscript.",
        "boundary": "Bulk/network association is not causal mechanism by itself.",
        "paths": "project/results/paper10_atlas/; project/results/paper12_network/summary.json; project/papers_hub_2026_0509/paper12_network.html",
    },
    {
        "branch": "neoantigen-clean-neo",
        "paper": "Neo/CROSS-Neo",
        "title": "Leakage-aware neoantigen prediction / vaccine platform",
        "state": "SEPARATE_UNIVERSE",
        "parent": "root",
        "history": "CROSS-Neo lockdown -> CLEAN-Neo blueprint -> BAR-Neo abstention/interpretability -> PAAD/THCA vaccine triage",
        "data": "ITSNdb, TESLA/CEDAR/IEDB-like public corpora audits, HLA class I benchmark rows, PAAD/THCA vaccine triage evidence, patient-level metadata requirements",
        "software": "Python benchmark harness, scikit-learn baselines, structure/ESMFold proxy features, leakage/public-overlap audits, BMA/abstention queues",
        "claim": "Current platform supports leakage-aware ranking, abstention, and patient triage; it does not prove definitive immunogenicity.",
        "meaning": "Separate methods/resource universe. Keep apart from THCA Paper 1-4 except as future platform context.",
        "boundary": "Do not imply thyroid vaccine readiness from PAAD/PDAC evidence; source-heldout performance remains heterogeneous.",
        "paths": "project/results/p_neo_bayesian_2026_05_09/NEOANTIGEN_METHOD_BLUEPRINT_2026_05_09.md; project/results/cross_neo_v1_lockdown/NEOANTIGEN_ALL_ALGORITHMS_ALL_TESTSETS_FAIR_COMPARISON_SUMMARY.md; project/results/clean_neobench_barneo_2026_05_09/",
    },
]


# Expanded 2026-05-10 registry: 14 numbered paper slots plus one separate methods repo.
# This supersedes the earlier strategic-only branch list above. The older list is
# kept in-file as historical context, but every repo-map page uses this registry.
PAPER_BRANCHES = [
    {
        "slot_type": "paper",
        "branch": "paper1-main",
        "paper": "Paper 1",
        "title": "DM1 molecular dark matter / 8-gene thyroid-lineage axis",
        "topic": "DM1 trunk / thyroid dedifferentiation",
        "opened": "2026-04-29",
        "last_update": "2026-05-10",
        "state": "SUBMIT_NOW",
        "parent": "root",
        "history": "v17/v18 Paper 1 trunk; Fig 8 v6-v14; v15-v18 deconvolution/spatial caveat; 0509 hub renewal.",
        "data": "TCGA-THCA RNA/HM450/clinical/SV, K2/PRJEB11591, Lee/GSE213647, GSE76039, Lu 2023/GSE250521, Pu 2021 single-cell, DepMap/PRISM reserve.",
        "software": "Python, pandas, scipy, statsmodels, scikit-learn, NuSVR, survival/meta-analysis scripts, matplotlib, HTML dossier builder.",
        "claim": "8-gene DM1/DM2 axis defines fusion-enriched, thyroid-lineage-silenced molecular dark matter with reviewer-defense evidence.",
        "meaning": "Trunk branch. Every other THCA branch either supports it, forks from it, or stays parked until Paper 1 is submitted.",
        "boundary": "Mechanism is correlative and motivates rather than confirms causality; protected manuscript prose remains author-keyboard only.",
        "paths": "project/results/p_deconv_2026_05_08/SUMMARY.md; project/HANDOFF_2026_05_09.md; project/papers_hub_2026_0509/paper1.html",
        "tokens": "paper1 fig8 deconv gpl570 reviewer_defense nature_cancer dm1 tds k2",
    },
    {
        "slot_type": "paper",
        "branch": "paper1b-braf-axis",
        "paper": "Paper 1B",
        "title": "BRAF-cPTC DM2 / DM1xTERT escape branch",
        "topic": "BRAF / TERT / MAPK subtype branch",
        "opened": "2026-05-09",
        "last_update": "2026-05-10",
        "state": "MERGE_OR_STANDALONE_DECISION",
        "parent": "paper1-main",
        "history": "BRAF Nature sprint -> takeover dossier -> survival/TERT/RAI/HM450/MAPK axis audit.",
        "data": "TCGA BRAF/RAS/TERT/clinical, Lee 2024, HM450, GSE76039, spatial BRAF inference, Paper 1 DM labels.",
        "software": "Python survival/effect-size scripts, Fisher tests, pandas/scipy/statsmodels, deconvolution, spatial scoring, HTML synthesis.",
        "claim": "BRAF-classic PTC is not generic DM1; DM2 and a DM1xTERT escape state carry the adverse arm.",
        "meaning": "High-value fork. It can strengthen Paper 1 or become standalone if it does not cannibalize the trunk.",
        "boundary": "Do not state generic DM1 is aggressive in BRAF-cPTC; the adverse subcase is DM1xTERT/escape framing.",
        "paths": "project/results/p2_braf_nature_sprint_2026_05_09/; project/papers_hub_2026_0509/paper1_paper2_braf_axis_takeover_dossier.html",
        "tokens": "braf tert paper1_paper2_braf dm2 mapk hm450 rai_response",
    },
    {
        "slot_type": "paper",
        "branch": "paper2-wsi-image",
        "paper": "Paper 2",
        "title": "H&E / WSI projection of DM1",
        "topic": "Digital pathology / image-DM1",
        "opened": "2026-05-04",
        "last_update": "2026-05-10",
        "state": "VIABLE_PILOT",
        "parent": "paper1-main",
        "history": "H&E->DM1 audit -> foundation CLAM/UNI/DINOv2 tests -> RunPod WSI recovery -> Path2Space bridge.",
        "data": "TCGA-THCA WSI slides, 50 embedded pilot slides, GSE230424/GSE250521/GSE248205 spatial/pathology resources, Paper 1 DM labels.",
        "software": "OpenSlide-style tiling, DINOv2/UNI embeddings, CLAM, sklearn LOSO models, RunPod A6000 jobs, spatial bridge scripts.",
        "claim": "H&E and spatial-pathology layers contain a moderate image-DM1 / thyroid-state signal in pilot settings.",
        "meaning": "Useful digital-pathology fork; pilot/methods paper unless external H&E validation or K2 slides arrive.",
        "boundary": "No clinical deployment claim; no leakage-prone slide/site inflation; small N and center/TSS confounding remain.",
        "paths": "project/results/p_wsi_dm1_2026_05_09/SUMMARY.md; project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/; project/papers_hub_2026_0509/paper2_image_dm1.html",
        "tokens": "wsi image_dm1 clam pathology path2space gse230424 gse250521 cv2",
    },
    {
        "slot_type": "paper",
        "branch": "paper2b-ht-immune",
        "paper": "Paper 2B",
        "title": "HT-overlap PTC immune/HLA expression context",
        "topic": "Hashimoto / HLA / immune context",
        "opened": "2026-05-04",
        "last_update": "2026-05-10",
        "state": "BIOLOGY_STRONG_ALLELE_WEAK",
        "parent": "paper2-wsi-image",
        "history": "HT signature transfer -> HLA two-paper synthesis -> spatial/single-cell immune-context validation.",
        "data": "TCGA HT signature transfer, GSE286332 PTC+HT, GSE248205, GSE250521, Lee/GSE213647, HLA expression modules.",
        "software": "Python module scoring, ssGSEA/singscore-style scoring, Fisher/OR tests, spatial score summaries, HLA execution board.",
        "claim": "HT-overlap/immune-active PTC has HLA-II/IFN/TLS/APC biology across expression and spatial resources.",
        "meaning": "Biology branch is real; allele-genetics branch is not ready.",
        "boundary": "Do not claim HLA allele association from cancer cohorts or n=9 vs n=9 GSE286332 allele tests.",
        "paths": "project/results/hla_deepdive_2026_05_08/; project/papers_hub_2026_0509/paper2_hla.html",
        "tokens": "ht hashimoto hla tls paper2_hla gse286332 aitd",
    },
    {
        "slot_type": "paper",
        "branch": "paper3-ici-readiness",
        "paper": "Paper 3",
        "title": "ICI vulnerability/readiness in molecularly dark thyroid cancer",
        "topic": "ICI / immune therapy readiness",
        "opened": "2026-05-04",
        "last_update": "2026-05-10",
        "state": "FROZEN_DESIGN",
        "parent": "paper1-main",
        "history": "Track A frozen -> Track B-lite public-cohort module scoring -> response-prediction boundary locked.",
        "data": "Public ICI cohorts, GSE76039/GSE33630/GSE60542/GSE126698 thyroid cohorts, GSE151179 RAI before/after, ATC/THCA immune modules.",
        "software": "Python module scoring, cohort pooling, response association tests, forest summaries, sign-consistency tables.",
        "claim": "Current data support immune vulnerability/readiness patterns, not thyroid ICI response prediction.",
        "meaning": "Future branch or Paper 11 reserve. It should not slow Paper 1.",
        "boundary": "No thyroid ICI response predictor claim without thyroid ICI-treated response data.",
        "paths": "project/results/paper3_ici_track_b_lite/phase2_summary.json; project/results/paper3_ici_track_b_lite/phase4_summary.json; project/papers_hub_2026_0509/paper3_ici_status_dossier_2026_05_09.html",
        "tokens": "paper3 ici checkpoint phase_C_ICI track18 immune",
    },
    {
        "slot_type": "paper",
        "branch": "paper4-gd-hla",
        "paper": "Paper 4",
        "title": "Korean GD / Pan-Asian HLA architecture",
        "topic": "GD / Pan-Asian HLA genetics",
        "opened": "2026-05-04",
        "last_update": "2026-05-10",
        "state": "UPGRADE_IF_DATA_ARRIVES",
        "parent": "paper2b-ht-immune",
        "history": "Pan-Asian literature extraction -> Chu/Shin/Chen anchors -> Korean replication readiness ledger.",
        "data": "Chu 2018, Shin 2019, Chen 2011 HLA allele rows, AFND/literature references, Korean replication cohort plan.",
        "software": "arcasHLA where applicable, Python meta-analysis, random-effects screening, BH-FDR tables, source extractability ledger.",
        "claim": "Pan-Asian GD HLA architecture is coherent around DPB1*05:01/B*46:01/A*02:07/C*01:02 signals.",
        "meaning": "High-upside branch if Korean adult GD germline NGS arrives.",
        "boundary": "Current result is not a definitive Korean GD discovery; de novo matched case-control data are needed.",
        "paths": "project/results/hla_deepdive_2026_05_08/track1_paper4_gd_panasian/; project/papers_hub_2026_0509/paper4_gd_hla.html",
        "tokens": "paper4 gd graves hla_ncomm panasian arcas",
    },
    {
        "slot_type": "paper",
        "branch": "paper5-npj-rai8",
        "paper": "Paper 5",
        "title": "npj-style 8-gene RAI / differentiation biomarker",
        "topic": "RAI biomarker / 8-gene clinical model",
        "opened": "2026-04-30",
        "last_update": "2026-05-09",
        "state": "LEGACY_SUPPORT",
        "parent": "paper1-main",
        "history": "v17 npj 8-gene biomarker pages -> robustness/clinical/decision-curve artifacts -> Paper 1 support layer.",
        "data": "TCGA-THCA, external GPL570 cohorts, K2/Lee validation, RAI/TDS panels, survival/clinical endpoints.",
        "software": "Python effect-size scripts, ROC/AUC, Cox/survival, calibration/DCA tables, figure atlas builder.",
        "claim": "A compact thyroid-differentiation/RAI-related score supports clinical readability of the DM1 trunk.",
        "meaning": "Support/release branch, not the current main manuscript.",
        "boundary": "Do not let legacy npj framing override current Paper 1 DM1 dark-matter claim hierarchy.",
        "paths": "project/results/v17*; project/papers_hub_2026_0509/paper5.html; project/papers_hub_2026_0509/portfolio_paper0.html",
        "tokens": "paper5 portfolio_paper0 npj v17 rai biomarker korean_k2",
    },
    {
        "slot_type": "paper",
        "branch": "paper6-dial",
        "paper": "Paper 6",
        "title": "DIAL / Bioinformatics model theorem branch",
        "topic": "Methods / DIAL framework",
        "opened": "2026-05-04",
        "last_update": "2026-05-09",
        "state": "METHODS_PARK",
        "parent": "paper5-npj-rai8",
        "history": "v15/v52 DIAL leak finding -> Paper B / Paper N method pages -> theorem-style reserve.",
        "data": "TCGA/K2/Lee-style tabular benchmarks, DIAL theorem artifacts, leakage audit files.",
        "software": "Python ML scripts, sklearn pipelines, leak diagnostics, theorem/demo HTML pages.",
        "claim": "DIAL-style methods and leak audits can become a methods paper if decoupled from Paper 1.",
        "meaning": "Methods side branch; useful but not paper-blocking for Paper 1.",
        "boundary": "Do not mix DIAL theorem claims into the biological manuscript unless only as a computational audit note.",
        "paths": "project/papers_hub_2026_0509/paper6.html; project/papers_hub_2026_0509/paperB.html; project/papers_hub_2026_0509/paperN.html",
        "tokens": "paper6 dial paperB paperN neurips bioinformatics leak",
    },
    {
        "slot_type": "paper",
        "branch": "paper7-agentic-framework",
        "paper": "Paper 7",
        "title": "Agentic research / dossier framework",
        "topic": "Research automation / agentic framework",
        "opened": "2026-05-04",
        "last_update": "2026-05-09",
        "state": "METHODS_DEMO",
        "parent": "paper6-dial",
        "history": "v18 agentic research framework -> hub/dossier automation -> visual review infrastructure.",
        "data": "Hub page inventories, figure/table/metric ledgers, generated dossier outputs, review routes.",
        "software": "Python HTML builders, static CSS, JSON/TSV manifests, automated audit/explanation layers.",
        "claim": "The research workflow can be packaged as an agentic review/dossier framework.",
        "meaning": "Infra/methods demo branch that supports the whole repo.",
        "boundary": "Framework page is operational tooling, not biological evidence.",
        "paths": "project/papers_hub_2026_0509/paper7.html; project/scripts/build_papers_hub_2026_0509.py",
        "tokens": "paper7 agentic research framework dossier hub_search reviewer_tour",
    },
    {
        "slot_type": "paper",
        "branch": "paper9-ic50-perturbation",
        "paper": "Paper 9",
        "title": "IC50 / PRISM / synthetic-lethality prioritization map",
        "topic": "Drug response / perturbation / IC50 proxy",
        "opened": "2026-05-06",
        "last_update": "2026-05-10",
        "state": "WETLAB_REQUIRED",
        "parent": "paper1-main",
        "history": "Paper 9 perturbation evidence -> synthetic lethality summary -> PRISM/DepMap actionability reserve.",
        "data": "DepMap dependency, PRISM drug-response LFC/AUC proxies, MAPK-inhibitor hits, SLC1A5/GLUD1/GLS metabolic candidates, spatial target overlays.",
        "software": "Python correlation/effect-size scripts, FDR enrichment, Spearman rho, Cohen's d, PRISM/DepMap matrix parsing.",
        "claim": "DM1-high/dedifferentiated state prioritizes MAPK-inhibitor sensitivity and glutamine/NAMPT/MYC candidate vulnerabilities.",
        "meaning": "IC50/drug-response map branch. Useful as actionability/reviewer reserve.",
        "boundary": "PRISM LFC is a drug-response proxy, not direct clinical IC50. No validated synthetic lethality without perturbation/rescue wet-lab.",
        "paths": "project/results/p_synthetic_lethality_2026_05_09/synthetic_lethality_summary.json; project/papers_hub_2026_0509/paper1_synthetic_lethality_dossier.html",
        "tokens": "paper9 synthetic prism depmap drug ic50 lethality perturbation vulnerability",
    },
    {
        "slot_type": "paper",
        "branch": "paper10-atlas",
        "paper": "Paper 10",
        "title": "DM1 synthetic-lethal / cell-line atlas",
        "topic": "Atlas / druggability resource",
        "opened": "2026-05-06",
        "last_update": "2026-05-10",
        "state": "SUPPORT_BRANCH",
        "parent": "paper9-ic50-perturbation",
        "history": "Paper 10 atlas pages -> cancer gene matrix -> vulnerability/druggability support.",
        "data": "Cell-line atlas, druggability candidates, cancer-gene matrix, PRISM/DepMap support rows.",
        "software": "Python atlas generation, matrix summaries, static HTML visual summaries.",
        "claim": "Atlas pages organize druggability and candidate vulnerability context for Paper 9/11.",
        "meaning": "Resource branch; useful as supplement/support rather than standalone during marathon.",
        "boundary": "Atlas co-occurrence is not target validation.",
        "paths": "project/results/paper10_atlas/; project/papers_hub_2026_0509/paper10_atlas.html",
        "tokens": "paper10 atlas cancer_gene_matrix druggability vulnerability",
    },
    {
        "slot_type": "paper",
        "branch": "paper11-pancancer",
        "paper": "Paper 11",
        "title": "Pan-cancer portable DM1 axis",
        "topic": "Pan-cancer DM1 portability",
        "opened": "2026-05-08",
        "last_update": "2026-05-10",
        "state": "NEXT_MAJOR_AFTER_P1",
        "parent": "paper1-main",
        "history": "Pan-cancer summary -> phase A-H modules -> Nature master dossier -> lineage mitigation.",
        "data": "TCGA pan-cancer, DepMap, PRISM, hallmark/pathway overlays, survival meta, lineage-specific mitigation, Korean multi-cohort support.",
        "software": "Python matrix scoring, lineage-wise correlations, Cox/survival tables, DepMap/PRISM overlays, matplotlib figures.",
        "claim": "The thyroid-derived DM1/dedifferentiation architecture ports to multiple cancer lineages and links to dependency/drug-response axes.",
        "meaning": "Best major branch after Paper 1.",
        "boundary": "Lineage-specificity and residual confounding must be presented up front.",
        "paths": "project/results/paper11_pancancer/summary.json; project/results/paper11_pancancer/nature_master_summary_v2.json; project/papers_hub_2026_0509/paper11_pancancer.html",
        "tokens": "paper11 pancancer nature_master mitigation lineage depmap prism",
    },
    {
        "slot_type": "paper",
        "branch": "paper12-network",
        "paper": "Paper 12",
        "title": "Candidate co-expression network architecture",
        "topic": "Network / module architecture",
        "opened": "2026-05-10",
        "last_update": "2026-05-10",
        "state": "ABSORB_INTO_P11",
        "parent": "paper11-pancancer",
        "history": "paper12 network summary -> candidate modules -> top edges -> support for Paper 11.",
        "data": "35 candidate genes, module summary, correlation matrix, top 137 network edges, cancer gene matrix support.",
        "software": "Python correlation networks, module summaries, edge ranking, static network figures.",
        "claim": "Network/module structure explains candidate architecture around LYN/NAMPT/TACSTD2/KCNN4/PAX8.",
        "meaning": "Support branch for Paper 11, not a separate near-term manuscript.",
        "boundary": "Bulk/network association is not causal mechanism by itself.",
        "paths": "project/results/paper12_network/summary.json; project/papers_hub_2026_0509/paper12_network.html",
        "tokens": "paper12 network module edges cancer_gene_matrix",
    },
    {
        "slot_type": "paper",
        "branch": "paper13-spatial-drug",
        "paper": "Paper 13",
        "title": "Spatial drug vulnerability / target-celltype overlap",
        "topic": "Spatial drug vulnerability",
        "opened": "2026-05-09",
        "last_update": "2026-05-10",
        "state": "VISUAL_STRONG_DATA_NEEDS_LOCK",
        "parent": "paper9-ic50-perturbation",
        "history": "spatial drug CV2 overlays -> target/celltype overlap -> match/mismatch decision boards.",
        "data": "GSE250521/GSE230424 spatial overlays, target gene expression, famous drug panels, validation ROI packs, stage/hotspot burden tables.",
        "software": "Spatial module scoring, target-state overlap, ROI gallery generation, CV2 visual summaries, static visual dossiers.",
        "claim": "Spatial maps can show where drug-target/state overlap and vulnerability hotspots occur across thyroid samples.",
        "meaning": "Excellent visual branch; needs locked quantitative controls before manuscript escalation.",
        "boundary": "Spatial co-localization is not drug efficacy; match/mismatch boards must stay explicit.",
        "paths": "project/papers_hub_2026_0509/spatial_drug_cv2_overlay_gallery.html; project/papers_hub_2026_0509/spatial_drug_target_celltype_overlap.html",
        "tokens": "spatial_drug target_celltype vulnerability cv2 roi selpercatinib encorafenib",
    },
    {
        "slot_type": "side_repo",
        "branch": "neoantigen-clean-neo",
        "paper": "Neo/CROSS-Neo",
        "title": "Leakage-aware neoantigen prediction / vaccine platform",
        "topic": "Neoantigen / vaccine methods",
        "opened": "2026-05-09",
        "last_update": "2026-05-10",
        "state": "SEPARATE_UNIVERSE",
        "parent": "root",
        "history": "CROSS-Neo lockdown -> CLEAN-Neo blueprint -> BAR-Neo abstention/interpretability -> MD/TCR audit.",
        "data": "ITSNdb, TESLA, CEDAR, dbPepNeo2, NeoDB, McPAS, IEDB, TCR/pMHC structures, MD audit outputs, PAAD/THCA vaccine triage.",
        "software": "CLEAN-Neo/CROSS-Neo harness, sklearn baselines, ESM2/structure proxy features, public-overlap audits, OpenMM/MD, BMA/abstention queues.",
        "claim": "Current platform supports leakage-aware ranking, abstention, and candidate triage; it does not prove definitive immunogenicity.",
        "meaning": "Separate methods/resource universe; keep apart from THCA Paper 1-4 except as future platform context.",
        "boundary": "Do not imply thyroid vaccine readiness from PAAD/PDAC evidence; source-heldout performance remains heterogeneous.",
        "paths": "project/results/p_neo_bayesian_2026_05_09/; project/results/clean_neobench_barneo_2026_05_09/; project/papers_hub_2026_0509/neoantigen_repo_2026_05_10.html",
        "tokens": "neo vaccine pmhc tcr barneo clean_neobench cross_neo md cancer_vaccine",
    },
]


SOFTWARE_METHODS = [
    ("RNA/module scoring", "pandas, numpy, scipy, ssGSEA/singscore-style z-score/module scoring", "Used for DM1/DM2, HLA/IFN/TLS, thyroid differentiation, MAPK output. Interpret as cohort-level molecular axis, not direct causality."),
    ("Methylation/HM450", "TCGA HM450 beta tables, scipy/statsmodels, effect sizes", "Used for promoter/mean-beta evidence. Interpret as correlative epigenetic silencing support."),
    ("Survival/meta-analysis", "Cox/log-rank style tables, fixed/random-effect summaries, confidence intervals", "Used for OS/PFI/prognosis. Check endpoint, event count, and cohort definition before claim use."),
    ("Deconvolution", "Lu 2023 scRNA pseudobulk, scipy NNLS, ridge-NNLS, LR-clip, sklearn NuSVR", "Used to test whether bulk DM1 signals survive cell-type composition adjustment. NuSVR is primary because NNLS collapses T cells."),
    ("Spatial transcriptomics", "GSE250521/GSE193581 spatial/scRNA score overlays and neighborhood summaries", "Used for anatomy/context/caveat support. Current spatial v15B-v18 is caveat/autopsy, not main Paper 1 proof."),
    ("WSI/image", "OpenSlide-style tiling, DINOv2 ViT-L fallback, CLAM/UNI when available, sklearn LOSO", "Used for H&E projection of DM1. Interpret as small pilot signal, not clinical pathology model."),
    ("DepMap/PRISM drug response", "DepMap dependency matrices, PRISM drug LFC/AUC-like response, Spearman/FDR/Cohen's d", "Used for IC50/drug-response map. Treat PRISM as public screen proxy; wet-lab validation needed."),
    ("Neoantigen benchmark", "CLEAN-Neo/CROSS-Neo harness, sklearn baselines, ESMFold/structure proxy, public-overlap audits", "Used for leakage-aware ranking and abstention. Interpret strict splits before pooled AUROC/AUPRC."),
    ("HTML dossier builder", "project/scripts/build_papers_hub_2026_0509.py, local HTML/CSS, generated JSON/TSV inventories", "Used to make the repo-style review surface. It does not create new biological evidence."),
]


def hrow(cells: Iterable[str], header: bool = False) -> str:
    tag = "th" if header else "td"
    return "<tr>" + "".join(f"<{tag}>{html.escape(str(c))}</{tag}>" for c in cells) + "</tr>"


def link_cell(path: str) -> str:
    safe = html.escape(path)
    if path.startswith("http") or path.startswith("#"):
        return f'<a href="{safe}">{safe}</a>'
    return f'<a href="{safe}">{safe}</a>'


def read_tsv_dicts(path: Path, limit: int | None = None) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if limit is not None:
        return rows[:limit]
    return rows


def build_renewal_block(meta: PageMeta, all_pages: list[PageMeta], idx: dict[str, list[str]]) -> str:
    terms = infer_terms(meta)
    group = infer_paper_group(meta)
    status_label, status_boundary, status_next = current_status_for_group(group)
    def_rows = "\n".join(
        hrow([t["key"], t["type"], t["definition"], t["use"]]) for t in terms
    )
    if not def_rows:
        def_rows = hrow(["No inferred term", "n/a", "Page has no obvious dataset token; see source audit below.", "Manual review required."])

    sibling_links = []
    prefix = meta.path.stem.split("_")[0].split("-")[0]
    for p in all_pages:
        if p.path.name == meta.path.name:
            continue
        if p.path.stem.startswith(prefix) or infer_paper_group(p) == group:
            sibling_links.append(p)
        if len(sibling_links) >= 8:
            break
    siblings = " ".join(
        f'<a class="renewal0509-chip" href="{html.escape(p.path.name)}">{html.escape(p.path.stem[:34])}</a>'
        for p in sibling_links
    )

    heading_preview = " → ".join([h[1] for h in meta.headings[:5]]) or "No h1-h3 headings detected."
    flow_cards = build_flow_cards(meta, limit=8)
    figure_cards = build_figure_cards(meta, idx=idx, limit=8)
    table_cards = build_table_cards(meta, limit=6)
    metric_strip = build_metric_strip(meta, limit=12)
    term_chips = " ".join(
        f'<span class="renewal0509-term">{html.escape(t["key"])}</span>' for t in terms[:10]
    )

    return f"""
<!-- RENEWAL_0509_BLOCK_START -->
<section class="renewal0509 renewal0509-top" id="renewal-0509-data-definitions">
  <style>
    .renewal0509{{--ink:#102033;--muted:#52606f;--paper:#f6f0e6;--card:#fffaf1;--line:#d4c6b3;--line2:#eadfcc;--red:#8f2d25;--blue:#244e73;--green:#426b50;--gold:#b58534;--coal:#17212f;font-family:'Newsreader','Noto Sans KR',serif;color:var(--ink);background:linear-gradient(135deg,#fffaf1 0%,#f1e2c9 100%);border:2px solid var(--coal);box-shadow:8px 8px 0 rgba(23,33,47,.14);margin:0 auto 28px;padding:0;line-height:1.6;max-width:1280px}}
    .renewal0509 *{{box-sizing:border-box}}
    .renewal0509 a{{color:var(--red);text-decoration:none;border-bottom:1px solid rgba(143,45,37,.28)}}
    .renewal0509 .renewal0509-hero{{padding:34px 38px 24px;position:relative;overflow:hidden;background:linear-gradient(145deg,#fffaf1 0%,#f0dfc4 100%)}}
    .renewal0509 .renewal0509-hero:after{{content:"";position:absolute;right:-72px;top:-72px;width:230px;height:230px;border:32px solid rgba(143,45,37,.10);border-radius:50%;pointer-events:none}}
    .renewal0509 .renewal0509-kicker,.renewal0509 .renewal0509-mini{{font-family:'JetBrains Mono',monospace;letter-spacing:.14em;text-transform:uppercase;font-weight:800}}
    .renewal0509 .renewal0509-kicker{{font-size:11px;color:var(--red);margin-bottom:12px}}
    .renewal0509 h2{{font-family:'Cormorant Garamond','Newsreader',serif;font-size:44px;line-height:1.05;margin:0 0 10px;color:var(--coal);max-width:1020px}}
    .renewal0509 h3{{font-family:'Cormorant Garamond','Newsreader',serif;font-size:30px;line-height:1.1;margin:26px 0 12px;color:var(--red)}}
    .renewal0509 h4{{font-family:'Cormorant Garamond','Newsreader',serif;font-size:24px;line-height:1.12;margin:4px 0 8px;color:var(--coal)}}
    .renewal0509 p{{margin:7px 0 11px}}
    .renewal0509 .renewal0509-lead{{font-size:18px;color:var(--muted);max-width:1040px;font-style:italic}}
    .renewal0509 .renewal0509-grid{{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));gap:10px;margin:20px 0 0}}
    .renewal0509 .renewal0509-stat{{background:#fff;border:1px solid var(--line2);border-radius:10px;padding:12px 14px;box-shadow:0 5px 16px rgba(23,33,47,.05)}}
    .renewal0509 .renewal0509-stat b{{display:block;font-family:'Cormorant Garamond',serif;font-size:34px;color:var(--red);line-height:1}}
    .renewal0509 .renewal0509-stat span{{font-family:'JetBrains Mono',monospace;font-size:10px;color:var(--muted);text-transform:uppercase;letter-spacing:.07em}}
    .renewal0509 .renewal0509-body{{padding:24px 34px 34px;background:rgba(255,250,241,.82)}}
    .renewal0509 .renewal0509-block{{background:rgba(255,255,255,.76);border:1px solid var(--line);border-radius:14px;padding:18px;margin:0 0 18px}}
    .renewal0509 .renewal0509-block-title{{display:flex;align-items:baseline;justify-content:space-between;gap:12px;margin-bottom:10px}}
    .renewal0509 .renewal0509-block-title span{{font-family:'JetBrains Mono',monospace;font-size:10px;color:var(--muted);text-transform:uppercase;letter-spacing:.08em}}
    .renewal0509 table{{width:100%;border-collapse:collapse;background:#fff;font-size:12px;margin:8px 0 12px}}
    .renewal0509 th,.renewal0509 td{{border:1px solid var(--line2);padding:8px 9px;vertical-align:top;text-align:left}}
    .renewal0509 th{{background:var(--coal);color:#fff;font-family:'JetBrains Mono',monospace;font-size:10px;text-transform:uppercase;letter-spacing:.05em}}
    .renewal0509 details{{margin-top:8px}}
    .renewal0509 summary{{cursor:pointer;font-family:'JetBrains Mono',monospace;font-size:12px;color:var(--red);font-weight:800}}
    .renewal0509-chip,.renewal0509-term{{display:inline-block;margin:3px 5px 3px 0;padding:5px 8px;border:1px solid var(--line);border-radius:999px;background:#fff;font-family:'JetBrains Mono',monospace;font-size:10px;color:var(--coal)!important}}
    .renewal0509-term{{border-radius:6px;border-color:#d9aaa3;background:#f7dedb;color:var(--red)!important;font-weight:800}}
    .renewal0509-metrics{{display:flex;flex-wrap:wrap;gap:7px;margin:14px 0 0}}
    .renewal0509-metric{{display:inline-block;background:#fff;border:1px solid var(--line);border-radius:6px;padding:6px 8px;font-family:'JetBrains Mono',monospace;font-size:10px;color:var(--coal);font-weight:800}}
    .renewal0509-metric.muted{{color:var(--muted);font-weight:600}}
    .renewal0509-callout{{border-left:5px solid var(--gold);background:#fff8e8;padding:14px 16px;margin:12px 0;border-radius:0 10px 10px 0}}
    .renewal0509-status{{display:grid;grid-template-columns:1fr 2fr 2fr;gap:10px;margin:14px 0 0}}
    .renewal0509-status div{{background:var(--coal);color:#f7efe0;border-radius:10px;padding:12px 14px}}
    .renewal0509-status b{{display:block;font-family:'JetBrains Mono',monospace;font-size:10px;letter-spacing:.08em;text-transform:uppercase;color:#ffd28a;margin-bottom:4px}}
    .renewal0509-status span{{display:block;font-size:12px;color:#f3f6fa}}
    .renewal0509-flow{{display:grid;grid-template-columns:repeat(auto-fit,minmax(230px,1fr));gap:10px}}
    .renewal0509-step{{background:#fff;border:1px solid var(--line2);border-radius:12px;padding:14px 14px 14px 58px;min-height:96px;position:relative}}
    .renewal0509-step b{{position:absolute;left:14px;top:14px;width:32px;height:32px;border-radius:50%;background:var(--coal);color:#fff;display:grid;place-items:center;font-family:'JetBrains Mono',monospace;font-size:12px}}
    .renewal0509-step em{{display:block;font-family:'JetBrains Mono',monospace;font-size:10px;color:var(--red);letter-spacing:.08em;font-style:normal;font-weight:800;margin-bottom:4px}}
    .renewal0509-step span{{display:block;font-size:13px;color:var(--ink)}}
    .renewal0509-figgrid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(290px,1fr));gap:14px}}
    .renewal0509-figcard{{background:#fff;border:1px solid var(--line);border-radius:14px;overflow:hidden;box-shadow:0 6px 18px rgba(23,33,47,.06)}}
    .renewal0509-figframe{{background:#f4eadb;border-bottom:1px solid var(--line2);aspect-ratio:16/9;display:grid;place-items:center;overflow:hidden}}
    .renewal0509-figframe img{{width:100%;height:100%;object-fit:contain;display:block}}
    .renewal0509-figbody{{padding:14px}}
    .renewal0509-mini{{font-size:10px;color:var(--blue);margin-bottom:4px}}
    .renewal0509-figbody p{{font-size:13px;color:var(--muted);line-height:1.45}}
    .renewal0509-figbody dl{{display:grid;grid-template-columns:96px 1fr;gap:6px 8px;margin:10px 0 0;font-size:12px}}
    .renewal0509-figbody dt{{font-family:'JetBrains Mono',monospace;font-size:9px;color:var(--red);text-transform:uppercase;letter-spacing:.06em;font-weight:800}}
    .renewal0509-figbody dd{{margin:0;color:var(--ink)}}
    .renewal0509-tablegrid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(240px,1fr));gap:10px}}
    .renewal0509-tablecard{{background:#fff;border:1px solid var(--line2);border-top:4px solid var(--blue);border-radius:10px;padding:12px}}
    .renewal0509-tablecard b{{font-family:'JetBrains Mono',monospace;font-size:10px;color:var(--blue);letter-spacing:.08em;text-transform:uppercase;display:block;margin-bottom:5px}}
    .renewal0509-tablecard strong{{display:block;color:var(--coal);font-size:14px;margin-bottom:5px}}
    .renewal0509-tablecard span{{display:block;color:var(--muted);font-size:12px;line-height:1.45}}
    .renewal0509-step-more{{border-color:var(--gold);background:#fff8e8}}
    @media(max-width:900px){{.renewal0509{{margin:0 10px 22px}}.renewal0509 .renewal0509-hero,.renewal0509 .renewal0509-body{{padding:22px 18px}}.renewal0509 .renewal0509-grid{{grid-template-columns:1fr 1fr}}.renewal0509-status{{grid-template-columns:1fr}}.renewal0509 h2{{font-size:34px}}}}
  </style>
  <div class="renewal0509-hero">
    <div class="renewal0509-kicker">papers_hub_2026_0509 · paper1-style dossier layer · generated {BUILD_DATE}</div>
    <h2>데이터 정의서 / Data Definition Sheet — {html.escape(meta.title)}</h2>
    <p class="renewal0509-lead">이 페이지는 원본 내용을 유지하면서, 첫 화면에 분석 데이터 정의서, 현재 claim boundary, 읽는 순서, figure/table 해설을 새로 붙인 2026-05-09 리뉴얼 버전입니다. 모든 그림은 “어떤 데이터가 들어갔는지 / 어떻게 읽는지 / 어디까지 주장 가능한지”를 같이 표시합니다.</p>
    <div class="renewal0509-grid">
      <div class="renewal0509-stat"><b>{len(meta.images)}</b><span>figures / visuals</span></div>
      <div class="renewal0509-stat"><b>{meta.table_count}</b><span>tables explained</span></div>
      <div class="renewal0509-stat"><b>{len(meta.headings)}</b><span>flow headings</span></div>
      <div class="renewal0509-stat"><b>{html.escape(group)}</b><span>page family</span></div>
    </div>
    <div class="renewal0509-metrics">{metric_strip}</div>
    <div class="renewal0509-status">
      <div><b>Current status</b><span>{html.escape(status_label)}</span></div>
      <div><b>Claim boundary</b><span>{html.escape(status_boundary)}</span></div>
      <div><b>Next use</b><span>{html.escape(status_next)}</span></div>
    </div>
  </div>
  <div class="renewal0509-body">
    <div class="renewal0509-block">
      <div class="renewal0509-block-title"><h3>0. Analysis Inputs</h3><span>data dictionary first</span></div>
      <p>{term_chips}</p>
      <table>
        {hrow(["Term", "Class", "Definition", "Where this page uses it"], header=True)}
        {def_rows}
      </table>
    </div>
    <div class="renewal0509-block">
      <div class="renewal0509-block-title"><h3>1. Reading Flow</h3><span>section logic</span></div>
      <div class="renewal0509-callout"><b>Flow summary:</b> {html.escape(heading_preview[:700])}</div>
      <div class="renewal0509-flow">{flow_cards}</div>
    </div>
    <div class="renewal0509-block">
      <div class="renewal0509-block-title"><h3>2. Figure Explanations</h3><span>preview; full gallery below</span></div>
      <div class="renewal0509-figgrid">{figure_cards}</div>
    </div>
    <div class="renewal0509-block">
      <div class="renewal0509-block-title"><h3>3. Table Explanations</h3><span>embedded data rows</span></div>
      <div class="renewal0509-tablegrid">{table_cards}</div>
    </div>
    <div class="renewal0509-block">
      <div class="renewal0509-block-title"><h3>4. Related Renewed Pages</h3><span>same family / same prefix</span></div>
      <p>{siblings or '<span class="renewal0509-chip">No immediate sibling detected</span>'}</p>
      <p style="font-family:'JetBrains Mono',monospace;font-size:10px;color:#52606f">Current situation: <a href="current_situation_2026_0509.html">current_situation_2026_0509.html</a> · Source page: project/papers_hub_2026_05_04/{html.escape(meta.path.name)} · Renewed page: project/papers_hub_2026_0509/{html.escape(meta.path.name)} · Voice-protected manuscript sections were not edited.</p>
    </div>
  </div>
</section>
<!-- RENEWAL_0509_BLOCK_END -->
"""


def build_annex(meta: PageMeta, idx: dict[str, list[str]], all_pages: list[PageMeta]) -> str:
    full_figure_cards = build_figure_cards(meta, idx=idx, limit=None)
    full_table_cards = build_table_cards(meta, limit=None)
    image_rows = []
    for i, img in enumerate(meta.images, 1):
        caption = meta.captions[i - 1] if i - 1 < len(meta.captions) else img.get("alt", "")
        matches = asset_matches(img.get("src", ""), idx)
        match_text = "; ".join(matches) if matches else "No same-basename result file found; use page asset path as provenance."
        image_rows.append(
            "<tr>"
            f"<td>{i}</td>"
            f"<td>{link_cell(img.get('src',''))}</td>"
            f"<td>{html.escape(normalize(caption)[:500] or 'No source caption/alt detected.')}</td>"
            f"<td>{html.escape(guess_visual_data(meta, img))}</td>"
            f"<td>{html.escape(match_text)}</td>"
            "</tr>"
        )
    if not image_rows:
        image_rows.append("<tr><td colspan='5'>No image tags detected in this page.</td></tr>")

    table_rows = []
    for i in range(1, meta.table_count + 1):
        caption = meta.table_captions[i - 1] if i - 1 < len(meta.table_captions) else f"Inline HTML table {i}"
        table_rows.append(
            "<tr>"
            f"<td>{i}</td>"
            f"<td>{html.escape(normalize(caption)[:500])}</td>"
            "<td>Rows are embedded in this HTML file. Treat the table body as the immediate data source; upstream TSV/JSON is listed when linked or named in the page.</td>"
            "</tr>"
        )
    if not table_rows:
        table_rows.append("<tr><td colspan='3'>No HTML tables detected in this page.</td></tr>")

    flow_rows = "\n".join(
        hrow([i, tag.upper(), txt[:240]]) for i, (tag, txt) in enumerate(meta.headings[:80], 1)
    )
    if not flow_rows:
        flow_rows = hrow(["-", "-", "No h1-h3 flow headings detected."])

    source_links = []
    for href in meta.links:
        if href.startswith("#") or href.startswith("javascript:"):
            continue
        if any(tok in href.lower() for tok in ["result", "asset", ".tsv", ".csv", ".json", ".md", ".png", ".pdf", ".html"]):
            source_links.append(href)
        if len(source_links) >= 40:
            break
    link_rows = "\n".join(hrow([i, href]) for i, href in enumerate(source_links, 1))
    if not link_rows:
        link_rows = hrow(["-", "No explicit source/result links detected beyond embedded figures/tables."])

    sibling_cards = []
    group = infer_paper_group(meta)
    for p in all_pages:
        if p.path.name == meta.path.name:
            continue
        if infer_paper_group(p) == group:
            sibling_cards.append(
                f'<a class="renewal0509-card" href="{html.escape(p.path.name)}"><b>{html.escape(p.path.stem[:52])}</b><span>{len(p.images)} figs · {p.table_count} tables</span></a>'
            )
        if len(sibling_cards) >= 12:
            break

    return f"""
<!-- RENEWAL_0509_ANNEX_START -->
<section class="renewal0509 renewal0509-annex" id="renewal-0509-figure-table-audit">
  <div class="renewal0509-kicker">v5 annex · figure/table/source audit</div>
  <h2>Figure · Table · Source Audit</h2>
  <p>아래 audit는 이 HTML 안의 모든 <code>img</code>와 <code>table</code>을 스캔해 만든 provenance layer입니다. 각 그림/표가 어떤 데이터 계열을 받았는지, 어떻게 읽어야 하는지, 어떤 로컬 path를 봐야 하는지, 어디까지가 자동 추정인지 분리합니다.</p>
  <div class="renewal0509-block">
    <div class="renewal0509-block-title"><h3>All Figure Explanation Cards</h3><span>visual-by-visual reading guide</span></div>
    <div class="renewal0509-figgrid">{full_figure_cards}</div>
  </div>
  <div class="renewal0509-block">
    <div class="renewal0509-block-title"><h3>All Table Explanation Cards</h3><span>table-by-table data guide</span></div>
    <div class="renewal0509-tablegrid">{full_table_cards}</div>
  </div>
  <details open>
    <summary>1. 모든 그림 / visual asset audit</summary>
    <table>
      {hrow(["#", "Asset path", "Caption / page explanation", "Data entered", "Matched local result/source"], header=True)}
      {''.join(image_rows)}
    </table>
  </details>
  <details open>
    <summary>2. 모든 표 / table audit</summary>
    <table>
      {hrow(["#", "Caption / heading", "Data entered"], header=True)}
      {''.join(table_rows)}
    </table>
  </details>
  <details open>
    <summary>3. 페이지 흐름 / section map</summary>
    <table>
      {hrow(["#", "Level", "Heading"], header=True)}
      {flow_rows}
    </table>
  </details>
  <details>
    <summary>4. 명시 source/result links</summary>
    <table>
      {hrow(["#", "Link"], header=True)}
      {link_rows}
    </table>
  </details>
  <details>
    <summary>5. Same-family renewed pages</summary>
    <div style="display:grid;grid-template-columns:repeat(auto-fit,minmax(220px,1fr));gap:10px;margin-top:10px">
      {''.join(sibling_cards) or '<span class="renewal0509-chip">No same-family pages detected</span>'}
    </div>
  </details>
  <style>
    .renewal0509-card{{display:block;background:#fff;border:1px solid #eadfcc;border-radius:8px;padding:10px 12px;color:#17212f!important;text-decoration:none!important}}
    .renewal0509-card b{{display:block;font-family:'JetBrains Mono',monospace;font-size:11px;color:#8f2d25}}
    .renewal0509-card span{{display:block;font-size:12px;color:#52606f;margin-top:4px}}
  </style>
</section>
<!-- RENEWAL_0509_ANNEX_END -->
"""


def inject_page(raw: str, meta: PageMeta, idx: dict[str, list[str]], all_pages: list[PageMeta]) -> str:
    body_block = build_renewal_block(meta, all_pages, idx)
    annex = build_annex(meta, idx, all_pages)
    if re.search(r"<body[^>]*>", raw, flags=re.I):
        raw = re.sub(r"(<body[^>]*>)", r"\1\n" + body_block, raw, count=1, flags=re.I)
    else:
        raw = body_block + raw
    if re.search(r"</body>", raw, flags=re.I):
        raw = re.sub(r"</body>", annex + "\n</body>", raw, count=1, flags=re.I)
    else:
        raw += annex
    raw = raw.replace("2026_05_04", "2026_0509")
    raw = raw.replace("8-Papers Hub", "2026-05-09 Renewed Hub")
    return raw


def copy_source_tree() -> None:
    if TARGET.exists():
        shutil.rmtree(TARGET)
    ignore = shutil.ignore_patterns("*.bak-*", "__pycache__", ".DS_Store")
    shutil.copytree(SOURCE, TARGET, ignore=ignore)


def collect_pages() -> list[PageMeta]:
    pages: list[PageMeta] = []
    for path in sorted(TARGET.glob("*.html")):
        raw = strip_old_blocks(path.read_text(encoding="utf-8", errors="ignore"))
        pages.append(parse_page(path, raw))
    return pages


def build_master_data_dictionary(pages: list[PageMeta]) -> str:
    term_rows = []
    for term in DATA_TERMS:
        used = [p.path.name for p in pages if term in infer_terms(p)]
        term_rows.append(
            "<tr>"
            f"<td>{html.escape(term['key'])}</td>"
            f"<td>{html.escape(term['type'])}</td>"
            f"<td>{html.escape(term['definition'])}</td>"
            f"<td>{html.escape(term['use'])}</td>"
            f"<td>{len(used)}</td>"
            "</tr>"
        )
    family_counts: dict[str, int] = {}
    for p in pages:
        family_counts[infer_paper_group(p)] = family_counts.get(infer_paper_group(p), 0) + 1
    family_rows = "\n".join(hrow([k, v]) for k, v in sorted(family_counts.items()))
    return f"""<!DOCTYPE html>
<html lang="ko"><head><meta charset="UTF-8"/><meta name="viewport" content="width=device-width,initial-scale=1.0"/>
<title>Data Definition Sheet 2026-05-09 · THCA Hub</title>{common_style()}</head><body>
<header class="master-hero"><div class="kicker">papers_hub_2026_0509</div><h1>데이터 정의서 / Data Definition Sheet</h1><p>모든 리뉴얼 페이지 앞에 삽입된 v5 data-definition dictionary의 원본입니다. 각 row는 페이지 텍스트/figure path 기반으로 자동 매칭됩니다.</p></header>
<main class="master-wrap">
<section class="panel"><h2>Definition Dictionary</h2><table>{hrow(["Term","Class","Definition","Use","Matched pages"], True)}{''.join(term_rows)}</table></section>
<section class="panel"><h2>Page Families</h2><table>{hrow(["Family","Page count"], True)}{family_rows}</table></section>
<section class="panel"><h2>Build Notes</h2><p>Source: <code>project/papers_hub_2026_05_04</code>. Target: <code>project/papers_hub_2026_0509</code>. Voice-protected manuscript sections were not edited.</p></section>
</main></body></html>"""


def build_current_situation_page(pages: list[PageMeta] | None = None) -> str:
    html_pages = len(pages) if pages is not None else 116
    audited = audit_pages(pages) if pages is not None else []
    figures_audited = sum(len(p.images) for p in audited) if pages is not None else 454
    tables_audited = sum(p.table_count for p in audited) if pages is not None else 839
    family_count = len({infer_paper_group(p) for p in audited}) if pages is not None else 18
    rows = [
        (
            "Paper 1",
            "Active manuscript spine",
            "Fig 8 mechanism layer complete through v18; v13/v14 are the load-bearing reviewer-defense layers.",
            "Use v2/v3/v5/v13/v14. Keep v15A/v15C as reserve. Treat GSE250521 v15B-v18 spatial as caveat/autopsy only.",
        ),
        (
            "Paper 1 voice guard",
            "Protected",
            "Hook, Aim, Discussion 3.1, Limitations 3.4, Cover Letter paragraph 1, and Q9 remain author-keyboard only.",
            "Web/data scaffolding is OK; no generated manuscript prose for protected anchors.",
        ),
        (
            "Paper 2 image-DM1",
            "Bounded pilot",
            "Old RAS-like/FVPTC AUC=1.000 is rejected as TSS artefact; UNI-final LOTO signal remains pilot-level only.",
            "Do not pitch as flagship image predictor until K2/external H&E validation and pathologist QC pass.",
        ),
        (
            "Paper 2 HT/HLA",
            "Separated claim track",
            "HT/PTC immune-expression ecology is separate from GD HLA genetics; cancer-cohort HLA genotype claims are not allowed.",
            "Use expression/spatial validation pages for Paper 2; move GD allele genetics to Paper 4.",
        ),
        (
            "Paper 3 ICI",
            "Vulnerability, not predictor",
            "The current status supports ICI vulnerability/readiness/prioritization, not thyroid ICI response prediction.",
            "Track A remains frozen. Track B requires explicit future start and raw thyroid ICI data to change claim level.",
        ),
        (
            "Paper 4 GD HLA",
            "Backlog / validation-gated",
            "Pan-Asian GD/AITD architecture is promising, but definitive Korean-led claim needs adult Korean GD germline HLA validation.",
            "Keep Paper 4 rows out of Paper 2 input claims.",
        ),
        (
            "Paper 5",
            "Paused reserve",
            "v17 npj 8-gene RAI biomarker package is ship-ready but paused to avoid same-data submission conflict.",
            "Reopen after Paper 1/Yu priority decision.",
        ),
        (
            "Paper 6 / Paper N / DIAL",
            "Forensic correction",
            "v52 LODO leak finding blocks downstream citation of earlier weak numbers; framework may survive as a methods story.",
            "Do not cite pre-v52 performance claims downstream.",
        ),
        (
            "Paper 7",
            "Framework / discipline",
            "Agentic research framework and marathon discipline are useful infrastructure, not manuscript sprint prose.",
            "Korean urgency phrases do not override voice-protected sections.",
        ),
        (
            "Paper 9-12",
            "Hypothesis / atlas layer",
            "Synthetic-lethality, pan-cancer transfer, and network pages are prioritization/atlas outputs, not validated clinical actionability.",
            "Keep caveats co-located: low thyroid cell-line n, low OS events, lineage-distance confounding, bulk co-expression only.",
        ),
        (
            "Neoantigen / vaccine",
            "Dossier ready, proof bounded",
            "Cancer vaccine dossier is strong as a traceable platform review, but local models support ranking/abstention, not definitive immunogenicity.",
            "PAAD vaccine-first is defensible only in resected/MRD low-burden framing; routine THCA vaccine-first is not defensible.",
        ),
    ]
    row_html = "\n".join(
        "<tr>"
        f"<td>{html.escape(a)}</td><td>{html.escape(b)}</td><td>{html.escape(c)}</td><td>{html.escape(d)}</td>"
        "</tr>"
        for a, b, c, d in rows
    )
    next_rows = [
        ("1", "Use this 2026_0509 hub for review, not the 2026_05_04 hub.", "The old hub is preserved as source history."),
        ("2", "For Paper 1: read `paper1.html` → `paper1_fig8_mechanism_dossier.html` → `paper1_reviewer_defense_dashboard.html`.", "This is the current manuscript-blocking path."),
        ("3", "For Paper 2/3/4: use the status dossier pages first, then older pages only as evidence archives.", "Prevents old NO-GO / backlog text from being mistaken for current claims."),
        ("4", "For future upgrades: replace automatic annex text with hand-curated page-specific source tables on the 10 highest-value pages.", "The current build is broad and traceable; the next pass should be depth-first."),
    ]
    next_html = "\n".join(
        f"<tr><td>{html.escape(a)}</td><td>{html.escape(b)}</td><td>{html.escape(c)}</td></tr>"
        for a, b, c in next_rows
    )
    return f"""<!DOCTYPE html>
<html lang="ko"><head><meta charset="UTF-8"/><meta name="viewport" content="width=device-width,initial-scale=1.0"/>
<title>Current Situation 2026-05-09 · THCA Hub</title>{common_style()}</head><body>
<section class="master-def">
  <div class="kicker">current situation · data definition first</div>
  <h2>데이터 정의서 / Data Definition Sheet</h2>
  <p>이 페이지는 분석 데이터를 새로 만들지 않습니다. 2026-05-09 현재 허브/결과물/리뷰어 방어 상태를 paper별로 정리하는 운영 대시보드입니다.</p>
  <details open><summary>Scope</summary><table>{hrow(["Input","Definition","Use"], True)}{hrow(["Hub pages","project/papers_hub_2026_0509/*.html","Current review surface"])}{hrow(["Result files","project/results/* plus page-linked TSV/JSON/MD","Evidence provenance"])}{hrow(["Voice-protected anchors","Hook/Aim/Disc 3.1/Limitations/Cover para 1/Q9","Do not generate manuscript prose"])}</table></details>
</section>
<header class="master-hero">
  <div class="kicker">THCA marathon status · 2026-05-09</div>
  <h1>Current Situation</h1>
  <p>지금 상태는 “전체 허브 리뉴얼 완료, Paper 1 우선, 나머지는 claim-boundary가 더 중요해진 상태”입니다. 아래 표는 페이지가 아니라 의사결정 표입니다.</p>
  <div class="stats">
    <div class="stat"><b>1</b><span>active main paper</span></div>
    <div class="stat"><b>v18</b><span>Paper 1 mechanism layer</span></div>
    <div class="stat"><b>{html_pages}</b><span>renewed html pages</span></div>
    <div class="stat"><b>{figures_audited}</b><span>figures audited</span></div>
    <div class="stat"><b>{tables_audited}</b><span>tables audited</span></div>
    <div class="stat"><b>{family_count}</b><span>page families</span></div>
    <div class="stat"><b>0</b><span>protected prose edits</span></div>
  </div>
</header>
<main class="master-wrap">
  <section class="panel"><h2>Situation Matrix</h2><table>{hrow(["Track","Current status","What is true now","Boundary / next move"], True)}{row_html}</table></section>
  <section class="panel"><h2>Immediate Operating Rules</h2><table>{hrow(["#","Rule","Reason"], True)}{next_html}</table></section>
  <section class="panel"><h2>Fast Links</h2><p><a href="index.html">Renewed index</a> · <a href="paper1.html">Paper 1 benchmark</a> · <a href="paper1_fig8_mechanism_dossier.html">Paper 1 Fig 8 dossier</a> · <a href="paper1_reviewer_defense_dashboard.html">Reviewer defense</a> · <a href="paper2_image_dm1_audit_dossier.html">Paper 2 audit</a> · <a href="paper3_ici_status_dossier_2026_05_09.html">Paper 3 ICI status</a> · <a href="hla_ncomm_upgrade_2026_05_09.html">HLA NComm upgrade</a> · <a href="cancer_vaccine_full_dossier.html">Cancer vaccine dossier</a></p></section>
  <section class="panel"><h2>Atlas Links</h2><p><a href="paper_repo_control_room_2026_05_10.html">Control room</a> · <a href="high_impact_board_2026_05_10.html">High-impact board</a> · <a href="impact_escalation_matrix_2026_05_10.html">Impact matrix</a> · <a href="killer_figure_storyboard_2026_05_10.html">Killer figures</a> · <a href="project_inventory_2026_05_10.html">Project inventory</a> · <a href="ic50_prediction_map_2026_05_10.html">IC50/PRISM map</a> · <a href="paper_repo_map_2026_05_10.html">Paper repo map</a> · <a href="paper_branch_manifest_2026_05_10.html">Branch manifest</a> · <a href="paper_date_topic_index_2026_05_10.html">Date/topic index</a> · <a href="software_method_ledger_2026_05_10.html">Software ledger</a> · <a href="deconvolution_repo_2026_05_10.html">Deconvolution repo</a> · <a href="neoantigen_repo_2026_05_10.html">Neoantigen repo</a> · <a href="hub_search_2026_0509.html">Search</a> · <a href="top22_review_pack_2026_0509.html">Top 22 review pack</a> · <a href="paper1_critical_path_2026_0509.html">Paper 1 critical path</a> · <a href="evidence_qc_dashboard_2026_0509.html">Evidence QC</a> · <a href="claim_boundary_ledger_2026_0509.html">Claim ledger</a> · <a href="figure_data_matrix_2026_0509.html">Figure data matrix</a> · <a href="figure_atlas_2026_0509.html">Figure atlas</a> · <a href="table_atlas_2026_0509.html">Table atlas</a> · <a href="flow_atlas_2026_0509.html">Flow atlas</a> · <a href="metric_atlas_2026_0509.html">Metric atlas</a> · <a href="source_atlas_2026_0509.html">Source atlas</a> · <a href="reviewer_tour_2026_0509.html">Reviewer tour</a> · <a href="family_index_2026_0509.html">Family dossiers</a></p></section>
</main></body></html>"""


def common_style() -> str:
    return """
<style>
*{box-sizing:border-box}body{margin:0;background:#0c1422;color:#e9eef7;font-family:'Noto Sans KR','Inter',Arial,sans-serif;line-height:1.58}
a{color:#ffd28a;text-decoration:none}.master-hero{padding:56px 30px;background:linear-gradient(135deg,#101d32,#182844 58%,#0c1422);border-bottom:1px solid #263750}
.master-hero h1{font-family:'Cormorant Garamond','Newsreader',serif;font-size:64px;line-height:1;margin:0 0 10px;color:#fff8e7}.master-hero p{max-width:980px;color:#cdd6e3;font-size:18px}.kicker{font-family:'JetBrains Mono',monospace;color:#ffd28a;letter-spacing:.22em;text-transform:uppercase;font-size:11px;font-weight:800;margin-bottom:12px}
.stats{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:10px;margin-top:20px}.stat{background:rgba(255,255,255,.05);border:1px solid rgba(255,210,138,.2);border-radius:10px;padding:13px}.stat b{display:block;color:#ffd28a;font-size:30px;font-family:'Cormorant Garamond',serif;line-height:1}.stat span{font-family:'JetBrains Mono',monospace;font-size:10px;color:#a8b4c4;text-transform:uppercase;letter-spacing:.08em}
.master-wrap{max-width:1320px;margin:0 auto;padding:28px 24px 80px}.panel{border:1px solid #263750;border-radius:12px;background:#101a2c;margin:0 0 20px;padding:20px}.panel h2{font-family:'Cormorant Garamond','Newsreader',serif;font-size:34px;color:#fff8e7;margin:0 0 8px}
.master-def{background:#fffaf1;color:#102033;border-bottom:2px solid #8f2d25;padding:20px 28px}.master-def h2{font-family:'Cormorant Garamond','Newsreader',serif;font-size:36px;line-height:1.05;margin:0 0 8px;color:#17212f}.master-def p{max-width:1100px}.master-def table{background:#fff;color:#102033}.master-def th{background:#17212f;color:#fff}.master-def td{border-color:#eadfcc}.master-def .kicker{color:#8f2d25}.master-def summary{cursor:pointer;font-family:'JetBrains Mono',monospace;font-size:12px;color:#8f2d25;font-weight:800;margin:8px 0}
table{width:100%;border-collapse:collapse;background:#0d1728;font-size:13px}th,td{border:1px solid #263750;padding:8px 10px;vertical-align:top;text-align:left}th{background:#16213a;color:#ffd28a;font-family:'JetBrains Mono',monospace;font-size:11px;text-transform:uppercase}.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(260px,1fr));gap:14px}.card{display:block;border:1px solid #263750;border-radius:10px;background:#0d1728;padding:14px;color:#e9eef7}.card b{display:block;color:#ffd28a}.card span{display:block;color:#a8b4c4;font-size:12px;margin-top:6px}.tag{font-family:'JetBrains Mono',monospace;font-size:10px;border:1px solid #31435e;border-radius:999px;padding:3px 7px;color:#a8b4c4;display:inline-block;margin:4px 4px 0 0}
code{font-family:'JetBrains Mono',monospace;color:#ffd28a}@media(max-width:760px){.master-hero h1{font-size:42px}.master-wrap{padding:18px 12px}table{font-size:12px}}
</style>
"""


def atlas_style() -> str:
    return """
<style>
:root{--ink:#102033;--muted:#52606f;--paper:#f6f0e6;--card:#fffaf1;--line:#d4c6b3;--line2:#eadfcc;--red:#8f2d25;--blue:#244e73;--green:#426b50;--gold:#b58534;--coal:#17212f}
*{box-sizing:border-box}body{margin:0;background:radial-gradient(circle at 10% 0%,rgba(181,133,52,.18),transparent 28rem),linear-gradient(135deg,#f8f2e8 0%,#efe3d1 55%,#f7efe0 100%);color:var(--ink);font-family:'Newsreader','Noto Sans KR',serif;line-height:1.62}
a{color:var(--red);text-decoration:none;border-bottom:1px solid rgba(143,45,37,.32)}.wrap{max-width:1340px;margin:0 auto;padding:30px 24px 80px}.topnav{display:flex;flex-wrap:wrap;gap:10px;align-items:center;justify-content:space-between;font-family:'JetBrains Mono',monospace;font-size:11px;color:var(--muted);margin-bottom:18px}.topnav a{border:1px solid var(--line);background:rgba(255,255,255,.65);border-radius:999px;padding:8px 10px}
.hero{border:2px solid var(--coal);background:linear-gradient(145deg,#fffaf1 0%,#f0dfc4 100%);padding:40px 42px;box-shadow:10px 10px 0 rgba(23,33,47,.16);position:relative;overflow:hidden}.hero:after{content:"";position:absolute;right:-80px;top:-80px;width:260px;height:260px;border:38px solid rgba(143,45,37,.12);border-radius:50%}.kicker,.mono{font-family:'JetBrains Mono',monospace}.kicker{letter-spacing:.18em;text-transform:uppercase;font-size:11px;color:var(--red);font-weight:800;margin-bottom:14px}h1{font-family:'Cormorant Garamond',serif;font-size:56px;line-height:1.02;margin:0 0 12px;color:var(--coal);max-width:980px}.subtitle{font-size:20px;color:var(--muted);max-width:1040px;font-style:italic}
.stats{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:10px;margin-top:20px}.stat{background:#fff;border:1px solid var(--line2);border-radius:10px;padding:13px}.stat b{display:block;font-family:'Cormorant Garamond',serif;color:var(--red);font-size:34px;line-height:1}.stat span{display:block;font-family:'JetBrains Mono',monospace;color:var(--muted);font-size:10px;letter-spacing:.08em;text-transform:uppercase}
section{background:rgba(255,250,241,.86);border:1px solid var(--line);border-radius:14px;padding:24px;margin:22px 0}h2{font-family:'Cormorant Garamond',serif;font-size:36px;line-height:1.08;margin:0 0 12px;color:var(--coal)}h3{font-family:'Cormorant Garamond',serif;font-size:27px;margin:20px 0 10px;color:var(--red)}
.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(300px,1fr));gap:14px}.card{background:#fff;border:1px solid var(--line2);border-radius:12px;padding:14px;box-shadow:0 5px 16px rgba(23,33,47,.05)}.card b{display:block;font-family:'JetBrains Mono',monospace;font-size:11px;letter-spacing:.06em;text-transform:uppercase;color:var(--red);margin-bottom:5px}.card strong{display:block;font-size:18px;color:var(--coal);line-height:1.16;margin-bottom:6px}.card span,.card p{display:block;font-size:13px;color:var(--muted);margin:5px 0}.renewal0509-metric{display:inline-block!important;background:#fff;border:1px solid var(--line);border-radius:6px;padding:5px 7px;font-family:'JetBrains Mono',monospace;font-size:10px;color:var(--coal);font-weight:800;margin:3px 4px 3px 0}.renewal0509-metric.muted{color:var(--muted);font-weight:600}
.figcard{display:grid;grid-template-columns:170px 1fr;gap:14px;min-height:168px}.thumb{border:1px solid var(--line2);border-radius:10px;background-color:#f4eadb;background-size:contain;background-position:center;background-repeat:no-repeat;min-height:142px}.dl{display:grid;grid-template-columns:96px 1fr;gap:5px 8px;font-size:12px;margin-top:8px}.dl dt{font-family:'JetBrains Mono',monospace;font-size:9px;color:var(--red);letter-spacing:.06em;text-transform:uppercase;font-weight:800}.dl dd{margin:0;color:var(--ink)}.flow{counter-reset:step;display:grid;gap:10px}.step{counter-increment:step;background:#fff;border:1px solid var(--line2);border-radius:12px;padding:13px 14px 13px 58px;position:relative}.step:before{content:counter(step);position:absolute;left:14px;top:13px;width:32px;height:32px;border-radius:50%;background:var(--coal);color:#fff;display:grid;place-items:center;font-family:'JetBrains Mono',monospace;font-size:12px;font-weight:800}.tag{display:inline-block;margin:3px 5px 3px 0;padding:4px 7px;border:1px solid var(--line);border-radius:999px;background:#fff;font-family:'JetBrains Mono',monospace;font-size:10px;color:var(--coal)}
.warning{border-left:5px solid var(--gold);background:#fff8e8;padding:13px 15px;border-radius:0 10px 10px 0;margin:14px 0}@media(max-width:820px){.wrap{padding:18px 12px}.hero{padding:28px 22px}h1{font-size:40px}.figcard{grid-template-columns:1fr}.thumb{min-height:210px}}
</style>
"""


def build_figure_atlas(pages: list[PageMeta], idx: dict[str, list[str]]) -> str:
    records = []
    for p in audit_pages(pages):
        group = infer_paper_group(p)
        for i, img in enumerate(p.images, 1):
            caption = p.captions[i - 1] if i - 1 < len(p.captions) else img.get("alt", "")
            caption = normalize(caption) or "No caption/alt detected in the source HTML."
            role, how, caveat = visual_reading_notes(p, img, caption)
            src = img.get("src", "")
            matches = asset_matches(src, idx)
            match_text = "; ".join(matches) if matches else "No same-basename result match listed."
            records.append((group, p, i, img, caption, role, how, caveat, match_text))
    cards = []
    for group, p, i, img, caption, role, how, caveat, match_text in records:
        src = html.escape(img.get("src", ""), quote=True)
        cards.append(
            '<article class="card figcard">'
            f'<a class="thumb" href="{html.escape(p.path.name)}" style="background-image:url({src})"></a>'
            "<div>"
            f"<b>{html.escape(group)} · Figure {i:02d}</b>"
            f'<strong><a href="{html.escape(p.path.name)}">{html.escape(Path(img.get("src", "")).name or p.path.name)}</a></strong>'
            f"<span>{html.escape(p.title[:150])}</span>"
            f"<p>{html.escape(caption[:520])}</p>"
            '<dl class="dl">'
            f"<dt>Data</dt><dd>{html.escape(guess_visual_data(p, img))}</dd>"
            f"<dt>Read</dt><dd>{html.escape(how)}</dd>"
            f"<dt>Boundary</dt><dd>{html.escape(caveat)}</dd>"
            f"<dt>Match</dt><dd>{html.escape(match_text[:240])}</dd>"
            "</dl></div></article>"
        )
    family_counts: dict[str, int] = {}
    for group, *_ in records:
        family_counts[group] = family_counts.get(group, 0) + 1
    stats = "".join(f'<span class="tag">{html.escape(k)}: {v}</span>' for k, v in sorted(family_counts.items()))
    return f"""<!DOCTYPE html><html lang="ko"><head><meta charset="UTF-8"/><meta name="viewport" content="width=device-width,initial-scale=1.0"/><title>Figure Atlas 2026-05-09</title>{atlas_style()}</head><body><div class="wrap">
<nav class="topnav"><a href="index.html">Renewed hub</a><a href="table_atlas_2026_0509.html">Table atlas</a><a href="flow_atlas_2026_0509.html">Flow atlas</a><a href="metric_atlas_2026_0509.html">Metric atlas</a><a href="source_atlas_2026_0509.html">Source atlas</a><a href="reviewer_tour_2026_0509.html">Reviewer tour</a><a href="family_index_2026_0509.html">Family dossiers</a></nav>
<header class="hero"><div class="kicker">papers_hub_2026_0509 · master figure atlas</div><h1>All Figure Explanations</h1><p class="subtitle">전체 리뉴얼 허브의 visual asset을 page/family별로 모은 그림 설명 atlas입니다. 각 card는 figure source, caption, data entered, how to read, boundary, source match를 같이 보여줍니다.</p><div class="stats"><div class="stat"><b>{len(records)}</b><span>figure records</span></div><div class="stat"><b>{len(family_counts)}</b><span>families</span></div><div class="stat"><b>0</b><span>protected prose edits</span></div></div></header>
<section><h2>Family Distribution</h2><p>{stats}</p></section><section><h2>Figure Cards</h2><div class="grid">{''.join(cards)}</div></section></div></body></html>"""


def build_table_atlas(pages: list[PageMeta]) -> str:
    records = []
    for p in audit_pages(pages):
        group = infer_paper_group(p)
        for i in range(1, p.table_count + 1):
            caption = p.table_captions[i - 1] if i - 1 < len(p.table_captions) else f"Inline HTML table {i}"
            records.append((group, p, i, normalize(caption)))
    cards = []
    for group, p, i, caption in records:
        terms = ", ".join(t["key"] for t in infer_terms(p)[:5])
        cards.append(
            '<article class="card">'
            f"<b>{html.escape(group)} · Table {i:02d}</b>"
            f'<strong><a href="{html.escape(p.path.name)}">{html.escape(caption[:220])}</a></strong>'
            f"<span>{html.escape(p.title[:160])}</span>"
            f"<p>Rows are embedded in the page HTML. Treat the table body as the immediate data source; linked TSV/JSON/MD paths remain binding where the source page names them.</p>"
            f"<p><span class=\"tag\">{html.escape(terms)}</span></p>"
            "</article>"
        )
    family_counts: dict[str, int] = {}
    for group, *_ in records:
        family_counts[group] = family_counts.get(group, 0) + 1
    stats = "".join(f'<span class="tag">{html.escape(k)}: {v}</span>' for k, v in sorted(family_counts.items()))
    return f"""<!DOCTYPE html><html lang="ko"><head><meta charset="UTF-8"/><meta name="viewport" content="width=device-width,initial-scale=1.0"/><title>Table Atlas 2026-05-09</title>{atlas_style()}</head><body><div class="wrap">
<nav class="topnav"><a href="index.html">Renewed hub</a><a href="figure_atlas_2026_0509.html">Figure atlas</a><a href="flow_atlas_2026_0509.html">Flow atlas</a><a href="metric_atlas_2026_0509.html">Metric atlas</a><a href="source_atlas_2026_0509.html">Source atlas</a><a href="reviewer_tour_2026_0509.html">Reviewer tour</a><a href="family_index_2026_0509.html">Family dossiers</a></nav>
<header class="hero"><div class="kicker">papers_hub_2026_0509 · master table atlas</div><h1>All Table Explanations</h1><p class="subtitle">전체 리뉴얼 허브의 HTML table을 한 곳에 모은 표 설명 atlas입니다. 각 card는 page, table caption, data dictionary terms, row provenance rule을 같이 보여줍니다.</p><div class="stats"><div class="stat"><b>{len(records)}</b><span>table records</span></div><div class="stat"><b>{len(family_counts)}</b><span>families</span></div><div class="stat"><b>0</b><span>protected prose edits</span></div></div></header>
<section><h2>Family Distribution</h2><p>{stats}</p></section><section><h2>Table Cards</h2><div class="grid">{''.join(cards)}</div></section></div></body></html>"""


def build_flow_atlas(pages: list[PageMeta]) -> str:
    cards = []
    for p in audit_pages(pages):
        group = infer_paper_group(p)
        status_label, status_boundary, status_next = current_status_for_group(group)
        flow = build_flow_cards(p, limit=10)
        metrics = build_metric_strip(p, limit=8)
        cards.append(
            '<article class="card">'
            f"<b>{html.escape(group)} · {len(p.headings)} headings</b>"
            f'<strong><a href="{html.escape(p.path.name)}">{html.escape(p.title[:160])}</a></strong>'
            f"<p>{metrics}</p>"
            f'<div class="warning"><b>{html.escape(status_label)}</b><br>{html.escape(status_boundary)}<br>{html.escape(status_next)}</div>'
            f'<div class="flow">{flow}</div>'
            "</article>"
        )
    return f"""<!DOCTYPE html><html lang="ko"><head><meta charset="UTF-8"/><meta name="viewport" content="width=device-width,initial-scale=1.0"/><title>Flow Atlas 2026-05-09</title>{atlas_style()}</head><body><div class="wrap">
<nav class="topnav"><a href="index.html">Renewed hub</a><a href="figure_atlas_2026_0509.html">Figure atlas</a><a href="table_atlas_2026_0509.html">Table atlas</a><a href="metric_atlas_2026_0509.html">Metric atlas</a><a href="source_atlas_2026_0509.html">Source atlas</a><a href="reviewer_tour_2026_0509.html">Reviewer tour</a><a href="family_index_2026_0509.html">Family dossiers</a></nav>
<header class="hero"><div class="kicker">papers_hub_2026_0509 · master flow atlas</div><h1>All Page Reading Flows</h1><p class="subtitle">각 페이지의 h1/h2/h3 흐름, key metrics, current claim boundary를 한 화면에서 비교하는 reading-flow atlas입니다.</p><div class="stats"><div class="stat"><b>{len(cards)}</b><span>page flows</span></div><div class="stat"><b>{sum(len(p.headings) for p in audit_pages(pages))}</b><span>headings mapped</span></div><div class="stat"><b>0</b><span>protected prose edits</span></div></div></header>
<section><h2>Flow Cards</h2><div class="grid">{''.join(cards)}</div></section></div></body></html>"""


def build_top22_review_pack(pages: list[PageMeta], idx: dict[str, list[str]]) -> str:
    source_pages = audit_pages(pages)
    top_pages = sorted(source_pages, key=lambda p: (-review_priority_score(p), p.path.name))[:22]

    page_cards = []
    for rank, p in enumerate(top_pages, 1):
        group = infer_paper_group(p)
        status_label, status_boundary, status_next = current_status_for_group(group)
        terms = " ".join(f'<span class="tag">{html.escape(t["key"])}</span>' for t in infer_terms(p)[:5])
        page_cards.append(
            '<article class="card">'
            f"<b>Rank {rank:02d} · score {review_priority_score(p)} · {html.escape(group)}</b>"
            f'<strong><a href="{html.escape(p.path.name)}">{html.escape(p.title[:180])}</a></strong>'
            f"<p>{terms}</p>"
            f"<p>{build_metric_strip(p, limit=10)}</p>"
            f'<div class="warning"><b>{html.escape(status_label)}</b><br>{html.escape(status_boundary)}<br>{html.escape(status_next)}</div>'
            "</article>"
        )

    figure_records = []
    for p in top_pages:
        for i, img in enumerate(p.images, 1):
            caption = p.captions[i - 1] if i - 1 < len(p.captions) else img.get("alt", "")
            caption = normalize(caption) or "No caption/alt detected."
            role, how, caveat = visual_reading_notes(p, img, caption)
            figure_records.append((review_priority_score(p), p, i, img, caption, role, how, caveat))
    figure_records = sorted(figure_records, key=lambda x: (-x[0], x[1].path.name, x[2]))[:22]
    fig_cards = []
    for _, p, i, img, caption, role, how, caveat in figure_records:
        src = html.escape(img.get("src", ""), quote=True)
        matches = asset_matches(img.get("src", ""), idx)
        match_text = "; ".join(matches[:3]) if matches else "No same-basename result match listed."
        fig_cards.append(
            '<article class="card figcard">'
            f'<a class="thumb" href="{html.escape(p.path.name)}" style="background-image:url({src})"></a>'
            "<div>"
            f"<b>{html.escape(infer_paper_group(p))} · Figure {i:02d} · {html.escape(role)}</b>"
            f'<strong><a href="{html.escape(p.path.name)}">{html.escape(Path(img.get("src", "")).name or p.path.name)}</a></strong>'
            f"<p>{html.escape(caption[:460])}</p>"
            f'<dl class="dl"><dt>Data</dt><dd>{html.escape(guess_visual_data(p, img))}</dd><dt>Read</dt><dd>{html.escape(how)}</dd><dt>Boundary</dt><dd>{html.escape(caveat)}</dd><dt>Match</dt><dd>{html.escape(match_text[:220])}</dd></dl>'
            "</div></article>"
        )

    table_records = []
    for p in top_pages:
        for i in range(1, p.table_count + 1):
            caption = p.table_captions[i - 1] if i - 1 < len(p.table_captions) else f"Inline HTML table {i}"
            table_records.append((review_priority_score(p), p, i, normalize(caption)))
    table_records = sorted(table_records, key=lambda x: (-x[0], x[1].path.name, x[2]))[:22]
    table_cards = []
    for _, p, i, caption in table_records:
        table_cards.append(
            '<article class="card">'
            f"<b>{html.escape(infer_paper_group(p))} · Table {i:02d}</b>"
            f'<strong><a href="{html.escape(p.path.name)}">{html.escape(caption[:240])}</a></strong>'
            f"<p>{html.escape(p.title[:180])}</p>"
            "<p>Embedded HTML rows are the immediate data source; linked TSV/JSON/MD paths remain binding when present.</p>"
            "</article>"
        )

    metric_cards = []
    for p in top_pages:
        metrics = extract_metric_tokens(p, limit=22)
        if not metrics:
            continue
        metric_html = " ".join(f'<span class="renewal0509-metric">{html.escape(m)}</span>' for m in metrics)
        metric_cards.append(
            '<article class="card">'
            f"<b>{html.escape(infer_paper_group(p))} · metrics</b>"
            f'<strong><a href="{html.escape(p.path.name)}">{html.escape(p.title[:170])}</a></strong>'
            f"<p>{metric_html}</p>"
            "</article>"
        )

    family_rows: dict[str, int] = {}
    for p in top_pages:
        family_rows[infer_paper_group(p)] = family_rows.get(infer_paper_group(p), 0) + 1
    family_tags = "".join(f'<span class="tag">{html.escape(k)}: {v}</span>' for k, v in sorted(family_rows.items()))

    return f"""<!DOCTYPE html><html lang="ko"><head><meta charset="UTF-8"/><meta name="viewport" content="width=device-width,initial-scale=1.0"/><title>Top 22 Review Pack 2026-05-09</title>{atlas_style()}</head><body><div class="wrap">
<nav class="topnav"><a href="index.html">Renewed hub</a><a href="reviewer_tour_2026_0509.html">Reviewer tour</a><a href="figure_atlas_2026_0509.html">Figure atlas</a><a href="table_atlas_2026_0509.html">Table atlas</a><a href="metric_atlas_2026_0509.html">Metric atlas</a><a href="source_atlas_2026_0509.html">Source atlas</a></nav>
<header class="hero"><div class="kicker">papers_hub_2026_0509 · top 22 review pack</div><h1>Top 22 Review Pack</h1><p class="subtitle">현재 상황에 맞춰 가장 먼저 봐야 할 22개 페이지, 22개 figure, 22개 table, metric 후보를 한 번에 모은 rapid review pack입니다. 점수는 paper priority, reviewer-defense relevance, figure/table density, heading flow, metric density로 계산했습니다.</p><div class="stats"><div class="stat"><b>22</b><span>priority pages</span></div><div class="stat"><b>{len(figure_records)}</b><span>priority figures</span></div><div class="stat"><b>{len(table_records)}</b><span>priority tables</span></div><div class="stat"><b>0</b><span>protected prose edits</span></div></div></header>
<section><h2>Family Mix</h2><p>{family_tags}</p></section>
<section><h2>Top 22 Pages</h2><div class="grid">{''.join(page_cards)}</div></section>
<section><h2>Top 22 Figures</h2><div class="grid">{''.join(fig_cards)}</div></section>
<section><h2>Top 22 Tables</h2><div class="grid">{''.join(table_cards)}</div></section>
<section><h2>Metric Pull List</h2><div class="grid">{''.join(metric_cards)}</div></section>
</div></body></html>"""


def build_search_manifest(pages: list[PageMeta], idx: dict[str, list[str]]) -> list[dict[str, object]]:
    records: list[dict[str, object]] = []
    for p in audit_pages(pages):
        group = infer_paper_group(p)
        status_label, status_boundary, status_next = current_status_for_group(group)
        records.append(
            {
                "page": p.path.name,
                "title": p.title,
                "family": group,
                "figures": len(p.images),
                "tables": p.table_count,
                "headings": len(p.headings),
                "metrics": extract_metric_tokens(p, limit=18),
                "terms": [t["key"] for t in infer_terms(p)],
                "sources": source_chips_for_page(p, idx, limit=12),
                "status": status_label,
                "boundary": status_boundary,
                "next": status_next,
                "text": p.text_sample[:1600],
                "score": review_priority_score(p),
            }
        )
    return records


def build_hub_search_page(records: list[dict[str, object]]) -> str:
    data = json.dumps(records, ensure_ascii=False)
    return f"""<!DOCTYPE html><html lang="ko"><head><meta charset="UTF-8"/><meta name="viewport" content="width=device-width,initial-scale=1.0"/><title>Hub Search 2026-05-09</title>{atlas_style()}<style>.searchbar{{display:grid;grid-template-columns:1fr 220px 140px;gap:10px;margin:18px 0}}input,select{{width:100%;border:1px solid var(--line);border-radius:8px;background:#fff;padding:11px 12px;font:14px 'Noto Sans KR',Arial,sans-serif;color:var(--ink)}}.result-count{{font-family:'JetBrains Mono',monospace;color:var(--red);font-weight:800}}@media(max-width:760px){{.searchbar{{grid-template-columns:1fr}}}}</style></head><body><div class="wrap">
<nav class="topnav"><a href="index.html">Renewed hub</a><a href="top22_review_pack_2026_0509.html">Top 22</a><a href="figure_atlas_2026_0509.html">Figure atlas</a><a href="claim_boundary_ledger_2026_0509.html">Claim ledger</a><a href="evidence_qc_dashboard_2026_0509.html">QC dashboard</a></nav>
<header class="hero"><div class="kicker">papers_hub_2026_0509 · static search</div><h1>Search The Renewed Hub</h1><p class="subtitle">page title, family, terms, metrics, status, source paths, and visible text sample을 검색합니다. 서버 없이 브라우저에서 바로 필터링됩니다.</p><div class="stats"><div class="stat"><b>{len(records)}</b><span>records</span></div><div class="stat"><b>0</b><span>protected prose edits</span></div></div></header>
<section><h2>Find Evidence</h2><div class="searchbar"><input id="q" placeholder="Search DM1, Fig 8, HLA, AUC, PRISM, source path..."/><select id="family"><option value="">All families</option></select><select id="sort"><option value="score">Priority</option><option value="figures">Figures</option><option value="tables">Tables</option><option value="title">Title</option></select></div><p class="result-count" id="count"></p><div class="grid" id="results"></div></section>
<script id="records" type="application/json">{data}</script>
<script>
const records = JSON.parse(document.getElementById('records').textContent);
const q = document.getElementById('q'), fam = document.getElementById('family'), sort = document.getElementById('sort'), results = document.getElementById('results'), count = document.getElementById('count');
[...new Set(records.map(r=>r.family))].sort().forEach(f=>{{ const o=document.createElement('option'); o.value=f; o.textContent=f; fam.appendChild(o); }});
function card(r){{return `<article class="card"><b>${{r.family}} · score ${{r.score}} · ${{r.figures}} figs · ${{r.tables}} tables</b><strong><a href="${{r.page}}">${{r.title}}</a></strong><p>${{r.status}} · ${{r.boundary}}</p><p>${{(r.metrics||[]).slice(0,10).map(x=>`<span class="renewal0509-metric">${{x}}</span>`).join(' ')}}</p><p>${{(r.terms||[]).slice(0,8).map(x=>`<span class="tag">${{x}}</span>`).join(' ')}}</p></article>`}}
function run(){{const needle=q.value.toLowerCase().trim(); let rows=records.filter(r=>(!fam.value||r.family===fam.value)); if(needle) rows=rows.filter(r=>JSON.stringify(r).toLowerCase().includes(needle)); rows.sort((a,b)=>sort.value==='title'?a.title.localeCompare(b.title):(b[sort.value]||0)-(a[sort.value]||0)); count.textContent=`${{rows.length}} / ${{records.length}} pages`; results.innerHTML=rows.slice(0,120).map(card).join('');}}
[q,fam,sort].forEach(el=>el.addEventListener('input',run)); run();
</script></div></body></html>"""


def build_evidence_qc_dashboard(pages: list[PageMeta], idx: dict[str, list[str]]) -> str:
    rows = []
    cards = []
    grade_counts: dict[str, int] = {}
    for p in sorted(audit_pages(pages), key=lambda x: (-review_priority_score(x), x.path.name)):
        score, grade = quality_score(p, idx)
        grade_counts[grade] = grade_counts.get(grade, 0) + 1
        items = quality_items(p, idx)
        checks = " ".join(f'<span class="tag">{k}: {"yes" if v else "no"}</span>' for k, v in items.items())
        rows.append(hrow([p.path.name, infer_paper_group(p), grade, score, len(p.images), p.table_count, len(p.headings), len(extract_metric_tokens(p, 12)), len(source_chips_for_page(p, idx, 12))]))
        cards.append(
            '<article class="card">'
            f"<b>Grade {grade} · score {score}</b>"
            f'<strong><a href="{html.escape(p.path.name)}">{html.escape(p.title[:170])}</a></strong>'
            f"<p>{checks}</p>"
            "</article>"
        )
    grade_tags = "".join(f'<span class="tag">Grade {k}: {v}</span>' for k, v in sorted(grade_counts.items()))
    return f"""<!DOCTYPE html><html lang="ko"><head><meta charset="UTF-8"/><meta name="viewport" content="width=device-width,initial-scale=1.0"/><title>Evidence QC Dashboard 2026-05-09</title>{atlas_style()}</head><body><div class="wrap">
<nav class="topnav"><a href="index.html">Renewed hub</a><a href="hub_search_2026_0509.html">Search</a><a href="claim_boundary_ledger_2026_0509.html">Claim ledger</a><a href="top22_review_pack_2026_0509.html">Top 22</a></nav>
<header class="hero"><div class="kicker">papers_hub_2026_0509 · evidence QC</div><h1>Evidence QC Dashboard</h1><p class="subtitle">각 페이지가 figures, tables, flow headings, metrics, source links, data definitions를 갖췄는지 자동 점검한 품질 대시보드입니다.</p><div class="stats"><div class="stat"><b>{len(rows)}</b><span>pages scored</span></div><div class="stat"><b>{grade_counts.get('A',0)}</b><span>grade A</span></div><div class="stat"><b>0</b><span>protected prose edits</span></div></div></header>
<section><h2>Grade Mix</h2><p>{grade_tags}</p></section><section><h2>QC Cards</h2><div class="grid">{''.join(cards)}</div></section><section><h2>QC Table</h2><table>{hrow(["Page","Family","Grade","Score","Figures","Tables","Headings","Metrics","Sources"], True)}{''.join(rows)}</table></section></div></body></html>"""


def build_claim_boundary_ledger(pages: list[PageMeta]) -> str:
    rows = []
    cards = []
    for p in sorted(audit_pages(pages), key=lambda x: (infer_paper_group(x), -review_priority_score(x))):
        group = infer_paper_group(p)
        status_label, status_boundary, status_next = current_status_for_group(group)
        rows.append(hrow([p.path.name, group, status_label, status_boundary, status_next, len(p.images), p.table_count, len(extract_metric_tokens(p, 12))]))
        cards.append(
            '<article class="card">'
            f"<b>{html.escape(group)} · {html.escape(status_label)}</b>"
            f'<strong><a href="{html.escape(p.path.name)}">{html.escape(p.title[:170])}</a></strong>'
            f"<p>{html.escape(status_boundary)}</p><p>{html.escape(status_next)}</p>"
            "</article>"
        )
    return f"""<!DOCTYPE html><html lang="ko"><head><meta charset="UTF-8"/><meta name="viewport" content="width=device-width,initial-scale=1.0"/><title>Claim Boundary Ledger 2026-05-09</title>{atlas_style()}</head><body><div class="wrap">
<nav class="topnav"><a href="index.html">Renewed hub</a><a href="hub_search_2026_0509.html">Search</a><a href="evidence_qc_dashboard_2026_0509.html">QC dashboard</a><a href="current_situation_2026_0509.html">Current situation</a></nav>
<header class="hero"><div class="kicker">papers_hub_2026_0509 · claim boundary ledger</div><h1>Claim Boundary Ledger</h1><p class="subtitle">모든 페이지를 현재 claim level과 boundary 기준으로 다시 묶은 ledger입니다. 오래된 페이지도 이 ledger 기준으로 읽어야 합니다.</p><div class="stats"><div class="stat"><b>{len(rows)}</b><span>pages</span></div><div class="stat"><b>0</b><span>protected prose edits</span></div></div></header>
<section><h2>Boundary Cards</h2><div class="grid">{''.join(cards)}</div></section><section><h2>Ledger Table</h2><table>{hrow(["Page","Family","Status","Claim boundary","Next use","Figures","Tables","Metrics"], True)}{''.join(rows)}</table></section></div></body></html>"""


def build_paper1_critical_path(pages: list[PageMeta]) -> str:
    paper1_pages = [p for p in audit_pages(pages) if infer_paper_group(p) == "Paper 1"]
    top = sorted(paper1_pages, key=lambda p: (-review_priority_score(p), p.path.name))
    buckets = [
        ("Start", ["paper1.html", "paper1_fig8_mechanism_dossier.html", "paper1_reviewer_defense_dashboard.html"]),
        ("Mechanism / Fig 8", ["fig8", "deconvolution", "spatial", "synthetic"]),
        ("Methods / Reproducibility", ["methods", "reproducibility", "gpl570"]),
        ("Decision / Journal", ["nature", "board", "master"]),
    ]
    sections = []
    for title, needles in buckets:
        cards = []
        for p in top:
            hay = p.path.name.lower()
            if any(n.lower() in hay for n in needles):
                cards.append(
                    '<article class="card">'
                    f"<b>{len(p.images)} figs · {p.table_count} tables · score {review_priority_score(p)}</b>"
                    f'<strong><a href="{html.escape(p.path.name)}">{html.escape(p.title[:170])}</a></strong>'
                    f"<p>{build_metric_strip(p, limit=8)}</p>"
                    "</article>"
                )
        sections.append(f"<section><h2>{html.escape(title)}</h2><div class=\"grid\">{''.join(cards) or '<article class=\"card\"><b>No matching page</b><p>No page matched this bucket.</p></article>'}</div></section>")
    return f"""<!DOCTYPE html><html lang="ko"><head><meta charset="UTF-8"/><meta name="viewport" content="width=device-width,initial-scale=1.0"/><title>Paper 1 Critical Path 2026-05-09</title>{atlas_style()}</head><body><div class="wrap">
<nav class="topnav"><a href="index.html">Renewed hub</a><a href="paper1.html">Paper 1</a><a href="top22_review_pack_2026_0509.html">Top 22</a><a href="claim_boundary_ledger_2026_0509.html">Claim ledger</a></nav>
<header class="hero"><div class="kicker">paper 1 · critical path · scaffolding only</div><h1>Paper 1 Critical Path</h1><p class="subtitle">Paper 1 review sequence만 묶은 operational page입니다. Hook/Aim/Discussion/Limitations/Q9/Cover paragraph voice-protected prose는 생성하지 않습니다.</p><div class="stats"><div class="stat"><b>{len(paper1_pages)}</b><span>Paper 1 pages</span></div><div class="stat"><b>{sum(len(p.images) for p in paper1_pages)}</b><span>figures</span></div><div class="stat"><b>{sum(p.table_count for p in paper1_pages)}</b><span>tables</span></div></div></header>
{''.join(sections)}</div></body></html>"""


def build_figure_data_matrix(pages: list[PageMeta], idx: dict[str, list[str]]) -> str:
    rows = []
    role_counts: dict[str, int] = {}
    group_counts: dict[str, int] = {}
    for p in audit_pages(pages):
        group = infer_paper_group(p)
        for i, img in enumerate(p.images, 1):
            caption = p.captions[i - 1] if i - 1 < len(p.captions) else img.get("alt", "")
            role, how, caveat = visual_reading_notes(p, img, caption)
            role_counts[role] = role_counts.get(role, 0) + 1
            group_counts[group] = group_counts.get(group, 0) + 1
            rows.append(hrow([p.path.name, group, i, Path(img.get("src","")).name, role, guess_visual_data(p, img), how, caveat, "; ".join(asset_matches(img.get("src",""), idx)[:3]) or "none"]))
    role_cards = "".join(f'<article class="card"><b>{v} figures</b><strong>{html.escape(k)}</strong></article>' for k, v in sorted(role_counts.items(), key=lambda x: -x[1]))
    group_tags = "".join(f'<span class="tag">{html.escape(k)}: {v}</span>' for k, v in sorted(group_counts.items()))
    return f"""<!DOCTYPE html><html lang="ko"><head><meta charset="UTF-8"/><meta name="viewport" content="width=device-width,initial-scale=1.0"/><title>Figure Data Matrix 2026-05-09</title>{atlas_style()}</head><body><div class="wrap">
<nav class="topnav"><a href="index.html">Renewed hub</a><a href="figure_atlas_2026_0509.html">Figure atlas</a><a href="hub_search_2026_0509.html">Search</a><a href="source_atlas_2026_0509.html">Source atlas</a></nav>
<header class="hero"><div class="kicker">papers_hub_2026_0509 · figure data matrix</div><h1>Figure Data Matrix</h1><p class="subtitle">모든 figure를 role, data entered, how to read, boundary, matched source path 기준으로 펼친 matrix입니다.</p><div class="stats"><div class="stat"><b>{len(rows)}</b><span>figures</span></div><div class="stat"><b>{len(role_counts)}</b><span>visual roles</span></div><div class="stat"><b>0</b><span>protected prose edits</span></div></div></header>
<section><h2>Family Distribution</h2><p>{group_tags}</p></section><section><h2>Visual Role Counts</h2><div class="grid">{role_cards}</div></section><section><h2>Matrix</h2><table>{hrow(["Page","Family","#","Asset","Role","Data entered","How to read","Boundary","Source match"], True)}{''.join(rows)}</table></section></div></body></html>"""


def build_metric_atlas(pages: list[PageMeta]) -> str:
    records = []
    for p in audit_pages(pages):
        metrics = extract_metric_tokens(p, limit=18)
        if not metrics:
            continue
        records.append((infer_paper_group(p), p, metrics))
    cards = []
    for group, p, metrics in records:
        metric_html = " ".join(f'<span class="renewal0509-metric">{html.escape(m)}</span>' for m in metrics)
        cards.append(
            '<article class="card">'
            f"<b>{html.escape(group)} · {len(metrics)} extracted tokens</b>"
            f'<strong><a href="{html.escape(p.path.name)}">{html.escape(p.title[:170])}</a></strong>'
            f"<p>{metric_html}</p>"
            "<p>자동 추출된 숫자 후보입니다. manuscript claim로 쓰기 전에는 해당 페이지의 figure/table card와 원본 source row를 같이 확인해야 합니다.</p>"
            "</article>"
        )
    return f"""<!DOCTYPE html><html lang="ko"><head><meta charset="UTF-8"/><meta name="viewport" content="width=device-width,initial-scale=1.0"/><title>Metric Atlas 2026-05-09</title>{atlas_style()}</head><body><div class="wrap">
<nav class="topnav"><a href="index.html">Renewed hub</a><a href="figure_atlas_2026_0509.html">Figure atlas</a><a href="table_atlas_2026_0509.html">Table atlas</a><a href="flow_atlas_2026_0509.html">Flow atlas</a><a href="source_atlas_2026_0509.html">Source atlas</a><a href="reviewer_tour_2026_0509.html">Reviewer tour</a></nav>
<header class="hero"><div class="kicker">papers_hub_2026_0509 · metric atlas</div><h1>Extracted Metrics</h1><p class="subtitle">각 페이지의 caption/text/table caption에서 HR, OR, AUC, rho, Cohen's d, p-value, percentages, ratios, confidence intervals 후보를 자동 추출한 숫자 atlas입니다.</p><div class="stats"><div class="stat"><b>{len(records)}</b><span>pages with metrics</span></div><div class="stat"><b>{sum(len(m) for _, _, m in records)}</b><span>metric tokens</span></div><div class="stat"><b>0</b><span>protected prose edits</span></div></div></header>
<section><h2>Metric Cards</h2><div class="grid">{''.join(cards)}</div></section></div></body></html>"""


def build_source_atlas(pages: list[PageMeta], idx: dict[str, list[str]]) -> str:
    cards = []
    source_count = 0
    for p in audit_pages(pages):
        group = infer_paper_group(p)
        links = []
        for href in p.links:
            if href.startswith("#") or href.startswith("javascript:"):
                continue
            if any(tok in href.lower() for tok in ["result", "asset", ".tsv", ".csv", ".json", ".md", ".png", ".pdf", ".html"]):
                links.append(href)
        image_matches = []
        for img in p.images[:8]:
            matches = asset_matches(img.get("src", ""), idx)
            if matches:
                image_matches.extend(matches[:2])
        items = []
        for value in links[:12] + image_matches[:12]:
            items.append(f'<span class="tag">{html.escape(value[:140])}</span>')
        if not items:
            items.append('<span class="tag">No explicit source/result link detected</span>')
        source_count += len(items)
        cards.append(
            '<article class="card">'
            f"<b>{html.escape(group)} · source links</b>"
            f'<strong><a href="{html.escape(p.path.name)}">{html.escape(p.title[:170])}</a></strong>'
            f"<p>{''.join(items)}</p>"
            "</article>"
        )
    return f"""<!DOCTYPE html><html lang="ko"><head><meta charset="UTF-8"/><meta name="viewport" content="width=device-width,initial-scale=1.0"/><title>Source Atlas 2026-05-09</title>{atlas_style()}</head><body><div class="wrap">
<nav class="topnav"><a href="index.html">Renewed hub</a><a href="figure_atlas_2026_0509.html">Figure atlas</a><a href="table_atlas_2026_0509.html">Table atlas</a><a href="flow_atlas_2026_0509.html">Flow atlas</a><a href="metric_atlas_2026_0509.html">Metric atlas</a><a href="reviewer_tour_2026_0509.html">Reviewer tour</a></nav>
<header class="hero"><div class="kicker">papers_hub_2026_0509 · source atlas</div><h1>Source And Result Paths</h1><p class="subtitle">각 페이지의 explicit source/result links와 figure basename으로 매칭된 local result paths를 모은 provenance atlas입니다.</p><div class="stats"><div class="stat"><b>{len(cards)}</b><span>pages scanned</span></div><div class="stat"><b>{source_count}</b><span>source chips</span></div><div class="stat"><b>0</b><span>protected prose edits</span></div></div></header>
<section><h2>Source Cards</h2><div class="grid">{''.join(cards)}</div></section></div></body></html>"""


def build_reviewer_tour(pages: list[PageMeta]) -> str:
    groups: dict[str, list[PageMeta]] = {}
    for p in audit_pages(pages):
        groups.setdefault(infer_paper_group(p), []).append(p)
    priority_order = [
        "Paper 1",
        "Paper 2",
        "Paper 3",
        "Paper 4",
        "Paper 6 / DIAL",
        "Paper 9",
        "Paper 10",
        "Paper 11",
        "Paper 12",
        "Neoantigen / Vaccine",
        "HLA / Paper 2-4",
        "Spatial / Trajectory",
        "Operations",
        "Other",
    ]
    cards = []
    step_no = 1
    for family in priority_order:
        ps = groups.get(family, [])
        if not ps:
            continue
        status_label, status_boundary, status_next = current_status_for_group(family)
        top_pages = sorted(ps, key=lambda p: (-(len(p.images) * 2 + p.table_count + len(p.headings)), p.path.name))[:8]
        links = " ".join(
            f'<a class="tag" href="{html.escape(p.path.name)}">{html.escape(p.path.stem[:46])}</a>'
            for p in top_pages
        )
        cards.append(
            '<article class="step">'
            f"<b>{step_no:02d}</b><strong>{html.escape(family)}</strong>"
            f"<p>{html.escape(status_label)} · {html.escape(status_boundary)} · {html.escape(status_next)}</p>"
            f"<p>{links}</p>"
            "</article>"
        )
        step_no += 1
    return f"""<!DOCTYPE html><html lang="ko"><head><meta charset="UTF-8"/><meta name="viewport" content="width=device-width,initial-scale=1.0"/><title>Reviewer Tour 2026-05-09</title>{atlas_style()}</head><body><div class="wrap">
<nav class="topnav"><a href="index.html">Renewed hub</a><a href="current_situation_2026_0509.html">Current situation</a><a href="figure_atlas_2026_0509.html">Figure atlas</a><a href="table_atlas_2026_0509.html">Table atlas</a><a href="metric_atlas_2026_0509.html">Metric atlas</a><a href="source_atlas_2026_0509.html">Source atlas</a></nav>
<header class="hero"><div class="kicker">papers_hub_2026_0509 · reviewer reading tour</div><h1>Reviewer Tour</h1><p class="subtitle">지금 상태에 맞춘 빠른 검토 순서입니다. Paper 1을 먼저 보고, 나머지는 claim boundary를 보면서 evidence archive로 읽도록 설계했습니다.</p><div class="stats"><div class="stat"><b>{len(cards)}</b><span>review steps</span></div><div class="stat"><b>{sum(len(v) for v in groups.values())}</b><span>source pages</span></div><div class="stat"><b>0</b><span>protected prose edits</span></div></div></header>
<section><h2>Route</h2><div class="flow">{''.join(cards)}</div></section></div></body></html>"""


def build_family_dossier(family: str, family_pages: list[PageMeta], idx: dict[str, list[str]]) -> str:
    status_label, status_boundary, status_next = current_status_for_group(family)
    page_cards = []
    for p in sorted(family_pages, key=lambda x: (-(len(x.images) + x.table_count), x.path.name)):
        terms = " ".join(f'<span class="tag">{html.escape(t["key"])}</span>' for t in infer_terms(p)[:5])
        page_cards.append(
            '<article class="card">'
            f'<b>{len(p.images)} figs · {p.table_count} tables · {len(p.headings)} headings</b>'
            f'<strong><a href="{html.escape(p.path.name)}">{html.escape(p.title[:170])}</a></strong>'
            f"<p>{terms}</p><p>{build_metric_strip(p, limit=8)}</p>"
            "</article>"
        )
    fig_records = []
    for p in family_pages:
        for i, img in enumerate(p.images[:4], 1):
            caption = p.captions[i - 1] if i - 1 < len(p.captions) else img.get("alt", "")
            caption = normalize(caption) or "No caption/alt detected."
            role, how, caveat = visual_reading_notes(p, img, caption)
            fig_records.append((p, i, img, caption, role, how, caveat))
    fig_cards = []
    for p, i, img, caption, role, how, caveat in fig_records[:36]:
        src = html.escape(img.get("src", ""), quote=True)
        fig_cards.append(
            '<article class="card figcard">'
            f'<a class="thumb" href="{html.escape(p.path.name)}" style="background-image:url({src})"></a>'
            "<div>"
            f"<b>Figure {i:02d} · {html.escape(role)}</b>"
            f'<strong><a href="{html.escape(p.path.name)}">{html.escape(Path(img.get("src", "")).name or p.path.name)}</a></strong>'
            f"<p>{html.escape(caption[:380])}</p>"
            f'<dl class="dl"><dt>Data</dt><dd>{html.escape(guess_visual_data(p, img))}</dd><dt>Read</dt><dd>{html.escape(how)}</dd><dt>Boundary</dt><dd>{html.escape(caveat)}</dd></dl>'
            "</div></article>"
        )
    table_cards = []
    for p in family_pages:
        for i in range(1, min(p.table_count, 4) + 1):
            caption = p.table_captions[i - 1] if i - 1 < len(p.table_captions) else f"Inline HTML table {i}"
            table_cards.append(
                '<article class="card">'
                f"<b>Table {i:02d} · {html.escape(p.path.name)}</b>"
                f'<strong><a href="{html.escape(p.path.name)}">{html.escape(normalize(caption)[:220])}</a></strong>'
                "<p>Embedded table body is the immediate source; linked TSV/JSON/MD remains binding when present.</p>"
                "</article>"
            )
    slug = slugify_family(family)
    return f"""<!DOCTYPE html><html lang="ko"><head><meta charset="UTF-8"/><meta name="viewport" content="width=device-width,initial-scale=1.0"/><title>{html.escape(family)} Dossier 2026-05-09</title>{atlas_style()}</head><body><div class="wrap">
<nav class="topnav"><a href="index.html">Renewed hub</a><a href="family_index_2026_0509.html">Family index</a><a href="figure_atlas_2026_0509.html">Figure atlas</a><a href="table_atlas_2026_0509.html">Table atlas</a><a href="metric_atlas_2026_0509.html">Metric atlas</a><a href="source_atlas_2026_0509.html">Source atlas</a><a href="reviewer_tour_2026_0509.html">Reviewer tour</a></nav>
<header class="hero"><div class="kicker">papers_hub_2026_0509 · family dossier</div><h1>{html.escape(family)}</h1><p class="subtitle">{html.escape(status_label)} · {html.escape(status_boundary)} · {html.escape(status_next)}</p><div class="stats"><div class="stat"><b>{len(family_pages)}</b><span>pages</span></div><div class="stat"><b>{sum(len(p.images) for p in family_pages)}</b><span>figures</span></div><div class="stat"><b>{sum(p.table_count for p in family_pages)}</b><span>tables</span></div><div class="stat"><b>{sum(len(p.headings) for p in family_pages)}</b><span>headings</span></div></div></header>
<section><h2>Pages In This Family</h2><div class="grid">{''.join(page_cards)}</div></section>
<section><h2>Representative Figure Explanations</h2><div class="grid">{''.join(fig_cards) or '<article class="card"><b>No figures</b><p>No image assets detected for this family.</p></article>'}</div></section>
<section><h2>Representative Table Explanations</h2><div class="grid">{''.join(table_cards) or '<article class="card"><b>No tables</b><p>No HTML tables detected for this family.</p></article>'}</div></section>
<section><h2>Source</h2><p><code>family_{html.escape(slug)}_dossier_2026_0509.html</code> generated from renewed hub page inventory. Voice-protected manuscript prose was not edited.</p></section>
</div></body></html>"""


def branch_matches_page(branch: dict[str, str], meta: PageMeta) -> bool:
    # Keep branch matching on stable identifiers. The body text contains
    # generated cross-links/data definitions that would make every branch match
    # every page after renewal injection.
    hay = f"{meta.path.name} {meta.title} {infer_paper_group(meta)}".lower()
    if branch["paper"] == infer_paper_group(meta):
        return True
    raw_tokens = branch.get("tokens", "")
    tokens = raw_tokens.split() if isinstance(raw_tokens, str) else list(raw_tokens)
    tokens = [tok for tok in tokens if len(tok) >= 3]
    if any(tok.lower() in hay for tok in tokens if tok):
        return True
    token_sets = {
        "paper1-main": ["paper1", "fig8", "deconv", "gpl570", "reviewer_defense"],
        "paper1b-braf-axis": ["braf", "tert", "paper1_paper2_braf"],
        "paper2-wsi-image": ["wsi", "image_dm1", "clam", "pathology"],
        "paper2b-ht-immune": ["ht", "hashimoto", "hla", "tls", "paper2_hla"],
        "paper3-ici-readiness": ["paper3", "ici", "checkpoint"],
        "paper4-gd-hla": ["paper4", "gd", "graves", "hla_ncomm"],
        "paper9-ic50-perturbation": ["paper9", "synthetic", "prism", "depmap", "drug", "ic50"],
        "paper11-pancancer": ["paper11", "pancancer", "nature_master"],
        "paper10-atlas": ["paper10", "atlas", "paper10_atlas", "cancer_gene_matrix"],
        "paper12-network": ["paper12", "network", "paper12_network", "module"],
        "paper13-spatial-drug": ["spatial_drug", "target_celltype", "vulnerability", "cv2"],
        "neoantigen-clean-neo": ["neo", "vaccine", "pmhc", "tcr", "barneo", "lumenix"],
    }
    fallback_tokens = [tok for tok in token_sets.get(branch["branch"], [branch["paper"].lower()]) if len(tok) >= 3]
    return any(tok in hay for tok in fallback_tokens)


def branch_pages(branch: dict[str, str], pages: list[PageMeta]) -> list[PageMeta]:
    return [p for p in audit_pages(pages) if branch_matches_page(branch, p)]


def repo_nav() -> str:
    return """<nav class="topnav"><a href="index.html">Renewed hub</a><a href="paper_repo_control_room_2026_05_10.html">Control room</a><a href="high_impact_board_2026_05_10.html">High-impact board</a><a href="impact_escalation_matrix_2026_05_10.html">Impact matrix</a><a href="killer_figure_storyboard_2026_05_10.html">Killer figures</a><a href="project_inventory_2026_05_10.html">Project inventory</a><a href="ic50_prediction_map_2026_05_10.html">IC50/PRISM map</a><a href="paper_repo_map_2026_05_10.html">Repo map</a><a href="paper_branch_manifest_2026_05_10.html">Branch manifest</a><a href="paper_date_topic_index_2026_05_10.html">Date/topic index</a><a href="software_method_ledger_2026_05_10.html">Software ledger</a><a href="deconvolution_repo_2026_05_10.html">Deconvolution repo</a><a href="neoantigen_repo_2026_05_10.html">Neoantigen repo</a></nav>"""


def build_repo_page_card(meta: PageMeta, idx: dict[str, list[str]], limit_sources: int = 6) -> str:
    srcs = " ".join(f'<span class="tag">{html.escape(s[:110])}</span>' for s in source_chips_for_page(meta, idx, limit_sources))
    return (
        '<article class="card">'
        f'<b>{html.escape(infer_paper_group(meta))} · {len(meta.images)} figs · {meta.table_count} tables · score {review_priority_score(meta)}</b>'
        f'<strong><a href="{html.escape(meta.path.name)}">{html.escape(meta.title[:170])}</a></strong>'
        f'<p>{build_metric_strip(meta, 8)}</p><p>{srcs}</p>'
        '</article>'
    )


def numbered_paper_branches() -> list[dict[str, str]]:
    return [b for b in PAPER_BRANCHES if b.get("slot_type") == "paper"]


def side_repo_branches() -> list[dict[str, str]]:
    return [b for b in PAPER_BRANCHES if b.get("slot_type") != "paper"]


def branch_by_id(branch_id: str) -> dict[str, str]:
    return next(b for b in PAPER_BRANCHES if b["branch"] == branch_id)


def branch_tag_link(branch: dict[str, str]) -> str:
    return f'<span class="tag">{html.escape(branch["paper"])} · {html.escape(branch["branch"])}</span>'


def branch_dossier_filename(branch: dict[str, str]) -> str:
    safe = re.sub(r"[^a-zA-Z0-9_-]+", "_", branch["branch"]).strip("_")
    return f"branch_dossier_{safe}_2026_05_10.html"


def branch_dossier_link(branch: dict[str, str]) -> str:
    href = branch_dossier_filename(branch)
    return f'<a href="{html.escape(href)}">{html.escape(href)}</a>'


REPO_LANES = [
    {
        "lane": "01 · 지금 집중",
        "meaning": "Paper 1을 제출 가능한 형태로 잠그는 lane입니다.",
        "branches": ["paper1-main"],
    },
    {
        "lane": "02 · IC50 / actionability map",
        "meaning": "버린 프로젝트가 아니라 Paper 9 중심의 drug-response/IC50 proxy lane입니다.",
        "branches": ["paper9-ic50-perturbation", "paper13-spatial-drug", "paper10-atlas"],
    },
    {
        "lane": "03 · 다음 큰 논문",
        "meaning": "Paper 1 이후 바로 큰 manuscript로 키울 수 있는 branch입니다.",
        "branches": ["paper1b-braf-axis", "paper11-pancancer"],
    },
    {
        "lane": "04 · 보조/응용 논문",
        "meaning": "그림과 응용성은 강하지만 claim boundary를 붙여서 운영해야 하는 branch입니다.",
        "branches": ["paper2-wsi-image", "paper2b-ht-immune", "paper12-network", "paper5-npj-rai8"],
    },
    {
        "lane": "05 · 데이터/실험 대기",
        "meaning": "새 cohort나 thyroid ICI outcome이 오면 급상승하는 branch입니다.",
        "branches": ["paper3-ici-readiness", "paper4-gd-hla"],
    },
    {
        "lane": "06 · 방법/별도 repo",
        "meaning": "생물학 paper가 아니라 method/infra/side repository로 분리해서 관리합니다.",
        "branches": ["paper6-dial", "paper7-agentic-framework", "neoantigen-clean-neo"],
    },
]


PRIMARY_PAGE_BY_BRANCH = {
    "paper1-main": "paper1.html",
    "paper1b-braf-axis": "paper1_paper2_braf_axis_takeover_dossier.html",
    "paper2-wsi-image": "paper2_image_dm1.html",
    "paper2b-ht-immune": "paper2_hla.html",
    "paper3-ici-readiness": "paper3_ici_status_dossier_2026_05_09.html",
    "paper4-gd-hla": "paper4_gd_hla.html",
    "paper5-npj-rai8": "paper5.html",
    "paper6-dial": "paper6.html",
    "paper7-agentic-framework": "paper7.html",
    "paper9-ic50-perturbation": "ic50_prediction_map_2026_05_10.html",
    "paper10-atlas": "paper10_atlas.html",
    "paper11-pancancer": "paper11_pancancer.html",
    "paper12-network": "paper12_network.html",
    "paper13-spatial-drug": "spatial_drug_cv2_overlay_gallery.html",
    "neoantigen-clean-neo": "neoantigen_repo_2026_05_10.html",
}


IMPACT_PROFILES = {
    "paper1-main": {
        "tier": "S",
        "route": "main manuscript / Nature Communications-tier target",
        "killer": "DM1 trunk map + 8-gene axis + Fig 8 mechanism defense + cross-cohort validation.",
        "blocker": "Author-keyboard manuscript completion and causal-language boundary.",
        "next": "Freeze Paper 1 figure order, keep mechanism language correlative, and use branch dossiers only as evidence scaffolding.",
    },
    "paper1b-braf-axis": {
        "tier": "A",
        "route": "merge-or-standalone mechanistic subtype paper",
        "killer": "BRAF/TERT/DM2 escape map with survival, RAI, HM450 and MAPK axis overlays.",
        "blocker": "Must not cannibalize or confuse Paper 1 trunk claim.",
        "next": "Make the branch decision explicit: Paper 1 supplement/reserve versus standalone BRAF-cPTC story.",
    },
    "paper2-wsi-image": {
        "tier": "A-",
        "route": "digital pathology pilot / methods-translational paper",
        "killer": "H&E-to-DM1 projection board with WSI tiles, CLAM/UNI/DINOv2 evidence, and leakage controls.",
        "blocker": "Small N, center/TSS leakage risk, and external slide availability.",
        "next": "Prioritize leakage-free external/held-out slide validation and concise figure-first narrative.",
    },
    "paper2b-ht-immune": {
        "tier": "B+",
        "route": "immune-context biology branch",
        "killer": "HT-overlap HLA-II/IFN/TLS/APC triangulation across bulk, spatial, and scRNA resources.",
        "blocker": "Allele-genetics evidence is weak; expression biology is stronger.",
        "next": "Separate HLA-expression biology from HLA-allele genetics in every figure caption.",
    },
    "paper3-ici-readiness": {
        "tier": "B",
        "route": "future ICI readiness/vulnerability paper",
        "killer": "Immune-readiness map with explicit no-thyroid-ICI-response boundary.",
        "blocker": "No thyroid ICI-treated response cohort.",
        "next": "Keep frozen until response data arrives; use as Paper 11/Paper 1 immune reserve only.",
    },
    "paper4-gd-hla": {
        "tier": "B",
        "route": "GD/Pan-Asian HLA upgrade if Korean data arrives",
        "killer": "Pan-Asian HLA architecture map anchored by DPB1*05:01/B*46:01/A*02:07/C*01:02.",
        "blocker": "No definitive matched Korean adult GD germline cohort yet.",
        "next": "Maintain as data-gated branch; do not oversell literature-meta as discovery.",
    },
    "paper5-npj-rai8": {
        "tier": "B",
        "route": "legacy support / clinical readability layer",
        "killer": "8-gene RAI/TDS clinical biomarker figure supporting Paper 1 readability.",
        "blocker": "Legacy framing can dilute current DM1 dark-matter story.",
        "next": "Use as supporting evidence, not the main claim frame.",
    },
    "paper6-dial": {
        "tier": "C+",
        "route": "methods reserve",
        "killer": "Leakage/theorem audit board showing why model discipline matters.",
        "blocker": "Not biological evidence and can distract from Paper 1.",
        "next": "Park unless a methods submission package is explicitly opened.",
    },
    "paper7-agentic-framework": {
        "tier": "B+",
        "route": "research automation / dossier framework",
        "killer": "End-to-end evidence repo: page inventory, figure/table/data ledgers, branch dossiers, no-abandoned-project audit.",
        "blocker": "Operational framework needs a clean use case and benchmark to be publishable.",
        "next": "Use it as infrastructure now; package later as a methods/demo paper.",
    },
    "paper9-ic50-perturbation": {
        "tier": "A",
        "route": "actionability / wet-lab proposal branch",
        "killer": "PRISM 1518-drug MAPK enrichment + top7 MAPK 7/7 + AZD-0364 FDR 5.88e-7 + SLC1A5/GLUD1/GLS rescue queue.",
        "blocker": "PRISM is proxy drug-response, not direct clinical IC50; wet-lab rescue needed.",
        "next": "Convert SLC1A5 > GLUD1 > GLS into executable perturbation-rescue experiments.",
    },
    "paper10-atlas": {
        "tier": "B+",
        "route": "resource/supplement branch for actionability",
        "killer": "Drug/target/cell-line atlas that makes Paper 9/Paper 11 searchable and reviewer-readable.",
        "blocker": "Atlas co-occurrence is not target validation.",
        "next": "Use as supplement/resource layer, not standalone biological proof.",
    },
    "paper11-pancancer": {
        "tier": "A",
        "route": "next major manuscript after Paper 1",
        "killer": "Pan-cancer portability figure linking DM1 axis to lineage, survival, DepMap, PRISM and hallmark layers.",
        "blocker": "Lineage specificity and residual confounding must be solved up front.",
        "next": "Promote after Paper 1 lock; keep lineage-wise caveat visible.",
    },
    "paper12-network": {
        "tier": "C+",
        "route": "support module for Paper 11",
        "killer": "Network/module architecture around LYN/NAMPT/TACSTD2/KCNN4/PAX8.",
        "blocker": "Bulk network association is not mechanism.",
        "next": "Absorb into Paper 11 as a mechanistic support panel.",
    },
    "paper13-spatial-drug": {
        "tier": "A-",
        "route": "visual actionability/spatial validation branch",
        "killer": "Spatial vulnerability overlays showing drug-target/state co-localization and mismatch boards.",
        "blocker": "Spatial co-localization is not drug efficacy; quantitative controls need locking.",
        "next": "Turn best overlays into a match/mismatch figure with explicit negative controls.",
    },
    "neoantigen-clean-neo": {
        "tier": "A",
        "route": "separate methods/resource universe",
        "killer": "Leakage-aware benchmark + abstention/interpretable neoantigen ranking + TCR/structure/MD audit extensions.",
        "blocker": "Do not mix with THCA Paper 1-4 claims; source-heldout performance remains heterogeneous.",
        "next": "Keep as separate repo/product/methods track with strict split and overlap audit front-loaded.",
    },
}


TIER_ORDER = {"S": 0, "A": 1, "A-": 2, "B+": 3, "B": 4, "C+": 5, "C": 6}


KILLER_FIGURE_ORDER = [
    "paper1-main",
    "paper9-ic50-perturbation",
    "paper11-pancancer",
    "paper2-wsi-image",
    "paper13-spatial-drug",
    "neoantigen-clean-neo",
    "paper1b-braf-axis",
]


def impact_profile(branch: dict[str, str]) -> dict[str, str]:
    return IMPACT_PROFILES.get(
        branch["branch"],
        {
            "tier": "B",
            "route": branch["state"],
            "killer": branch["claim"],
            "blocker": branch["boundary"],
            "next": branch["meaning"],
        },
    )


def primary_page_link(branch: dict[str, str]) -> str:
    href = PRIMARY_PAGE_BY_BRANCH.get(branch["branch"], "paper_branch_manifest_2026_05_10.html")
    return f'<a href="{html.escape(href)}">{html.escape(href)}</a>'


def killer_panels_for_branch(branch: dict[str, str]) -> list[tuple[str, str, str, str, str, str]]:
    """Panel label, title, data, software/method, use, boundary."""
    bid = branch["branch"]
    if bid == "paper1-main":
        return [
            ("A", "DM1 trunk map", "TCGA-THCA DM labels, 8-gene score, driver/clinical strata", "Python module scoring and effect-size scripts", "First-page identity panel for Paper 1", "Mechanism language remains correlative."),
            ("B", "Cross-cohort validation strip", "K2/PRJEB11591, Lee/GSE213647, GSE76039, GPL570/Korean validation rows", "forest/effect-size summaries", "Shows portability beyond TCGA", "Each cohort keeps its own assay/filter boundary."),
            ("C", "Fig 8 mechanism defense", "HM450, thyroid-lineage genes, methylation/expression support", "methylation/effect-size audits", "Supports lineage-silencing layer", "Motivates rather than proves causality."),
            ("D", "Deconvolution residual map", "Lu 2023 scRNA pseudobulk, TCGA bulk, nu-SVR residualization", "NuSVR plus NNLS-family audits", "Separates composition from retained molecular axis", "Cell-type adjustment is sensitivity analysis."),
            ("E", "Spatial caveat board", "GSE250521/GSE193581 spatial checks and v16/v17 caveats", "spatial scoring and residualization", "Shows what validates and what does not", "Negative rows stay visible."),
            ("F", "Reviewer defense ledger", "claim ledger, figure/table matrix, source paths", "HTML/TSV/JSON dossier builder", "Makes evidence reviewer-readable", "No protected manuscript prose generated."),
        ]
    if bid == "paper9-ic50-perturbation":
        return [
            ("A", "PRISM MAPK response map", "PRISM 1518-drug response matrix and MAPK annotations", "Spearman/FDR/Cohen's d", "Main actionability signal", "Proxy drug response, not clinical IC50."),
            ("B", "Top7 MAPK 7/7 enrichment", "synthetic_lethality_summary.json MAPK top-rank counts", "drug-class enrichment", "Shows pharmacologic coherence", "Prioritization, not efficacy."),
            ("C", "AZD-0364 anchor", "AZD delta LFC, d=-0.594, FDR=5.88e-7", "effect-size audit", "Named wet-lab anchor", "Needs dose-response validation."),
            ("D", "DepMap genetic overlay", "MYC/NAMPT/SLC1A5/GLUD1/GLS dependency rows", "CRISPR dependency effect sizes", "Separates genetic candidates from drug screen", "Pan-essentiality controls required."),
            ("E", "SLC1A5 > GLUD1 > GLS queue", "paper9_wetlab_perturbation_matrix.tsv", "validation design table", "Executable perturbation/rescue plan", "No synthetic-lethal claim until rescue passes."),
            ("F", "Spatial mismatch/caveat panel", "Paper 13 overlays and v16/v17 negative spatial rows", "spatial target-state overlap", "Keeps visual actionability honest", "Co-localization is not efficacy."),
        ]
    if bid == "paper11-pancancer":
        return [
            ("A", "Pan-cancer portability map", "TCGA pan-cancer expression and DM1 score", "lineage-wise matrix scoring", "Shows whether thyroid-derived axis travels", "Lineage effects front-loaded."),
            ("B", "Lineage mitigation panel", "lineage residual/confounding rows", "stratified models and residual checks", "Prevents overbroad claim", "Do not pool away lineage specificity."),
            ("C", "Survival/meta panel", "survival endpoints and hazard/effect-size summaries", "Cox/meta-analysis tables", "Clinical relevance layer", "Endpoint/event counts govern claim."),
            ("D", "Dependency/drug overlay", "DepMap and PRISM overlays", "correlation/effect-size scripts", "Links DM1 state to actionability", "Screen data is prioritization."),
            ("E", "Hallmark/pathway panel", "hallmark/pathway module overlays", "module scoring", "Biological readability", "Pathway correlation is not causality."),
            ("F", "Promotion gate", "Paper 11 dossier, source ledger, lineage caveat", "HTML dossier builder", "Next-major-manuscript decision", "Promote after Paper 1 lock."),
        ]
    if bid == "paper2-wsi-image":
        return [
            ("A", "WSI tile atlas", "TCGA-THCA WSI and embedded pilot slides", "OpenSlide-style tiling", "Visual entry point for image-DM1", "No deployment claim."),
            ("B", "Foundation embedding panel", "DINOv2/UNI/CLAM features", "embedding extraction and sklearn models", "Shows signal source", "Small N and center leakage controlled."),
            ("C", "LOSO validation panel", "held-out slide/site rows", "LOSO models and AUC/effect sizes", "Core leakage-control proof", "Avoid site/TSS inflation."),
            ("D", "Path2Space bridge", "GSE230424/GSE250521/GSE248205 pathology/spatial resources", "spatial-pathology bridge scripts", "Connects morphology to molecular state", "Pilot evidence only."),
            ("E", "Failure controls", "audit tables and mismatch cases", "QC dashboards", "Prevents overfit story", "Negative panels stay co-located."),
            ("F", "Clinical boundary", "branch dossier/source ledger", "HTML audit builder", "States what can be claimed", "No clinical pathology model yet."),
        ]
    if bid == "paper13-spatial-drug":
        return [
            ("A", "Drug overlay atlas", "GSE250521/GSE230424 spatial overlays and famous drug panels", "spatial module scoring", "Visual actionability hook", "Overlay is not efficacy."),
            ("B", "Target-state overlap", "target genes, DM1/TDS/RAI state modules", "target-state overlap scoring", "Shows alignment of target and state", "Needs quantitative controls."),
            ("C", "Hotspot zoom strips", "ROI/hotspot images and contours", "CV2 ROI gallery generation", "Makes spatial story inspectable", "ROI choice transparent."),
            ("D", "Match/mismatch board", "positive overlays and negative cases", "match/mismatch decision tables", "Reviewer-honest control", "Mismatch cases not hidden."),
            ("E", "Quant control panel", "stage/hotspot burden and detection controls", "residualization/QC tables", "Locks credibility", "Must pass before escalation."),
            ("F", "Wet-lab translation", "candidate targets and tissue pockets", "ROI validation planning", "Links to IHC/RNAscope/perturbation", "Proposal until tested."),
        ]
    if bid == "neoantigen-clean-neo":
        return [
            ("A", "Leakage contract", "ITSNdb/TESLA/CEDAR/dbPepNeo2/NeoDB overlap audits", "strict split/public-overlap checks", "Defines benchmark credibility", "Never mix leaked and strict rows."),
            ("B", "Model comparison", "CLEAN-Neo/CROSS-Neo and competitor rows", "benchmark harness", "Ranks against baselines", "Source-heldout heterogeneity remains."),
            ("C", "Abstention/calibration", "BMA queues and uncertainty rows", "calibration/abstention analysis", "Usability under uncertainty", "Abstention is not biology."),
            ("D", "Structure/TCR layer", "ESM2/structure proxy, TCR/pMHC/MD outputs", "structure proxy/MD pipeline", "Mechanistic interpretability", "Proxy is not experimental binding."),
            ("E", "Patient triage", "PAAD/THCA vaccine triage evidence", "patient metadata gating", "Separates plausible use cases", "No broad thyroid vaccine claim."),
            ("F", "Separate repo boundary", "Neo dossier and source paths", "HTML/source ledger builder", "Keeps methods universe separate", "Do not merge with Paper 1-4 claims."),
        ]
    if bid == "paper1b-braf-axis":
        return [
            ("A", "BRAF/TERT subtype map", "TCGA BRAF/RAS/TERT/clinical strata", "Fisher/effect-size scripts", "Defines subtype fork", "Do not generalize to all DM1."),
            ("B", "DM2/escape state", "Paper 1 DM labels plus BRAF axis rows", "module scoring", "Shows trunk relation", "Avoid cannibalizing Paper 1."),
            ("C", "Survival/RAI panel", "clinical and RAI-response rows", "survival/effect-size tables", "Clinical readability", "Endpoint counts govern claim."),
            ("D", "HM450/MAPK panel", "HM450, MAPK output, lineage genes", "methylation/pathway scoring", "Mechanistic subtype readability", "Correlative only."),
            ("E", "External support", "Lee 2024 and spatial BRAF inference", "external validation scripts", "External support layer", "Inference is not genotyping."),
            ("F", "Merge gate", "takeover dossier and source ledger", "HTML synthesis", "Manuscript routing decision", "Do not slow Paper 1."),
        ]
    profile = impact_profile(branch)
    return [
        ("A", "Branch identity", branch["data"], branch["software"], profile["route"], branch["boundary"]),
        ("B", "Killer figure", profile["killer"], branch["software"], profile["next"], profile["blocker"]),
        ("C", "Source ledger", branch["paths"], "HTML/TSV/JSON inventories", "Trace every claim to a local file", "Source files remain binding."),
    ]


RESULT_DIR_BRANCH_RULES = [
    ("paper1-main", ["00_qc", "01_spatial_score", "02_stage_trend", "bulk_rnaseq", "dark_matter", "dm1_robustness", "external_expression", "gpl570", "landa", "ncomm_push", "p1_driver", "p_deconv", "v17", "v5", "v6_scrna", "v8p1"]),
    ("paper1b-braf-axis", ["braf_nature", "braf_multimodal", "tert_recovery"]),
    ("paper2-wsi-image", ["03_pathology", "image_dm1", "p2_pillar", "p2_power", "pathology", "spatial_full", "wsi_dm1"]),
    ("paper2b-ht-immune", ["gse286332", "hashimoto", "hla_deepdive", "hla_two_paper", "ht_biology", "htptc", "p5_8gene_vs_hla", "paper2_ht"]),
    ("paper3-ici-readiness", ["paper3_ici", "track18_ici", "phase_c_ici", "p3_p9_full_execution"]),
    ("paper4-gd-hla", ["d4p1", "gd_hla", "graves", "hla_deepdive", "paper4_gd"]),
    ("paper5-npj-rai8", ["d7p3", "d8b", "d8c", "v17p", "ultimate", "unified_model"]),
    ("paper6-dial", ["dial", "v10_aaai", "v15_neurips"]),
    ("paper7-agentic-framework", ["agentic", "paper_portfolio", "portfolio", "pantheonos", "transfer_packages"]),
    ("paper9-ic50-perturbation", ["drug_discovery", "drug_phaseb", "p_synthetic_lethality", "paper9", "perturbation", "repurposing", "sl_first_pass", "synthetic_lethality", "v13_drug", "v14_drug", "v7_repurposing"]),
    ("paper10-atlas", ["paper10_atlas", "cancer_gene_matrix"]),
    ("paper11-pancancer", ["pancancer", "paper11"]),
    ("paper12-network", ["network", "paper12"]),
    ("paper13-spatial-drug", ["spatial_drug", "target_celltype", "vulnerability", "cv2", "spatial_cnv", "spatial_proteomics"]),
    ("neoantigen-clean-neo", ["barneo", "cancer_vaccine", "clean_neobench", "cross_neo", "neo_bayesian", "neo_product", "neo_review", "p_cancer_vaccine", "vaccine"]),
]

FOUNDATION_RESULT_NEEDLES = [
    "audit_",
    "figs",
    "figures",
    "figures_for_advisor",
    "geo_search",
    "high_impact_topic_pilots",
    "json",
    "methylation",
    "microarray",
    "ml",
    "proteogenomic",
    "tables",
    "terminology",
]


def owner_tags(owner_ids: list[str]) -> str:
    if not owner_ids:
        return '<span class="tag">operations-parking-lot</span>'
    tags = []
    for owner_id in owner_ids:
        b = branch_by_id(owner_id)
        tags.append(f'<span class="tag">{html.escape(b["paper"])} · {html.escape(owner_id)}</span>')
    return " ".join(tags)


def result_dir_inventory() -> list[dict[str, object]]:
    base = ROOT / "project" / "results"
    if not base.exists():
        return []
    rows: list[dict[str, object]] = []
    for d in sorted((p for p in base.iterdir() if p.is_dir()), key=lambda p: p.name.lower()):
        low = d.name.lower()
        owners: list[str] = []
        for branch_id, needles in RESULT_DIR_BRANCH_RULES:
            if any(n in low for n in needles):
                owners.append(branch_id)
        owners = list(dict.fromkeys(owners))
        if owners:
            disposition = "OWNED_BRANCH"
            policy = "active_or_supporting_branch"
        elif any(n in low for n in FOUNDATION_RESULT_NEEDLES):
            disposition = "FOUNDATION_ARCHIVE"
            policy = "retained_support_archive"
        else:
            disposition = "TRIAGE_PARKING_LOT"
            policy = "not_abandoned_review_queue"
        try:
            child_dirs = sum(1 for c in d.iterdir() if c.is_dir())
            child_files = sum(1 for c in d.iterdir() if c.is_file())
        except OSError:
            child_dirs = 0
            child_files = 0
        summary = ""
        for candidate in ["SUMMARY.md", "summary.json", "README.md", "MANIFEST.json", "phase2_summary.json"]:
            if (d / candidate).exists():
                summary = candidate
                break
        rows.append(
            {
                "path": f"project/results/{d.name}",
                "name": d.name,
                "owners": owners,
                "disposition": disposition,
                "policy": policy,
                "child_dirs": child_dirs,
                "child_files": child_files,
                "summary": summary,
            }
        )
    return rows


def write_project_inventory_files() -> None:
    rows = result_dir_inventory()
    serializable_rows = [{**row, "owners": ";".join(row["owners"])} for row in rows]
    (TARGET / "project_inventory_2026_05_10.json").write_text(
        json.dumps(serializable_rows, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    cols = ["path", "owners", "disposition", "policy", "child_dirs", "child_files", "summary"]
    lines = ["\t".join(cols)]
    for row in serializable_rows:
        lines.append("\t".join(str(row.get(c, "")).replace("\t", " ").replace("\n", " ") for c in cols))
    (TARGET / "project_inventory_2026_05_10.tsv").write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_project_inventory_page(pages: list[PageMeta], idx: dict[str, list[str]]) -> str:
    del idx
    rows = result_dir_inventory()
    active = [r for r in rows if r["disposition"] == "OWNED_BRANCH"]
    archives = [r for r in rows if r["disposition"] == "FOUNDATION_ARCHIVE"]
    parking = [r for r in rows if r["disposition"] == "TRIAGE_PARKING_LOT"]
    ic50_owner_ids = {"paper9-ic50-perturbation", "paper10-atlas", "paper13-spatial-drug"}
    ic50_dirs = [r for r in rows if ic50_owner_ids.intersection(set(r["owners"]))]
    ic50_pages = branch_pages(branch_by_id("paper9-ic50-perturbation"), pages)

    ic50_dir_rows = "".join(
        "<tr>"
        f"<td>{html.escape(str(r['path']))}</td>"
        f"<td>{owner_tags(list(r['owners']))}</td>"
        f"<td>{html.escape(str(r['disposition']))}</td>"
        f"<td>{html.escape(str(r['summary'] or ''))}</td>"
        "</tr>"
        for r in ic50_dirs
    )
    all_rows = "".join(
        "<tr>"
        f"<td><code>{html.escape(str(r['path']))}</code></td>"
        f"<td>{owner_tags(list(r['owners']))}</td>"
        f"<td>{html.escape(str(r['disposition']))}</td>"
        f"<td>{html.escape(str(r['policy']))}</td>"
        f"<td>{r['child_dirs']}</td><td>{r['child_files']}</td>"
        f"<td>{html.escape(str(r['summary'] or ''))}</td>"
        "</tr>"
        for r in rows
    )
    parking_rows = "".join(
        hrow([r["path"], r["disposition"], r["policy"], f"{r['child_dirs']} dirs / {r['child_files']} files", r["summary"] or ""])
        for r in parking[:80]
    )
    top_ic50_pages = "".join(
        build_repo_page_card(p, {}) for p in sorted(ic50_pages, key=lambda p: -review_priority_score(p))[:12]
    )

    return f"""<!DOCTYPE html><html lang="ko"><head><meta charset="UTF-8"/><meta name="viewport" content="width=device-width,initial-scale=1.0"/><title>Project Inventory 2026-05-10</title>{atlas_style()}</head><body><div class="wrap">
{repo_nav()}
<header class="hero"><div class="kicker">no abandoned project audit · 2026-05-10</div><h1>Project Inventory / Parking Lot</h1><p class="subtitle">IC50/PRISM map은 Paper 9 active branch입니다. 나머지 결과 디렉터리도 owned branch, foundation archive, triage parking lot으로 보존해서 버려진 프로젝트가 없게 관리합니다.</p><div class="stats"><div class="stat"><b>{len(rows)}</b><span>result dirs scanned</span></div><div class="stat"><b>{len(active)}</b><span>owned branch dirs</span></div><div class="stat"><b>{len(archives)}</b><span>support archives</span></div><div class="stat"><b>{len(parking)}</b><span>triage parking lot</span></div><div class="stat"><b>0</b><span>abandoned by policy</span></div><div class="stat"><b>{len(ic50_dirs)}</b><span>IC50/actionability dirs</span></div></div></header>
<section><h2>IC50 Is Active</h2><div class="warning"><b>Owner:</b> Paper 9 · paper9-ic50-perturbation. <b>Primary page:</b> <a href="ic50_prediction_map_2026_05_10.html">ic50_prediction_map_2026_05_10.html</a>. PRISM/DepMap은 direct clinical IC50이 아니라 public drug-response/dependency proxy입니다.</div><table>{hrow(["Result dir","Owner","Disposition","Summary file"], True)}{ic50_dir_rows}</table></section>
<section><h2>IC50 Evidence Pages</h2><div class="grid">{top_ic50_pages}</div></section>
<section><h2>Parking-Lot Rule</h2><table>{hrow(["Bucket","Meaning","Action"], True)}{hrow(["OWNED_BRANCH","이미 14 paper slot 또는 side repo에 연결됨","branch primary page에서 관리"])}{hrow(["FOUNDATION_ARCHIVE","figure/table/audit/raw support 성격","삭제하지 않고 evidence archive로 보존"])}{hrow(["TRIAGE_PARKING_LOT","branch owner가 아직 약하거나 여러 paper에 걸침","다음 audit 때 owner를 승격하거나 archive로 확정"])}</table></section>
<section><h2>Triage Parking Lot</h2><table>{hrow(["Path","Disposition","Policy","Immediate contents","Summary"], True)}{parking_rows or hrow(["none","none","none","none","none"])}</table></section>
<section><h2>All Result Directories</h2><table>{hrow(["Path","Owner","Disposition","Policy","Dirs","Files","Summary"], True)}{all_rows}</table></section>
<section><h2>Machine Files</h2><p><a href="project_inventory_2026_05_10.tsv">project_inventory_2026_05_10.tsv</a> · <a href="project_inventory_2026_05_10.json">project_inventory_2026_05_10.json</a></p></section>
</div></body></html>"""


def branch_representative_figures(branch: dict[str, str], ps: list[PageMeta], idx: dict[str, list[str]], limit: int = 12) -> str:
    cards: list[str] = []
    for page in sorted(ps, key=lambda p: -review_priority_score(p)):
        for i, img in enumerate(page.images, 1):
            if len(cards) >= limit:
                return "".join(cards)
            caption = page.captions[i - 1] if i - 1 < len(page.captions) else img.get("alt", "")
            caption = normalize(caption) or "No caption/alt detected."
            role, how, caveat = visual_reading_notes(page, img, caption)
            source_matches = "; ".join(asset_matches(img.get("src", ""), idx)[:3]) or guess_visual_data(page, img)
            src = html.escape(img.get("src", ""), quote=True)
            cards.append(
                '<article class="card figcard">'
                f'<a class="thumb" href="{html.escape(page.path.name)}" style="background-image:url({src})"></a>'
                '<div>'
                f'<b>{html.escape(branch["paper"])} figure · {html.escape(role)}</b>'
                f'<strong><a href="{html.escape(page.path.name)}">{html.escape(Path(img.get("src", "")).name or page.path.name)}</a></strong>'
                f'<p>{html.escape(caption[:430])}</p>'
                f'<dl class="dl"><dt>Data</dt><dd>{html.escape(source_matches)}</dd><dt>Read</dt><dd>{html.escape(how)}</dd><dt>Boundary</dt><dd>{html.escape(caveat)}</dd></dl>'
                '</div></article>'
            )
    return "".join(cards) or '<article class="card"><b>No detected figures</b><p>This branch is tracked by tables, source ledgers, and dossier text, but no image assets were detected in matched pages.</p></article>'


def branch_representative_tables(branch: dict[str, str], ps: list[PageMeta], limit: int = 12) -> str:
    del branch
    cards: list[str] = []
    for page in sorted(ps, key=lambda p: (-(p.table_count + len(p.images)), p.path.name)):
        for i in range(1, min(page.table_count, 6) + 1):
            if len(cards) >= limit:
                return "".join(cards)
            caption = page.table_captions[i - 1] if i - 1 < len(page.table_captions) else f"Inline HTML table {i}"
            cards.append(
                '<article class="card">'
                f'<b>Table {i:02d} · {html.escape(page.path.name)}</b>'
                f'<strong><a href="{html.escape(page.path.name)}">{html.escape(normalize(caption)[:240])}</a></strong>'
                '<p><b>Data entered:</b> Embedded HTML table body and any linked TSV/JSON/MD on the page.</p>'
                '<p><b>How to read:</b> Treat rows as the immediate evidence matrix for the page claim; check directionality and cohort labels before reuse.</p>'
                '<p><b>Boundary:</b> Table presence is provenance, not automatic causal proof.</p>'
                '</article>'
            )
    return "".join(cards) or '<article class="card"><b>No detected tables</b><p>No embedded HTML tables were found in matched pages; use source path and figure/data ledgers for provenance.</p></article>'


def branch_overview_figures(branch: dict[str, str], ps: list[PageMeta]) -> str:
    pages_n = len(ps)
    figs_n = sum(len(p.images) for p in ps)
    tables_n = sum(p.table_count for p in ps)
    rows = [
        ("Overview Fig 1", "Branch repo map", f"{branch['paper']} registry row, parent branch, opened/last_update, state, {pages_n} matched pages", "Explains where this paper lives in the portfolio and what history it inherits.", "Organization figure only; evidence still lives in source pages."),
        ("Overview Fig 2", "Data provenance map", branch["data"], "Shows which cohorts, public screens, scores, or resources feed this branch.", "Dataset availability does not imply every claim is validated."),
        ("Overview Fig 3", "Software/method map", branch["software"], "Shows which computational stack produced the branch outputs.", "Software describes execution, not biological causality."),
        ("Overview Fig 4", "Claim and boundary map", branch["claim"], "States the strongest allowed claim and prevents overreach.", branch["boundary"]),
        ("Overview Fig 5", "Evidence density map", f"{pages_n} pages, {figs_n} figures, {tables_n} tables", "Shows whether a branch is figure-heavy, table-heavy, or still thin.", "Counts are audit metrics, not journal readiness by themselves."),
        ("Overview Fig 6", "Next-use map", branch["meaning"], "Shows whether to submit now, reserve, merge, park, or validate experimentally.", "Strategic disposition can change after new data."),
    ]
    return "".join(
        '<article class="card">'
        f'<b>{html.escape(label)}</b><strong>{html.escape(title)}</strong>'
        f'<p><b>Data entered:</b> {html.escape(data)}</p>'
        f'<p><b>What it explains:</b> {html.escape(read)}</p>'
        f'<p><b>Boundary:</b> {html.escape(boundary)}</p>'
        '</article>'
        for label, title, data, read, boundary in rows
    )


def branch_flow_cards(branch: dict[str, str]) -> str:
    steps = [
        ("Data intake", branch["data"], "Make sure every cohort/resource is named before reading claims."),
        ("Score/model layer", branch["software"], "Identify whether this is module scoring, survival, deconvolution, WSI, PRISM, HLA, or neoantigen modeling."),
        ("Evidence pages", branch["paths"], "Use linked pages and result paths as the source of truth."),
        ("Representative visuals", "Figure/table audit cards", "Read what data entered each visual and what the visual can/cannot support."),
        ("Claim", branch["claim"], "Use the strongest allowed claim only inside the branch boundary."),
        ("Disposition", f"{branch['state']} · {branch['meaning']}", "Decide submit, merge, reserve, validate, or park."),
    ]
    return "".join(
        '<article class="step">'
        f'<strong>{html.escape(title)}</strong>'
        f'<p>{html.escape(text)}</p>'
        f'<p><span class="tag">{html.escape(use[:120])}</span></p>'
        '</article>'
        for title, text, use in steps
    )


def build_branch_dossier(branch: dict[str, str], pages: list[PageMeta], idx: dict[str, list[str]]) -> str:
    ps = branch_pages(branch, pages)
    profile = impact_profile(branch)
    page_cards = "".join(build_repo_page_card(p, idx) for p in sorted(ps, key=lambda p: -review_priority_score(p))[:30])
    term_counts: dict[str, int] = {}
    for page in ps:
        for term in infer_terms(page):
            term_counts[term["key"]] = term_counts.get(term["key"], 0) + 1
    term_rows = "".join(
        hrow([term, n, "Detected from matched page text/title/assets", "Use only after checking page-specific source paths"])
        for term, n in sorted(term_counts.items(), key=lambda x: (-x[1], x[0]))[:14]
    ) or hrow(["Branch registry data", len(ps), branch["data"], branch["boundary"]])
    source_rows = "".join(
        hrow([p.path.name, infer_paper_group(p), len(p.images), p.table_count, "; ".join(source_chips_for_page(p, idx, 4)) or "page-local content"])
        for p in sorted(ps, key=lambda p: -review_priority_score(p))[:20]
    ) or hrow(["No matched page", branch["paper"], 0, 0, branch["paths"]])
    table_explanation_rows = "".join(hrow(r) for r in [
        ["Registry table", "branch title/topic/state/parent/date", "Portfolio ownership and history", "Not biological evidence."],
        ["Data definition table", "branch data + inferred page terms", "Defines what data entered the branch", "Dataset names must be source-checked before manuscript use."],
        ["Representative table cards", "Embedded HTML tables from matched pages", "Shows result rows, metrics, and claim support", "Directionality/covariate labels must be read per page."],
        ["Source path table", "Top matched pages plus linked files", "Tells exactly where to click for provenance", "Generated inventory can lag if source files change before rebuild."],
    ])
    return f"""<!DOCTYPE html><html lang="ko"><head><meta charset="UTF-8"/><meta name="viewport" content="width=device-width,initial-scale=1.0"/><title>{html.escape(branch['paper'])} Branch Dossier 2026-05-10</title>{atlas_style()}</head><body><div class="wrap">
{repo_nav()}
<header class="hero"><div class="kicker">high-impact branch dossier · {html.escape(branch['branch'])}</div><h1>{html.escape(branch['paper'])} · {html.escape(branch['title'])}</h1><p class="subtitle">{html.escape(branch['meaning'])}</p><div class="stats"><div class="stat"><b>{html.escape(branch['state'])}</b><span>state</span></div><div class="stat"><b>{len(ps)}</b><span>matched pages</span></div><div class="stat"><b>{sum(len(p.images) for p in ps)}</b><span>figures explained</span></div><div class="stat"><b>{sum(p.table_count for p in ps)}</b><span>tables explained</span></div><div class="stat"><b>{html.escape(branch['last_update'])}</b><span>last update</span></div><div class="stat"><b>{html.escape(branch['parent'])}</b><span>parent</span></div></div></header>
<section><h2>Data Definition First</h2><div class="warning"><b>Claim:</b> {html.escape(branch['claim'])}<br/><b>Boundary:</b> {html.escape(branch['boundary'])}</div><table>{hrow(["Data/software block","Definition or source","Where it is used","Interpretation rule"], True)}{hrow(["Data", branch["data"], "Branch evidence and linked pages", branch["boundary"]])}{hrow(["Software", branch["software"], "Analysis execution and visualization", "Software output is not causal proof by itself."])}{hrow(["Source paths", branch["paths"], "Primary local provenance", "Use linked files as source of truth."])}</table></section>
<section><h2>Impact Profile</h2><table>{hrow(["Impact tier","Target route","Killer figure","Current blocker","Next 48h action"], True)}{hrow([profile["tier"], profile["route"], profile["killer"], profile["blocker"], profile["next"]])}</table></section>
<section><h2>Whole-Branch Flow</h2><div class="flow">{branch_flow_cards(branch)}</div></section>
<section><h2>Representative Overview Figures</h2><div class="grid">{branch_overview_figures(branch, ps)}</div></section>
<section><h2>Representative Evidence Figures</h2><div class="grid">{branch_representative_figures(branch, ps, idx)}</div></section>
<section><h2>Representative Table Explanations</h2><div class="grid">{branch_representative_tables(branch, ps)}</div></section>
<section><h2>Detected Data Terms</h2><table>{hrow(["Term","Matched pages","Definition source","Interpretation rule"], True)}{term_rows}</table></section>
<section><h2>Table Logic Overview</h2><table>{hrow(["Table type","Data entered","Where used","Boundary"], True)}{table_explanation_rows}</table></section>
<section><h2>Source Page Ledger</h2><table>{hrow(["Page","Family","Figures","Tables","Source chips"], True)}{source_rows}</table></section>
<section><h2>Tracked Evidence Pages</h2><div class="grid">{page_cards or '<article class="card"><b>No matched source pages</b><p>Use branch registry/source paths until a dedicated page exists.</p></article>'}</div></section>
</div></body></html>"""


def build_high_impact_board(pages: list[PageMeta], idx: dict[str, list[str]]) -> str:
    source_pages = audit_pages(pages)
    cards = []
    rows = []
    for branch in PAPER_BRANCHES:
        ps = branch_pages(branch, source_pages)
        score = sum(review_priority_score(p) for p in ps)
        rows.append((branch, ps, score))
    rows.sort(key=lambda item: (TIER_ORDER.get(impact_profile(item[0])["tier"], 9), -item[2], item[0]["paper"]))
    tier_rows = []
    killer_rows = []
    for branch, ps, score in rows:
        primary = PRIMARY_PAGE_BY_BRANCH.get(branch["branch"], "paper_branch_manifest_2026_05_10.html")
        profile = impact_profile(branch)
        tier_rows.append(
            "<tr>"
            f"<td><span class=\"tag\">{html.escape(profile['tier'])}</span></td>"
            f"<td>{html.escape(branch['paper'])}</td>"
            f"<td>{html.escape(profile['route'])}</td>"
            f"<td>{len(ps)}</td><td>{sum(len(p.images) for p in ps)}</td><td>{sum(p.table_count for p in ps)}</td>"
            f"<td>{html.escape(profile['blocker'])}</td>"
            f"<td><a href=\"{html.escape(branch_dossier_filename(branch))}\">dossier</a></td>"
            "</tr>"
        )
        killer_rows.append(
            "<tr>"
            f"<td>{html.escape(branch['paper'])}</td>"
            f"<td>{html.escape(profile['killer'])}</td>"
            f"<td>{html.escape(branch['data'])}</td>"
            f"<td>{html.escape(profile['next'])}</td>"
            f"<td>{html.escape(branch['boundary'])}</td>"
            "</tr>"
        )
        cards.append(
            '<article class="card">'
            f'<b>{html.escape(branch["paper"])} · tier {html.escape(profile["tier"])} · score {score} · {html.escape(branch["state"])}</b>'
            f'<strong>{html.escape(branch["title"])}</strong>'
            f'<p><span class="tag">{len(ps)} pages</span><span class="tag">{sum(len(p.images) for p in ps)} figs</span><span class="tag">{sum(p.table_count for p in ps)} tables</span><span class="tag">parent: {html.escape(branch["parent"])}</span></p>'
            f'<p><b>Target route:</b> {html.escape(profile["route"])}</p>'
            f'<p><b>Killer figure:</b> {html.escape(profile["killer"])}</p>'
            f'<p><b>Blocker:</b> {html.escape(profile["blocker"])}</p>'
            f'<p><b>Next 48h:</b> {html.escape(profile["next"])}</p>'
            f'<p><a href="{html.escape(branch_dossier_filename(branch))}">branch dossier</a> · <a href="{html.escape(primary)}">primary page</a></p>'
            '</article>'
        )
    lane_cards = []
    for lane in REPO_LANES:
        links = " ".join(
            f'<a class="tag" href="{html.escape(branch_dossier_filename(branch_by_id(branch_id)))}">{html.escape(branch_by_id(branch_id)["paper"])}</a>'
            for branch_id in lane["branches"]
        )
        lane_cards.append(
            '<article class="step">'
            f'<strong>{html.escape(lane["lane"])}</strong>'
            f'<p>{html.escape(lane["meaning"])}</p><p>{links}</p>'
            '</article>'
        )
    source_rows = "".join(hrow(r) for r in [
        ["Branch dossiers", "14 numbered paper slots + Neo side repo", "Each has data definition, flow, representative figures/tables, source ledger", "Generated review surface; source files remain binding."],
        ["IC50 map", "PRISM/DepMap/Paper 9 candidate matrices", "Actionability branch with decision board and wet-lab queue", "Proxy until dose-response validation."],
        ["Project inventory", "project/results top-level directory audit", "No-abandoned-project owner/parking-lot ledger", "Parking lot needs future owner decision."],
        ["Control room", "PAPER_BRANCHES registry", "Top-level repo-style map", "Does not edit protected manuscript prose."],
    ])
    return f"""<!DOCTYPE html><html lang="ko"><head><meta charset="UTF-8"/><meta name="viewport" content="width=device-width,initial-scale=1.0"/><title>High-Impact Board 2026-05-10</title>{atlas_style()}</head><body><div class="wrap">
{repo_nav()}
<header class="hero"><div class="kicker">portfolio escalation board · 2026-05-10</div><h1>High-Impact Board</h1><p class="subtitle">가능성에서 멈추지 않도록 모든 branch를 dossier 단위로 승격한 board입니다. 각 paper는 데이터 정의서, 전체 흐름, 대표 figure/table 설명, claim boundary, source ledger를 가집니다.</p><div class="stats"><div class="stat"><b>{len(numbered_paper_branches())}</b><span>paper dossiers</span></div><div class="stat"><b>{len(side_repo_branches())}</b><span>side repos</span></div><div class="stat"><b>{len(source_pages)}</b><span>source pages</span></div><div class="stat"><b>{sum(len(p.images) for p in source_pages)}</b><span>figures audited</span></div><div class="stat"><b>{sum(p.table_count for p in source_pages)}</b><span>tables audited</span></div><div class="stat"><b>0</b><span>abandoned by policy</span></div></div></header>
<section><h2>Escalation Flow</h2><div class="flow">{''.join(lane_cards)}</div></section>
<section><h2>Impact Tier Matrix</h2><table>{hrow(["Tier","Branch","Target route","Pages","Figures","Tables","Blocker","Dossier"], True)}{''.join(tier_rows)}</table></section>
<section><h2>Killer Figure Roadmap</h2><table>{hrow(["Branch","Killer figure","Data entered","Next 48h action","Boundary"], True)}{''.join(killer_rows)}</table></section>
<section><h2>Source Contract</h2><table>{hrow(["Layer","Data entered","Where used","Boundary"], True)}{source_rows}</table></section>
<section><h2>Branch Dossiers</h2><div class="grid">{''.join(cards)}</div></section>
</div></body></html>"""


def build_impact_escalation_matrix(pages: list[PageMeta], idx: dict[str, list[str]]) -> str:
    del idx
    source_pages = audit_pages(pages)
    rows = []
    for branch in PAPER_BRANCHES:
        ps = branch_pages(branch, source_pages)
        score = sum(review_priority_score(p) for p in ps)
        rows.append((branch, ps, score, impact_profile(branch)))
    rows.sort(key=lambda item: (TIER_ORDER.get(item[3]["tier"], 9), -item[2], item[0]["paper"]))

    tier_counts: dict[str, int] = {}
    for _, _, _, profile in rows:
        tier_counts[profile["tier"]] = tier_counts.get(profile["tier"], 0) + 1

    headline_rows = []
    blocker_rows = []
    action_cards = []
    figure_sequence_rows = []
    for i, (branch, ps, score, profile) in enumerate(rows, 1):
        dossier = branch_dossier_filename(branch)
        primary = PRIMARY_PAGE_BY_BRANCH.get(branch["branch"], "paper_branch_manifest_2026_05_10.html")
        headline_rows.append(
            "<tr>"
            f"<td>{i}</td><td><span class=\"tag\">{html.escape(profile['tier'])}</span></td>"
            f"<td>{html.escape(branch['paper'])}</td><td>{html.escape(profile['route'])}</td>"
            f"<td>{score}</td><td>{len(ps)}</td><td>{sum(len(p.images) for p in ps)}</td><td>{sum(p.table_count for p in ps)}</td>"
            f"<td><a href=\"{html.escape(dossier)}\">dossier</a> · <a href=\"{html.escape(primary)}\">primary</a></td>"
            "</tr>"
        )
        blocker_rows.append(
            "<tr>"
            f"<td>{html.escape(branch['paper'])}</td>"
            f"<td>{html.escape(profile['blocker'])}</td>"
            f"<td>{html.escape(profile['next'])}</td>"
            f"<td>{html.escape(branch['boundary'])}</td>"
            "</tr>"
        )
        figure_sequence_rows.append(
            "<tr>"
            f"<td>{i}</td><td>{html.escape(branch['paper'])}</td>"
            f"<td>{html.escape(profile['killer'])}</td>"
            f"<td>{html.escape(branch['data'])}</td>"
            f"<td>{html.escape(branch['software'])}</td>"
            "</tr>"
        )
        if profile["tier"] in {"S", "A", "A-"}:
            action_cards.append(
                '<article class="step">'
                f'<strong>{html.escape(branch["paper"])} · tier {html.escape(profile["tier"])}</strong>'
                f'<p>{html.escape(profile["next"])}</p>'
                f'<p><span class="tag">{html.escape(profile["route"])}</span><a class="tag" href="{html.escape(dossier)}">dossier</a></p>'
                '</article>'
            )

    route_rows = "".join(hrow(r) for r in [
        ["Submit now", "Paper 1", "Evidence already dense; manuscript voice is bottleneck.", "Do not generate protected prose."],
        ["Scale next", "Paper 1B, Paper 11, Neo/CROSS-Neo", "High upside once Paper 1 is locked or kept separate.", "Keep claim boundaries explicit."],
        ["Validate experimentally", "Paper 9, Paper 13", "Actionability signal is strong enough to design wet-lab/spatial validation.", "No clinical IC50/drug efficacy claim yet."],
        ["Park/support", "Paper 5, Paper 6, Paper 7, Paper 10, Paper 12", "Useful infrastructure/supplement/resource branches.", "Do not let side branches dilute the main manuscript."],
        ["Data gate", "Paper 3, Paper 4", "High upside only if missing outcome/genotype data arrives.", "Freeze until required data exists."],
    ])

    return f"""<!DOCTYPE html><html lang="ko"><head><meta charset="UTF-8"/><meta name="viewport" content="width=device-width,initial-scale=1.0"/><title>Impact Escalation Matrix 2026-05-10</title>{atlas_style()}</head><body><div class="wrap">
{repo_nav()}
<header class="hero"><div class="kicker">impact escalation matrix · 2026-05-10</div><h1>Impact Escalation Matrix</h1><p class="subtitle">각 branch를 논문/제안서 관점에서 다시 정렬했습니다. 핵심은 tier, target route, killer figure, blocker, next 48h action입니다.</p><div class="stats"><div class="stat"><b>{tier_counts.get('S', 0)}</b><span>S tier</span></div><div class="stat"><b>{tier_counts.get('A', 0) + tier_counts.get('A-', 0)}</b><span>A/A- tier</span></div><div class="stat"><b>{sum(1 for _, _, _, p in rows if p['tier'].startswith('B'))}</b><span>B tier</span></div><div class="stat"><b>{len(rows)}</b><span>branches</span></div><div class="stat"><b>{sum(len(p.images) for p in source_pages)}</b><span>figures audited</span></div><div class="stat"><b>{sum(p.table_count for p in source_pages)}</b><span>tables audited</span></div></div></header>
<section><h2>Route Decision Matrix</h2><table>{hrow(["Route","Branches","Why","Boundary"], True)}{route_rows}</table></section>
<section><h2>Next 48h Action Queue</h2><div class="flow">{''.join(action_cards)}</div></section>
<section><h2>Ranked Escalation Table</h2><table>{hrow(["Rank","Tier","Branch","Target route","Impact score","Pages","Figures","Tables","Links"], True)}{''.join(headline_rows)}</table></section>
<section><h2>Killer Figure Sequence</h2><table>{hrow(["#","Branch","Killer figure","Data entered","Software/method"], True)}{''.join(figure_sequence_rows)}</table></section>
<section><h2>Blockers And Kill-Switches</h2><table>{hrow(["Branch","Current blocker","Next action","Do-not-overclaim boundary"], True)}{''.join(blocker_rows)}</table></section>
</div></body></html>"""


def build_killer_figure_storyboard(pages: list[PageMeta], idx: dict[str, list[str]]) -> str:
    source_pages = audit_pages(pages)
    selected: list[tuple[dict[str, str], list[PageMeta], dict[str, str], list[tuple[str, str, str, str, str, str]]]] = []
    for branch_id in KILLER_FIGURE_ORDER:
        branch = branch_by_id(branch_id)
        ps = branch_pages(branch, source_pages)
        selected.append((branch, ps, impact_profile(branch), killer_panels_for_branch(branch)))

    total_panels = sum(len(panels) for _, _, _, panels in selected)
    total_pages = sum(len(ps) for _, ps, _, _ in selected)
    total_figs = sum(len(p.images) for _, ps, _, _ in selected for p in ps)
    total_tables = sum(p.table_count for _, ps, _, _ in selected for p in ps)

    contract_rows = "".join(
        hrow(r)
        for r in [
            ["Panel label", "A-F panel identity per branch", "Turns each project into a figure package", "Do not merge unrelated branches in one panel."],
            ["Data entered", "Exact cohort, assay, matrix, or result file type", "Lets reviewer trace every visual claim", "Dataset naming must be checked against source file before manuscript use."],
            ["Software/method", "Scoring, residualization, model, audit, or HTML builder", "Shows how the output was produced", "Software output is not causal proof."],
            ["Narrative job", "What the panel should prove or clarify", "Prevents decorative figures", "Narrative stays inside evidence boundary."],
            ["Boundary", "Known caveat or kill-switch", "Keeps high impact honest", "Negative/caveat panels remain visible."],
            ["Source pages", "Top matched hub pages and linked assets", "Fast click-through to figures/tables", "Generated hub is a map; local result files remain source of truth."],
        ]
    )

    package_cards = []
    package_sections = []
    for rank, (branch, ps, profile, panels) in enumerate(selected, 1):
        dossier = branch_dossier_filename(branch)
        primary = PRIMARY_PAGE_BY_BRANCH.get(branch["branch"], "paper_branch_manifest_2026_05_10.html")
        panel_rows = "".join(
            hrow([label, title, data, method, use, boundary])
            for label, title, data, method, use, boundary in panels
        )
        panel_flow = "".join(
            '<article class="step">'
            f'<strong>Panel {html.escape(label)} · {html.escape(title)}</strong>'
            f'<p><b>Data:</b> {html.escape(data)}</p>'
            f'<p><b>Method:</b> {html.escape(method)}</p>'
            f'<p><b>Job:</b> {html.escape(use)}</p>'
            f'<p><span class="tag">{html.escape(boundary)}</span></p>'
            '</article>'
            for label, title, data, method, use, boundary in panels
        )
        source_cards = "".join(
            build_repo_page_card(p, idx)
            for p in sorted(ps, key=lambda p: -review_priority_score(p))[:6]
        ) or '<article class="card"><b>No matched source page yet</b><p>Use branch registry paths until a dedicated page is generated.</p></article>'
        source_rows = "".join(
            hrow(
                [
                    p.path.name,
                    infer_paper_group(p),
                    len(p.images),
                    p.table_count,
                    "; ".join(source_chips_for_page(p, idx, 4)) or "page-local content",
                ]
            )
            for p in sorted(ps, key=lambda p: -review_priority_score(p))[:10]
        ) or hrow(["No matched page", branch["paper"], 0, 0, branch["paths"]])

        package_cards.append(
            '<article class="card">'
            f'<b>#{rank} · {html.escape(branch["paper"])} · tier {html.escape(profile["tier"])}</b>'
            f'<strong>{html.escape(branch["title"])}</strong>'
            f'<p><span class="tag">{len(panels)} figure panels</span><span class="tag">{len(ps)} pages</span><span class="tag">{sum(len(p.images) for p in ps)} figs</span><span class="tag">{sum(p.table_count for p in ps)} tables</span></p>'
            f'<p><b>Killer figure:</b> {html.escape(profile["killer"])}</p>'
            f'<p><b>Next action:</b> {html.escape(profile["next"])}</p>'
            f'<p><a href="#pkg-{html.escape(branch["branch"])}">storyboard</a> · <a href="{html.escape(dossier)}">dossier</a> · <a href="{html.escape(primary)}">primary</a></p>'
            '</article>'
        )

        package_sections.append(
            f"""<section id="pkg-{html.escape(branch['branch'])}">
<h2>{rank:02d}. {html.escape(branch['paper'])} Killer Figure Package</h2>
<div class="warning"><b>Target route:</b> {html.escape(profile['route'])}<br/><b>Killer figure:</b> {html.escape(profile['killer'])}<br/><b>Boundary:</b> {html.escape(profile['blocker'])}</div>
<table>{hrow(["Branch field","Value","How to use"], True)}{hrow(["Data definition", branch["data"], "Input material for figure panels and result tables."])}{hrow(["Software/method ledger", branch["software"], "Methods box and reproducibility path."])}{hrow(["Git-style history", branch["history"], "Explains why this branch exists and what it inherits."])}{hrow(["Source paths", branch["paths"], "Click-through/local provenance before manuscript claim."])}</table>
<h3>Panel Flow</h3><div class="flow">{panel_flow}</div>
<h3>Panel Assembly Table</h3><table>{hrow(["Panel","Figure panel title","Data entered","Software/method","Narrative job","Boundary"], True)}{panel_rows}</table>
<h3>Representative Source Pages</h3><div class="grid">{source_cards}</div>
<h3>Source Ledger</h3><table>{hrow(["Page","Family","Figures","Tables","Source chips"], True)}{source_rows}</table>
</section>"""
        )

    priority_rows = "".join(
        hrow(
            [
                str(i),
                branch["paper"],
                profile["tier"],
                profile["route"],
                profile["killer"],
                profile["next"],
                branch["boundary"],
            ]
        )
        for i, (branch, _, profile, _) in enumerate(selected, 1)
    )

    return f"""<!DOCTYPE html><html lang="ko"><head><meta charset="UTF-8"/><meta name="viewport" content="width=device-width,initial-scale=1.0"/><title>Killer Figure Storyboard 2026-05-10</title>{atlas_style()}</head><body><div class="wrap">
{repo_nav()}
<header class="hero"><div class="kicker">high-impact figure architecture · 2026-05-10</div><h1>Killer Figure Storyboard</h1><p class="subtitle">Paper 1 수준의 figure-heavy dossier를 다른 핵심 branch에도 적용하기 위한 패널 설계도입니다. 각 패널은 데이터, 소프트웨어/방법, 해석 역할, caveat, source page를 같이 묶습니다.</p><div class="stats"><div class="stat"><b>{len(selected)}</b><span>priority packages</span></div><div class="stat"><b>{total_panels}</b><span>defined panels</span></div><div class="stat"><b>{total_pages}</b><span>matched source pages</span></div><div class="stat"><b>{total_figs}</b><span>source figures</span></div><div class="stat"><b>{total_tables}</b><span>source tables</span></div><div class="stat"><b>0</b><span>protected prose edits</span></div></div></header>
<section><h2>How To Use This Page</h2><table>{hrow(["Field","Data entered","Why it matters","Boundary"], True)}{contract_rows}</table></section>
<section><h2>Priority Figure Queue</h2><table>{hrow(["Rank","Branch","Tier","Target route","Killer figure","Next action","Boundary"], True)}{priority_rows}</table></section>
<section><h2>Package Cards</h2><div class="grid">{''.join(package_cards)}</div></section>
{''.join(package_sections)}
</div></body></html>"""


def build_paper_repo_control_room(pages: list[PageMeta], idx: dict[str, list[str]]) -> str:
    source_pages = audit_pages(pages)
    counts: dict[str, tuple[int, int, int]] = {}
    for branch in PAPER_BRANCHES:
        ps = branch_pages(branch, source_pages)
        counts[branch["branch"]] = (len(ps), sum(len(p.images) for p in ps), sum(p.table_count for p in ps))

    lane_sections = []
    for lane in REPO_LANES:
        mini_cards = []
        for branch_id in lane["branches"]:
            b = branch_by_id(branch_id)
            pc, fc, tc = counts[b["branch"]]
            mini_cards.append(
                '<article class="card">'
                f'<b>{html.escape(b["paper"])} · {html.escape(b["state"])}</b>'
                f'<strong>{html.escape(b["title"])}</strong>'
                f'<p><span class="tag">{html.escape(b["topic"])}</span><span class="tag">parent: {html.escape(b["parent"])}</span><span class="tag">last: {html.escape(b["last_update"])}</span></p>'
                f'<p>{html.escape(b["meaning"])}</p>'
                f'<p><b>Boundary:</b> {html.escape(b["boundary"])}</p>'
                f'<p><span class="tag">{pc} pages</span><span class="tag">{fc} figs</span><span class="tag">{tc} tables</span> {primary_page_link(b)} · {branch_dossier_link(b)}</p>'
                '</article>'
            )
        lane_sections.append(
            '<section>'
            f'<h2>{html.escape(lane["lane"])}</h2>'
            f'<p>{html.escape(lane["meaning"])}</p>'
            f'<div class="grid">{"".join(mini_cards)}</div>'
            '</section>'
        )

    table_rows = []
    for b in PAPER_BRANCHES:
        pc, fc, tc = counts[b["branch"]]
        table_rows.append(
            "<tr>"
            f"<td>{html.escape(b['paper'])}</td>"
            f"<td>{html.escape(b['topic'])}</td>"
            f"<td>{html.escape(b['state'])}</td>"
            f"<td>{html.escape(b['parent'])}</td>"
            f"<td>{html.escape(b['opened'])}</td>"
            f"<td>{html.escape(b['last_update'])}</td>"
            f"<td>{pc}</td><td>{fc}</td><td>{tc}</td>"
            f"<td>{primary_page_link(b)} · {branch_dossier_link(b)}</td>"
            "</tr>"
        )

    date_rows = []
    for date in sorted({b.get("last_update", "") for b in PAPER_BRANCHES}, reverse=True):
        bs = [b for b in PAPER_BRANCHES if b.get("last_update", "") == date]
        date_rows.append(f"<tr><td>{html.escape(date)}</td><td>{' '.join(branch_tag_link(b) for b in bs)}</td></tr>")

    graph_rows = []
    for b in PAPER_BRANCHES:
        graph_rows.append(hrow([b["parent"], b["branch"], b["paper"], b["state"], b["history"]]))

    next_rows = [
        ("1", "Paper 1만 submit lane으로 고정", "새 분석보다 author-voice manuscript completion이 병목입니다."),
        ("2", "IC50/PRISM은 Paper 9 active branch", "버린 프로젝트가 아니라 actionability map이며 wet-lab queue까지 연결합니다."),
        ("3", "Paper 1B와 Paper 11은 다음 큰 후보", "둘 다 upside가 크지만 Paper 1을 방해하면 안 됩니다."),
        ("4", "Project inventory는 no-abandoned-project ledger", "owner가 약한 결과물도 parking-lot으로 보존하고 다음 audit에서 승격합니다."),
        ("5", "Paper 3/4는 데이터 gate", "thyroid ICI outcome, Korean adult GD NGS가 오기 전까지는 frozen/gated입니다."),
        ("6", "Neo/CROSS-Neo는 별도 repo", "THCA Paper 1-4와 섞지 않고 methods/resource로 관리합니다."),
    ]
    next_html = "".join(hrow(r) for r in next_rows)

    return f"""<!DOCTYPE html><html lang="ko"><head><meta charset="UTF-8"/><meta name="viewport" content="width=device-width,initial-scale=1.0"/><title>Paper Repo Control Room 2026-05-10</title>{atlas_style()}</head><body><div class="wrap">
{repo_nav()}
<header class="hero"><div class="kicker">repo control room · 2026-05-10</div><h1>Paper Repo Control Room</h1><p class="subtitle">흩어진 paper 후보를 14개 numbered paper slots + 1 side repo로 정리한 최상위 지도입니다. 여기서 날짜, 주제, parent branch, status, primary page를 먼저 보고 들어가면 됩니다.</p><div class="stats"><div class="stat"><b>{len(numbered_paper_branches())}</b><span>paper slots</span></div><div class="stat"><b>{len(side_repo_branches())}</b><span>side repos</span></div><div class="stat"><b>{len(REPO_LANES)}</b><span>lanes</span></div><div class="stat"><b>{len(source_pages)}</b><span>tracked pages</span></div><div class="stat"><b>{sum(len(p.images) for p in source_pages)}</b><span>figures audited</span></div><div class="stat"><b>{sum(p.table_count for p in source_pages)}</b><span>tables audited</span></div></div></header>
<section><h2>Read First</h2><table>{hrow(["#","Rule","Why"], True)}{next_html}</table></section>
{''.join(lane_sections)}
<section><h2>Date Buckets</h2><table>{hrow(["Last update","Branches"], True)}{''.join(date_rows)}</table></section>
<section><h2>Registry Table</h2><table>{hrow(["Paper","Topic","State","Parent","Opened","Last update","Pages","Figures","Tables","Primary page"], True)}{''.join(table_rows)}</table></section>
<section><h2>Parent / History Graph</h2><table>{hrow(["Parent","Branch","Paper","State","History"], True)}{''.join(graph_rows)}</table></section>
<section><h2>Machine Registry</h2><p><a href="high_impact_board_2026_05_10.html">High-impact board</a> · <a href="paper_branch_registry_2026_05_10.tsv">TSV registry</a> · <a href="paper_branch_registry_2026_05_10.json">JSON registry</a> · <a href="project_inventory_2026_05_10.html">Project inventory / parking lot</a></p></section>
</div></body></html>"""


def build_paper_repo_map(pages: list[PageMeta], idx: dict[str, list[str]]) -> str:
    source_pages = audit_pages(pages)
    branch_cards = []
    graph_rows = []
    for branch in PAPER_BRANCHES:
        ps = branch_pages(branch, source_pages)
        graph_rows.append(hrow([
            branch["paper"],
            branch["branch"],
            branch.get("slot_type", "paper"),
            branch.get("topic", ""),
            branch.get("opened", ""),
            branch.get("last_update", ""),
            branch["parent"],
            branch["state"],
            len(ps),
            sum(len(p.images) for p in ps),
            sum(p.table_count for p in ps),
            branch["history"],
        ]))
        branch_cards.append(
            '<article class="card">'
            f'<b>{html.escape(branch["paper"])} · {html.escape(branch.get("topic", ""))}</b>'
            f'<strong>{html.escape(branch["branch"])} / {html.escape(branch["paper"])}</strong>'
            f'<p>{html.escape(branch["title"])}</p><p><span class="tag">opened: {html.escape(branch.get("opened", ""))}</span><span class="tag">last: {html.escape(branch.get("last_update", ""))}</span><span class="tag">{html.escape(branch["state"])}</span><span class="tag">{len(ps)} pages</span><span class="tag">{sum(len(p.images) for p in ps)} figs</span><span class="tag">{sum(p.table_count for p in ps)} tables</span></p>'
            f'<p><span class="tag">parent: {html.escape(branch["parent"])}</span><span class="tag">slot: {html.escape(branch.get("slot_type", "paper"))}</span></p>'
            f'<p><b>Claim:</b> {html.escape(branch["claim"])}</p>'
            f'<p><b>Data:</b> {html.escape(branch["data"])}</p>'
            f'<p><b>Software:</b> {html.escape(branch["software"])}</p>'
            f'<p><b>Boundary:</b> {html.escape(branch["boundary"])}</p>'
            '</article>'
        )
    date_rows = []
    for date in sorted({b.get("last_update", "") for b in PAPER_BRANCHES}, reverse=True):
        bs = [b for b in PAPER_BRANCHES if b.get("last_update", "") == date]
        date_rows.append(f"<tr><td>{html.escape(date)}</td><td>{' '.join(branch_tag_link(b) for b in bs)}</td></tr>")
    topic_rows = []
    for topic in sorted({b.get("topic", "") for b in PAPER_BRANCHES}):
        bs = [b for b in PAPER_BRANCHES if b.get("topic", "") == topic]
        topic_rows.append(f"<tr><td>{html.escape(topic)}</td><td>{' '.join(branch_tag_link(b) for b in bs)}</td></tr>")
    main_count = len(numbered_paper_branches())
    side_count = len(side_repo_branches())
    return f"""<!DOCTYPE html><html lang="ko"><head><meta charset="UTF-8"/><meta name="viewport" content="width=device-width,initial-scale=1.0"/><title>Paper Repo Map 2026-05-10</title>{atlas_style()}</head><body><div class="wrap">
{repo_nav()}
<header class="hero"><div class="kicker">papers_hub_2026_0509 · repo-style portfolio map · 2026-05-10</div><h1>Paper Repository Map</h1><p class="subtitle">논문 포트폴리오를 git repo처럼 나눈 지도입니다. 14개 numbered paper slots는 main registry, Neo/CROSS-Neo는 separate side repo입니다. 각 branch는 opened/last_update, topic, parent, history, tracked files, data/software/boundary를 가집니다.</p><div class="stats"><div class="stat"><b>{main_count}</b><span>numbered paper slots</span></div><div class="stat"><b>{side_count}</b><span>side repos</span></div><div class="stat"><b>{len(source_pages)}</b><span>tracked pages</span></div><div class="stat"><b>{sum(len(p.images) for p in source_pages)}</b><span>tracked figures</span></div><div class="stat"><b>{sum(p.table_count for p in source_pages)}</b><span>tracked tables</span></div><div class="stat"><b>0</b><span>protected prose edits</span></div></div></header>
<section><h2>Date Index</h2><table>{hrow(["Last update","Branches"], True)}{''.join(date_rows)}</table></section>
<section><h2>Topic Index</h2><table>{hrow(["Topic","Branches"], True)}{''.join(topic_rows)}</table></section>
<section><h2>Repository Graph</h2><table>{hrow(["Paper","Branch","Slot","Topic","Opened","Last update","Parent","State","Pages","Figures","Tables","History"], True)}{''.join(graph_rows)}</table></section>
<section><h2>Branch Cards</h2><div class="grid">{''.join(branch_cards)}</div></section>
<section><h2>Interpretation Rule</h2><div class="warning"><b>Paper 1 is trunk.</b> BRAF, WSI, HT/HLA, ICI, IC50/PRISM, pan-cancer, deconvolution, spatial drug, and network branches can reference trunk history, but each keeps its own claim boundary and data provenance. Neo/CROSS-Neo is a separate methods repo.</div></section>
</div></body></html>"""


def build_paper_branch_manifest(pages: list[PageMeta], idx: dict[str, list[str]]) -> str:
    rows = []
    cards = []
    for branch in PAPER_BRANCHES:
        ps = branch_pages(branch, pages)
        links = " ".join(f'<a class="tag" href="{html.escape(p.path.name)}">{html.escape(p.path.stem[:42])}</a>' for p in sorted(ps, key=lambda p: -review_priority_score(p))[:12])
        rows.append(hrow([
            branch["paper"],
            branch["branch"],
            branch.get("slot_type", "paper"),
            branch.get("topic", ""),
            branch.get("opened", ""),
            branch.get("last_update", ""),
            branch["state"],
            branch["parent"],
            len(ps),
            sum(len(p.images) for p in ps),
            sum(p.table_count for p in ps),
            branch["data"],
            branch["software"],
            branch["boundary"],
        ]))
        cards.append(
            '<article class="card">'
            f'<b>{html.escape(branch["paper"])} · {html.escape(branch.get("topic", ""))}</b><strong>{html.escape(branch["title"])}</strong>'
            f'<p><span class="tag">opened: {html.escape(branch.get("opened", ""))}</span><span class="tag">last: {html.escape(branch.get("last_update", ""))}</span><span class="tag">{html.escape(branch["state"])}</span><span class="tag">parent: {html.escape(branch["parent"])}</span></p>'
            f'<p><b>Meaning:</b> {html.escape(branch["meaning"])}</p><p><b>Tracked pages:</b> {links or "none"}</p>'
            '</article>'
        )
    main_count = len(numbered_paper_branches())
    side_count = len(side_repo_branches())
    return f"""<!DOCTYPE html><html lang="ko"><head><meta charset="UTF-8"/><meta name="viewport" content="width=device-width,initial-scale=1.0"/><title>Paper Branch Manifest 2026-05-10</title>{atlas_style()}</head><body><div class="wrap">
{repo_nav()}
<header class="hero"><div class="kicker">repo map · branch manifest</div><h1>Paper Branch Manifest</h1><p class="subtitle">각 paper branch가 어떤 날짜, 주제, 데이터, 소프트웨어, claim boundary, tracked pages를 갖는지 보여주는 manifest입니다. Main registry는 14 numbered paper slots입니다.</p><div class="stats"><div class="stat"><b>{main_count}</b><span>numbered paper slots</span></div><div class="stat"><b>{side_count}</b><span>side repos</span></div><div class="stat"><b>data</b><span>driven claims</span></div><div class="stat"><b>0</b><span>protected prose edits</span></div></div></header>
<section><h2>Branch Cards</h2><div class="grid">{''.join(cards)}</div></section>
<section><h2>Manifest Table</h2><table>{hrow(["Paper","Branch","Slot","Topic","Opened","Last update","State","Parent","Pages","Figures","Tables","Data","Software","Boundary"], True)}{''.join(rows)}</table></section>
</div></body></html>"""


def build_paper_date_topic_index(pages: list[PageMeta], idx: dict[str, list[str]]) -> str:
    source_pages = audit_pages(pages)
    rows = []
    for branch in PAPER_BRANCHES:
        ps = branch_pages(branch, source_pages)
        rows.append(hrow([
            branch.get("last_update", ""),
            branch.get("opened", ""),
            branch.get("topic", ""),
            branch.get("paper", ""),
            branch.get("branch", ""),
            branch.get("state", ""),
            branch.get("parent", ""),
            len(ps),
            sum(len(p.images) for p in ps),
            sum(p.table_count for p in ps),
        ]))
    date_sections = []
    for date in sorted({b.get("last_update", "") for b in PAPER_BRANCHES}, reverse=True):
        bs = [b for b in PAPER_BRANCHES if b.get("last_update", "") == date]
        date_sections.append(
            '<article class="card">'
            f'<b>last update · {html.escape(date)}</b>'
            f'<p>{" ".join(branch_tag_link(b) for b in bs)}</p>'
            '</article>'
        )
    topic_sections = []
    for topic in sorted({b.get("topic", "") for b in PAPER_BRANCHES}):
        bs = [b for b in PAPER_BRANCHES if b.get("topic", "") == topic]
        topic_sections.append(
            '<article class="card">'
            f'<b>topic</b><strong>{html.escape(topic)}</strong>'
            f'<p>{" ".join(branch_tag_link(b) for b in bs)}</p>'
            '</article>'
        )
    return f"""<!DOCTYPE html><html lang="ko"><head><meta charset="UTF-8"/><meta name="viewport" content="width=device-width,initial-scale=1.0"/><title>Paper Date Topic Index 2026-05-10</title>{atlas_style()}</head><body><div class="wrap">
{repo_nav()}
<header class="hero"><div class="kicker">repo map · date/topic index</div><h1>Paper Date / Topic Index</h1><p class="subtitle">14개 paper slot과 Neo side repo를 날짜별, 주제별, parent/history별로 보는 index입니다. repo 관리 관점에서는 이 표가 branch registry의 source of truth입니다.</p><div class="stats"><div class="stat"><b>{len(numbered_paper_branches())}</b><span>numbered paper slots</span></div><div class="stat"><b>{len(side_repo_branches())}</b><span>side repos</span></div><div class="stat"><b>{len({b.get("topic", "") for b in PAPER_BRANCHES})}</b><span>topics</span></div><div class="stat"><b>{len({b.get("last_update", "") for b in PAPER_BRANCHES})}</b><span>update dates</span></div></div></header>
<section><h2>Machine Registry Files</h2><p><a href="paper_branch_registry_2026_05_10.tsv">paper_branch_registry_2026_05_10.tsv</a> · <a href="paper_branch_registry_2026_05_10.json">paper_branch_registry_2026_05_10.json</a></p></section>
<section><h2>Date Buckets</h2><div class="grid">{''.join(date_sections)}</div></section>
<section><h2>Topic Buckets</h2><div class="grid">{''.join(topic_sections)}</div></section>
<section><h2>Sortable Registry Table</h2><table>{hrow(["Last update","Opened","Topic","Paper","Branch","State","Parent","Pages","Figures","Tables"], True)}{''.join(rows)}</table></section>
</div></body></html>"""


def build_software_method_ledger(pages: list[PageMeta]) -> str:
    cards = "".join(
        '<article class="card">'
        f'<b>method ledger</b><strong>{html.escape(name)}</strong><p><b>Software:</b> {html.escape(software)}</p><p><b>How to interpret:</b> {html.escape(read)}</p>'
        '</article>'
        for name, software, read in SOFTWARE_METHODS
    )
    rows = "".join(hrow([b["branch"], b["data"], b["software"], b["boundary"]]) for b in PAPER_BRANCHES)
    return f"""<!DOCTYPE html><html lang="ko"><head><meta charset="UTF-8"/><meta name="viewport" content="width=device-width,initial-scale=1.0"/><title>Software Method Ledger 2026-05-10</title>{atlas_style()}</head><body><div class="wrap">
{repo_nav()}
<header class="hero"><div class="kicker">repo map · software/method ledger</div><h1>Software And Method Ledger</h1><p class="subtitle">각 분석에서 어떤 데이터와 소프트웨어를 썼고, 그 결과를 어떻게 해석해야 하는지 정리한 ledger입니다.</p><div class="stats"><div class="stat"><b>{len(SOFTWARE_METHODS)}</b><span>method blocks</span></div><div class="stat"><b>{len(PAPER_BRANCHES)}</b><span>branches</span></div><div class="stat"><b>0</b><span>protected prose edits</span></div></div></header>
<section><h2>Method Cards</h2><div class="grid">{cards}</div></section>
<section><h2>Branch Method Matrix</h2><table>{hrow(["Branch","Data","Software","Interpretation boundary"], True)}{rows}</table></section>
</div></body></html>"""


def build_ic50_prediction_map(pages: list[PageMeta], idx: dict[str, list[str]]) -> str:
    branch = next(b for b in PAPER_BRANCHES if b["branch"] == "paper9-ic50-perturbation")
    ps = branch_pages(branch, pages)
    cards = "".join(build_repo_page_card(p, idx) for p in sorted(ps, key=lambda p: -review_priority_score(p)))
    summary_path = ROOT / "project" / "results" / "p_synthetic_lethality_2026_05_09" / "synthetic_lethality_summary.json"
    candidate_path = ROOT / "project" / "results" / "p_synthetic_lethality_2026_05_09" / "synthetic_lethality_candidate_matrix.tsv"
    wetlab_path = ROOT / "project" / "results" / "paper9_perturbation_2026_05_06" / "paper9_wetlab_perturbation_matrix.tsv"
    summary = json.loads(summary_path.read_text(encoding="utf-8")) if summary_path.exists() else {}
    prism = summary.get("prism", {})
    depmap = summary.get("depmap", {})
    metabolic = summary.get("metabolic", {})
    candidates = read_tsv_dicts(candidate_path)
    wetlab = read_tsv_dicts(wetlab_path)

    def fmt_num(value: object, digits: int = 3) -> str:
        if value is None or value == "":
            return "n/a"
        try:
            value_f = float(value)
        except (TypeError, ValueError):
            return str(value)
        if abs(value_f) < 0.001 and value_f != 0:
            return f"{value_f:.2e}"
        return f"{value_f:.{digits}f}"

    metric_rows = "".join(hrow(r) for r in [
        ["PRISM universe", f"{prism.get('n_drugs', 'n/a')} drugs; FDR hits {prism.get('fdr_hits', 'n/a')}; canonical MAPK top7 {prism.get('top7_canonical_mapk', 'n/a')}; top15 {prism.get('top15_canonical_mapk', 'n/a')}", "This is the strongest actionability signal. Treat as drug-response proxy, not clinical IC50."],
        ["AZD-0364 anchor", f"delta LFC {fmt_num(prism.get('azd_delta_lfc'))}; Cohen's d {fmt_num(prism.get('azd_d'))}; FDR {fmt_num(prism.get('azd_fdr'))}", "Use as a named MAPK-inhibitor anchor for reviewer reserve and wet-lab selection."],
        ["DM1 x MAPK drug response", f"rho all {fmt_num(prism.get('rho_all'))}, p {fmt_num(prism.get('rho_all_p'))}; rho excluding thyroid {fmt_num(prism.get('rho_exclude_thyroid'))}, p {fmt_num(prism.get('rho_exclude_thyroid_p'))}", "Negative rho means higher DM1 score tracks stronger drug-response signal in the public screen direction."],
        ["DepMap MYC/NAMPT", f"MYC d {fmt_num(depmap.get('myc_d'))}, p {fmt_num(depmap.get('myc_p'))}; NAMPT d {fmt_num(depmap.get('nampt_d'))}, p {fmt_num(depmap.get('nampt_p'))}", "Functional prioritization only; common-essential/pan-essential risk needs controls."],
        ["Metabolic candidates", "; ".join(f"{gene} rho {fmt_num(vals.get('dependency_rho'))}, p {fmt_num(vals.get('dependency_p'))}" for gene, vals in metabolic.items()), "SLC1A5 > GLUD1 > GLS for thyroid wet-lab queue; strict thyroid n=11 remains underpowered."],
    ])
    candidate_rows = "".join(
        hrow([
            row.get("tier", ""),
            row.get("candidate_program", ""),
            row.get("perturbation", ""),
            row.get("headline_metrics", ""),
            row.get("disposition", ""),
            row.get("required_next_validation", ""),
        ])
        for row in candidates
    )
    wetlab_rows = "".join(
        hrow([
            row.get("model_state", ""),
            row.get("perturbation", ""),
            row.get("candidate", ""),
            row.get("readouts", ""),
            row.get("go_criterion", ""),
            row.get("falsification", ""),
        ])
        for row in wetlab
    )
    decision_rows = "".join(hrow(r) for r in [
        ["GO", "MAPK inhibitor PRISM axis", "Keep as Paper 9 actionability map and Paper 1 reviewer reserve.", "PRISM response proxy; not direct IC50 or clinical response."],
        ["GO", "SLC1A5 / GLUD1 / GLS wet-lab queue", "Turn into perturbation-rescue experiment plan.", "Needs matched thyroid lineage-high vs lineage-silenced models."],
        ["HOLD", "MYC / NAMPT", "Keep as functional prioritization and comparator controls.", "Pan-essentiality/common-essential risk must be separated."],
        ["NO HEADLINE", "Spatial anti-correlation rescue", "Use only as honest negative/caveat.", "v16/v17 spatial tests did not support mechanism."],
    ])
    flow_cards = "".join(
        '<article class="step">'
        f'<strong>{html.escape(title)}</strong>'
        f'<p>{html.escape(text)}</p>'
        f'<p><span class="tag">{html.escape(data)}</span><span class="tag">{html.escape(software)}</span></p>'
        '</article>'
        for title, text, data, software in [
            ("1. Define state", "DM1-high / lineage-silenced state is the index state; lineage-high or DM1-low is the comparator.", "DM1 score, TDS/RAI genes, cell-line metadata", "Python/pandas module scoring"),
            ("2. Screen drug response", "PRISM drug-response matrix is scanned for compounds whose response tracks the DM1/MAPK axis.", "PRISM 1518-drug matrix", "Spearman, FDR, Cohen's d"),
            ("3. Name the pharmacology axis", "The top signal is not random drugs: canonical MAPK inhibitors dominate the top ranked set.", "MAPK/MEK/RAF/ERK drug annotations", "Drug-class enrichment"),
            ("4. Cross-check dependency", "DepMap and metabolic dependencies give orthogonal genetic-priority candidates.", "DepMap CRISPR/dependency rows", "High-vs-low effect sizes"),
            ("5. Translate to wet-lab", "SLC1A5, GLUD1, and GLS become perturbation-rescue experiments with explicit falsification rules.", "wetlab perturbation matrix", "validation design table"),
            ("6. Lock boundary", "The current evidence supports actionability prioritization, not validated clinical IC50 or synthetic lethality.", "claim ledger", "repo-style boundary control"),
        ]
    )
    overview_figures = "".join(
        '<article class="card">'
        f'<b>{html.escape(fid)}</b>'
        f'<strong>{html.escape(title)}</strong>'
        f'<p><b>Data entered:</b> {html.escape(data)}</p>'
        f'<p><b>What it explains:</b> {html.escape(read)}</p>'
        f'<p><b>Boundary:</b> {html.escape(boundary)}</p>'
        '</article>'
        for fid, title, data, read, boundary in [
            ("Figure 1", "Paper 9 branch overview map", "Paper 9 registry, Paper 10 atlas, Paper 13 spatial-drug branch, project inventory owner tags", "Shows IC50/PRISM as an active actionability branch connected to atlas and spatial validation layers.", "Branch map is project organization, not a biological result."),
            ("Figure 2", "PRISM MAPK-inhibitor signal", "synthetic_lethality_summary.json: n_drugs, FDR hits, top7/top15 MAPK enrichment, AZD-0364 metrics", "Shows why the strongest current drug-response story is MAPK-axis sensitivity in the public screen.", "PRISM is a public response proxy; dose-response IC50 is still wet-lab work."),
            ("Figure 3", "DM1 score versus drug-response direction", "PRISM rho_all / rho_exclude_thyroid values and MAPK LFC direction", "Shows that the signal remains when thyroid lines are excluded, reducing one-lineage artifact concern.", "Correlation is not causal drug mechanism."),
            ("Figure 4", "DepMap dependency overlay", "MYC, NAMPT, SLC1A5, GLUD1, GLS dependency effect sizes and p-values", "Separates pharmacologic MAPK signal from genetic/metabolic candidate prioritization.", "MYC/NAMPT have common-essential risk and need controls."),
            ("Figure 5", "Wet-lab translation queue", "paper9_wetlab_perturbation_matrix.tsv: model state, perturbation, readouts, go criteria, falsification", "Turns the map into testable experiments instead of leaving it as a vague possibility.", "No wet-lab result is claimed until rescue/falsification passes."),
            ("Figure 6", "Negative/caveat pocket", "spatial v16/v17 negative MAPK anti-correlation rows from candidate matrix boundary lines", "Keeps failed spatial-rescue evidence visible so the page is reviewer-honest.", "Use as caveat, not as support."),
        ]
    )
    extracted_figures = []
    for page in sorted(ps, key=lambda p: -review_priority_score(p)):
        for i, img in enumerate(page.images, 1):
            if len(extracted_figures) >= 12:
                break
            caption = page.captions[i - 1] if i - 1 < len(page.captions) else img.get("alt", "")
            caption = normalize(caption) or "No caption/alt detected."
            role, how, caveat = visual_reading_notes(page, img, caption)
            src = html.escape(img.get("src", ""), quote=True)
            extracted_figures.append(
                '<article class="card figcard">'
                f'<a class="thumb" href="{html.escape(page.path.name)}" style="background-image:url({src})"></a>'
                '<div>'
                f'<b>Evidence figure · {html.escape(role)}</b>'
                f'<strong><a href="{html.escape(page.path.name)}">{html.escape(Path(img.get("src", "")).name or page.path.name)}</a></strong>'
                f'<p>{html.escape(caption[:420])}</p>'
                f'<dl class="dl"><dt>Data</dt><dd>{html.escape(guess_visual_data(page, img))}</dd><dt>Read</dt><dd>{html.escape(how)}</dd><dt>Boundary</dt><dd>{html.escape(caveat)}</dd></dl>'
                '</div></article>'
            )
        if len(extracted_figures) >= 12:
            break
    extracted_fig_html = "".join(extracted_figures) or '<article class="card"><b>No source-page figures</b><p>The IC50 map still has structured overview figures and source tables above.</p></article>'
    table_explanations = "".join(
        '<article class="card">'
        f'<b>{html.escape(tid)}</b>'
        f'<strong>{html.escape(title)}</strong>'
        f'<p><b>Rows/columns:</b> {html.escape(rows_desc)}</p>'
        f'<p><b>Use:</b> {html.escape(use)}</p>'
        f'<p><b>Boundary:</b> {html.escape(boundary)}</p>'
        '</article>'
        for tid, title, rows_desc, use, boundary in [
            ("Table 1", "Source-of-truth metrics", "Layer, current result, interpretation; values pulled from synthetic_lethality_summary.json", "Front-door quantitative summary for PRISM/DepMap/metabolic claims.", "Use only values linked to the JSON source."),
            ("Table 2", "Decision board", "Disposition, program, use now, boundary", "Tells which programs are GO, HOLD, or negative caveat.", "Does not replace wet-lab validation."),
            ("Table 3", "Candidate matrix", "Tier, program, perturbation, headline metrics, disposition, next validation", "Paper 9 program registry with source paths and next experiments.", "Candidate status is explicit; no validated synthetic-lethal claim."),
            ("Table 4", "Wet-lab queue", "Model state, perturbation, candidate, readouts, go criterion, falsification", "Converts computational map into executable validation design.", "Only a design until experiments are run."),
            ("Table 5", "Source paths", "File path and how the page uses it", "Lets every metric/table be traced back to TSV/JSON/MD evidence.", "If source file changes, rebuild the page."),
        ]
    )
    data_ledger_rows = "".join(hrow(r) for r in [
        ["DM1 score / TDS axis", "Paper 1 trunk labels and thyroid-lineage/dedifferentiation scores", "Defines high/low state for drug and dependency comparisons", "Do not reinterpret as causal without perturbation"],
        ["PRISM drug response", "project/results/p_synthetic_lethality_2026_05_09/synthetic_lethality_summary.json", "Main MAPK-inhibitor actionability map and AZD-0364 anchor", "Proxy response; not direct clinical IC50"],
        ["DepMap dependency", "same JSON plus Paper 11/DepMap support rows", "MYC/NAMPT and metabolic genetic-priority overlay", "Common-essential and lineage confounding controls required"],
        ["Candidate matrix", "synthetic_lethality_candidate_matrix.tsv", "Program-level GO/HOLD/caveat decisions", "Disposition is not validation"],
        ["Wet-lab matrix", "paper9_wetlab_perturbation_matrix.tsv", "Executable perturbation-rescue plan", "No biological claim until experiment passes go criteria"],
        ["Spatial drug branch", "Paper 13 spatial-drug pages and v16/v17 caveat rows", "Visual/contextual branch for where targets may overlap", "Co-localization is not drug efficacy"],
    ])
    return f"""<!DOCTYPE html><html lang="ko"><head><meta charset="UTF-8"/><meta name="viewport" content="width=device-width,initial-scale=1.0"/><title>IC50 Prediction Map 2026-05-10</title>{atlas_style()}</head><body><div class="wrap">
{repo_nav()}
<header class="hero"><div class="kicker">paper9 active branch · IC50/PRISM/DepMap actionability map</div><h1>IC50 Prediction Map</h1><p class="subtitle">이 페이지는 버려진 아이디어가 아니라 Paper 9의 active actionability branch입니다. 현재 강한 부분은 PRISM MAPK-inhibitor drug-response proxy이고, wet-lab으로 넘길 부분은 SLC1A5/GLUD1/GLS perturbation-rescue queue입니다.</p><div class="stats"><div class="stat"><b>{prism.get('n_drugs', '1518')}</b><span>PRISM drugs</span></div><div class="stat"><b>{prism.get('top7_canonical_mapk', '7/7')}</b><span>top7 MAPK</span></div><div class="stat"><b>{fmt_num(prism.get('azd_fdr'), 2)}</b><span>AZD-0364 FDR</span></div><div class="stat"><b>{len(ps)}</b><span>tracked pages</span></div><div class="stat"><b>{sum(len(p.images) for p in ps)}</b><span>figures</span></div><div class="stat"><b>{sum(p.table_count for p in ps)}</b><span>tables</span></div></div></header>
<section><h2>Claim Boundary</h2><div class="warning">{html.escape(branch["boundary"])} <b>Operational rule:</b> call it IC50 map only as shorthand; the data source is public drug-response/dependency proxy until dose-response wet-lab validation exists.</div></section>
<section><h2>Whole-Analysis Flow</h2><div class="flow">{flow_cards}</div></section>
<section><h2>Data-To-Use Ledger</h2><table>{hrow(["Data block","Source","Where it is used","Interpretation boundary"], True)}{data_ledger_rows}</table></section>
<section><h2>Representative Overview Figures</h2><div class="grid">{overview_figures}</div></section>
<section><h2>Source Of Truth Metrics</h2><table>{hrow(["Layer","Current result","How to interpret"], True)}{metric_rows}</table></section>
<section><h2>Decision Board</h2><table>{hrow(["Disposition","Program","Use now","Boundary"], True)}{decision_rows}</table></section>
<section><h2>Candidate Matrix</h2><table>{hrow(["Tier","Candidate program","Perturbation","Headline metrics","Disposition","Required next validation"], True)}{candidate_rows}</table></section>
<section><h2>Wet-Lab Queue</h2><table>{hrow(["Model state","Perturbation","Candidate","Readouts","Go criterion","Falsification"], True)}{wetlab_rows}</table></section>
<section><h2>Representative Table Explanations</h2><div class="grid">{table_explanations}</div></section>
<section><h2>Representative Evidence Figure Explanations</h2><div class="grid">{extracted_fig_html}</div></section>
<section><h2>Tracked Evidence Pages</h2><div class="grid">{cards}</div></section>
<section><h2>Source Paths</h2><table>{hrow(["File","Use"], True)}{hrow([str(summary_path.relative_to(ROOT)), "Metric JSON for PRISM/DepMap/metabolic headline values"])}{hrow([str(candidate_path.relative_to(ROOT)), "Program-level candidate matrix and disposition"])}{hrow([str(wetlab_path.relative_to(ROOT)), "Concrete validation queue and falsification criteria"])}</table></section>
</div></body></html>"""


def build_deconvolution_repo(pages: list[PageMeta], idx: dict[str, list[str]]) -> str:
    branch = next(b for b in PAPER_BRANCHES if b["branch"] == "paper1-main")
    ps = [p for p in audit_pages(pages) if any(tok in f"{p.path.name} {p.title}".lower() for tok in ["deconv", "spatial", "fig8", "paper1"])]
    cards = "".join(build_repo_page_card(p, idx) for p in sorted(ps, key=lambda p: -review_priority_score(p))[:50])
    rows = "".join(hrow(r) for r in [
        ["Reference", "Lu 2023/GSE193581 author_celltype pseudobulk, 67,678 cells, 8 cell types", "Cell-type composition reference, not patient-level truth."],
        ["Bulk", "TCGA-THCA log2(TPM+1), canonical DM labels and 8-gene RAI score", "Primary Paper 1 cohort."],
        ["Algorithms", "NNLS, Ridge-NNLS, LR-clip, sklearn NuSVR", "nu-SVR primary because NNLS-family collapses T cell to zero."],
        ["Interpretation", "47% nu-SVR retention after full residualization; v13/v14 MAPK/TDS16 lock", "Bulk signal is partly composition, partly retained molecular axis."],
        ["Boundary", "GSE250521 spatial v15B-v18 caveat/autopsy only", "Do not make it the main mechanism proof."],
    ])
    return f"""<!DOCTYPE html><html lang="ko"><head><meta charset="UTF-8"/><meta name="viewport" content="width=device-width,initial-scale=1.0"/><title>Deconvolution Repo 2026-05-10</title>{atlas_style()}</head><body><div class="wrap">
{repo_nav()}
<header class="hero"><div class="kicker">paper1 branch · deconvolution history</div><h1>Deconvolution Repository</h1><p class="subtitle">bulk deconvolution, Fig 8 mechanism, spatial caveat, v13/v14 reviewer defense를 하나의 repo history로 묶은 페이지입니다.</p><div class="stats"><div class="stat"><b>{len(ps)}</b><span>tracked pages</span></div><div class="stat"><b>{sum(len(p.images) for p in ps)}</b><span>figures</span></div><div class="stat"><b>{sum(p.table_count for p in ps)}</b><span>tables</span></div></div></header>
<section><h2>Method / Interpretation Ledger</h2><table>{hrow(["Layer","Data or software","How to interpret"], True)}{rows}</table></section>
<section><h2>Tracked Evidence Pages</h2><div class="grid">{cards}</div></section>
<section><h2>Branch Boundary</h2><div class="warning">{html.escape(branch["boundary"])}</div></section>
</div></body></html>"""


def build_neoantigen_repo(pages: list[PageMeta], idx: dict[str, list[str]]) -> str:
    branch = next(b for b in PAPER_BRANCHES if b["branch"] == "neoantigen-clean-neo")
    ps = branch_pages(branch, pages)
    cards = "".join(build_repo_page_card(p, idx) for p in sorted(ps, key=lambda p: -review_priority_score(p))[:80])
    rows = "".join(hrow(r) for r in [
        ["Benchmark contract", "exact peptide-HLA, near peptide, source-protein, study/patient, HLA/supertype, public-overlap audits", "Use strict split rows before pooled claims."],
        ["Current best rows", "rule_gated_C_QK_structure diagnostic AUPRC up to 0.647; locked internal rows lower", "Diagnostic and locked rows must not be mixed."],
        ["CLEAN-Neo", "delta sequence, presentation, structure, retrieval, TCR latent, patient MIL, calibration/abstention", "Architecture plan; not yet definitive SOTA proof."],
        ["Vaccine triage", "PAAD resected/MRD low-burden plausible; THCA routine vaccine-first not defensible", "Patient-selection-first, not broad vaccine claim."],
    ])
    return f"""<!DOCTYPE html><html lang="ko"><head><meta charset="UTF-8"/><meta name="viewport" content="width=device-width,initial-scale=1.0"/><title>Neoantigen Repo 2026-05-10</title>{atlas_style()}</head><body><div class="wrap">
{repo_nav()}
<header class="hero"><div class="kicker">separate universe · neoantigen/vaccine repo</div><h1>Neoantigen Repository</h1><p class="subtitle">CROSS-Neo, CLEAN-Neo, BAR-Neo, cancer vaccine triage를 THCA Paper 1-4와 분리된 methods/resource repo로 정리합니다.</p><div class="stats"><div class="stat"><b>{len(ps)}</b><span>tracked pages</span></div><div class="stat"><b>{sum(len(p.images) for p in ps)}</b><span>figures</span></div><div class="stat"><b>{sum(p.table_count for p in ps)}</b><span>tables</span></div></div></header>
<section><h2>Method / Interpretation Ledger</h2><table>{hrow(["Layer","Data or software","How to interpret"], True)}{rows}</table></section>
<section><h2>Tracked Evidence Pages</h2><div class="grid">{cards}</div></section>
<section><h2>Branch Boundary</h2><div class="warning">{html.escape(branch["boundary"])}</div></section>
</div></body></html>"""


def build_family_index(pages: list[PageMeta]) -> str:
    groups: dict[str, list[PageMeta]] = {}
    for p in audit_pages(pages):
        groups.setdefault(infer_paper_group(p), []).append(p)
    cards = []
    for family, ps in sorted(groups.items()):
        slug = slugify_family(family)
        status_label, status_boundary, status_next = current_status_for_group(family)
        cards.append(
            '<article class="card">'
            f'<b>{len(ps)} pages · {sum(len(p.images) for p in ps)} figs · {sum(p.table_count for p in ps)} tables</b>'
            f'<strong><a href="family_{html.escape(slug)}_dossier_2026_0509.html">{html.escape(family)}</a></strong>'
            f"<p>{html.escape(status_label)} · {html.escape(status_boundary)} · {html.escape(status_next)}</p>"
            "</article>"
        )
    return f"""<!DOCTYPE html><html lang="ko"><head><meta charset="UTF-8"/><meta name="viewport" content="width=device-width,initial-scale=1.0"/><title>Family Dossiers 2026-05-09</title>{atlas_style()}</head><body><div class="wrap">
<nav class="topnav"><a href="index.html">Renewed hub</a><a href="figure_atlas_2026_0509.html">Figure atlas</a><a href="table_atlas_2026_0509.html">Table atlas</a><a href="flow_atlas_2026_0509.html">Flow atlas</a><a href="metric_atlas_2026_0509.html">Metric atlas</a><a href="source_atlas_2026_0509.html">Source atlas</a><a href="reviewer_tour_2026_0509.html">Reviewer tour</a></nav>
<header class="hero"><div class="kicker">papers_hub_2026_0509 · paper-family dossiers</div><h1>Family Dossiers</h1><p class="subtitle">각 paper family별로 pages, figures, tables, current claim boundary를 묶은 review entry point입니다.</p><div class="stats"><div class="stat"><b>{len(groups)}</b><span>families</span></div><div class="stat"><b>{sum(len(v) for v in groups.values())}</b><span>source pages</span></div><div class="stat"><b>0</b><span>protected prose edits</span></div></div></header>
<section><h2>Dossier Links</h2><div class="grid">{''.join(cards)}</div></section></div></body></html>"""


def write_paper_branch_registry_files(pages: list[PageMeta]) -> None:
    source_pages = audit_pages(pages)
    rows = []
    for branch in PAPER_BRANCHES:
        ps = branch_pages(branch, source_pages)
        rows.append(
            {
                "paper": branch.get("paper", ""),
                "branch": branch.get("branch", ""),
                "slot_type": branch.get("slot_type", ""),
                "topic": branch.get("topic", ""),
                "opened": branch.get("opened", ""),
                "last_update": branch.get("last_update", ""),
                "state": branch.get("state", ""),
                "parent": branch.get("parent", ""),
                "tracked_pages": len(ps),
                "tracked_figures": sum(len(p.images) for p in ps),
                "tracked_tables": sum(p.table_count for p in ps),
                "title": branch.get("title", ""),
                "claim": branch.get("claim", ""),
                "data": branch.get("data", ""),
                "software": branch.get("software", ""),
                "boundary": branch.get("boundary", ""),
                "history": branch.get("history", ""),
                "paths": branch.get("paths", ""),
            }
        )
    (TARGET / "paper_branch_registry_2026_05_10.json").write_text(
        json.dumps(rows, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    cols = list(rows[0].keys()) if rows else []
    lines = ["\t".join(cols)]
    for row in rows:
        lines.append("\t".join(str(row.get(c, "")).replace("\t", " ").replace("\n", " ") for c in cols))
    (TARGET / "paper_branch_registry_2026_05_10.tsv").write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_atlas_pages(pages: list[PageMeta], idx: dict[str, list[str]]) -> None:
    source_pages = audit_pages(pages)
    (TARGET / "figure_atlas_2026_0509.html").write_text(build_figure_atlas(source_pages, idx), encoding="utf-8")
    (TARGET / "table_atlas_2026_0509.html").write_text(build_table_atlas(source_pages), encoding="utf-8")
    (TARGET / "flow_atlas_2026_0509.html").write_text(build_flow_atlas(source_pages), encoding="utf-8")
    (TARGET / "metric_atlas_2026_0509.html").write_text(build_metric_atlas(source_pages), encoding="utf-8")
    (TARGET / "source_atlas_2026_0509.html").write_text(build_source_atlas(source_pages, idx), encoding="utf-8")
    (TARGET / "reviewer_tour_2026_0509.html").write_text(build_reviewer_tour(source_pages), encoding="utf-8")
    (TARGET / "top22_review_pack_2026_0509.html").write_text(build_top22_review_pack(source_pages, idx), encoding="utf-8")
    manifest = build_search_manifest(source_pages, idx)
    (TARGET / "search_manifest_2026_0509.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    (TARGET / "hub_search_2026_0509.html").write_text(build_hub_search_page(manifest), encoding="utf-8")
    (TARGET / "evidence_qc_dashboard_2026_0509.html").write_text(build_evidence_qc_dashboard(source_pages, idx), encoding="utf-8")
    (TARGET / "claim_boundary_ledger_2026_0509.html").write_text(build_claim_boundary_ledger(source_pages), encoding="utf-8")
    (TARGET / "paper1_critical_path_2026_0509.html").write_text(build_paper1_critical_path(source_pages), encoding="utf-8")
    (TARGET / "figure_data_matrix_2026_0509.html").write_text(build_figure_data_matrix(source_pages, idx), encoding="utf-8")
    (TARGET / "paper_repo_control_room_2026_05_10.html").write_text(build_paper_repo_control_room(source_pages, idx), encoding="utf-8")
    (TARGET / "high_impact_board_2026_05_10.html").write_text(build_high_impact_board(source_pages, idx), encoding="utf-8")
    (TARGET / "impact_escalation_matrix_2026_05_10.html").write_text(build_impact_escalation_matrix(source_pages, idx), encoding="utf-8")
    (TARGET / "killer_figure_storyboard_2026_05_10.html").write_text(build_killer_figure_storyboard(source_pages, idx), encoding="utf-8")
    (TARGET / "project_inventory_2026_05_10.html").write_text(build_project_inventory_page(source_pages, idx), encoding="utf-8")
    (TARGET / "paper_repo_map_2026_05_10.html").write_text(build_paper_repo_map(source_pages, idx), encoding="utf-8")
    (TARGET / "paper_branch_manifest_2026_05_10.html").write_text(build_paper_branch_manifest(source_pages, idx), encoding="utf-8")
    (TARGET / "paper_date_topic_index_2026_05_10.html").write_text(build_paper_date_topic_index(source_pages, idx), encoding="utf-8")
    (TARGET / "software_method_ledger_2026_05_10.html").write_text(build_software_method_ledger(source_pages), encoding="utf-8")
    (TARGET / "ic50_prediction_map_2026_05_10.html").write_text(build_ic50_prediction_map(source_pages, idx), encoding="utf-8")
    (TARGET / "deconvolution_repo_2026_05_10.html").write_text(build_deconvolution_repo(source_pages, idx), encoding="utf-8")
    (TARGET / "neoantigen_repo_2026_05_10.html").write_text(build_neoantigen_repo(source_pages, idx), encoding="utf-8")
    for branch in PAPER_BRANCHES:
        (TARGET / branch_dossier_filename(branch)).write_text(build_branch_dossier(branch, source_pages, idx), encoding="utf-8")
    write_paper_branch_registry_files(source_pages)
    write_project_inventory_files()
    groups: dict[str, list[PageMeta]] = {}
    for p in source_pages:
        groups.setdefault(infer_paper_group(p), []).append(p)
    for family, ps in groups.items():
        slug = slugify_family(family)
        (TARGET / f"family_{slug}_dossier_2026_0509.html").write_text(
            build_family_dossier(family, ps, idx), encoding="utf-8"
        )
    (TARGET / "family_index_2026_0509.html").write_text(build_family_index(source_pages), encoding="utf-8")


def build_master_index(pages: list[PageMeta]) -> str:
    # `pages` still contains the transformed source index at this point; this
    # function overwrites it with a generated master index. Account for the
    # generated page shape so the hero stats match the final deployed files.
    audited = audit_pages(pages)
    total_images = sum(len(p.images) for p in audited)
    total_tables = sum(p.table_count for p in audited)
    family_counts: dict[str, int] = {}
    for p in audited:
        family_counts[infer_paper_group(p)] = family_counts.get(infer_paper_group(p), 0) + 1

    cards = []
    for p in pages:
        group = infer_paper_group(p)
        title = p.title or p.path.stem
        terms = ", ".join(t["key"] for t in infer_terms(p)[:4])
        cards.append(
            f'<a class="card" href="{html.escape(p.path.name)}"><b>{html.escape(title[:92])}</b>'
            f'<span>{html.escape(group)} · {len(p.images)} figures · {p.table_count} tables · {len(p.headings)} headings</span>'
            f'<span>{html.escape(terms)}</span></a>'
        )

    family_rows = "\n".join(hrow([k, v]) for k, v in sorted(family_counts.items()))
    status_rows = "\n".join(
        hrow([family, family_counts[family], *current_status_for_group(family)])
        for family in sorted(family_counts)
    )
    definition_rows = "\n".join(
        hrow([t["key"], t["type"], t["definition"], t["use"]]) for t in DATA_TERMS
    )
    pages_by_size = sorted(pages, key=lambda p: p.path.stat().st_size, reverse=True)[:25]
    heavy_rows = "\n".join(
        hrow([p.path.name, f"{p.path.stat().st_size/1024:.1f} KB", len(p.images), p.table_count, infer_paper_group(p)])
        for p in pages_by_size
    )

    return f"""<!DOCTYPE html>
<html lang="ko"><head><meta charset="UTF-8"/><meta name="viewport" content="width=device-width,initial-scale=1.0"/>
<title>THCA Papers Hub 2026-05-09 · Renewed</title>{common_style()}</head><body>
<section class="master-def">
  <div class="kicker">data definition first · applies to every renewed page</div>
  <h2>데이터 정의서 / Data Definition Sheet</h2>
  <p>이 허브는 모든 페이지 앞에 같은 방식의 v5 정의서를 붙였습니다. 아래는 index용 master 정의서이고, 각 하위 페이지에는 그 페이지가 실제 참조하는 cohort/score/statistic subset과 figure/table provenance audit가 별도로 붙어 있습니다.</p>
  <details open><summary>Master term dictionary</summary><table>{hrow(["Term","Class","Definition","Use"], True)}{definition_rows}</table></details>
</section>
<header class="master-hero">
  <div class="kicker">THCA papers hub · renewed 2026-05-09 · 14-slot repo registry</div>
  <h1>14 Paper Slots · All Pages Renewed</h1>
  <p>Paper 1 수준을 전체 허브에 강제 적용한 리뉴얼 빌드입니다. 14개 numbered paper slot과 Neo/CROSS-Neo side repo를 git-style branch registry로 관리하고, 모든 HTML 페이지 첫머리에는 데이터 정의서, 하단에는 figure/table/source audit를 붙였습니다.</p>
  <div class="stats">
    <div class="stat"><b>{len(numbered_paper_branches())}</b><span>paper slots</span></div>
    <div class="stat"><b>{len(side_repo_branches())}</b><span>side repos</span></div>
    <div class="stat"><b>{len(pages)}</b><span>renewed html pages</span></div>
    <div class="stat"><b>{total_images}</b><span>figures audited</span></div>
    <div class="stat"><b>{total_tables}</b><span>tables audited</span></div>
    <div class="stat"><b>{len(DATA_TERMS)}</b><span>definition terms</span></div>
    <div class="stat"><b>{len(family_counts)}</b><span>page families</span></div>
    <div class="stat"><b>v5</b><span>data dictionary</span></div>
  </div>
</header>
<main class="master-wrap">
  <section class="panel">
    <h2>Start Here</h2>
    <p><a href="paper_repo_control_room_2026_05_10.html">Control room</a> · <a href="high_impact_board_2026_05_10.html">High-impact board</a> · <a href="impact_escalation_matrix_2026_05_10.html">Impact matrix</a> · <a href="killer_figure_storyboard_2026_05_10.html">Killer figures</a> · <a href="project_inventory_2026_05_10.html">Project inventory</a> · <a href="ic50_prediction_map_2026_05_10.html">IC50/PRISM map</a> · <a href="paper_repo_map_2026_05_10.html">Paper repo map</a> · <a href="paper_branch_manifest_2026_05_10.html">Branch manifest</a> · <a href="paper_date_topic_index_2026_05_10.html">Date/topic index</a> · <a href="software_method_ledger_2026_05_10.html">Software ledger</a> · <a href="deconvolution_repo_2026_05_10.html">Deconvolution repo</a> · <a href="neoantigen_repo_2026_05_10.html">Neoantigen repo</a> · <a href="hub_search_2026_0509.html">Search</a> · <a href="current_situation_2026_0509.html">Current situation</a> · <a href="top22_review_pack_2026_0509.html">Top 22 review pack</a> · <a href="paper1_critical_path_2026_0509.html">Paper 1 critical path</a> · <a href="evidence_qc_dashboard_2026_0509.html">Evidence QC</a> · <a href="claim_boundary_ledger_2026_0509.html">Claim ledger</a> · <a href="figure_data_matrix_2026_0509.html">Figure data matrix</a> · <a href="reviewer_tour_2026_0509.html">Reviewer tour</a> · <a href="figure_atlas_2026_0509.html">Figure atlas</a> · <a href="table_atlas_2026_0509.html">Table atlas</a> · <a href="flow_atlas_2026_0509.html">Flow atlas</a> · <a href="metric_atlas_2026_0509.html">Metric atlas</a> · <a href="source_atlas_2026_0509.html">Source atlas</a> · <a href="family_index_2026_0509.html">Family dossiers</a> · <a href="data_definitions_2026_0509.html">Data definitions master</a> · <a href="page_inventory_2026_0509.tsv">TSV inventory</a> · <a href="page_inventory_2026_0509.json">JSON inventory</a> · <a href="paper1.html">Paper 1 benchmark page</a></p>
    <table>{hrow(["Family","Renewed pages"], True)}{family_rows}</table>
  </section>
  <section class="panel">
    <h2>Current Claim Boundaries</h2>
    <table>{hrow(["Family","Pages","Current status","Claim boundary","Next use"], True)}{status_rows}</table>
  </section>
  <section class="panel">
    <h2>Heavy Pages First</h2>
    <table>{hrow(["Page","Size","Figures","Tables","Family"], True)}{heavy_rows}</table>
  </section>
  <section class="panel">
    <h2>All Renewed Pages</h2>
    <div class="grid">{''.join(cards)}</div>
  </section>
  <section class="panel">
    <h2>Build Source</h2>
    <p><code>project/scripts/build_papers_hub_2026_0509.py</code> copied <code>project/papers_hub_2026_05_04</code> into <code>project/papers_hub_2026_0509</code>, removed stale v4 definition blocks, injected v5 per-page data definitions, and appended visual/table/source audits. Deployed mirror: <code>/var/www/papers/papers_hub_2026_0509</code>.</p>
  </section>
</main></body></html>"""


def write_inventory(pages: list[PageMeta]) -> None:
    rows = []
    for p in pages:
        rows.append(
            {
                "page": p.path.name,
                "title": p.title,
                "family": infer_paper_group(p),
                "images": len(p.images),
                "tables": p.table_count,
                "headings": len(p.headings),
                "definitions": [t["key"] for t in infer_terms(p)],
                "bytes": p.path.stat().st_size,
            }
        )
    (TARGET / "page_inventory_2026_0509.json").write_text(
        json.dumps(rows, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    tsv = ["page\ttitle\tfamily\timages\ttables\theadings\tdefinitions\tbytes"]
    for r in rows:
        tsv.append(
            "\t".join(
                [
                    r["page"],
                    r["title"].replace("\t", " "),
                    r["family"],
                    str(r["images"]),
                    str(r["tables"]),
                    str(r["headings"]),
                    "; ".join(r["definitions"]),
                    str(r["bytes"]),
                ]
            )
        )
    (TARGET / "page_inventory_2026_0509.tsv").write_text("\n".join(tsv) + "\n", encoding="utf-8")


def deploy_target() -> None:
    if DEPLOY.exists():
        shutil.rmtree(DEPLOY)
    shutil.copytree(TARGET, DEPLOY)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build the renewed 2026-05-09 papers hub.")
    parser.add_argument(
        "--no-deploy",
        action="store_true",
        help="Build project/papers_hub_2026_0509 locally without copying to /var/www.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if not SOURCE.exists():
        raise SystemExit(f"Missing source hub: {SOURCE}")
    copy_source_tree()
    idx = result_inventory()
    pages = collect_pages()
    # Write transformed pages except index; index gets a new master landing page.
    for meta in pages:
        raw = strip_old_blocks(meta.path.read_text(encoding="utf-8", errors="ignore"))
        meta2 = parse_page(meta.path, raw)
        meta.path.write_text(inject_page(raw, meta2, idx, pages), encoding="utf-8")
    pages = collect_pages()
    (TARGET / "data_definitions_2026_0509.html").write_text(build_master_data_dictionary(pages), encoding="utf-8")
    (TARGET / "current_situation_2026_0509.html").write_text(build_current_situation_page(), encoding="utf-8")
    pages = collect_pages()
    write_atlas_pages(pages, idx)
    pages = collect_pages()
    (TARGET / "index.html").write_text(build_master_index(pages), encoding="utf-8")
    pages = collect_pages()
    (TARGET / "current_situation_2026_0509.html").write_text(build_current_situation_page(pages), encoding="utf-8")
    pages = collect_pages()
    write_inventory(pages)
    if args.no_deploy:
        deployed = False
    else:
        deploy_target()
        deployed = True
    audited = audit_pages(pages)
    print(f"built {TARGET}")
    if deployed:
        print(f"deployed {DEPLOY}")
    else:
        print("deploy skipped (--no-deploy)")
    print(f"html pages: {len(pages)}")
    print(f"figures audited: {sum(len(p.images) for p in audited)}")
    print(f"tables audited: {sum(p.table_count for p in audited)}")


if __name__ == "__main__":
    main()
