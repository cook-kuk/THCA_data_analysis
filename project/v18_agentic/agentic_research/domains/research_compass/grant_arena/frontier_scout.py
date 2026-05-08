"""frontier_scout — extract methods + keywords from top researchers' papers.

Concept: keep a watchlist of high-impact researchers (Hernán, Pearl, Topol,
Mukherjee, Lipton, Hoffman, Aerts, Esteva, ...) and their handles, pull their
recent papers, and extract two distilled signals:

    1. methods    — named techniques (target trial emulation, causal forest,
                    foundation model, IPW, doubly robust, MAE, contrastive,
                    spatial transcriptomics, ...)
    2. keywords   — topic terms (oncology, ICI response, multi-omics,
                    digital twin, surgical AI, dedifferentiation, ...)

These flow back into the arena as priors for FrontierTechAgent, the
method_index, and the BRIC trend agent.

Offline-first. Live wiring (PubMed esearch, arXiv API, bioRxiv,
Semantic Scholar Graph API) is v0.3.
"""
from __future__ import annotations

import re
from collections import Counter
from dataclasses import dataclass, field
from typing import Iterable, Sequence


@dataclass
class ResearcherHandle:
    name: str
    pubmed_id: str = ""        # Author Identifier (NCBI ORCID-linked)
    arxiv_id: str = ""
    semantic_scholar_id: str = ""
    google_scholar_id: str = ""
    notes: str = ""


@dataclass
class ScoutedPaper:
    title: str
    venue: str
    year: int
    abstract: str = ""
    authors: list[str] = field(default_factory=list)
    methods: list[str] = field(default_factory=list)
    keywords: list[str] = field(default_factory=list)
    source: str = ""           # "pubmed" | "arxiv" | "biorxiv" | "semantic_scholar" | "stub"
    url: str = ""


# Lexicons used to distill method/keyword signals from titles+abstracts.
# Hand-curated; ordered roughly by specificity. Extend per domain.
_METHOD_LEXICON = (
    # causal / stats
    "target trial emulation", "causal forest", "doubly robust", "instrumental variable",
    "propensity score", "inverse probability weighting", "g-computation",
    "synthetic control", "regression discontinuity", "difference in differences",
    # ml / dl architectures
    "foundation model", "transformer", "diffusion model", "graph neural network",
    "contrastive learning", "self-supervised", "masked autoencoder",
    "mixture of experts", "state space model", "mamba", "retrieval augmented",
    # genomics / imaging
    "spatial transcriptomics", "single cell", "scRNA-seq", "ATAC-seq",
    "multi-omics", "multiplex immunofluorescence", "digital pathology",
    "whole slide image", "CLAM", "MIL",
    # surgical / physical AI
    "surgical phase recognition", "instrument tracking", "surgical foundation",
    "physical AI", "world model",
    # quantum
    "tensor network", "matrix product state", "quantum kernel",
    "variational quantum", "quantum machine learning",
)

_KEYWORD_LEXICON = (
    "oncology", "thyroid cancer", "papillary thyroid", "active surveillance",
    "dedifferentiation", "RAI resistance", "ICI response", "tumor microenvironment",
    "EMT", "CTC", "liquid biopsy", "digital twin", "counterfactual",
    "biomarker", "drug discovery", "synthetic lethality", "clinical decision",
    "robotic surgery", "endocrine surgery", "pathology", "spatial",
)


def _extract_phrases(text: str, lexicon: Sequence[str]) -> list[str]:
    if not text:
        return []
    t = text.lower()
    found = [p for p in lexicon if p.lower() in t]
    return list(dict.fromkeys(found))  # dedupe, preserve order


def annotate(paper: ScoutedPaper) -> ScoutedPaper:
    """Run lexicon-based method+keyword extraction on a paper in-place."""
    text = f"{paper.title}\n{paper.abstract}"
    if not paper.methods:
        paper.methods = _extract_phrases(text, _METHOD_LEXICON)
    if not paper.keywords:
        paper.keywords = _extract_phrases(text, _KEYWORD_LEXICON)
    return paper


async def scout_researcher(
    handle: ResearcherHandle,
    *,
    days: int = 365,
    live: bool = False,
    stub_papers: Sequence[ScoutedPaper] = (),
) -> list[ScoutedPaper]:
    """Pull recent papers for one researcher and annotate.

    Live mode (v0.3) hits PubMed esearch, arXiv API, Semantic Scholar Graph
    API; offline returns stub_papers (after annotation).
    """
    if not live:
        return [annotate(p) for p in stub_papers]
    raise NotImplementedError("Live researcher scout deferred to v0.3")


async def scout_watchlist(
    handles: Sequence[ResearcherHandle],
    *,
    days: int = 365,
    live: bool = False,
    stub_by_name: dict[str, Sequence[ScoutedPaper]] | None = None,
) -> list[ScoutedPaper]:
    """Run scout_researcher across a whole watchlist; flatten + annotate."""
    out: list[ScoutedPaper] = []
    stubs = stub_by_name or {}
    for h in handles:
        out.extend(await scout_researcher(
            h, days=days, live=live, stub_papers=stubs.get(h.name, ()),
        ))
    return out


def aggregate_signals(papers: Sequence[ScoutedPaper]) -> dict[str, Counter]:
    """Collapse a paper list into method / keyword / venue frequency tables.

    Use this to update FrontierTechAgent priors or bump method_index entries.
    """
    methods: Counter[str] = Counter()
    keywords: Counter[str] = Counter()
    venues: Counter[str] = Counter()
    for p in papers:
        methods.update(p.methods)
        keywords.update(p.keywords)
        if p.venue:
            venues[p.venue] += 1
    return {"methods": methods, "keywords": keywords, "venues": venues}


def default_watchlist() -> list[ResearcherHandle]:
    """Seed list of high-leverage researchers for biomedical+causal+AI work."""
    return [
        ResearcherHandle("Miguel Hernán", notes="target trial emulation, causal inference"),
        ResearcherHandle("Judea Pearl", notes="causal hierarchy, do-calculus"),
        ResearcherHandle("Eric Topol", notes="medical AI, multimodal foundation"),
        ResearcherHandle("Pranav Rajpurkar", notes="medical foundation models"),
        ResearcherHandle("James Zou", notes="ML for biomedicine, fairness"),
        ResearcherHandle("Hugo Aerts", notes="radiomics, oncology AI"),
        ResearcherHandle("Faisal Mahmood", notes="digital pathology foundation"),
        ResearcherHandle("Andrew Beck", notes="computational pathology"),
        ResearcherHandle("Daphne Koller", notes="ML, Insitro, drug discovery"),
        ResearcherHandle("Regina Barzilay", notes="biomedical NLP, drug design"),
    ]


def merge_into_method_index_seed(
    aggregated: Counter,
    *,
    min_count: int = 2,
) -> list[dict]:
    """Convert top-N method counts into MethodEntry-shaped dicts.

    Caller can pass these into MethodIndex.extend([MethodEntry(**d) for d in seeds]).
    """
    seeds = []
    for name, count in aggregated.most_common():
        if count < min_count:
            break
        seeds.append({
            "name": name,
            "paper_title": f"frontier-scout aggregate (n={count})",
            "year": 0,
            "domains": [],
            "keywords": [name],
        })
    return seeds


def _is_recent(paper: ScoutedPaper, current_year: int, window: int = 3) -> bool:
    return paper.year >= current_year - window


def hot_methods(papers: Sequence[ScoutedPaper], *, current_year: int, top_k: int = 10) -> list[tuple[str, int]]:
    """Methods that show up most in the last `window` years across the watchlist."""
    recent = [p for p in papers if _is_recent(p, current_year)]
    counts: Counter[str] = Counter()
    for p in recent:
        counts.update(p.methods)
    return counts.most_common(top_k)
