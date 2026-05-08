"""Smoke tests for the v19 research_compass extension."""
from __future__ import annotations

import asyncio

import pytest

from agentic_research.domains.research_compass import (
    DatasetEntry,
    JournalItem,
    MethodIndex,
    ResearchCompass,
    ResearchProfile,
    load_cv,
    rank_topics,
)
from agentic_research.domains.research_compass.method_index import default_seed_entries


_CV_TEXT = """\
# Test User

## Skills
- RNA-seq
- HLA imputation

## Datasets
- TCGA-THCA
- GSE286332
"""


def test_load_cv_parses_sections():
    p = load_cv(_CV_TEXT, name="tester")
    assert p.name == "tester"
    assert "RNA-seq" in p.skills
    assert "TCGA-THCA" in p.datasets_used
    assert "GSE286332" in p.keywords  # picked up by GSE regex


def test_method_index_search_ranks_overlap():
    idx = MethodIndex(default_seed_entries())
    hits = idx.search(["RNA-seq", "batch effect"], top_k=3)
    assert hits, "expected at least one hit"
    top_name = hits[0][0].name
    assert top_name in {"ComBat-seq", "DESeq2"}, f"unexpected top method: {top_name}"


def test_rank_topics_orders_by_score_and_tags_tiers():
    profile = ResearchProfile(
        name="t",
        skills=["RNA-seq"],
        prior_topics=["thyroid"],
        methods_used=["ComBat-seq"],
        datasets_used=["TCGA-THCA"],
        keywords=["GSE286332"],
    )
    items = [
        JournalItem(title="thyroid RNA-seq paper", journal="X", date="2026-04-01",
                    keywords=["thyroid", "RNA-seq", "PTC"]),
        JournalItem(title="microbiome paper", journal="Y", date="2026-03-01",
                    keywords=["microbiome", "soil"]),
    ]
    ds = [
        DatasetEntry(accession="GSE286332", source="GEO", assay="RNA-seq",
                     keywords=["thyroid", "PTC"], access="open"),
        DatasetEntry(accession="UKB", source="UKB", assay="exome",
                     keywords=["burden"], access="controlled"),
    ]
    idx = MethodIndex(default_seed_entries())
    cards = rank_topics(profile, items, idx, ds, top_k=2)
    assert len(cards) == 2
    assert cards[0].title.startswith("thyroid"), "thyroid item should outrank microbiome"
    assert cards[0].score > cards[1].score
    assert cards[0].tier in {"HIGH", "MEDIUM", "LOW"}


def test_pipeline_offline_runs_end_to_end():
    compass = ResearchCompass(
        live=False,
        stub_journal_items=[
            JournalItem(title="thyroid foundation model", journal="Nature", date="2026-04-01",
                        keywords=["thyroid", "foundation model", "RNA-seq"]),
        ],
        stub_datasets=[
            DatasetEntry(accession="GSE286332", source="GEO", assay="RNA-seq",
                         keywords=["thyroid"], access="open"),
        ],
        top_k=3,
    )
    result = asyncio.run(compass.run(_CV_TEXT))
    assert result.profile.name
    assert result.orchestration_summary["ok"] == 3
    assert len(result.hypotheses) >= 1
    assert result.hypotheses[0].seed_papers
