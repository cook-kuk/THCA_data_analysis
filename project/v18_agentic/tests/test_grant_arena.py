"""Smoke tests for grant_arena."""
from __future__ import annotations

import asyncio

from agentic_research.domains.research_compass import (
    DatasetEntry, MethodIndex, ResearchProfile,
)
from agentic_research.domains.research_compass.method_index import default_seed_entries
from agentic_research.domains.research_compass.grant_arena import (
    GrantArena, parse_rfp, score_topic, screen_collision_risk,
    score_dataset_reuse, is_jcr_q1, is_top_ai_venue, BRICItem,
)
from agentic_research.domains.research_compass.grant_arena.frontier_scout import (
    ScoutedPaper, aggregate_signals, annotate, hot_methods,
)


def test_parse_rfp_extracts_keywords_and_budget():
    rfp = parse_rfp("Advanced AI, 디지털 헬스, 양자, 30억 3년형")
    assert rfp.budget_scale == "30억"
    assert rfp.duration_years == 3
    assert "Advanced AI" in rfp.preferred_keywords
    assert "디지털 헬스" in rfp.preferred_keywords
    assert "양자" in rfp.preferred_keywords


def test_scoring_10_factor_and_tiering():
    s = score_topic(
        novelty=0.9, feasibility=0.7, pi_fit=0.85, funder_fit=0.8,
        publishability=0.75, patentability=0.7, clinical_impact=0.7,
        scalability=0.6, pilot_speed=0.6,
        patent_collision_risk=0.1, overclaim_risk=0.1,
    )
    assert s.total > 0.55
    assert s.tier == "HIGH"
    assert "novelty" in s.components


def test_biodata_reuse_score_known_dataset():
    s = score_dataset_reuse("TCGA", disease_relevance=0.9)
    assert 0 < s.composite <= 1.0
    assert s.accession == "TCGA"


def test_jcr_and_top_venue_filters():
    # Default floor is "TOP" (NC and above) per user preference
    assert is_jcr_q1("Nature Medicine")
    assert is_jcr_q1("nature communications")
    assert is_jcr_q1("Nature Medicine").tier == "TOP"
    # Below-NC Q1 venues are filtered OUT by default
    assert is_jcr_q1("Bioinformatics") is None
    assert is_jcr_q1("Bioinformatics", min_tier="Q1") is not None
    # Predatory / unknown stay rejected at any tier
    assert is_jcr_q1("Some Predatory Journal", min_tier="Q1") is None
    assert is_top_ai_venue("NeurIPS")
    assert is_top_ai_venue("icml")
    assert not is_top_ai_venue("Random Workshop")


def test_patent_priorart_returns_hits_for_overlapping_topic():
    hits = screen_collision_risk(
        "Oncology digital twin for thyroid",
        ["digital twin", "oncology"],
    )
    assert hits, "expected at least one prior-art hit from stub corpus"
    assert hits[0].similarity > 0


def test_frontier_scout_extracts_methods_and_aggregates():
    p = ScoutedPaper(
        title="Target trial emulation in oncology",
        venue="JAMA Oncology", year=2024,
        abstract="We review target trial emulation and propensity score methods.",
    )
    annotate(p)
    assert "target trial emulation" in p.methods
    assert "propensity score" in p.methods
    assert "oncology" in p.keywords
    sigs = aggregate_signals([p])
    assert sigs["methods"]["target trial emulation"] == 1


def test_frontier_scout_hot_methods_recent_window():
    papers = [
        ScoutedPaper(title="foundation model paper", venue="Nature", year=2025,
                     methods=["foundation model"]),
        ScoutedPaper(title="old paper", venue="Old Journal", year=2010,
                     methods=["foundation model"]),
    ]
    hot = hot_methods(papers, current_year=2026, top_k=5)
    assert hot, "expected at least one hot method"
    names = [m for m, _ in hot]
    assert "foundation model" in names


def test_arena_runs_end_to_end_offline():
    profile = ResearchProfile(
        name="t",
        skills=["agentic AI", "spatial transcriptomics"],
        prior_topics=["thyroid cancer"],
        methods_used=["ComBat-seq", "target trial emulation"],
        datasets_used=["TCGA-THCA"],
    )
    rfp = parse_rfp("Advanced AI, 디지털 헬스, 정밀의료, 30억 3년")
    datasets = [
        DatasetEntry(accession="TCGA-THCA", source="GDC", assay="RNA-seq",
                     keywords=["thyroid", "PTC"], access="open"),
        DatasetEntry(accession="GSE286332", source="GEO", assay="RNA-seq",
                     keywords=["thyroid", "Korean"], access="open"),
    ]
    method_index = MethodIndex(default_seed_entries())
    bric = [BRICItem(title="Korean PTC", korean_pi="Yu", affiliation="SNUBH",
                     journal="Thyroid", year=2025, keywords=["thyroid", "active surveillance"])]

    arena = GrantArena(topics_per_generator=2, evolve_top_k=2)
    result = arena.run(profile=profile, rfp=rfp, datasets=datasets,
                       method_index=method_index, bric_items=bric)

    assert result.topics
    assert result.ranking
    assert result.final_scores
    top = result.ranking[0]
    assert top.title in result.final_scores
    assert result.final_scores[top.title].total >= result.final_scores[result.ranking[-1].title].total
    assert result.battle_log, "expected at least one Elo battle"
