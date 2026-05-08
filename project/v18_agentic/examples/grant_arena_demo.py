"""Offline demo for grant_arena.

Feeds Yu-Cook stub CV + a tiny RFP + stub datasets + stub BRIC items + stub
frontier-scout papers into the arena. Prints the top-5 ranked grant topics
and their critique log.

Run:
    cd project/v18_agentic && python3 -m examples.grant_arena_demo
"""
from __future__ import annotations

import asyncio

from agentic_research.domains.research_compass import (
    DatasetEntry, JournalItem, MethodIndex, ResearchProfile, load_cv,
)
from agentic_research.domains.research_compass.method_index import default_seed_entries
from agentic_research.domains.research_compass.grant_arena import (
    GrantArena, parse_rfp,
)
from agentic_research.domains.research_compass.grant_arena.frontier_scout import (
    ResearcherHandle, ScoutedPaper, aggregate_signals, default_watchlist,
    merge_into_method_index_seed, scout_watchlist,
)
from agentic_research.domains.research_compass.grant_arena.mcp_servers.bric_trend import BRICItem
from agentic_research.domains.research_compass.method_index import MethodEntry


SAMPLE_CV = """\
# Seungho Cook

## Skills
- agentic AI
- multi-omics integration
- spatial transcriptomics
- digital pathology
- quantum machine learning

## Research
- thyroid cancer molecular dark matter
- papillary thyroid cancer active surveillance
- HLA immunogenomics

## Methods
- ComBat-seq
- arcasHLA
- CLAM
- target trial emulation

## Datasets
- TCGA-THCA
- GSE286332
- PRJEB11591
"""

SAMPLE_RFP = """\
삼성미래기술육성사업 Technology 트랙 30억 3년형.
지정테마: Advanced AI, 차세대 로봇, 디지털 헬스, 정밀의료, 바이오 융합, 양자.
평가기준: 창의성, 원천성, 도전성, 차별성, 실행가능성.
"""


def _stub_datasets() -> list[DatasetEntry]:
    return [
        DatasetEntry(accession="TCGA-THCA", source="GDC", n_samples=496, assay="RNA-seq",
                     keywords=["thyroid", "PTC", "TCGA"], access="open"),
        DatasetEntry(accession="GSE286332", source="GEO", n_samples=18, assay="RNA-seq",
                     keywords=["thyroid", "Hashimoto", "Korean"], access="open"),
        DatasetEntry(accession="DepMap", source="DepMap", n_samples=1100, assay="CRISPR",
                     keywords=["dependency", "synthetic lethality"], access="open"),
        DatasetEntry(accession="UKB-Showcase", source="UKB", n_samples=500_000, assay="exome",
                     keywords=["GeneBass", "burden"], access="controlled"),
    ]


def _stub_bric() -> list[BRICItem]:
    return [
        BRICItem(title="Spatial atlas of Korean PTC progression",
                 korean_pi="Kim, J.", affiliation="SNU", journal="Nat Commun", year=2025,
                 keywords=["thyroid cancer", "spatial transcriptomics"]),
        BRICItem(title="Active surveillance long-term outcomes (MAeSTro-EXP)",
                 korean_pi="Yu, H.W.", affiliation="SNUBH", journal="Thyroid", year=2025,
                 keywords=["active surveillance", "papillary thyroid", "Korean cohort"]),
    ]


def _stub_scout_papers() -> dict[str, list[ScoutedPaper]]:
    return {
        "Miguel Hernán": [
            ScoutedPaper(
                title="Target trial emulation in oncology: a tutorial",
                venue="JAMA Oncology", year=2024,
                abstract="We review target trial emulation and propensity score methods for "
                         "comparative effectiveness in oncology.",
            ),
        ],
        "Eric Topol": [
            ScoutedPaper(
                title="Foundation models for medicine: opportunities and risks",
                venue="Nature Medicine", year=2025,
                abstract="A foundation model approach for digital pathology and multi-omics integration "
                         "with self-supervised contrastive learning.",
            ),
        ],
        "Faisal Mahmood": [
            ScoutedPaper(
                title="A pan-cancer foundation model on whole slide images",
                venue="Nature", year=2025,
                abstract="We pretrain a vision transformer with masked autoencoder and contrastive "
                         "learning on millions of whole slide image tiles.",
            ),
        ],
    }


async def main():
    profile = load_cv(SAMPLE_CV)
    rfp = parse_rfp(SAMPLE_RFP, program="삼성미래기술육성사업", track="Technology")
    datasets = _stub_datasets()
    bric_items = _stub_bric()

    handles = [h for h in default_watchlist() if h.name in {"Miguel Hernán", "Eric Topol", "Faisal Mahmood"}]
    scouted = await scout_watchlist(handles, live=False, stub_by_name=_stub_scout_papers())
    signals = aggregate_signals(scouted)

    method_index = MethodIndex(default_seed_entries())
    new_seeds = merge_into_method_index_seed(signals["methods"], min_count=1)
    method_index.extend([MethodEntry(**s) for s in new_seeds])

    arena = GrantArena(topics_per_generator=2, evolve_top_k=3)
    result = arena.run(
        profile=profile, rfp=rfp, datasets=datasets,
        method_index=method_index, bric_items=bric_items,
    )

    print(f"\nProfile: {profile.name}  | skills={len(profile.skills)} | methods={len(profile.methods_used)}")
    print(f"RFP: {rfp.program} {rfp.track} {rfp.budget_scale} | preferred={rfp.preferred_keywords}")
    print(f"Frontier-scout: {len(scouted)} papers | top methods: {signals['methods'].most_common(5)}")
    print(f"Method index size: {len(method_index)} (seed + {len(new_seeds)} from scout)")
    print(f"Arena: {len(result.topics)} topics, {len(result.battle_log)} battles")

    print("\n=== Top 5 ranked grant topics ===")
    for i, t in enumerate(result.ranking[:5], 1):
        sc = result.final_scores[t.title]
        print(f"\n[{i}] {sc.tier} score={sc.total} elo={t.elo:.0f}  by={t.proposed_by}")
        print(f"    {t.title}")
        print(f"    keywords: {t.keywords}")
        crits = result.critiques.get(t.title, [])
        if crits:
            print(f"    critiques: {len(crits)}")
            for c in crits[:2]:
                print(f"      - [{c.by}] {c.attack}")


if __name__ == "__main__":
    asyncio.run(main())
