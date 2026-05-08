"""Offline demo for the v19 research_compass extension.

Runs end-to-end with stub data — no network. Demonstrates:
  CV (markdown text)
    -> ResearchProfile
    -> tier-2 fan-out (journal_feed | data_registry | method_index)
    -> ranked TopicHypothesis cards

Run:
    cd project/v18_agentic && python3 -m examples.research_compass_demo
"""
from __future__ import annotations

import asyncio
from pathlib import Path

from agentic_research.domains.research_compass import (
    DatasetEntry,
    JournalItem,
    ResearchCompass,
)


SAMPLE_CV = """\
# Seungho Cook — CV (excerpt)

## Skills
- RNA-seq analysis (bulk + single-cell)
- HLA imputation from RNA-seq
- Pathology image deep learning (CLAM, attention MIL)
- Cohort harmonization (LODO ComBat)

## Research
- thyroid cancer molecular dark matter (DM1)
- Hashimoto-overlap PTC (PTC+HT)
- Korean cohort HLA forensics

## Methods
- ComBat-seq
- arcasHLA
- CLAM
- DESeq2

## Datasets
- TCGA-THCA
- GSE286332
- PRJEB11591
"""


def _stub_journal() -> list[JournalItem]:
    return [
        JournalItem(
            title="Pan-cancer dedifferentiation axis predicts ICI response",
            journal="Nat Med", date="2026-04-12",
            keywords=["dedifferentiation", "ICI", "thyroid", "pan-cancer", "RNA-seq"],
            source="rss",
        ),
        JournalItem(
            title="HLA class II haplotype distributions in East Asian thyroid cancer",
            journal="BRIC-한빛사", date="2026-03-22",
            keywords=["HLA", "thyroid cancer", "Korean", "class II"],
            source="bric",
        ),
        JournalItem(
            title="Foundation models for H&E whole-slide images: benchmarks",
            journal="Nature", date="2025-11-08",
            keywords=["pathology", "WSI", "foundation model", "H&E", "deep learning"],
            source="rss",
        ),
        JournalItem(
            title="A new soil microbiome atlas of the Andean highlands",
            journal="ISME J", date="2026-02-01",
            keywords=["microbiome", "soil", "16S", "ecology"],
            source="rss",
        ),
    ]


def _stub_datasets() -> list[DatasetEntry]:
    return [
        DatasetEntry(accession="GSE286332", source="GEO", n_samples=18, assay="RNA-seq",
                     keywords=["thyroid", "PTC", "Hashimoto", "Korean"], access="open"),
        DatasetEntry(accession="thca_tcga_pub", source="cBioPortal", n_samples=496, assay="RNA-seq",
                     keywords=["thyroid", "TCGA", "PTC"], access="open"),
        DatasetEntry(accession="PRJEB11591", source="SRA", n_samples=260, assay="RNA-seq",
                     keywords=["thyroid", "Korean", "HLA"], access="open"),
        DatasetEntry(accession="UKB-Showcase", source="UKB", n_samples=500_000, assay="exome",
                     keywords=["GeneBass", "burden", "phenotype"], access="controlled"),
    ]


async def main():
    compass = ResearchCompass(
        live=False,
        stub_journal_items=_stub_journal(),
        stub_datasets=_stub_datasets(),
        top_k=5,
    )
    result = await compass.run(SAMPLE_CV)

    print(f"\nProfile: {result.profile.name}")
    print(f"  skills:        {result.profile.skills}")
    print(f"  prior_topics:  {result.profile.prior_topics}")
    print(f"  methods_used:  {result.profile.methods_used}")
    print(f"  datasets_used: {result.profile.datasets_used}")
    print(f"  keywords:      {result.profile.keywords}")
    print(f"\nFetched: {len(result.journal_items)} papers, {len(result.datasets)} datasets, {len(result.method_index)} methods indexed")
    print(f"Orchestration: {result.orchestration_summary['ok']}/{result.orchestration_summary['n_tasks']} tasks ok")

    print("\n=== Top hypotheses ===")
    for i, h in enumerate(result.hypotheses, 1):
        ds = h.suggested_dataset.accession if h.suggested_dataset else "—"
        mt = h.suggested_method.name if h.suggested_method else "—"
        print(f"\n[{i}] {h.tier} (score={h.score})  {h.title}")
        print(f"     dataset: {ds} | method: {mt}")
        print(f"     why: {h.why_this_researcher}")
        print(f"     components: {h.components}")


if __name__ == "__main__":
    asyncio.run(main())
