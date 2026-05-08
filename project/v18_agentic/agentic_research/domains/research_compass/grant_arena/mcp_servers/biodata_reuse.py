"""biodata_reuse_mcp — score how heavily a public biomedical dataset is reused.

Scoring follows the brainstorm spec:

    BioData Reuse Score
      = 0.30 × log(1 + dataset_citing_papers)
      + 0.20 × high_impact_journal_usage
      + 0.15 × top_AI_conference_usage
      + 0.15 × disease_relevance
      + 0.10 × multimodal_value
      + 0.10 × accessibility/reproducibility

For the offline scaffold we hand-curate ~10 canonical datasets with rough
estimates from public sources. v0.3 wires PubMed/Semantic-Scholar mention
counts via live MCP tools.
"""
from __future__ import annotations

import math
from dataclasses import dataclass


@dataclass
class DatasetReuseScore:
    accession: str
    raw_citations: int
    high_impact_usage: float       # 0..1
    top_ai_usage: float            # 0..1
    disease_relevance: float       # 0..1 (set per-query)
    multimodal_value: float        # 0..1
    accessibility: float           # 0..1
    composite: float = 0.0

    def recompute(self) -> "DatasetReuseScore":
        cite_term = math.log1p(self.raw_citations) / math.log1p(20000)  # normalize
        self.composite = round(
            0.30 * cite_term
            + 0.20 * self.high_impact_usage
            + 0.15 * self.top_ai_usage
            + 0.15 * self.disease_relevance
            + 0.10 * self.multimodal_value
            + 0.10 * self.accessibility,
            3,
        )
        return self


def default_registry() -> dict[str, DatasetReuseScore]:
    """Curated seed registry. Citation counts are order-of-magnitude estimates."""
    seed = [
        DatasetReuseScore("TCGA",          15000, 1.0, 0.9, 0.8, 1.0, 1.0),
        DatasetReuseScore("GEO",           50000, 0.9, 0.7, 0.7, 0.9, 1.0),
        DatasetReuseScore("GTEx",           7000, 1.0, 0.8, 0.6, 0.7, 1.0),
        DatasetReuseScore("DepMap",         3000, 1.0, 0.9, 0.9, 0.8, 1.0),
        DatasetReuseScore("CPTAC",          2000, 1.0, 0.7, 0.8, 0.9, 0.9),
        DatasetReuseScore("TCIA",           4000, 0.8, 0.7, 0.7, 0.8, 1.0),
        DatasetReuseScore("HTAN",            500, 1.0, 0.8, 0.9, 1.0, 0.9),
        DatasetReuseScore("CELLxGENE",       800, 0.9, 0.9, 0.7, 0.9, 1.0),
        DatasetReuseScore("LINCS-L1000",    1500, 0.8, 0.9, 0.7, 0.8, 1.0),
        DatasetReuseScore("UKBiobank",     18000, 1.0, 0.8, 0.6, 0.9, 0.5),
        DatasetReuseScore("MIMIC-IV",       6000, 0.9, 1.0, 0.4, 0.7, 0.7),
    ]
    return {d.accession: d.recompute() for d in seed}


def score_dataset_reuse(accession: str, *, disease_relevance: float = 0.5) -> DatasetReuseScore:
    """Look up a dataset, override disease_relevance for the current query."""
    reg = default_registry()
    if accession not in reg:
        return DatasetReuseScore(
            accession=accession, raw_citations=0,
            high_impact_usage=0.3, top_ai_usage=0.2,
            disease_relevance=disease_relevance,
            multimodal_value=0.3, accessibility=0.5,
        ).recompute()
    d = reg[accession]
    d.disease_relevance = disease_relevance
    return d.recompute()
