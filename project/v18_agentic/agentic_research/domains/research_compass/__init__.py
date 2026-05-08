"""research_compass — CV-conditioned topic recommendation for working researchers.

Five-stage pipeline (phylo.bio composite-scoring × feynman.is multi-agent ×
v18 TieredParallelOrchestrator):

  1. cv_profile     researcher CV  -> ResearchProfile (skills, prior topics, methods, datasets)
  2. journal_feed   recent papers  -> JournalItem[]   (BRIC 한빛사, bioRxiv, top-tier RSS)
  3. method_index   method papers  -> MethodEntry[]   (keyword-indexed)
  4. data_registry  open data      -> DatasetEntry[]  (GEO/SRA/cBioPortal/UKB metadata)
  5. topic_ranker   join + score   -> TopicHypothesis[] (genetics-driven ranking, phylo style)

Offline-first: every stage accepts injected stub data so demos and tests run
without network. Live fetchers are behind ``live=True`` flags.
"""
from __future__ import annotations

from .cv_profile import ResearchProfile, load_cv
from .journal_feed import JournalItem, fetch_journal_feed
from .method_index import MethodEntry, MethodIndex
from .data_registry import DatasetEntry, query_datasets
from .topic_ranker import TopicHypothesis, rank_topics
from .pipeline import ResearchCompass

__all__ = [
    "ResearchProfile",
    "load_cv",
    "JournalItem",
    "fetch_journal_feed",
    "MethodEntry",
    "MethodIndex",
    "DatasetEntry",
    "query_datasets",
    "TopicHypothesis",
    "rank_topics",
    "ResearchCompass",
]
