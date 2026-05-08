"""Five-stage orchestration that reuses v18's TieredParallelOrchestrator.

Tier 1: cv_profile (single, sequential — feeds tier 2 query terms)
Tier 2: journal_feed | data_registry | method_index seeding (parallel fan-out)
Tier 3: topic_ranker (joins tier-2 outputs)

This is the same Phase-B shape used in the v17p35 thyroid sprint, just
applied to topic discovery instead of biology validation.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Sequence

from ...patterns.tiered_dag import Tier, TieredParallelOrchestrator
from ...core.orchestrator import Task

from .cv_profile import ResearchProfile, load_cv
from .data_registry import DatasetEntry, query_datasets
from .journal_feed import JournalItem, fetch_journal_feed
from .method_index import MethodEntry, MethodIndex, default_seed_entries
from .topic_ranker import TopicHypothesis, rank_topics


@dataclass
class CompassResult:
    profile: ResearchProfile
    journal_items: list[JournalItem] = field(default_factory=list)
    datasets: list[DatasetEntry] = field(default_factory=list)
    method_index: MethodIndex | None = None
    hypotheses: list[TopicHypothesis] = field(default_factory=list)
    orchestration_summary: dict | None = None


class ResearchCompass:
    """End-to-end pipeline: CV in, ranked TopicHypothesis cards out."""

    def __init__(
        self,
        *,
        live: bool = False,
        method_entries: Sequence[MethodEntry] | None = None,
        stub_journal_items: Sequence[JournalItem] | None = None,
        stub_datasets: Sequence[DatasetEntry] | None = None,
        days_lookback: int = 90,
        top_k: int = 10,
    ):
        self.live = live
        self.method_entries = list(method_entries) if method_entries else default_seed_entries()
        self.stub_journal_items = list(stub_journal_items or [])
        self.stub_datasets = list(stub_datasets or [])
        self.days_lookback = days_lookback
        self.top_k = top_k

    async def run(self, cv: str) -> CompassResult:
        profile = load_cv(cv)
        terms = sorted(profile.all_terms())

        async def _journal():
            return await fetch_journal_feed(
                days=self.days_lookback,
                query=" OR ".join(terms[:8]),
                live=self.live,
                stub_items=self.stub_journal_items,
            )

        async def _data():
            return await query_datasets(
                terms=terms,
                live=self.live,
                stub_entries=self.stub_datasets,
            )

        async def _methods():
            return MethodIndex(self.method_entries)

        orch = TieredParallelOrchestrator([
            Tier(name="fetch", tasks=[
                Task("journal_feed", _journal),
                Task("data_registry", _data),
                Task("method_index", _methods),
            ]),
        ])
        summary = await orch.run()

        # TieredParallelOrchestrator drops raw payloads from its summary, so we
        # re-await the three coroutines to grab the actual values. Cheap in
        # offline mode (stubs); in live mode, swap to a custom DAG that keeps
        # results, or memoize the fetchers.
        journal_items = await _journal()
        datasets = await _data()
        method_index = await _methods()

        hypotheses = rank_topics(
            profile=profile,
            journal_items=journal_items,
            method_index=method_index,
            datasets=datasets,
            top_k=self.top_k,
        )

        return CompassResult(
            profile=profile,
            journal_items=journal_items,
            datasets=datasets,
            method_index=method_index,
            hypotheses=hypotheses,
            orchestration_summary=summary,
        )
