"""Stage 5 — composite-scored TopicHypothesis ranking.

Phylo.bio uses a 4-factor weighted score for genetic targets. We adapt that
shape for research-topic recommendation:

    score = w_competence * competence_overlap        # CV ↔ method+dataset terms
          + w_novelty    * recency_signal            # journal-feed freshness
          + w_feasibility* dataset_access            # open vs controlled
          + w_method_fit * method_index_overlap      # CV methods ↔ matched methods

Defaults: 0.35 / 0.25 / 0.25 / 0.15 (mirrors phylo's 35/25/25/15 split).

Outputs ``TopicHypothesis`` cards: title, why-this-researcher, suggested
dataset, suggested method, citations, tier (HIGH/MEDIUM/LOW).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Sequence

from .cv_profile import ResearchProfile
from .data_registry import DatasetEntry
from .journal_feed import JournalItem
from .method_index import MethodEntry, MethodIndex


@dataclass
class TopicHypothesis:
    title: str
    tier: str                          # "HIGH" | "MEDIUM" | "LOW"
    score: float
    why_this_researcher: str
    suggested_dataset: DatasetEntry | None = None
    suggested_method: MethodEntry | None = None
    seed_papers: list[JournalItem] = field(default_factory=list)
    components: dict[str, float] = field(default_factory=dict)


_DEFAULT_WEIGHTS = {
    "competence":  0.35,
    "novelty":     0.25,
    "feasibility": 0.25,
    "method_fit":  0.15,
}


def _overlap(a: set[str], b: set[str]) -> float:
    if not a or not b:
        return 0.0
    return len(a & b) / max(1, len(a))


def _tier_for(score: float) -> str:
    if score >= 0.55:
        return "HIGH"
    if score >= 0.30:
        return "MEDIUM"
    return "LOW"


def rank_topics(
    profile: ResearchProfile,
    journal_items: Sequence[JournalItem],
    method_index: MethodIndex,
    datasets: Sequence[DatasetEntry],
    *,
    weights: dict[str, float] | None = None,
    top_k: int = 10,
) -> list[TopicHypothesis]:
    """Join the four signals and produce ranked TopicHypothesis cards.

    One card per (journal_item × best_dataset × best_method) join, scored by
    the four-component weighted sum.
    """
    w = {**_DEFAULT_WEIGHTS, **(weights or {})}
    profile_terms = profile.all_terms()
    cards: list[TopicHypothesis] = []

    for j in journal_items:
        j_terms = {k.lower() for k in j.keywords if k}
        competence = _overlap(profile_terms, j_terms)

        ds_scored = [(d, _overlap(j_terms | profile_terms, d.all_terms())) for d in datasets]
        ds_scored.sort(key=lambda x: x[1], reverse=True)
        best_ds, best_ds_score = ds_scored[0] if ds_scored else (None, 0.0)
        feasibility = best_ds_score * (1.0 if (best_ds and best_ds.access == "open") else 0.5)

        method_hits = method_index.search(profile.methods_used + list(j_terms), top_k=1)
        best_method, method_fit = method_hits[0] if method_hits else (None, 0.0)

        novelty = 1.0 if j.date >= "2024-01-01" else 0.5

        score = (
            w["competence"]  * competence
            + w["novelty"]     * novelty
            + w["feasibility"] * feasibility
            + w["method_fit"]  * method_fit
        )

        why = (
            f"Overlap {len(profile_terms & j_terms)} terms with your CV; "
            f"dataset {best_ds.accession if best_ds else 'n/a'} matches {best_ds_score:.2f}; "
            f"method {best_method.name if best_method else 'n/a'} fit {method_fit:.2f}"
        )

        cards.append(TopicHypothesis(
            title=j.title,
            tier=_tier_for(score),
            score=round(score, 3),
            why_this_researcher=why,
            suggested_dataset=best_ds,
            suggested_method=best_method,
            seed_papers=[j],
            components={"competence": round(competence, 3), "novelty": round(novelty, 3),
                        "feasibility": round(feasibility, 3), "method_fit": round(method_fit, 3)},
        ))

    cards.sort(key=lambda c: c.score, reverse=True)
    return cards[:top_k]
