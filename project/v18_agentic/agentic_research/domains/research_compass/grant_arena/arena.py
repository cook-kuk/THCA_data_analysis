"""Competitive arena: generate → battle → critique → evolve → judge.

Workflow (deterministic, offline):

    1. Each GeneratorAgent emits N candidate topics.
    2. Topics enter a pairwise Elo tournament; "winner" of A vs B is the
       one with the higher composite score (from scoring.score_topic).
    3. CriticAgents apply factor-delta penalties/bonuses in sequence.
    4. EvolutionAgent (built-in) mutates the top-K topics by combining
       keywords/datasets/methods across the winners.
    5. Final ranking by composite score; HIGH/MEDIUM/LOW tiers.

Pure Python; no LLM calls. v0.3 swaps the judge and evolution mutator for
LLM-backed implementations.
"""
from __future__ import annotations

import itertools
from dataclasses import dataclass, field
from typing import Sequence

from ..cv_profile import ResearchProfile
from ..data_registry import DatasetEntry
from ..method_index import MethodIndex
from .agents.critics import (
    BudgetCritic, CriticAgent, Critique, FeasibilityCritic,
    NoveltyCritic, PatentCritic, Reviewer2Agent,
)
from .agents.generators import (
    BioDataReuseAgent, CancerBiologyAgent, ClinicalNeedAgent, CrazyIdeaAgent,
    FrontierTechAgent, FundingStrategyAgent, GeneratorAgent, GrantTopic,
    KoreanTrendAgent, PatentWhiteSpaceAgent,
)
from .mcp_servers.bric_trend import BRICItem
from .mcp_servers.rfp_reader import RFP
from .scoring import OpportunityScore, ScoringWeights, score_topic


@dataclass
class BattleLog:
    pair: tuple[str, str]                    # (topic_a.title, topic_b.title)
    winner: str
    scores: tuple[float, float]


@dataclass
class ArenaResult:
    topics: list[GrantTopic] = field(default_factory=list)
    final_scores: dict[str, OpportunityScore] = field(default_factory=dict)
    battle_log: list[BattleLog] = field(default_factory=list)
    critiques: dict[str, list[Critique]] = field(default_factory=dict)
    ranking: list[GrantTopic] = field(default_factory=list)


def _topic_score(topic: GrantTopic, weights: ScoringWeights | None = None) -> OpportunityScore:
    return score_topic(weights=weights, **topic.factor_inputs)


def _elo_update(a: GrantTopic, b: GrantTopic, a_won: bool, k: float = 32.0) -> None:
    expected_a = 1.0 / (1.0 + 10 ** ((b.elo - a.elo) / 400))
    score_a = 1.0 if a_won else 0.0
    a.elo += k * (score_a - expected_a)
    b.elo += k * ((1 - score_a) - (1 - expected_a))


class GrantArena:
    """End-to-end orchestrator for the competitive grant-topic search."""

    def __init__(
        self,
        *,
        generators: Sequence[GeneratorAgent] | None = None,
        critics: Sequence[CriticAgent] | None = None,
        weights: ScoringWeights | None = None,
        topics_per_generator: int = 5,
        evolve_top_k: int = 3,
    ):
        self.generators = list(generators) if generators is not None else self._default_generators()
        self.critics = list(critics) if critics is not None else self._default_critics()
        self.weights = weights
        self.topics_per_generator = topics_per_generator
        self.evolve_top_k = evolve_top_k

    @staticmethod
    def _default_generators() -> list[GeneratorAgent]:
        return [
            ClinicalNeedAgent(), FrontierTechAgent(), CancerBiologyAgent(),
            BioDataReuseAgent(), KoreanTrendAgent(), PatentWhiteSpaceAgent(),
            FundingStrategyAgent(), CrazyIdeaAgent(),
        ]

    @staticmethod
    def _default_critics() -> list[CriticAgent]:
        return [FeasibilityCritic(), NoveltyCritic(), PatentCritic(),
                Reviewer2Agent(), BudgetCritic()]

    def run(
        self,
        *,
        profile: ResearchProfile,
        rfp: RFP,
        datasets: Sequence[DatasetEntry],
        method_index: MethodIndex,
        bric_items: Sequence[BRICItem] = (),
    ) -> ArenaResult:
        result = ArenaResult()

        for gen in self.generators:
            result.topics.extend(gen.generate(
                profile=profile, rfp=rfp, datasets=datasets,
                method_index=method_index, bric_items=bric_items,
                n=self.topics_per_generator,
            ))

        for a, b in itertools.combinations(result.topics, 2):
            sa, sb = _topic_score(a, self.weights).total, _topic_score(b, self.weights).total
            a_won = sa >= sb
            _elo_update(a, b, a_won)
            result.battle_log.append(BattleLog(
                pair=(a.title, b.title),
                winner=a.title if a_won else b.title,
                scores=(round(sa, 3), round(sb, 3)),
            ))

        for topic in result.topics:
            result.critiques[topic.title] = [c.critique(topic) for c in self.critics]

        survivors = sorted(result.topics, key=lambda t: _topic_score(t, self.weights).total, reverse=True)
        for i, parent in enumerate(survivors[:self.evolve_top_k]):
            for j, partner in enumerate(survivors[:self.evolve_top_k]):
                if i >= j:
                    continue
                child = self._mutate(parent, partner)
                result.topics.append(child)

        for topic in result.topics:
            result.final_scores[topic.title] = _topic_score(topic, self.weights)
        result.ranking = sorted(
            result.topics,
            key=lambda t: result.final_scores[t.title].total,
            reverse=True,
        )
        return result

    @staticmethod
    def _mutate(parent: GrantTopic, partner: GrantTopic) -> GrantTopic:
        merged_kw = list(dict.fromkeys(parent.keywords + partner.keywords))[:6]
        title = f"[evolved] {parent.title} + {partner.proposed_by}"
        factors = {
            k: round((parent.factor_inputs.get(k, 0.5) + partner.factor_inputs.get(k, 0.5)) / 2, 3)
            for k in set(parent.factor_inputs) | set(partner.factor_inputs)
        }
        factors["novelty"] = min(1.0, factors.get("novelty", 0.5) + 0.05)
        return GrantTopic(
            title=title, proposed_by="evolution",
            keywords=merged_kw,
            suggested_dataset=parent.suggested_dataset or partner.suggested_dataset,
            suggested_method=parent.suggested_method or partner.suggested_method,
            rationale=f"Crossover: {parent.proposed_by} × {partner.proposed_by}",
            factor_inputs=factors,
        )
