"""Pattern 4 -- External Reviewer Loop.

In the v3 -> v17p35 sprint we ran a recurring loop:
  1. Local agent produces an artifact (analysis, manuscript draft, figure).
  2. An external reviewer (a different LLM, or a peer) reads the compressed
     dump and returns critique + ranked attack vectors.
  3. The local agent ingests the critique, generates a defense list, and
     re-enters the work loop with the highest-leverage gaps as the new goal.

This pattern is the single highest-leverage way to escape local optima.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Awaitable, Callable

from ..core.agent import Goal


@dataclass
class Critique:
    reviewer_id: str
    attack_vectors: list[dict[str, Any]] = field(default_factory=list)
    raw_text: str = ""

    def top_attacks(self, k: int = 3) -> list[dict[str, Any]]:
        return sorted(
            self.attack_vectors, key=lambda a: a.get("severity", 0), reverse=True
        )[:k]


@dataclass
class DefenseEntry:
    attack_id: str
    defense_text: str
    artifact_paths: list[str] = field(default_factory=list)


class ReviewerLoop:
    """Run an artifact through one or more external reviewers and collect defenses."""

    def __init__(
        self,
        reviewers: list[Callable[[dict[str, Any]], Awaitable[Critique]]],
        defense_builder: Callable[[Critique, dict[str, Any]], Awaitable[list[DefenseEntry]]],
    ):
        self.reviewers = reviewers
        self.defense_builder = defense_builder

    async def review(
        self, artifacts: dict[str, Any]
    ) -> tuple[list[Critique], list[DefenseEntry]]:
        import asyncio

        critiques = await asyncio.gather(*[r(artifacts) for r in self.reviewers])
        defenses: list[DefenseEntry] = []
        for c in critiques:
            d = await self.defense_builder(c, artifacts)
            defenses.extend(d)
        return list(critiques), defenses

    @staticmethod
    def goal_from_attacks(attacks: list[dict[str, Any]], parent: Goal) -> Goal:
        """Promote the top unaddressed attack into the next iteration's goal."""
        descriptor = "; ".join(a.get("title", a.get("description", "")) for a in attacks)
        return Goal(
            description=f"Address reviewer attacks: {descriptor}",
            domain_prior=parent.domain_prior,
            success_criteria=[
                f"defense for: {a.get('title')}" for a in attacks
            ],
            max_iterations=parent.max_iterations,
        )
