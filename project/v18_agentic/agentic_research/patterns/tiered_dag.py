"""Pattern 2 -- Tiered Parallel DAG.

Wraps `core.orchestrator.TieredDAG` with conveniences for the common case:
declare tiers as lists, run them, get a summary back. This is the v17p35
Phase B shape: 13 tasks, 4 tiers, ~2.5 hours wall time vs ~9 hours serial.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Awaitable, Callable

from ..core.orchestrator import Task, TaskResult, TieredDAG


@dataclass
class Tier:
    name: str
    tasks: list[Task]


class TieredParallelOrchestrator:
    """Sugar around TieredDAG for the common 'declare tiers, run' workflow."""

    def __init__(self, tiers: list[Tier]):
        self.tiers = tiers
        flat: list[Task] = []
        for idx, tier in enumerate(tiers, start=1):
            for t in tier.tasks:
                t.tier = idx
                flat.append(t)
        self._dag = TieredDAG(flat)

    async def run(self) -> dict[str, Any]:
        results = await self._dag.run()
        summary = self._dag.summary(results)
        summary["per_tier_names"] = {
            i + 1: tier.name for i, tier in enumerate(self.tiers)
        }
        summary["results"] = [
            {
                "name": r.name,
                "tier": r.tier,
                "status": r.status,
                "elapsed_seconds": r.elapsed_seconds,
                "error": r.error,
            }
            for r in results
        ]
        return summary

    @staticmethod
    def task(name: str, func: Callable[[], Awaitable[Any]], timeout_s: int = 1800) -> Task:
        return Task(name=name, func=func, timeout_s=timeout_s)
