"""TieredDAG — phase-tiered parallel task orchestrator.

This is the v17p35 Phase B pattern: 13 tasks split into 4 tiers, with each
tier awaiting the previous before fanning out. Inside a tier, tasks run
concurrently via asyncio.gather.
"""
from __future__ import annotations

import asyncio
import time
from dataclasses import dataclass, field
from typing import Any, Awaitable, Callable


@dataclass
class Task:
    name: str
    func: Callable[[], Awaitable[Any]]
    tier: int = 1
    timeout_s: int = 1800
    notes: str = ""


@dataclass
class TaskResult:
    name: str
    tier: int
    status: str  # "ok" | "error" | "timeout"
    result: Any = None
    error: str | None = None
    elapsed_seconds: float = 0.0


class TieredDAG:
    """Run tasks in tiers; within a tier they run in parallel."""

    def __init__(self, tasks: list[Task]):
        self.tasks = tasks

    async def _run_one(self, task: Task) -> TaskResult:
        t0 = time.time()
        try:
            res = await asyncio.wait_for(task.func(), timeout=task.timeout_s)
            return TaskResult(name=task.name, tier=task.tier, status="ok", result=res, elapsed_seconds=time.time() - t0)
        except asyncio.TimeoutError:
            return TaskResult(name=task.name, tier=task.tier, status="timeout", elapsed_seconds=time.time() - t0)
        except Exception as e:  # noqa: BLE001
            return TaskResult(name=task.name, tier=task.tier, status="error", error=f"{type(e).__name__}: {e}", elapsed_seconds=time.time() - t0)

    async def run(self) -> list[TaskResult]:
        results: list[TaskResult] = []
        tiers = sorted({t.tier for t in self.tasks})
        for tier in tiers:
            tier_tasks = [t for t in self.tasks if t.tier == tier]
            tier_results = await asyncio.gather(*[self._run_one(t) for t in tier_tasks])
            results.extend(tier_results)
        return results

    def summary(self, results: list[TaskResult]) -> dict[str, Any]:
        return {
            "n_tasks": len(results),
            "ok": sum(1 for r in results if r.status == "ok"),
            "error": sum(1 for r in results if r.status == "error"),
            "timeout": sum(1 for r in results if r.status == "timeout"),
            "by_tier": {
                tier: {
                    "ok": sum(1 for r in results if r.tier == tier and r.status == "ok"),
                    "error": sum(1 for r in results if r.tier == tier and r.status != "ok"),
                }
                for tier in sorted({r.tier for r in results})
            },
            "total_seconds": sum(r.elapsed_seconds for r in results),
            "wall_seconds": max((r.elapsed_seconds for r in results), default=0.0),
        }
