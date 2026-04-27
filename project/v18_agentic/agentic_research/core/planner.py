"""Planner — decomposes a Goal into tiered Tasks.

The base implementation is rule-based and intentionally dumb: subclasses
swap in LLM-driven decomposition. The point is the seam, not the algorithm.
"""
from __future__ import annotations

from typing import Any

from .agent import Goal
from .orchestrator import Task


class Planner:
    def plan(self, goal: Goal) -> list[Task]:
        """Override in subclasses. Default: a single tier-1 placeholder task."""

        async def noop():
            return {"goal": goal.description}

        return [Task(name="placeholder", func=noop, tier=1)]
