"""Pattern 1 -- Fix-Amplify-Synthesis cycle.

The shape that drove v3 -> v17p35:
  Phase N      ("Fix")       -- targeted intervention on a specific defect
  Phase N+1    ("Amplify")    -- extend the fix's surface area, run downstream
  Phase N+2    ("Synthesis")  -- audit, write up, declare done or re-enter

The agent body alternates between these phases; the verifier decides which
to dispatch next.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Awaitable, Callable

from ..core.agent import BaseAgent, AgentResult, Goal


class Phase(str, Enum):
    FIX = "fix"
    AMPLIFY = "amplify"
    SYNTHESIS = "synthesis"


@dataclass
class CycleStep:
    phase: Phase
    handler: Callable[[Goal, dict[str, Any]], Awaitable[dict[str, Any]]]


class FixAmplifyAgent(BaseAgent):
    """Agent that cycles fix -> amplify -> synthesis until the verifier passes.

    Concrete usage: pass three handler coroutines, one per phase. Each handler
    receives (goal, last_artifacts) and returns new artifacts. The verifier
    decides whether to break or re-enter the cycle.
    """

    def __init__(
        self,
        goal: Goal,
        fix_handler: Callable[[Goal, dict[str, Any]], Awaitable[dict[str, Any]]],
        amplify_handler: Callable[[Goal, dict[str, Any]], Awaitable[dict[str, Any]]],
        synthesis_handler: Callable[[Goal, dict[str, Any]], Awaitable[dict[str, Any]]],
        memory=None,
        verifier=None,
    ):
        super().__init__(goal=goal, memory=memory, verifier=verifier)
        self._handlers = {
            Phase.FIX: fix_handler,
            Phase.AMPLIFY: amplify_handler,
            Phase.SYNTHESIS: synthesis_handler,
        }
        self._phase_seq = [Phase.FIX, Phase.AMPLIFY, Phase.SYNTHESIS]
        self._last_artifacts: dict[str, Any] = {}

    async def step(self, iteration: int) -> dict[str, Any]:
        phase = self._phase_seq[iteration % len(self._phase_seq)]
        handler = self._handlers[phase]
        new = await handler(self.goal, self._last_artifacts)
        new["__phase"] = phase.value
        self._last_artifacts = {**self._last_artifacts, **new}
        return new
