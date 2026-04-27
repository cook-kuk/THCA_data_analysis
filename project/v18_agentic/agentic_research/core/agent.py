"""BaseAgent — the six-component spine.

Every agent has: Goal, Planner, Executor, Memory, Tools, Verifier. The base
class wires them together; subclasses (FixAmplifyAgent, ReviewerLoop, etc.)
override the loop shape but reuse the components.

Distilled from the six manual roles played during the v3->v17p35 thyroid
sprint:
  1. Goal owner   -- what does "done" look like?
  2. Planner      -- goal -> task list with a DAG
  3. Executor     -- tasks -> side effects (files, API calls, models)
  4. Memory       -- persistent project state across iterations
  5. Tools        -- bash, fetch, LLM, MCP -- the things the executor calls
  6. Verifier     -- was the iteration's output good enough? if not, re-plan
"""
from __future__ import annotations

import time
from dataclasses import dataclass, field
from typing import Any, Awaitable, Callable, Protocol


@dataclass
class Goal:
    """What the agent is trying to achieve, in plain text + structured prior."""

    description: str
    domain_prior: dict[str, Any] = field(default_factory=dict)
    success_criteria: list[str] = field(default_factory=list)
    max_iterations: int = 5


@dataclass
class AgentResult:
    success: bool
    iterations: int
    artifacts: dict[str, Any] = field(default_factory=dict)
    audit_trail: list[dict[str, Any]] = field(default_factory=list)
    elapsed_seconds: float = 0.0
    notes: list[str] = field(default_factory=list)


class Tool(Protocol):
    name: str

    async def __call__(self, *args, **kwargs) -> Any: ...  # noqa: D401


class BaseAgent:
    """Skeleton agent. Override `step` for the work loop."""

    def __init__(
        self,
        goal: Goal,
        planner=None,
        executor: Callable[[Any], Awaitable[Any]] | None = None,
        memory=None,
        tools: list[Tool] | None = None,
        verifier=None,
    ):
        self.goal = goal
        self.planner = planner
        self.executor = executor
        self.memory = memory
        self.tools = tools or []
        self.verifier = verifier

    async def step(self, iteration: int) -> dict[str, Any]:
        """Override in subclasses. Returns one iteration's artifacts + diagnostics."""
        raise NotImplementedError

    async def run(self) -> AgentResult:
        t0 = time.time()
        result = AgentResult(success=False, iterations=0)
        for i in range(self.goal.max_iterations):
            artifacts = await self.step(i)
            result.iterations = i + 1
            result.audit_trail.append({"iteration": i, "artifacts_keys": list(artifacts.keys())})
            if self.memory is not None:
                self.memory.record(f"iter_{i}", artifacts)
            verdict = None
            if self.verifier is not None:
                verdict = await self.verifier.verify(artifacts, self.goal)
                result.audit_trail[-1]["verdict"] = verdict.summary()
            if verdict is not None and verdict.passed:
                result.success = True
                result.artifacts.update(artifacts)
                break
            result.artifacts.update(artifacts)
        result.elapsed_seconds = time.time() - t0
        return result
