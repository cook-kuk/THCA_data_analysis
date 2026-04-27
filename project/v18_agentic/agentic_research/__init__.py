"""agentic_research -- patterns for multi-phase research orchestration.

Distilled from a 6-month thyroid-cancer transcriptomic project (v3 -> v17p35)
that ran ~80 phase-level sprints across discovery, validation, and paper
preparation. The framework codifies five recurring patterns and the six
agent components those patterns assume.

Quick start
-----------
>>> import asyncio
>>> from agentic_research import Goal, FixAmplifyAgent
>>> # see examples/thyroid_replay.py for a working demo
"""
from .core.agent import BaseAgent, AgentResult, Goal
from .core.orchestrator import TieredDAG, Task, TaskResult
from .core.memory import ProjectMemory, Fact
from .core.verifier import Verifier, VerificationResult
from .core.planner import Planner

from .patterns.fix_amplify import FixAmplifyAgent
from .patterns.tiered_dag import TieredParallelOrchestrator
from .patterns.dual_dump import DualDumpRecorder
from .patterns.reviewer_loop import ReviewerLoop
from .patterns.failure_catalog import FailureCatalog

__version__ = "0.1.0"
__all__ = [
    "BaseAgent",
    "AgentResult",
    "Goal",
    "TieredDAG",
    "Task",
    "TaskResult",
    "ProjectMemory",
    "Fact",
    "Verifier",
    "VerificationResult",
    "Planner",
    "FixAmplifyAgent",
    "TieredParallelOrchestrator",
    "DualDumpRecorder",
    "ReviewerLoop",
    "FailureCatalog",
]
