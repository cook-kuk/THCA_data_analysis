"""Smoke tests — package imports + each pattern runs end-to-end on mocks."""
from __future__ import annotations

import asyncio
import tempfile
from pathlib import Path

import pytest

from agentic_research import (
    DualDumpRecorder,
    Fact,
    FailureCatalog,
    FixAmplifyAgent,
    Goal,
    ProjectMemory,
    ReviewerLoop,
    Task,
    TieredParallelOrchestrator,
    Verifier,
)
from agentic_research.patterns.failure_catalog import FailureRecord
from agentic_research.patterns.reviewer_loop import Critique, DefenseEntry
from agentic_research.patterns.tiered_dag import Tier


def test_imports_clean():
    assert Goal and FixAmplifyAgent and TieredParallelOrchestrator
    assert ReviewerLoop and FailureCatalog and DualDumpRecorder
    assert ProjectMemory and Verifier


def test_fix_amplify_agent_runs():
    async def fix(g, last):
        return {"step": "fix"}

    async def amp(g, last):
        return {"step": "amplify"}

    async def syn(g, last):
        return {"step": "synthesis"}

    agent = FixAmplifyAgent(
        goal=Goal(description="t", max_iterations=2),
        fix_handler=fix, amplify_handler=amp, synthesis_handler=syn,
        verifier=Verifier(),
    )
    res = asyncio.run(agent.run())
    assert res.iterations >= 1


def test_tiered_dag_runs_in_parallel():
    counter = {"calls": 0}

    async def f():
        counter["calls"] += 1
        return {"ok": True}

    orch = TieredParallelOrchestrator([
        Tier(name="t1", tasks=[Task("a", f), Task("b", f)]),
        Tier(name="t2", tasks=[Task("c", f)]),
    ])
    summary = asyncio.run(orch.run())
    assert summary["ok"] == 3
    assert counter["calls"] == 3


def test_project_memory_records_and_recalls():
    with tempfile.TemporaryDirectory() as tmp:
        mem = ProjectMemory(tmp)
        mem.record("iter_0", {"hello": "world"})
        mem.fact(Fact(topic="auc", claim="0.954", confidence="HIGH"))
        mem.fact(Fact(topic="auc", claim="0.81", confidence="LOW"))
        # _facts_cache should be invalidated when adding facts
        high_only = mem.recall("auc", min_confidence="HIGH")
        assert len(high_only) == 1
        all_facts = mem.recall("auc", min_confidence="LOW")
        assert len(all_facts) == 2


def test_dual_dump_writes_both_files():
    with tempfile.TemporaryDirectory() as tmp:
        rec = DualDumpRecorder(tmp)
        paths = rec.dump("phase_X", {"key": "value"}, "# Summary\n")
        assert paths.raw.exists()
        assert paths.compressed.exists()
        assert paths.compressed.read_text().startswith("# Summary")


def test_failure_catalog_round_trip():
    with tempfile.TemporaryDirectory() as tmp:
        cat = FailureCatalog(Path(tmp) / "fail.jsonl")
        cat.log(FailureRecord(
            failure_id="X1",
            title="bad assumption",
            what_was_tried="x",
            why_it_failed="y",
            severity="HIGH",
        ))
        cat.log(FailureRecord(
            failure_id="X2",
            title="another",
            what_was_tried="x",
            why_it_failed="y",
            severity="LOW",
        ))
        assert len(cat.list()) == 2
        assert len(cat.list(severity="HIGH")) == 1


def test_reviewer_loop_returns_attacks_and_defenses():
    async def rev(arts):
        return Critique(reviewer_id="r1", attack_vectors=[
            {"id": "A", "title": "x", "severity": 5},
        ])

    async def defb(crit, arts):
        return [DefenseEntry(attack_id=a["id"], defense_text="d") for a in crit.top_attacks(3)]

    loop = ReviewerLoop(reviewers=[rev], defense_builder=defb)
    crits, defs = asyncio.run(loop.review({}))
    assert len(crits) == 1
    assert len(defs) == 1
