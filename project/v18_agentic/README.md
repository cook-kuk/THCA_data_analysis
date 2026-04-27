# agentic-research

Reusable patterns for multi-phase research orchestration, distilled from a
6-month thyroid-cancer transcriptomic project that ran roughly 80 phase-level
sprints (`v3 -> v17p35`) across discovery, validation, and paper preparation.

This is not a framework that schedules LLM calls for you. It is a small set
of primitives — `BaseAgent`, `TieredDAG`, `ProjectMemory`, `FixAmplifyAgent`,
`ReviewerLoop`, `DualDumpRecorder`, `FailureCatalog` — that codify the five
patterns that kept that project moving when manual orchestration would have
collapsed.

## Why this exists

When the work is "produce a defensible thyroid-cancer paper," you don't need
a long-running autonomous agent. You need:

1. A way to alternate between targeted fixes, broad amplification, and
   honest synthesis without losing the audit trail (**fix-amplify cycle**).
2. A way to fan out a sprint into a few tiers of parallel tasks and survive
   when one of them fails (**tiered parallel DAG**).
3. A way to record both the raw artifact (replayable) and a compressed
   markdown summary (skim-able) for every phase (**dual dump**).
4. A way to run an external reviewer over your output and ingest the
   critique as the next iteration's goal (**reviewer loop**).
5. A way to log failed approaches with enough context that future-you doesn't
   repeat them (**failure catalog**).

These five patterns map onto six classical agent components: **Goal,
Planner, Executor, Memory, Tools, Verifier**.

## Install

```bash
pip install -e .
# optional demo extras
pip install -e .[demo]
```

## Quick start

```python
import asyncio
from agentic_research import Goal, FixAmplifyAgent, Verifier

async def fix(goal, last):       return {"step": "fix",       "claim": "found a defect"}
async def amplify(goal, last):   return {"step": "amplify",   "claim": "extended to 3 cohorts"}
async def synthesis(goal, last): return {"step": "synthesis", "claim": "wrote the section"}

agent = FixAmplifyAgent(
    goal=Goal(description="Demo", max_iterations=3),
    fix_handler=fix, amplify_handler=amplify, synthesis_handler=synthesis,
    verifier=Verifier(),
)
print(asyncio.run(agent.run()))
```

## Live demo

```bash
python -m agentic_research.examples.thyroid_replay --max-iter 3
```

Walks one full agent through all five patterns using mock handlers. No
external services required.

## Layout

```
agentic_research/
  core/
    agent.py          BaseAgent, Goal, AgentResult, Tool protocol
    orchestrator.py   TieredDAG, Task, TaskResult
    memory.py         ProjectMemory + Fact store
    verifier.py       Verifier + check protocol
    planner.py        Planner skeleton
  patterns/
    fix_amplify.py    FixAmplifyAgent (Pattern 1)
    tiered_dag.py     TieredParallelOrchestrator (Pattern 2)
    dual_dump.py      DualDumpRecorder (Pattern 3)
    reviewer_loop.py  ReviewerLoop, Critique, DefenseEntry (Pattern 4)
    failure_catalog.py FailureCatalog, FailureRecord (Pattern 5)
  executors/
    async_fetch.py    aiohttp multi-URL fetcher with audit
    parallel_pool.py  multiprocessing pool wrapper
    claude_code.py    Claude Code CLI subprocess dispatcher
  examples/
    thyroid_replay.py 200-line reproduction of the v17 sprint shape
    quick_start.py    5-minute hello world
docs/
  tutorial_01_what_is_agentic.md
  tutorial_02_six_components.md
  tutorial_03_five_patterns.md
  references.md
  related_work_table.md
```

## Tutorials

1. [What is agentic? (Anthropic-aligned definition)](docs/tutorial_01_what_is_agentic.md)
2. [The six components of an agent](docs/tutorial_02_six_components.md)
3. [The five patterns](docs/tutorial_03_five_patterns.md)

## Status

Alpha. The interfaces are stable enough to build on but not stable enough to
claim semver compliance. Pinned to `0.1.x` until at least one external user
reports back.

## License

MIT.

## Acknowledgement

This package is a distillation, not an invention. The five patterns are
rediscoveries of practices that show up in the agent literature
(Reflexion, ReAct, multi-agent research systems) and in the practical
multi-phase project shown in `examples/thyroid_replay.py`. See
`docs/references.md` for primary sources.
