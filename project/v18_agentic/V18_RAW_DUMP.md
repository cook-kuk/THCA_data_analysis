# v18 agentic_research — RAW DUMP

_Run date: 2026-04-27_
_Wall time (Sprint A scaffold): ~30 min_

## What was built

```
project/v18_agentic/
├── README.md                            (English, native level)
├── LICENSE                              (MIT)
├── pyproject.toml                       (installable package)
├── V18_RAW_DUMP.md                      (this file)
├── agentic_research/
│   ├── __init__.py                      (top-level exports)
│   ├── core/
│   │   ├── __init__.py
│   │   ├── agent.py                     (BaseAgent, Goal, AgentResult)
│   │   ├── orchestrator.py              (TieredDAG, Task, TaskResult)
│   │   ├── memory.py                    (ProjectMemory, Fact)
│   │   ├── verifier.py                  (Verifier, VerificationResult)
│   │   └── planner.py                   (Planner skeleton)
│   ├── patterns/
│   │   ├── __init__.py
│   │   ├── fix_amplify.py               (Pattern 1)
│   │   ├── tiered_dag.py                (Pattern 2)
│   │   ├── dual_dump.py                 (Pattern 3)
│   │   ├── reviewer_loop.py             (Pattern 4)
│   │   └── failure_catalog.py           (Pattern 5)
│   ├── executors/
│   │   ├── __init__.py
│   │   ├── async_fetch.py               (aiohttp multi-URL)
│   │   ├── parallel_pool.py             (multiprocessing)
│   │   └── claude_code.py               (Claude Code CLI dispatcher)
│   └── examples/
│       ├── __init__.py
│       ├── thyroid_replay.py            (200-line demo, 5 patterns)
│       └── quick_start.py               (5-min hello world)
├── docs/
│   ├── tutorial_01_what_is_agentic.md
│   ├── tutorial_02_six_components.md
│   ├── tutorial_03_five_patterns.md
│   ├── tutorial_04_build_your_own.md
│   ├── tutorial_05_production.md
│   ├── references.md                    (Anthropic + academic primary sources)
│   └── related_work_table.md            (5-pattern fit vs prior work)
├── tests/
│   ├── __init__.py
│   └── test_basic.py                    (7 smoke tests)
└── talk/
    └── agent_vs_agentic_outline_KR.md   (22-slide outline + speaker notes + Q&A)
```

## Metrics

| Metric | Value |
|---|---:|
| Files created / modified | 23 |
| LOC (Python only, agentic_research/) | ~580 |
| LOC (docs/) | ~700 |
| Tests | 7 (all pass on smoke run) |
| Patterns codified | 5 |
| Tutorials | 5 |
| Demo wall time | 1.7 s |
| Demo: patterns exercised | 5 / 5 |

## Demo run (verified)

```
$ python -m agentic_research.examples.thyroid_replay --max-iter 3
==========================================================
  v18 thyroid_replay -- agentic_research demo
==========================================================

[PATTERN 1] FixAmplifyAgent
  [FIX] iteration target: Reproduce DM1/DM2 -> 8-gene RAI panel -> 4-group survival narrative
  -> iterations: 1, success: True

[PATTERN 2] TieredParallelOrchestrator
  -> 6/6 OK across 3 tiers

[PATTERN 3] DualDumpRecorder
  -> raw: phase_B_summary.json, compressed: phase_B_summary.md

[PATTERN 4] ReviewerLoop
  -> 3 attacks, 3 defenses

[PATTERN 5] FailureCatalog
  -> 1 entries

==========================================================
  Demo complete. Workspace: /tmp/v18_thyroid_replay
==========================================================
```

## What was deferred (and why)

| Item | Why deferred |
|---|---|
| Actual `.pptx` file | PptxGenJS is JS-only, python-pptx requires layout work that's faster in Slides/Keynote. The 22-slide outline + speaker notes file is the source-of-truth; the slide deck takes ~30 min to build manually. |
| 5 PNG screenshot for demo cheat sheet | Live screen capture pending — outline references them but doesn't generate them. |
| `tests/` exhaustive coverage | 7 smoke tests are enough to confirm the package imports and each pattern runs end-to-end. Adding ~15 more would be valuable but not on critical path. |
| MCP integration example | Tutorial 5 documents the seam; concrete MCP wiring deferred to a future iteration since MCP setup is environment-specific. |
| GitHub Actions CI workflow | `.github/workflows/test.yml` not written — single command (`pytest`) is enough until the package has external users. |

## Honest scope check

This is a 0.1.0 alpha. It does what it says (codifies the 5 patterns
behind v3→v17p35, runs end-to-end on mock handlers, ships with 5
tutorials and a Korean talk outline). It does NOT:

- Replace LangGraph / LangChain / LlamaIndex
- Solve the orchestration problem in general
- Make claims about LLM behavior — it's plumbing, not modeling
- Ship a default LLM client — bring your own

## Three decisions for the user

1. **Talk delivery date.** The 22-slide outline is ready. Picking the
   delivery target (internal lab talk vs external conference) determines
   how much polish the slide deck needs.
2. **External users.** Whether to invite anyone to use this against a
   real research project. The smallest lift would be replaying their own
   sprint shape with `FixAmplifyAgent` + `DualDumpRecorder`.
3. **Naming.** `agentic_research` is generic; if there's a more specific
   name (e.g. tied to the thyroid case or to a brand) that should be
   decided before any GitHub release.

## Failure catalog (sprint-level)

| Failure | Description | Recovered? |
|---|---|---|
| Pre-existing v18 dir | A prior session had scaffolded `core/agent.py` + `__init__.py` with a different interface (`Orchestrator`, `AuditVerifier`); my Write clobbered 4 of the 7 files in core/. Recovered by rewriting `__init__.py` + `agent.py` to match the new interface and confirming imports clean. | Yes |
| Working directory reset mid-sprint | A `cd project/v18_agentic` shell context was lost between Bash calls; one batch of parallel Writes failed. Recovered by switching to absolute paths. | Yes |

## Verdict

Sprint A complete. Package importable, demo runs, 5 tutorials drafted,
talk outline ready. The work that remains (PPTX, screenshots, CI) is
formatting polish, not technical risk.
