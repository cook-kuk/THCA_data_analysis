# Tutorial 2 — The six components

Every agent in this package has the same six components. They show up in
`core/agent.py`'s `BaseAgent` constructor.

| # | Component  | What it owns                              | This package's class    |
|--:|-----------|-------------------------------------------|-------------------------|
| 1 | Goal      | "Done" definition + success criteria      | `Goal` (dataclass)      |
| 2 | Planner   | Goal → tiered task list                   | `Planner`               |
| 3 | Executor  | Task → side effects                       | `executors/*`           |
| 4 | Memory    | Cross-iteration state + audit             | `ProjectMemory`, `Fact` |
| 5 | Tools     | Things the executor calls                 | Implementer's choice    |
| 6 | Verifier  | "Was the artifact good enough?"           | `Verifier`              |

## How they wire up

```python
agent = FixAmplifyAgent(
    goal=Goal(description="...", success_criteria=["AUC > 0.95"]),
    fix_handler=...,        # implements step "fix"     using executor + tools
    amplify_handler=...,    # implements step "amplify"
    synthesis_handler=...,  # implements step "synthesis"
    memory=ProjectMemory("./workspace"),
    verifier=Verifier(checks=[my_auc_check]),
)
result = await agent.run()
```

The base class does the loop:

```
for i in range(goal.max_iterations):
    artifacts = await self.step(i)         # ← uses Planner / Executor / Tools
    self.memory.record(...)                # ← Memory writes
    verdict = await self.verifier.verify(artifacts, self.goal)
    if verdict.passed: break
```

## What each component should and shouldn't do

### Goal

**Should** carry the success criteria as plain English plus structured priors
(domain knowledge, paths, panel names). **Should not** carry mutable state
— that's Memory's job.

### Planner

**Should** decompose goals into tiered tasks. **Should not** execute the
tasks — that's the Executor.

### Executor

**Should** know how to actually run a task (subprocess, HTTP call, MP pool).
**Should not** decide what to run — that's the Planner.

### Memory

**Should** persist three things: per-iteration audit, structured facts
(claims, AUCs, p-values), and pointers to large external artifacts.
**Should not** become a database — keep it boring.

### Tools

**Should** be small, composable, individually testable. **Should not** be
where business logic lives.

### Verifier

**Should** make a binary "good enough?" decision and explain why if it
says no. **Should not** be optional — without a verifier the agent is just
a sequential pipeline.

## Mapping to the v17 thyroid sprint

| Component | What it was during v3 → v17p35 |
|---|---|
| Goal      | "Defensible npj submission" |
| Planner   | The user, picking phase boundaries by hand |
| Executor  | Multiple Claude Code sessions across terminals |
| Memory    | The repo + `MEMORY.md` file in `.claude/` |
| Tools     | bash, fetch, plot, statsmodels, lifelines, scipy |
| Verifier  | The user reading the dual dump + an external LLM |

Most of the "manual orchestration" the user did was acting as Planner and
Verifier. The point of this package is to let those two roles be
componentized so future projects don't pay the manual cost again.
