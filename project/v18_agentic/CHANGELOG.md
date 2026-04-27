# Changelog

All notable changes to this project will be documented in this file.

## [0.1.0] — 2026-04-27

### Added
- `BaseAgent`, `Goal`, `AgentResult` core abstractions.
- `TieredDAG` and `TieredParallelOrchestrator` for tier-level parallel execution.
- `ProjectMemory` and `Fact` for cross-iteration memory + structured fact recall.
- `Verifier` and `VerificationResult` with pluggable check protocol.
- `Planner` skeleton.
- Five pattern modules: `FixAmplifyAgent`, `TieredParallelOrchestrator`,
  `DualDumpRecorder`, `ReviewerLoop`, `FailureCatalog`.
- Three executor implementations: `AsyncFetchExecutor`, `ParallelPoolExecutor`,
  `ClaudeCodeExecutor`.
- 200-line `examples/thyroid_replay.py` demo exercising all five patterns.
- 5-minute `examples/quick_start.py` hello-world.
- 7 smoke tests in `tests/test_basic.py`.
- 5 tutorials covering: what is agentic, six components, five patterns,
  build-your-own, and production concerns.
- 22-slide Korean talk outline + speaker notes + demo cheat sheet.
- References + related work table mapping patterns to prior literature.
- GitHub Actions test workflow.

### Notes
- Alpha release. Interfaces may change before 0.2.0.
- No PyPI release yet — install via `pip install -e .`.
