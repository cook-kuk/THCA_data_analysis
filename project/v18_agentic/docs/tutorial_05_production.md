# Tutorial 5 — Production concerns

Once an agent leaves the demo, four things matter: cost, failure modes,
observability, and integration with the rest of the world.

## Cost

Most of the cost in agentic systems is context tokens, not compute. The
biggest savings come from:

1. **Compress before re-feeding.** The Dual Dump pattern produces a
   compressed markdown summary precisely so the next iteration reads the
   summary, not the raw artifact.
2. **Recall, don't re-derive.** `ProjectMemory.recall(topic)` is cheaper
   than re-running the analysis. Use it.
3. **Cap iterations.** `Goal.max_iterations` is not optional. Pick a
   number, tune by observation, never disable.

For Anthropic-API-backed implementations:

- Use prompt caching for the system prompt + tool definitions.
- Send tool results as `cache_control: ephemeral` blocks where applicable.
- A long-running agent that doesn't use prompt caching is leaving
  60-90% of its cost on the table.

## Failure modes

### Mode 1 — Loop confidence

The agent thinks it's done (verifier returned `passed=True`) but it's
actually wrong. Mitigation: **multiple verifiers**, including at least one
that didn't see the iteration. The Reviewer Loop is the canonical fix.

### Mode 2 — Context bloat

Memory grows monotonically; eventually the next iteration doesn't fit in
the model's window. Mitigation: **bounded memory** with explicit eviction.
`ProjectMemory.recall(topic, min_confidence="HIGH")` is bounded; reading
the entire `audit/` directory is not.

### Mode 3 — Untracked side effects

The agent wrote a file, made an HTTP call, edited a config, and now you
can't replay. Mitigation: **dual dump every external action**. If it
wrote bytes anywhere, log the path + hash + timestamp.

## Observability

A production agent needs three log streams:

1. **Audit** (per-iteration artifacts) — `ProjectMemory.audit_dir`.
2. **Trace** (per-tool-call trace) — wrap your `Tool` implementations to
   emit a structured event per call.
3. **Verdict** (per-iteration verifier output) — already in
   `AgentResult.audit_trail`.

If you can't answer "what did the agent do at iteration 3?" in 10 seconds,
fix the logging.

## MCP integration

For Claude-based deployments, exposing your tools via the Model Context
Protocol (MCP) is the right primitive. The package's `Tool` protocol is
intentionally narrow so it can be implemented by either a local Python
function or an MCP-backed call.

```python
class MCPTool:
    name: str = "mcp_query"
    async def __call__(self, query: str) -> dict:
        return await mcp_client.call("query", {"q": query})
```

## Kill switches

Anything that runs more than a minute in production needs:

1. A wall-time budget (use `asyncio.wait_for` at the top of `run()`).
2. A cost budget (track tokens; abort if exceeded).
3. A signal handler (SIGTERM → flush memory → exit).

The `BaseAgent.run()` loop is intentionally trivial so you can wrap it
with whatever budget enforcement your runtime needs.

## Anti-patterns

- **No verifier.** A loop without a verifier is a pipeline.
- **Verifier that never says no.** A verifier that always returns
  `passed=True` is no verifier at all. Test it on bad artifacts.
- **Hidden global state.** Module-level globals make replay impossible.
  Pass state through `Goal.domain_prior` and `ProjectMemory`.
- **Ad-hoc retry loops.** If your code has `for _ in range(3): try: ...
  except: continue`, you wanted the failure catalog instead.
