# Tutorial 1 — What is "agentic"?

The word "agent" is overloaded. This tutorial pins it down for the purposes
of this package.

## Anthropic-aligned definition

From Anthropic's "Building effective agents" (2024):

- **Workflow** — code paths predetermined by the developer; the LLM picks
  among fixed branches.
- **Agent** — the LLM dynamically directs its own tool use and control
  flow.

That distinction is what we mean by *agentic*. A pipeline that runs three
LLM calls in a fixed sequence is not agentic. A loop that lets the LLM
decide whether to re-enter the loop, what tool to call next, or what
sub-goal to dispatch is.

## The five-rung spectrum

| Rung | Name              | Who decides control flow? | Example                          |
|-----:|-------------------|---------------------------|----------------------------------|
|    1 | Workflow          | Developer                 | LangChain SequentialChain        |
|    2 | Tool-augmented    | Developer + LLM picks tool| Function-calling                 |
|    3 | Single agent      | LLM decides loop          | ReAct, AutoGPT                   |
|    4 | Multi-agent       | Orchestrator agent        | Multi-agent research system      |
|    5 | Autonomous agent  | Self-deployed             | Long-running ops agents          |

This package targets rungs 3 and 4. Rung 5 needs production guardrails
(rate limits, money budgets, kill switches) outside the scope here.

## Why the distinction matters

The cost / latency / reliability tradeoffs change at each rung:

- Rung 1 is cheap, fast, predictable. Use it when you can.
- Rung 3 is expensive and non-deterministic. Use it only when the control
  flow really cannot be predetermined — e.g. when the next step depends on
  an artifact you don't have yet.
- Rung 5 will lose money if you let it.

The five patterns in this package are designed for rung 3-4 work where the
human still owns the goal and the verifier, but the LLM owns the
intermediate planning.

## Anti-patterns

- **Rung-3 framework for a rung-1 problem.** If your control flow is
  knowable, write the loop. Don't ask an LLM to decide.
- **Rung-4 multi-agent for a rung-3 problem.** Two agents that mostly do
  what one agent could do are a context-window leak.
- **Hidden non-determinism.** If the LLM picks the next step, log the
  decision. The audit trail is what makes the work defensible.
