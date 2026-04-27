# References

Primary sources behind the patterns in this package.

## Anthropic primary sources

| Year | Title | URL / DOI | Used here for |
|---|---|---|---|
| 2024 | Building effective agents | anthropic.com/research/building-effective-agents | Definition of "agent" vs "workflow" (Tutorial 1) |
| 2025 | How we built our multi-agent research system | anthropic.com/engineering/multi-agent-research-system | Multi-agent orchestration shape (Tutorial 1, rung 4) |
| 2025 | Model Context Protocol specification | modelcontextprotocol.io | Tool protocol (Tutorial 5, MCP integration) |

## Academic primary sources

| Year | Authors | Title | Venue | Used here for |
|---|---|---|---|---|
| 2022 | Yao S. et al. | ReAct: Synergizing Reasoning and Acting in Language Models | ICLR 2023 (preprint 2022) | Core single-agent loop shape |
| 2023 | Shinn N. et al. | Reflexion: Language Agents with Verbal Reinforcement Learning | NeurIPS 2023 | Reviewer loop pattern |
| 2023 | Wang G. et al. | Voyager: An Open-Ended Embodied Agent | NeurIPS 2023 (preprint) | Skill library / Memory durability |
| 2023 | Madaan A. et al. | Self-Refine: Iterative Refinement with Self-Feedback | NeurIPS 2023 | Fix-Amplify cycle precursor |
| 2023 | Park J.S. et al. | Generative Agents: Interactive Simulacra of Human Behavior | UIST 2023 | Memory + retrieval architecture |

## Tooling references

| Tool | What we use it for |
|---|---|
| `aiohttp` | `executors/async_fetch.py` — multi-URL concurrent retrieval |
| `asyncio` | `core/orchestrator.py` — tier-level parallelism |
| `concurrent.futures.ProcessPoolExecutor` | `executors/parallel_pool.py` — CPU-bound parallel |
| LangGraph (referenced, not depended on) | Comparable orchestrator with LLM-decided edges |

## Where this package sits

This is *not* an LLM framework — it is a thin orchestration layer that
sits on top of one. The intent is the same as ReAct + Reflexion +
Self-Refine combined with a multi-phase project manager, but expressed as
plain Python primitives that a small team can read in ~500 LOC.
