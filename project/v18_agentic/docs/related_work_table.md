# Related work — fit table

How each prior work maps onto the five patterns in this package.

| Prior work | Fix-Amplify | Tiered DAG | Dual Dump | Reviewer Loop | Failure Catalog |
|---|---|---|---|---|---|
| ReAct (Yao 2022)              | partial | — | — | — | — |
| Reflexion (Shinn 2023)        | partial | — | — | **direct fit** | partial |
| Self-Refine (Madaan 2023)     | **direct fit** | — | — | partial | — |
| Voyager (Wang 2023)           | partial | — | — | — | **direct fit** |
| Generative Agents (Park 2023) | — | — | partial | — | — |
| LangGraph (Harrison 2024)     | partial | **direct fit** | — | — | — |
| Multi-agent research (Anthropic 2025) | — | **direct fit** | partial | partial | — |

## What is novel here

The package is not novel as a research contribution — every individual
pattern has at least one prior reference. What is novel is the
*combination*: a small set of primitives that, together, support a
multi-month research project where the human owns Goal + Verifier and the
agents own Planner + Executor + Memory.

The thyroid replay (`examples/thyroid_replay.py`) is the empirical
existence proof that the combination works on a real published-paper
target.
