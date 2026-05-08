# cancer_discovery_os — Yu-Cook Lab integration layer

> One module that knows about every Yu-Cook asset (papers, cohorts, datasets,
> discoveries, platforms, grants) and routes free-text queries to the right
> subsystem (research_compass, grant_arena, frontier_scout, Biomni).

## What it integrates

| Yu-Cook asset                            | Where it shows up                            |
|------------------------------------------|----------------------------------------------|
| Paper 1 (DM1), Paper 2 (image→DM1)       | `manifest.default_yu_cook_manifest()`        |
| Paper 3 (ICI), Paper 4 (Korean GD HLA)   | manifest entries with status flags           |
| Paper 11 (pan-cancer prognostic)         | manifest entries                             |
| Korean K2 / BABA / CTC-EMT / MAeSTro     | manifest cohort entries                      |
| GSE286332 / proteogenomic / RAI-dediff   | manifest dataset+discovery entries           |
| v18 framework / research_compass / grant_arena / neoantigen hub | manifest platform entries |
| Causal K-Thyro 30억 grant design         | manifest grant entry                         |

| External capability                      | Adapter                                       |
|------------------------------------------|-----------------------------------------------|
| Stanford SNAP Biomni (CRISPR/scRNA/...)  | `biomni_adapter.BiomniAdapter` (stub now)     |
| Top-researcher paper monitoring          | `grant_arena.frontier_scout`                  |

## Runner dispatch

```python
from agentic_research.domains.cancer_discovery_os import CancerDiscoveryOS, OSQuery

os = CancerDiscoveryOS()
os.handle(OSQuery("what assets do we have on PTC?"))         # → asset lookup
os.handle(OSQuery("plan a CRISPR screen for thyroid dedifferentiation"))  # → BiomniAdapter (stub)
os.handle(OSQuery("Samsung grant topic for Yu-Cook combo", payload={...}))  # → GrantArena
os.handle(OSQuery("recommend a research direction", payload={"cv": "..."}))  # → ResearchCompass
```

Intent classification is rule-based for v0; v0.3 swaps in an LLM
classifier and wires real Biomni `A1.go(...)` calls.

## What's *not* in this module

- LLM API calls (compass / arena / Biomni all run offline)
- Biomni installation (8GB+ datalake, code execution sandbox)
- Real PubMed / arXiv / bioRxiv fetchers
- Any UI / dashboard

These all become v0.3 after paper ship.

## Marathon-mode contract

- 2026-05-04 to 2026-06-13: Paper 1 (DM1) bioRxiv, Paper 2 (image→DM1),
  Paper 3 (ICI Track A submit), Paper 11 (pan-cancer) take priority.
- This module is *scaffolding* — interfaces, manifest, adapters. No
  manuscript prose, no voice-protected sections, no new analyses.
- Do NOT enable `BiomniAdapter.available = True` and start running
  `agent.go(...)` until paper ship is complete.
