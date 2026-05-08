# grant_arena — competitive multi-agent layer on research_compass

> Distilled from a long brainstorm with another LLM. Layered on top of
> `research_compass` v0 to turn topic *recommendation* into topic
> *discovery → battle → critique → evolution → ranked proposal seeds*.

## Layout

```
grant_arena/
├── scoring.py                 10-factor composite (mega-grant tuned)
├── arena.py                   Generate → Elo battle → critique → evolve → judge
├── frontier_scout.py          Watchlist top researchers, extract method/keyword signals
├── agents/
│   ├── generators.py          8 generators (Clinical, Frontier, Biology, BioData,
│   │                          Korean, Patent, Funding, Crazy)
│   └── critics.py             5 critics (Feasibility, Novelty, Patent, Reviewer-2, Budget)
└── mcp_servers/               7 MCP-shaped tool stubs (offline; live in v0.3)
    ├── rfp_reader.py
    ├── biodata_reuse.py       BioData Reuse Score (TCGA, GEO, DepMap, ...)
    ├── jcr_paper.py           Q1 biomed allowlist
    ├── ai_topvenue.py         NeurIPS/ICML/ICLR/AAAI/CVPR/ACL/...
    ├── bric_trend.py          BRIC 한빛사 fetcher slot
    └── patent_priorart.py     Prior-art collision risk screening (NOT legal FTO)
```

## How it sits on research_compass

| research_compass v0           | grant_arena v0.2 extension                         |
|-------------------------------|----------------------------------------------------|
| `cv_profile.ResearchProfile`  | feeds generators directly                          |
| `journal_feed`                | augmented by `bric_trend` + `frontier_scout`       |
| `method_index`                | augmented by `frontier_scout.merge_into_method_index_seed` |
| `data_registry.DatasetEntry`  | scored by `biodata_reuse` MCP                      |
| `topic_ranker` (4-factor)     | superseded by `scoring` (10-factor) + `arena`      |

## 10-factor composite

```
score = + 0.18 novelty
        + 0.15 feasibility
        + 0.15 pi_fit
        + 0.12 funder_fit
        + 0.10 publishability
        + 0.10 patentability
        + 0.10 clinical_impact
        + 0.05 scalability
        + 0.05 pilot_speed
        - 0.15 patent_collision_risk
        - 0.10 overclaim_risk
```

Tiers: HIGH ≥0.55, MEDIUM ≥0.30, LOW <0.30.

## Frontier-scout idea

The watchlist captures top researchers (Hernán, Pearl, Topol, Rajpurkar,
Aerts, Mahmood, Koller, Barzilay, ...). For each, scout their last 1-3
years of papers, run lexicon-based method/keyword extraction, and feed the
counts back into:

- `MethodIndex` — bumps for methods that surge across the watchlist
- `FrontierTechAgent` priors — biases topic generation toward what's *currently hot*
- `bric_trend` cross-check — flag Korean PIs working on the same hot methods

This is what makes the system *learn what's working globally* instead of
just remixing static seed entries.

## Status

- Offline scaffold: complete, tested.
- Live MCP wiring (PubMed esearch, arXiv API, bioRxiv, Semantic Scholar Graph,
  KIPRIS, WIPO, Google Patents, JCR/Clarivate): deferred to v0.3.
- LLM-backed judge / generator / mutator: deferred to v0.3.

## Marathon-mode note

This is scaffolding/infra layered on top of the v18 framework. It does not
generate manuscript prose, does not touch voice-protected sections, and
does not require new prospective data collection. Live wiring + LLM
generators are deferred to post-2026-06-13.
