# research_compass — CV-conditioned topic recommendation

> "phylo.bio composite scoring × feynman.is multi-agent × v18 TieredParallelOrchestrator"

The `research_compass` domain extends `agentic_research` from a research
*execution* framework into a research *targeting* framework. Given a
researcher's CV plus live signals from the literature and open data, it
returns ranked topic hypotheses tailored to that specific researcher's
competence, available datasets, and method fit.

## Five-stage pipeline

```
Stage 1  cv_profile      CV markdown      ->  ResearchProfile (skills, prior topics, methods, datasets)
Stage 2  journal_feed    last N days      ->  JournalItem[]   (BRIC 한빛사, bioRxiv, top-tier RSS)   ─┐
Stage 3  method_index    method library   ->  MethodEntry[]   (keyword-indexed; embedding-ready)     │ Tier 2 fan-out
Stage 4  data_registry   open data        ->  DatasetEntry[]  (GEO, SRA, cBioPortal, UKB)            ─┘
Stage 5  topic_ranker    join + score     ->  TopicHypothesis[] (HIGH/MEDIUM/LOW tiers, phylo style)
```

Stages 2–4 fan out via `TieredParallelOrchestrator`. Stage 5 joins.

## Composite score

Mirrors phylo.bio's 35/25/25/15 split, retargeted to topic recommendation:

| Component   | Weight | What it measures |
|-------------|--------|------------------|
| competence  | 0.35   | Jaccard overlap of CV terms with paper keywords |
| novelty     | 0.25   | Recency signal (≥2024-01-01 → 1.0, else 0.5) |
| feasibility | 0.25   | Best-matching dataset's overlap × access (open=1, controlled=0.5) |
| method_fit  | 0.15   | Best-matching method's keyword overlap with CV methods + paper |

Tiers: `score ≥ 0.55 → HIGH`, `≥ 0.30 → MEDIUM`, `< 0.30 → LOW`.

## Quick start

```bash
cd project/v18_agentic
python3 -m examples.research_compass_demo     # offline, stub data
pytest tests/test_research_compass.py -q      # 4 unit tests
```

```python
from agentic_research.domains.research_compass import ResearchCompass

compass = ResearchCompass(live=False, stub_journal_items=..., stub_datasets=...)
result = await compass.run(cv_text)
for h in result.hypotheses:
    print(h.tier, h.score, h.title, "->", h.suggested_dataset.accession)
```

## Offline-first design

Every fetcher accepts a `live=False` default and a `stub_*` injection slot.
This keeps tests/demo deterministic and fast (~0.03s) and lets you build
the pipeline against real data sources incrementally without breaking the
contract.

Live wiring TODOs (sorted by effort):

1. **bioRxiv API** (easy) — `GET /details/biorxiv/{from}/{to}` + client-side keyword filter.
2. **BRIC 한빛사 scraper** (easy) — listing page is static HTML.
3. **GEO esearch + esummary** (easy) — NCBI Entrez, no key needed for low rate.
4. **cBioPortal `/api/studies`** (easy) — public REST.
5. **Top-tier RSS fan-out** (easy) — `feedparser` per journal.
6. **SRA esearch** (medium) — accession metadata is fine; pulling FASTQ isn't this layer's job.
7. **UKB Showcase field metadata** (medium) — open metadata only; participant data stays DUA-gated.
8. **Sentence-transformer method index** (medium) — swap `MethodIndex.search()` for embedding cosine while keeping the API.
9. **CV → ResearchProfile via LLM** (medium) — current regex-based parser is a strong baseline; LLM extraction adds nuance for unstructured CVs.
10. **Reviewer loop integration** (medium) — wrap the ranker in `ReviewerLoop` so a critic agent attacks each hypothesis card before output.

## How this relates to the v18 framework

`research_compass` is the first concrete `domains/` module. It demonstrates
that the v18 spine (Goal/Planner/Executor/Memory/Tools/Verifier) and the 5
patterns (FixAmplify, TieredDAG, DualDump, ReviewerLoop, FailureCatalog)
are domain-agnostic — the same patterns that drove the v3→v17p35 thyroid
sprint can drive a topic-discovery agent, with no changes to `core/` or
`patterns/`.

## Marathon-mode note (2026-05-08)

This is scaffolding/infrastructure under marathon-mode rules (5/4–6/13);
it is not paper-blocking and does not generate manuscript prose. Live
fetchers and an LLM-backed CV parser are deferred to post-6/13 unless an
in-flight Paper (1/2/3/11) requires them.
