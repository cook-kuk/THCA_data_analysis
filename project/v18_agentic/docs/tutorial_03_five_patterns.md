# Tutorial 3 — The five patterns

Each pattern comes from a recurring shape in the v3 → v17p35 thyroid
sprint. Each is a small composable class — pick the ones you need.

---

## Pattern 1 — Fix-Amplify-Synthesis cycle

**Class:** `FixAmplifyAgent` (`patterns/fix_amplify.py`)

The shape:
- **Fix** — targeted intervention on a specific defect
- **Amplify** — extend the fix's surface area
- **Synthesis** — audit, write up, declare done or re-enter

```python
agent = FixAmplifyAgent(
    goal=Goal(description="...", max_iterations=6),
    fix_handler=my_fix,       # async (goal, last) -> dict
    amplify_handler=my_amp,
    synthesis_handler=my_syn,
    verifier=Verifier(checks=[...]),
)
result = await agent.run()
```

**Why it works:** prevents two failure modes. (a) Endless fixes without
amplification (perfectionism on one cohort). (b) Endless amplification
without synthesis (no one writes the paper).

---

## Pattern 2 — Tiered Parallel DAG

**Class:** `TieredParallelOrchestrator` (`patterns/tiered_dag.py`)

```python
orchestrator = TieredParallelOrchestrator([
    Tier(name="discovery", tasks=[
        Task("clustering", clustering_fn),
        Task("panel_screen", panel_fn),
    ]),
    Tier(name="validation", tasks=[
        Task("external", external_fn),
        Task("survival", survival_fn),
    ]),
])
summary = await orchestrator.run()
```

**Why it works:** within-tier parallelism gives you wall-time wins; tier
boundaries respect data dependencies. v17p35 Phase B ran 13 tasks across
4 tiers in ~2.5h vs ~9h serial.

---

## Pattern 3 — Dual Dump

**Class:** `DualDumpRecorder` (`patterns/dual_dump.py`)

```python
recorder = DualDumpRecorder("./workspace/dumps")
recorder.dump(
    name="phase_B_summary",
    raw_payload=results_dict,                 # full machine-readable payload
    compressed_md="# Phase B\n\n...",         # 60-second skim
)
```

**Why it works:** reviewers (human or LLM) read the compressed markdown.
Replays read the raw JSON/TSV. The two stay in lockstep because they're
written together.

---

## Pattern 4 — Reviewer Loop

**Class:** `ReviewerLoop` (`patterns/reviewer_loop.py`)

```python
loop = ReviewerLoop(
    reviewers=[mock_reviewer_LLM],
    defense_builder=my_defense_handler,
)
critiques, defenses = await loop.review(artifacts)

# Promote top attacks into the next iteration's goal
next_goal = ReviewerLoop.goal_from_attacks(
    critiques[0].top_attacks(3), parent=current_goal
)
```

**Why it works:** the highest-leverage way to escape local optima is to
hand the artifact to an external mind that doesn't share your assumptions.
v17p35 ran this loop on every Phase B output.

---

## Pattern 5 — Honest Failure Catalog

**Class:** `FailureCatalog` (`patterns/failure_catalog.py`)

```python
catalog = FailureCatalog("./workspace/failures.jsonl")
catalog.log(FailureRecord(
    failure_id="V1_TERT_RECOVERY",
    title="v1 TERT recovery missed thca_tcga_pub",
    what_was_tried="cBio: thca_tcga + thca_tcga_pan_can_atlas_2018",
    why_it_failed="Both use WES MAF; TERT promoter not captured by WES baits",
    next_action="Iterate ALL cBio thyroid studies",
    severity="HIGH",
))
print(catalog.render_markdown())
```

**Why it works:** the catalog becomes the limitations section of the paper
and the test fixture for "don't repeat this" agents in the next sprint.
v17p35 logged 11 failures; 7 became reviewer-defense entries.

---

## Combining patterns

The thyroid replay (`examples/thyroid_replay.py`) uses all five in one
script:

1. `FixAmplifyAgent` for the discovery loop.
2. `TieredParallelOrchestrator` for Phase B-style fan-out.
3. `DualDumpRecorder` to write the audit pair.
4. `ReviewerLoop` to ingest external critique.
5. `FailureCatalog` to record what didn't work.

Run it:

```bash
python -m agentic_research.examples.thyroid_replay --max-iter 3
```

Output is written to `/tmp/v18_thyroid_replay/`. Read the `dumps/*.md`
files first — that's where the five-pattern flow becomes legible.
