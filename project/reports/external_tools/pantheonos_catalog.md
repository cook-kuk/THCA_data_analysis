# PantheonOS port catalog (THCA repo)

**Source**: `aristoteleo/PantheonOS` (BSD-2-Clause, v0.5.3 2026-04-22) + `aristoteleo/pantheon-cli`
**Pulled**: 2026-05-08 (commit @ main)
**Staged at**: `project/external/pantheonos/`
**Authority**: `v19_marathon_override_pantheonos_2026_05_08` (one-shot, voice-protected sections still locked)

PantheonOS is a multi-agent framework + 1,000+ skill marketplace for biomedical data science. We are NOT installing the runtime (NATS-based, heavy). We extract the **portable artifacts** — skills (code recipes), agent prompts (orchestration templates), and self-contained analysis modules — and adapt to our cohort.

---

## 1. Inventory pulled

### 1.1 `examples/` trajectories (8 — staged subset)

| Trajectory | Type | Reusable artifact | Our local path |
|---|---|---|---|
| `code_distillation` | Code-opt agent | (skipped — meta) | — |
| `evolution_batch_correction` | GA-evolved Harmony | algorithm + figures | (not staged; reviewer-reserve) |
| `evolution_gene_panel` | RL-based MERFISH panel | self-contained Python | `project/external/pantheonos/evolution_gene_panel/` |
| `evolution_topact` | Spatial cell-type caller | TopACT classifier + GA wrapper | `project/external/pantheonos/evolution_topact/` |
| `fastq_processing` | FASTQ → counts | (we already use STAR; skipped) | — |
| `paper_reporter` (v1) | Bio-agent for paper drafting | bioagent.md prompt | `project/external/pantheonos/paper_reporter_v2/paper_reporter_v1.py` |
| `paper_reporter_v2` | **Literature search agent** (not a manuscript writer — clarification) | Tasks/Search agent with duckduckgo+web_crawl | `project/external/pantheonos/paper_reporter_v2/paper_reporter_v2.py` |
| `single_cell_spatial_analysis` | sc/spatial multi-agent team | 6 agent prompts + 2 skills + omics_expert_team.py | `project/external/pantheonos/single_cell_spatial_analysis/` |

### 1.2 `single_cell_spatial_analysis` sub-cases (6, prompt-only)

`als_motor_cortex`, `chd_heart_spatial` (fetal heart MERFISH 2D+3D), `gene_panel_selection_immune`, `mouse_embryo_e9.5`, `pbmc3k`, `stero_seq_embryo`. None ported — these are demos for non-thyroid tissues. The `chd_heart_spatial/detailed_prompts/sc_mapping.md` is worth reading as a worked example of the MOSCOT mapping skill.

### 1.3 Reusable skills (portable code recipes)

| Skill | What it does | Deps | Status |
|---|---|---|---|
| `skill_sc_spatial_mapping.md` | MOSCOT optimal-transport sc → spatial mapping; gene + celltype imputation | `pip install moscot` | Ready (paper9 candidate) |
| `skill_3d_viz.md` | PyVista 3D point-cloud + GIF rotation for spatial expression / celltype | `pip install pyvista` | Ready (Visium has 2D coords; 3D needs serial sections) |

### 1.4 Reusable agent prompts (orchestration templates)

`agent_leader.md` (orchestrator), `agent_analysis_expert.md` (scverse QC→PCA→UMAP→Leiden→DEG→annotation playbook), `agent_biologist.md` (hypothesis generator), `agent_reporter.md` (PDF report writer), `agent_browser_use.md` (web search), `agent_system_manager.md` (env install).

→ These are **prompt templates we can crib** for our v18 agentic_research framework, not a runtime to run.

### 1.5 Self-contained Python modules

| Module | Lines | What it is |
|---|---|---|
| `evolution_topact/topact_classifier.py` | ~3.7K | TopACT spatial cell-type classifier (independent of evolution wrapper) |
| `evolution_topact/topact_spatial.py` | ~17.7K | Core spatial logic |
| `evolution_topact/run_evolution.py` | ~7.4K | GA driver (requires PantheonOS runtime) |
| `evolution_gene_panel/rl_gene_panel.py` | ~25.5K | RL panel optimizer (PyTorch) |

→ TopACT classifier itself is portable. The GA/RL drivers need PantheonOS — skip those, use the bare classifiers.

---

## 2. Applicability to our papers

| Paper | Status | Best PantheonOS port | Marathon judgment |
|---|---|---|---|
| Paper 1 (DM1 dark matter) | manuscript-write marathon | `paper_reporter_v2` (literature mining for Discussion related-work) | ✅ infra, allowed |
| Paper 2 (H&E → DM1) | PASS_LAUNCH 2026-05-08 | `evolution_gene_panel` for IO MERFISH design once we move to a panel-validation phase | Backlog (not paper-blocking) |
| Paper 3 (ICI vulnerability) | Track A FROZEN, B-lite done | (none directly) | Track B still BLOCKED |
| Paper 9 (paper9-perturbation, current branch, 100 spatial figs) | active | TopACT classifier on Visium + MOSCOT skill if we get matched sc | Override-permitted scaffolding |
| Paper 11 (pan-cancer) | active | `paper_reporter_v2` (related work for Nat Commun reach) | ✅ infra |
| Paper 4 (GD HLA backlog) | gated 4/4 | (none) | — |

---

## 3. Recommended ports (ranked)

1. **paper_reporter_v2 → literature mining helper** (Paper 1 + Paper 11). Repurpose: feed each paper's theme, retrieve related Nat-Cancer/JCI Insight cites. Output → Discussion bib enrichment. **Voice-protected sections untouched.** Implementation: thin wrapper that uses our existing `requests`/`scholarly`/`crossref` instead of PantheonOS runtime — same prompt, our infra.

2. **MOSCOT skill → DM1 sc-spatial bridge**. We DO have matched data:
   - sc reference: `project/results/v17_lu2023/GSE193581_hvg_adata.h5ad` (Lu 2023, 67,678 cells × 2000 HVG, `author_celltype`)
   - spatial: GSE250521 Visium 15+ samples at `project/data/processed/GSE250521/`
   - alt sc: `/data/thca/scrna/processed/classical_baseline_F12.h5ad` (17.9k, mutation-stratified)

3. **TopACT classifier on Visium → paper9 alternative cell-type call**. Compare to existing niche clustering (S_F31, S_F47). Self-contained, no PantheonOS runtime needed. **Stub script staged; user picks Visium sample to run.**

4. **PyVista 3D viz → only if we acquire serial Visium sections**. Currently 2D — defer.

5. **evolution_gene_panel → Paper 2 panel design future**. Not today.

---

## 4. What's been done vs what's pending

✅ Done (this session):
- 31 source files staged at `project/external/pantheonos/`
- Catalog written (this file)

🟡 Stubbed (user input needed before running, 4 scripts):
- `scripts/external/pantheonos_lit_search.py` — Crossref-based lit search, no PantheonOS deps. Needs `--theme` per paper. Test: `--theme "DM1 thyroid dark matter" --out project/reports/lit_search/dm1.md`
- `scripts/external/pantheonos_lr_paper9.py` — Squidpy LR analysis (Gallery #6). Needs Visium .h5ad + cluster column. Reqs: `pip install squidpy omnipath`
- `scripts/external/pantheonos_topact_visium.py` — TopACT spatial cell-type. Needs Visium + sc reference (sc not currently available for THCA). Imports staged TopACT classifier
- `scripts/external/pantheonos_moscot_stub.py` — MOSCOT OT mapping. Reviewer-reserve only — no matched sc+spatial pair yet. Reqs: `pip install moscot`

❌ Not ported:
- PantheonOS runtime (NATS, multi-agent chatroom) — too heavy, not paper-blocking
- evolution_batch_correction, evolution_gene_panel runtime — reviewer-reserve only
- 6 sub-case demos — non-thyroid tissues

---

## 5. Marathon footnote

This port is authorized as a one-shot infra burst per `v19_marathon_override_pantheonos_2026_05_08`. Default ("새 분석 = paper-blocking only") resumes after this session. Voice-protected sections (Hook/Aim/Disc 3.1/Limitations/Cover Para 1/Q9) are NOT touched by anything in this catalog.
