---
title: "Paper 2 advisor brief — CHANGELOG"
date: 2026-05-01
---

# CHANGELOG — `p2_advisor_discussion.html`

This file tracks substantive changes to the advisor brief.

## [unreleased] — 2026-05-01 (autonomous infra batch)

### Added — companion docs (4 files)
- `DATA_SOURCES_INDEX.md` — 21 figures + 5 tables → underlying TSV/JSON path mapping (26/26 source paths verified ✅).
- `PROMPT_DECISION_LOG.md` — audit log of v4 audit integration: 1차 wrong scope (DM1 fusion+epigenetic) → user correction → 2차 HLA-only filter → 5 propagation actions checklist.
- `README.md` — directory navigation, Pillar status table, DM1↔DM2 boundary clarification, PDF re-generation command.
- `CHANGELOG.md` — this file.

### Added — PDF deliverable
- `p2_advisor_discussion.pdf` — 56-page printable / email-able. Headless Chromium (`/home/seungho/.cache/ms-playwright/chromium-1208/chrome-linux64/chrome`) `--virtual-time-budget=15000` allows Plotly JS to render. 28 numeric labels + Korean text fully embedded. 2.6 MB.

### Pillar I — PARTIAL → STRONG
- Hero tagline + byline date → `2026-05-03` (Pillar I STRONG 격상).
- Pillar grid: status `Partial` → `★ Strong (격상)`.
- § 3.1 entire rewrite — 8 sub-sections (3.1.1–3.1.8 + Methods + Discussion paste-ready paragraphs):
  - Figure 2 ★★★ — qualitative single-allele scatter → 6-allele 3-arm forest (paper-defining)
  - Figure 2b — random-effects DerSimonian-Laird pooled OR + Cochran I²
  - Figure 2c — Korean sub-cohort heterogeneity (K2 / Lee / GSE286332-PTC)
  - Figure 2d — 4-scenario sensitivity DPB1*05:01
  - Table 1.0 — sub-cohort heterogeneity matrix
  - Table 1.1 — 4-scenario sensitivity matrix
  - Table 1.2 — Pillar I 5-row honest disclosure
  - § 3.1.6 Methods (paste-ready, 99.78% verbatim with `methods_paragraph.md`)
  - § 3.1.7 Discussion (paste-ready, 97.81% verbatim with `discussion_paragraph.md`, +1 sentence advisor-friendly Cochran I² enhancement)

### Citation correction (Chen 2018 → Chu 2018)
- All occurrences of "Chen 2018" → "Chu X et al. 2018 J Med Genet 55:685–692, doi:10.1136/jmedgenet-2017-105146 (PMC 6161647)".
- Verified via PMC 6161647 web fetch 2026-05-03.
- Footer source documents updated.

### § 3.6 Mediation — autocorrelation residualization expansion
- 3.6.1 sub-section added — P5 8-gene × HLA-II residualization Table 2.0 (TCGA + GSE286332 dual-cohort).
- 3.6.2 sub-section added — pooled meta 2-cohort 8-gene RAI Cohen d = −1.78 [−2.00, −1.56], τ²=0.

### Discussion Q&A
- Added Q9 — citation propagation (Chen → Chu) advisor decision.
- Added Q10 — venue ladder reach (Cell Rep Med default vs Nat Commun reach with Pillar I STRONG).
- Added Q11 — DQB1*02:01 0/874 freq (biological signal vs arcasHLA typing artifact).
- Added Q12 — haplotype-level (DRB1-DQA1-DQB1 trio + DPB1 LD) future work.

### Gap 3 RESOLVED
- § 4 Gap 3 styled as resolved (green border, strikethrough title).
- Gap 3 was: "4-cohort random-effects forest meta — Pillar I 격상 조건". Resolved via Chu 2018 published summary statistics integration (no Taiwan/Japanese cohorts needed).

### Scenario B 확률 50% → 55% ↑
- Decision rule callout updated to reflect Pillar I STRONG strengthens default Cell Rep Med reach.

### Footer enhancements
- Companion docs (DATA_SOURCES_INDEX, PROMPT_DECISION_LOG, DEPRECATED) embedded in footer.
- "22 figures (4 신규), 6 paper-defining (red box)".

### Stats
- HTML: 1322 → 1586 라인 (+264 라인, +20%)
- File: 95 kB → 128 kB (+33 kB)
- Plotly figures: 18 → 21 (4 신규: Fig 2 재작성 + Fig 2b/c/d 추가)
- Paper-defining figures: 5 → 6 (Fig 2 ★★★ 추가)
- Tables: 2 → 5 (Tables 1.0, 1.1, 1.2, 2.0 신규)
- Discussion Q: 8 → 12 (Q9–Q12 신규)
- All 54 hard-coded data points (Pillar I figures) verified ✅ vs source TSV.

## [pre-2026-05-01] — original brief

- 18 Plotly figures, 5 paper-defining (Fig 6, 9, 12, 14, 16)
- 5 Pillars + Mediation, 8 Discussion Q
- Pillar I status: PARTIAL (qualitative direction only, quantitative meta deferred to Gap 3)
- Citation: "Chen 2018 Han Chinese GD" (incorrect — actual paper Chu 2018)

## Companion infra changes (outside `manuscript_p2_brief/`)

### `project/results/d4p1_panasian_meta/`
- Added `DEPRECATED.md` — partial deprecation (3 files deprecated, 5 KEEP including `korean_PTC_pool_n908.tsv` used by canonical script + manuscript S9).

### `project/notebooks_or_scripts/`
- `v17_D4P1_forest_meta.py` — added DEPRECATED warning header + Python `DeprecationWarning` runtime hook.
- `v17_paper2_pillar1_forest.py` — extended docstring with citation history + canonical status + 4 downstream consumers + memory pointer.

### `project/results/` — 7 new READMEs
- `p3_gse286332/README.md` (Pillar II)
- `d4p2_tcga_hashimoto_signature/README.md` (Pillar III TCGA arm)
- `d5p6_bcr_repertoire/README.md` (Pillar IV)
- `d6p7_dm1_subcluster/README.md` (Pillar V — cross-paper boundary noted)
- `d8c_dm1_subB_x_K2_NBNR/README.md` (Pillar V Korean transfer)
- `d3p5_pdm1_gradient/README.md` (§ 3.6 Mediation)
- `d8b_korean_replication/README.md` (Pillar III Korean arm)

Each README: YAML frontmatter (dir / pillar / status / generator), files inventory + brief usage cross-ref, key claims (정량 + p-value), companion link.

---

*Maintainer: Seungho Cook (kukshomr@gmail.com). Update on each substantive brief revision.*
