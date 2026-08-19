---
title: Manuscript file inventory
audit_date: 2026-07-30
auditor: manuscript-audit-agent
---

# 00 — Manuscript file inventory

## Canonical submission files (v2 draft, latest)

| File | Size | Modified | Role |
|---|---|---|---|
| `manuscript_v8/NATURE_COMMUNICATIONS_FULL_DRAFT_v2_2026_07_08.md` | 42,757 B | 2026-07-08 | **PRIMARY CANONICAL DRAFT** — 8-figure NC target, integrates A2 interaction + IHC 3-plex + Lee replication. Use this as submission source. |
| `manuscript_v8/REVIEWER_FIRST_NC_REBUILD_2026_07_08.md` | 11,062 B | 2026-07-08 | Figure architecture specification — 5-figure reviewer-first plan. CONFLICTS with v2 8-figure structure (see BLOCKERS). |
| `manuscript_v8/SUBMISSION_CONTROL/DM1_PROJECT_STATE.md` | — | 2026-07-30 | Canonical project state register |
| `manuscript_v8/SUBMISSION_CONTROL/AUDIT_LOCKED_RESULTS.md` | — | 2026-07-30 | Immutable quantitative ledger |
| `manuscript_v8/SUBMISSION_CONTROL/CLAIM_LANGUAGE_MATRIX.md` | — | 2026-07-30 | Preferred/prohibited phrasing rules |
| `manuscript_v8/SUBMISSION_CONTROL/VOICE_PROTECTED_SLOTS.md` | — | 2026-07-30 | Author-only sections (V1–V6) |
| `manuscript_v8/SUBMISSION_CONTROL/NATURE_COMMUNICATIONS_CHECKLIST.md` | — | 2026-07-30 | NC submission requirements |
| `manuscript_v8/STATUS_INSILICO_FINALIZATION_2026_07_30.md` | 12,381 B | 2026-07-30 | Current project status in Korean/English |

## Section files (older v1 basis — partially superseded by v2)

| File | Size | Modified | Role |
|---|---|---|---|
| `manuscript_v8/01_abstract.md` | 5,990 B | 2026-05-08 | Abstract — v1 basis, does NOT reflect v2 interaction/IHC findings |
| `manuscript_v8/04_results.md` | 16,595 B | 2026-05-13 | Results — v1 basis (6-section), superseded by v2 sections §3-§9 |
| `manuscript_v8/05_figure_captions.md` | 36,019 B | 2026-05-13 | Cell Press 8-fig v3 captions — CRM-anchor build, not current NC target |
| `manuscript_v8/05_figure_captions_NC.md` | 27,379 B | 2026-05-27 | NC-architecture captions (6-main-figure) — does NOT match v2 8-figure structure |
| `manuscript_v8/06_discussion.md` | 9,899 B | 2026-05-08 | Discussion — v1 basis |
| `manuscript_v8/07_star_methods.md` | 16,767 B | 2026-05-13 | STAR methods — cohort edits done, 2x TODO placeholders remain |
| `manuscript_v8/08_cover_letter.md` | 3,764 B | 2026-05-08 | Cover letter — voice-protected slot unfilled |
| `manuscript_v8/09_reviewer_qa.md` | 15,107 B | 2026-05-13 | Reviewer Q&A (Q1-Q14) — Q9 voice-protected |
| `manuscript_v8/13_supplementary_tables.md` | 7,331 B | 2026-05-04 | Supplementary table shells |

## Supporting documents

| File | Size | Modified | Role |
|---|---|---|---|
| `manuscript_v8/NATURE_COMMUNICATIONS_FULL_DRAFT_2026_07_08.md` | 25,637 B | 2026-07-08 | v1 draft (5-figure) — superseded by v2; retain as reference |
| `manuscript_v8/NATURE_CANCER_MANUSCRIPT_ARCHITECTURE_2026_07_08.md` | 11,613 B | 2026-07-08 | NC/NM scope comparison — planning document only |
| `manuscript_v8/EDITORIAL_REVIEW_RUTHLESS_CUT_2026_06_05.md` | 21,976 B | 2026-06-05 | Reviewer ruthless-cut editorial audit |
| `manuscript_v8/MANUSCRIPT_REORG_6FIG_PLAN_2026_06_05.md` | 31,619 B | 2026-06-05 | 6-figure reorganization plan |
| `manuscript_v8/DM1_MASTER_BRIEF_FOR_GPT_AND_MEETING_2026_06_25.md` | 28,269 B | 2026-06-25 | Meeting brief |
| `manuscript_v8/email_reply_minsu_20260713.md` | 9,290 B | 2026-07-13 | Meeting correspondence |
| `manuscript_v8/kakao_minsu_20260714.md` | 2,813 B | 2026-07-15 | Kakao correspondence |

## Figure assets

### In `manuscript_v8/figures/` (only 2 files)

| File | Role |
|---|---|
| `figures/Fig7_dm1_mechanism.pdf/png` | Legacy Fig 7 mechanism panel |
| `figures/Fig8_epigenetic.pdf/png` | Legacy Fig 8 epigenetic panel |

### In `project/results/manuscript_v8_nc_main/` (main figure build directory)

| File(s) | Role |
|---|---|
| `Fig1_discovery_axis.pdf/png` + annotated | Main Fig 1 |
| `Fig2_fusion_mechanism.pdf/png` + annotated | Main Fig 2 |
| `Fig3_epigenetic.pdf/png` + annotated | Main Fig 3 |
| `Fig4_sc_validation.pdf/png` + annotated | Main Fig 4 |
| `Fig5_survival_portability.pdf/png` + annotated | Main Fig 5 |
| `Fig6_reflex_translation.pdf/png` + annotated | Main Fig 6 |
| `ED1–ED15_*.png` | Extended Data placeholders (15 files) |
| `EV_master_forest_annotated.png`, `EV_pergene_matrix_annotated.png` | External validation atlas |
| `External_Validation_Atlas.pdf/png/annotated` | Full validation atlas |

**MISSING from nc_main: Fig 7 and Fig 8** — these are in `project/dm1_story_web/public/figures/` and `dist/figures/` as:
- `fig_deep1_head_to_head.png` → v2 Fig 7 panel A
- `fig_deep2_prognostic_predictive.png` → v2 Fig 7 panel B
- `fig_a2_replication.png` → v2 Fig 7 panel C
- `fig_deep4_multiomics.png` → v2 Fig 7 panel D
- `fig_ihc3_killer.png` → v2 Fig 8
- `fig_deep5_snubh_simulation.png` → v2 Fig 8 panel E

**ALSO MISSING: Fig 3C (beta-expression scatter for TPO/DIO1/TSHR/TG)** — confirmed not built.
**ALSO MISSING: Fig 3D (per-driver-class beta bar)** — rebuild noted as needed.

### In `project/results/figures/` (additional assets)

| File | Role |
|---|---|
| `F3_driver_neutrality/F3_driver_neutrality.pdf/png` | Driver neutrality panel → Fig 1 |
| `F4_pangenome_robustness/F4_pangenome_robustness.pdf/png` | ARI ladder → Fig 1 |
| `suppl/SF1–SF6_*.pdf/png` | Supplementary figures |

## Archive / superseded files

| File | Role |
|---|---|
| `manuscript_v8/_archive_2026_05_03_sprint_draft.md` | Archived sprint draft |
| `manuscript_v8/10_full_manuscript_compiled.md` | Compiled v1 manuscript |
| `manuscript_v8/14_self_verification_report.md` | v1 self-verification |
| `manuscript_v8/_for_web_claude_*.md` | Review preparation notes |
| `manuscript_v8/PAPER2_*.md` | Paper 2 related planning documents |
| `manuscript_v8/cross_contamination_check_paper1.md` | Cross-contamination audit |
| `manuscript_v8/assets/p1_onepage_audit/` | One-page audit figures and TSVs |

## Critical structural observation

**The project currently has THREE parallel figure architectures in active use:**
1. v2 draft uses 8 main figures (Fig 1-6 from nc_main + Fig 7-8 from dm1_story_web)
2. `05_figure_captions_NC.md` describes a 6-figure NC architecture
3. `REVIEWER_FIRST_NC_REBUILD_2026_07_08.md` specifies a 5-figure reviewer-first plan

No single resolved, unified figure architecture exists. This is a FATAL BLOCKER for submission.
