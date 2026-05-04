# v15 NeurIPS 2025 — Editorial update (2026-05-04, post-v52 audit + closure battery)

**Companion document for** `v15_neurips_SUBMIT_FINAL.pdf`
**Date:** 2026-05-04
**Author:** Seungho Cook
**Source paper:** v15 NeurIPS 2025 submission, audit_report 2026-04-27 (12 pages, 345 KB, anonymized)
**Hub view:** http://40.82.129.113/papers_hub_2026_05_04/paper6.html

> **Why this note exists.** The v15 NeurIPS submission was complete and submitted on 2026-04-27 with the original v5.1 THCA DIAL = 0.494 / 4-of-5 batch_entangled finding. A subsequent forensic re-audit (v52 LODO ComBat, 2026-04-29) showed that the 0.494 figure was a normalization-leak artifact, not a true cross-cohort generalization signal. This document records the post-submission state. The PDF itself is not modified (submission artifact).

---

## 1. Submission state (frozen 2026-04-27)

| field | value |
|---|---|
| Submission file | `v15_neurips_SUBMIT_FINAL.pdf` (also `v15_NEURIPS_DASHBOARD.pdf`) |
| Date | 2026-04-27 (audit_report final compile) |
| Pages | 12 |
| Size | 345 KB anonymized |
| Audit components on disk | claims_audit · prior_work_positioning_table · openreview_form_fields · reproducibility_checklist_filled · reproducibility_verification · self_review_simulation · rebuttal_prep_v1 · camera_ready_prep_plan · deadline_calendar · countdown_apr27_to_may06 · anonymous_code.zip · v15_NEURIPS_SUBMIT_DUMP |
| Status | Submitted; no withdrawal |

## 2. Post-submission discovery — v52 LODO leak finding (2026-04-29)

| event | date | finding |
|---|---|---|
| Pre-submit audit | 2026-04-27 | audit_report 2026-04-27 final compile; v5.1 DIAL = 0.494 cite as cross-cohort generalization signal |
| v17 audit phase | 2026-04-28 | normalization order reviewed — train/test split 가 ComBat 적용 후 였을 가능성 의심 |
| v52 forensic re-audit | 2026-04-29 | 5/5 THCA classifier under LODO ComBat (Leave-One-Domain-Out, fit per fold on train only) → all DIAL = 0.000 / true_biology |
| Halt decision | 2026-04-29 | downstream paper (Paper 1, Paper 2A/B) 들에서 v5.1 0.494 number cite 금지 |
| Memory anchor | 2026-04-30 | `v52_lodo_finding` 작성 |

**Key result:** The v5.1 THCA DIAL = 0.494 number reported in the v15 NeurIPS submission was a **normalization-leak artifact**, not a true cross-cohort generalization signal. Under proper LODO ComBat (cohort-pair-aware batch correction) all 5/5 THCA classifiers return DIAL = 0.000. The DIAL **methodology** (defined in submission Sections 3–4) is not invalidated; the THCA-specific **empirical headline** is.

## 3. What survives in the submission text

| section | survives v52 | comment |
|---|---|---|
| Theorem 1 (phase transition) | ✅ unchanged | Theoretical result |
| Theorem 2 (extension) | ✅ unchanged | Theoretical result |
| Theorem 2 proof audit | ✅ unchanged | Audit doc preserved |
| Methodology framework (Sections 3–4) | ✅ unchanged | DIAL definition + computation |
| Synthetic benchmark | ✅ unchanged | Controlled simulation |
| Quantum bound proposition | ✅ unchanged | Theoretical |
| **THCA-specific cancer headline (Section 5/6)** | ❌ **superseded** | v5.1 → v52: DIAL 0.494 → 0.000 |
| Cross-cohort generalization claim for THCA | ❌ superseded | Replaced by "consistent with Theorem 1 prediction under controlled cohort separation" |

## 4. Reproducibility — original vs corrected

| protocol | result | reproducibility path |
|---|---|---|
| v5.1 (leaky ComBat: fit on pooled X, then split) | DIAL = 0.494 across 4-of-5 classifiers | `results/v5/dial_results.tsv` (preserved) |
| v52 (proper LODO ComBat: fit per fold on train only) | DIAL = 0.000 across 5/5 classifiers | `results/v5p2_fix/v5p2_dial_proper_lodo.tsv` |
| Permutation test (domain label shuffle) | DIAL distribution centered at 0 | `results/v5p2_fix/permutation_dial.tsv` |
| kBET acceptance rate | v5.1: ~0.45 / v52: ~0.92 | `results/v5p2_fix/kbet_post.tsv` |
| iLISI score (cohort) | v5.1: ~1.4 / v52: ~1.9 | same |

**Reviewer instruction:** Both `results/v5/` (original v5.1) and `results/v5p2_fix/` (corrected v52) trees are preserved on disk. The corrected v52 result (DIAL = 0.000) is consistent with Theorem 1's phase-transition prediction under strict cohort separation.

## 5. Cross-paper references (THCA project as of 2026-05-04)

| Paper | Topic | Status | Hub |
|---|---|---|---|
| Paper 1 | 8-gene RAI/DM1 molecular axis | active manuscript_v8; bioRxiv target 2026-06-13 | http://40.82.129.113/papers_hub_2026_05_04/paper1.html |
| Paper 2A | H&E → DM1 closure battery | NO-GO verdict (D 폐기); methods supplement asset | http://40.82.129.113/papers_hub_2026_05_04/paper2a.html |
| Paper 2B | Hashimoto-overlap PTC mechanism | active brief; Pillar I-VI documented | http://40.82.129.113/papers_hub_2026_05_04/paper2b.html |
| Paper 3 | ICI vulnerability dark thyroid cancer | Track A frozen (chmod 444); Track B gated | http://40.82.129.113/papers_hub_2026_05_04/paper3.html |
| Paper 4 | Korean GD HLA / Pan-Asian | backlog; gated by Bundang FFPE access | http://40.82.129.113/papers_hub_2026_05_04/paper4.html |
| Paper 5 | npj Precision Oncology 8-gene biomarker | ship-ready | http://40.82.129.113/papers_hub_2026_05_04/paper5.html |
| **Paper 6** | **this v15 NeurIPS DIAL paper as v52 forensic re-audit case study** | **methods paper candidate (Briefings in Bioinformatics)** | **http://40.82.129.113/papers_hub_2026_05_04/paper6.html** |
| Paper 7 | v18 agentic_research framework | Korean Bioinformatics Society talk + methods paper | http://40.82.129.113/papers_hub_2026_05_04/paper7.html |

## 6. Position relative to recent multimodal AI

GigaTIME (Valanarasu et al., **Cell** 189[2]:386–400.e19, January 2026, DOI 10.1016/j.cell.2025.11.016) translates H&E → virtual multiplex-IF using ~40 × 10⁶ paired cells across 14,256 patients and 24 cancer types. This work is complementary to our DIAL methodology (different task class) and provides a useful task-class boundary reference for the closure battery NO-GO documented in Paper 2A. See `project/reports/2026_05_04_gigatime_manuscript_use_only.md` (commit `a1fcac3`) for the manuscript-use-only utility.

## 7. Recommended forward path

| action | who | when |
|---|---|---|
| Withdraw v15 from NeurIPS if review is in-progress (assess based on submission cycle) | author | TBD |
| Repackage as Paper 6 methods paper for **Briefings in Bioinformatics** (best fit, IF ~6.8) or **GigaScience** (open data, IF ~6.0) | author | post-marathon (post-2026-06-13) |
| Co-publish with Paper 7 (v18 agentic_research framework as Pattern 1 archetype) | author | concurrent |
| Author Tier 1 needed data for methods paper | author / Yu / external collab | 6-12 months |
| File pre-registration (OSF) for any future cross-cohort transcriptional classifier study | author | template ready |

## 8. Anchor files (current state)

- `project/submission/v15_neurips/audit_report.md` — original NeurIPS audit (2026-04-27)
- `project/submission/v15_neurips/v15_neurips_SUBMIT_FINAL.pdf` — submitted artifact (frozen)
- `project/submission/v15_neurips/v15_neurips_update_note_2026_05_04.md` — **this file**
- `project/reports/2026_05_04_paper1_dm1_full_molecular_only_lock.md` — Paper 1 manuscript-safe envelope (commit `f66a6bb`)
- `project/reports/2026_05_04_image_dm1_final_nogo_decision.md` — closure decision authority (commit `1d16a4e`)
- `project/reports/2026_05_04_gigatime_manuscript_use_only.md` — GigaTIME utility (commit `a1fcac3`)
- Memory anchors: `v52_lodo_finding`, `v17_marathon_mode_post_pillar1`, `v17_sprint_vs_marathon_violation`

---

*Generated 2026-05-04. The v15 NeurIPS PDF itself is not modified (submission artifact). This document is a companion update note for reviewers + cross-paper navigation. For the broader THCA project state, navigate to the paper6 hub view.*
