# Full repo cleanup plan — post-Codex multi-session marathon
**Date:** 2026-05-06
**Trigger:** Codex sessions 89–94 produced ~159 untracked + 26 modified files spanning Paper 1 / 2 / 3 / 4 / 9.
**Mode:** ship-prep cleanup only — Claude does not author manuscript prose; voice-protected sections require explicit user approval before commit.
**Author:** Seungho Cook

---

## 0. Hard rules for this cleanup

- Claude does NOT write Hook / Aim / Discussion §3.1 / Limitations / Cover letter Para 1 / Reviewer Q9.
- Claude MAY git-stage and commit user-/Codex-authored prose if explicitly approved per-group.
- Per-paper commit groups so user can ship Paper 1 today and gate Papers 2/3/4/9 later if desired.
- 131 MB GEO download + 16 MB Korean HLA reference panel → gitignore (reproducible public refs).

---

## 1. Gitignore additions (recommended first)

| Path | Reason | Size |
|---|---|---|
| `project/results/aggressive_sprint_2026_05_06/geo/` | GEO download cache (GSE151179 family soft + GPL annot) | 131 MB |
| `project/external_refs/kim_korean_hla_ref/KOR_REF/` | Kim Korean HLA imputation reference (Beagle/PLINK binaries; public) | ~14 MB |
| `project/external_refs/kim_korean_hla_ref/HapMap3_CHB_JPT/` | HapMap3 reference (public) | ~250 KB |
| `project/external_refs/kim_korean_hla_ref/KOR_REF_1.0.zip` | zipped redundancy of above | 1.3 MB |
| `project/external_refs/korean_ngs_hla_controls/baek2021_plos/*.docx` | Baek 2021 Suppl. docx (publisher-hosted) | ~20 KB |
| `project/results/p_gpl570_validation/raw/` | already gitignored (line 75) | — |

**Action:** append the paths above to `.gitignore` and verify with `git check-ignore`.

Optionally KEEP committed (small, derivative): `project/external_refs/kim_korean_hla_ref/README.txt`, `project/external_refs/korean_ngs_hla_controls/baek2021_plos/parsed/*.txt`, `project/external_refs/korean_ngs_hla_controls/baek2021_plos/article.html` (227 KB — case by case).

---

## 2. Commit groups

Groups are ordered to allow stop-after-Paper-1 if desired. Each group is a single git commit unless noted.

### G1 — Paper 1 GPL570 validation pack *(my prior pass, still pending)*

**Files**
- `project/notebooks_or_scripts/p_gpl570_validation_first_pass.py`
- `project/scripts/build_gpl570_validation_html.py`
- `project/results/p_gpl570_validation/{sample_metadata.tsv, probe_to_gene_panel.tsv, GSE*_expression_gene_log.tsv.gz, GSE*_scores.tsv, gpl570_score_tests.tsv, gpl570_meta_effect_summary.tsv, gpl570_*.png}`
- `project/papers_hub_2026_05_04/gpl570_validation_pack.html`
- `project/reports/2026_05_04_gpl570_validation_access_check.md`
- `project/reports/2026_05_04_gpl570_external_validation_report.md`
- `project/reports/2026_05_04_gpl570_validation_FULL_VIEW.md`
- `.gitignore` (add `project/results/p_gpl570_validation/raw/`)

**Excluded:** `project/results/p_gpl570_validation/raw/` (~50 MB Series Matrices, gitignored).

### G2 — Paper 1 supporting reports + hub assets (Codex output)

**Files**
- `project/reports/2026_05_04_paper1_supp_methods_he_negative_feasibility.md`
- `project/reports/2026_05_04_paper1_venue_justification_gigatime.md`
- `project/reports/2026_05_04_session_recap_paper1_strategy_plus_landa_replication.md`
- `project/reports/2026_05_04_voice_hook_fact_brief.md`
- `project/reports/2026_05_06_paper1_web_integration_full_result.md`
- `project/papers_hub_2026_05_04/assets/paper1/` (10+ figures, ~10 MB total — fig_2x4_RAI_DM1_per_stage.png is 7.5 MB)
- `project/papers_hub_2026_05_04/paper1_nature_cancer_board.html`

**Note:** `paper1.html` itself is in the `M` set (existing tracked file modified by Codex). If approved, include in this group.

### G3 — Voice-protected manuscript_v8 prose *(user-authored via Codex — needs explicit OK)*

**Files (modified, voice-protected)**
- `project/manuscript_v8/00_title_candidates.md` (+14/−110)
- `project/manuscript_v8/01_abstract.md` (+1/−1)
- `project/manuscript_v8/03_introduction.md` (+5/−55)
- `project/manuscript_v8/04_intro_1_1_hook.md` (+14, no deletions)
- `project/manuscript_v8/04_results.md` (+10/−48)
- `project/manuscript_v8/05_figure_captions.md` (+62/−22)
- `project/manuscript_v8/06_discussion.md` (+22/−81)
- `project/manuscript_v8/07_star_methods.md` (+14/−32)
- `project/manuscript_v8/08_cover_letter.md` (+17/−34)
- `project/manuscript_v8/09_reviewer_qa.md` (+11/−25)
- `project/manuscript_v8/13_supplementary_tables.md` (+6/−6)

**New files**
- `project/manuscript_v8/05_supp_figure_captions_v17_dm.md`
- `project/manuscript_v8/10_full_manuscript_compiled.md`
- `project/manuscript_v8/14_self_verification_report.md`

**Hard requirement:** user explicit OK before commit. Claude did not author any of this prose.

### G4 — Paper 2 HLA work (Codex sessions 89–94)

**Scripts**
- `project/notebooks_or_scripts/paper2_hla_6allele_exploratory_stats.py`
- `project/notebooks_or_scripts/paper2_hla_aggressive_web_audit.py`
- `project/notebooks_or_scripts/paper2_hla_allele_vs_ngs_baseline.py`
- `project/notebooks_or_scripts/paper2_kim_hla_reference_panel_validation.py`
- `project/notebooks_or_scripts/paper2_paper4_hla_validation_and_meta_update.py`

**Reports** (8 total, dated 2026-05-06)
- `2026_05_06_paper2_6allele_exploratory_stats_report.md`
- `2026_05_06_paper2_6allele_forest_source_table_report.md`
- `2026_05_06_paper2_allele_vs_ngs_baseline_validation.md`
- `2026_05_06_paper2_hla_professor_figure_page.md`
- `2026_05_06_paper2_hla_web_integration_full_result.md`
- `2026_05_06_paper2_kim2014_hla_reference_panel_validation.md`
- `2026_05_06_paper2_korean_baseline_extraction_report.md`
- `2026_05_06_paper2_paper4_hla_professor_pages_full.md`
- `2026_05_06_hla_more_data_registry.md`
- `2026_05_06_hla_validation_and_screening_meta_update.md`
- `2026_05_06_korean_ngs_hla_control_bundle.md`

**Web/hub**
- `project/paper2-hla/index.html`
- `project/papers_hub_2026_05_04/paper2_hla.html`
- `project/papers_hub_2026_05_04/paper2_master_view.html`
- `project/papers_hub_2026_05_04/assets/paper2_hla/` (figures)

**Older P2 reports (May 4)**
- `project/reports/2026_05_04_HLA_strengthening_situation.md`
- `project/reports/2026_05_04_hla_beta_advisor_packet.md`
- `project/reports/2026_05_04_hla_beta_yu_send_package.md`
- `project/reports/2026_05_04_hla_strengthening_beta_plan.md`
- `project/reports/2026_05_04_paper1_paper2_scope_split_audit.md`
- `project/papers_hub_2026_05_04/hla_advisor_packet.html`
- `project/papers_hub_2026_05_04/hla_beta_plan.html`
- `project/papers_hub_2026_05_04/hla_situation.html`
- `project/papers_hub_2026_05_04/ht_vs_gd_audit.html`
- `project/papers_hub_2026_05_04/yu_send_package.html`
- `project/papers_hub_2026_05_04/pillar1_v2_spec.html`

**External refs** (committed metadata only after gitignore)
- `project/external_refs/kim_korean_hla_ref/README.txt`
- `project/external_refs/korean_ngs_hla_controls/baek2021_plos/parsed/{S9.txt, S10.txt}`

**manuscript_v8 internal P2 work**
- `_LEE_2014_FILL_DECISION_LOCK_2026_05_04.md`, `_PAPER2_MASTER_VIEW.md`, `PAPER2_PILLAR1_V2_FOREST_META_SPEC_2026_05_04.md`, `PAPER2_TERMINOLOGY_CORRECTION_PLAN_2026_05_04.md`, `_paper2_isolated_resume_response.md`, `_session_origin_paper2_isolated.md`, `_terminology_scan_hits.tsv`, `_terminology_scanner.py`

**manuscript_p2_brief**
- M `paper2_brief.html`, `paper2_brief.pdf`, `p2_advisor_discussion.html`
- `SITUATION_OVERVIEW.md`, `YU_MEETING_REVIEW_PACKET_2026_05_04.md`, `YU_MEETING_TALK_CARD_2026_05_04.md`

### G5 — Paper 4 GD HLA (Codex)

**Files**
- `project/notebooks_or_scripts/paper4_gd_hla_professor_assets.py`
- `project/paper4-hla/index.html`
- `project/papers_hub_2026_05_04/assets/paper4_hla/` (figures)
- `project/reports/2026_05_06_paper4_gd_hla_web_integration_full_result.md`

(Paper 2/4 cross-cutting reports already in G4: `2026_05_06_paper2_paper4_hla_professor_pages_full.md`)

### G6 — Paper 3 ICI (Codex)

**Files**
- `project/reports/2026_05_06_paper3_ici_trackA_plus_strengthening.md`
- `project/reports/2026_05_06_ici_vs_synthetic_lethality_boundary_map.md` *(also touches Paper 9)*

### G7 — Paper 9 synthetic lethality (Ruppin-style) (Codex)

**Files**
- `project/paper9-perturbation/index.html`
- `project/notebooks_or_scripts/post_pod_B_wsi_analysis.py`, `post_pod_C_alphafold_analysis.py`, `post_pod_D_meta3cohort.py`
- `project/scripts/p9_sl_first_pass.py`, `p9_sprint1_gls_validation.py`, `p9_sprint2_metabolic_prioritization.py`, `p9_final_web_packaging.py`, `p3_p9_full_public_execution.py`, `render_paper9_plan.py`, `render_paper9_results.py`
- `project/reports/2026_05_04_paper9_sl_first_pass_results.md`
- `project/reports/2026_05_04_paper9_synthetic_lethality_ruppin_style_plan.md`
- `project/reports/2026_05_06_paper9_synthetic_lethality_ruppin_style_strengthening.md`
- `project/papers_hub_2026_05_04/paper9_first_pass.pdf`, `paper9_plan.pdf`
- `project/results/paper9_sl_first_pass/`, `project/results/p3_p9_full_execution/`, `project/reports/p3_p9_full_execution/`

### G8 — Aggressive feasibility sprint (today)

**Files**
- `project/notebooks_or_scripts/aggressive_sprint_2026_05_06.py`
- `project/results/aggressive_sprint_2026_05_06/*.tsv`, `*.json`, `*.md` (excluding `geo/`)
- `project/reports/2026_05_06_advance_paper_feasibility_map.md`
- Top-level (move into project/reports/?): `dark_matter_feasibility_report.md`, `dark_matter_go_no_go_criteria.md`, `professor_summary_kr.md`

**Excluded:** `project/results/aggressive_sprint_2026_05_06/geo/` — 131 MB GEO download (gitignored).

### G9 — Pod B/C/D recovery + reports (May 4)

**Files**
- `project/notebooks_or_scripts/poll_and_recover_pod_D.sh`
- `project/reports/2026_05_04_pod_C_alphafold_raw.md`, `2026_05_04_pod_D_k2_star_postprocess.md`, `2026_05_04_pod_cd_completion_report.md`, `2026_05_04_bib_m3m4_low_batch_completion.md`
- `project/results/d7p3_k2_starred/`

### G10 — Hub web infrastructure

**Modified**
- `project/papers_hub_2026_05_04/index.html`
- `project/papers_hub_2026_05_04/paper1.html` (overlap with G2 — keep here unless G2 chosen alone)
- `project/8_papers_hub_2026_05_04.html` (−1198 lines — major restructure)
- `project/three_papers_index.html`, `.pdf`
- `project/submission/papers_overview.html`, `.pdf`
- `project/scripts/serve_secure.py`

**New**
- `project/papers_hub_2026_05_04/master_view.html`, `master_situation_2026_05_04.html`
- `project/papers_hub_2026_05_04/intro_strategy_audit.html`
- `project/papers_hub_2026_05_04/paper11_pancancer.html`, `.pdf`
- `project/papers_hub_2026_05_04/paper12_network.html`, `.pdf`
- `project/papers_hub_2026_05_04/paperB.html`, `paperG.html`, `paperN.html`
- `project/scripts/build_paper10_atlas_html.py`, `build_paper11_12_html.py`, `p10_atlas_figures.py`, `p11_p12_extended_atlas.py`
- `project/scripts/nature_cancer_feasibility_analysis.py`, `nature_cancer_new_data_analysis.py`, `nature_cancer_scrna_epithelial_filter.py`
- `project/results/paper10_atlas/`, `paper11_pancancer/`, `paper12_network/`, `nature_cancer_feasibility/`
- `project/results/p_external_expression_validation/`

### G11 — Strategy / session MDs + pod state

**Top-level**
- `GIGATIME_BUNDLE_2026_05_04.md`, `GIGATIME_STRATEGY_2026_05_04.md`
- `MARATHON_REENTRY_FINAL_2026_05_04.md`
- `PAPERS_5_OVERVIEW_2026_05_04.md`
- `POD_CD_RUNNING_STATE_2026_05_04.md`, `POD_CLEANUP_AND_PAPER1_LOCK_2026_05_04.md`
- `project/CODEX_PROMPT_8_PAPERS.md`

**Reports**
- `project/reports/2026_05_04_MASTER_SITUATION_VIEW.md`
- `project/reports/2026_05_04_intro_v1_strategy_audit.md`
- `project/reports/2026_05_04_voice_hook_audit_workflow.md`
- `project/reports/2026_05_04_paper2a_*` (3 files — Paper 2 reentry/methods/reviewer-Q)
- `project/reports/2026_05_04_supp_methods_closure_negative_feasibility.md`
- `project/results/terminology_correction_2026_05_04/SESSION_REPORT.md` + M files

**Misc**
- `thyroid_public_dataset_registry.tsv`
- `project_external_st/results/extra/g1_external_tile_metadata.tsv.gz`
- `project_external_st/src/12_gpu/`

### G12 — This cleanup plan

- `project/reports/2026_05_06_FULL_CLEANUP_PLAN.md` (this file)

---

## 3. Excluded from any commit (gitignored or out-of-scope)

| Path | Reason |
|---|---|
| `project/results/aggressive_sprint_2026_05_06/geo/` | 131 MB GEO download |
| `project/external_refs/kim_korean_hla_ref/{KOR_REF, HapMap3_CHB_JPT, *.zip}` | ~16 MB public ref |
| `project/external_refs/korean_ngs_hla_controls/baek2021_plos/*.docx` | publisher-hosted Baek 2021 supp |
| `project/results/p_landa_2016/raw/` | already gitignored |
| `project/results/p_gpl570_validation/raw/` | already gitignored after G1 |

---

## 4. Recommended sequencing

| Choice | What you say | Effect |
|---|---|---|
| **Ship Paper 1 only, hold rest** | `commit G1+G2 only` | Paper 1 GPL570 validation + Paper 1 hub assets shipped; everything else stays untracked. Safe minimum. |
| **Ship Paper 1 + manuscript prose** | `commit G1+G2+G3` | Adds the voice-protected manuscript_v8 edits (user-authored via Codex). |
| **Ship everything except voice-protected** | `commit G1, G2, G4, G5, G6, G7, G8, G9, G10, G11, G12` | All Codex outputs shipped; manuscript_v8 prose held for explicit approval. |
| **Full ship (transparent)** | `commit all groups` | Every group committed in order. Voice-protected prose included with explicit user authorization. |
| **Custom selection** | e.g. `commit G1, G3, G7` | Cherry-pick groups. |

Each commit will be made separately so history reads as a clean sequence by topic.

---

## 5. Pre-flight checks Claude will run before staging anything

1. Apply gitignore additions in §1 first; verify with `git check-ignore -v` for the listed large paths.
2. Confirm zero binary file >10 MB enters any staging area (defensive scan).
3. Confirm `project/results/aggressive_sprint_2026_05_06/geo/` is excluded from `git add`.
4. Confirm no `__pycache__` / `.pyc` slips in.
5. For G3 (voice-protected) only: re-display the diff stat lines and require an explicit second confirmation.

---

## 6. What Claude will NOT do

- Will not write Hook / Aim / Discussion §3.1 / Limitations / Cover letter Para 1 / Reviewer Q9 prose.
- Will not lightly polish or paraphrase voice-protected sentences.
- Will not commit Papers 2/3/4/9 without an explicit per-group approval.
- Will not push or open PRs unless told.

---

## 7. Approval phrases (next user input)

- `commit G1+G2 only`
- `commit G1+G2+G3`
- `commit G1, G2, G4, G5, G6, G7, G8, G9, G10, G11, G12`
- `commit all groups`
- Custom: `commit Gn, Gm, ...`
- Or: `audit X first` / `do not commit, change Y`

---

**Plan ready. No commits executed. Awaiting choice.**
