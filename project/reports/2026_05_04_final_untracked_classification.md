# Final untracked classification — marathon re-entry cleanup

**Date:** 2026-05-04 · **Owner:** Seungho Cook
**Trigger:** marathon re-entry final cleanup before Paper 1 voice-hook entry.
**Rule:** Claude does NOT delete files. Claude does NOT commit data, embeddings, deprecated-angle scripts, or runtime-state JSON.

---

## A — commit as report (this session, `f66a6bb` follow-up commit)

| path | size | nature |
|---|---|---|
| `SHARE_2026_05_04_SESSION_RECAP.md` | 9.9 KB | parallel agent's clean session recap (8 commits, 4-paper status). Safe shareable single-file. |
| `project/reports/2026_05_04_runpod_user_stop_reminder.md` | new | this session — user action authority for Pod C/D stop |
| `project/reports/2026_05_04_final_untracked_classification.md` | new | this file |
| `SITUATION_BACKGROUND_ACTIVITY_2026_05_04.md` | already committed in `f66a6bb` | (skip — no re-commit needed) |

---

## B — leave untracked / deprecated image-axis (do NOT commit, do NOT run, do NOT delete by Claude)

These are direct implementations of the closed image-DM1 / H&E pathology angle, or output processors for pods that ran in that workflow. Marathon scope forbids re-execution. They remain on disk for reference under the §5 re-entry conditions of `2026_05_04_image_dm1_final_nogo_decision.md` (none currently met).

| path | size | role |
|---|---|---|
| `project_external_st/src/12_gpu/setup.sh` | 1.3 KB | GPU env setup |
| `project_external_st/src/12_gpu/extract_external_tiles.py` | 4.7 KB | external-ST tile extraction |
| `project_external_st/src/12_gpu/g1_embed.py` | 4.5 KB | G1 ResNet50 embedding |
| `project_external_st/src/12_gpu/g1_ridge_loso.py` | 5.6 KB | G1 LOSO Ridge train |
| `project_external_st/src/12_gpu/g2_slide_regress.py` | 3.6 KB | G2 slide-level regression |
| `project_external_st/src/12_gpu/g2_tcga_he_pipeline.sh` | 2.5 KB | G2 TCGA H&E pipeline driver |
| `project_external_st/src/12_gpu/g2_tile_embed.py` | 5.0 KB | G2 tile embed |
| `project_external_st/src/12_gpu/g3_celldart.py` | 5.4 KB | G3 CellDART deconvolution |
| `project_external_st/src/12_gpu/README.md` | 2.9 KB | G1/G2/G3 doc |
| `project_external_st/results/extra/g1_external_tile_metadata.tsv.gz` | 116 KB | G1 metadata output |
| `project/notebooks_or_scripts/post_pod_B_wsi_analysis.py` | 4.7 KB | Pod B (image-DM1 WSI) output processor |
| `project/notebooks_or_scripts/post_pod_C_alphafold_analysis.py` | 3.8 KB | Pod C (AlphaFold) output processor — non-image-DM1, but parallel-pod workflow scope |
| `project/notebooks_or_scripts/post_pod_D_meta3cohort.py` | 8.0 KB | Pod D (K2 STAR meta-3-cohort) output processor — same |

→ **Marathon-mode action: none.** Post-marathon: user decides keep or delete.

---

## C — user decision (review before commit)

| path | size / change | reason for hold |
|---|---|---|
| `project/manuscript_p2_brief/p2_advisor_discussion.html` | 2-line modification | Paper 2 brief — paused paper, user reviews edits |
| `project/manuscript_p2_brief/paper2_brief.html` | 10-line modification | same |
| `project/manuscript_p2_brief/paper2_brief.pdf` | binary regen | same |
| `project/results/terminology_correction_2026_05_04/before_after_diff.md` | +49 lines | Paper 2 Task B audit append (parallel agent) |
| `project/results/terminology_correction_2026_05_04/files_modified.txt` | +5 lines | same |
| `project/results/terminology_correction_2026_05_04/SESSION_REPORT.md` | new | parallel agent — user reviews/commits with Paper 2 work |
| `project/submission/papers_overview.html` | 4-line modification | top-level papers index — user reviews |
| `project/submission/papers_overview.pdf` | binary regen | same |
| `project/three_papers_index.html` | 14-line modification | three-paper index — user reviews |
| `project/three_papers_index.pdf` | binary regen (large diff) | same |
| `project/manuscript_v8/_session_origin_paper2_isolated.md` | new | parallel agent — Paper 2 origin trace; user reviews |
| `POD_CD_RUNNING_STATE_2026_05_04.md` | new (root) | parallel agent's pod-state report; superseded by `runpod_user_stop_reminder.md` written this session, but user may want both for record |
| `POD_CLEANUP_AND_PAPER1_LOCK_2026_05_04.md` | new (root) | this-session deliverable from earlier; user may commit at next block |
| `PAPERS_5_OVERVIEW_2026_05_04.md` | new (root) | parallel agent — 5-paper schema snapshot; user reviews framing |
| `project/reports/2026_05_04_voice_hook_fact_brief.md` | new (this session, **uncommitted by spec STEP 5**) | factual bullets only, no Hook prose; ready for voice-hook author session |

---

## D — never commit (data / runtime / regenerable)

| pattern | reason |
|---|---|
| `project/data/raw/*`, `project/data/processed/*`, `project/data_raw`, `project/data_processed` | gitignored already; raw + per-sample h5ad + tile PNGs |
| `project/results/03_pathology_poc/embeddings_*.npz` | gitignored already; regen via `embed_resnet50.py` |
| `project/results/03_pathology_poc/simple_features.npy` | already tracked from prior work; do not re-add |
| `project/.runpod_ssh.json`, `project/.runpod_pods.json` | gitignored; runtime IP/port + pod ID dictionary |
| `*.h5ad`, `*.fastq.gz`, `*.bam`, `*.bai`, `*.cram`, `*.alignment.p` | gitignored; bioinformatics intermediates |
| `azure_gpu_phaseA/phaseA_gpu_pkg.tar.gz` | gitignored; 614 MB regen via `package_for_gpu.sh` |
| `project/results/figures_for_advisor/fig_2x4_RAI_DM1_per_stage.png` etc. | gitignored heavy raster |

---

## Summary by class

| class | count | total size |
|---|---|---|
| A — committed this session | 3 (after dedup) | ~25 KB |
| B — leave untracked deprecated | 13 | ~52 KB |
| C — user decision | 14+ | mixed |
| D — never commit | many | gitignored |

---

## Cross-references

- `2026_05_04_runpod_user_stop_reminder.md` — Pod C/D user action
- `2026_05_04_paper1_dm1_full_molecular_only_lock.md` — Paper 1 manuscript-safe envelope (commit `f66a6bb`)
- `2026_05_04_paper1_dm1_full_conflict_audit.md` — parallel auditor verification (commit `f66a6bb`)
- `2026_05_04_image_dm1_final_nogo_decision.md` — closure decision authority (commit `1d16a4e`)
- `2026_05_04_paper2_post_image_dm1_nogo_status.md` — Paper 2 scope cleanup (commit `1d16a4e`)
- `2026_05_04_marathon_state_after_image_dm1_nogo.md` — marathon snapshot (commit `1d16a4e`)
- `2026_05_04_voice_hook_fact_brief.md` — voice-hook fact briefing (this session, uncommitted, author-ready)
