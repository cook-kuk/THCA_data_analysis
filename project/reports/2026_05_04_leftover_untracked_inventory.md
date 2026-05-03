# Leftover untracked / modified file inventory — post closure pivot

**Date:** 2026-05-04 · **Owner:** Seungho Cook
**Trigger:** Closure pivot commits `1d16a4e` + `c42bd6a` left several files in working tree that need decision. This memo classifies each.

**Rule:** No file deleted by Claude. No data/tile/embedding committed. No script committed without user explicit approval.

---

## Classification key

- **A** — commit as internal report (safe; Claude-staged this commit)
- **B** — leave untracked / gitignore (deprecated angle or runtime state; user may inspect/delete manually)
- **C** — needs manual user delete (only if user decides to clean)
- **D** — needs explicit user decision before any action

---

## Files

| path | size | nature | classification | reason |
|---|---|---|---|---|
| `PAPER1_DM1_FULL_2026_05_04.md` | 19,655 B | modified (+49 lines from 5/4 audit by parallel agent — adds top-of-file `INTERNAL WORKING DOCUMENT` banner + per-section `[MANUSCRIPT-SAFE candidate]/[INTERNAL-ONLY]/[DEPRECATED]` tags + new "Remaining risks (NOT all risks resolved)" section listing N1, Q3, H&E NO-GO, no wet-lab, narrow pan-cancer, voice-protected pending) | **D** | Modification is conservative (banner + tags, not prose). Recommend user review the diff (`git diff PAPER1_DM1_FULL_2026_05_04.md`) and commit if accepting. Claude does NOT auto-commit because PAPER1_DM1_FULL is voice-territory adjacent. |
| `SITUATION_BACKGROUND_ACTIVITY_2026_05_04.md` | 11,096 B | new — dated report describing `runpod_dispatch_v2.sh` loop, Pod C+D dispatch state, 9-file working-tree decision menu. Provides Q1/Q2/Q3 user-decision frame. | **A** | Report-shaped; provides historical record of background activity. Safe to commit. |
| `project/reports/2026_05_04_paper1_dm1_full_conflict_audit.md` | 22,828 B | new — claim-by-claim source verification (16/16 source TSVs verified) + identifies §9/§11 RunPod recommendations as superseded by closure NO-GO. Authority for the molecular-only lock. | **A** | Report. Safe to commit as parallel auditor's contribution. Referenced by the lock memo (`2026_05_04_paper1_dm1_full_molecular_only_lock.md` §1). |
| `project/notebooks_or_scripts/post_pod_B_wsi_analysis.py` | 4,722 B | new — Pod B output processor (WSI pathology embeddings → Hashimoto classifier) | **B** | Pod B = image-DM1 / pathology AI territory. Deprecated angle per closure NO-GO. Leaving untracked rather than committing (commit would normalize the deprecated workflow). User may keep on disk for inspection or delete manually. |
| `project/notebooks_or_scripts/post_pod_C_alphafold_analysis.py` | 3,848 B | new — Pod C output processor (AlphaFold structures) | **B** | AlphaFold is unrelated to image-DM1 closure but also outside Claude's authorized work scope this session. Leave untracked, user owns. |
| `project/notebooks_or_scripts/post_pod_D_meta3cohort.py` | 8,177 B | new — Pod D output processor (K2 STAR / meta 3-cohort) | **B** | Same as C. User owns. |
| `project_external_st/src/12_gpu/` (dir, 8 files: `setup.sh`, `extract_external_tiles.py`, `g1_embed.py`, `g1_ridge_loso.py`, `g2_slide_regress.py`, `g2_tcga_he_pipeline.sh`, `g2_tile_embed.py`, `g3_celldart.py`, `README.md`) | ~36 KB total | new — H&E-DM1 G1/G2/G3 GPU pipeline implementations + TCGA HE pipeline + CellDART | **B** | Direct image-DM1 pipeline implementations. Closure battery rejected this angle (`2026_05_04_image_dm1_final_nogo_decision.md`). Committing would normalize the deprecated workflow against marathon spec. Leave untracked. Do NOT delete (preserves option for §5 re-entry conditions if ever met). |
| `project_external_st/results/extra/g1_external_tile_metadata.tsv.gz` | 118,120 B | new — external-ST tile metadata for the G1 GPU pipeline | **B** | Output of deprecated angle scripts. Leave untracked. |
| `project/manuscript_p2_brief/p2_advisor_discussion.html` | 137,948 B | modified (2-line change) | **D** | Paper 2 brief HTML — author/other-agent territory. Claude does not modify Paper 2 brief HTML in this session. User reviews. |
| `project/results/terminology_correction_2026_05_04/before_after_diff.md` | 7,964 B | modified (+49 lines, Paper 2 Task B audit append) | **D** | Paper 2 Task B audit append by parallel agent. User reviews / commits with related Paper 2 work. |
| `project/results/terminology_correction_2026_05_04/files_modified.txt` | 603 B | modified (+5 lines) | **D** | Paper 2 Task B file-list update. Same as above. |
| `project/submission/papers_overview.html` | 32,150 B | modified (4-line change) | **D** | Top-level papers-overview HTML — author/other-agent territory. User reviews. |
| `project/results/03_pathology_poc/tile_metadata.tsv.gz` | (already tracked) | committed in `1d16a4e` | n/a | Closure outputs already in repo at clean state. Listed for spec completeness. |
| `project/results/03_pathology_poc/tile_metadata_resid.tsv.gz` | (already tracked) | committed in `1d16a4e` | n/a | Same. |
| `project/.runpod_ssh.json` (not in working tree status anymore) | (gitignored) | runtime state | n/a | Added to `.gitignore` in `1d16a4e`. Live SSH IP/port dictionary; do not commit. |

---

## Summary by classification

### A — committed by this session (with Claude's selective stage)

- `project/reports/2026_05_04_paper1_dm1_full_molecular_only_lock.md` (Claude-written, this session)
- `project/reports/2026_05_04_leftover_untracked_inventory.md` (this file)
- `project/reports/2026_05_04_paper1_dm1_full_conflict_audit.md` (parallel auditor; safe to include)
- `SITUATION_BACKGROUND_ACTIVITY_2026_05_04.md` (background activity report; safe to include)

### B — leave untracked (deprecated angle or non-Paper-1 workflow)

- `project/notebooks_or_scripts/post_pod_B_wsi_analysis.py`
- `project/notebooks_or_scripts/post_pod_C_alphafold_analysis.py`
- `project/notebooks_or_scripts/post_pod_D_meta3cohort.py`
- `project_external_st/src/12_gpu/*` (8 files)
- `project_external_st/results/extra/g1_external_tile_metadata.tsv.gz`

### C — none recommended for delete by Claude

(Claude does not delete. User may decide later whether to remove `12_gpu/`, `post_pod_*.py`, and `g1_external_tile_metadata.tsv.gz` once the marathon completes and §5 re-entry conditions are clearly not pursued.)

### D — needs user review / decision

- `PAPER1_DM1_FULL_2026_05_04.md` modification (banner + per-section tags + risks section by parallel agent — recommend accepting and committing on user review)
- `project/manuscript_p2_brief/p2_advisor_discussion.html` (2-line change)
- `project/results/terminology_correction_2026_05_04/before_after_diff.md` (Paper 2 Task B append)
- `project/results/terminology_correction_2026_05_04/files_modified.txt` (Paper 2 Task B append)
- `project/submission/papers_overview.html` (4-line change)

---

## Recommendation

1. **Accept this commit** (A items only) — safe, no data, no scripts, no manuscript prose modification.
2. **Review D items separately** at next marathon block (especially `PAPER1_DM1_FULL` which has the parallel agent's banner+tags ready).
3. **Defer B items** (do not commit, do not delete) — they document the deprecated angle's implementation but should not be normalized into the repo while the angle is closed under marathon mode.
4. **Do not run any B-classified scripts** during the marathon (image-DM1 / G1/G2/G3 pipelines are blocked).

---

## Cross-references

- `2026_05_04_image_dm1_final_nogo_decision.md` (decision authority)
- `2026_05_04_paper1_dm1_full_molecular_only_lock.md` (this session's lock)
- `2026_05_04_paper1_dm1_full_conflict_audit.md` (parallel auditor)
- `2026_05_04_marathon_state_after_image_dm1_nogo.md`
- `1d16a4e decision: drop image-dm1 pathology angle after closure battery`
- `c42bd6a fix: close paper1 marathon audit low batch L1-L6` (parallel)
- `778b255 decision: archive image-dm1 closure battery no-go` (parallel)
