# Pod cleanup + Paper 1 DM1_FULL molecular-only lock — session report

**Date:** 2026-05-04 · **Owner:** Seungho Cook · **Session commit:** `f66a6bb`
**Files written this session:** 4 (all under `project/reports/` + 1 root memo)
**Files NOT modified this session:** PAPER1_DM1_FULL (parallel agent has banner+tags ready), Paper 2 brief HTML, papers_overview.html, all `[B]` deprecated-angle scripts

---

## 1. Background process status

### Local
- No active Python compute (process scan empty for `runpod|g1_|g2_|g3_|embed|ridge|tcga|WSI|CellDART|spaceranger|train_loso|closure|negative_controls|extract_tiles|compute_resid|simple_features|post_pod|12_gpu`)
- Dispatch loop `/tmp/runpod_dispatch_v2.sh` 자연 종료 at **01:27:38** (last iter=6, dispatched 2/3)

### Remote (RunPod, observed read-only — NOT killed per spec)
| Pod | endpoint | role | shutdown |
|---|---|---|---|
| Pod C | 216.81.151.3:10960 | AlphaFold | 24hr forced auto-shutdown set |
| Pod D | 157.157.221.29:21123 | K2 STAR meta3cohort | 24hr forced auto-shutdown set |

→ **Self-terminate within 24 h.** No Claude action. User decides if early stop needed.

---

## 2. Molecular-only lock verdict

**Paper 1 = molecular DM1 axis full GO.**
**Image-DM1 / RunPod G1-G2-G3 / TCGA WSI / IP angle full STOP.**

### Manuscript-safe envelope (16 source-verified claims)

| # | claim | placement |
|---|---|---|
| 1 | TCGA bulk DM1 vs THYROID_NONOVERLAP r=−0.885 (p=5.1×10⁻¹⁸⁸, n=561) | **Main** validation |
| 2 | TCGA DM1 cluster vs DM1_like score (DM1=0.00, DM2=−0.66, p≪0.001, n=157) | **Main** anchor |
| 3 | PFI dichotomized HR=2.04 [1.15–3.61] p=0.015; DFI multivar HR=1.41 p=0.025 (n=461) | **Main** KM + Cox forest |
| 4 | TF mechanism — FOXE1 −1.21 / NKX2-1 −0.65; STAT3 +1.55 / FOSL1 +0.98 / DNMT1 +0.74 (n=261) | **Main** mechanism panel |
| 5 | TROP2/TACSTD2 Δmean +4.33, FDR=1.2×10⁻³² (n=261) | **Main** volcano + drug overlay (tumor-level only) |
| 6 | Tumor vs normal — DM1 p=1.8×10⁻²⁰, TROP2 1.4× tumor>normal (n=517) | **Main** sub-panel |
| 7 | Hallmark GSEA — IL6_JAK_STAT3 +8.0 (p=6×10⁻¹⁵), OXPHOS −5.5 (p=1×10⁻⁷) | **Main/Supp** mechanism |
| 8 | TF network coordination (FOXE1-PAX8 r=+0.55, etc.) | **Supp** |
| 9 | Bootstrap 95% CI 12-slide [−0.997, −0.916] | **Supp** |
| 10 | 8-gene drop-one-out r ∈ [−0.987, −0.981] | **Supp** referee Q |
| 11 | Moran's I per slide mean 0.36, 26/28 perm p<0.05 | **Supp** spatial structure |
| 12 | Bivariate Moran (DM1 × NONOVERLAP) all 12 negative | **Supp** spatial anti-corr |
| 13 | Margin-distance gradient — 23/28 negative ρ | **Supp** |
| 14 | Pan-cancer THCA r=−0.892 outlier (33 cancers) | **Supp** thyroid-specific |
| 15 | N1 dark-matter (BRAF-neg ∩ RAS-neg, n=156) HR=1.20 NS | **Supp** honest negative |
| 16 | Harmony 28-slide UMAP (ST QC backbone) | **Supp** Methods |

### Excluded / deprecated (must NOT appear in submission text)

- H&E → DM1 prediction at tile level
- TCGA WSI external validation of DM1
- RunPod G1/G2/G3 pipeline
- Image-DM1 IP / patent claim
- "All risks resolved" / "TDS overlap, inflammation 모두 해소"
- "Cancer Cell 도전권 20-30%" / "Nat Cancer 50-65% reach" / venue-probability framing
- "DM1-high spots are TROP2-high" / spot-level colocalization (Q3 ρ=−0.016)
- "H&E-inferable molecular subtype"

### Required cautious language

| context | preferred phrasing |
|---|---|
| outcome association | `supports independent prognostic association` (NOT `proves`) |
| TF mechanism | `consistent with a 4-step framework` (NOT `demonstrates`) |
| TROP2 | `tumor-level therapeutic vulnerability` / `tumor-population vulnerability` |
| spatial transcriptomics | `supportive spatial validation` |
| wet-lab status | `functional validation remains required` (must appear in Limitations) |
| H&E negative | "We pre-tested H&E-based prediction... and found no learnable signal beyond random panels..." |

### Forbidden language (STOP if seen in draft)

`proves` · `establishes` · `all risks resolved` / `모든 risk 해소` · `Cancer Cell-ready` / `Cancer Cell 도전권` / `Nat Cancer reach` · venue-probability percentages · `H&E-inferable` / `H&E predicts DM1` · `pathology-AI triage of DM1` · `TROP2 spatially colocalizes with DM1-high spots` / `DM1-high spots are TROP2-high` · `tile-level DM1 inference` / `WSI-validated DM1` · `RunPod G1/G2/G3` (in active-claim form) · `morphology-derived DM1`

---

## 3. Leftover inventory classification

| classification | count | items |
|---|---|---|
| **A** committed (this session, `f66a6bb`) | 4 | molecular_only_lock · leftover_inventory · conflict_audit · SITUATION_BACKGROUND_ACTIVITY |
| **B** leave untracked (deprecated angle, do NOT commit) | 11 | post_pod_B/C/D + 8 files in `project_external_st/src/12_gpu/` + `g1_external_tile_metadata.tsv.gz` |
| **C** Claude does not delete | 0 | (none) |
| **D** user decision | 5+ | PAPER1_DM1_FULL diff (banner+tags by parallel agent), `p2_advisor_discussion.html`, `terminology_correction_2026_05_04/*` (P2 Task B append), `papers_overview.html`, `manuscript_v8/_session_origin_paper2_isolated.md` |

### B-item details (image-DM1 deprecated angle — leave on disk, do NOT commit, do NOT run)

```
project/notebooks_or_scripts/post_pod_B_wsi_analysis.py        # Pod B = image-DM1
project/notebooks_or_scripts/post_pod_C_alphafold_analysis.py  # AlphaFold output processor
project/notebooks_or_scripts/post_pod_D_meta3cohort.py         # K2 STAR output processor
project_external_st/src/12_gpu/setup.sh
project_external_st/src/12_gpu/extract_external_tiles.py
project_external_st/src/12_gpu/g1_embed.py
project_external_st/src/12_gpu/g1_ridge_loso.py
project_external_st/src/12_gpu/g2_slide_regress.py
project_external_st/src/12_gpu/g2_tcga_he_pipeline.sh
project_external_st/src/12_gpu/g2_tile_embed.py
project_external_st/src/12_gpu/g3_celldart.py
project_external_st/src/12_gpu/README.md
project_external_st/results/extra/g1_external_tile_metadata.tsv.gz
```

→ marathon 종료 후 5 re-entry conditions (`2026_05_04_image_dm1_final_nogo_decision.md` §5) 평가 시 reference 가능. Marathon 동안에는 실행/commit 금지.

---

## 4. Commit hash

```
f66a6bbd4a8af60753d2e7e887af650446563ada
docs: lock paper1 molecular-only interpretation after dm1 full audit

4 files changed, 708 insertions(+)
- project/reports/2026_05_04_paper1_dm1_full_molecular_only_lock.md  (lock authority)
- project/reports/2026_05_04_leftover_untracked_inventory.md         (A/B/C/D class)
- project/reports/2026_05_04_paper1_dm1_full_conflict_audit.md       (parallel agent's verification)
- SITUATION_BACKGROUND_ACTIVITY_2026_05_04.md                        (background activity report)
```

### Recent commit chain (today)

```
f66a6bb docs: lock paper1 molecular-only interpretation after dm1 full audit   ← THIS
c42bd6a fix: close paper1 marathon audit low batch L1-L6                       (parallel)
1d16a4e decision: drop image-dm1 pathology angle after closure battery
778b255 decision: archive image-dm1 closure battery no-go                      (parallel)
e4bb222 chore: resolve bib medium issues and runpod stragglers                 (parallel)
25d693c chore: archive spatial and pathology feasibility artifacts             (parallel)
ef023dd fix: close paper1 marathon audit high medium issues                    (parallel)
c55d582 freeze: paper3 ici track a design bundle                               (parallel)
```

---

## 5. Remaining user decisions

1. **`PAPER1_DM1_FULL_2026_05_04.md` diff (D)** — parallel agent has added top-of-file `INTERNAL WORKING DOCUMENT` banner + per-section `[MANUSCRIPT-SAFE candidate]` / `[INTERNAL-ONLY]` / `[DEPRECATED]` tags + new "Remaining risks (NOT 'all risks resolved')" section. Recommend accept + commit:
   ```bash
   git diff PAPER1_DM1_FULL_2026_05_04.md   # review
   git add PAPER1_DM1_FULL_2026_05_04.md
   git commit -m "docs: tag paper1 dm1_full bundle as internal molecular-only"
   ```

2. **Paper 2 brief modifications (D)** — bundle reviewed before commit:
   - `project/manuscript_p2_brief/p2_advisor_discussion.html` (2 lines)
   - `project/results/terminology_correction_2026_05_04/before_after_diff.md` (P2 Task B append)
   - `project/results/terminology_correction_2026_05_04/files_modified.txt` (P2 Task B)
   - `project/results/terminology_correction_2026_05_04/SESSION_REPORT.md` (new)
   - `project/submission/papers_overview.html` (4 lines)

3. **`project/manuscript_v8/_session_origin_paper2_isolated.md`** — new file by parallel agent. Read 후 commit/discard.

4. **B items (deprecated angle scripts) — keep on disk?** — marathon 끝나고 5 re-entry conditions 평가 시 reference. 안 쓸 거면 marathon 종료 후 manual delete.

5. **RunPod Pod C / Pod D** — 24h auto-shutdown OK인지, 즉시 종료 원하시는지.

---

## 6. Recommended next action

1. **voice-hook ready** → Paper 1 voice-protected (Hook / Aim / Discussion 3.1 / Limitations / Cover Para 1 / Q9) — **본인 키보드 only** (Claude X per `v17_sprint_vs_marathon_violation` memory)
2. **cosmetic cleanup 먼저** → 위 D items 5개 review + 묶어서 commit → 그 후 voice-hook
3. **safe-to-defer Claude work**: Paper 1 supplementary Methods scaffolding for closure battery negative-feasibility entry (Claude OK, voice-protected 아님)

---

## scp 명령 (Windows PowerShell)

```powershell
# 이 단일 요약 파일:
scp seungho@40.82.129.113:/home/seungho/personal/THCA_data_analysis/POD_CLEANUP_AND_PAPER1_LOCK_2026_05_04.md .

# 본문 lock 보고서들 (commit f66a6bb 안에 모두 포함, 별도 scp 필요시):
scp seungho@40.82.129.113:/home/seungho/personal/THCA_data_analysis/project/reports/2026_05_04_paper1_dm1_full_molecular_only_lock.md .
scp seungho@40.82.129.113:/home/seungho/personal/THCA_data_analysis/project/reports/2026_05_04_leftover_untracked_inventory.md .
scp seungho@40.82.129.113:/home/seungho/personal/THCA_data_analysis/project/reports/2026_05_04_paper1_dm1_full_conflict_audit.md .
scp seungho@40.82.129.113:/home/seungho/personal/THCA_data_analysis/SITUATION_BACKGROUND_ACTIVITY_2026_05_04.md .
```

---

## 7. Cross-reference index

- `project/reports/2026_05_04_paper1_dm1_full_molecular_only_lock.md` — Paper 1 manuscript-safe envelope (this session, authority)
- `project/reports/2026_05_04_leftover_untracked_inventory.md` — A/B/C/D classification
- `project/reports/2026_05_04_paper1_dm1_full_conflict_audit.md` — parallel agent's claim-by-claim verification
- `project/reports/2026_05_04_image_dm1_final_nogo_decision.md` — closure battery decision authority
- `project/reports/2026_05_04_paper2_post_image_dm1_nogo_status.md` — Paper 2 scope cleanup
- `project/reports/2026_05_04_marathon_state_after_image_dm1_nogo.md` — marathon snapshot
- `CLOSURE_BATTERY_2026_05_04.md` — top-level closure summary
- `pathology_dm1_phaseA_cpu_verdict_2026_05_04.md` — Phase A NO-GO baseline
- `pathology_dm1_closure_battery_2026_05_04.md` — full forensic battery
- `STATUS_PAPER1_2026_05_02.md`, `STATUS_PAPER2_2026_05_02.md`, `STATUS_PAPER3_2026_05_02.md`
- Memory: `v17_marathon_mode_post_pillar1`, `v17_sprint_vs_marathon_violation`, `v18_paper2_HT_isolated`, `v19_paper3_GD_gating`, `v17_npj_ship_status`, `v17_dark_matter_pivot_2026_04_29`

---

*Generated 2026-05-04 by Claude (Opus 4.7) under marathon-mode discipline. No new analysis, no manuscript prose modification, no RunPod/GPU dispatch, no new download. RunPod Pod C/D observed read-only and left to self-terminate.*
