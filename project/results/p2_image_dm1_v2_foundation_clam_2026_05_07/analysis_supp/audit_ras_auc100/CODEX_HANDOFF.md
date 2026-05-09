# Codex handoff — Paper 2 image-DM1 audit (2026-05-09)

You are continuing a Paper 2 image-pathology audit started by Claude Code on 2026-05-09. This document is self-contained — read it top to bottom, then pick up from §6 "Open tasks".

## 1. TL;DR — what this audit found

A user noticed that Paper 2's TCGA-THCA image-DM1 model (UNI foundation model + CLAM, n=59 slides) reported subgroup AUC = 1.000 in RAS-like / FVPTC cases (n=16, n_pos=3) and asked whether it was a data-split artefact. A 4-phase audit (12 tests + 4 corrective baselines, all CPU, all <30 min total) established:

1. **The 1.000 is a TSS batch-effect artefact.** Pooled leave-one-TSS-out (LOTO) cross-validation: overall AUC drops to 0.602; within-RAS-like AUC drops to **0.308** (worse than random). A logistic regression on TSS dummies alone — with zero image features — already achieves within-RAS-like AUC = 1.000 because the 16 RAS-like slides split deterministically by TSS center (DJ/FK = 3 DM1, EM = 10 DM2).

2. **The image channel adds zero net gain over clinical metadata on TCGA.** Multimodal LR (CLAM_prob + histology + sex + TSS dummies) gives OOF AUC = 0.756, which is *below* the clinical-only baseline of 0.768.

3. **Sex stratification fails.** Female AUC 0.83 [0.70, 0.94] vs Male AUC **0.35** [0.00, 0.73] — model is worse than random in n=13 males. ComBat partially fixes this (Male median 0.75 across 4 init seeds).

4. **The 4 corrective baselines (S1 ComBat, S2 TSS-balanced split, S3 multimodal, S5 EM-out) cannot recover the 1.000.** TCGA-only honest ceiling: overall AUC 0.65–0.80, RAS-like AUC 0.70–0.90 (highly seed-dependent with n_pos=3).

5. **The fundamental fix is more data, not more analysis.** Paper 2 requires K2 H&E (or FFPE multi-site) external validation with TSS-balanced design before any image-DM1 main-text claim is defensible.

## 2. File map — where everything lives

Working directory: `/home/seungho/personal/THCA_data_analysis/project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/`

### Audit scripts (all CPU, all reusable)
- `scripts/audit_clam_train_lib.py` — shared CLAM training library
- `scripts/audit_ras_auc100.py` — Phase 1: C1-C5 (subgroup overlap, perm null, shortcut probe, fold composition)
- `scripts/audit_stratified_retrain.py` — Phase 1+: Strategy A/B/C (multi-seed KFold, StratifiedKFold(3), LOO on 16 RAS-like)
- `scripts/audit_phase2_confounds.py` — Phase 2: E1-E5 (sex, TSS, profile, gap, clinical baseline)
- `scripts/audit_phase3_decisive.py` — Phase 3: E6-E8 (init variance, label-shuffle null, leave-one-TSS-out)
- `scripts/audit_phase4_solutions.py` — Phase 4: S1-S5 (TSS-ComBat, TSS-balanced split, multimodal LR, EM-out subsample)
- `scripts/audit_make_figures.py` + `audit_make_figures_phase23.py` — figures A-G
- `scripts/audit_make_report_v2.py` — final report builder

### Audit deliverables (under `analysis_supp/audit_ras_auc100/`)
- `AUDIT_REPORT_v2.md` (247 lines) — final synthesis report
- `METHODS_PAPER_OUTLINE.md` — Path 2 (Methods paper sequel to DIAL)
- `K2_HE_PUSH_EMAIL_DRAFT.md` — Path 3 push email (Korean + English)
- `figures/fig_A..G_*.{png,pdf}` — 7 figures (raw 16-slide distribution, perm null, fold imbalance, retrain summary, sex+TSS confound, clinical baseline, phase 3 decisive)
- `*.tsv`, `*.json` — 17 raw-evidence files
- `merged_preds_subgroups.tsv` — convenience: 59-slide OOF preds + subgroup labels in one place

### Paper 2 dossier (modified by audit)
- `REPORT_YU_REVIEW_PAPER2_2026_05_08.html` — header now contains red audit caveat box
- `manuscript_draft/02_methods_results.md` — methods §2 has new "Post-hoc audit" paragraph

### Inputs (don't modify)
- `phase2_tcga_clam/features/*.pt` — 88 UNI features (200 tiles × 1024-d each)
- `phase2_tcga_clam/slide_manifest.tsv` — 90 slides with file_id, submitter_id, dm
- `phase2_tcga_clam/clam_per_slide_predictions.tsv` — original 5-fold OOF preds (overall AUC 0.746)
- `/home/seungho/personal/THCA_data_analysis/project/metadata/sample_master_v3.tsv` — TCGA-THCA case metadata (histology_subtype, molecular_subtype, sex, age, ajcc_stage_group)

## 3. How to re-run anything

```bash
cd /home/seungho/personal/THCA_data_analysis
# Phase 1 (no retraining, 1 min)
python3 project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/scripts/audit_ras_auc100.py

# Phase 1+ (multi-seed retrain, ~15 min CPU)
python3 project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/scripts/audit_stratified_retrain.py

# Phase 2 (no retraining, <10 s)
python3 project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/scripts/audit_phase2_confounds.py

# Phase 3 (E6 + E7 + E8 retrain, ~17 min CPU)
python3 project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/scripts/audit_phase3_decisive.py

# Phase 4 (S1 ComBat + S2 TSS-balanced + S3 LR + S5 EM-out, ~17 min CPU)
python3 project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/scripts/audit_phase4_solutions.py

# Build figures + final report
python3 project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/scripts/audit_make_figures.py
python3 project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/scripts/audit_make_figures_phase23.py
python3 project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/scripts/audit_make_report_v2.py
```

## 4. Project constraints (read before acting)

- **Marathon mode (active 2026-04-30 → 2026-06-13).** Voice-protected sections (Hook · Aim · Disc 3.1 · Limitations · Cover Para 1 · Q9) are author-keyboard only — Codex must NOT generate prose for those. Scaffolding/infra (this audit, methods paper outline, scripts, figures) is allowed. New analyses must be paper-blocking only.
- **Disk layout:** `project/results/...` and `project/data/...` are bind-mounted to `/data/thca/...` (512G SSD). Don't write large outputs (>100 MB) to other paths.
- **Outreach is author-only.** Email drafts (`K2_HE_PUSH_EMAIL_DRAFT.md`) are for the user to edit and send, never auto-send.
- **TCGA-THCA WSI**: only 88 cases have UNI features. Don't fabricate or extrapolate.
- **Voice-protected**: don't write prose for Hook / Aim / Disc 3.1 / Limitations / Cover Para 1 / Q9. Methods, supplementary, captions, infrastructure are fine.

## 5. Three forward paths (already drafted)

The user has been presented with three high-impact forward paths:

| Path | Goal | Status | Blocker |
|---|---|---|---|
| **Path 1** — Multi-cohort multimodal nomogram | NC / Nat Cancer ship | requires K2 H&E + TCGA full WSI | K2 timing unknown |
| **Path 2** — Methods paper (DIAL sequel, this audit packaged) | NC Methods | outline drafted | needs LUAD/BRCA replication post-marathon |
| **Path 3** — Paper 1 absorption | safe NC | requires Paper 2 demoted to Paper 1 supp | user/Yu agreement needed |

User said "다 해야지 다" — execute all three. Path 2 outline is `METHODS_PAPER_OUTLINE.md`. Path 3 has been started by adding audit caveats to Paper 2 dossier; full demotion to Paper 1 supp is the natural next step.

## 6. Open tasks for Codex (concrete, ranked by ROI)

### High-ROI / can ship before marathon ends (2026-06-13)

**T1. Generalize the audit framework to LUAD + BRCA TCGA cohorts.**
- Why: Methods paper (Path 2) needs ≥2 cohorts beyond thyroid for cross-cohort generalization claim.
- Inputs needed: per-slide UNI features for TCGA-LUAD and TCGA-BRCA (foundation-model embedding standard pipeline). Likely need to download via GDC + run UNI inference. ~50 GB raw WSI per cohort, ~1 hr per cohort on RTX A6000.
- Deliverable: `audit_*_LUAD.py`, `audit_*_BRCA.py` (mirror thyroid scripts) + `METHODS_PAPER_TABLES.md` summarizing 12-test outputs across 3 cohorts.
- Scripts in `scripts/audit_*.py` are parameterized — should mostly just need new feature_dir paths.

**T2. Build a single combined HTML dossier (`AUDIT_DOSSIER.html`) for advisor review.**
- Why: user style preference (`feedback_dossier_pattern.md`) — sticky TOC + traceable numbers + co-located honest caveats + decision matrix. User said "대박" 2026-05-08 to that pattern.
- Inputs: 7 figures + AUDIT_REPORT_v2.md + METHODS_PAPER_OUTLINE.md + K2_HE_PUSH_EMAIL_DRAFT.md.
- Layout: §1 verdict box → §2 phase 1 confound landscape → §3 phase 2/3 decisive evidence → §4 phase 4 solutions → §5 forward paths → §6 supplementary.
- Reference style: existing `papers_hub_2026_05_04/*.html` files have established patterns; `feedback_dossier_pattern.md` memory describes the dossier convention.

**T3. Add the audit results to Paper 1 supplementary as a TSS-confound case study.**
- Why: Path 3 (Paper 1 absorption). Currently the audit is in Paper 2 dossier only.
- Action: write a 2-page Paper 1 supp section ("Image-DM1 pilot — TSS confound disclosure") covering data, audit, conclusion, recommendation for K2 validation. Place in `project/manuscript_v8/sup/` or wherever Paper 1 supp materials live (check `manuscript_v8/` directory structure first).
- Constraint: this is supp text, not main-text Hook/Aim/Disc 3.1/Limitations — auto-generation allowed.

### Medium-ROI / needs K2 timing first

**T4. K2 H&E retrieval status check.**
- Action: ask user for the recipient name + email of the K2 cohort PI / SNUH-GMI contact.
- Then: fill in `[받는 분 호칭]` / `[recipient]` placeholders in `K2_HE_PUSH_EMAIL_DRAFT.md` and present the final draft to user for sending. Do NOT send.

**T5. If K2 timing is "available within 3 months":** start Path 1 by drafting the multi-cohort multimodal training pipeline (clinical + image + RNA → DM1 risk score → OS Cox HR).

**T6. If K2 timing is "unavailable / >12 months":** lock in Path 2 + Path 3 as the actual plan; archive Path 1 as future work.

### Low-ROI but useful for completeness

**T7. Replicate the audit on the original ViT-L (n=59, AUC 0.746) and UNI (n=54, AUC 0.874) checkpoints.**
- Currently the audit is on the ViT-L OOF predictions (`clam_per_slide_predictions.tsv`). The UNI-final result (`bootstrap_auc_uni.json`) reports overall 0.874 but I haven't audited it.
- Action: re-run audit with UNI per-slide predictions (file path TBD — search `phase2_tcga_clam_UNI/` if it exists).

**T8. Tile-budget ablation under TSS-LOTO.**
- Already done (50 / 100 / 150 / 200 tiles, in-sample only) in `manuscript_draft/02_methods_results.md`. Re-run as out-of-fold under LOTO to make the tile-budget claim defensible.

## 7. Reproducibility / sanity-check

To verify the audit before adding new work, re-run Phase 1:
```bash
python3 project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/scripts/audit_ras_auc100.py
```
Expected output (key numbers):
- C1 RAS_like ∩ FVPTC = 16/16
- C2 lowest DM1 prob = 0.549, highest DM2 prob = 0.532, gap = 0.017
- C3 RAS-like permutation p-two-sided = 0.0030
- C4 histology-only LR AUC = 0.680
- C5 fold composition: folds 1, 2, 4 have 0 RAS-like positives in test

Phase 3 LOTO key numbers (E8):
- Pooled overall LOTO AUC = 0.602
- Pooled RAS-like LOTO AUC = 0.308
- Per-TSS LOTO: BJ=0.93, DE=1.00, DJ=0.67, EM=0.95, ET=0.50, FK=1.00, others n/a

If any of these don't match within ±0.01, the dataset has changed since this handoff was written.

## 8. Quick-glance summary table (for Codex's first response back to user)

| Audit number | Value | Meaning |
|---|---:|---|
| Original RAS-like AUC | 1.000 | Headline before audit |
| LOTO RAS-like AUC | 0.308 | Worse than random when TSS withheld |
| TSS-only LR within-RAS | 1.000 | Image not needed for "perfect" ranking |
| Multimodal LR AUC | 0.756 | < clinical-only 0.768 (image gain = 0) |
| Male AUC (raw) | 0.35 | Sex-stratified failure |
| Honest TCGA-only ceiling | 0.65–0.80 | What can defensibly be reported |

## 9. Last user input

User asked: "이거 codex에서 이어가게 하고싶어" (= "I want to continue this in Codex").

Most likely next user instruction: "T2 (HTML dossier) 만들어줘" or "T1 (LUAD/BRCA generalization) 시작해" or "T3 (Paper 1 supp) 추가해". Confirm with user which task to pick up before starting.
