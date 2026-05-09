# Methods paper outline — TCGA-thyroid pathology AI confound audit framework

**Working title:** "Site-of-origin confounds in TCGA pathology AI for thyroid molecular subtyping: a 12-test audit framework with three corrective baselines"

**Target venue:** Nature Communications (Methods) or Nat Comm Medicine. DIAL paper sequel (`v5_broken/v5_dial_paper.tex`) sits well as the methodological lineage — DIAL audits batch-effect leakage in *bulk RNA* classifiers; this audits the same pathology in *tile-level WSI* classifiers.

**Status:** Drafted 2026-05-09 from p2 audit data. Marathon mode = scaffolding/infra allowed. Full draft after 2026-06-13 marathon ends.

---

## 1. Pitch (one paragraph)

Pathology foundation models (UNI, CONCH, GigaPath) report increasingly impressive per-slide AUCs on TCGA cancer subtype tasks, but the field's standard validation recipe — random KFold cross-validation on a single cohort — is structurally insensitive to tissue-source-site (TSS) batch effects baked into TCGA's multi-center sampling. We use the published Paper 2 image-DM1 result on TCGA-THCA (n=59, foundation-model UNI + CLAM, headline subgroup AUC = 1.000 in RAS-like cases) as a worked example to (i) introduce a 12-test audit framework, (ii) show that the apparent 1.000 collapses to 0.308 under leave-one-TSS-out cross-validation, (iii) demonstrate that a 3-feature logistic regression on (histology + sex + TSS-dummy) already attains 0.768 — exceeding the foundation-model AUC, and (iv) propose three corrective baselines (TSS-ComBat at the tile level, TSS-balanced StratifiedKFold, and a multimodal clinical+image LR) that re-establish honest AUC ranges. The framework runs in <30 min on CPU using pre-extracted features and applies to any TCGA pathology task.

## 2. Why this matters now

- TCGA pathology AI has accelerated since 2023 (Howard 2021 *Nat Comm* original site-confound paper; UNI/CONCH 2024; CHIEF 2024; GigaPath 2024 — all published with random KFold on single cohorts).
- Most papers do **not** report leave-one-TSS-out, do **not** disentangle clinical-only baselines, do **not** report sex-stratified AUC.
- The result: a generation of "foundation model > classical CNN" claims that may be 50% TSS-batch in origin.
- DIAL (this lab, 2024 prep) demonstrated the analogous problem in bulk-RNA cancer-subtype classifiers and provided a label-flip diagnostic. **The image-AI version of DIAL has not been published.**

## 3. The 12-test audit framework

| # | Test | Detects | Compute |
|---|---|---|---|
| C1 | Subgroup overlap matrix | redundant subgroup definitions | <1 s |
| C2 | Raw OOF prediction distribution within subgroup | per-sample separation margin | <1 s |
| C3 | Permutation null AUC | model ranking ≠ random baseline | 10 s |
| C4 | Histology-only LR baseline | trivial-feature shortcut | <1 s |
| C5 | Per-fold class composition | unstratified-fold imbalance | <1 s |
| E1 | Sex-stratified AUC | sex-asymmetric failure | <1 s |
| E2 | TSS × subtype crosstab | TSS-confound prior | <1 s |
| E3 | Subgroup case-level profile | identifying anchor cases | <1 s |
| E4 | Histology × subtype gap | disentanglement feasibility | <1 s |
| E5 | Clinical-covariate baseline (LR) | image net gain over clinical | 5 s |
| E6 | Init-only seed variance (split fixed) | stochastic instability | 5 min |
| E7 | Label-shuffle null retraining | training-level memorisation null | 5 min |
| E8 | Leave-one-TSS-out retraining | TSS batch dependence | 8 min |

**The diagnostic verdict** is read off as a 4-cell quadrant (analogous to DIAL):
- (high OOF AUC, high LOTO AUC) = real signal
- (high OOF AUC, low LOTO AUC) = TSS-batch artefact
- (low OOF AUC, low LOTO AUC) = no signal
- (low OOF AUC, high LOTO AUC) = unlikely (suggests training instability)

## 4. Three corrective baselines

| Baseline | What it does | Expected effect |
|---|---|---|
| S1 — TSS-ComBat at the tile level | per-feature mean-centering by TSS group before MIL aggregation | removes additive batch component; reduces apparent AUC to honest ceiling |
| S2 — StratifiedKFold on (label × TSS-bucket) | guarantees TSS distribution is matched between train/val of each fold | preserves TSS-balanced AUC; combined with S1 gives honest range |
| S3 — Multimodal clinical + image LR | quantifies image net gain over (histology + sex + TSS dummies) | reveals when the image channel is redundant with clinical metadata |

## 5. Worked example — Paper 2 / TCGA-THCA / UNI+CLAM

**Headline before audit:** RAS-like AUC = 1.000 (95% boot CI [1.000, 1.000]); FVPTC AUC = 1.000; overall AUC = 0.746 (ViT-L) / 0.874 (UNI).

**After 12-test audit:**
- LOTO RAS-like AUC = 0.308 (was 1.000)
- LOTO overall AUC = 0.602 (was 0.746)
- TSS-only LR within-RAS-like AUC = 1.000 (image not needed)
- Multimodal LR AUC = 0.756 < clinical-only 0.768
- Sex Male AUC = 0.35; ComBat correction lifts to 0.75 (median over inits)

**Honest ceiling:** TCGA-only image-DM1 AUC ≈ 0.65–0.80 with clinical-or-better performance only achievable through multimodal integration; image alone provides 0% net gain on this cohort.

## 6. Generalization claims

The framework is task-agnostic. Apply to any TCGA pathology classifier:
- Need: per-slide UNI/CONCH/GigaPath features (200 tiles × 1024-d sufficient)
- Need: case-level metadata (TSS, sex, histology, molecular subtype)
- Output: 12 audit numbers + 3 correction baselines
- Time: <30 min CPU for n<200, <2 hr GPU for n>1000

## 7. Limitations

- 4-cell verdict requires both OOF AUC and LOTO AUC; LOTO needs ≥10 TSS centers in the cohort.
- ComBat at tile level is a strong simplification; full ComBat (with per-tile covariance) is more thorough but rarely warranted at n<500.
- Doesn't address slide-level scanner artefacts directly (color-norm preprocessing recommended in addition).

## 8. Deliverables (expected)

- Methods paper with 4 figures + 1 supplementary 12-test panel
- GitHub repo: `claude-code/audit-pathology-ai` with reusable scripts:
  - `audit_phase1.py` (C1-C5)
  - `audit_phase2_confounds.py` (E1-E5)
  - `audit_phase3_decisive.py` (E6-E8)
  - `audit_phase4_solutions.py` (S1-S3)
  - `audit_make_figures.py` + `audit_make_report.py`
- All scripts already exist at `project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/scripts/audit_*.py` — ready to package.

## 9. Co-author plan

- Cook (lead, methods + audit framework)
- Yu (senior, mentorship + clinical interpretation)
- Cho (RNA-seq cohort context, paper 1 link)
- TBD pathology consultant

## 10. Timeline (post-marathon, 2026-06-13 onward)

- Week 1–2: package GitHub repo + dockerize
- Week 3–4: extend audit to 2 more TCGA cohorts (LUAD + BRCA) for generalization
- Week 5–6: write manuscript draft
- Week 7: internal review (Yu)
- Week 8: submission to Nat Comm Methods

## 11. Why this is shippable (defensive scoping)

- All compute already done for thyroid worked example (audit phase 1-4 complete 2026-05-09).
- LUAD + BRCA replication is mechanical (same scripts, swap input).
- DIAL paper provides the methodological precedent — reviewers familiar with the framing.
- Paper 2 itself doesn't have to be retracted or revised before this is submitted; this paper *frames* the Paper 2 result as a worked example of the field-wide problem.
