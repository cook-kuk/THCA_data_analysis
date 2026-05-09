---
title: "Supplementary Text ST1 - Image-DM1 pilot TSS-confound disclosure"
date: 2026-05-09
status: v1 supplementary case study; Paper 1 absorption path; not main-text prose
source_audit: project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/analysis_supp/audit_ras_auc100/
---

# Supplementary Text ST1. Image-DM1 pilot - TSS-confound disclosure

## Purpose and scope

This supplementary case study documents the post-hoc audit of the Paper 2 image-DM1 pilot and defines how the result should be used in Paper 1. The audit is included here as a transparent negative-control example: a visually plausible TCGA pathology signal can be inflated by tissue-source-site (TSS) structure when sample size is small and subgroup positives are sparse. This section does not support a main-text claim that hematoxylin-and-eosin whole-slide images independently predict DM1. It supports the opposite boundary: TCGA-only image-DM1 remains a pilot result that requires K2 H&E or a multi-site FFPE validation cohort before it can be used as a primary claim.

## Data

The audited model used TCGA-THCA whole-slide image features extracted with a foundation pathology encoder and trained with CLAM attention multiple-instance learning. The primary audited out-of-fold prediction table contains 59 slides with DM1/DM2 labels, per-slide probabilities, and fold assignments. The original report highlighted a RAS-like/FVPTC subgroup AUC of 1.000 in 16 slides, including 3 DM1-positive and 13 DM2-negative cases. Metadata were joined from `project/metadata/sample_master_v3.tsv` to recover histology subtype, molecular subtype, sex, age, stage, and TCGA TSS code. Source files are:

| Item | Path |
|---|---|
| Original OOF predictions | `project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/phase2_tcga_clam/clam_per_slide_predictions.tsv` |
| Slide manifest | `project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/phase2_tcga_clam/slide_manifest.tsv` |
| Audit root | `project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/analysis_supp/audit_ras_auc100/` |
| Dossier | `project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/analysis_supp/audit_ras_auc100/AUDIT_DOSSIER.html` |
| Hub copy | `project/papers_hub_2026_05_04/paper2_image_dm1_audit_dossier.html` |

## Audit design

The audit used four layers. Phase 1 ran no-retraining checks: subgroup overlap, raw prediction distribution, permutation null, histology-only shortcut baseline, and fold composition. Phase 2 quantified confounders: sex-stratified performance, TSS distribution, case-level profiles, histology by molecular-subtype gaps, and a clinical-only logistic-regression baseline. Phase 3 performed decisive retraining controls: fixed-split initialization variance, label-shuffle null, and leave-one-TSS-out (LOTO) cross-validation. Phase 4 tested corrective baselines: TSS-ComBat feature adjustment, TSS-balanced splitting, multimodal logistic regression, and an EM-out sensitivity analysis.

The central criterion was whether the image signal survived TSS holdout. Random K-fold and stratified folds can place the same tissue-source-site centers into training and test sets. LOTO removes that shortcut by withholding all slides from one TSS at a time and pooling the held-out predictions.

## Results

The RAS-like and FVPTC subgroup labels were identical in the 59-slide subset: RAS_like intersection FVPTC = 16/16, with no off-diagonal cPTC/RAS-like or FVPTC/BRAF-like cases available for disentanglement. The apparent perfect ranking in the 16-slide subgroup was narrow: the lowest DM1 probability was 0.549 and the highest DM2 probability was 0.532, giving a separation gap of 0.017. The permutation p value for the RAS-like AUC was 0.0030, but the bootstrap confidence interval of [1.000, 1.000] reflected deterministic resampling of the same perfect ordering, not center-generalized performance.

The key audit findings were:

| Test | Value | Interpretation | Source |
|---|---:|---|---|
| Original RAS-like/FVPTC AUC | 1.000 | Pre-audit headline in n=16 | `AUDIT_SUMMARY.json` |
| RAS-like/FVPTC overlap | 16/16 | One subgroup labeled two ways | `AUDIT_SUMMARY.json` |
| RAS-like separation gap | 0.017 | One near-boundary rank supports the perfect AUC | `c2_prob_distribution.tsv` |
| Folds with zero RAS-like positives | 3/5 | Original K-fold split is poorly balanced for the subgroup | `c5_RAS_like_fold_composition.tsv` |
| Male AUC | 0.350 | Sex-stratified failure in n=13 males | `PHASE2_AUDIT_SUMMARY.json` |
| Clinical-only LR AUC | 0.768 | Histology + sex + TSS beats image-inclusive multimodal LR | `phase4/s3_multimodal_LR.tsv` |
| Multimodal LR AUC | 0.756 | CLAM probability adds no net gain over metadata | `phase4/s3_multimodal_LR.tsv` |
| TSS-only LR within RAS-like | 1.000 | Perfect ranking can be reproduced without image features | `tss_only_within_raslike_lr.tsv` |
| Pooled LOTO overall AUC | 0.602 | Center holdout weakens overall performance | `PHASE3_AUDIT_SUMMARY.json` |
| Pooled LOTO RAS-like AUC | 0.308 | Center holdout reverses the subgroup signal | `PHASE3_AUDIT_SUMMARY.json` |
| UNI-final OOF AUC | 0.874 | Foundation encoder improves overall ranking in n=54 | `analysis_supp/bootstrap_auc_uni.json` |
| UNI-final LOTO overall AUC | 0.852 | UNI retains center-holdout overall performance | `analysis_supp/audit_uni_loto/UNI_LOTO_SUMMARY.json` |
| UNI-final LOTO RAS-like AUC | 0.718 | Honest subgroup estimate is no longer perfect and remains underpowered | `analysis_supp/audit_uni_loto/UNI_LOTO_SUMMARY.json` |

The TSS distribution explained why the subgroup looked perfect under ordinary splits. All 3 RAS-like DM1 cases came from DJ or FK, while all 10 RAS-like cases from EM were DM2. A logistic regression using only TSS dummy variables, with no image features, achieved within-RAS-like AUC=1.000. Under pooled LOTO for the audited ViT-L feature table, the overall AUC dropped to 0.602 and the RAS-like AUC dropped to 0.308, which is worse than random. A later GPU-run LOTO audit of the UNI-final feature table refined rather than erased this conclusion: UNI retained strong overall center-holdout performance (LOTO AUC 0.852), but the RAS-like subgroup estimate was 0.718 rather than 1.000 and still rested on only 3 positive RAS-like cases.

## Corrective baselines

No corrective baseline recovered the original 1.000 in a center-generalized manner. TSS-ComBat reduced the apparent batch component and partially repaired the male failure, but the overall median AUC remained in the 0.63-0.71 range across tested initialization seeds. TSS-balanced splitting preserved high RAS-like AUC when training and test folds shared centers, which is useful as a matched-center ceiling but does not answer the center-holdout question. Multimodal logistic regression quantified the image channel directly: clinical-only AUC was 0.768, CLAM-only AUC in the corrective baseline was 0.691, and multimodal AUC was 0.756. Removing the 10 EM RAS-like DM2 anchor slides reduced the RAS-like median to approximately 0.72-0.78 depending on initialization.

## Conclusion for Paper 1

For Paper 1, the image-DM1 result should be treated as a cautionary supplementary pilot, not as a supporting pillar. The appropriate use is a TSS-confound disclosure showing that TCGA single-cohort pathology AI can appear stronger than it is when subgroup labels, histology, sex, and acquisition center are entangled. The UNI-final LOTO result also preserves a positive future path: the foundation-model image channel may carry real overall DM1 signal, but the RAS-like/FVPTC subgroup claim should be reported as underpowered and center-sensitive until K2 H&E or a multi-site FFPE cohort validates it externally. The finding does not weaken the Paper 1 molecular DM1 axis, which is based on RNA, methylation, fusion, survival, and cross-cohort evidence; it only blocks a premature image-derived extension of that axis.

## Recommendation

The Paper 2 image-DM1 line should proceed only after external validation. The minimum defensible design is K2 H&E or a multi-site FFPE cohort with a locked TSS-balanced split, clinical-only baseline, multimodal image-gain estimate, sex-stratified AUC, and leave-one-site-out or leave-one-center-out validation. Until those data are available, report the TCGA image model as an audited pilot: UNI shows promising center-holdout overall performance, but the original perfect RAS-like/FVPTC subgroup headline is not defensible.
