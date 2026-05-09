# Methods paper tables - cross-cohort audit status

Status: 2026-05-09 scaffold. THCA worked example is complete. LUAD and BRCA are staged but blocked because per-slide UNI features, cohort manifests, and audit metadata are not present locally.

## Table 1. Cohort execution status

| Cohort | Proposed pathology task | Features | Manifest | Metadata | OOF predictions | Audit status |
|---|---|---:|---:|---:|---:|---|
| TCGA-THCA | DM1 vs DM2 image pilot | present | present | present | present | complete, 12 tests + 4 corrective baselines |
| TCGA-LUAD | task to lock after metadata assembly (candidate: driver or subtype task) | missing | missing | missing | missing | scaffold only |
| TCGA-BRCA | task to lock after metadata assembly (candidate: PAM50/ER-status task) | missing | missing | missing | missing | scaffold only |

Status files:

| Cohort | Status JSON | Status Markdown |
|---|---|---|
| TCGA-LUAD | `analysis_supp/audit_multicohort_generalization/TCGA_LUAD_input_status.json` | `analysis_supp/audit_multicohort_generalization/TCGA_LUAD_input_status.md` |
| TCGA-BRCA | `analysis_supp/audit_multicohort_generalization/TCGA_BRCA_input_status.json` | `analysis_supp/audit_multicohort_generalization/TCGA_BRCA_input_status.md` |

## Table 2. THCA worked-example audit outputs

| Audit | Metric | Value | Evidence |
|---|---|---:|---|
| C1 | RAS_like intersection FVPTC | 16/16 | `AUDIT_SUMMARY.json` |
| C2 | RAS-like separation gap | 0.017 | `c2_prob_distribution.tsv` |
| C3 | RAS-like permutation p, two-sided | 0.0030 | `c3_permutation_test.json` |
| C4 | Histology-only LR AUC | 0.680 | `c4_shortcut_baselines.json` |
| C5 | RAS-like test folds with zero positives | 3/5 | `c5_RAS_like_fold_composition.tsv` |
| E1 | Female AUC | 0.833 | `PHASE2_AUDIT_SUMMARY.json` |
| E1 | Male AUC | 0.350 | `PHASE2_AUDIT_SUMMARY.json` |
| E2 | RAS-like EM cases | 10/16, all DM2 | `e2_tss_among_RASlike.tsv` |
| E4 | Off-diagonal cPTC/RAS-like + FVPTC/BRAF-like cases | 0 | `PHASE2_AUDIT_SUMMARY.json` |
| E5 | Clinical-only LR AUC | 0.768 | `PHASE2_AUDIT_SUMMARY.json` |
| E6 | Fixed-split init RAS-like AUC range | 0.615-1.000 | `e6_init_variance.tsv` |
| E7 | Label-shuffle overall AUC range | 0.534-0.648 | `e7_label_shuffle_null.tsv` |
| E8 | Pooled LOTO overall AUC | 0.602 | `PHASE3_AUDIT_SUMMARY.json` |
| E8 | Pooled LOTO RAS-like AUC | 0.308 | `PHASE3_AUDIT_SUMMARY.json` |
| S1 | TSS-ComBat overall AUC range | 0.626-0.707 | `phase4/s1_tss_combat_per_init.tsv` |
| S2 | TSS-balanced split overall AUC range | 0.731-0.844 | `phase4/s2_tss_balanced_split.tsv` |
| S3 | Clinical-only vs multimodal LR | 0.768 vs 0.756 | `phase4/s3_multimodal_LR.tsv` |
| S5 | EM-out RAS-like AUC range | 0.667-0.778 | `phase4/s5_em_out.tsv` |

## Table 2b. UNI-final audit

This table audits both the available UNI-final OOF prediction file and the GPU-run UNI leave-one-TSS-out retraining audit.

| Metric | Value | Evidence |
|---|---:|---|
| UNI-final n slides | 54 | `analysis_supp/audit_uni_predictions/UNI_PREDICTION_ONLY_AUDIT_SUMMARY.json` |
| UNI-final overall AUC | 0.874 | `analysis_supp/bootstrap_auc_uni.json`; `analysis_supp/audit_uni_predictions/UNI_PREDICTION_ONLY_AUDIT_SUMMARY.json` |
| UNI-final RAS-like AUC | 0.795 | `analysis_supp/audit_uni_predictions/UNI_PREDICTION_ONLY_AUDIT_SUMMARY.json` |
| UNI-final BRAF-like AUC | 0.891 | `analysis_supp/audit_uni_predictions/UNI_PREDICTION_ONLY_AUDIT_SUMMARY.json` |
| UNI-final Female AUC | 0.893 | `analysis_supp/audit_uni_predictions/UNI_PREDICTION_ONLY_AUDIT_SUMMARY.json` |
| UNI-final Male AUC | 0.812 | `analysis_supp/audit_uni_predictions/UNI_PREDICTION_ONLY_AUDIT_SUMMARY.json` |
| UNI-final clinical-only OOF AUC | 0.735 | `analysis_supp/audit_uni_predictions/UNI_PREDICTION_ONLY_AUDIT_SUMMARY.json` |
| TSS-only within-RAS-like AUC | 1.000 | `analysis_supp/audit_uni_predictions/UNI_PREDICTION_ONLY_AUDIT_SUMMARY.json` |
| UNI-final LOTO overall AUC | 0.852 | `analysis_supp/audit_uni_loto/UNI_LOTO_SUMMARY.json` |
| UNI-final LOTO RAS-like AUC | 0.718 | `analysis_supp/audit_uni_loto/UNI_LOTO_SUMMARY.json` |
| UNI-final LOTO BRAF-like AUC | 0.855 | `analysis_supp/audit_uni_loto/UNI_LOTO_SUMMARY.json` |
| UNI-final LOTO Female AUC | 0.862 | `analysis_supp/audit_uni_loto/UNI_LOTO_SUMMARY.json` |
| UNI-final LOTO Male AUC | 0.875 | `analysis_supp/audit_uni_loto/UNI_LOTO_SUMMARY.json` |

Interpretation: the ViT-L RAS-like 1.000 headline was a TSS artefact, but the UNI-final model is not dead. UNI retains strong center-holdout overall performance (LOTO AUC 0.852), while the RAS-like subgroup falls to an honest underpowered estimate (0.718 rather than 1.000).

## Table 3. Cross-cohort placeholder matrix

| Test | TCGA-THCA | TCGA-LUAD | TCGA-BRCA |
|---|---:|---:|---:|
| C1 subgroup overlap | complete | pending inputs | pending inputs |
| C2 raw subgroup distribution | complete | pending inputs | pending inputs |
| C3 permutation null | complete | pending inputs | pending inputs |
| C4 trivial-feature baseline | complete | pending inputs | pending inputs |
| C5 fold composition | complete | pending inputs | pending inputs |
| E1 sex-stratified AUC | complete | pending inputs | pending inputs |
| E2 TSS x subtype crosstab | complete | pending inputs | pending inputs |
| E3 anchor-case profile | complete | pending inputs | pending inputs |
| E4 histology x subtype gap | complete | pending inputs | pending inputs |
| E5 clinical-covariate baseline | complete | pending inputs | pending inputs |
| E6 init-only variance | complete | pending inputs | pending inputs |
| E7 label-shuffle null | complete | pending inputs | pending inputs |
| E8 leave-one-TSS-out | complete | pending inputs | pending inputs |
| S1 TSS-ComBat | complete | pending inputs | pending inputs |
| S2 TSS-balanced split | complete | pending inputs | pending inputs |
| S3 multimodal LR | complete | pending inputs | pending inputs |
| S5 dominant-site-out sensitivity | complete | pending inputs | pending inputs |

## Table 4. Expected input schema for LUAD/BRCA

| File | Required columns | Notes |
|---|---|---|
| `slide_manifest.tsv` | `file_id`, `submitter_id`, `label`, `tss` | One row per WSI with one binary task label. |
| `features/*.pt` | one tensor file per `file_id` | Expected shape: tiles x feature_dim, e.g. 200 x 1024 for UNI. |
| cohort metadata TSV | `submitter_id`, `sex`, `histology`, `subtype`, `tss` | Extra columns can include age, stage, mutation, PAM50, ER/PR/HER2. |
| `clam_per_slide_predictions.tsv` | `slide`, `fold`, `label`, `prob_pos`, `n_tiles` | Created after first CLAM OOF run. |

## Commands

```bash
python3 project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/scripts/audit_tcga_luad.py
python3 project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/scripts/audit_tcga_brca.py
```

Both commands currently write input-status reports and stop at blockers. Once features and manifests are present, the same wrapper names should become the cohort-specific 12-test entry points.
