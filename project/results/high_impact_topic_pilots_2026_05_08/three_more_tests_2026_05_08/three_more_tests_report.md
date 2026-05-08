# Three more high-impact tests (2026-05-08)

## Test 1 — T00 CNV robustness

Question: does the cBioPortal ARMDRIVER signature survive simple morphology/purity/clinical stratification, and does it add signal over covariates?

| subset | n DM2 | n DM1 | DM2 rate | DM1 rate | OR | p |
|---|---:|---:|---:|---:|---:|---:|
| full_Class6_DM2_vs_DM1 | 48 | 76 | 0.354 | 0.013 | 41.129 | 1.607e-07 |
| follicular_pct<median(70.000) | 7 | 42 | 0.571 | 0.000 | NA | 1.652e-04 |
| stage=Stage I | 23 | 50 | 0.348 | 0.020 | 26.133 | 2.609e-04 |
| purity>=median(0.685) | 30 | 21 | 0.400 | 0.000 | NA | 6.133e-04 |
| histology_subtype=cPTC | 18 | 55 | 0.278 | 0.018 | 20.769 | 0.003 |
| BRAFV600E_RAS=Ras-like | 41 | 13 | 0.390 | 0.000 | NA | 0.006 |
| purity<median(0.685) | 11 | 40 | 0.364 | 0.025 | 22.286 | 0.006 |
| histology_subtype=FVPTC | 27 | 14 | 0.407 | 0.000 | NA | 0.007 |
| stage=Stage III | 7 | 9 | 0.571 | 0.000 | NA | 0.019 |
| follicular_pct>=median(70.000) | 36 | 20 | 0.333 | 0.050 | 9.500 | 0.020 |
| ARM_SCNA_CLUSTER=SomeSCNA | 10 | 12 | 0.400 | 0.083 | 7.333 | 0.135 |
| ARM_SCNA_CLUSTER=Quiet | 24 | 55 | 0.000 | 0.000 | NA | 1.000 |

Cross-validated covariates-only AUC: 0.935
Cross-validated covariates+CNV-signature AUC: 0.947
Mean delta AUC: 0.012

## Test 2 — Evidence-grounded retrieval

Leave-one-slide-out KNN retrieval; no target-specific regression. `fusion_a06` uses 60% morphology and 40% stage text.

| target | best mode | best rho | best AUROC q75 | morph rho | text rho | fusion_a06 rho |
|---|---|---:|---:|---:|---:|---:|
| CAF_ECM_score | fusion_a08 | 0.090 | 0.563 | 0.072 | -0.115 | 0.088 |
| DM1_like_score | morphology_only | 0.173 | 0.670 | 0.173 | -0.038 | 0.163 |
| Hypoxia_score | morphology_only | 0.108 | 0.547 | 0.108 | -0.005 | 0.059 |
| RAI_8_score | morphology_only | 0.173 | 0.535 | 0.173 | -0.038 | 0.163 |
| TDS_like_score | morphology_only | 0.142 | 0.528 | 0.142 | -0.049 | 0.121 |
| TROP2 | fusion_a06 | 0.532 | 0.825 | 0.021 | 0.464 | 0.532 |

## Test 3 — Haiku-style counterfactual stage text

Holding morphology fixed, change only the stage text in fusion retrieval and summarize predicted target shifts.

| target | PTC-PT | LPTC-PT | ATC-PT | LPTC-PTC | ATC-PTC |
|---|---:|---:|---:|---:|---:|
| CAF_ECM_score | 4.412e-04 | -0.077 | 0.084 | -0.078 | 0.083 |
| DM1_like_score | 0.024 | -0.083 | 0.029 | -0.107 | 0.004 |
| Hypoxia_score | -0.021 | 0.021 | -0.005 | 0.042 | 0.016 |
| RAI_8_score | -0.024 | 0.083 | -0.029 | 0.107 | -0.004 |
| TDS_like_score | -3.448e-04 | 0.100 | 0.001 | 0.101 | 0.001 |
| TROP2 | 4.670 | 5.422 | -0.759 | 0.752 | -5.428 |

## Bottom line

- T00 remains the main biology flagship if the signature improves DM2 classification beyond covariates and remains enriched in key strata.
- T12 is strengthened only if retrieval, not just regression, recovers TROP2 with morphology+text.
- Counterfactual shifts should be presented as hypothesis-generating, matching Haiku's own caution.

## Files
- `project/results/high_impact_topic_pilots_2026_05_08/three_more_tests_2026_05_08/test1_cnv_stratified_robustness.tsv`
- `project/results/high_impact_topic_pilots_2026_05_08/three_more_tests_2026_05_08/test1_cnv_covariate_auc.tsv`
- `project/results/high_impact_topic_pilots_2026_05_08/three_more_tests_2026_05_08/test2_knn_retrieval_metrics.tsv`
- `project/results/high_impact_topic_pilots_2026_05_08/three_more_tests_2026_05_08/test2_knn_retrieval_predictions.tsv.gz`
- `project/results/high_impact_topic_pilots_2026_05_08/three_more_tests_2026_05_08/test3_counterfactual_stage_summary.tsv`
- `project/results/high_impact_topic_pilots_2026_05_08/three_more_tests_2026_05_08/test3_counterfactual_stage_predictions.tsv.gz`
