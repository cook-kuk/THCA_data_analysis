# CROSS-Neo v6 wetlab algorithm comparison

Generated: 2026-05-10T22:36:45

Result source: `smoke_confirmatory_threshold_v6`

## Bottom line

- Best comparator in this result sheet: `BigMHC_EL` (AUPRC 0.866, AUROC 0.664, top5 precision 1.000).
- Countable candidates: 24/24
- PASS calls: 17

## Algorithm metrics

| algorithm                           | score_column                        |   n |   hits |   hit_rate |    AUROC |    AUPRC |   top1_hits |   top1_precision |   top2_hits |   top2_precision |   top3_hits |   top3_precision |   top5_hits |   top5_precision |   top10_hits |   top10_precision |   top12_hits |   top12_precision |   top24_hits |   top24_precision |
|:------------------------------------|:------------------------------------|----:|-------:|-----------:|---------:|---------:|------------:|-----------------:|------------:|-----------------:|------------:|-----------------:|------------:|-----------------:|-------------:|------------------:|-------------:|------------------:|-------------:|------------------:|
| BigMHC_EL                           | bigmhc_el_score                     |  24 |     17 |   0.708333 | 0.663866 | 0.865599 |           1 |                1 |           2 |              1   |           3 |         1        |           5 |              1   |            8 |               0.8 |            9 |          0.75     |           17 |          0.708333 |
| CROSS_Neo_integrated_immunogenicity | immunogenicity_discovery_score      |  24 |     17 |   0.708333 | 0.596639 | 0.790492 |           1 |                1 |           2 |              1   |           2 |         0.666667 |           4 |              0.8 |            8 |               0.8 |            9 |          0.75     |           17 |          0.708333 |
| CROSS_Neo_claim_safe_integrated     | immunogenicity_claim_safe_score     |  24 |     17 |   0.708333 | 0.588235 | 0.788084 |           1 |                1 |           2 |              1   |           2 |         0.666667 |           4 |              0.8 |            8 |               0.8 |            9 |          0.75     |           17 |          0.708333 |
| CROSS_Neo_stress_guarded            | stress_guarded_discovery_score      |  24 |     17 |   0.708333 | 0.521008 | 0.776958 |           1 |                1 |           2 |              1   |           3 |         1        |           4 |              0.8 |            7 |               0.7 |            9 |          0.75     |           17 |          0.708333 |
| CROSS_Neo_BMA_v2                    | bma_v2_discovery_score              |  24 |     17 |   0.708333 | 0.453782 | 0.766804 |           1 |                1 |           2 |              1   |           3 |         1        |           4 |              0.8 |            6 |               0.6 |            7 |          0.583333 |           17 |          0.708333 |
| CROSS_Neo_finetuned_priority        | finetuned_experiment_priority_score |  24 |     17 |   0.708333 | 0.420168 | 0.681142 |           0 |                0 |           1 |              0.5 |           2 |         0.666667 |           3 |              0.6 |            7 |               0.7 |            8 |          0.666667 |           17 |          0.708333 |
| BigMHC_IM                           | bigmhc_im_score                     |  24 |     17 |   0.708333 | 0.453782 | 0.672325 |           0 |                0 |           1 |              0.5 |           1 |         0.333333 |           3 |              0.6 |            7 |               0.7 |            8 |          0.666667 |           17 |          0.708333 |

## Hit/fail decomposition

| call                             |   n |
|:---------------------------------|----:|
| hit_specificity_control          |   4 |
| hit_clean_antigen                |   3 |
| fail_presentation_or_specificity |   3 |
| hit_positive_qc                  |   3 |
| hit_model_boundary_resolution    |   3 |
| hit_tcr_md_mechanism             |   2 |
| fail_label_rescue                |   2 |
| hit_label_rescue                 |   2 |
| fail_assay_qc                    |   1 |
| fail_model_boundary              |   1 |

## Endpoint unlocks

| plate_v4_arm          |   success_count |   threshold_successes |   confirmatory_threshold_successes | exploratory_unlock_status   | confirmatory_unlock_status   | claim_layer_after_unlock                      |
|:----------------------|----------------:|----------------------:|-----------------------------------:|:----------------------------|:-----------------------------|:----------------------------------------------|
| B_mechanism_TCR_MD    |               2 |                     1 |                                  2 | UNLOCKED                    | UNLOCKED                     | orthogonal TCR/MD mechanistic support         |
| A_clean_discovery     |               3 |                     2 |                                  3 | UNLOCKED                    | UNLOCKED                     | wetlab-supported clean antigen discovery lane |
| C_label_rescue        |               2 |                     1 |                                  2 | UNLOCKED                    | UNLOCKED                     | label-noise or assay-mismatch evidence        |
| E_positive_QC_control |               3 |                     3 |                                  3 | UNLOCKED                    | UNLOCKED                     | assay QC only                                 |
| D_specificity_moat    |               4 |                     3 |                                  4 | UNLOCKED                    | UNLOCKED                     | specificity-aware control moat                |
| F_model_boundary      |               3 |                     2 |                                  3 | UNLOCKED                    | UNLOCKED                     | model-boundary clarification                  |

## Claim boundary

The default run uses a smoke-confirmatory result sheet, not real prospective wetlab data. When actual assay calls are entered, rerun this same script with `--candidate-calls <interpreter output>/candidate_calls_v6.tsv --result-source actual_wetlab_v6`.
