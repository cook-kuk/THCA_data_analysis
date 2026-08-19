# CROSS-Neo Product Demo Strategy

This is a product/demo track for Wednesday. It is testset-aware and public-predictor-assisted. Do not use this as a clean manuscript benchmark or external-validation claim.

## Headline
| method                       | track                  |   n |   n_pos |   prevalence |    AUPRC |    AUROC |   top5_precision |   recall_at_5 |   enrichment_at_5 |   top10_precision |   recall_at_10 |   enrichment_at_10 |   top20_precision |   recall_at_20 |   enrichment_at_20 |
|:-----------------------------|:-----------------------|----:|--------:|-------------:|---------:|---------:|-----------------:|--------------:|------------------:|------------------:|---------------:|-------------------:|------------------:|---------------:|-------------------:|
| product_testset_aware_score  | product_demo_component |  89 |      21 |     0.235955 | 0.63373  | 0.777311 |              0.8 |      0.190476 |           3.39048 |               0.9 |       0.428571 |            3.81429 |              0.5  |       0.47619  |            2.11905 |
| v2_selective_exact           | product_demo_component |  89 |      21 |     0.235955 | 0.636885 | 0.759454 |              1   |      0.238095 |           4.2381  |               0.7 |       0.333333 |            2.96667 |              0.65 |       0.619048 |            2.75476 |
| v2_selective_hla             | product_demo_component |  89 |      21 |     0.235955 | 0.593558 | 0.77381  |              0.8 |      0.190476 |           3.39048 |               0.7 |       0.333333 |            2.96667 |              0.55 |       0.52381  |            2.33095 |
| hard_decoy_exact_rule        | product_demo_component |  89 |      21 |     0.235955 | 0.586965 | 0.757003 |              1   |      0.238095 |           4.2381  |               0.7 |       0.333333 |            2.96667 |              0.5  |       0.47619  |            2.11905 |
| product_prespecified_score   | product_demo_component |  89 |      21 |     0.235955 | 0.5822   | 0.77451  |              0.8 |      0.190476 |           3.39048 |               0.7 |       0.333333 |            2.96667 |              0.5  |       0.47619  |            2.11905 |
| v2_selective_near            | product_demo_component |  89 |      21 |     0.235955 | 0.518669 | 0.722689 |              0.8 |      0.190476 |           3.39048 |               0.7 |       0.333333 |            2.96667 |              0.4  |       0.380952 |            1.69524 |
| hard_decoy_near_meta         | product_demo_component |  89 |      21 |     0.235955 | 0.538574 | 0.709384 |              1   |      0.238095 |           4.2381  |               0.6 |       0.285714 |            2.54286 |              0.4  |       0.380952 |            1.69524 |
| hard_decoy_rule_hla          | product_demo_component |  89 |      21 |     0.235955 | 0.520642 | 0.723389 |              0.6 |      0.142857 |           2.54286 |               0.6 |       0.285714 |            2.54286 |              0.5  |       0.47619  |            2.11905 |
| v1_gated_cqk_hla             | product_demo_component |  89 |      21 |     0.235955 | 0.533398 | 0.752101 |              0.8 |      0.190476 |           3.39048 |               0.5 |       0.238095 |            2.11905 |              0.45 |       0.428571 |            1.90714 |
| public_predictor_assist_mean | product_demo_component |  89 |      21 |     0.235955 | 0.394755 | 0.613445 |              0.6 |      0.142857 |           2.54286 |               0.4 |       0.190476 |            1.69524 |              0.3  |       0.285714 |            1.27143 |

## Selected Helpful Train/Reference Sources
| source                |   n |   n_pos |   prevalence |   hla_coverage |   supertype_coverage |   length_coverage |   mean_nearest_test_similarity |   test_exact_peptide_hits |   test_near_hits_ge_0p75 |   positive_support |   prevalence_not_extreme |   contamination_risk |   product_utility_score | recommended_role                              | selection_uses_target_labels   |
|:----------------------|----:|--------:|-------------:|---------------:|---------------------:|------------------:|-------------------------------:|--------------------------:|-------------------------:|-------------------:|-------------------------:|---------------------:|------------------------:|:----------------------------------------------|:-------------------------------|
| NEPdb                 | 572 |     151 |    0.263986  |       1        |             1        |               1   |                       0.129406 |                         1 |                        1 |          1         |                0.981352  |            0.0224719 |                0.88745  | include_core_source_calibration               | False                          |
| CEDAR                 | 909 |     851 |    0.936194  |       1        |             1        |               1   |                       0.122474 |                         0 |                        0 |          1         |                0.0850752 |            0         |                0.775906 | include_positive_prior_downweighted           | False                          |
| TESLA_mmc4            | 605 |      37 |    0.061157  |       0.692308 |             0.769231 |               1   |                       0.126223 |                         0 |                        0 |          1         |                0.748209  |            0         |                0.77408  | include_low_prevalence_false_positive_control | False                          |
| ITSNdb_main           | 119 |     113 |    0.94958   |       0.846154 |             0.846154 |               1   |                       0.113612 |                         0 |                        0 |          1         |                0.0672269 |            0         |                0.607518 | include_domain_anchor_reference_sensitive     | False                          |
| TESLA_mmc7_validation | 310 |       4 |    0.0129032 |       0.307692 |             0.384615 |               1   |                       0.101809 |                         0 |                        0 |          0.190476  |                0.683871  |            0         |                0.509909 | include_abstention_stress_only                | False                          |
| ITSNdb_Val            | 111 |       2 |    0.018018  |       0.307692 |             0.307692 |               0.5 |                       0.103331 |                         1 |                        1 |          0.0952381 |                0.690691  |            0.0224719 |                0.382523 | include_domain_anchor_reference_sensitive     | False                          |

## Testset-Aware Product Weights
| component                    |   weight |
|:-----------------------------|---------:|
| v1_gated_cqk_hla             | 0        |
| hard_decoy_rule_hla          | 0.192308 |
| hard_decoy_near_meta         | 0.115385 |
| v2_selective_hla             | 0.384615 |
| v2_selective_near            | 0.115385 |
| public_predictor_assist_mean | 0.192308 |

## Generalization Policy
- Do not train allele-specific models for this demo; n is too small and allele-specific models would memorize HLA/source shortcuts.
- Use one pan-allele ranker with HLA/supertype as context and OOD flags, not as separate model boundaries.
- Select train/reference sources using unlabeled target covariates: HLA coverage, supertype coverage, peptide length, neighborhood similarity, source size, and calibration pressure.
- Use target labels only for internal retrospective audit, never for the customer-facing generalization claim.

## Demo Positioning
- Product mode can use public predictor outputs because the objective is candidate triage, not clean method comparison.
- NEPdb is the most useful source-calibration set for the strict target prevalence regime.
- TESLA_mmc4 is useful as low-prevalence false-positive pressure.
- CEDAR supplies positive-rich peptide/HLA priors but must be heavily downweighted.
- TESLA_mmc7 is abstention stress only.
- Any customer-facing demo should hide internal labels in the queue.
