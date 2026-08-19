# CROSS-Neo v0 Decision Report

This is an internal/locked split report. It does not claim external validation or quantum advantage.

## Best Internal Model

- split: `hla_stratified_group_5fold`
- feature group: `C_counterfactual`
- model: `rf_secondary`
- AUPRC: 0.388
- AUROC: 0.628
- top10 precision: 0.400
- enrichment@10: 1.695

## Why CROSS-Neo Exists

- `Structure_LR` is clean but shallow.
- `Wave8` is strong but reference-sensitive.
- `QK-NoAnchor` is promising but small-n internal only.
- CROSS-Neo separates retrieval evidence, structure geometry, counterfactual peptide/HLA encoding, quantum fixed features, and OOD abstention so gains can be audited.

## Strict/Internal Top Models

| split_name                       | feature_group                 | model          |   n |   n_pos |    AUPRC |    AUROC |   top10_precision |   enrichment_at_10 |
|:---------------------------------|:------------------------------|:---------------|----:|--------:|---------:|---------:|------------------:|-------------------:|
| repeated_stratified_5x5_internal | C_counterfactual              | rf_secondary   | 445 |     105 | 0.352871 | 0.609916 |               0.5 |            2.11905 |
| repeated_stratified_5x5_internal | A_structure_baseline_features | rf_secondary   | 445 |     105 | 0.286703 | 0.575994 |               0   |            0       |
| repeated_stratified_5x5_internal | F_CROSS_all                   | rf_secondary   | 445 |     105 | 0.285343 | 0.58591  |               0.3 |            1.27143 |
| repeated_stratified_5x5_internal | C_counterfactual              | elastic_net_lr | 445 |     105 | 0.271118 | 0.581261 |               0   |            0       |
| repeated_stratified_5x5_internal | A_structure_baseline_features | elastic_net_lr | 445 |     105 | 0.24624  | 0.534454 |               0   |            0       |
| repeated_stratified_5x5_internal | F_CROSS_all                   | elastic_net_lr | 445 |     105 | 0.243632 | 0.537647 |               0.1 |            0.42381 |
| repeated_stratified_5x5_internal | E_quantum_fixed               | elastic_net_lr | 445 |     105 | 0.238189 | 0.481176 |               0.4 |            1.69524 |
| repeated_stratified_5x5_internal | B_retrieval_only_clean_flags  | elastic_net_lr | 445 |     105 | 0.23553  | 0.496429 |               0.4 |            1.69524 |
| repeated_stratified_5x5_internal | B_retrieval_only_clean_flags  | rf_secondary   | 445 |     105 | 0.234566 | 0.497997 |               0.1 |            0.42381 |
| repeated_stratified_5x5_internal | E_quantum_fixed               | rf_secondary   | 445 |     105 | 0.222928 | 0.489636 |               0.1 |            0.42381 |

## HLA Robustness Snapshot

| split_name                 | feature_group                 | model          |   n |   n_pos |    AUPRC |    AUROC |   top10_precision |   enrichment_at_10 |
|:---------------------------|:------------------------------|:---------------|----:|--------:|---------:|---------:|------------------:|-------------------:|
| hla_stratified_group_5fold | C_counterfactual              | rf_secondary   |  89 |      21 | 0.388228 | 0.628151 |               0.4 |            1.69524 |
| hla_stratified_group_5fold | A_structure_baseline_features | rf_secondary   |  89 |      21 | 0.372024 | 0.544118 |               0.3 |            1.27143 |
| hla_supertype_heldout      | C_counterfactual              | rf_secondary   |  72 |      15 | 0.368256 | 0.691228 |               0.4 |            1.92    |
| hla_supertype_heldout      | A_structure_baseline_features | rf_secondary   |  72 |      15 | 0.354303 | 0.499415 |               0.3 |            1.44    |
| hla_stratified_group_5fold | C_counterfactual              | elastic_net_lr |  89 |      21 | 0.294583 | 0.582633 |               0.1 |            0.42381 |

## Quantum Fallback

| split_name                       | branch                 |   n |   n_pos |    AUPRC |    AUROC |
|:---------------------------------|:-----------------------|----:|--------:|---------:|---------:|
| near_peptide_cluster_holdout     | qk_quantum_only_gamma1 |  89 |      21 | 0.459695 | 0.678571 |
| hla_supertype_heldout            | qk_quantum_only_gamma1 |  72 |      15 | 0.449836 | 0.635088 |
| hla_stratified_group_5fold       | qk_quantum_only_gamma1 |  89 |      21 | 0.426975 | 0.653361 |
| hla_stratified_group_5fold       | qk_no_anchor_gamma1    |  89 |      21 | 0.42504  | 0.665966 |
| repeated_stratified_5x5_internal | qk_no_anchor_gamma1    | 445 |     105 | 0.375286 | 0.632941 |
| hla_supertype_heldout            | qk_no_anchor_gamma1    |  72 |      15 | 0.371807 | 0.671345 |
| repeated_stratified_5x5_internal | qk_quantum_only_gamma1 | 445 |     105 | 0.349639 | 0.619832 |
| near_peptide_cluster_holdout     | qk_no_anchor_gamma1    |  89 |      21 | 0.34541  | 0.659664 |
| repeated_stratified_5x5_internal | qk_clean_no_tcr_gamma1 | 445 |     105 | 0.312195 | 0.621653 |
| hla_stratified_group_5fold       | qk_clean_no_tcr_gamma1 |  89 |      21 | 0.30988  | 0.614146 |

## Fixed Late Fusion

| split_name                       | fusion                                              | status                    |   n |   n_pos |   prevalence |    AUPRC |    AUROC |    Brier |   top5_precision |   top10_precision |
|:---------------------------------|:----------------------------------------------------|:--------------------------|----:|--------:|-------------:|---------:|---------:|---------:|-----------------:|------------------:|
| near_peptide_cluster_holdout     | counterfactual_rf_plus_qk_quantum_only_gamma1_w0.25 | exploratory               |  89 |      21 |     0.235955 | 0.484806 | 0.689776 | 0.209057 |              0.8 |               0.6 |
| hla_stratified_group_5fold       | counterfactual_rf_plus_qk_quantum_only_gamma1_w0.75 | exploratory               |  89 |      21 |     0.235955 | 0.474182 | 0.712885 | 0.195652 |              0.6 |               0.4 |
| hla_stratified_group_5fold       | counterfactual_rf_plus_qk_no_anchor_gamma1_w0.5     | prespecified_equal_weight |  89 |      21 |     0.235955 | 0.472913 | 0.691877 | 0.197007 |              0.8 |               0.5 |
| hla_supertype_heldout            | counterfactual_rf_plus_qk_quantum_only_gamma1_w0.75 | exploratory               |  72 |      15 |     0.208333 | 0.470344 | 0.721637 | 0.188199 |              0.6 |               0.4 |
| hla_supertype_heldout            | counterfactual_rf_plus_qk_quantum_only_gamma1_w0.25 | exploratory               |  72 |      15 |     0.208333 | 0.456347 | 0.646784 | 0.226131 |              0.6 |               0.5 |
| hla_supertype_heldout            | counterfactual_rf_plus_qk_quantum_only_gamma1_w0.5  | prespecified_equal_weight |  72 |      15 |     0.208333 | 0.454616 | 0.679532 | 0.202305 |              0.6 |               0.4 |
| hla_stratified_group_5fold       | counterfactual_rf_plus_qk_no_anchor_gamma1_w0.75    | exploratory               |  89 |      21 |     0.235955 | 0.4524   | 0.679972 | 0.193009 |              0.8 |               0.5 |
| hla_stratified_group_5fold       | counterfactual_rf_plus_qk_quantum_only_gamma1_w0.5  | prespecified_equal_weight |  89 |      21 |     0.235955 | 0.442783 | 0.693978 | 0.204889 |              0.4 |               0.4 |
| hla_stratified_group_5fold       | counterfactual_rf_plus_qk_quantum_only_gamma1_w0.25 | exploratory               |  89 |      21 |     0.235955 | 0.436284 | 0.670868 | 0.223052 |              0.6 |               0.5 |
| repeated_stratified_5x5_internal | counterfactual_rf_plus_qk_quantum_only_gamma1_w0.5  | prespecified_equal_weight | 445 |     105 |     0.235955 | 0.432176 | 0.64465  | 0.203766 |              1   |               0.7 |
| hla_stratified_group_5fold       | counterfactual_rf_plus_qk_no_anchor_gamma1_w0.25    | exploratory               |  89 |      21 |     0.235955 | 0.428628 | 0.668067 | 0.207332 |              0.6 |               0.4 |
| repeated_stratified_5x5_internal | counterfactual_rf_plus_qk_quantum_only_gamma1_w0.75 | exploratory               | 445 |     105 |     0.235955 | 0.428313 | 0.634146 | 0.19543  |              1   |               1   |

## Abstention

|   coverage |   kept_n |    AUPRC |   top10_precision |   enrichment_at_10 |
|-----------:|---------:|---------:|------------------:|-------------------:|
|        1   |       89 | 0.388228 |               0.4 |           1.69524  |
|        0.9 |       80 | 0.343207 |               0.4 |           1.88235  |
|        0.8 |       71 | 0.364602 |               0.4 |           1.775    |
|        0.7 |       62 | 0.307415 |               0.2 |           0.953846 |
|        0.6 |       53 | 0.272091 |               0.2 |           0.963636 |
|        0.5 |       44 | 0.346775 |               0.5 |           2.44444  |
|        0.4 |       36 | 0.249722 |               0.2 |           1.44     |
|        0.3 |       27 | 0.145262 |               0.1 |           0.9      |

## Leakage Controls

- Public predictor scores (`MHCflurry`, `NetMHCpan`, `BigMHC`, `PRIME`, `MixMHCpred`, `NetMHCstabpan`) were not used as model features.
- Retrieval features were recomputed train-fold only.
- Scaling/calibration/model fitting occurred inside each train fold.
- Fixed quantum gamma was used; no test-fold gamma search.
- Source-window, patient, and time splits are reported as unavailable when metadata are missing.

## Decision

Decision: **HOLD**.

Promote only after the same feature groups survive a true external/time/study split and a public-corpus overlap audit.
