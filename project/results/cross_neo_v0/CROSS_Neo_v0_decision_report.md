# CROSS-Neo v0 Decision Report

This is an internal/locked split report. It does not claim external validation or quantum advantage.

## Best Eligible Internal Model

Eligibility for this headline row requires `n >= 50`, so tiny heldout artifacts do not dominate the decision.

- split: `exact_peptide_hla_holdout`
- feature group: `C_counterfactual`
- model: `rf_secondary`
- n: 89, positives: 21, prevalence: 0.236
- AUPRC: 0.525
- AUROC: 0.715
- top10 precision: 0.600
- enrichment@10: 2.543

## Raw Best Flag

Raw maximum AUPRC is `E_quantum_fixed / xgboost_depth2 / study_heldout` with n=9, AUPRC=0.811, AUROC=0.600. Treat this as descriptive only when n is small.

## Why CROSS-Neo Exists

- `Structure_LR` is clean but shallow.
- `Wave8` is strong but reference-sensitive.
- `QK-NoAnchor` is promising but small-n internal only.
- CROSS-Neo separates retrieval evidence, structure geometry, counterfactual peptide/HLA encoding, quantum fixed features, and OOD abstention so gains can be audited.

## Strict/Internal Top Models

| split_name                       | feature_group                 | model                 |   n |   n_pos |    AUPRC |    AUROC |   top10_precision |   enrichment_at_10 |
|:---------------------------------|:------------------------------|:----------------------|----:|--------:|---------:|---------:|------------------:|-------------------:|
| repeated_stratified_5x5_internal | C_counterfactual              | elastic_net_lr        | 445 |     105 | 0.460206 | 0.646246 |               0.8 |            3.39048 |
| repeated_stratified_5x5_internal | F_CROSS_all                   | elastic_net_lr        | 445 |     105 | 0.42474  | 0.663109 |               0.7 |            2.96667 |
| repeated_stratified_5x5_internal | C_counterfactual              | rf_secondary          | 445 |     105 | 0.407385 | 0.687619 |               0.5 |            2.11905 |
| repeated_stratified_5x5_internal | C_counterfactual              | xgboost_depth2        | 445 |     105 | 0.383647 | 0.682325 |               0.4 |            1.69524 |
| repeated_stratified_5x5_internal | C_counterfactual              | calibrated_linear_svm | 445 |     105 | 0.351517 | 0.603277 |               0.4 |            1.69524 |
| repeated_stratified_5x5_internal | F_CROSS_all                   | xgboost_depth2        | 445 |     105 | 0.33565  | 0.633669 |               0.5 |            2.11905 |
| repeated_stratified_5x5_internal | F_CROSS_all                   | calibrated_linear_svm | 445 |     105 | 0.326984 | 0.613585 |               0.5 |            2.11905 |
| repeated_stratified_5x5_internal | F_CROSS_all                   | rf_secondary          | 445 |     105 | 0.316061 | 0.628683 |               0.3 |            1.27143 |
| repeated_stratified_5x5_internal | A_structure_baseline_features | rf_secondary          | 445 |     105 | 0.286703 | 0.575994 |               0   |            0       |
| repeated_stratified_5x5_internal | A_structure_baseline_features | xgboost_depth2        | 445 |     105 | 0.274559 | 0.54056  |               0.1 |            0.42381 |

## Strict No-Reference / Near-Cluster Holdout

| split_name                   | feature_group    | model                 |   n |   n_pos |   prevalence |    AUPRC |    AUROC |   top10_precision |   enrichment_at_10 |
|:-----------------------------|:-----------------|:----------------------|----:|--------:|-------------:|---------:|---------:|------------------:|-------------------:|
| near_peptide_cluster_holdout | C_counterfactual | rf_secondary          |  89 |      21 |     0.235955 | 0.441926 | 0.671569 |               0.5 |            2.11905 |
| near_peptide_cluster_holdout | C_counterfactual | elastic_net_lr        |  89 |      21 |     0.235955 | 0.434058 | 0.630952 |               0.6 |            2.54286 |
| near_peptide_cluster_holdout | F_CROSS_all      | elastic_net_lr        |  89 |      21 |     0.235955 | 0.429078 | 0.620448 |               0.4 |            1.69524 |
| near_peptide_cluster_holdout | C_counterfactual | calibrated_linear_svm |  89 |      21 |     0.235955 | 0.368858 | 0.542017 |               0.4 |            1.69524 |
| near_peptide_cluster_holdout | F_CROSS_all      | rf_secondary          |  89 |      21 |     0.235955 | 0.341079 | 0.566527 |               0.4 |            1.69524 |
| near_peptide_cluster_holdout | E_quantum_fixed  | elastic_net_lr        |  89 |      21 |     0.235955 | 0.338528 | 0.576331 |               0.3 |            1.27143 |
| near_peptide_cluster_holdout | F_CROSS_all      | xgboost_depth2        |  89 |      21 |     0.235955 | 0.329462 | 0.592437 |               0.4 |            1.69524 |
| near_peptide_cluster_holdout | E_quantum_fixed  | xgboost_depth2        |  89 |      21 |     0.235955 | 0.328673 | 0.560224 |               0.4 |            1.69524 |

## HLA Robustness Snapshot

| split_name                 | feature_group    | model          |   n |   n_pos |    AUPRC |    AUROC |   top10_precision |   enrichment_at_10 |
|:---------------------------|:-----------------|:---------------|----:|--------:|---------:|---------:|------------------:|-------------------:|
| hla_stratified_group_5fold | C_counterfactual | rf_secondary   |  89 |      21 | 0.50549  | 0.717087 |               0.6 |            2.54286 |
| hla_supertype_heldout      | C_counterfactual | elastic_net_lr |  72 |      15 | 0.496912 | 0.738012 |               0.5 |            2.4     |
| hla_stratified_group_5fold | C_counterfactual | elastic_net_lr |  89 |      21 | 0.492459 | 0.679972 |               0.6 |            2.54286 |
| hla_supertype_heldout      | C_counterfactual | rf_secondary   |  72 |      15 | 0.477884 | 0.747368 |               0.5 |            2.4     |
| hla_stratified_group_5fold | F_CROSS_all      | elastic_net_lr |  89 |      21 | 0.427589 | 0.65056  |               0.5 |            2.11905 |

## Quantum Fallback

| split_name                       | branch                 |   n |   n_pos |    AUPRC |    AUROC |
|:---------------------------------|:-----------------------|----:|--------:|---------:|---------:|
| near_peptide_cluster_holdout     | qk_quantum_only_gamma1 |  89 |      21 | 0.459695 | 0.678571 |
| exact_peptide_hla_holdout        | qk_no_anchor_gamma1    |  89 |      21 | 0.452113 | 0.703081 |
| hla_supertype_heldout            | qk_quantum_only_gamma1 |  72 |      15 | 0.449836 | 0.635088 |
| exact_peptide_hla_holdout        | qk_clean_no_tcr_gamma1 |  89 |      21 | 0.44444  | 0.720588 |
| exact_peptide_hla_holdout        | qk_quantum_only_gamma1 |  89 |      21 | 0.432427 | 0.654062 |
| hla_stratified_group_5fold       | qk_quantum_only_gamma1 |  89 |      21 | 0.426975 | 0.653361 |
| hla_stratified_group_5fold       | qk_no_anchor_gamma1    |  89 |      21 | 0.42504  | 0.665966 |
| repeated_stratified_5x5_internal | qk_no_anchor_gamma1    | 445 |     105 | 0.375286 | 0.632941 |
| hla_supertype_heldout            | qk_no_anchor_gamma1    |  72 |      15 | 0.371807 | 0.671345 |
| repeated_stratified_5x5_internal | qk_quantum_only_gamma1 | 445 |     105 | 0.349639 | 0.619832 |

## Fixed Late Fusion

| split_name                 | fusion                                              | status                    |   n |   n_pos |   prevalence |    AUPRC |    AUROC |    Brier |   top5_precision |   top10_precision |
|:---------------------------|:----------------------------------------------------|:--------------------------|----:|--------:|-------------:|---------:|---------:|---------:|-----------------:|------------------:|
| study_heldout              | counterfactual_rf_plus_qk_no_anchor_gamma1_w0.75    | exploratory               |   9 |       5 |     0.555556 | 0.644444 | 0.4      | 0.314121 |              0.4 |          0.555556 |
| study_heldout              | counterfactual_rf_plus_qk_quantum_only_gamma1_w0.5  | prespecified_equal_weight |   9 |       5 |     0.555556 | 0.596825 | 0.3      | 0.31676  |              0.4 |          0.555556 |
| exact_peptide_hla_holdout  | counterfactual_rf_plus_qk_no_anchor_gamma1_w0.75    | exploratory               |  89 |      21 |     0.235955 | 0.596809 | 0.768207 | 0.178298 |              0.8 |          0.7      |
| study_heldout              | counterfactual_rf_plus_qk_quantum_only_gamma1_w0.75 | exploratory               |   9 |       5 |     0.555556 | 0.591111 | 0.3      | 0.312124 |              0.4 |          0.555556 |
| study_heldout              | counterfactual_rf_plus_qk_clean_no_tcr_gamma1_w0.75 | exploratory               |   9 |       5 |     0.555556 | 0.591111 | 0.3      | 0.317647 |              0.4 |          0.555556 |
| exact_peptide_hla_holdout  | counterfactual_rf_plus_qk_no_anchor_gamma1_w0.5     | prespecified_equal_weight |  89 |      21 |     0.235955 | 0.586965 | 0.757003 | 0.184277 |              1   |          0.7      |
| study_heldout              | counterfactual_rf_plus_qk_clean_no_tcr_gamma1_w0.5  | prespecified_equal_weight |   9 |       5 |     0.555556 | 0.563492 | 0.2      | 0.326366 |              0.2 |          0.555556 |
| hla_stratified_group_5fold | counterfactual_rf_plus_qk_quantum_only_gamma1_w0.75 | exploratory               |  89 |      21 |     0.235955 | 0.555039 | 0.761204 | 0.18385  |              0.6 |          0.6      |
| exact_peptide_hla_holdout  | counterfactual_rf_plus_qk_quantum_only_gamma1_w0.75 | exploratory               |  89 |      21 |     0.235955 | 0.553634 | 0.766106 | 0.181363 |              0.8 |          0.8      |
| hla_supertype_heldout      | counterfactual_rf_plus_qk_quantum_only_gamma1_w0.75 | exploratory               |  72 |      15 |     0.208333 | 0.546529 | 0.768421 | 0.178696 |              0.6 |          0.6      |
| exact_peptide_hla_holdout  | counterfactual_rf_plus_qk_clean_no_tcr_gamma1_w0.5  | prespecified_equal_weight |  89 |      21 |     0.235955 | 0.545582 | 0.765406 | 0.186388 |              0.6 |          0.6      |
| exact_peptide_hla_holdout  | counterfactual_rf_plus_qk_clean_no_tcr_gamma1_w0.75 | exploratory               |  89 |      21 |     0.235955 | 0.5452   | 0.757703 | 0.179308 |              0.8 |          0.6      |

## Source-Heldout Stress

| heldout_study         | method                                   |   n |   n_pos |   prevalence |     AUPRC |    AUROC |    Brier |   top5_precision |   recall_at_5 |   enrichment_at_5 |   top10_precision |   recall_at_10 |   enrichment_at_10 |
|:----------------------|:-----------------------------------------|----:|--------:|-------------:|----------:|---------:|---------:|-----------------:|--------------:|------------------:|------------------:|---------------:|-------------------:|
| CEDAR                 | sourceheld_counterfactual_rf             | 913 |     851 |    0.932092  | 0.918082  | 0.416341 | 0.3709   |              0.8 |    0.00470035 |          0.858284 |               0.9 |      0.0105758 |            0.96557 |
| NEPdb                 | sourceheld_counterfactual_rf             | 886 |     354 |    0.399549  | 0.40534   | 0.49521  | 0.261428 |              0   |    0          |          0        |               0   |      0         |            0       |
| TESLA_mmc4            | sourceheld_counterfactual_rf             | 610 |      37 |    0.0606557 | 0.0650109 | 0.549267 | 0.24215  |              0   |    0          |          0        |               0   |      0         |            0       |
| TESLA_mmc7_validation | sourceheld_counterfactual_rf             | 319 |       6 |    0.0188088 | 0.0378096 | 0.655485 | 0.212548 |              0   |    0          |          0        |               0   |      0         |            0       |
| CEDAR                 | sourceheld_prespecified_late_fusion_w0.5 | 913 |     851 |    0.932092  | 0.921816  | 0.439028 | 0.413353 |              1   |    0.00587544 |          1.07286  |               1   |      0.0117509 |            1.07286 |
| NEPdb                 | sourceheld_prespecified_late_fusion_w0.5 | 886 |     354 |    0.399549  | 0.395509  | 0.479297 | 0.293726 |              0   |    0          |          0        |               0   |      0         |            0       |
| TESLA_mmc4            | sourceheld_prespecified_late_fusion_w0.5 | 610 |      37 |    0.0606557 | 0.0551107 | 0.464931 | 0.31171  |              0   |    0          |          0        |               0   |      0         |            0       |
| TESLA_mmc7_validation | sourceheld_prespecified_late_fusion_w0.5 | 319 |       6 |    0.0188088 | 0.0168969 | 0.294462 | 0.241671 |              0   |    0          |          0        |               0   |      0         |            0       |
| CEDAR                 | sourceheld_qk_compact_gamma1             | 913 |     851 |    0.932092  | 0.926744  | 0.461317 | 0.527638 |              0.8 |    0.00470035 |          0.858284 |               0.9 |      0.0105758 |            0.96557 |
| NEPdb                 | sourceheld_qk_compact_gamma1             | 886 |     354 |    0.399549  | 0.406426  | 0.482637 | 0.401408 |              0.4 |    0.00564972 |          1.00113  |               0.4 |      0.0112994 |            1.00113 |
| TESLA_mmc4            | sourceheld_qk_compact_gamma1             | 610 |      37 |    0.0606557 | 0.0516625 | 0.437904 | 0.476714 |              0   |    0          |          0        |               0   |      0         |            0       |
| TESLA_mmc7_validation | sourceheld_qk_compact_gamma1             | 319 |       6 |    0.0188088 | 0.0145965 | 0.180511 | 0.350284 |              0   |    0          |          0        |               0   |      0         |            0       |

## Public Overlap Status

| source     | status                       | note                                                                          |
|:-----------|:-----------------------------|:------------------------------------------------------------------------------|
| MHCflurry  | unresolved_public_pretrained | IEDB/MS-ligand/affinity training corpus not row-audited locally               |
| NetMHCpan  | unresolved_public_pretrained | BA/EL public training corpus not row-audited locally                          |
| BigMHC     | unresolved_public_pretrained | presentation/immunogenicity release exists but not downloaded into this audit |
| PRIME      | unresolved_public_pretrained | public immunogenicity training set not row-audited locally                    |
| MixMHCpred | unresolved_public_pretrained | public ligand training corpus not row-audited locally                         |
| CEDAR      | local_train_pool_present     | exact/near overlap against local train pool auditable                         |
| TESLA      | local_train_pool_present     | exact/near overlap against local train pool auditable                         |
| NEPdb      | local_train_pool_present     | exact/near overlap against local train pool auditable                         |

## Abstention

|   coverage |   kept_n |    AUPRC |   top10_precision |   enrichment_at_10 |
|-----------:|---------:|---------:|------------------:|-------------------:|
|        1   |       89 | 0.441926 |               0.5 |            2.11905 |
|        0.9 |       80 | 0.359444 |               0.5 |            2.35294 |
|        0.8 |       71 | 0.261842 |               0.2 |            1.09231 |
|        0.7 |       62 | 0.24344  |               0.4 |            2.48    |
|        0.6 |       53 | 0.479054 |               0.3 |            1.9875  |
|        0.5 |       44 | 0.103918 |               0   |            0       |
|        0.4 |       36 | 0.150677 |               0.1 |            0.72    |
|        0.3 |       27 | 0.212632 |               0.2 |            1.35    |

## Leakage Controls

- Public predictor scores (`MHCflurry`, `NetMHCpan`, `BigMHC`, `PRIME`, `MixMHCpred`, `NetMHCstabpan`) were not used as model features.
- Retrieval features were recomputed train-fold only.
- Scaling/calibration/model fitting occurred inside each train fold.
- Fixed quantum gamma was used; no test-fold gamma search.
- Source-window, patient, and time splits are reported as unavailable when metadata are missing.

## Decision

Decision: **HOLD**.

Reason: eligible internal/HLA metrics improved over the shallow structure baseline, but source-heldout NEPdb/TESLA performance collapses at top-k. Promote only after the same feature groups survive a true external/time/study split and a public-corpus overlap audit.
