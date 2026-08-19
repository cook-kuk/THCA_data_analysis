# CROSS-Neo v1 Hard-Decoy/Focal Repair

Hard decoys are now a bounded auxiliary contrastive confidence signal, not a standalone motif model. Public predictor scores are not used.

## HLA-Stratified Locked
| split_name                 | model                              |   n |   n_pos |   prevalence |    AUPRC |    AUROC |    Brier |   top5_precision |   recall_at_5 |   enrichment_at_5 |   top10_precision |   recall_at_10 |   enrichment_at_10 |
|:---------------------------|:-----------------------------------|----:|--------:|-------------:|---------:|---------:|---------:|-----------------:|--------------:|------------------:|------------------:|---------------:|-------------------:|
| hla_stratified_group_5fold | hard_decoy_rule_aux_C_QK_no_anchor |  89 |      21 |     0.235955 | 0.520642 | 0.723389 | 0.1811   |              0.6 |      0.142857 |           2.54286 |               0.6 |       0.285714 |            2.54286 |
| hla_stratified_group_5fold | hard_decoy_nested_meta_C_QK_aux    |  89 |      21 |     0.235955 | 0.465173 | 0.712885 | 0.172868 |              0.6 |      0.142857 |           2.54286 |               0.5 |       0.238095 |            2.11905 |

## Near-Peptide Locked
| split_name                   | model                              |   n |   n_pos |   prevalence |    AUPRC |    AUROC |    Brier |   top5_precision |   recall_at_5 |   enrichment_at_5 |   top10_precision |   recall_at_10 |   enrichment_at_10 |
|:-----------------------------|:-----------------------------------|----:|--------:|-------------:|---------:|---------:|---------:|-----------------:|--------------:|------------------:|------------------:|---------------:|-------------------:|
| near_peptide_cluster_holdout | hard_decoy_nested_meta_C_QK_aux    |  89 |      21 |     0.235955 | 0.538574 | 0.709384 | 0.178976 |              1   |     0.238095  |           4.2381  |               0.6 |       0.285714 |            2.54286 |
| near_peptide_cluster_holdout | hard_decoy_rule_aux_C_QK_no_anchor |  89 |      21 |     0.235955 | 0.42672  | 0.695378 | 0.187007 |              0.4 |     0.0952381 |           1.69524 |               0.5 |       0.238095 |            2.11905 |

## Source-Stress
| split_name                           | model                                  |   n |   n_pos |   prevalence |     AUPRC |    AUROC |    Brier |   top5_precision |   recall_at_5 |   enrichment_at_5 |   top10_precision |   recall_at_10 |   enrichment_at_10 |
|:-------------------------------------|:---------------------------------------|----:|--------:|-------------:|----------:|---------:|---------:|-----------------:|--------------:|------------------:|------------------:|---------------:|-------------------:|
| source_heldout_CEDAR                 | hard_decoy_sequence_only_source_stress | 909 |     851 |    0.936194  | 0.929679  | 0.462114 | 0.313464 |              0.8 |    0.00470035 |          0.854524 |               0.8 |     0.00940071 |           0.854524 |
| source_heldout_NEPdb                 | hard_decoy_sequence_only_source_stress | 572 |     151 |    0.263986  | 0.425611  | 0.713439 | 0.222431 |              0.8 |    0.0264901  |          3.03046  |               0.7 |     0.0463576  |           2.65166  |
| source_heldout_TESLA_mmc4            | hard_decoy_sequence_only_source_stress | 605 |      37 |    0.061157  | 0.118611  | 0.576561 | 0.195559 |              0.4 |    0.0540541  |          6.54054  |               0.3 |     0.0810811  |           4.90541  |
| source_heldout_TESLA_mmc7_validation | hard_decoy_sequence_only_source_stress | 310 |       4 |    0.0129032 | 0.0825419 | 0.70098  | 0.190101 |              0.2 |    0.25       |         15.5      |               0.1 |     0.25       |           7.75     |

## Decoy Manifest
| split_name                           | fold_id               |   train_n |   train_pos |   hard_decoy_n |   pattern_reliability |
|:-------------------------------------|:----------------------|----------:|------------:|---------------:|----------------------:|
| repeated_stratified_5x5_internal     | rs_00                 |        71 |          16 |            125 |              0.75     |
| repeated_stratified_5x5_internal     | rs_01                 |        71 |          17 |            127 |              0.75     |
| repeated_stratified_5x5_internal     | rs_02                 |        71 |          17 |            135 |              0.75     |
| repeated_stratified_5x5_internal     | rs_03                 |        71 |          17 |            132 |              0.75     |
| repeated_stratified_5x5_internal     | rs_04                 |        72 |          17 |            136 |              0.75     |
| repeated_stratified_5x5_internal     | rs_05                 |        71 |          16 |            126 |              0.75     |
| repeated_stratified_5x5_internal     | rs_06                 |        71 |          17 |            129 |              0.75     |
| repeated_stratified_5x5_internal     | rs_07                 |        71 |          17 |            130 |              0.75     |
| repeated_stratified_5x5_internal     | rs_08                 |        71 |          17 |            130 |              0.75     |
| repeated_stratified_5x5_internal     | rs_09                 |        72 |          17 |            136 |              0.75     |
| repeated_stratified_5x5_internal     | rs_10                 |        71 |          16 |            126 |              0.75     |
| repeated_stratified_5x5_internal     | rs_11                 |        71 |          17 |            126 |              0.75     |
| repeated_stratified_5x5_internal     | rs_12                 |        71 |          17 |            132 |              0.75     |
| repeated_stratified_5x5_internal     | rs_13                 |        71 |          17 |            131 |              0.75     |
| repeated_stratified_5x5_internal     | rs_14                 |        72 |          17 |            134 |              0.75     |
| repeated_stratified_5x5_internal     | rs_15                 |        71 |          16 |            126 |              0.75     |
| repeated_stratified_5x5_internal     | rs_16                 |        71 |          17 |            125 |              0.75     |
| repeated_stratified_5x5_internal     | rs_17                 |        71 |          17 |            129 |              0.75     |
| repeated_stratified_5x5_internal     | rs_18                 |        71 |          17 |            130 |              0.75     |
| repeated_stratified_5x5_internal     | rs_19                 |        72 |          17 |            133 |              0.75     |
| repeated_stratified_5x5_internal     | rs_20                 |        71 |          16 |            125 |              0.75     |
| repeated_stratified_5x5_internal     | rs_21                 |        71 |          17 |            130 |              0.75     |
| repeated_stratified_5x5_internal     | rs_22                 |        71 |          17 |            129 |              0.75     |
| repeated_stratified_5x5_internal     | rs_23                 |        71 |          17 |            125 |              0.75     |
| repeated_stratified_5x5_internal     | rs_24                 |        72 |          17 |            135 |              0.75     |
| hla_stratified_group_5fold           | hla_00                |        70 |          19 |            151 |              0.75     |
| hla_stratified_group_5fold           | hla_01                |        56 |          11 |             82 |              0.75     |
| hla_stratified_group_5fold           | hla_02                |        73 |          19 |            149 |              0.75     |
| hla_stratified_group_5fold           | hla_03                |        79 |          17 |            134 |              0.75     |
| hla_stratified_group_5fold           | hla_04                |        78 |          18 |            142 |              0.75     |
| hla_supertype_heldout                | super_00              |        84 |          19 |            142 |              0.75     |
| hla_supertype_heldout                | super_01              |        56 |          11 |             82 |              0.75     |
| hla_supertype_heldout                | super_02              |        74 |          20 |            158 |              0.75     |
| hla_supertype_heldout                | super_09              |        70 |          19 |            151 |              0.75     |
| near_peptide_cluster_holdout         | near_00               |        71 |          17 |            133 |              0.75     |
| near_peptide_cluster_holdout         | near_01               |        70 |          16 |            128 |              0.75     |
| near_peptide_cluster_holdout         | near_02               |        72 |          17 |            133 |              0.75     |
| near_peptide_cluster_holdout         | near_03               |        72 |          17 |            130 |              0.75     |
| near_peptide_cluster_holdout         | near_04               |        71 |          17 |            130 |              0.75     |
| exact_peptide_hla_holdout            | pmhc_00               |        71 |          17 |            132 |              0.75     |
| exact_peptide_hla_holdout            | pmhc_01               |        71 |          17 |            130 |              0.75     |
| exact_peptide_hla_holdout            | pmhc_02               |        71 |          16 |            126 |              0.75     |
| exact_peptide_hla_holdout            | pmhc_03               |        72 |          17 |            136 |              0.75     |
| exact_peptide_hla_holdout            | pmhc_04               |        71 |          17 |            128 |              0.75     |
| study_heldout                        | study_00              |        80 |          16 |            125 |              0.75     |
| study_heldout                        | study_01              |         9 |           5 |             39 |              0.396889 |
| source_heldout_CEDAR                 | CEDAR                 |      1487 |         192 |           1491 |              0.75     |
| source_heldout_NEPdb                 | NEPdb                 |      1824 |         892 |           6842 |              0.75     |
| source_heldout_TESLA_mmc4            | TESLA_mmc4            |      1791 |        1006 |           7734 |              0.75     |
| source_heldout_TESLA_mmc7_validation | TESLA_mmc7_validation |      2086 |        1039 |           7952 |              0.75     |

## Interpretation
- Use the rule auxiliary model if it improves locked AUPRC/top-k without hurting near-peptide holdout.
- Use sequence-only source stress only as a business triage diagnostic, not external validation.
- Fold safety: decoys and motif statistics are created from train rows only.

## Fold-Averaged Positive Pattern
|   position_1based | top_positive_enriched_residues   |   mean_positive_entropy |   mean_background_entropy |
|------------------:|:---------------------------------|------------------------:|--------------------------:|
|                 1 | E,D,M,N,P                        |                 4.05244 |                   3.96327 |
|                 2 | C,D,E,G,H                        |                 3.876   |                   3.50293 |
|                 3 | C,E,R,D,L                        |                 3.91918 |                   4.07932 |
|                 4 | W,V,E,G,K                        |                 3.99452 |                   4.06835 |
|                 5 | M,Y,W,T,G                        |                 4.03795 |                   4.07305 |
|                 6 | D,Q,V,H,E                        |                 3.91121 |                   4.01068 |
|                 7 | R,W,I,E,L                        |                 4.05987 |                   4.11413 |
|                 8 | C,H,N,W,L                        |                 4.0577  |                   4.04552 |
|                 9 | C,D,H,N,P                        |                 4.0121  |                   3.74089 |
|                10 | A,C,D,E,G                        |                 4.21183 |                   4.0873  |
|                11 | A,C,D,E,F                        |                 4.32193 |                   4.32193 |
|                12 | A,C,D,E,F                        |                 4.32193 |                   4.32193 |
|                13 | A,C,D,E,F                        |                 4.32193 |                   4.32193 |
|                14 | A,C,D,E,F                        |                 4.32193 |                   4.32193 |
|                15 | A,C,D,E,F                        |                 4.32193 |                   4.32193 |
