# CROSS-Neo v1 Decoy/Focal Positive-Pattern Report

This branch implements the earlier decoy-negative idea in a fold-safe way. Decoys and positive PWM/log-odds features are generated from each train fold only.

## Metrics
| split_name                           | model                        |   n |   n_pos |   prevalence |     AUPRC |    AUROC |     Brier |   top5_precision |   recall_at_5 |   enrichment_at_5 |   top10_precision |   recall_at_10 |   enrichment_at_10 |
|:-------------------------------------|:-----------------------------|----:|--------:|-------------:|----------:|---------:|----------:|-----------------:|--------------:|------------------:|------------------:|---------------:|-------------------:|
| source_heldout_CEDAR                 | decoy_focal_positive_pattern | 909 |     851 |    0.936194  | 0.932369  | 0.508935 | 0.916048  |              1   |    0.00587544 |          1.06816  |          0.9      |     0.0105758  |            0.96134 |
| study_heldout                        | decoy_focal_positive_pattern |   9 |       5 |    0.555556  | 0.71      | 0.75     | 0.51542   |              0.8 |    0.8        |          1.44     |          0.555556 |     1          |            1       |
| source_heldout_NEPdb                 | decoy_focal_positive_pattern | 572 |     151 |    0.263986  | 0.395733  | 0.658964 | 0.261779  |              1   |    0.0331126  |          3.78808  |          0.6      |     0.0397351  |            2.27285 |
| hla_supertype_heldout                | decoy_focal_positive_pattern |  72 |      15 |    0.208333  | 0.371964  | 0.597661 | 0.171646  |              0.4 |    0.133333   |          1.92     |          0.3      |     0.2        |            1.44    |
| near_peptide_cluster_holdout         | decoy_focal_positive_pattern |  89 |      21 |    0.235955  | 0.284199  | 0.586835 | 0.212771  |              0.2 |    0.047619   |          0.847619 |          0.1      |     0.047619   |            0.42381 |
| hla_stratified_group_5fold           | decoy_focal_positive_pattern |  89 |      21 |    0.235955  | 0.258583  | 0.544818 | 0.212995  |              0.2 |    0.047619   |          0.847619 |          0.1      |     0.047619   |            0.42381 |
| repeated_stratified_5x5_internal     | decoy_focal_positive_pattern | 445 |     105 |    0.235955  | 0.210893  | 0.443838 | 0.222558  |              0   |    0          |          0        |          0.1      |     0.00952381 |            0.42381 |
| exact_peptide_hla_holdout            | decoy_focal_positive_pattern |  89 |      21 |    0.235955  | 0.202504  | 0.415966 | 0.225639  |              0.2 |    0.047619   |          0.847619 |          0.1      |     0.047619   |            0.42381 |
| source_heldout_TESLA_mmc7_validation | decoy_focal_positive_pattern | 310 |       4 |    0.0129032 | 0.186724  | 0.961601 | 0.0225295 |              0.2 |    0.25       |         15.5      |          0.2      |     0.5        |           15.5     |
| source_heldout_TESLA_mmc4            | decoy_focal_positive_pattern | 605 |      37 |    0.061157  | 0.0717574 | 0.523268 | 0.0628048 |              0   |    0          |          0        |          0        |     0          |            0       |

## Decoy Manifest Summary
|        | split_name                       | fold_id   |   train_n |   train_pos |   decoy_n |   features_n |
|:-------|:---------------------------------|:----------|----------:|------------:|----------:|-------------:|
| count  | 45                               | 45        |  45       |    45       |  45       |           45 |
| unique | 6                                | 45        | nan       |   nan       | nan       |          nan |
| top    | repeated_stratified_5x5_internal | rs_00     | nan       |   nan       | nan       |          nan |
| freq   | 25                               | 1         | nan       |   nan       | nan       |          nan |
| mean   | nan                              | nan       |  71.3778  |    16.8222  |  67.2889  |          210 |
| std    | nan                              | nan       |   4.34474 |     1.54168 |   6.16671 |            0 |
| min    | nan                              | nan       |  56       |    11       |  44       |          210 |
| 25%    | nan                              | nan       |  71       |    17       |  68       |          210 |
| 50%    | nan                              | nan       |  71       |    17       |  68       |          210 |
| 75%    | nan                              | nan       |  72       |    17       |  68       |          210 |
| max    | nan                              | nan       |  84       |    20       |  80       |          210 |

## Descriptive Strict Positive Pattern
|   position_1based | top_positive_enriched_residues   | top_logodds                   |   positive_entropy |   background_entropy |
|------------------:|:---------------------------------|:------------------------------|-------------------:|---------------------:|
|                 1 | E,D,M,N,P                        | 1.161,1.161,1.161,1.161,1.161 |            3.9284  |              3.91258 |
|                 2 | C,D,E,G,H                        | 1.161,1.161,1.161,1.161,1.161 |            3.64796 |              3.36036 |
|                 3 | C,E,D,V,R                        | 1.161,1.161,0.782,0.650,0.650 |            3.70077 |              4.06355 |
|                 4 | W,G,V,E,H                        | 1.161,0.709,0.650,0.525,0.314 |            3.83922 |              4.04997 |
|                 5 | M,Y,W,G,E                        | 1.161,1.161,1.161,0.960,0.825 |            3.90143 |              4.0487  |
|                 6 | D,Q,E,H,L                        | 1.161,1.161,0.709,0.573,0.498 |            3.71569 |              3.98171 |
|                 7 | R,W,L,I,Y                        | 1.161,1.161,0.515,0.373,0.373 |            3.94401 |              4.10568 |
|                 8 | N,T,W,H,V                        | 0.825,0.650,0.650,0.650,0.573 |            3.93406 |              4.01889 |
|                 9 | C,D,Q,N,W                        | 1.161,1.161,1.161,1.161,1.161 |            3.8733  |              3.64511 |
|                10 | A,C,D,E,G                        | 0.431,0.431,0.431,0.431,0.431 |            4.07103 |              3.89354 |
|                11 | A,C,D,E,F                        | 0.000,0.000,0.000,0.000,0.000 |            4.32193 |              4.32193 |
|                12 | A,C,D,E,F                        | 0.000,0.000,0.000,0.000,0.000 |            4.32193 |              4.32193 |
|                13 | A,C,D,E,F                        | 0.000,0.000,0.000,0.000,0.000 |            4.32193 |              4.32193 |
|                14 | A,C,D,E,F                        | 0.000,0.000,0.000,0.000,0.000 |            4.32193 |              4.32193 |
|                15 | A,C,D,E,F                        | 0.000,0.000,0.000,0.000,0.000 |            4.32193 |              4.32193 |

Interpretation: this is a focal-style reweighted logistic ranker over positive-pattern and decoy features; public predictor scores are not used.
