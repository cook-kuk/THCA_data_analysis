# GSE230424 Path2Space-style H&E-to-thyroid spatial-axis screen

## Verdict

- Dataset: 4 Visium thyroid slides from GSE230424 with shipped H&E JPEGs, matrices, barcodes, features, and tissue positions.
- GEO sample records expose P1-P4 only; Supplementary Table S5 recovers sample labels: P1/P2 = PTC+HT, P3/P4 = HT.
- Primary tests remain label-free LOSO and sample-centered because disease groups are only n=2 slides per group.
- Spots modeled: 15,489; samples: 4.
- Top H&E-predictable axis: **DM1_low_RAI_score** with pooled rho **0.626**, sample-centered rho **0.644**, domain-centered rho **0.910**.
- QC-only is stronger for this axis (sample-centered rho **0.732**; coord+QC **0.723**), so the conservative claim is morphology/QC-aligned tissue state plus a residual H&E component.
- Coord+QC residual alignment for top axis: rho **0.283**; within-sample permutation p = **0.0010**.
- DM1/low-RAI axis sample-centered rho: **0.644**; AP/TLS composite sample-centered rho: **0.048**.

## Interpretation Boundary

This is useful as an external spatial thyroid control for Paper 2's image-to-spatial-RNA direction. The strongest raw prediction tracks QC/tissue-density structure, but a smaller H&E-aligned residual remains after coordinate+QC adjustment. Disease labels were recovered from Supplementary Table S5, but with only two PTC+HT and two HT slides this should not be framed as a robust disease-group validation. It is not Paper 1 causal mechanism evidence.

## Model Comparison

| target                     | model            |     n |   pooled_rho |   pooled_p |   sample_centered_rho |   sample_centered_p |   median_sample_rho |   samples_rho_gt_0_2 |   n_features |
|:---------------------------|:-----------------|------:|-------------:|-----------:|----------------------:|--------------------:|--------------------:|---------------------:|-------------:|
| DM1_low_RAI_score          | HE_tile_features | 15489 |       0.6257 |     0      |                0.6444 |                   0 |              0.586  |                    4 |           52 |
| DM1_low_RAI_score          | Coord_only       | 15489 |       0.3213 |     0      |                0.3064 |                   0 |              0.4068 |                    3 |            2 |
| DM1_low_RAI_score          | QC_only          | 15489 |       0.7375 |     0      |                0.7321 |                   0 |              0.6629 |                    4 |            3 |
| DM1_low_RAI_score          | Coord_QC         | 15489 |       0.722  |     0      |                0.7227 |                   0 |              0.7191 |                    4 |            5 |
| RAI8_lineage_score         | HE_tile_features | 15489 |       0.6257 |     0      |                0.6444 |                   0 |              0.586  |                    4 |           52 |
| RAI8_lineage_score         | Coord_only       | 15489 |       0.3213 |     0      |                0.3064 |                   0 |              0.4068 |                    3 |            2 |
| RAI8_lineage_score         | QC_only          | 15489 |       0.7375 |     0      |                0.7321 |                   0 |              0.6629 |                    4 |            3 |
| RAI8_lineage_score         | Coord_QC         | 15489 |       0.722  |     0      |                0.7227 |                   0 |              0.7191 |                    4 |            5 |
| MAPK_output_score          | HE_tile_features | 15489 |       0.4346 |     0      |                0.2963 |                   0 |              0.2516 |                    2 |           52 |
| MAPK_output_score          | Coord_only       | 15489 |      -0.1582 |     0      |               -0.128  |                   0 |             -0.1845 |                    1 |            2 |
| MAPK_output_score          | QC_only          | 15489 |       0.6413 |     0      |                0.5355 |                   0 |              0.5299 |                    4 |            3 |
| MAPK_output_score          | Coord_QC         | 15489 |       0.5969 |     0      |                0.5027 |                   0 |              0.5086 |                    4 |            5 |
| AP_TLS_composite_score     | HE_tile_features | 15489 |      -0.0607 |     0      |                0.0479 |                   0 |              0.1369 |                    1 |           52 |
| AP_TLS_composite_score     | Coord_only       | 15489 |       0.3429 |     0      |                0.3502 |                   0 |              0.348  |                    3 |            2 |
| AP_TLS_composite_score     | QC_only          | 15489 |       0.5159 |     0      |                0.6301 |                   0 |              0.6075 |                    3 |            3 |
| AP_TLS_composite_score     | Coord_QC         | 15489 |       0.5611 |     0      |                0.6329 |                   0 |              0.6673 |                    4 |            5 |
| HLA_II_AP_score            | HE_tile_features | 15489 |       0.2276 |     0      |                0.199  |                   0 |              0.2085 |                    2 |           52 |
| HLA_II_AP_score            | Coord_only       | 15489 |       0.3418 |     0      |                0.3726 |                   0 |              0.4252 |                    3 |            2 |
| HLA_II_AP_score            | QC_only          | 15489 |       0.5533 |     0      |                0.6163 |                   0 |              0.6182 |                    3 |            3 |
| HLA_II_AP_score            | Coord_QC         | 15489 |       0.6124 |     0      |                0.6549 |                   0 |              0.6864 |                    4 |            5 |
| B_TLS_score                | HE_tile_features | 15489 |      -0.384  |     0      |               -0.175  |                   0 |              0.0513 |                    0 |           52 |
| B_TLS_score                | Coord_only       | 15489 |       0.2389 |     0      |                0.2633 |                   0 |              0.267  |                    2 |            2 |
| B_TLS_score                | QC_only          | 15489 |       0.2516 |     0      |                0.4098 |                   0 |              0.4661 |                    3 |            3 |
| B_TLS_score                | Coord_QC         | 15489 |       0.246  |     0      |                0.369  |                   0 |              0.4065 |                    3 |            5 |
| CD36_SPP1_macrophage_score | HE_tile_features | 15489 |       0.0292 |     0.0003 |                0.2085 |                   0 |              0.2209 |                    3 |           52 |
| CD36_SPP1_macrophage_score | Coord_only       | 15489 |       0.1889 |     0      |                0.2162 |                   0 |              0.2067 |                    3 |            2 |
| CD36_SPP1_macrophage_score | QC_only          | 15489 |       0.5905 |     0      |                0.6554 |                   0 |              0.6039 |                    4 |            3 |
| CD36_SPP1_macrophage_score | Coord_QC         | 15489 |       0.6012 |     0      |                0.667  |                   0 |              0.6294 |                    4 |            5 |
| Tumor_ZCCHC12_score        | HE_tile_features | 15489 |       0.2253 |     0      |                0.2327 |                   0 |              0.1926 |                    2 |           52 |
| Tumor_ZCCHC12_score        | Coord_only       | 15489 |      -0.2504 |     0      |               -0.1427 |                   0 |             -0.1361 |                    0 |            2 |
| Tumor_ZCCHC12_score        | QC_only          | 15489 |       0.6372 |     0      |                0.511  |                   0 |              0.4642 |                    4 |            3 |
| Tumor_ZCCHC12_score        | Coord_QC         | 15489 |       0.5999 |     0      |                0.4988 |                   0 |              0.4847 |                    4 |            5 |

## Residual Alignment

| target                     | within_sample_adjustment   |   residual_alignment_rho |      p |     n |
|:---------------------------|:---------------------------|-------------------------:|-------:|------:|
| DM1_low_RAI_score          | sample_mean_only           |                   0.6444 | 0      | 15489 |
| DM1_low_RAI_score          | coord_only                 |                   0.567  | 0      | 15489 |
| DM1_low_RAI_score          | qc_only                    |                   0.3316 | 0      | 15489 |
| DM1_low_RAI_score          | coord_qc                   |                   0.2829 | 0      | 15489 |
| RAI8_lineage_score         | sample_mean_only           |                   0.6444 | 0      | 15489 |
| RAI8_lineage_score         | coord_only                 |                   0.567  | 0      | 15489 |
| RAI8_lineage_score         | qc_only                    |                   0.3316 | 0      | 15489 |
| RAI8_lineage_score         | coord_qc                   |                   0.2829 | 0      | 15489 |
| MAPK_output_score          | sample_mean_only           |                   0.2963 | 0      | 15489 |
| MAPK_output_score          | coord_only                 |                   0.2629 | 0      | 15489 |
| MAPK_output_score          | qc_only                    |                   0.055  | 0      | 15489 |
| MAPK_output_score          | coord_qc                   |                   0.0687 | 0      | 15489 |
| AP_TLS_composite_score     | sample_mean_only           |                   0.0479 | 0      | 15489 |
| AP_TLS_composite_score     | coord_only                 |                   0.0336 | 0      | 15489 |
| AP_TLS_composite_score     | qc_only                    |                   0.056  | 0      | 15489 |
| AP_TLS_composite_score     | coord_qc                   |                   0.0552 | 0      | 15489 |
| HLA_II_AP_score            | sample_mean_only           |                   0.199  | 0      | 15489 |
| HLA_II_AP_score            | coord_only                 |                   0.173  | 0      | 15489 |
| HLA_II_AP_score            | qc_only                    |                   0.0722 | 0      | 15489 |
| HLA_II_AP_score            | coord_qc                   |                   0.0688 | 0      | 15489 |
| B_TLS_score                | sample_mean_only           |                  -0.175  | 0      | 15489 |
| B_TLS_score                | coord_only                 |                  -0.1728 | 0      | 15489 |
| B_TLS_score                | qc_only                    |                   0.0294 | 0.0003 | 15489 |
| B_TLS_score                | coord_qc                   |                   0.0258 | 0.0013 | 15489 |
| CD36_SPP1_macrophage_score | sample_mean_only           |                   0.2085 | 0      | 15489 |
| CD36_SPP1_macrophage_score | coord_only                 |                   0.1783 | 0      | 15489 |
| CD36_SPP1_macrophage_score | qc_only                    |                   0.0161 | 0.0453 | 15489 |
| CD36_SPP1_macrophage_score | coord_qc                   |                   0.0336 | 0      | 15489 |
| Tumor_ZCCHC12_score        | sample_mean_only           |                   0.2327 | 0      | 15489 |
| Tumor_ZCCHC12_score        | coord_only                 |                   0.197  | 0      | 15489 |
| Tumor_ZCCHC12_score        | qc_only                    |                   0.0252 | 0.0017 | 15489 |
| Tumor_ZCCHC12_score        | coord_qc                   |                   0.0474 | 0      | 15489 |

## Domain Aggregation

| target                     |   n_domains |   pooled_domain_rho |   sample_centered_domain_rho |      p |
|:---------------------------|------------:|--------------------:|-----------------------------:|-------:|
| DM1_low_RAI_score          |          39 |              0.7824 |                       0.9099 | 0      |
| RAI8_lineage_score         |          39 |              0.7824 |                       0.9099 | 0      |
| MAPK_output_score          |          39 |              0.6899 |                       0.504  | 0      |
| AP_TLS_composite_score     |          39 |             -0.2504 |                      -0.0526 | 0.1242 |
| HLA_II_AP_score            |          39 |              0.4136 |                       0.2814 | 0.0089 |
| B_TLS_score                |          39 |             -0.6472 |                      -0.32   | 0      |
| CD36_SPP1_macrophage_score |          39 |             -0.1715 |                       0.4134 | 0.2966 |
| Tumor_ZCCHC12_score        |          39 |              0.2263 |                       0.3411 | 0.1659 |

## Permutation Alignment

| target                     |   observed_sample_centered_rho |   n_permutations |   empirical_p |
|:---------------------------|-------------------------------:|-----------------:|--------------:|
| DM1_low_RAI_score          |                         0.6444 |             1000 |         0.001 |
| RAI8_lineage_score         |                         0.6444 |             1000 |         0.001 |
| MAPK_output_score          |                         0.2963 |             1000 |         0.001 |
| AP_TLS_composite_score     |                         0.0479 |             1000 |         0.001 |
| HLA_II_AP_score            |                         0.199  |             1000 |         0.001 |
| B_TLS_score                |                        -0.175  |             1000 |         0.001 |
| CD36_SPP1_macrophage_score |                         0.2085 |             1000 |         0.001 |
| Tumor_ZCCHC12_score        |                         0.2327 |             1000 |         0.001 |

## Sample Context

| sample   | disease_group   |   RAI8_lineage_score |   MAPK_output_score |   HLA_II_AP_score |   B_TLS_score |   T_cell_score |   CD36_SPP1_macrophage_score |   Tumor_ZCCHC12_score |   DM1_low_RAI_score |   AP_TLS_composite_score |   total_counts |   n_genes_by_counts |   pct_counts_mt |
|:---------|:----------------|---------------------:|--------------------:|------------------:|--------------:|---------------:|-----------------------------:|----------------------:|--------------------:|-------------------------:|---------------:|--------------------:|----------------:|
| P1       | PTC+HT          |               0.1418 |              0.0685 |            0.038  |        0.181  |         0.0874 |                       0.0635 |               -0.0112 |             -0.1418 |                   0.1022 |        6006.85 |             2036.02 |          5.3295 |
| P2       | PTC+HT          |              -0.0286 |             -0.0247 |            0.0002 |        0.1493 |         0.0175 |                       0.1805 |                0.1086 |              0.0286 |                   0.0557 |        4539.17 |             1662.24 |          2.5275 |
| P3       | HT              |              -0.3473 |             -0.1586 |           -0.179  |       -0.0423 |         0.0037 |                      -0.1694 |               -0.2497 |              0.3473 |                  -0.0725 |        3210.83 |             1166.36 |          4.1937 |
| P4       | HT              |               0.2128 |              0.1022 |            0.1271 |       -0.2074 |        -0.0842 |                      -0.0243 |                0.154  |             -0.2128 |                  -0.0548 |        4306.64 |             1769.11 |          2.7499 |
