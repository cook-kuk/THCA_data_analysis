# GSE230424 residual-target controls

## Verdict

- Spots: 15,489; samples: 4.
- Strict test: residualize observed spatial targets within each sample, then train leave-one-sample-out H&E models to predict the residual targets.
- DM1/low-RAI residual after coordinate+QC adjustment: H&E sample-centered rho **0.232**, median sample rho **0.195**, permutation p **0.0010**.
- MAPK residual target rho 0.061; HLA-II/AP 0.186; Tumor/ZCCHC12 0.098.

## Interpretation

This strengthens the conservative GSE230424 claim: the major DM1/RAI signal is QC/tissue-density aligned, but H&E still predicts a smaller residual component after coordinate+QC removal. This remains Paper 2 image-to-spatial-RNA support, not Paper 1 causal mechanism evidence.

## Model Comparison

| target                     | residual_target_adjustment   | model            |     n |   pooled_rho |   pooled_p |   sample_centered_rho |   sample_centered_p |   median_sample_rho |   samples_rho_gt_0_1 |
|:---------------------------|:-----------------------------|:-----------------|------:|-------------:|-----------:|----------------------:|--------------------:|--------------------:|---------------------:|
| DM1_low_RAI_score          | sample_mean_only             | HE_tile_features | 15489 |       0.6231 |     0      |                0.6894 |              0      |              0.651  |                    4 |
| DM1_low_RAI_score          | sample_mean_only             | Coord_only       | 15489 |       0.3187 |     0      |                0.3109 |              0      |              0.4423 |                    3 |
| DM1_low_RAI_score          | sample_mean_only             | QC_only          | 15489 |       0.6303 |     0      |                0.723  |              0      |              0.6421 |                    4 |
| DM1_low_RAI_score          | sample_mean_only             | Coord_QC         | 15489 |       0.6172 |     0      |                0.7003 |              0      |              0.7089 |                    4 |
| DM1_low_RAI_score          | coord_only                   | HE_tile_features | 15489 |       0.4904 |     0      |                0.5345 |              0      |              0.4566 |                    4 |
| DM1_low_RAI_score          | coord_only                   | Coord_only       | 15489 |      -0.0189 |     0.0189 |               -0.0271 |              0.0007 |             -0.0054 |                    0 |
| DM1_low_RAI_score          | coord_only                   | QC_only          | 15489 |       0.5408 |     0      |                0.6341 |              0      |              0.5275 |                    4 |
| DM1_low_RAI_score          | coord_only                   | Coord_QC         | 15489 |       0.545  |     0      |                0.6314 |              0      |              0.4931 |                    4 |
| DM1_low_RAI_score          | qc_only                      | HE_tile_features | 15489 |       0.3012 |     0      |                0.3295 |              0      |              0.2707 |                    4 |
| DM1_low_RAI_score          | qc_only                      | Coord_only       | 15489 |       0.2038 |     0      |                0.1997 |              0      |              0.3416 |                    3 |
| DM1_low_RAI_score          | qc_only                      | QC_only          | 15489 |      -0.0214 |     0.0078 |               -0.0115 |              0.1539 |             -0.0574 |                    0 |
| DM1_low_RAI_score          | qc_only                      | Coord_QC         | 15489 |       0.2112 |     0      |                0.2089 |              0      |              0.3755 |                    3 |
| DM1_low_RAI_score          | coord_qc                     | HE_tile_features | 15489 |       0.2162 |     0      |                0.2323 |              0      |              0.1947 |                    4 |
| DM1_low_RAI_score          | coord_qc                     | Coord_only       | 15489 |      -0.0053 |     0.5094 |                0.005  |              0.5364 |              0.0041 |                    0 |
| DM1_low_RAI_score          | coord_qc                     | QC_only          | 15489 |      -0.0428 |     0      |               -0.0477 |              0      |             -0.0196 |                    0 |
| DM1_low_RAI_score          | coord_qc                     | Coord_QC         | 15489 |      -0.0305 |     0.0001 |               -0.028  |              0.0005 |             -0.0118 |                    0 |
| RAI8_lineage_score         | sample_mean_only             | HE_tile_features | 15489 |       0.6231 |     0      |                0.6894 |              0      |              0.651  |                    4 |
| RAI8_lineage_score         | sample_mean_only             | Coord_only       | 15489 |       0.3187 |     0      |                0.3109 |              0      |              0.4423 |                    3 |
| RAI8_lineage_score         | sample_mean_only             | QC_only          | 15489 |       0.6303 |     0      |                0.723  |              0      |              0.6421 |                    4 |
| RAI8_lineage_score         | sample_mean_only             | Coord_QC         | 15489 |       0.6172 |     0      |                0.7003 |              0      |              0.7089 |                    4 |
| RAI8_lineage_score         | coord_only                   | HE_tile_features | 15489 |       0.4904 |     0      |                0.5345 |              0      |              0.4566 |                    4 |
| RAI8_lineage_score         | coord_only                   | Coord_only       | 15489 |      -0.0189 |     0.0189 |               -0.0271 |              0.0007 |             -0.0054 |                    0 |
| RAI8_lineage_score         | coord_only                   | QC_only          | 15489 |       0.5408 |     0      |                0.6341 |              0      |              0.5275 |                    4 |
| RAI8_lineage_score         | coord_only                   | Coord_QC         | 15489 |       0.545  |     0      |                0.6314 |              0      |              0.4931 |                    4 |
| RAI8_lineage_score         | qc_only                      | HE_tile_features | 15489 |       0.3012 |     0      |                0.3295 |              0      |              0.2707 |                    4 |
| RAI8_lineage_score         | qc_only                      | Coord_only       | 15489 |       0.2038 |     0      |                0.1997 |              0      |              0.3416 |                    3 |
| RAI8_lineage_score         | qc_only                      | QC_only          | 15489 |      -0.0214 |     0.0078 |               -0.0115 |              0.1539 |             -0.0574 |                    0 |
| RAI8_lineage_score         | qc_only                      | Coord_QC         | 15489 |       0.2112 |     0      |                0.2089 |              0      |              0.3755 |                    3 |
| RAI8_lineage_score         | coord_qc                     | HE_tile_features | 15489 |       0.2162 |     0      |                0.2323 |              0      |              0.1947 |                    4 |
| RAI8_lineage_score         | coord_qc                     | Coord_only       | 15489 |      -0.0053 |     0.5094 |                0.005  |              0.5364 |              0.0041 |                    0 |
| RAI8_lineage_score         | coord_qc                     | QC_only          | 15489 |      -0.0428 |     0      |               -0.0477 |              0      |             -0.0196 |                    0 |
| RAI8_lineage_score         | coord_qc                     | Coord_QC         | 15489 |      -0.0305 |     0.0001 |               -0.028  |              0.0005 |             -0.0118 |                    0 |
| MAPK_output_score          | sample_mean_only             | HE_tile_features | 15489 |       0.2853 |     0      |                0.3256 |              0      |              0.3198 |                    3 |
| MAPK_output_score          | sample_mean_only             | Coord_only       | 15489 |      -0.2254 |     0      |               -0.2075 |              0      |             -0.216  |                    0 |
| MAPK_output_score          | sample_mean_only             | QC_only          | 15489 |       0.3885 |     0      |                0.5126 |              0      |              0.4967 |                    4 |
| MAPK_output_score          | sample_mean_only             | Coord_QC         | 15489 |       0.3426 |     0      |                0.4504 |              0      |              0.4634 |                    4 |
| MAPK_output_score          | coord_only                   | HE_tile_features | 15489 |       0.2584 |     0      |                0.2946 |              0      |              0.2581 |                    3 |
| MAPK_output_score          | coord_only                   | Coord_only       | 15489 |       0.013  |     0.1057 |                0.0252 |              0.0017 |              0.0216 |                    0 |
| MAPK_output_score          | coord_only                   | QC_only          | 15489 |       0.3582 |     0      |                0.4739 |              0      |              0.4277 |                    4 |
| MAPK_output_score          | coord_only                   | Coord_QC         | 15489 |       0.3636 |     0      |                0.4619 |              0      |              0.433  |                    4 |
| MAPK_output_score          | qc_only                      | HE_tile_features | 15489 |      -0.0121 |     0.1324 |                0.033  |              0      |              0.0289 |                    1 |
| MAPK_output_score          | qc_only                      | Coord_only       | 15489 |      -0.0575 |     0      |               -0.0598 |              0      |             -0.086  |                    0 |
| MAPK_output_score          | qc_only                      | QC_only          | 15489 |      -0.0006 |     0.94   |               -0.0482 |              0      |              0.0048 |                    0 |
| MAPK_output_score          | qc_only                      | Coord_QC         | 15489 |      -0.0537 |     0      |               -0.0544 |              0      |             -0.0708 |                    0 |
| MAPK_output_score          | coord_qc                     | HE_tile_features | 15489 |       0.0407 |     0      |                0.0607 |              0      |              0.0424 |                    1 |
| MAPK_output_score          | coord_qc                     | Coord_only       | 15489 |       0.0195 |     0.0154 |                0.0184 |              0.0223 |              0.0242 |                    1 |
| MAPK_output_score          | coord_qc                     | QC_only          | 15489 |      -0.0102 |     0.2038 |               -0.0273 |              0.0007 |             -0.004  |                    0 |
| MAPK_output_score          | coord_qc                     | Coord_QC         | 15489 |       0.0009 |     0.9123 |               -0.0093 |              0.2446 |             -0.0142 |                    0 |
| HLA_II_AP_score            | sample_mean_only             | HE_tile_features | 15489 |       0.1417 |     0      |                0.1842 |              0      |              0.1838 |                    2 |
| HLA_II_AP_score            | sample_mean_only             | Coord_only       | 15489 |       0.3686 |     0      |                0.3728 |              0      |              0.4097 |                    4 |
| HLA_II_AP_score            | sample_mean_only             | QC_only          | 15489 |       0.4589 |     0      |                0.6147 |              0      |              0.6048 |                    4 |
| HLA_II_AP_score            | sample_mean_only             | Coord_QC         | 15489 |       0.5583 |     0      |                0.6568 |              0      |              0.6772 |                    4 |
| HLA_II_AP_score            | coord_only                   | HE_tile_features | 15489 |       0.0103 |     0.1998 |                0.0259 |              0.0012 |              0.111  |                    2 |
| HLA_II_AP_score            | coord_only                   | Coord_only       | 15489 |       0.0319 |     0.0001 |                0.0276 |              0.0006 |             -0.0084 |                    0 |
| HLA_II_AP_score            | coord_only                   | QC_only          | 15489 |       0.4425 |     0      |                0.5836 |              0      |              0.544  |                    4 |
| HLA_II_AP_score            | coord_only                   | Coord_QC         | 15489 |       0.4477 |     0      |                0.5722 |              0      |              0.5757 |                    4 |
| HLA_II_AP_score            | qc_only                      | HE_tile_features | 15489 |       0.2391 |     0      |                0.2471 |              0      |              0.2751 |                    4 |
| HLA_II_AP_score            | qc_only                      | Coord_only       | 15489 |       0.2412 |     0      |                0.2464 |              0      |              0.3115 |                    3 |
| HLA_II_AP_score            | qc_only                      | QC_only          | 15489 |       0.0238 |     0.0031 |               -0.0141 |              0.0785 |             -0.0082 |                    0 |
| HLA_II_AP_score            | qc_only                      | Coord_QC         | 15489 |       0.2238 |     0      |                0.2397 |              0      |              0.3041 |                    3 |
| HLA_II_AP_score            | coord_qc                     | HE_tile_features | 15489 |       0.1861 |     0      |                0.1861 |              0      |              0.2239 |                    3 |
| HLA_II_AP_score            | coord_qc                     | Coord_only       | 15489 |      -0.024  |     0.0028 |               -0.0304 |              0.0002 |             -0.016  |                    0 |
| HLA_II_AP_score            | coord_qc                     | QC_only          | 15489 |       0.0137 |     0.0871 |                0.0161 |              0.0448 |              0.0262 |                    0 |
| HLA_II_AP_score            | coord_qc                     | Coord_QC         | 15489 |       0.0009 |     0.9147 |               -0.0083 |              0.2988 |              0.0034 |                    0 |
| CD36_SPP1_macrophage_score | sample_mean_only             | HE_tile_features | 15489 |       0.2054 |     0      |                0.3171 |              0      |              0.2784 |                    3 |
| CD36_SPP1_macrophage_score | sample_mean_only             | Coord_only       | 15489 |       0.2357 |     0      |                0.2444 |              0      |              0.2391 |                    4 |
| CD36_SPP1_macrophage_score | sample_mean_only             | QC_only          | 15489 |       0.5268 |     0      |                0.6498 |              0      |              0.5715 |                    4 |
| CD36_SPP1_macrophage_score | sample_mean_only             | Coord_QC         | 15489 |       0.5928 |     0      |                0.6797 |              0      |              0.6297 |                    4 |
| CD36_SPP1_macrophage_score | coord_only                   | HE_tile_features | 15489 |       0.1732 |     0      |                0.2537 |              0      |              0.2114 |                    3 |
| CD36_SPP1_macrophage_score | coord_only                   | Coord_only       | 15489 |       0.0089 |     0.2676 |               -0.0062 |              0.4396 |             -0.0148 |                    0 |
| CD36_SPP1_macrophage_score | coord_only                   | QC_only          | 15489 |       0.5221 |     0      |                0.6311 |              0      |              0.5403 |                    4 |
| CD36_SPP1_macrophage_score | coord_only                   | Coord_QC         | 15489 |       0.5145 |     0      |                0.6146 |              0      |              0.5581 |                    4 |
| CD36_SPP1_macrophage_score | qc_only                      | HE_tile_features | 15489 |       0.104  |     0      |                0.1179 |              0      |              0.1277 |                    3 |
| CD36_SPP1_macrophage_score | qc_only                      | Coord_only       | 15489 |       0.1393 |     0      |                0.1419 |              0      |              0.1837 |                    3 |
| CD36_SPP1_macrophage_score | qc_only                      | QC_only          | 15489 |       0.0184 |     0.022  |                0.0095 |              0.2364 |              0.0184 |                    0 |
| CD36_SPP1_macrophage_score | qc_only                      | Coord_QC         | 15489 |       0.1349 |     0      |                0.1364 |              0      |              0.1743 |                    3 |
| CD36_SPP1_macrophage_score | coord_qc                     | HE_tile_features | 15489 |       0.0692 |     0      |                0.0847 |              0      |              0.0835 |                    1 |
| CD36_SPP1_macrophage_score | coord_qc                     | Coord_only       | 15489 |      -0.0037 |     0.6487 |               -0.0072 |              0.3669 |              0.0041 |                    0 |
| CD36_SPP1_macrophage_score | coord_qc                     | QC_only          | 15489 |       0.0102 |     0.2052 |                0.0355 |              0      |              0.0251 |                    0 |
| CD36_SPP1_macrophage_score | coord_qc                     | Coord_QC         | 15489 |       0.0026 |     0.7498 |                0.0161 |              0.0447 |              0.0294 |                    0 |
| Tumor_ZCCHC12_score        | sample_mean_only             | HE_tile_features | 15489 |       0.3147 |     0      |                0.3035 |              0      |              0.2922 |                    4 |
| Tumor_ZCCHC12_score        | sample_mean_only             | Coord_only       | 15489 |      -0.08   |     0      |               -0.0938 |              0      |             -0.135  |                    0 |
| Tumor_ZCCHC12_score        | sample_mean_only             | QC_only          | 15489 |       0.3863 |     0      |                0.5336 |              0      |              0.4586 |                    4 |
| Tumor_ZCCHC12_score        | sample_mean_only             | Coord_QC         | 15489 |       0.4043 |     0      |                0.5066 |              0      |              0.4817 |                    4 |
| Tumor_ZCCHC12_score        | coord_only                   | HE_tile_features | 15489 |       0.3517 |     0      |                0.3449 |              0      |              0.3133 |                    4 |
| Tumor_ZCCHC12_score        | coord_only                   | Coord_only       | 15489 |       0.0034 |     0.6738 |                0.006  |              0.4525 |              0.0322 |                    0 |
| Tumor_ZCCHC12_score        | coord_only                   | QC_only          | 15489 |       0.3977 |     0      |                0.5111 |              0      |              0.4635 |                    4 |
| Tumor_ZCCHC12_score        | coord_only                   | Coord_QC         | 15489 |       0.3898 |     0      |                0.486  |              0      |              0.4184 |                    4 |
| Tumor_ZCCHC12_score        | qc_only                      | HE_tile_features | 15489 |      -0.0401 |     0      |                0.0319 |              0.0001 |              0.0313 |                    0 |
| Tumor_ZCCHC12_score        | qc_only                      | Coord_only       | 15489 |       0.1412 |     0      |                0.1133 |              0      |              0.1218 |                    2 |
| Tumor_ZCCHC12_score        | qc_only                      | QC_only          | 15489 |       0.01   |     0.2142 |               -0.049  |              0      |             -0.004  |                    0 |
| Tumor_ZCCHC12_score        | qc_only                      | Coord_QC         | 15489 |       0.1348 |     0      |                0.1158 |              0      |              0.1369 |                    2 |
| Tumor_ZCCHC12_score        | coord_qc                     | HE_tile_features | 15489 |       0.0734 |     0      |                0.0982 |              0      |              0.0941 |                    2 |
| Tumor_ZCCHC12_score        | coord_qc                     | Coord_only       | 15489 |      -0.0298 |     0.0002 |               -0.0549 |              0      |             -0.0774 |                    0 |
| Tumor_ZCCHC12_score        | coord_qc                     | QC_only          | 15489 |      -0.016  |     0.0471 |               -0.0453 |              0      |             -0.0219 |                    0 |
| Tumor_ZCCHC12_score        | coord_qc                     | Coord_QC         | 15489 |      -0.0248 |     0.002  |               -0.0464 |              0      |             -0.0445 |                    0 |

## Permutation

| target                     |   coord_qc_residual_target_he_rho |   n_permutations |   empirical_p |
|:---------------------------|----------------------------------:|-----------------:|--------------:|
| DM1_low_RAI_score          |                            0.2162 |             1000 |         0.001 |
| RAI8_lineage_score         |                            0.2162 |             1000 |         0.001 |
| MAPK_output_score          |                            0.0407 |             1000 |         0.001 |
| HLA_II_AP_score            |                            0.1861 |             1000 |         0.001 |
| CD36_SPP1_macrophage_score |                            0.0692 |             1000 |         0.001 |
| Tumor_ZCCHC12_score        |                            0.0734 |             1000 |         0.001 |
