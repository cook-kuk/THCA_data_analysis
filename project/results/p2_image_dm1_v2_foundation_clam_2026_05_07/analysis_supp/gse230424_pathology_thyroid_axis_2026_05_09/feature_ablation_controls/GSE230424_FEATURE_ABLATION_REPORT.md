# GSE230424 feature-family ablation controls

## Verdict

- DM1/low-RAI raw target all-feature rho: **0.644**.
- Best raw feature family: **radius96_only** rho **0.645**.
- DM1/low-RAI coord+QC residual target all-feature rho: **0.232**.
- Best residual feature family: **radius96_only** rho **0.246**.
- Density/texture residual rho: 0.174; RGB-distribution residual rho: 0.143.

## Interpretation

The DM1/RAI signal is not a single-feature artifact, but density/texture and stain/color summaries carry much of the raw signal. After coordinate+QC residualization, the retained signal is smaller and feature-family dependent; this supports a cautious Paper 2 image-to-spatial-RNA claim with explicit QC/stain caveats.

## Full Results

| target              | target_mode              | feature_set      |     n |   n_features |   pooled_rho |   sample_centered_rho |   median_sample_rho |   samples_rho_gt_0_1 |
|:--------------------|:-------------------------|:-----------------|------:|-------------:|-------------:|----------------------:|--------------------:|---------------------:|
| DM1_low_RAI_score   | raw_smoothed             | all_HE_features  | 15489 |           52 |       0.6257 |                0.6444 |              0.586  |                    4 |
| DM1_low_RAI_score   | raw_smoothed             | radius48_only    | 15489 |           26 |       0.5532 |                0.5296 |              0.46   |                    4 |
| DM1_low_RAI_score   | raw_smoothed             | radius96_only    | 15489 |           26 |       0.6263 |                0.6451 |              0.5886 |                    4 |
| DM1_low_RAI_score   | raw_smoothed             | rgb_distribution | 15489 |           24 |       0.4826 |                0.417  |              0.4105 |                    3 |
| DM1_low_RAI_score   | raw_smoothed             | hsv_gray_summary | 15489 |           16 |       0.6179 |                0.6178 |              0.5259 |                    4 |
| DM1_low_RAI_score   | raw_smoothed             | density_texture  | 15489 |           10 |       0.5185 |                0.5151 |              0.4962 |                    3 |
| DM1_low_RAI_score   | raw_smoothed             | stain_proxy_only | 15489 |            4 |       0.4348 |                0.2834 |              0.3347 |                    3 |
| DM1_low_RAI_score   | raw_smoothed             | color_means_only | 15489 |           14 |       0.4885 |                0.4065 |              0.4202 |                    3 |
| DM1_low_RAI_score   | coord_qc_residual_target | all_HE_features  | 15489 |           52 |       0.2162 |                0.2323 |              0.1947 |                    4 |
| DM1_low_RAI_score   | coord_qc_residual_target | radius48_only    | 15489 |           26 |       0.2038 |                0.2159 |              0.1933 |                    4 |
| DM1_low_RAI_score   | coord_qc_residual_target | radius96_only    | 15489 |           26 |       0.2293 |                0.2459 |              0.2104 |                    4 |
| DM1_low_RAI_score   | coord_qc_residual_target | rgb_distribution | 15489 |           24 |       0.143  |                0.1426 |              0.1218 |                    2 |
| DM1_low_RAI_score   | coord_qc_residual_target | hsv_gray_summary | 15489 |           16 |       0.2094 |                0.2388 |              0.213  |                    4 |
| DM1_low_RAI_score   | coord_qc_residual_target | density_texture  | 15489 |           10 |       0.1771 |                0.1737 |              0.1653 |                    3 |
| DM1_low_RAI_score   | coord_qc_residual_target | stain_proxy_only | 15489 |            4 |       0.1575 |                0.1671 |              0.1493 |                    3 |
| DM1_low_RAI_score   | coord_qc_residual_target | color_means_only | 15489 |           14 |       0.1773 |                0.1967 |              0.1656 |                    4 |
| MAPK_output_score   | raw_smoothed             | all_HE_features  | 15489 |           52 |       0.4346 |                0.2963 |              0.2516 |                    3 |
| MAPK_output_score   | raw_smoothed             | radius48_only    | 15489 |           26 |       0.3943 |                0.2262 |              0.1859 |                    2 |
| MAPK_output_score   | raw_smoothed             | radius96_only    | 15489 |           26 |       0.4342 |                0.2853 |              0.243  |                    3 |
| MAPK_output_score   | raw_smoothed             | rgb_distribution | 15489 |           24 |       0.389  |                0.2147 |              0.1542 |                    3 |
| MAPK_output_score   | raw_smoothed             | hsv_gray_summary | 15489 |           16 |       0.4341 |                0.2805 |              0.2154 |                    3 |
| MAPK_output_score   | raw_smoothed             | density_texture  | 15489 |           10 |       0.286  |                0.1566 |              0.1009 |                    2 |
| MAPK_output_score   | raw_smoothed             | stain_proxy_only | 15489 |            4 |       0.2207 |                0.0313 |              0.0892 |                    2 |
| MAPK_output_score   | raw_smoothed             | color_means_only | 15489 |           14 |       0.3283 |                0.1191 |              0.1211 |                    2 |
| MAPK_output_score   | coord_qc_residual_target | all_HE_features  | 15489 |           52 |       0.0407 |                0.0607 |              0.0424 |                    1 |
| MAPK_output_score   | coord_qc_residual_target | radius48_only    | 15489 |           26 |       0.02   |                0.0375 |              0.0267 |                    1 |
| MAPK_output_score   | coord_qc_residual_target | radius96_only    | 15489 |           26 |       0.0323 |                0.0567 |              0.0316 |                    1 |
| MAPK_output_score   | coord_qc_residual_target | rgb_distribution | 15489 |           24 |       0.0548 |                0.0699 |              0.0718 |                    1 |
| MAPK_output_score   | coord_qc_residual_target | hsv_gray_summary | 15489 |           16 |       0.0271 |                0.0576 |              0.0556 |                    1 |
| MAPK_output_score   | coord_qc_residual_target | density_texture  | 15489 |           10 |      -0.0074 |               -0.004  |              0.0188 |                    0 |
| MAPK_output_score   | coord_qc_residual_target | stain_proxy_only | 15489 |            4 |       0.0156 |                0.0266 |              0.0396 |                    0 |
| MAPK_output_score   | coord_qc_residual_target | color_means_only | 15489 |           14 |       0.0099 |                0.0301 |              0.0153 |                    0 |
| HLA_II_AP_score     | raw_smoothed             | all_HE_features  | 15489 |           52 |       0.2276 |                0.199  |              0.2085 |                    4 |
| HLA_II_AP_score     | raw_smoothed             | radius48_only    | 15489 |           26 |       0.1135 |                0.0616 |              0.0556 |                    2 |
| HLA_II_AP_score     | raw_smoothed             | radius96_only    | 15489 |           26 |       0.2359 |                0.2084 |              0.2175 |                    4 |
| HLA_II_AP_score     | raw_smoothed             | rgb_distribution | 15489 |           24 |       0.0629 |                0.0264 |              0.0743 |                    2 |
| HLA_II_AP_score     | raw_smoothed             | hsv_gray_summary | 15489 |           16 |       0.1237 |                0.0876 |              0.1094 |                    2 |
| HLA_II_AP_score     | raw_smoothed             | density_texture  | 15489 |           10 |      -0.0086 |               -0.038  |             -0.0323 |                    2 |
| HLA_II_AP_score     | raw_smoothed             | stain_proxy_only | 15489 |            4 |      -0.3827 |               -0.44   |             -0.4579 |                    1 |
| HLA_II_AP_score     | raw_smoothed             | color_means_only | 15489 |           14 |      -0.2442 |               -0.3018 |             -0.409  |                    1 |
| HLA_II_AP_score     | coord_qc_residual_target | all_HE_features  | 15489 |           52 |       0.1861 |                0.1861 |              0.2239 |                    3 |
| HLA_II_AP_score     | coord_qc_residual_target | radius48_only    | 15489 |           26 |       0.1466 |                0.1452 |              0.1587 |                    3 |
| HLA_II_AP_score     | coord_qc_residual_target | radius96_only    | 15489 |           26 |       0.1926 |                0.1947 |              0.2262 |                    3 |
| HLA_II_AP_score     | coord_qc_residual_target | rgb_distribution | 15489 |           24 |       0.1824 |                0.1918 |              0.231  |                    3 |
| HLA_II_AP_score     | coord_qc_residual_target | hsv_gray_summary | 15489 |           16 |       0.1571 |                0.1649 |              0.2079 |                    3 |
| HLA_II_AP_score     | coord_qc_residual_target | density_texture  | 15489 |           10 |       0.0739 |                0.0642 |              0.0771 |                    1 |
| HLA_II_AP_score     | coord_qc_residual_target | stain_proxy_only | 15489 |            4 |       0.1443 |                0.149  |              0.1384 |                    2 |
| HLA_II_AP_score     | coord_qc_residual_target | color_means_only | 15489 |           14 |       0.151  |                0.1393 |              0.1747 |                    3 |
| Tumor_ZCCHC12_score | raw_smoothed             | all_HE_features  | 15489 |           52 |       0.2253 |                0.2327 |              0.1926 |                    2 |
| Tumor_ZCCHC12_score | raw_smoothed             | radius48_only    | 15489 |           26 |       0.237  |                0.179  |              0.1421 |                    2 |
| Tumor_ZCCHC12_score | raw_smoothed             | radius96_only    | 15489 |           26 |       0.2373 |                0.2355 |              0.1993 |                    3 |
| Tumor_ZCCHC12_score | raw_smoothed             | rgb_distribution | 15489 |           24 |       0.0472 |                0.101  |              0.1097 |                    2 |
| Tumor_ZCCHC12_score | raw_smoothed             | hsv_gray_summary | 15489 |           16 |       0.2599 |                0.2577 |              0.2542 |                    4 |
| Tumor_ZCCHC12_score | raw_smoothed             | density_texture  | 15489 |           10 |       0.1792 |                0.1842 |              0.1098 |                    2 |
| Tumor_ZCCHC12_score | raw_smoothed             | stain_proxy_only | 15489 |            4 |      -0.1189 |               -0.1208 |             -0.0947 |                    1 |
| Tumor_ZCCHC12_score | raw_smoothed             | color_means_only | 15489 |           14 |       0.0307 |               -0.0805 |             -0.1247 |                    1 |
| Tumor_ZCCHC12_score | coord_qc_residual_target | all_HE_features  | 15489 |           52 |       0.0734 |                0.0982 |              0.0941 |                    2 |
| Tumor_ZCCHC12_score | coord_qc_residual_target | radius48_only    | 15489 |           26 |       0.082  |                0.1069 |              0.1008 |                    2 |
| Tumor_ZCCHC12_score | coord_qc_residual_target | radius96_only    | 15489 |           26 |       0.0752 |                0.1056 |              0.102  |                    2 |
| Tumor_ZCCHC12_score | coord_qc_residual_target | rgb_distribution | 15489 |           24 |       0.0486 |                0.0595 |              0.0627 |                    0 |
| Tumor_ZCCHC12_score | coord_qc_residual_target | hsv_gray_summary | 15489 |           16 |       0.088  |                0.106  |              0.1053 |                    2 |
| Tumor_ZCCHC12_score | coord_qc_residual_target | density_texture  | 15489 |           10 |       0.0816 |                0.0757 |              0.0953 |                    2 |
| Tumor_ZCCHC12_score | coord_qc_residual_target | stain_proxy_only | 15489 |            4 |       0.0265 |                0.0359 |              0.0518 |                    1 |
| Tumor_ZCCHC12_score | coord_qc_residual_target | color_means_only | 15489 |           14 |       0.0485 |                0.05   |              0.0483 |                    1 |
