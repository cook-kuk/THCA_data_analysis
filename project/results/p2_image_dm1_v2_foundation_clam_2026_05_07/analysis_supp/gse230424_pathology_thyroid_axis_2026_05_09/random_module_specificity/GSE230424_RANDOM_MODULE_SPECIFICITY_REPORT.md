# GSE230424 random-module specificity controls

## Verdict

- Matched random modules: 250 expression/detection-matched 8-gene modules.
- Actual DM1/RAI raw H&E rho: **0.644**; random-module percentile **99.2%**; empirical upper-tail p = **0.0120**.
- Actual DM1/RAI coord+QC residual-target H&E rho: **0.232**; random-module percentile **100.0%**; empirical upper-tail p = **0.0040**.
- Panel leave-one-out raw rho range: **0.628 to 0.654**.
- Panel leave-one-out residual-target rho range: **0.220 to 0.238**.

## Interpretation

The control asks whether the H&E-to-DM1/RAI result is just another spatially smooth random gene module. A high percentile supports specificity; a lower percentile would force a broader tissue-state framing. Leave-one-out stability tests whether one extreme thyroid-lineage gene drives the result.

## Actual And Leave-One-Out Results

| target             | mode                     |     n |   pooled_rho |   sample_centered_rho |   sample_centered_p |   median_sample_rho |   samples_rho_gt_0_2 | class      |
|:-------------------|:-------------------------|------:|-------------:|----------------------:|--------------------:|--------------------:|---------------------:|:-----------|
| actual_DM1_low_RAI | raw_smoothed             | 15489 |       0.6257 |                0.6444 |                   0 |              0.586  |                    4 | actual_dm1 |
| loo_drop_DIO1      | raw_smoothed             | 15489 |       0.6386 |                0.6523 |                   0 |              0.5753 |                    4 | panel_loo  |
| loo_drop_FOXE1     | raw_smoothed             | 15489 |       0.6303 |                0.6337 |                   0 |              0.5774 |                    4 | panel_loo  |
| loo_drop_NKX2-1    | raw_smoothed             | 15489 |       0.6337 |                0.6406 |                   0 |              0.5771 |                    4 | panel_loo  |
| loo_drop_PAX8      | raw_smoothed             | 15489 |       0.6306 |                0.6433 |                   0 |              0.5796 |                    4 | panel_loo  |
| loo_drop_SLC5A5    | raw_smoothed             | 15489 |       0.6456 |                0.6542 |                   0 |              0.5945 |                    4 | panel_loo  |
| loo_drop_TG        | raw_smoothed             | 15489 |       0.6147 |                0.633  |                   0 |              0.6048 |                    4 | panel_loo  |
| loo_drop_TPO       | raw_smoothed             | 15489 |       0.5802 |                0.6352 |                   0 |              0.5843 |                    4 | panel_loo  |
| loo_drop_TSHR      | raw_smoothed             | 15489 |       0.5961 |                0.6276 |                   0 |              0.5678 |                    4 | panel_loo  |
| actual_DM1_low_RAI | coord_qc_residual_target | 15489 |       0.2162 |                0.2323 |                   0 |              0.1947 |                    2 | actual_dm1 |
| loo_drop_DIO1      | coord_qc_residual_target | 15489 |       0.2215 |                0.2362 |                   0 |              0.1985 |                    2 | panel_loo  |
| loo_drop_FOXE1     | coord_qc_residual_target | 15489 |       0.2135 |                0.2307 |                   0 |              0.1903 |                    2 | panel_loo  |
| loo_drop_NKX2-1    | coord_qc_residual_target | 15489 |       0.2121 |                0.227  |                   0 |              0.1874 |                    2 | panel_loo  |
| loo_drop_PAX8      | coord_qc_residual_target | 15489 |       0.2164 |                0.2306 |                   0 |              0.1942 |                    2 | panel_loo  |
| loo_drop_SLC5A5    | coord_qc_residual_target | 15489 |       0.2221 |                0.238  |                   0 |              0.1988 |                    2 | panel_loo  |
| loo_drop_TG        | coord_qc_residual_target | 15489 |       0.205  |                0.2211 |                   0 |              0.1923 |                    2 | panel_loo  |
| loo_drop_TPO       | coord_qc_residual_target | 15489 |       0.206  |                0.2242 |                   0 |              0.1904 |                    2 | panel_loo  |
| loo_drop_TSHR      | coord_qc_residual_target | 15489 |       0.2043 |                0.2196 |                   0 |              0.1847 |                    2 | panel_loo  |

## Random Module Null Summary

| mode                     |   null_mean |   null_median |   null_p90 |   null_p95 |   actual |
|:-------------------------|------------:|--------------:|-----------:|-----------:|---------:|
| raw_smoothed             |      0.5388 |        0.5517 |     0.6111 |     0.6234 |   0.6444 |
| coord_qc_residual_target |      0.1338 |        0.1386 |     0.1842 |     0.1971 |   0.2323 |

## Top Random Raw Modules

| target                  | mode         |     n |   pooled_rho |   sample_centered_rho |   sample_centered_p |   median_sample_rho |   samples_rho_gt_0_2 | class          |
|:------------------------|:-------------|------:|-------------:|----------------------:|--------------------:|--------------------:|---------------------:|:---------------|
| null_random_matched_168 | raw_smoothed | 15489 |       0.6182 |                0.6539 |                   0 |              0.5734 |                    4 | random_matched |
| null_random_matched_124 | raw_smoothed | 15489 |       0.6841 |                0.6509 |                   0 |              0.5923 |                    4 | random_matched |
| null_random_matched_237 | raw_smoothed | 15489 |       0.6854 |                0.6397 |                   0 |              0.5482 |                    4 | random_matched |
| null_random_matched_031 | raw_smoothed | 15489 |       0.6623 |                0.6379 |                   0 |              0.5045 |                    4 | random_matched |
| null_random_matched_105 | raw_smoothed | 15489 |       0.6645 |                0.6375 |                   0 |              0.5002 |                    4 | random_matched |
| null_random_matched_244 | raw_smoothed | 15489 |       0.6523 |                0.6326 |                   0 |              0.5609 |                    4 | random_matched |
| null_random_matched_172 | raw_smoothed | 15489 |       0.6481 |                0.632  |                   0 |              0.5401 |                    4 | random_matched |
| null_random_matched_087 | raw_smoothed | 15489 |       0.657  |                0.6317 |                   0 |              0.5418 |                    4 | random_matched |
| null_random_matched_096 | raw_smoothed | 15489 |       0.6814 |                0.6314 |                   0 |              0.5505 |                    4 | random_matched |
| null_random_matched_042 | raw_smoothed | 15489 |       0.6373 |                0.6283 |                   0 |              0.5618 |                    4 | random_matched |
| null_random_matched_145 | raw_smoothed | 15489 |       0.6318 |                0.628  |                   0 |              0.4939 |                    4 | random_matched |
| null_random_matched_011 | raw_smoothed | 15489 |       0.6263 |                0.6257 |                   0 |              0.5099 |                    4 | random_matched |

## Top Random Residual Modules

| target                  | mode                     |     n |   pooled_rho |   sample_centered_rho |   sample_centered_p |   median_sample_rho |   samples_rho_gt_0_2 | class          |
|:------------------------|:-------------------------|------:|-------------:|----------------------:|--------------------:|--------------------:|---------------------:|:---------------|
| null_random_matched_169 | coord_qc_residual_target | 15489 |       0.2277 |                0.231  |                   0 |              0.2489 |                    3 | random_matched |
| null_random_matched_096 | coord_qc_residual_target | 15489 |       0.2158 |                0.2255 |                   0 |              0.1971 |                    2 | random_matched |
| null_random_matched_124 | coord_qc_residual_target | 15489 |       0.208  |                0.2235 |                   0 |              0.2156 |                    3 | random_matched |
| null_random_matched_244 | coord_qc_residual_target | 15489 |       0.2086 |                0.221  |                   0 |              0.1979 |                    2 | random_matched |
| null_random_matched_008 | coord_qc_residual_target | 15489 |       0.2088 |                0.2176 |                   0 |              0.201  |                    2 | random_matched |
| null_random_matched_165 | coord_qc_residual_target | 15489 |       0.1947 |                0.2149 |                   0 |              0.2112 |                    3 | random_matched |
| null_random_matched_145 | coord_qc_residual_target | 15489 |       0.2237 |                0.214  |                   0 |              0.2127 |                    3 | random_matched |
| null_random_matched_105 | coord_qc_residual_target | 15489 |       0.1939 |                0.2006 |                   0 |              0.2039 |                    2 | random_matched |
| null_random_matched_113 | coord_qc_residual_target | 15489 |       0.1847 |                0.2002 |                   0 |              0.1828 |                    1 | random_matched |
| null_random_matched_238 | coord_qc_residual_target | 15489 |       0.1859 |                0.1993 |                   0 |              0.2036 |                    2 | random_matched |
| null_random_matched_172 | coord_qc_residual_target | 15489 |       0.1876 |                0.199  |                   0 |              0.1959 |                    2 | random_matched |
| null_random_matched_130 | coord_qc_residual_target | 15489 |       0.19   |                0.1982 |                   0 |              0.1997 |                    2 | random_matched |

## Random Module Gene Map

See `gse230424_random_matched_modules.tsv` for the full matched gene map.
