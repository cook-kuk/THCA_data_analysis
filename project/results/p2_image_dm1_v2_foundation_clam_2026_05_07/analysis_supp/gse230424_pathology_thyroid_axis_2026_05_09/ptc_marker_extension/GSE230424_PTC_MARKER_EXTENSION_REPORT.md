# GSE230424 PTC-marker extension

## Verdict

- PTC marker genes used: 57 from article Supplementary Table S6.
- H&E sample-centered rho for the PTC-marker spatial score: **0.159**.
- QC-only sample-centered rho: **0.415**; coord+QC rho: **0.367**.
- Strict coordinate+QC residual-target H&E rho: **0.071**.
- PTC+HT minus HT sample-mean PTC marker score: **0.162** (n=2/group, descriptive only).

## Interpretation

The article's own PTC-specific spatial marker program is H&E-predictable, but QC/tissue-density baselines remain strong. Treat this as a descriptive extension of Paper 2 image-to-spatial-RNA support, not a robust disease-group validation.

## Model Comparison

| target                    | model            |     n |   pooled_rho |   pooled_p |   sample_centered_rho |   sample_centered_p |   median_sample_rho |   samples_rho_gt_0_2 |
|:--------------------------|:-----------------|------:|-------------:|-----------:|----------------------:|--------------------:|--------------------:|---------------------:|
| PTC_specific_marker_score | HE_tile_features | 15489 |       0.0371 |          0 |                0.1592 |                   0 |              0.1307 |                    1 |
| PTC_specific_marker_score | Coord_only       | 15489 |      -0.1497 |          0 |               -0.081  |                   0 |             -0.0251 |                    0 |
| PTC_specific_marker_score | QC_only          | 15489 |       0.5399 |          0 |                0.4149 |                   0 |              0.427  |                    4 |
| PTC_specific_marker_score | Coord_QC         | 15489 |       0.4142 |          0 |                0.3672 |                   0 |              0.5262 |                    3 |

## Residual Target

| target                    | residual_target_adjustment   | model            |     n |   pooled_rho |   sample_centered_rho |      p |
|:--------------------------|:-----------------------------|:-----------------|------:|-------------:|----------------------:|-------:|
| PTC_specific_marker_score | coord_qc                     | HE_tile_features | 15489 |       0.0143 |                0.0708 | 0.0745 |
