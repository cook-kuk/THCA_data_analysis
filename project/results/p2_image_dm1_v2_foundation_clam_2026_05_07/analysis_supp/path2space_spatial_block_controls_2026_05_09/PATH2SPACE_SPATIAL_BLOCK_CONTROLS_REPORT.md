# Path2Space-inspired spatial-block controls

## Verdict

- GSE250521 block-centered rho: **0.380**; top-decile block-stratified permutation p = **0.0010**.
- GSE230424 block-centered rho: **0.497**; top-decile block-stratified permutation p = **0.0010**.
- GSE230424 coord+QC residual block-centered rho: **0.198**; top-decile block-stratified permutation p = **0.0010**.

## Interpretation

The block-centered test removes broad coordinate-domain means before computing correlation. The block-stratified hotspot null shuffles predicted hotspot labels only within spatial blocks, preserving broad predicted hotspot density by region. Residual signal under these controls argues against a purely broad-domain artifact.

## Summary

| evidence                                    |     n |   n_groups |   n_spatial_blocks |   sample_centered_rho |   sample_centered_p |   spatial_block_centered_rho |   spatial_block_centered_p |   top10_precision |   top10_precision_lift |   top10_tp_over_expected |   top10_or |   top10_fisher_p |   top10_domain_stratified_perm_p |   domain_stratified_null_tp_mean |   domain_stratified_null_tp_p95 |   observed_tp |
|:--------------------------------------------|------:|-----------:|-------------------:|----------------------:|--------------------:|-----------------------------:|---------------------------:|------------------:|-----------------------:|-------------------------:|-----------:|-----------------:|---------------------------------:|---------------------------------:|--------------------------------:|--------------:|
| GSE250521 UNI DM1/RAI smoothed              |  3200 |         16 |                 64 |                0.4404 |                   0 |                       0.3799 |                          0 |            0.35   |                 3.5    |                   3.5    |     6.9172 |                0 |                            0.001 |                           52.717 |                              62 |           112 |
| GSE230424 H&E DM1/low-RAI smoothed          | 15489 |          4 |                 43 |                0.6444 |                   0 |                       0.4971 |                          0 |            0.3327 |                 3.3224 |                   3.3224 |     6.2153 |                0 |                            0.001 |                          265.716 |                             287 |           516 |
| GSE230424 H&E DM1/low-RAI coord+QC residual | 15489 |          4 |                 43 |                0.2323 |                   0 |                       0.198  |                          0 |            0.1631 |                 1.629  |                   1.629  |     1.8981 |                0 |                            0.001 |                          176.249 |                             193 |           253 |

## Per-Slide/Sample Detail

| evidence                                    | group             |    n |   n_blocks |   sample_rho |   block_centered_rho |   top10_precision |   top10_precision_lift |   top10_or |
|:--------------------------------------------|:------------------|-----:|-----------:|-------------:|---------------------:|------------------:|-----------------------:|-----------:|
| GSE250521 UNI DM1/RAI smoothed              | GSM7980860_N-1    |  200 |          4 |       0.5598 |               0.3533 |            0.45   |                 4.5    |    12.5702 |
| GSE250521 UNI DM1/RAI smoothed              | GSM7980861_N-2    |  200 |          4 |       0.4073 |               0.433  |            0.4    |                 4      |     9.3333 |
| GSE250521 UNI DM1/RAI smoothed              | GSM7980862_N-3    |  200 |          4 |       0.2627 |               0.2384 |            0.2    |                 2      |     2.5625 |
| GSE250521 UNI DM1/RAI smoothed              | GSM7980863_N-4    |  200 |          4 |       0.3214 |               0.3177 |            0.35   |                 3.5    |     6.9172 |
| GSE250521 UNI DM1/RAI smoothed              | GSM7980864_PTC-1  |  200 |          4 |       0.466  |               0.4847 |            0.15   |                 1.5    |     1.692  |
| GSE250521 UNI DM1/RAI smoothed              | GSM7980865_PTC-2  |  200 |          4 |       0.6057 |               0.6229 |            0.55   |                 5.5    |    23.2222 |
| GSE250521 UNI DM1/RAI smoothed              | GSM7980866_PTC-3  |  200 |          4 |       0.1922 |               0.2647 |            0.6    |                 6      |    32.25   |
| GSE250521 UNI DM1/RAI smoothed              | GSM7980867_PTC-4  |  200 |          4 |       0.3382 |               0.417  |            0.5    |                 5      |    17      |
| GSE250521 UNI DM1/RAI smoothed              | GSM7980868_LPTC-1 |  200 |          4 |       0.4815 |               0.2625 |            0.6    |                 6      |    32.25   |
| GSE250521 UNI DM1/RAI smoothed              | GSM7980869_LPTC-2 |  200 |          4 |       0.6565 |               0.5146 |            0      |                 0      |     0      |
| GSE250521 UNI DM1/RAI smoothed              | GSM7980870_LPTC-3 |  200 |          4 |       0.5242 |               0.4935 |            0.1    |                 1      |     1      |
| GSE250521 UNI DM1/RAI smoothed              | GSM7980871_LPTC-4 |  200 |          4 |       0.3146 |               0.4318 |            0.35   |                 3.5    |     6.9172 |
| GSE250521 UNI DM1/RAI smoothed              | GSM7980872_ATC-1  |  200 |          4 |       0.0298 |               0.0491 |            0.3    |                 3      |     5.0816 |
| GSE250521 UNI DM1/RAI smoothed              | GSM7980873_ATC-2  |  200 |          4 |      -0.0488 |               0.0055 |            0.05   |                 0.5    |     0.446  |
| GSE250521 UNI DM1/RAI smoothed              | GSM7980874_ATC-3  |  200 |          4 |       0.0861 |               0.1579 |            0.15   |                 1.5    |     1.692  |
| GSE250521 UNI DM1/RAI smoothed              | GSM7980875_ATC-4  |  200 |          4 |       0.5313 |               0.4833 |            0.85   |                 8.5    |   334.333  |
| GSE230424 H&E DM1/low-RAI smoothed          | P1                | 3647 |         10 |       0.5638 |               0.4386 |            0.2959 |                 2.9565 |     4.9463 |
| GSE230424 H&E DM1/low-RAI smoothed          | P2                | 3154 |          9 |       0.6082 |               0.4187 |            0.3133 |                 3.127  |     5.5104 |
| GSE230424 H&E DM1/low-RAI smoothed          | P3                | 4062 |         12 |       0.7871 |               0.6894 |            0.3194 |                 3.1878 |     5.7233 |
| GSE230424 H&E DM1/low-RAI smoothed          | P4                | 4626 |         12 |       0.445  |               0.3837 |            0.3866 |                 3.8628 |     8.6087 |
| GSE230424 H&E DM1/low-RAI coord+QC residual | P1                | 3647 |         10 |       0.1667 |               0.1661 |            0.2    |                 1.9984 |     2.5599 |
| GSE230424 H&E DM1/low-RAI coord+QC residual | P2                | 3154 |          9 |       0.3602 |               0.3532 |            0.1139 |                 1.1371 |     1.1746 |
| GSE230424 H&E DM1/low-RAI coord+QC residual | P3                | 4062 |         12 |       0.1603 |               0.0839 |            0.0442 |                 0.4414 |     0.3885 |
| GSE230424 H&E DM1/low-RAI coord+QC residual | P4                | 4626 |         12 |       0.2226 |               0.2019 |            0.2721 |                 2.719  |     4.2448 |
