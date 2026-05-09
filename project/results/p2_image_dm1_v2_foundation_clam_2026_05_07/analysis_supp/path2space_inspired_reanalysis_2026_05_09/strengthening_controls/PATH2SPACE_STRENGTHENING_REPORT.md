# Path2Space strengthening controls

## Verdict

- UNI H&E DM1/RAI slide-centered rho: **0.440**.
- Coordinate-only DM1/RAI slide-centered rho: 0.133; QC-only: 0.476; coordinate+QC: 0.474.
- UNI minus coordinate-only delta: **+0.308**; UNI minus coordinate+QC delta: -0.033.
- After within-slide QC residualization, UNI DM1/RAI residual alignment rho remains **0.303**; after coordinate+QC residualization it remains 0.298.
- ResNet50 H&E DM1/RAI slide-centered rho: 0.298.
- Spatial-domain DM1/RAI pooled rho: **0.437** across 96 coordinate domains.
- Within-slide permutation empirical p for DM1/RAI: **0.0010** (1000 permutations).

## Interpretation

The image-derived DM1/RAI signal survives slide-centering and beats coordinate-only plus ResNet50 baselines, but coordinate+QC is a strong measured-ST baseline. Therefore the safe claim is not that UNI uniquely explains all local RNA structure; it is that pathology carries a reproducible DM1/RAI spatial signal, a component remains after QC residualization, and domain averaging strengthens it. The boundary remains unchanged: this strengthens Paper 2 image-to-spatial-RNA framing, not Paper 1 MAPK causal mechanism.

## Model Comparison

| target            | model          |    n |   pooled_rho |   pooled_p |   slide_centered_rho |   slide_centered_p |   median_slide_rho |   slides_rho_gt_0_2 |   features |
|:------------------|:---------------|-----:|-------------:|-----------:|---------------------:|-------------------:|-------------------:|--------------------:|-----------:|
| DM1_like_score    | UNI_HE         | 3200 |       0.392  |     0      |               0.4404 |             0      |             0.3727 |                  12 |       1024 |
| DM1_like_score    | ResNet50_HE    | 3200 |       0.2701 |     0      |               0.2979 |             0      |             0.2472 |                   9 |       2048 |
| DM1_like_score    | Coord_only     | 3200 |       0.109  |     0      |               0.1325 |             0      |             0.1103 |                   7 |          2 |
| DM1_like_score    | QC_only        | 3200 |       0.4234 |     0      |               0.4763 |             0      |             0.478  |                  14 |          3 |
| DM1_like_score    | Coord_QC       | 3200 |       0.4281 |     0      |               0.4736 |             0      |             0.4402 |                  15 |          5 |
| DM1_like_score    | Stage_only     | 3200 |      -0.0396 |     0.0251 |               0.0527 |             0.0029 |           nan      |                   0 |          4 |
| DM1_like_score    | Coord_QC_stage | 3200 |       0.4196 |     0      |               0.4697 |             0      |             0.4392 |                  15 |          9 |
| TDS_like_score    | UNI_HE         | 3200 |       0.3807 |     0      |               0.4099 |             0      |             0.3379 |                  12 |       1024 |
| TDS_like_score    | ResNet50_HE    | 3200 |       0.252  |     0      |               0.2603 |             0      |             0.1985 |                   8 |       2048 |
| TDS_like_score    | Coord_only     | 3200 |       0.1244 |     0      |               0.1509 |             0      |             0.1804 |                   7 |          2 |
| TDS_like_score    | QC_only        | 3200 |       0.3858 |     0      |               0.4289 |             0      |             0.4105 |                  13 |          3 |
| TDS_like_score    | Coord_QC       | 3200 |       0.3928 |     0      |               0.4289 |             0      |             0.4119 |                  13 |          5 |
| TDS_like_score    | Stage_only     | 3200 |      -0.0635 |     0.0003 |              -0.026  |             0.141  |           nan      |                   0 |          4 |
| TDS_like_score    | Coord_QC_stage | 3200 |       0.3883 |     0      |               0.4246 |             0      |             0.4118 |                  13 |          9 |
| MAPK_output_score | UNI_HE         | 3200 |       0.2983 |     0      |               0.3283 |             0      |             0.3065 |                  12 |       1024 |
| MAPK_output_score | ResNet50_HE    | 3200 |       0.279  |     0      |               0.2989 |             0      |             0.2513 |                  12 |       2048 |
| MAPK_output_score | Coord_only     | 3200 |       0.0401 |     0.0235 |               0.0262 |             0.1381 |             0.037  |                   3 |          2 |
| MAPK_output_score | QC_only        | 3200 |       0.5136 |     0      |               0.6145 |             0      |             0.5884 |                  16 |          3 |
| MAPK_output_score | Coord_QC       | 3200 |       0.5101 |     0      |               0.6063 |             0      |             0.5841 |                  16 |          5 |
| MAPK_output_score | Stage_only     | 3200 |       0.0488 |     0.0057 |               0.0029 |             0.8681 |           nan      |                   0 |          4 |
| MAPK_output_score | Coord_QC_stage | 3200 |       0.47   |     0      |               0.6061 |             0      |             0.5835 |                  16 |          9 |
| CAF_ECM_score     | UNI_HE         | 3200 |       0.3184 |     0      |               0.3913 |             0      |             0.2895 |                   9 |       1024 |
| CAF_ECM_score     | ResNet50_HE    | 3200 |       0.1971 |     0      |               0.218  |             0      |             0.0876 |                   7 |       2048 |
| CAF_ECM_score     | Coord_only     | 3200 |      -0.134  |     0      |              -0.1375 |             0      |            -0.1338 |                   5 |          2 |
| CAF_ECM_score     | QC_only        | 3200 |       0.2042 |     0      |               0.2491 |             0      |             0.2853 |                  10 |          3 |
| CAF_ECM_score     | Coord_QC       | 3200 |       0.1277 |     0      |               0.1547 |             0      |             0.2188 |                   9 |          5 |
| CAF_ECM_score     | Stage_only     | 3200 |      -0.0411 |     0.0201 |               0.0063 |             0.7233 |           nan      |                   0 |          4 |
| CAF_ECM_score     | Coord_QC_stage | 3200 |       0.1136 |     0      |               0.1707 |             0      |             0.2314 |                  10 |          9 |
| Epithelial_score  | UNI_HE         | 3200 |       0.2437 |     0      |               0.3035 |             0      |             0.2736 |                  12 |       1024 |
| Epithelial_score  | ResNet50_HE    | 3200 |       0.2424 |     0      |               0.2797 |             0      |             0.2417 |                  10 |       2048 |
| Epithelial_score  | Coord_only     | 3200 |       0.0435 |     0.0139 |               0.0319 |             0.0715 |             0.0014 |                   4 |          2 |
| Epithelial_score  | QC_only        | 3200 |       0.5077 |     0      |               0.6137 |             0      |             0.5513 |                  15 |          3 |
| Epithelial_score  | Coord_QC       | 3200 |       0.5049 |     0      |               0.6024 |             0      |             0.5465 |                  15 |          5 |
| Epithelial_score  | Stage_only     | 3200 |       0.0535 |     0.0025 |               0.0339 |             0.0555 |           nan      |                   0 |          4 |
| Epithelial_score  | Coord_QC_stage | 3200 |       0.4732 |     0      |               0.6025 |             0      |             0.5457 |                  15 |          9 |

## Domain Summary

| target            |   n_domains |   pooled_domain_rho |   pooled_domain_p |   slide_centered_domain_rho |   slide_centered_domain_p |   median_slide_domain_rho |   slides_domain_rho_gt_0_2 |
|:------------------|------------:|--------------------:|------------------:|----------------------------:|--------------------------:|--------------------------:|---------------------------:|
| DM1_like_score    |          96 |              0.4368 |            0      |                      0.5811 |                         0 |                    0.5429 |                         10 |
| TDS_like_score    |          96 |              0.4063 |            0      |                      0.5617 |                         0 |                    0.3429 |                         10 |
| MAPK_output_score |          96 |              0.3958 |            0.0001 |                      0.5267 |                         0 |                    0.6    |                         13 |
| CAF_ECM_score     |          96 |              0.3749 |            0.0002 |                      0.6581 |                         0 |                    0.6286 |                         11 |
| Epithelial_score  |          96 |              0.2888 |            0.0043 |                      0.5189 |                         0 |                    0.5143 |                         12 |

## Permutation Null

| target            |   observed_slide_centered_rho |   null_mean |   null_sd |   null_p95 |   null_p99 |   empirical_p_greater_equal |   n_permutations |
|:------------------|------------------------------:|------------:|----------:|-----------:|-----------:|----------------------------:|-----------------:|
| DM1_like_score    |                        0.4404 |      0.0038 |    0.0185 |     0.035  |     0.0451 |                       0.001 |             1000 |
| TDS_like_score    |                        0.4099 |      0.0039 |    0.017  |     0.0303 |     0.0423 |                       0.001 |             1000 |
| MAPK_output_score |                        0.3283 |      0.0015 |    0.018  |     0.0311 |     0.0434 |                       0.001 |             1000 |

## Residual Alignment

| target            | within_slide_adjustment   |   residual_alignment_rho |   p |    n |
|:------------------|:--------------------------|-------------------------:|----:|-----:|
| DM1_like_score    | slide_mean_only           |                   0.4404 |   0 | 3200 |
| DM1_like_score    | coord_only                |                   0.4382 |   0 | 3200 |
| DM1_like_score    | qc_only                   |                   0.3027 |   0 | 3200 |
| DM1_like_score    | coord_qc                  |                   0.2975 |   0 | 3200 |
| TDS_like_score    | slide_mean_only           |                   0.4099 |   0 | 3200 |
| TDS_like_score    | coord_only                |                   0.4054 |   0 | 3200 |
| TDS_like_score    | qc_only                   |                   0.2774 |   0 | 3200 |
| TDS_like_score    | coord_qc                  |                   0.2733 |   0 | 3200 |
| MAPK_output_score | slide_mean_only           |                   0.3283 |   0 | 3200 |
| MAPK_output_score | coord_only                |                   0.3381 |   0 | 3200 |
| MAPK_output_score | qc_only                   |                   0.1307 |   0 | 3200 |
| MAPK_output_score | coord_qc                  |                   0.1513 |   0 | 3200 |
| CAF_ECM_score     | slide_mean_only           |                   0.3913 |   0 | 3200 |
| CAF_ECM_score     | coord_only                |                   0.3188 |   0 | 3200 |
| CAF_ECM_score     | qc_only                   |                   0.3431 |   0 | 3200 |
| CAF_ECM_score     | coord_qc                  |                   0.2801 |   0 | 3200 |
| Epithelial_score  | slide_mean_only           |                   0.3035 |   0 | 3200 |
| Epithelial_score  | coord_only                |                   0.333  |   0 | 3200 |
| Epithelial_score  | qc_only                   |                   0.1196 |   0 | 3200 |
| Epithelial_score  | coord_qc                  |                   0.153  |   0 | 3200 |
