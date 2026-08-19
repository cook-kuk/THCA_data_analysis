# GSE250521 stage-generalization controls

## Verdict

- Existing leave-one-slide-out DM1/RAI stage-centered rho: **0.392**.
- Existing leave-one-slide-out DM1/RAI by-stage pooled rho: PT **0.246**, PTC **0.487**, LPTC **0.471**, ATC **0.191**.
- Leave-one-stage-out raw DM1/RAI pooled rho: **0.335**; stage-centered rho **0.401**.
- Leave-one-stage-out coord+QC residual-target DM1/RAI pooled rho: **0.140**; stage-centered rho **0.170**.

## Interpretation

The existing leave-one-slide-out signal is not only a global stage mean because stage-centered evaluation remains positive. The raw signal also transfers in leave-one-stage-out evaluation, especially into PTC and LPTC. The caveat is the stricter coord+QC residual target: it stays positive but is weaker, and ATC is the least stable stage.

Safe wording: use GSE250521 as a local image-to-DM1/RAI spatial result that is not just stage mean structure, with stage heterogeneity and weaker residual-stage transfer disclosed.

## Existing Leave-One-Slide-Out Summary

| target            | evaluation    |   pooled_rho |   pooled_p |   slide_centered_rho |   slide_centered_p |   stage_centered_rho |   stage_centered_p |    n |
|:------------------|:--------------|-------------:|-----------:|---------------------:|-------------------:|---------------------:|-------------------:|-----:|
| DM1_like_score    | existing_loso |       0.392  |          0 |               0.4404 |                  0 |               0.3919 |                  0 | 3200 |
| TDS_like_score    | existing_loso |       0.3807 |          0 |               0.4099 |                  0 |               0.3753 |                  0 | 3200 |
| MAPK_output_score | existing_loso |       0.2983 |          0 |               0.3283 |                  0 |               0.2916 |                  0 | 3200 |

## Existing Leave-One-Slide-Out By Stage

| target            | evaluation    | stage   |   n_spots |   pooled_within_stage_rho |   pooled_within_stage_p |   slide_centered_within_stage_rho |   slide_centered_within_stage_p |   median_slide_rho |   slides_rho_gt_0_2 |
|:------------------|:--------------|:--------|----------:|--------------------------:|------------------------:|----------------------------------:|--------------------------------:|-------------------:|--------------------:|
| DM1_like_score    | existing_loso | PT      |       800 |                    0.2462 |                       0 |                            0.4018 |                               0 |             0.3643 |                   4 |
| DM1_like_score    | existing_loso | PTC     |       800 |                    0.4868 |                       0 |                            0.4651 |                               0 |             0.4021 |                   3 |
| DM1_like_score    | existing_loso | LPTC    |       800 |                    0.4711 |                       0 |                            0.5319 |                               0 |             0.5029 |                   4 |
| DM1_like_score    | existing_loso | ATC     |       800 |                    0.1908 |                       0 |                            0.2383 |                               0 |             0.0579 |                   1 |
| TDS_like_score    | existing_loso | PT      |       800 |                    0.2226 |                       0 |                            0.3786 |                               0 |             0.3379 |                   4 |
| TDS_like_score    | existing_loso | PTC     |       800 |                    0.4699 |                       0 |                            0.4256 |                               0 |             0.3256 |                   3 |
| TDS_like_score    | existing_loso | LPTC    |       800 |                    0.4654 |                       0 |                            0.5097 |                               0 |             0.4908 |                   4 |
| TDS_like_score    | existing_loso | ATC     |       800 |                    0.2468 |                       0 |                            0.2265 |                               0 |             0.0986 |                   1 |
| MAPK_output_score | existing_loso | PT      |       800 |                    0.1653 |                       0 |                            0.1717 |                               0 |             0.1752 |                   2 |
| MAPK_output_score | existing_loso | PTC     |       800 |                    0.3199 |                       0 |                            0.2956 |                               0 |             0.3114 |                   4 |
| MAPK_output_score | existing_loso | LPTC    |       800 |                    0.3772 |                       0 |                            0.4908 |                               0 |             0.4993 |                   4 |
| MAPK_output_score | existing_loso | ATC     |       800 |                    0.2217 |                       0 |                            0.2216 |                               0 |             0.1988 |                   2 |

## Leave-One-Stage-Out Summary

| target            | mode                     | evaluation          |   pooled_rho |   pooled_p |   stage_centered_rho |   stage_centered_p |    n |
|:------------------|:-------------------------|:--------------------|-------------:|-----------:|---------------------:|-------------------:|-----:|
| DM1_like_score    | raw_smoothed             | leave_one_stage_out |       0.3353 |     0      |               0.4012 |             0      | 3200 |
| DM1_like_score    | coord_qc_residual_target | leave_one_stage_out |       0.1405 |     0      |               0.17   |             0      | 3200 |
| TDS_like_score    | raw_smoothed             | leave_one_stage_out |       0.3218 |     0      |               0.3769 |             0      | 3200 |
| TDS_like_score    | coord_qc_residual_target | leave_one_stage_out |       0.1365 |     0      |               0.1559 |             0      | 3200 |
| MAPK_output_score | raw_smoothed             | leave_one_stage_out |       0.2952 |     0      |               0.2785 |             0      | 3200 |
| MAPK_output_score | coord_qc_residual_target | leave_one_stage_out |       0.0612 |     0.0005 |               0.0714 |             0.0001 | 3200 |

## Leave-One-Stage-Out By Held Stage

| target            | mode                     | evaluation          | held_out_stage   |   n_spots |    rho |      p |   slide_centered_rho |   slide_centered_p |
|:------------------|:-------------------------|:--------------------|:-----------------|----------:|-------:|-------:|---------------------:|-------------------:|
| DM1_like_score    | raw_smoothed             | leave_one_stage_out | PT               |       800 | 0.2505 | 0      |               0.404  |             0      |
| DM1_like_score    | raw_smoothed             | leave_one_stage_out | PTC              |       800 | 0.4572 | 0      |               0.4082 |             0      |
| DM1_like_score    | raw_smoothed             | leave_one_stage_out | LPTC             |       800 | 0.4785 | 0      |               0.5572 |             0      |
| DM1_like_score    | raw_smoothed             | leave_one_stage_out | ATC              |       800 | 0.214  | 0      |               0.1972 |             0      |
| DM1_like_score    | coord_qc_residual_target | leave_one_stage_out | PT               |       800 | 0.1667 | 0      |               0.1909 |             0      |
| DM1_like_score    | coord_qc_residual_target | leave_one_stage_out | PTC              |       800 | 0.1965 | 0      |               0.1859 |             0      |
| DM1_like_score    | coord_qc_residual_target | leave_one_stage_out | LPTC             |       800 | 0.1512 | 0      |               0.2199 |             0      |
| DM1_like_score    | coord_qc_residual_target | leave_one_stage_out | ATC              |       800 | 0.0972 | 0.0059 |               0.1168 |             0.0009 |
| TDS_like_score    | raw_smoothed             | leave_one_stage_out | PT               |       800 | 0.2417 | 0      |               0.4136 |             0      |
| TDS_like_score    | raw_smoothed             | leave_one_stage_out | PTC              |       800 | 0.4239 | 0      |               0.3588 |             0      |
| TDS_like_score    | raw_smoothed             | leave_one_stage_out | LPTC             |       800 | 0.4466 | 0      |               0.5237 |             0      |
| TDS_like_score    | raw_smoothed             | leave_one_stage_out | ATC              |       800 | 0.2276 | 0      |               0.1919 |             0      |
| TDS_like_score    | coord_qc_residual_target | leave_one_stage_out | PT               |       800 | 0.1569 | 0      |               0.1995 |             0      |
| TDS_like_score    | coord_qc_residual_target | leave_one_stage_out | PTC              |       800 | 0.1833 | 0      |               0.1736 |             0      |
| TDS_like_score    | coord_qc_residual_target | leave_one_stage_out | LPTC             |       800 | 0.1357 | 0.0001 |               0.2001 |             0      |
| TDS_like_score    | coord_qc_residual_target | leave_one_stage_out | ATC              |       800 | 0.082  | 0.0204 |               0.0852 |             0.0159 |
| MAPK_output_score | raw_smoothed             | leave_one_stage_out | PT               |       800 | 0.1152 | 0.0011 |               0.0847 |             0.0166 |
| MAPK_output_score | raw_smoothed             | leave_one_stage_out | PTC              |       800 | 0.2888 | 0      |               0.285  |             0      |
| MAPK_output_score | raw_smoothed             | leave_one_stage_out | LPTC             |       800 | 0.4212 | 0      |               0.4533 |             0      |
| MAPK_output_score | raw_smoothed             | leave_one_stage_out | ATC              |       800 | 0.1804 | 0      |               0.146  |             0      |
| MAPK_output_score | coord_qc_residual_target | leave_one_stage_out | PT               |       800 | 0.0443 | 0.2107 |               0.0418 |             0.238  |
| MAPK_output_score | coord_qc_residual_target | leave_one_stage_out | PTC              |       800 | 0.0882 | 0.0126 |               0.0939 |             0.0079 |
| MAPK_output_score | coord_qc_residual_target | leave_one_stage_out | LPTC             |       800 | 0.0746 | 0.035  |               0.0744 |             0.0355 |
| MAPK_output_score | coord_qc_residual_target | leave_one_stage_out | ATC              |       800 | 0.0617 | 0.0811 |               0.0739 |             0.0365 |
