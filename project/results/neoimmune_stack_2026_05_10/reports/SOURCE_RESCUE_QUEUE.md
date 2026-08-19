# Source Rescue Queue

This queue is the practical handoff for source families where public reconstruction is still thin.

## Reading rule
- Positive labels with high rescue priority are the candidates most worth discussing with wet-lab collaborators.
- Negative labels can still be useful as failure cases, but they are not rescue targets.

| dataset_source                           |   rows |   positives |   top_priority |   median_priority |   mean_clean |   mean_public_best |   mean_delta |
|:-----------------------------------------|-------:|------------:|---------------:|------------------:|-------------:|-------------------:|-------------:|
| CEDAR                                    |    863 |         816 |       0.8539   |          0.660238 |     0.874479 |           0.837567 |    0.0434579 |
| dbPepNeo2_MHCI                           |    344 |         344 |       0.83892  |          0.656049 |     0.892926 |           0.745615 |    0.147311  |
| TESLA_mmc4                               |    273 |          63 |       0.752232 |          0.6237   |     0.938864 |           0.923305 |    0.0156216 |
| TESLA_mmc7_validation                    |     41 |           8 |       0.685411 |          0.507685 |     0.888607 |           0.90787  |   -0.0189757 |
| dbPepNeo2_MHCII                          |      9 |           9 |       0.609821 |          0.525121 |     0.823037 |         nan        |  nan         |
| IMPROVE_Neoepitopes_CEDAR_benchmark_data |     40 |          40 |       0.487411 |          0.445179 |     0.736836 |         nan        |  nan         |

## Priority rule
- High `clean_science_score` plus positive `local_rescue_delta_vs_public` means the candidate is rescued by local branch evidence beyond public predictors.
- High `model_disagreement_score` means the candidate is a good audit case even when wet-lab validation is not immediate.

## Claim boundary
- This is a retrospective prioritization queue only.
- No candidate is a validated vaccine target without orthogonal assay support.
