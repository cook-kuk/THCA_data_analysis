# BioDarwin RF anchor expansion v1

Generated: 2026-05-11T01:03:02

## Best validation-like
| split                    | algorithm            |   n |   positives |   positive_rate |    AUPRC |    AUROC |   top5_precision |   top5_hits |   top10_precision |   top10_hits |   top24_precision |   top24_hits |
|:-------------------------|:---------------------|----:|------------:|----------------:|---------:|---------:|-----------------:|------------:|------------------:|-------------:|------------------:|-------------:|
| academic_validation_like | Current_mode_bank_v4 | 430 |          11 |       0.0255814 | 0.886307 | 0.972445 |                1 |           5 |               0.9 |            9 |          0.416667 |           10 |

## Best industrial locked
| split                | algorithm        |   n |   positives |   positive_rate |    AUPRC |    AUROC |   top5_precision |   top5_hits |   top10_precision |   top10_hits |   top24_precision |   top24_hits |
|:---------------------|:-----------------|----:|------------:|----------------:|---------:|---------:|-----------------:|------------:|------------------:|-------------:|------------------:|-------------:|
| industrial_locked_v0 | Public_anchor_v3 |  25 |           9 |            0.36 | 0.627979 | 0.638889 |              0.8 |           4 |               0.5 |            5 |             0.375 |            9 |

## Focus metrics
| split                    | algorithm              |   n |   positives |   positive_rate |    AUPRC |    AUROC |   top5_precision |   top5_hits |   top10_precision |   top10_hits |   top24_precision |   top24_hits |
|:-------------------------|:-----------------------|----:|------------:|----------------:|---------:|---------:|-----------------:|------------:|------------------:|-------------:|------------------:|-------------:|
| academic_validation_like | RF_biophys             | 430 |          11 |       0.0255814 | 0.323505 | 0.890432 |              0.4 |           2 |               0.3 |            3 |          0.166667 |            4 |
| academic_validation_like | RF_anchor_heuristic_v1 | 430 |          11 |       0.0255814 | 0.401705 | 0.944456 |              0.4 |           2 |               0.3 |            3 |          0.291667 |            7 |
| academic_validation_like | RF_anchor_stack_v1     | 430 |          11 |       0.0255814 | 0.343092 | 0.939032 |              0.4 |           2 |               0.4 |            4 |          0.25     |            6 |
| academic_validation_like | Current_mode_bank_v4   | 430 |          11 |       0.0255814 | 0.886307 | 0.972445 |              1   |           5 |               0.9 |            9 |          0.416667 |           10 |
| academic_validation_like | Public_anchor_v3       | 430 |          11 |       0.0255814 | 0.31884  | 0.797787 |              0.6 |           3 |               0.4 |            4 |          0.25     |            6 |
| academic_validation_like | BigMHC_IM              | 430 |          11 |       0.0255814 | 0.333705 | 0.827728 |              0.6 |           3 |               0.4 |            4 |          0.25     |            6 |
| industrial_locked_v0     | RF_biophys             |  25 |           9 |       0.36      | 0.316413 | 0.354167 |              0.2 |           1 |               0.3 |            3 |          0.375    |            9 |
| industrial_locked_v0     | RF_anchor_heuristic_v1 |  25 |           9 |       0.36      | 0.402111 | 0.506944 |              0.4 |           2 |               0.4 |            4 |          0.375    |            9 |
| industrial_locked_v0     | RF_anchor_stack_v1     |  25 |           9 |       0.36      | 0.345172 | 0.423611 |              0.2 |           1 |               0.3 |            3 |          0.333333 |            8 |
| industrial_locked_v0     | Current_mode_bank_v4   |  25 |           9 |       0.36      | 0.627979 | 0.638889 |              0.8 |           4 |               0.5 |            5 |          0.375    |            9 |
| industrial_locked_v0     | Public_anchor_v3       |  25 |           9 |       0.36      | 0.627979 | 0.638889 |              0.8 |           4 |               0.5 |            5 |          0.375    |            9 |
| industrial_locked_v0     | BigMHC_IM              |  25 |           9 |       0.36      | 0.578121 | 0.611111 |              0.8 |           4 |               0.4 |            4 |          0.375    |            9 |

## Claim boundary
- RF is used as an AUPRC anchor and routing feature.
- This is an experiment driver, not a locked SOTA claim.
- Overlap-aware interpretation still applies.