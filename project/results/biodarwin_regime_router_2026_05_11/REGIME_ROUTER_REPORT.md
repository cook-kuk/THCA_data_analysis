# BioDarwin regime router v1

Generated: 2026-05-11T01:09:11

## Best validation-like
| split                    | algorithm            |   n |   positives |   positive_rate |    AUPRC |    AUROC |   top5_precision |   top5_hits |   top10_precision |   top10_hits |   top24_precision |   top24_hits |
|:-------------------------|:---------------------|----:|------------:|----------------:|---------:|---------:|-----------------:|------------:|------------------:|-------------:|------------------:|-------------:|
| academic_validation_like | Current_mode_bank_v4 | 430 |          11 |       0.0255814 | 0.886307 | 0.972445 |                1 |           5 |               0.9 |            9 |          0.416667 |           10 |

## Best industrial locked
| split                | algorithm        |   n |   positives |   positive_rate |    AUPRC |    AUROC |   top5_precision |   top5_hits |   top10_precision |   top10_hits |   top24_precision |   top24_hits |
|:---------------------|:-----------------|----:|------------:|----------------:|---------:|---------:|-----------------:|------------:|------------------:|-------------:|------------------:|-------------:|
| industrial_locked_v0 | Regime_router_v1 |  25 |           9 |            0.36 | 0.629766 | 0.645833 |              0.8 |           4 |               0.4 |            4 |             0.375 |            9 |

## Claim boundary
- Internal-feature rows use current mode bank.
- External rows use public anchor + BigMHC with RF disagreement penalty.
- This is a router, not a locked SOTA claim.