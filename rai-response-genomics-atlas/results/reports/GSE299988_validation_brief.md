# GSE299988 — RAI panel validation brief

- n samples merged: 14
- panel genes available per sample: 8
- label group counts: {'avid': 5, 'refractory': 5}

## Tests

| comparison                      |   n_test |   n_ref |   median_test |   median_ref |   cohen_d_test_minus_ref |    mwu_p |   AUC_test_low_panel |
|:--------------------------------|---------:|--------:|--------------:|-------------:|-------------------------:|---------:|---------------------:|
| tumor_vs_non_neoplastic         |       10 |       4 |     0.0147973 |     0.618194 |                   -1.608 | 0.001998 |                  1   |
| refractory_vs_avid_within_tumor |        5 |       5 |     0.125956  |    -0.298239 |                    1.128 | 0.150794 |                  0.2 |

## Wording rule

Any claim derived from this report must declare the dataset tier (Tier 2 RAI avidity for GSE151179 / GSE299988; Tier 1 for Boucai 2023 if available). 
Use 'associated with' or 'stratifies', never 'predicts RAI response'.