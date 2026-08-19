# GSE151179 — RAI panel validation brief

- n samples merged: 52
- panel genes available per sample: 8
- label group counts: {'avid': 4, 'refractory': 35}

## Tests

| comparison                            |   n_test |   n_ref |   median_test |   median_ref |   cohen_d_test_minus_ref |         mwu_p |   AUC_test_low_panel |   n_groups | labels                              |     stat |          p |
|:--------------------------------------|---------:|--------:|--------------:|-------------:|-------------------------:|--------------:|---------------------:|-----------:|:------------------------------------|---------:|-----------:|
| tumor_vs_non_neoplastic               |       39 |      13 |     -0.162736 |    0.850811  |                   -1.772 |   9.45208e-07 |             0.95858  |        nan | nan                                 | nan      | nan        |
| refractory_vs_avid_within_tumor       |       35 |       4 |     -0.168386 |    0.101794  |                   -0.244 |   0.487897    |             0.614286 |        nan | nan                                 | nan      | nan        |
| uptake_no_vs_yes_within_tumor         |       20 |      19 |     -0.201447 |   -0.0771816 |                   -0.266 |   0.683703    |             0.539474 |        nan | nan                                 | nan      | nan        |
| persistence_vs_remission_within_tumor |       35 |       4 |     -0.168386 |    0.101794  |                   -0.244 |   0.487897    |             0.614286 |        nan | nan                                 | nan      | nan        |
| kruskal_by_lesion_class_tumor_only    |      nan |     nan |    nan        |  nan         |                  nan     | nan           |           nan        |          4 | brafv600e:15,fusion:9,ptert:3,wt:11 |   2.7307 |   0.435035 |

## Wording rule

Any claim derived from this report must declare the dataset tier (Tier 2 RAI avidity for GSE151179 / GSE299988; Tier 1 for Boucai 2023 if available). 
Use 'associated with' or 'stratifies', never 'predicts RAI response'.