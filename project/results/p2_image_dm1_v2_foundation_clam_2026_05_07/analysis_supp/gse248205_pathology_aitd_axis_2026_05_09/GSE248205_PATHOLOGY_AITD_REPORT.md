# GSE248205 H&E-to-AITD spatial-axis scout

## Verdict

- Spots: 16,985 across 8 thyroid Visium samples (Control n=2, GD n=3, HT n=3).
- H&E tile features -> AP/TLS composite: pooled rho **0.590**, sample-centered rho **-0.090**.
- AP/TLS coordinate-only centered rho 0.096; group-only centered rho -0.014.
- Domain-level AP/TLS sample-centered rho **-0.164**.
- HLA-II/AP centered rho -0.046; B/TLS -0.103; Thyrocyte -0.011.

## Interpretation

This is a fast handcrafted-feature screen, not a UNI/foundation-model result. The pooled association is driven by sample/group structure and does not survive the sample-centered test. Treat GSE248205 as a caveat-only negative control for local H&E-to-AITD spatial prediction, not as positive Paper 2 evidence and not as Paper 1 mechanism support.

## Model comparison

| target           | model            |     n |   pooled_rho |   pooled_p |   sample_centered_rho |   sample_centered_p |   median_sample_rho |   samples_rho_gt_0_2 |   n_features |
|:-----------------|:-----------------|------:|-------------:|-----------:|----------------------:|--------------------:|--------------------:|---------------------:|-------------:|
| AP_TLS_composite | HE_tile_features | 16985 |       0.5902 |     0      |               -0.0903 |              0      |             -0.0473 |                    0 |           52 |
| AP_TLS_composite | Coord_only       | 16985 |       0.0022 |     0.7703 |                0.0958 |              0      |              0.0611 |                    3 |            2 |
| AP_TLS_composite | Group_only       | 16985 |       0.7914 |     0      |               -0.0137 |              0.0735 |            nan      |                    0 |            3 |
| AP_TLS_composite | Coord_group      | 16985 |       0.7612 |     0      |                0.0731 |              0      |              0.0361 |                    3 |            5 |
| HLA_II_AP        | HE_tile_features | 16985 |       0.5829 |     0      |               -0.0457 |              0      |             -0.0637 |                    1 |           52 |
| HLA_II_AP        | Coord_only       | 16985 |       0.0973 |     0      |                0.1828 |              0      |              0.1718 |                    4 |            2 |
| HLA_II_AP        | Group_only       | 16985 |       0.5593 |     0      |                0.0214 |              0.0054 |            nan      |                    0 |            3 |
| HLA_II_AP        | Coord_group      | 16985 |       0.5707 |     0      |                0.158  |              0      |              0.1806 |                    3 |            5 |
| B_TLS            | HE_tile_features | 16985 |       0.5023 |     0      |               -0.1033 |              0      |             -0.0905 |                    0 |           52 |
| B_TLS            | Coord_only       | 16985 |      -0.0272 |     0.0004 |                0.0248 |              0.0012 |             -0.0711 |                    3 |            2 |
| B_TLS            | Group_only       | 16985 |       0.8265 |     0      |               -0.0146 |              0.0579 |            nan      |                    0 |            3 |
| B_TLS            | Coord_group      | 16985 |       0.7875 |     0      |               -0.0133 |              0.0831 |             -0.0845 |                    2 |            5 |
| CD74_MIF_axis    | HE_tile_features | 16985 |       0.5394 |     0      |               -0.0553 |              0      |             -0.0341 |                    0 |           52 |
| CD74_MIF_axis    | Coord_only       | 16985 |       0.1168 |     0      |                0.1136 |              0      |              0.1053 |                    4 |            2 |
| CD74_MIF_axis    | Group_only       | 16985 |       0.514  |     0      |               -0.0537 |              0      |            nan      |                    0 |            3 |
| CD74_MIF_axis    | Coord_group      | 16985 |       0.5432 |     0      |                0.1078 |              0      |              0.0996 |                    4 |            5 |
| Thyrocyte        | HE_tile_features | 16985 |       0.3976 |     0      |               -0.0109 |              0.1559 |              0.0279 |                    1 |           52 |
| Thyrocyte        | Coord_only       | 16985 |      -0.176  |     0      |               -0.0566 |              0      |             -0.0004 |                    2 |            2 |
| Thyrocyte        | Group_only       | 16985 |       0.5229 |     0      |               -0.0364 |              0      |            nan      |                    0 |            3 |
| Thyrocyte        | Coord_group      | 16985 |       0.5891 |     0      |               -0.0401 |              0      |              0.086  |                    2 |            5 |

## Domain aggregation

| target           |   n_domains |   pooled_domain_rho |   sample_centered_domain_rho |      p |
|:-----------------|------------:|--------------------:|-----------------------------:|-------:|
| AP_TLS_composite |          77 |              0.6374 |                      -0.1644 | 0      |
| HLA_II_AP        |          77 |              0.6264 |                      -0.1417 | 0      |
| B_TLS            |          77 |              0.5162 |                      -0.2022 | 0      |
| CD74_MIF_axis    |          77 |              0.5646 |                      -0.0679 | 0      |
| Thyrocyte        |          77 |              0.4411 |                      -0.0999 | 0.0001 |
