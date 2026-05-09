# UNI-final leave-one-TSS-out audit

| Metric | Value |
|---|---:|
| n slides | 54 |
| Original UNI OOF AUC | 0.874 |
| UNI LOTO pooled overall AUC | 0.852 |
| UNI LOTO pooled RAS-like AUC | 0.718 |
| UNI LOTO pooled BRAF-like AUC | 0.855 |
| UNI LOTO Female AUC | 0.862 |
| UNI LOTO Male AUC | 0.875 |

## Per-TSS held-out folds

| tss_held_out   |   n_test |   n_pos |   n_ras_test |   n_ras_pos_test |   tss_holdout_auc |   seconds |
|:---------------|---------:|--------:|-------------:|-----------------:|------------------:|----------:|
| BJ             |       12 |       3 |            1 |                0 |          0.777778 |      21.1 |
| CE             |        4 |       4 |            0 |                0 |        nan        |      21.7 |
| DE             |        3 |       1 |            1 |                0 |          1        |      21.9 |
| DJ             |        7 |       6 |            2 |                2 |          0.833333 |      20.9 |
| EL             |        4 |       4 |            0 |                0 |        nan        |      21.5 |
| EM             |       11 |       0 |           10 |                0 |        nan        |      20.2 |
| ET             |        4 |       2 |            0 |                0 |          0.5      |      21.7 |
| FK             |        3 |       2 |            1 |                1 |          1        |      11.4 |
| FY             |        3 |       0 |            1 |                0 |        nan        |       3.7 |
| J8             |        1 |       1 |            0 |                0 |        nan        |       3.5 |
| KS             |        2 |       2 |            0 |                0 |        nan        |       4.2 |

## Interpretation

This is the center-holdout counterpart to the UNI-final OOF result. It should be cited before making any UNI image-DM1 claim from TCGA.
