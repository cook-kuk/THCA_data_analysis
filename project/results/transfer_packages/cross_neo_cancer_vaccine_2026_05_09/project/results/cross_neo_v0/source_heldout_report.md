# Source-Heldout Stress Test

This is not external validation. It holds out local train-pool sources to test source shortcut sensitivity.

| heldout_study         | method                                   |   n |   n_pos |   prevalence |     AUPRC |    AUROC |    Brier |   top5_precision |   recall_at_5 |   enrichment_at_5 |   top10_precision |   recall_at_10 |   enrichment_at_10 |
|:----------------------|:-----------------------------------------|----:|--------:|-------------:|----------:|---------:|---------:|-----------------:|--------------:|------------------:|------------------:|---------------:|-------------------:|
| CEDAR                 | sourceheld_counterfactual_rf             | 913 |     851 |    0.932092  | 0.932467  | 0.449812 | 0.36596  |              1   |    0.00587544 |          1.07286  |               1   |      0.0117509 |            1.07286 |
| NEPdb                 | sourceheld_counterfactual_rf             | 886 |     354 |    0.399549  | 0.406253  | 0.499825 | 0.261804 |              0   |    0          |          0        |               0   |      0         |            0       |
| TESLA_mmc4            | sourceheld_counterfactual_rf             | 610 |      37 |    0.0606557 | 0.0690419 | 0.567756 | 0.242893 |              0   |    0          |          0        |               0   |      0         |            0       |
| TESLA_mmc7_validation | sourceheld_counterfactual_rf             | 319 |       6 |    0.0188088 | 0.0575599 | 0.747604 | 0.213933 |              0   |    0          |          0        |               0   |      0         |            0       |
| CEDAR                 | sourceheld_prespecified_late_fusion_w0.5 | 913 |     851 |    0.932092  | 0.922543  | 0.444335 | 0.410916 |              1   |    0.00587544 |          1.07286  |               0.9 |      0.0105758 |            0.96557 |
| NEPdb                 | sourceheld_prespecified_late_fusion_w0.5 | 886 |     354 |    0.399549  | 0.394412  | 0.477481 | 0.293803 |              0   |    0          |          0        |               0   |      0         |            0       |
| TESLA_mmc4            | sourceheld_prespecified_late_fusion_w0.5 | 610 |      37 |    0.0606557 | 0.0583761 | 0.472006 | 0.311812 |              0   |    0          |          0        |               0   |      0         |            0       |
| TESLA_mmc7_validation | sourceheld_prespecified_late_fusion_w0.5 | 319 |       6 |    0.0188088 | 0.0174647 | 0.317359 | 0.24176  |              0   |    0          |          0        |               0   |      0         |            0       |
| CEDAR                 | sourceheld_qk_compact_gamma1             | 913 |     851 |    0.932092  | 0.926744  | 0.461317 | 0.527638 |              0.8 |    0.00470035 |          0.858284 |               0.9 |      0.0105758 |            0.96557 |
| NEPdb                 | sourceheld_qk_compact_gamma1             | 886 |     354 |    0.399549  | 0.406426  | 0.482637 | 0.401408 |              0.4 |    0.00564972 |          1.00113  |               0.4 |      0.0112994 |            1.00113 |
| TESLA_mmc4            | sourceheld_qk_compact_gamma1             | 610 |      37 |    0.0606557 | 0.0516625 | 0.437904 | 0.476714 |              0   |    0          |          0        |               0   |      0         |            0       |
| TESLA_mmc7_validation | sourceheld_qk_compact_gamma1             | 319 |       6 |    0.0188088 | 0.0145965 | 0.180511 | 0.350284 |              0   |    0          |          0        |               0   |      0         |            0       |
