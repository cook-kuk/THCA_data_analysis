# ESM2-650m Fast RunPod Evaluation

| split_name                   | model_name                              |   n |   n_pos |   prevalence |    AUPRC |    AUROC |   top10_precision |   top20_precision |
|:-----------------------------|:----------------------------------------|----:|--------:|-------------:|---------:|---------:|------------------:|------------------:|
| exact_peptide_hla_holdout    | v2_multimodal_esm2_650m_lr_fast         |  89 |      21 |     0.235955 | 0.522969 | 0.741597 |               0.7 |              0.5  |
| exact_peptide_hla_holdout    | v2_cf_esm2_650m_lr_fast                 |  89 |      21 |     0.235955 | 0.515206 | 0.730392 |               0.7 |              0.6  |
| exact_peptide_hla_holdout    | v2_groupdro_proxy_cf_esm2_650m_lr_fast  |  89 |      21 |     0.235955 | 0.484957 | 0.717787 |               0.6 |              0.55 |
| hla_stratified_group_5fold   | v2_multimodal_esm2_650m_lr_fast         |  89 |      21 |     0.235955 | 0.476056 | 0.730392 |               0.6 |              0.45 |
| hla_stratified_group_5fold   | v2_cf_esm2_650m_lr_fast                 |  89 |      21 |     0.235955 | 0.46951  | 0.70098  |               0.5 |              0.5  |
| hla_stratified_group_5fold   | v2_groupdro_proxy_cf_esm2_650m_lr_fast  |  89 |      21 |     0.235955 | 0.437989 | 0.685574 |               0.4 |              0.4  |
| hla_supertype_heldout        | v2_multimodal_esm2_650m_lr_fast         |  72 |      15 |     0.208333 | 0.514148 | 0.780117 |               0.6 |              0.45 |
| hla_supertype_heldout        | v2_cf_esm2_650m_lr_fast                 |  72 |      15 |     0.208333 | 0.469068 | 0.729825 |               0.4 |              0.3  |
| hla_supertype_heldout        | v2_source_balanced_cf_esm2_650m_lr_fast |  72 |      15 |     0.208333 | 0.448622 | 0.723977 |               0.4 |              0.4  |
| near_peptide_cluster_holdout | v2_cf_esm2_650m_lr_fast                 |  89 |      21 |     0.235955 | 0.420486 | 0.637255 |               0.7 |              0.4  |
| near_peptide_cluster_holdout | v2_multimodal_esm2_650m_lr_fast         |  89 |      21 |     0.235955 | 0.393627 | 0.632353 |               0.5 |              0.4  |
| near_peptide_cluster_holdout | v2_groupdro_proxy_cf_esm2_650m_lr_fast  |  89 |      21 |     0.235955 | 0.381249 | 0.620448 |               0.6 |              0.4  |

Claim boundary: frozen ESM2-650M features, train-fold-only scaling/model fitting, no public predictor scores.
