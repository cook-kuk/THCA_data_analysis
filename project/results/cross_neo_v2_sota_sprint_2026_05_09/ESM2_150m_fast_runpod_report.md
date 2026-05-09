# ESM2-150m Fast RunPod Evaluation

| split_name                   | model_name                              |   n |   n_pos |   prevalence |    AUPRC |    AUROC |   top10_precision |   top20_precision |
|:-----------------------------|:----------------------------------------|----:|--------:|-------------:|---------:|---------:|------------------:|------------------:|
| exact_peptide_hla_holdout    | v2_cf_esm2_150m_lr_fast                 |  89 |      21 |     0.235955 | 0.540566 | 0.735294 |               0.8 |              0.5  |
| exact_peptide_hla_holdout    | v2_multimodal_esm2_150m_lr_fast         |  89 |      21 |     0.235955 | 0.525202 | 0.741597 |               0.7 |              0.5  |
| exact_peptide_hla_holdout    | v2_source_balanced_cf_esm2_150m_lr_fast |  89 |      21 |     0.235955 | 0.513314 | 0.72479  |               0.7 |              0.45 |
| hla_stratified_group_5fold   | v2_multimodal_esm2_150m_lr_fast         |  89 |      21 |     0.235955 | 0.475045 | 0.742997 |               0.5 |              0.5  |
| hla_stratified_group_5fold   | v2_cf_esm2_150m_lr_fast                 |  89 |      21 |     0.235955 | 0.471511 | 0.72479  |               0.5 |              0.45 |
| hla_stratified_group_5fold   | v2_groupdro_proxy_cf_esm2_150m_lr_fast  |  89 |      21 |     0.235955 | 0.447968 | 0.701681 |               0.4 |              0.4  |
| hla_supertype_heldout        | v2_cf_esm2_150m_lr_fast                 |  72 |      15 |     0.208333 | 0.48111  | 0.71345  |               0.4 |              0.4  |
| hla_supertype_heldout        | v2_multimodal_esm2_150m_lr_fast         |  72 |      15 |     0.208333 | 0.470509 | 0.74269  |               0.5 |              0.4  |
| hla_supertype_heldout        | v2_groupdro_proxy_cf_esm2_150m_lr_fast  |  72 |      15 |     0.208333 | 0.456731 | 0.700585 |               0.4 |              0.35 |
| near_peptide_cluster_holdout | v2_cf_esm2_150m_lr_fast                 |  89 |      21 |     0.235955 | 0.365435 | 0.590336 |               0.6 |              0.4  |
| near_peptide_cluster_holdout | v2_multimodal_esm2_150m_lr_fast         |  89 |      21 |     0.235955 | 0.363116 | 0.602241 |               0.5 |              0.35 |
| near_peptide_cluster_holdout | v2_esm2_150m_lr_fast                    |  89 |      21 |     0.235955 | 0.356577 | 0.582633 |               0.3 |              0.35 |

Claim boundary: frozen ESM2-150M features, train-fold-only scaling/model fitting, no public predictor scores.
