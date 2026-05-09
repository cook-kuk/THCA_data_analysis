# ESM2-35m Fast RunPod Evaluation

| split_name                   | model_name                             |   n |   n_pos |   prevalence |    AUPRC |    AUROC |   top10_precision |   top20_precision |
|:-----------------------------|:---------------------------------------|----:|--------:|-------------:|---------:|---------:|------------------:|------------------:|
| exact_peptide_hla_holdout    | v2_multimodal_esm2_35m_lr_fast         |  89 |      21 |     0.235955 | 0.505271 | 0.721989 |               0.6 |              0.65 |
| exact_peptide_hla_holdout    | v2_esm2_35m_lr_fast                    |  89 |      21 |     0.235955 | 0.494123 | 0.677871 |               0.5 |              0.55 |
| exact_peptide_hla_holdout    | v2_cf_esm2_35m_lr_fast                 |  89 |      21 |     0.235955 | 0.476705 | 0.70098  |               0.6 |              0.6  |
| hla_stratified_group_5fold   | v2_multimodal_esm2_35m_lr_fast         |  89 |      21 |     0.235955 | 0.466458 | 0.716387 |               0.6 |              0.55 |
| hla_stratified_group_5fold   | v2_groupdro_proxy_cf_esm2_35m_lr_fast  |  89 |      21 |     0.235955 | 0.454769 | 0.668768 |               0.4 |              0.3  |
| hla_stratified_group_5fold   | v2_source_balanced_cf_esm2_35m_lr_fast |  89 |      21 |     0.235955 | 0.439766 | 0.663866 |               0.4 |              0.3  |
| hla_supertype_heldout        | v2_multimodal_esm2_35m_lr_fast         |  72 |      15 |     0.208333 | 0.450586 | 0.74152  |               0.5 |              0.45 |
| hla_supertype_heldout        | v2_cf_esm2_35m_lr_fast                 |  72 |      15 |     0.208333 | 0.435867 | 0.688889 |               0.3 |              0.35 |
| hla_supertype_heldout        | v2_groupdro_proxy_cf_esm2_35m_lr_fast  |  72 |      15 |     0.208333 | 0.410468 | 0.678363 |               0.3 |              0.35 |
| near_peptide_cluster_holdout | v2_cf_esm2_35m_lr_fast                 |  89 |      21 |     0.235955 | 0.401278 | 0.584734 |               0.6 |              0.45 |
| near_peptide_cluster_holdout | v2_multimodal_esm2_35m_lr_fast         |  89 |      21 |     0.235955 | 0.401144 | 0.605042 |               0.5 |              0.35 |
| near_peptide_cluster_holdout | v2_esm2_35m_lr_fast                    |  89 |      21 |     0.235955 | 0.364912 | 0.565126 |               0.3 |              0.2  |

Claim boundary: frozen ESM2-35M features, train-fold-only scaling/model fitting, no public predictor scores.
