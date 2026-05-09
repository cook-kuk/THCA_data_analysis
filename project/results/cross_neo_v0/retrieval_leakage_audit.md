# Retrieval Leakage Audit

Retrieval features were computed using only each fold's strict-set train rows.

| split_name                       | retrieval_leakage_risk   |   n |
|:---------------------------------|:-------------------------|----:|
| exact_peptide_hla_holdout        | clean_no_reference       |  87 |
| exact_peptide_hla_holdout        | near_hit                 |   2 |
| hla_stratified_group_5fold       | clean_no_reference       |  89 |
| hla_supertype_heldout            | clean_no_reference       |  72 |
| near_peptide_cluster_holdout     | clean_no_reference       |  89 |
| repeated_stratified_5x5_internal | clean_no_reference       | 437 |
| repeated_stratified_5x5_internal | near_hit                 |   8 |
| study_heldout                    | clean_no_reference       |  89 |
