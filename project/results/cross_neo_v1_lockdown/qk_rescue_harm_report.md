# CROSS-Neo v1 QK Rescue/Harm Report

Definitions: a rescued positive is a positive outside anchor top-10 or below the fold median under `anchor_rf` but top-10 under QK/fusion. A harmed negative is a negative outside anchor top-10 but moved into top-10 by QK/fusion.

## Summary

| split_name                           | comparator               | event            |   n |   median_rank_pct_delta |   median_score_delta |
|:-------------------------------------|:-------------------------|:-----------------|----:|------------------------:|---------------------:|
| exact_peptide_hla_holdout            | best_foldsafe_fusion     | harmed_negative  |   9 |              0.222222   |            0.131843  |
| exact_peptide_hla_holdout            | best_foldsafe_fusion     | rescued_positive |   2 |              0.401961   |            0.13648   |
| exact_peptide_hla_holdout            | best_foldsafe_fusion     | stable_positive  |  15 |              0.0555556  |            0.069784  |
| exact_peptide_hla_holdout            | best_foldsafe_fusion     | unstable_case    |  15 |             -0.277778   |           -0.0772914 |
| exact_peptide_hla_holdout            | qk_no_anchor             | harmed_negative  |  13 |              0.444444   |            0.231129  |
| exact_peptide_hla_holdout            | qk_no_anchor             | rescued_positive |   4 |              0.611111   |            0.280836  |
| exact_peptide_hla_holdout            | qk_no_anchor             | stable_positive  |  11 |              0.0588235  |            0.249472  |
| exact_peptide_hla_holdout            | qk_no_anchor             | unstable_case    |  26 |             -0.352941   |           -0.155581  |
| exact_peptide_hla_holdout            | qk_quantum_only          | harmed_negative  |  14 |              0.527778   |            0.310198  |
| exact_peptide_hla_holdout            | qk_quantum_only          | rescued_positive |   4 |              0.602941   |            0.384949  |
| exact_peptide_hla_holdout            | qk_quantum_only          | stable_positive  |  12 |             -0.112745   |            0.210974  |
| exact_peptide_hla_holdout            | qk_quantum_only          | unstable_case    |  23 |             -0.333333   |           -0.134495  |
| hla_stratified_group_5fold           | best_foldsafe_fusion     | harmed_negative  |   5 |              0.263158   |            0.0794904 |
| hla_stratified_group_5fold           | best_foldsafe_fusion     | rescued_positive |   6 |              0.316667   |            0.120167  |
| hla_stratified_group_5fold           | best_foldsafe_fusion     | stable_positive  |  11 |              0          |            0.0435041 |
| hla_stratified_group_5fold           | best_foldsafe_fusion     | unstable_case    |  15 |             -0.272727   |           -0.0259041 |
| hla_stratified_group_5fold           | qk_no_anchor             | harmed_negative  |  11 |              0.421053   |            0.226074  |
| hla_stratified_group_5fold           | qk_no_anchor             | rescued_positive |   7 |              0.4        |            0.332158  |
| hla_stratified_group_5fold           | qk_no_anchor             | stable_positive  |   9 |             -0.0909091  |            0.0849248 |
| hla_stratified_group_5fold           | qk_no_anchor             | unstable_case    |  32 |             -0.324561   |           -0.0538215 |
| hla_stratified_group_5fold           | qk_quantum_only          | harmed_negative  |  15 |              0.5        |            0.40798   |
| hla_stratified_group_5fold           | qk_quantum_only          | rescued_positive |   5 |              0.5        |            0.485559  |
| hla_stratified_group_5fold           | qk_quantum_only          | stable_positive  |   9 |             -0.030303   |            0.175521  |
| hla_stratified_group_5fold           | qk_quantum_only          | unstable_case    |  26 |             -0.366029   |           -0.0672937 |
| hla_supertype_heldout                | best_foldsafe_fusion     | harmed_negative  |   8 |              0.357895   |            0.209649  |
| hla_supertype_heldout                | best_foldsafe_fusion     | rescued_positive |   3 |              0.368421   |            0.337076  |
| hla_supertype_heldout                | best_foldsafe_fusion     | stable_positive  |   8 |              0.0818182  |            0.167007  |
| hla_supertype_heldout                | best_foldsafe_fusion     | unstable_case    |  19 |             -0.315789   |           -0.137923  |
| hla_supertype_heldout                | qk_no_anchor             | harmed_negative  |  10 |              0.394737   |            0.28472   |
| hla_supertype_heldout                | qk_no_anchor             | rescued_positive |   4 |              0.409091   |            0.362047  |
| hla_supertype_heldout                | qk_no_anchor             | stable_positive  |   6 |             -0.127273   |            0.106071  |
| hla_supertype_heldout                | qk_no_anchor             | unstable_case    |  27 |             -0.333333   |           -0.050262  |
| hla_supertype_heldout                | qk_quantum_only          | harmed_negative  |  14 |              0.533333   |            0.424556  |
| hla_supertype_heldout                | qk_quantum_only          | rescued_positive |   2 |              0.475279   |            0.480668  |
| hla_supertype_heldout                | qk_quantum_only          | stable_positive  |   5 |              0.0909091  |            0.317224  |
| hla_supertype_heldout                | qk_quantum_only          | unstable_case    |  26 |             -0.366029   |           -0.052836  |
| near_peptide_cluster_holdout         | best_foldsafe_fusion     | harmed_negative  |   9 |              0.222222   |            0.0868855 |
| near_peptide_cluster_holdout         | best_foldsafe_fusion     | rescued_positive |   2 |              0.205882   |            0.120236  |
| near_peptide_cluster_holdout         | best_foldsafe_fusion     | stable_positive  |  13 |              0.0588235  |            0.0395243 |
| near_peptide_cluster_holdout         | best_foldsafe_fusion     | unstable_case    |   6 |             -0.270468   |           -0.0133931 |
| near_peptide_cluster_holdout         | qk_no_anchor             | harmed_negative  |  13 |              0.529412   |            0.344629  |
| near_peptide_cluster_holdout         | qk_no_anchor             | rescued_positive |   5 |              0.5        |            0.305568  |
| near_peptide_cluster_holdout         | qk_no_anchor             | stable_positive  |  11 |              0.0555556  |            0.235604  |
| near_peptide_cluster_holdout         | qk_no_anchor             | unstable_case    |  21 |             -0.444444   |           -0.115829  |
| near_peptide_cluster_holdout         | qk_quantum_only          | harmed_negative  |  15 |              0.411765   |            0.324983  |
| near_peptide_cluster_holdout         | qk_quantum_only          | rescued_positive |   4 |              0.5        |            0.449458  |
| near_peptide_cluster_holdout         | qk_quantum_only          | stable_positive  |  12 |             -0.0820433  |            0.18114   |
| near_peptide_cluster_holdout         | qk_quantum_only          | unstable_case    |  22 |             -0.428105   |           -0.129272  |
| repeated_stratified_5x5_internal     | best_foldsafe_fusion     | harmed_negative  |  60 |              0.400327   |            0.149378  |
| repeated_stratified_5x5_internal     | best_foldsafe_fusion     | rescued_positive |  10 |              0.277778   |            0.111633  |
| repeated_stratified_5x5_internal     | best_foldsafe_fusion     | stable_positive  |  71 |              0          |            0.135029  |
| repeated_stratified_5x5_internal     | best_foldsafe_fusion     | unstable_case    |  85 |             -0.333333   |           -0.172098  |
| repeated_stratified_5x5_internal     | qk_no_anchor             | harmed_negative  |  87 |              0.5        |            0.29718   |
| repeated_stratified_5x5_internal     | qk_no_anchor             | rescued_positive |  15 |              0.444444   |            0.302141  |
| repeated_stratified_5x5_internal     | qk_no_anchor             | stable_positive  |  59 |              0          |            0.164214  |
| repeated_stratified_5x5_internal     | qk_no_anchor             | unstable_case    | 133 |             -0.444444   |           -0.146911  |
| repeated_stratified_5x5_internal     | qk_quantum_only          | harmed_negative  |  88 |              0.485294   |            0.296673  |
| repeated_stratified_5x5_internal     | qk_quantum_only          | rescued_positive |  16 |              0.5        |            0.403318  |
| repeated_stratified_5x5_internal     | qk_quantum_only          | stable_positive  |  52 |              0.0555556  |            0.204058  |
| repeated_stratified_5x5_internal     | qk_quantum_only          | unstable_case    | 139 |             -0.444444   |           -0.121922  |
| source_heldout_CEDAR                 | best_foldsafe_fusion     | rescued_positive |   8 |              0.0671067  |            0.275761  |
| source_heldout_CEDAR                 | best_foldsafe_fusion     | stable_positive  |   2 |              0.00110011 |            0.23627   |
| source_heldout_CEDAR                 | best_foldsafe_fusion     | unstable_case    | 516 |             -0.259626   |           -0.0773464 |
| source_heldout_CEDAR                 | source_qk_compact_gamma1 | harmed_negative  |   1 |              0.706271   |            0.641855  |
| source_heldout_CEDAR                 | source_qk_compact_gamma1 | rescued_positive |   9 |              0.792079   |            0.656444  |
| source_heldout_CEDAR                 | source_qk_compact_gamma1 | unstable_case    | 513 |             -0.253025   |           -0.180895  |
| source_heldout_NEPdb                 | best_foldsafe_fusion     | harmed_negative  |   7 |              0.0340909  |            0.150072  |
| source_heldout_NEPdb                 | best_foldsafe_fusion     | unstable_case    | 318 |             -0.26792    |           -0.102338  |
| source_heldout_NEPdb                 | source_qk_compact_gamma1 | harmed_negative  |   6 |              0.787587   |            0.599002  |
| source_heldout_NEPdb                 | source_qk_compact_gamma1 | rescued_positive |   4 |              0.524476   |            0.524715  |
| source_heldout_NEPdb                 | source_qk_compact_gamma1 | unstable_case    | 328 |             -0.267045   |           -0.125178  |
| source_heldout_TESLA_mmc4            | best_foldsafe_fusion     | harmed_negative  |   6 |              0.0297521  |            0.197374  |
| source_heldout_TESLA_mmc4            | best_foldsafe_fusion     | unstable_case    | 344 |              0.255372   |            0.0471736 |
| source_heldout_TESLA_mmc4            | source_qk_compact_gamma1 | harmed_negative  |  10 |              0.42562    |            0.501623  |
| source_heldout_TESLA_mmc4            | source_qk_compact_gamma1 | unstable_case    | 354 |             -0.254545   |            0.103453  |
| source_heldout_TESLA_mmc7_validation | best_foldsafe_fusion     | harmed_negative  |   9 |              0.0741935  |            0.235223  |
| source_heldout_TESLA_mmc7_validation | best_foldsafe_fusion     | unstable_case    | 153 |             -0.251613   |           -0.0826962 |
| source_heldout_TESLA_mmc7_validation | source_qk_compact_gamma1 | harmed_negative  |   9 |              0.351613   |            0.523996  |
| source_heldout_TESLA_mmc7_validation | source_qk_compact_gamma1 | unstable_case    | 152 |             -0.254839   |           -0.209717  |
| study_heldout                        | best_foldsafe_fusion     | rescued_positive |   3 |              0          |           -0.0193685 |
| study_heldout                        | best_foldsafe_fusion     | stable_positive  |   2 |              0          |            0.0549896 |
| study_heldout                        | best_foldsafe_fusion     | unstable_case    |   1 |              0.333333   |            0.0659464 |
| study_heldout                        | qk_no_anchor             | rescued_positive |   3 |              0.111111   |            0.0665052 |
| study_heldout                        | qk_no_anchor             | stable_positive  |   2 |             -0.333333   |            0.0347745 |
| study_heldout                        | qk_no_anchor             | unstable_case    |   3 |              0.333333   |            0.354709  |
| study_heldout                        | qk_quantum_only          | rescued_positive |   3 |              0.222222   |            0.159083  |
| study_heldout                        | qk_quantum_only          | stable_positive  |   2 |             -0.444444   |           -0.110197  |
| study_heldout                        | qk_quantum_only          | unstable_case    |   3 |              0.333333   |            0.259191  |

## Interpretation

- QK remains a bounded fallback/fusion component, not a standalone main claim.
- Harmed-negative rows should be inspected before promoting any QK-heavy fusion.
