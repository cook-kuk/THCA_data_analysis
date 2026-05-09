# Neoantigen All-Algorithm / All-Testset Fair Comparison

Workbook: `project/results/cross_neo_v1_lockdown/NEOANTIGEN_ALL_ALGORITHMS_ALL_TESTSETS_FAIR_COMPARISON.xlsx`

This workbook is strict: it includes legacy/exploratory/public-comparator rows, but flags claim status and leakage-control tier so they are not mixed with reviewer-safe locked rows.

## Counts

- metric rows raw: 9240
- metric rows deduplicated for ranking: 4351
- source tables scanned: 81
- loaded source tables: 80

## Best Primary Locked Rows

| testset                      | algorithm                                           |   n |   n_pos |   prevalence |    AUPRC |    AUROC |   top5_precision |   top10_precision | claim_status                  |
|:-----------------------------|:----------------------------------------------------|----:|--------:|-------------:|---------:|---------:|-----------------:|------------------:|:------------------------------|
| exact_peptide_hla_holdout    | rule_gated_C_QK_structure                           |  89 |      21 |     0.235955 | 0.647073 | 0.778711 |              1   |               0.7 | diagnostic                    |
| exact_peptide_hla_holdout    | rule_gated_C_QK_structure                           |  88 |      20 |     0.227273 | 0.623016 | 0.768382 |              1   |               0.6 | diagnostic                    |
| exact_peptide_hla_holdout    | rule_gated_C_QK_structure                           |  87 |      19 |     0.218391 | 0.622607 | 0.772446 |              0.8 |               0.6 | diagnostic                    |
| hla_stratified_group_5fold   | nested_rf_qk_quantum_only_train_selected            |  89 |      21 |     0.235955 | 0.555039 | 0.761204 |              0.6 |               0.6 | reviewer_safe_internal_locked |
| hla_stratified_group_5fold   | counterfactual_rf_plus_qk_quantum_only_gamma1_w0.75 |  89 |      21 |     0.235955 | 0.555039 | 0.761204 |              0.6 |               0.6 | exploratory_descriptive_only  |
| hla_stratified_group_5fold   | prespecified_equal_weight_C_QK_no_anchor            |  88 |      20 |     0.227273 | 0.537125 | 0.757353 |              0.8 |               0.5 | diagnostic                    |
| hla_supertype_heldout        | nested_learned_gate_C_QK_structure                  |  72 |      15 |     0.208333 | 0.611634 | 0.816374 |              0.8 |               0.5 | diagnostic                    |
| hla_supertype_heldout        | counterfactual_rf_plus_qk_quantum_only_gamma1_w0.75 |  72 |      15 |     0.208333 | 0.546529 | 0.768421 |              0.6 |               0.6 | exploratory_descriptive_only  |
| hla_supertype_heldout        | prespecified_lr_qk_quantum_only_w0.5                |  72 |      15 |     0.208333 | 0.544922 | 0.767251 |              0.6 |               0.6 | reviewer_safe_internal_locked |
| near_peptide_cluster_holdout | counterfactual_rf_plus_qk_quantum_only_gamma1_w0.25 |  89 |      21 |     0.235955 | 0.518892 | 0.710084 |              0.8 |               0.7 | exploratory_descriptive_only  |
| near_peptide_cluster_holdout | rule_gate_rf_qk_fallback_train_selected             |  89 |      21 |     0.235955 | 0.518669 | 0.722689 |              0.8 |               0.7 | reviewer_safe_internal_locked |
| near_peptide_cluster_holdout | counterfactual_rf_plus_qk_quantum_only_gamma1_w0.75 |  89 |      21 |     0.235955 | 0.518669 | 0.722689 |              0.8 |               0.7 | exploratory_descriptive_only  |

## Use Rules

- Public pretrained predictors are comparator-only unless row-level public overlap is resolved.
- QK standalone rows are diagnostic/fallback only; no quantum advantage claim.
- Source-heldout failures remain visible and should not be filtered out.
- Primary ranking should use AUPRC and top-k before AUROC.
