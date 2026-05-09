# CROSS-Neo v1 Lockdown Decision Report

This is an internal locked-split and source-heldout stress report. It does not claim external validation or quantum advantage.

## A. Executive Decision

Decision: **HOLD**.

Fold-safe fusion improves AUPRC or top-10 over `anchor_rf` on 4 of 4 primary locked splits. Source-heldout collapse remains: True. Public overlap unresolved: True. Fusion weight instability flagged: True.

## B. Best Conservative Candidate

`nested_rf_qk_quantum_only_train_selected` with mean primary AUPRC=0.519, mean top10=0.625, split coverage=4.

Conservative candidates were restricted to `anchor_rf`, `anchor_lr`, prespecified w0.5 fusion, nested fold-safe fusion, rule-gated fold-safe fusion, and source top-k rescue variants that do not use heldout labels for fitting.

## C. Best Exploratory Candidate

| split_name                 | fusion                                              | status      |   n |   n_pos |   prevalence |    AUPRC |    AUROC |    Brier |   top5_precision |   top10_precision |
|:---------------------------|:----------------------------------------------------|:------------|----:|--------:|-------------:|---------:|---------:|---------:|-----------------:|------------------:|
| study_heldout              | counterfactual_rf_plus_qk_no_anchor_gamma1_w0.75    | exploratory |   9 |       5 |     0.555556 | 0.644444 | 0.4      | 0.314121 |              0.4 |          0.555556 |
| exact_peptide_hla_holdout  | counterfactual_rf_plus_qk_no_anchor_gamma1_w0.75    | exploratory |  89 |      21 |     0.235955 | 0.596809 | 0.768207 | 0.178298 |              0.8 |          0.7      |
| study_heldout              | counterfactual_rf_plus_qk_quantum_only_gamma1_w0.75 | exploratory |   9 |       5 |     0.555556 | 0.591111 | 0.3      | 0.312124 |              0.4 |          0.555556 |
| study_heldout              | counterfactual_rf_plus_qk_clean_no_tcr_gamma1_w0.75 | exploratory |   9 |       5 |     0.555556 | 0.591111 | 0.3      | 0.317647 |              0.4 |          0.555556 |
| hla_stratified_group_5fold | counterfactual_rf_plus_qk_quantum_only_gamma1_w0.75 | exploratory |  89 |      21 |     0.235955 | 0.555039 | 0.761204 | 0.18385  |              0.6 |          0.6      |
| exact_peptide_hla_holdout  | counterfactual_rf_plus_qk_quantum_only_gamma1_w0.75 | exploratory |  89 |      21 |     0.235955 | 0.553634 | 0.766106 | 0.181363 |              0.8 |          0.8      |
| hla_supertype_heldout      | counterfactual_rf_plus_qk_quantum_only_gamma1_w0.75 | exploratory |  72 |      15 |     0.208333 | 0.546529 | 0.768421 | 0.178696 |              0.6 |          0.6      |
| exact_peptide_hla_holdout  | counterfactual_rf_plus_qk_clean_no_tcr_gamma1_w0.75 | exploratory |  89 |      21 |     0.235955 | 0.5452   | 0.757703 | 0.179308 |              0.8 |          0.6      |

Exploratory rows are descriptive only; fixed w0.75 is not promoted unless selected by inner train folds.

## D. Whether v1 Beats v0 Anchor

| split_name                   | method                                   |    AUPRC |   anchor_AUPRC |   top10_precision |   anchor_top10 | improves_AUPRC   | improves_top10   | passes   |
|:-----------------------------|:-----------------------------------------|---------:|---------------:|------------------:|---------------:|:-----------------|:-----------------|:---------|
| exact_peptide_hla_holdout    | prespecified_rf_qk_no_anchor_w0.5        | 0.586965 |       0.525436 |               0.7 |            0.6 | True             | True             | True     |
| hla_stratified_group_5fold   | nested_rf_qk_quantum_only_train_selected | 0.555039 |       0.50549  |               0.6 |            0.6 | True             | False            | True     |
| hla_supertype_heldout        | prespecified_lr_qk_quantum_only_w0.5     | 0.544922 |       0.477884 |               0.6 |            0.5 | True             | True             | True     |
| near_peptide_cluster_holdout | rule_gate_rf_qk_fallback_train_selected  | 0.518669 |       0.441926 |               0.7 |            0.5 | True             | True             | True     |

## E. Why Source-Heldout Collapses

| heldout_study         |   n |   positives |   prevalence |   positive_top10 |     AUPRC |    AUROC |   top10_precision |   median_positive_rank | cause_prevalence_too_low   | cause_source_shift   | cause_score_miscalibration   | cause_positive_low_tail   | cause_label_definition_mismatch_possible   |
|:----------------------|----:|------------:|-------------:|-----------------:|----------:|---------:|------------------:|-----------------------:|:---------------------------|:---------------------|:-----------------------------|:--------------------------|:-------------------------------------------|
| CEDAR                 | 913 |         851 |    0.932092  |               10 | 0.921816  | 0.439028 |                 1 |                  461   | False                      | True                 | True                         | True                      | False                                      |
| NEPdb                 | 886 |         354 |    0.399549  |                0 | 0.395509  | 0.479297 |                 0 |                  439.5 | False                      | True                 | True                         | False                     | True                                       |
| TESLA_mmc4            | 610 |          37 |    0.0606557 |                0 | 0.0551107 | 0.464931 |                 0 |                  360   | True                       | True                 | True                         | True                      | True                                       |
| TESLA_mmc7_validation | 319 |           6 |    0.0188088 |                0 | 0.0168969 | 0.294462 |                 0 |                  237.5 | True                       | True                 | True                         | True                      | True                                       |

The collapse is most consistent with source distribution/label-definition shift plus low-prevalence top-k brittleness in TESLA. NEPdb/TESLA positives often land outside the high-score tail; public overlap remains unresolved.

## F. QK Interpretation

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

QK should remain a bounded fallback/fusion component. It is useful in selected locked splits but can move negatives into top-k, so QK-heavy claims are not reviewer-safe.

## G. Public Overlap

| source        | status                    | local_file_found_or_missing                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        | action_needed                            |
|:--------------|:--------------------------|:-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|:-----------------------------------------|
| MHCflurry     | local_file_found_unparsed | project/data/external_benchmarks/ITSNdb/Colab/MHCFlurry.py;project/data/external_benchmarks/ITSNdb/Colab/MHCFlurry_Colab.ipynb;project/results/p_neo_bayesian_2026_05_09/wave3_algorithm_sweep/mhcflurry_input.csv;project/results/p_neo_bayesian_2026_05_09/wave3_algorithm_sweep/02_run_mhcflurry.py;project/results/p_neo_bayesian_2026_05_09/wave3_algorithm_sweep/mhcflurry_scored.tsv;project/results/p_neo_bayesian_2026_05_09/wave5b/mhcflurry_raw.parquet                                                 | download/parse official training corpus  |
| NetMHCpan     | local_file_found_unparsed | project/data/external_benchmarks/ITSNdb/R/netMHCpan.R;project/data/external_benchmarks/ITSNdb/R/netMHCpan_cohorts.R;project/data/external_benchmarks/ITSNdb/man/Install_netMHCPan.Rd;project/data/external_benchmarks/ITSNdb/man/RunNetMHCPan.Rd;project/data/external_benchmarks/ITSNdb/man/RunNetMHCPan_peptides.Rd;project/results/p_neo_bayesian_2026_05_09/wave11/wave11_predictions/NetMHCpan_4.1__itsndb.tsv                                                                                                | download/parse official training corpus  |
| BigMHC        | local_file_found_unparsed | project/results/p_neo_bayesian_2026_05_09/wave3_algorithm_sweep/bigmhc_input.csv;project/results/p_neo_bayesian_2026_05_09/wave3_algorithm_sweep/03_run_bigmhc.py;project/results/p_neo_bayesian_2026_05_09/wave3_algorithm_sweep/bigmhc_im_output.csv;project/results/p_neo_bayesian_2026_05_09/wave3_algorithm_sweep/bigmhc_scored.tsv;project/results/p_neo_bayesian_2026_05_09/wave3_algorithm_sweep/03b_run_bigmhc_el.py;project/results/p_neo_bayesian_2026_05_09/wave3_algorithm_sweep/bigmhc_el_output.csv | download/parse official training corpus  |
| PRIME         | local_file_found_unparsed | project/data/external_benchmarks/ITSNdb/man/Install_PRIME.Rd;project/data/external_benchmarks/ITSNdb/man/RunPRIME.Rd;project/results/hla_deepdive_2026_05_08/track2_k2_baseline/figures/F14_LD_Dprime_heatmaps.pdf;project/results/hla_deepdive_2026_05_08/track2_k2_baseline/figures/F14_LD_Dprime_heatmaps.png;project/results/p_neo_bayesian_2026_05_09/wave3_algorithm_sweep/prime_allele_index.tsv;project/results/p_neo_bayesian_2026_05_09/wave3_algorithm_sweep/04_run_prime.py                            | download/parse official training corpus  |
| MixMHCpred    | local_file_found_unparsed | project/data/external_benchmarks/ITSNdb/R/mixMHCpred.R;project/data/external_benchmarks/ITSNdb/man/RunMixMHCpred.Rd                                                                                                                                                                                                                                                                                                                                                                                                | download/parse official training corpus  |
| NetMHCstabpan | local_file_found_unparsed | project/results/p_neo_bayesian_2026_05_09/wave9/run_netmhcstabpan.py;project/results/p_neo_bayesian_2026_05_09/wave9/predictions_netmhcstabpan.tsv;project/results/p_neo_bayesian_2026_05_09/wave9/auroc_netmhcstabpan.tsv;project/results/p_neo_bayesian_2026_05_09/wave9/_meta_netmhcstabpan.json                                                                                                                                                                                                                | download/parse official training corpus  |
| IEDB          | local_audit_available     | project/data/external_benchmarks/ITSNdb/R/IEDB_CIImm.R                                                                                                                                                                                                                                                                                                                                                                                                                                                             | row-level audit completed for local rows |
| CEDAR         | local_audit_available     | project/results/cross_neo_v0/master_table.tsv                                                                                                                                                                                                                                                                                                                                                                                                                                                                      | row-level audit completed for local rows |
| NEPdb         | local_audit_available     | project/results/cross_neo_v0/master_table.tsv                                                                                                                                                                                                                                                                                                                                                                                                                                                                      | row-level audit completed for local rows |
| TESLA         | local_audit_available     | project/results/cross_neo_v0/master_table.tsv                                                                                                                                                                                                                                                                                                                                                                                                                                                                      | row-level audit completed for local rows |

Unresolved public comparators/downloads:

| source        | expected_url_or_search_term                                   | local_target_path                             |
|:--------------|:--------------------------------------------------------------|:----------------------------------------------|
| MHCflurry     | MHCflurry training data presentation affinity IEDB MS ligands | project/data/external/mhcflurry_training/     |
| NetMHCpan     | NetMHCpan 4.1 training data BA EL peptide HLA                 | project/data/external/netmhcpan_training/     |
| BigMHC        | BigMHC Mendeley el_train im_train peptide HLA immunogenicity  | project/data/external/bigmhc/                 |
| PRIME         | PRIME neoantigen immunogenicity training data peptide HLA     | project/data/external/prime/                  |
| MixMHCpred    | MixMHCpred training ligands peptide HLA allele                | project/data/external/mixmhcpred/             |
| NetMHCstabpan | NetMHCstabpan training data peptide HLA stability             | project/data/external/netmhcstabpan_training/ |

## H. Reviewer-Safe Claim

"CROSS-Neo v1 is an internal locked-split prioritization framework that combines mutant-WT counterfactual encoding with fold-safe quantum-kernel fallback. It improves internal/HLA-stratified ranking and top-k enrichment but remains source-shift limited; therefore, it is presented as a stress-tested candidate rather than an externally validated neoantigen predictor."

## I. Forbidden Claims

- No external validation.
- No quantum advantage.
- No clinical vaccine selection claim.
- No public pretrained comparator cleanliness claim unless the overlap audit resolves it.

## J. Next Experimental Validation

- Independent external time split.
- Peptide-HLA near-neighbor clean benchmark.
- Wet-lab or literature-backed top-k candidate validation.
- Public overlap audit completion for MHCflurry, NetMHCpan, BigMHC, PRIME, MixMHCpred, NetMHCstabpan, and IEDB.
