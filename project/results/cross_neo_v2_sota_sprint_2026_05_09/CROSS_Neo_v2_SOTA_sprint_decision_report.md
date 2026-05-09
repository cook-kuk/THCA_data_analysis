# CROSS-Neo v2 SOTA Sprint Decision Report

Executive verdict: **PROMOTE_TO_EVALUATION_BENCHMARK_ONLY**

## Best Reviewer-Safe Model

v2_rule_moe_cf_plm_structure across primary locked splits: mean AUPRC=0.519, mean top10=0.550.

## Best Exploratory/Diagnostic Model

| split_name           | model_name                  |    AUPRC |   top10_precision | claim_status   |
|:---------------------|:----------------------------|---------:|------------------:|:---------------|
| source_heldout_CEDAR | v2_plm_lr_frozen_hash_pilot | 0.957745 |                 1 | diagnostic     |

## v2 vs v1

| split                        | anchor_rf                          | v1_best_foldsafe_fusion                                           | v2_best_reviewer_safe                     | v2_beats_anchor   | v2_beats_v1_fusion   |
|:-----------------------------|:-----------------------------------|:------------------------------------------------------------------|:------------------------------------------|:------------------|:---------------------|
| exact_peptide_hla_holdout    | anchor_rf AUPRC=0.525, top10=0.600 | prespecified_rf_qk_no_anchor_w0.5 AUPRC=0.587, top10=0.700        | v2_cf_plm_rf AUPRC=0.561, top10=0.600     | True              | False                |
| near_peptide_cluster_holdout | anchor_rf AUPRC=0.442, top10=0.500 | rule_gate_rf_qk_fallback_train_selected AUPRC=0.519, top10=0.700  | v2_multimodal_lr AUPRC=0.436, top10=0.500 | False             | False                |
| hla_stratified_group_5fold   | anchor_rf AUPRC=0.505, top10=0.600 | nested_rf_qk_quantum_only_train_selected AUPRC=0.555, top10=0.600 | v2_cf_plm_lr AUPRC=0.563, top10=0.500     | True              | True                 |
| hla_supertype_heldout        | anchor_rf AUPRC=0.478, top10=0.500 | prespecified_lr_qk_quantum_only_w0.5 AUPRC=0.545, top10=0.600     | v2_multimodal_lr AUPRC=0.585, top10=0.600 | True              | True                 |

Primary locked splits where v2 beats both anchor_rf and v1 best fold-safe fusion by AUPRC or top10: 2/4.

## Source-Heldout Rescue

| split_name                           | model_name                        |   n |   n_pos |   prevalence |     AUPRC |   top10_precision |   top20_precision | claim_status                  |
|:-------------------------------------|:----------------------------------|----:|--------:|-------------:|----------:|------------------:|------------------:|:------------------------------|
| low_prevalence_stress_split          | v2_multimodal_lr                  | 915 |      41 |    0.0448087 | 0.0615867 |               0.1 |              0.1  | reviewer_safe_internal_locked |
| low_prevalence_stress_split          | v2_cf_plm_lr                      | 915 |      41 |    0.0448087 | 0.0615672 |               0.1 |              0.1  | reviewer_safe_internal_locked |
| low_prevalence_stress_split          | v2_plm_lr_frozen_hash_pilot       | 915 |      41 |    0.0448087 | 0.0605782 |               0.1 |              0.1  | diagnostic                    |
| source_heldout_CEDAR                 | v2_plm_lr_frozen_hash_pilot       | 909 |     851 |    0.936194  | 0.957745  |               1   |              1    | diagnostic                    |
| source_heldout_CEDAR                 | v2_multimodal_lr                  | 909 |     851 |    0.936194  | 0.946896  |               1   |              1    | reviewer_safe_internal_locked |
| source_heldout_CEDAR                 | v2_cf_plm_lr                      | 909 |     851 |    0.936194  | 0.946877  |               1   |              1    | reviewer_safe_internal_locked |
| source_heldout_NEPdb                 | v2_groupdro_proxy_cf_lr           | 572 |     151 |    0.263986  | 0.281349  |               0.5 |              0.45 | reviewer_safe_internal_locked |
| source_heldout_NEPdb                 | source_qk_compact_gamma1          | 886 |     354 |    0.399549  | 0.406426  |               0.4 |              0.4  | diagnostic                    |
| source_heldout_NEPdb                 | v2_multimodal_lr                  | 572 |     151 |    0.263986  | 0.30528   |               0.3 |              0.4  | reviewer_safe_internal_locked |
| source_heldout_TESLA_mmc4            | v2_groupdro_proxy_cf_lr           | 605 |      37 |    0.061157  | 0.077088  |               0.1 |              0.15 | reviewer_safe_internal_locked |
| source_heldout_TESLA_mmc4            | v2_class_balanced_cf_plm_hgb      | 605 |      37 |    0.061157  | 0.0841894 |               0.1 |              0.05 | reviewer_safe_internal_locked |
| source_heldout_TESLA_mmc4            | v2_cf_plm_rf_hla_ranknorm         | 605 |      37 |    0.061157  | 0.0792775 |               0.1 |              0.05 | reviewer_safe_internal_locked |
| source_heldout_TESLA_mmc7_validation | anchor_rf                         | 319 |       6 |    0.0188088 | 0.0378096 |               0   |              0    | reviewer_safe_internal_locked |
| source_heldout_TESLA_mmc7_validation | v2_counterfactual_lr_hla_ranknorm | 310 |       4 |    0.0129032 | 0.0332845 |               0   |              0    | reviewer_safe_internal_locked |
| source_heldout_TESLA_mmc7_validation | v2_cf_plm_rf_hla_ranknorm         | 310 |       4 |    0.0129032 | 0.0277166 |               0   |              0    | reviewer_safe_internal_locked |

Caveat: source-heldout rows are descriptive stress tests. Imported v1 diagnostic rows can have different test-set cardinalities from v2-generated source splits, so do not use them as direct source-level winner claims unless n/test definitions match.

NEPdb nonzero top10: True. TESLA nonzero top20: True.

## Public Overlap

Public comparator training corpora are still unresolved locally; public predictor scores remain excluded from features.

## QK Interpretation

QK remains a bounded fallback/diagnostic branch. It may help individual locked/source splits, but the claim boundary forbids quantum advantage language and requires rescue/harm accounting before promotion.

## Exact Claims Allowed

- CROSS-Neo 2.0 is an internal locked-split prioritization/evaluation framework.
- It performs contamination-aware overlap auditing and source-heldout stress testing.
- Reviewer-safe claims should emphasize AUPRC, top-k precision, enrichment, calibration, and failure modes.

## Exact Claims Forbidden

- External validation.
- Clinical vaccine selection or patient-ready utility.
- Quantum advantage.
- SOTA over public predictors until public overlap and clean comparator benchmarking are resolved.

## Figures Generated

- fig1_problem_and_audit_design.png
- fig2_cross_neo_v2_architecture.png
- fig3_locked_split_performance.png
- fig4_source_heldout_rescue.png
- fig5_public_overlap_audit.png
- fig6_abstention_curve.png
- fig7_case_level_rescue_harm.png
- fig8_claim_boundary_decision_tree.png

## Next 72-Hour Action Plan

1. Replace hashed PLM pilot with ESM2/ProtT5 embeddings under the same split-safe feature interface.
2. Complete local public corpus downloads for MHCflurry/NetMHCpan/BigMHC/PRIME/MixMHCpred/NetMHCstabpan and rerun overlap audit.
3. Run source-heldout rescue with true source-balanced/GroupDRO deep objective if source top-k remains weak.
4. Pre-register wet-lab top-k candidates only after overlap status is assigned.

## Submission Target Recommendation

If public overlap remains unresolved or source-heldout remains unstable, target NeurIPS Datasets & Benchmarks / MLCB benchmark-method framing first. Upgrade to Nature Machine Intelligence/Nature Communications only after independent clean comparator or wet-lab evidence.

## Reproducibility

- Output workbook: `/home/seungho/personal/THCA_data_analysis/project/results/cross_neo_v2_sota_sprint_2026_05_09/CROSS_Neo_v2_all_results_summary.xlsx`
- Environment: `{'python': '3.12.3', 'platform': 'Linux-6.17.0-1011-azure-x86_64-with-glibc2.39', 'seed': 20260509, 'repo': '/home/seungho/personal/THCA_data_analysis'}`
- Full rerun: `bash project/scripts/cross_neo_v2/run_cross_neo_v2_sota_sprint_all.sh`
