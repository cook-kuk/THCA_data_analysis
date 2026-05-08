# Spatial CNV reviewer stress battery - 2026-05-08

Purpose: attack the spatial CNV ecosystem claim with slide-level, QC/composition, condition-stratified, spatial-block, and influence analyses.

## Verdict

### arm_level_expression_cnv
- QC/composition predicts territory within slide: median apparent AUROC **0.97**.
- QC/composition leave-one-slide-out AUROC: **0.601**.
- TROP2_raw: **fails_or_reverses**; median slide beta -0.00349, positive slides 6/12, condition+ 2/3, block p 0.00249, LOSO sign changes 0.
- TACSTD2_z: **directionally_inconsistent**; median slide beta 0.00982, positive slides 7/12, condition+ 1/3, block p 0.00249, LOSO sign changes 0.
- Macrophage_TAM_z: **survives_strong**; median slide beta 0.107, positive slides 10/12, condition+ 3/3, block p 0.00249, LOSO sign changes 0.
- Tcell_cytotoxic_z: **survives_moderate**; median slide beta 0.0469, positive slides 9/12, condition+ 2/3, block p 0.00249, LOSO sign changes 0.
- TLS_B_z: **survives_moderate**; median slide beta 0.0394, positive slides 8/12, condition+ 2/3, block p 0.0698, LOSO sign changes 0.
### gene_bin_infercnv_lite
- QC/composition predicts territory within slide: median apparent AUROC **0.985**.
- QC/composition leave-one-slide-out AUROC: **0.621**.
- TROP2_raw: **directionally_inconsistent**; median slide beta 0.135, positive slides 7/12, condition+ 2/3, block p 0.00249, LOSO sign changes 0.
- TACSTD2_z: **directionally_inconsistent**; median slide beta 0.0369, positive slides 7/12, condition+ 1/3, block p 0.00249, LOSO sign changes 0.
- Macrophage_TAM_z: **survives_strong**; median slide beta 0.0916, positive slides 10/12, condition+ 3/3, block p 0.00249, LOSO sign changes 0.
- Tcell_cytotoxic_z: **survives_moderate**; median slide beta 0.0643, positive slides 9/12, condition+ 2/3, block p 0.00249, LOSO sign changes 0.
- TLS_B_z: **survives_moderate**; median slide beta 0.08, positive slides 10/12, condition+ 2/3, block p 0.0299, LOSO sign changes 0.

Interpretation: treat TROP2/TACSTD2 as QC-sensitive exploratory phenotypes unless they survive slide-level and block-level tests. The most defensible biology is whatever survives as positive across slides, conditions, block permutation, and leave-one-slide influence.

## Scorecard

| Label set | Outcome | Verdict | median beta | +slides | +conditions | block beta | block p | LOSO sign changes |
|---|---|---|---:|---:|---:|---:|---:|---:|
| arm_level_expression_cnv | Macrophage_TAM_z | survives_strong | 0.107 | 10/12 | 3/3 | 0.423 | 0.00249 | 0 |
| arm_level_expression_cnv | TACSTD2_z | directionally_inconsistent | 0.00982 | 7/12 | 1/3 | -1.05 | 0.00249 | 0 |
| arm_level_expression_cnv | TLS_B_z | survives_moderate | 0.0394 | 8/12 | 2/3 | 0.199 | 0.0698 | 0 |
| arm_level_expression_cnv | TROP2_raw | fails_or_reverses | -0.00349 | 6/12 | 2/3 | -1.55 | 0.00249 | 0 |
| arm_level_expression_cnv | Tcell_cytotoxic_z | survives_moderate | 0.0469 | 9/12 | 2/3 | 0.803 | 0.00249 | 0 |
| gene_bin_infercnv_lite | Macrophage_TAM_z | survives_strong | 0.0916 | 10/12 | 3/3 | 0.465 | 0.00249 | 0 |
| gene_bin_infercnv_lite | TACSTD2_z | directionally_inconsistent | 0.0369 | 7/12 | 1/3 | -1 | 0.00249 | 0 |
| gene_bin_infercnv_lite | TLS_B_z | survives_moderate | 0.08 | 10/12 | 2/3 | 0.217 | 0.0299 | 0 |
| gene_bin_infercnv_lite | TROP2_raw | directionally_inconsistent | 0.135 | 7/12 | 2/3 | -1.5 | 0.00249 | 0 |
| gene_bin_infercnv_lite | Tcell_cytotoxic_z | survives_moderate | 0.0643 | 9/12 | 2/3 | 0.837 | 0.00249 | 0 |

## Raw Slide High-Low Summary

| Label set | Outcome | Scope | median delta | +slides | sign p | Wilcoxon p |
|---|---|---|---:|---:|---:|---:|
| arm_level_expression_cnv | Macrophage_TAM_z | all_cancer_slides | 0.346 | 10/12 | 0.0386 | 0.0269 |
| arm_level_expression_cnv | TACSTD2_z | all_cancer_slides | 1.05 | 10/12 | 0.0386 | 0.0269 |
| arm_level_expression_cnv | TLS_B_z | all_cancer_slides | 0.176 | 9/12 | 0.146 | 0.301 |
| arm_level_expression_cnv | TROP2_raw | all_cancer_slides | 4.44 | 10/12 | 0.0386 | 0.0122 |
| arm_level_expression_cnv | Tcell_cytotoxic_z | all_cancer_slides | 0.084 | 10/12 | 0.0386 | 0.0425 |
| gene_bin_infercnv_lite | Macrophage_TAM_z | all_cancer_slides | 0.392 | 10/12 | 0.0386 | 0.0269 |
| gene_bin_infercnv_lite | TACSTD2_z | all_cancer_slides | 1.13 | 10/12 | 0.0386 | 0.021 |
| gene_bin_infercnv_lite | TLS_B_z | all_cancer_slides | 0.153 | 9/12 | 0.146 | 0.301 |
| gene_bin_infercnv_lite | TROP2_raw | all_cancer_slides | 4.51 | 10/12 | 0.0386 | 0.00928 |
| gene_bin_infercnv_lite | Tcell_cytotoxic_z | all_cancer_slides | 0.0859 | 10/12 | 0.0386 | 0.0269 |

## QC/Composition Territory Predictability

| Label set | Mode | n | AUROC | high rate |
|---|---|---:|---:|---:|
| arm_level_expression_cnv | leave_one_slide_out | 43180 | 0.601 | 0.332 |
| gene_bin_infercnv_lite | leave_one_slide_out | 43180 | 0.621 | 0.318 |
| arm_level_expression_cnv | within_slide_apparent | 43180 | 0.97 | 0.329 |
| gene_bin_infercnv_lite | within_slide_apparent | 43180 | 0.985 | 0.298 |

## Output files
- `/home/seungho/personal/THCA_data_analysis/project/results/high_impact_topic_pilots_2026_05_08/spatial_cnv_reviewer_stress_battery_2026_05_08/condition_pooled_residual_effects.tsv`
- `/home/seungho/personal/THCA_data_analysis/project/results/high_impact_topic_pilots_2026_05_08/spatial_cnv_reviewer_stress_battery_2026_05_08/loso_influence.tsv`
- `/home/seungho/personal/THCA_data_analysis/project/results/high_impact_topic_pilots_2026_05_08/spatial_cnv_reviewer_stress_battery_2026_05_08/qc_composition_predicts_territory.tsv`
- `/home/seungho/personal/THCA_data_analysis/project/results/high_impact_topic_pilots_2026_05_08/spatial_cnv_reviewer_stress_battery_2026_05_08/reviewer_stress_scorecard.tsv`
- `/home/seungho/personal/THCA_data_analysis/project/results/high_impact_topic_pilots_2026_05_08/spatial_cnv_reviewer_stress_battery_2026_05_08/slide_high_low_effects.tsv`
- `/home/seungho/personal/THCA_data_analysis/project/results/high_impact_topic_pilots_2026_05_08/spatial_cnv_reviewer_stress_battery_2026_05_08/slide_high_low_summary.tsv`
- `/home/seungho/personal/THCA_data_analysis/project/results/high_impact_topic_pilots_2026_05_08/spatial_cnv_reviewer_stress_battery_2026_05_08/slide_residual_effects.tsv`
- `/home/seungho/personal/THCA_data_analysis/project/results/high_impact_topic_pilots_2026_05_08/spatial_cnv_reviewer_stress_battery_2026_05_08/slide_residual_summary.tsv`
- `/home/seungho/personal/THCA_data_analysis/project/results/high_impact_topic_pilots_2026_05_08/spatial_cnv_reviewer_stress_battery_2026_05_08/spatial_block_permutation.tsv`