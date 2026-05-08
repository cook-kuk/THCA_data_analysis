# Spatial CNV covariate ladder - 2026-05-08

Purpose: quantify how territory-outcome associations change after adding slide, QC, and broad tissue-composition covariates.

## Verdict

### arm_level_expression_cnv
- TROP2_raw: sample_only: beta 1.95, r 0.203; sample_qc: beta -0.377, r -0.0469; sample_qc_composition: beta -0.382, r -0.0477.
- TACSTD2_z: sample_only: beta 0.424, r 0.0989; sample_qc: beta -0.357, r -0.0864; sample_qc_composition: beta -0.36, r -0.0875.
- Macrophage_TAM_z: sample_only: beta 0.27, r 0.121; sample_qc: beta 0.181, r 0.076; sample_qc_composition: beta 0.176, r 0.0742.
- TLS_B_z: sample_only: beta 0.0804, r 0.0429; sample_qc: beta 0.0406, r 0.0208; sample_qc_composition: beta 0.04, r 0.0205.
- Tcell_cytotoxic_z: sample_only: beta 0.338, r 0.16; sample_qc: beta 0.193, r 0.0865; sample_qc_composition: beta 0.19, r 0.0855.
### gene_bin_infercnv_lite
- TROP2_raw: sample_only: beta 2.14, r 0.221; sample_qc: beta -0.348, r -0.0424; sample_qc_composition: beta -0.354, r -0.0434.
- TACSTD2_z: sample_only: beta 0.474, r 0.11; sample_qc: beta -0.366, r -0.0867; sample_qc_composition: beta -0.37, r -0.0879.
- Macrophage_TAM_z: sample_only: beta 0.291, r 0.128; sample_qc: beta 0.207, r 0.0853; sample_qc_composition: beta 0.203, r 0.0839.
- TLS_B_z: sample_only: beta 0.0824, r 0.0435; sample_qc: beta 0.0503, r 0.0253; sample_qc_composition: beta 0.0498, r 0.0251.
- Tcell_cytotoxic_z: sample_only: beta 0.357, r 0.167; sample_qc: beta 0.21, r 0.0923; sample_qc_composition: beta 0.208, r 0.0915.

Interpretation: TROP2/TACSTD2 is a weak claim after QC/composition adjustment. The immune interface, especially TAM and cytotoxic/TLS axes, is the more defensible spatial ecosystem claim.

## Full Ladder

| Label set | Covariates | Outcome | n | beta | r | p |
|---|---|---|---:|---:|---:|---:|
| arm_level_expression_cnv | sample_only | TACSTD2_z | 43180 | 0.424 | 0.0989 | 2.59e-94 |
| arm_level_expression_cnv | sample_only | TROP2_raw | 43180 | 1.95 | 0.203 | 0 |
| arm_level_expression_cnv | sample_only | Macrophage_TAM_z | 43180 | 0.27 | 0.121 | 1.92e-139 |
| arm_level_expression_cnv | sample_only | TLS_B_z | 43180 | 0.0804 | 0.0429 | 4.99e-19 |
| arm_level_expression_cnv | sample_only | Tcell_cytotoxic_z | 43180 | 0.338 | 0.16 | 7.67e-246 |
| arm_level_expression_cnv | sample_only | RAI_thyroid_z | 43180 | 0.154 | 0.138 | 3.95e-182 |
| arm_level_expression_cnv | sample_only | EMT_stress_z | 43180 | 0.179 | 0.162 | 3.41e-251 |
| arm_level_expression_cnv | sample_qc | TACSTD2_z | 43180 | -0.357 | -0.0864 | 2.75e-72 |
| arm_level_expression_cnv | sample_qc | TROP2_raw | 43180 | -0.377 | -0.0469 | 1.85e-22 |
| arm_level_expression_cnv | sample_qc | Macrophage_TAM_z | 43180 | 0.181 | 0.076 | 2.84e-56 |
| arm_level_expression_cnv | sample_qc | TLS_B_z | 43180 | 0.0406 | 0.0208 | 1.49e-05 |
| arm_level_expression_cnv | sample_qc | Tcell_cytotoxic_z | 43180 | 0.193 | 0.0865 | 1.79e-72 |
| arm_level_expression_cnv | sample_qc | RAI_thyroid_z | 43180 | -0.058 | -0.0551 | 2.28e-30 |
| arm_level_expression_cnv | sample_qc | EMT_stress_z | 43180 | -0.043 | -0.0428 | 5.79e-19 |
| arm_level_expression_cnv | sample_qc_composition | TACSTD2_z | 43180 | -0.36 | -0.0875 | 4.50e-74 |
| arm_level_expression_cnv | sample_qc_composition | TROP2_raw | 43180 | -0.382 | -0.0477 | 3.40e-23 |
| arm_level_expression_cnv | sample_qc_composition | Macrophage_TAM_z | 43180 | 0.176 | 0.0742 | 8.04e-54 |
| arm_level_expression_cnv | sample_qc_composition | TLS_B_z | 43180 | 0.04 | 0.0205 | 1.99e-05 |
| arm_level_expression_cnv | sample_qc_composition | Tcell_cytotoxic_z | 43180 | 0.19 | 0.0855 | 8.50e-71 |
| arm_level_expression_cnv | sample_qc_composition | RAI_thyroid_z | 43180 | -0.0563 | -0.0535 | 9.97e-29 |
| arm_level_expression_cnv | sample_qc_composition | EMT_stress_z | 43180 | -0.0446 | -0.0448 | 1.36e-20 |
| gene_bin_infercnv_lite | sample_only | TACSTD2_z | 43180 | 0.474 | 0.11 | 2.55e-115 |
| gene_bin_infercnv_lite | sample_only | TROP2_raw | 43180 | 2.14 | 0.221 | 0 |
| gene_bin_infercnv_lite | sample_only | Macrophage_TAM_z | 43180 | 0.291 | 0.128 | 1.54e-157 |
| gene_bin_infercnv_lite | sample_only | TLS_B_z | 43180 | 0.0824 | 0.0435 | 1.63e-19 |
| gene_bin_infercnv_lite | sample_only | Tcell_cytotoxic_z | 43180 | 0.357 | 0.167 | 4.21e-268 |
| gene_bin_infercnv_lite | sample_only | RAI_thyroid_z | 43180 | 0.163 | 0.144 | 4.04e-200 |
| gene_bin_infercnv_lite | sample_only | EMT_stress_z | 43180 | 0.195 | 0.175 | 2.44e-293 |
| gene_bin_infercnv_lite | sample_qc | TACSTD2_z | 43180 | -0.366 | -0.0867 | 7.53e-73 |
| gene_bin_infercnv_lite | sample_qc | TROP2_raw | 43180 | -0.348 | -0.0424 | 1.17e-18 |
| gene_bin_infercnv_lite | sample_qc | Macrophage_TAM_z | 43180 | 0.207 | 0.0853 | 1.36e-70 |
| gene_bin_infercnv_lite | sample_qc | TLS_B_z | 43180 | 0.0503 | 0.0253 | 1.43e-07 |
| gene_bin_infercnv_lite | sample_qc | Tcell_cytotoxic_z | 43180 | 0.21 | 0.0923 | 2.50e-82 |
| gene_bin_infercnv_lite | sample_qc | RAI_thyroid_z | 43180 | -0.0641 | -0.0597 | 2.35e-35 |
| gene_bin_infercnv_lite | sample_qc | EMT_stress_z | 43180 | -0.0374 | -0.0365 | 3.34e-14 |
| gene_bin_infercnv_lite | sample_qc_composition | TACSTD2_z | 43180 | -0.37 | -0.0879 | 7.05e-75 |
| gene_bin_infercnv_lite | sample_qc_composition | TROP2_raw | 43180 | -0.354 | -0.0434 | 1.93e-19 |
| gene_bin_infercnv_lite | sample_qc_composition | Macrophage_TAM_z | 43180 | 0.203 | 0.0839 | 2.90e-68 |
| gene_bin_infercnv_lite | sample_qc_composition | TLS_B_z | 43180 | 0.0498 | 0.0251 | 1.84e-07 |
| gene_bin_infercnv_lite | sample_qc_composition | Tcell_cytotoxic_z | 43180 | 0.208 | 0.0915 | 5.17e-81 |
| gene_bin_infercnv_lite | sample_qc_composition | RAI_thyroid_z | 43180 | -0.0626 | -0.0584 | 6.82e-34 |
| gene_bin_infercnv_lite | sample_qc_composition | EMT_stress_z | 43180 | -0.0384 | -0.0378 | 4.16e-15 |

## Output files
- `/home/seungho/personal/THCA_data_analysis/project/results/high_impact_topic_pilots_2026_05_08/spatial_cnv_covariate_ladder_2026_05_08/spatial_cnv_covariate_ladder.tsv`