# Spatial CNV territory QC comparison - 2026-05-08

Purpose: compare whether the original arm-level expression-CNV territory and the new gene-bin inferCNV-lite territory survive QC/slide confound checks.

## Verdict

### arm_level_expression_cnv
- Slide-median log-total-count delta: **2.01**.
- Slide-median n-gene delta: **2.26e+03**.
- Residualized TROP2_raw: beta **-0.377**, r **-0.0469**, p **1.85e-22**.
- Residualized TACSTD2_z: beta **-0.357**, r **-0.0864**, p **2.75e-72**.
- Residualized Macrophage_TAM_z: beta **0.181**, r **0.076**, p **2.84e-56**.
- Residualized TLS_B_z: beta **0.0406**, r **0.0208**, p **1.49e-05**.
- Residualized Tcell_cytotoxic_z: beta **0.193**, r **0.0865**, p **1.79e-72**.
### gene_bin_infercnv_lite
- Slide-median log-total-count delta: **2.04**.
- Slide-median n-gene delta: **2.31e+03**.
- Residualized TROP2_raw: beta **-0.348**, r **-0.0424**, p **1.17e-18**.
- Residualized TACSTD2_z: beta **-0.366**, r **-0.0867**, p **7.53e-73**.
- Residualized Macrophage_TAM_z: beta **0.207**, r **0.0853**, p **1.36e-70**.
- Residualized TLS_B_z: beta **0.0503**, r **0.0253**, p **1.43e-07**.
- Residualized Tcell_cytotoxic_z: beta **0.21**, r **0.0923**, p **2.50e-82**.

Interpretation: high-CNV expression territories are strongly coupled to library depth/n-gene counts. TROP2/TACSTD2 enrichment is therefore not robust as a direct CNV-territory consequence after QC residualization; TAM/T cell/TLS immune ecosystem signals are more defensible.

## QC High-Low Summary

| Label set | Scope | Axis | delta high-low | d | p |
|---|---|---|---:|---:|---:|
| arm_level_expression_cnv | pooled_cancer_spots | log_total_counts | 1.68 | 1.26 | 0 |
| arm_level_expression_cnv | pooled_cancer_spots | n_genes_by_counts | 1.41e+03 | 0.951 | 0 |
| arm_level_expression_cnv | pooled_cancer_spots | pct_mito | -0.213 | -0.0619 | 6.71e-52 |
| gene_bin_infercnv_lite | pooled_cancer_spots | log_total_counts | 1.71 | 1.28 | 0 |
| gene_bin_infercnv_lite | pooled_cancer_spots | n_genes_by_counts | 1.47e+03 | 0.999 | 0 |
| gene_bin_infercnv_lite | pooled_cancer_spots | pct_mito | -0.171 | -0.0496 | 2.22e-51 |
| arm_level_expression_cnv | slide_median_effect | log_total_counts | 2.01 | 4.09 | NA |
| arm_level_expression_cnv | slide_median_effect | n_genes_by_counts | 2.26e+03 | 3.09 | NA |
| arm_level_expression_cnv | slide_median_effect | pct_mito | -0.489 | -0.586 | NA |
| gene_bin_infercnv_lite | slide_median_effect | log_total_counts | 2.04 | 4.07 | NA |
| gene_bin_infercnv_lite | slide_median_effect | n_genes_by_counts | 2.31e+03 | 3.55 | NA |
| gene_bin_infercnv_lite | slide_median_effect | pct_mito | -0.505 | -0.49 | NA |

## Residualized Ecosystem

| Label set | Axis | n | raw delta | raw d | residual beta | residual r | p |
|---|---|---:|---:|---:|---:|---:|---:|
| arm_level_expression_cnv | TACSTD2_z | 43180 | 0.409 | 0.179 | -0.357 | -0.0864 | 2.75e-72 |
| arm_level_expression_cnv | TROP2_raw | 43180 | 2.1 | 0.412 | -0.377 | -0.0469 | 1.85e-22 |
| arm_level_expression_cnv | Macrophage_TAM_z | 43180 | 0.256 | 0.231 | 0.181 | 0.076 | 2.84e-56 |
| arm_level_expression_cnv | TLS_B_z | 43180 | 0.125 | 0.133 | 0.0406 | 0.0208 | 1.49e-05 |
| arm_level_expression_cnv | Tcell_cytotoxic_z | 43180 | 0.323 | 0.241 | 0.193 | 0.0865 | 1.79e-72 |
| arm_level_expression_cnv | RAI_thyroid_z | 43180 | 0.206 | 0.238 | -0.058 | -0.0551 | 2.28e-30 |
| arm_level_expression_cnv | EMT_stress_z | 43180 | 0.168 | 0.308 | -0.043 | -0.0428 | 5.79e-19 |
| arm_level_expression_cnv | DM1_like_score | 2400 | -0.0735 | -0.17 | 0.0919 | 0.0979 | 1.55e-06 |
| arm_level_expression_cnv | RAI_8_score | 2400 | 0.0735 | 0.17 | -0.0919 | -0.0979 | 1.55e-06 |
| arm_level_expression_cnv | n_fib | 2400 | 7.76 | 0.227 | -2.47 | -0.0366 | 0.0733 |
| gene_bin_infercnv_lite | TACSTD2_z | 43180 | 0.493 | 0.216 | -0.366 | -0.0867 | 7.53e-73 |
| gene_bin_infercnv_lite | TROP2_raw | 43180 | 2.36 | 0.465 | -0.348 | -0.0424 | 1.17e-18 |
| gene_bin_infercnv_lite | Macrophage_TAM_z | 43180 | 0.263 | 0.237 | 0.207 | 0.0853 | 1.36e-70 |
| gene_bin_infercnv_lite | TLS_B_z | 43180 | 0.135 | 0.144 | 0.0503 | 0.0253 | 1.43e-07 |
| gene_bin_infercnv_lite | Tcell_cytotoxic_z | 43180 | 0.393 | 0.294 | 0.21 | 0.0923 | 2.50e-82 |
| gene_bin_infercnv_lite | RAI_thyroid_z | 43180 | 0.221 | 0.255 | -0.0641 | -0.0597 | 2.35e-35 |
| gene_bin_infercnv_lite | EMT_stress_z | 43180 | 0.163 | 0.298 | -0.0374 | -0.0365 | 3.34e-14 |
| gene_bin_infercnv_lite | DM1_like_score | 2400 | -0.0828 | -0.191 | 0.0916 | 0.0951 | 3.05e-06 |
| gene_bin_infercnv_lite | RAI_8_score | 2400 | 0.0828 | 0.191 | -0.0916 | -0.0951 | 3.05e-06 |
| gene_bin_infercnv_lite | n_fib | 2400 | 7.59 | 0.223 | -3.21 | -0.0464 | 0.0229 |

## Output files
- `/home/seungho/personal/THCA_data_analysis/project/results/high_impact_topic_pilots_2026_05_08/spatial_cnv_territory_qc_compare_2026_05_08/territory_ecosystem_residualized_compare.tsv`
- `/home/seungho/personal/THCA_data_analysis/project/results/high_impact_topic_pilots_2026_05_08/spatial_cnv_territory_qc_compare_2026_05_08/territory_qc_highlow_compare.tsv`