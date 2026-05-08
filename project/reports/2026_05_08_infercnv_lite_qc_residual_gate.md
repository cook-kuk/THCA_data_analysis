# inferCNV-lite QC/residual gate - 2026-05-08

Purpose: test whether inferCNV-lite T00-high territories are mostly library-size/quality artifacts, and whether ecosystem enrichment survives slide and QC residualization.

## Verdict

- Cancer spots evaluated: **43180**.
- Slide-median high-low log-total-count delta: **2.04**.
- Slide-median high-low n-gene delta: **2.31e+03**.
- Slide-median high-low pct-mito delta: **-0.505**.
- Residualized TROP2_raw association with T00-high after slide + QC: beta **-0.348**, r **-0.0424**, p **1.17e-18**.
- Residualized TACSTD2_z association: beta **-0.366**, r **-0.0867**, p **7.53e-73**.
- Residualized Macrophage_TAM_z association: beta **0.207**, r **0.0853**, p **1.36e-70**.
- Interpretation: if residualized ecosystem betas remain positive, the territory ecosystem claim is not just a slide/QC artifact. QC deltas should still be disclosed because expression-CNV methods are inherently sensitive to transcriptome quality and cell composition.

## QC High-Low Checks

| Scope | Axis | delta high-low | d | p |
|---|---|---:|---:|---:|
| pooled_cancer_spots | log_total_counts | 1.71 | 1.28 | 0 |
| pooled_cancer_spots | n_genes_by_counts | 1.47e+03 | 0.999 | 0 |
| pooled_cancer_spots | pct_mito | -0.171 | -0.0496 | 2.22e-51 |
| slide_median_effect | log_total_counts | 2.04 | 4.07 | NA |
| slide_median_effect | n_genes_by_counts | 2.31e+03 | 3.55 | NA |
| slide_median_effect | pct_mito | -0.505 | -0.49 | NA |

## Residualized Ecosystem Associations

| Axis | n | raw delta | raw d | residual beta | residual r | p |
|---|---:|---:|---:|---:|---:|---:|
| TACSTD2_z | 43180 | 0.493 | 0.216 | -0.366 | -0.0867 | 7.53e-73 |
| TROP2_raw | 43180 | 2.36 | 0.465 | -0.348 | -0.0424 | 1.17e-18 |
| Macrophage_TAM_z | 43180 | 0.263 | 0.237 | 0.207 | 0.0853 | 1.36e-70 |
| TLS_B_z | 43180 | 0.135 | 0.144 | 0.0503 | 0.0253 | 1.43e-07 |
| Tcell_cytotoxic_z | 43180 | 0.393 | 0.294 | 0.21 | 0.0923 | 2.50e-82 |
| RAI_thyroid_z | 43180 | 0.221 | 0.255 | -0.0641 | -0.0597 | 2.35e-35 |
| EMT_stress_z | 43180 | 0.163 | 0.298 | -0.0374 | -0.0365 | 3.34e-14 |
| DM1_like_score | 2400 | -0.0828 | -0.191 | 0.0916 | 0.0951 | 3.05e-06 |
| RAI_8_score | 2400 | 0.0828 | 0.191 | -0.0916 | -0.0951 | 3.05e-06 |
| n_fib | 2400 | 7.59 | 0.223 | -3.21 | -0.0464 | 0.0229 |

## Output files
- `/home/seungho/personal/THCA_data_analysis/project/results/high_impact_topic_pilots_2026_05_08/infercnv_lite_qc_residual_gate_2026_05_08/infercnv_lite_ecosystem_residualized_for_slide_qc.tsv`
- `/home/seungho/personal/THCA_data_analysis/project/results/high_impact_topic_pilots_2026_05_08/infercnv_lite_qc_residual_gate_2026_05_08/infercnv_lite_qc_highlow.tsv`
- `/home/seungho/personal/THCA_data_analysis/project/results/high_impact_topic_pilots_2026_05_08/infercnv_lite_qc_residual_gate_2026_05_08/spot_qc.tsv.gz`