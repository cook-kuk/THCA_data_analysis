# Fast gene-order inferCNV-lite gate - 2026-05-08

This is a fast gene-order/bin-level expression-CNV analysis using normal thyroid Visium spots as reference. It is a practical bridge between arm-level expression-CNV territories and formal allele-specific CNV callers.

## Verdict

- Cancer-slide median KNN coherence z for inferCNV-lite T00 territories: **57.9**.
- Cancer-slide median bin-bootstrap ARI: **0.859**.
- Pooled cancer concordance with prior arm-level territory: ARI **0.706**, signature rho **1**, high-territory AUROC **0.648**.
- Pooled all-slide concordance with prior arm-level territory: ARI **0.68**, signature rho **1**.
- Interpretation: gene-order inferCNV-lite supports the spatial expression-CNV territory story strongly enough for a main/supplementary figure, while still requiring cautious wording.

## Data
- Spots analyzed: **57997**.
- Slides analyzed: **16**.
- Cancer spots: **43180**.

## Ecosystem Enrichment: InferCNV-lite T00-high vs T00-low

| Scope | Axis | delta high-low | d | p |
|---|---|---:|---:|---:|
| pooled_cancer_spots | DM1_like_score | -0.183 | -0.439 | 1.84e-23 |
| slide_median_effect | DM1_like_score | -0.288 | -0.836 | NA |
| pooled_cancer_spots | Macrophage_TAM_z | 0.38 | 0.339 | 7.63e-290 |
| slide_median_effect | Macrophage_TAM_z | 0.392 | 0.386 | NA |
| pooled_cancer_spots | RAI_thyroid_z | 0.431 | 0.501 | 0 |
| slide_median_effect | RAI_thyroid_z | 0.314 | 0.605 | NA |
| pooled_cancer_spots | TACSTD2_z | 0.818 | 0.357 | 0 |
| slide_median_effect | TACSTD2_z | 1.13 | 0.536 | NA |
| pooled_cancer_spots | TLS_B_z | 0.235 | 0.252 | 1.30e-274 |
| slide_median_effect | TLS_B_z | 0.153 | 0.364 | NA |
| pooled_cancer_spots | TROP2_raw | 3.5 | 0.709 | 0 |
| slide_median_effect | TROP2_raw | 4.51 | 1.09 | NA |

## Pathology Bridge

| Mode | Target | Metric | Value | n |
|---|---|---|---:|---:|
| UNI_image_only | infercnv_t00_high | AUROC | 0.506 | 3200 |
| cell_morph_only | infercnv_t00_high | AUROC | 0.546 | 3200 |
| stage_text_only | infercnv_t00_high | AUROC | 0.492 | 3200 |
| UNI_plus_cell_morph | infercnv_t00_high | AUROC | 0.507 | 3200 |
| UNI_plus_stage_text | infercnv_t00_high | AUROC | 0.499 | 3200 |
| UNI_cell_stage | infercnv_t00_high | AUROC | 0.501 | 3200 |
| UNI_image_only | infercnv_aneuploidy_high | AUROC | 0.465 | 3200 |
| cell_morph_only | infercnv_aneuploidy_high | AUROC | 0.475 | 3200 |
| stage_text_only | infercnv_aneuploidy_high | AUROC | 0.483 | 3200 |
| UNI_plus_cell_morph | infercnv_aneuploidy_high | AUROC | 0.466 | 3200 |
| UNI_plus_stage_text | infercnv_aneuploidy_high | AUROC | 0.466 | 3200 |
| UNI_cell_stage | infercnv_aneuploidy_high | AUROC | 0.467 | 3200 |
| UNI_image_only | infercnv_t00_signature | spearman | 0.0862 | 3200 |
| cell_morph_only | infercnv_t00_signature | spearman | -0.188 | 3200 |
| stage_text_only | infercnv_t00_signature | spearman | -0.279 | 3200 |
| UNI_plus_cell_morph | infercnv_t00_signature | spearman | 0.075 | 3200 |
| UNI_plus_stage_text | infercnv_t00_signature | spearman | -0.0132 | 3200 |
| UNI_cell_stage | infercnv_t00_signature | spearman | -0.0221 | 3200 |
| UNI_image_only | infercnv_aneuploidy_score | spearman | 0.465 | 3200 |
| cell_morph_only | infercnv_aneuploidy_score | spearman | -0.229 | 3200 |
| stage_text_only | infercnv_aneuploidy_score | spearman | 0.393 | 3200 |
| UNI_plus_cell_morph | infercnv_aneuploidy_score | spearman | 0.447 | 3200 |
| UNI_plus_stage_text | infercnv_aneuploidy_score | spearman | 0.435 | 3200 |
| UNI_cell_stage | infercnv_aneuploidy_score | spearman | 0.419 | 3200 |

## Manuscript Use
Use this as a sensitivity layer: gene-order expression-CNV territories are spatially coherent and concordant with the previous arm-level territory, but they are not allele-specific clone calls. The safest high-impact framing remains TCGA/cBio genomic CNV residual class as the anchor, spatial transcriptomics as the territory/ecosystem map, and pathology/UNI as the downstream phenotype reader.

## Output files
- `/home/seungho/personal/THCA_data_analysis/project/results/high_impact_topic_pilots_2026_05_08/infercnv_lite_fast_2026_05_08/infercnv_lite_fast_bin_names.tsv`
- `/home/seungho/personal/THCA_data_analysis/project/results/high_impact_topic_pilots_2026_05_08/infercnv_lite_fast_2026_05_08/infercnv_lite_fast_clone_summary.tsv`
- `/home/seungho/personal/THCA_data_analysis/project/results/high_impact_topic_pilots_2026_05_08/infercnv_lite_fast_2026_05_08/infercnv_lite_fast_coherence.tsv`
- `/home/seungho/personal/THCA_data_analysis/project/results/high_impact_topic_pilots_2026_05_08/infercnv_lite_fast_2026_05_08/infercnv_lite_fast_concordance_to_arm_territory.tsv`
- `/home/seungho/personal/THCA_data_analysis/project/results/high_impact_topic_pilots_2026_05_08/infercnv_lite_fast_2026_05_08/infercnv_lite_fast_ecosystem_enrichment.tsv`
- `/home/seungho/personal/THCA_data_analysis/project/results/high_impact_topic_pilots_2026_05_08/infercnv_lite_fast_2026_05_08/infercnv_lite_fast_per_spot.tsv.gz`
- `/home/seungho/personal/THCA_data_analysis/project/results/high_impact_topic_pilots_2026_05_08/infercnv_lite_fast_2026_05_08/infercnv_lite_fast_stability.tsv`
- `/home/seungho/personal/THCA_data_analysis/project/results/high_impact_topic_pilots_2026_05_08/infercnv_lite_fast_2026_05_08/pathology_to_infercnv_lite_fast_metrics.tsv`
- `/home/seungho/personal/THCA_data_analysis/project/results/high_impact_topic_pilots_2026_05_08/infercnv_lite_fast_2026_05_08/pathology_to_infercnv_lite_fast_predictions.tsv.gz`