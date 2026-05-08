# Pathology to ST-CNV territory bridge - 2026-05-08

Ground truth: spatial transcriptome-derived T00 CNV territory from the territory gate. Prediction: leave-one-slide-out from UNI tile embeddings and image-derived cell morphology features.

## Verdict

- UNI image-only to ST-CNV high territory: AUROC **0.534**.
- UNI image-only to continuous ST-CNV signature: rho **0.0862**.
- Cell morphology-only AUROC: **0.555**.
- Stage/text-only AUROC: **0.462**.
- Best overall mode: **cell_morph_only**, AUROC **0.555**.
- Interpretation: pathology image alone does not recover ST-CNV territory strongly; this is an important negative boundary.

## Data
- Tile embeddings evaluated: **3200** spots.
- Slides: **16**.
- T00-high territory prevalence in evaluated tiles: **35.5%**.

## Main Metrics

| Mode | Target | Metric | Value | n |
|---|---|---|---:|---:|
| UNI_image_only | t00_territory_high | AUROC | 0.534 | 3200 |
| cell_morph_only | t00_territory_high | AUROC | 0.555 | 3200 |
| stage_text_only | t00_territory_high | AUROC | 0.462 | 3200 |
| UNI_plus_cell_morph | t00_territory_high | AUROC | 0.537 | 3200 |
| UNI_image_only | t00_spatial_cnv_signature | spearman | 0.0862 | 3200 |
| cell_morph_only | t00_spatial_cnv_signature | spearman | -0.188 | 3200 |
| stage_text_only | t00_spatial_cnv_signature | spearman | -0.279 | 3200 |
| UNI_plus_cell_morph | t00_spatial_cnv_signature | spearman | 0.075 | 3200 |
| UNI_image_only | TACSTD2_z | spearman | 0.507 | 3200 |
| cell_morph_only | TACSTD2_z | spearman | -0.136 | 3200 |
| stage_text_only | TACSTD2_z | spearman | 0.371 | 3200 |
| UNI_plus_cell_morph | TACSTD2_z | spearman | 0.499 | 3200 |
| UNI_image_only | Macrophage_TAM_z | spearman | 0.475 | 3200 |
| cell_morph_only | Macrophage_TAM_z | spearman | 0.0383 | 3200 |
| stage_text_only | Macrophage_TAM_z | spearman | 0.47 | 3200 |
| UNI_plus_cell_morph | Macrophage_TAM_z | spearman | 0.47 | 3200 |

## Paper Meaning

If image-only stays modest, the paper should not claim H&E can replace ST-CNV mapping. The stronger high-IF framing is: spatial transcriptomics reveals CNV-like territories; pathology can partly read ecosystem programs such as TROP2/TAM/TLS, but direct ST/CNV measurement remains needed for the CNV state itself.

## Output files

- `/home/seungho/personal/THCA_data_analysis/project/results/high_impact_topic_pilots_2026_05_08/pathology_to_st_cnv_territory_fast_2026_05_08/pathology_to_st_cnv_metrics.tsv`
- `/home/seungho/personal/THCA_data_analysis/project/results/high_impact_topic_pilots_2026_05_08/pathology_to_st_cnv_territory_fast_2026_05_08/pathology_to_st_cnv_per_slide_metrics.tsv`
- `/home/seungho/personal/THCA_data_analysis/project/results/high_impact_topic_pilots_2026_05_08/pathology_to_st_cnv_territory_fast_2026_05_08/pathology_to_st_cnv_predictions.tsv.gz`