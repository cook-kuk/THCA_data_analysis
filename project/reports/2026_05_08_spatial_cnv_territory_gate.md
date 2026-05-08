# Spatial-CNV territory gate - 2026-05-08

Purpose: decide whether the H1/H2 paper can claim spatially coherent CNV-like territories, not just spot-wise arm-expression noise.

## Gate verdict

- Median same-territory KNN z across cancer slides: **58.7**.
- Median high-vs-low territory CNV-score delta: **0.0224**.
- Median leave-arm-out ARI across cancer slides: **0.87**.
- Interpretation: spatial territories are real enough for a figure and a subclaim, but still expression-CNV proxy, not allele-specific clone phylogeny.

## Territory coherence by slide

| Slide | Condition | high fraction | high-low CNV delta | KNN same obs | null | z | perm p |
|---|---|---:|---:|---:|---:|---:|---:|
| GSM7980860_N-1 | N | 0.153 | 0.0282 | 0.744 | 0.384 | 76.3 | 0.002 |
| GSM7980861_N-2 | N | 0.5 | 0.0238 | 0.676 | 0.384 | 70.5 | 0.002 |
| GSM7980862_N-3 | N | 0.437 | 0.01 | 0.561 | 0.377 | 44 | 0.002 |
| GSM7980863_N-4 | N | 0.478 | 0.0158 | 0.506 | 0.394 | 23.8 | 0.002 |
| GSM7980864_PTC-1 | PTC | 0.329 | 0.0233 | 0.724 | 0.335 | 68.2 | 0.002 |
| GSM7980865_PTC-2 | PTC | 0.47 | 0.0112 | 0.574 | 0.378 | 51.9 | 0.002 |
| GSM7980866_PTC-3 | PTC | 0.49 | 0.00607 | 0.6 | 0.416 | 46.3 | 0.002 |
| GSM7980867_PTC-4 | PTC | 0.251 | 0.0167 | 0.67 | 0.386 | 76 | 0.002 |
| GSM7980868_LPTC-1 | LPTC | 0.403 | 0.0228 | 0.622 | 0.37 | 37.4 | 0.002 |
| GSM7980869_LPTC-2 | LPTC | 0.276 | 0.0351 | 0.581 | 0.342 | 59.4 | 0.002 |
| GSM7980870_LPTC-3 | LPTC | 0.233 | 0.0277 | 0.691 | 0.36 | 83.8 | 0.002 |
| GSM7980871_LPTC-4 | LPTC | 0.213 | 0.0303 | 0.69 | 0.494 | 58 | 0.002 |
| GSM7980872_ATC-1 | ATC | 0.446 | 0.0648 | 0.679 | 0.364 | 80.2 | 0.002 |
| GSM7980873_ATC-2 | ATC | 0.33 | 0.0221 | 0.783 | 0.385 | 78.1 | 0.002 |
| GSM7980874_ATC-3 | ATC | 0.0401 | 0.00634 | 0.607 | 0.489 | 17.4 | 0.002 |
| GSM7980875_ATC-4 | ATC | 0.363 | 0.0202 | 0.61 | 0.4 | 38.6 | 0.002 |

## Ecosystem enrichment in high-CNV territory

| Scope | Axis | delta high-low | d | p |
|---|---|---:|---:|---:|
| pooled_cancer_spots | Macrophage_TAM_z | 0.38 | 0.338 | 9.20e-296 |
| slide_median_effect | Macrophage_TAM_z | 0.346 | 0.35 | NA |
| pooled_cancer_spots | CAF_ECM_z | 0.104 | 0.098 | 1.30e-10 |
| slide_median_effect | CAF_ECM_z | 0.038 | 0.0416 | NA |
| pooled_cancer_spots | TLS_B_z | 0.225 | 0.241 | 8.58e-260 |
| slide_median_effect | TLS_B_z | 0.176 | 0.346 | NA |
| pooled_cancer_spots | Tcell_cytotoxic_z | 0.36 | 0.246 | 3.31e-198 |
| slide_median_effect | Tcell_cytotoxic_z | 0.084 | 0.183 | NA |
| pooled_cancer_spots | TACSTD2_z | 0.703 | 0.306 | 2.23e-243 |
| slide_median_effect | TACSTD2_z | 1.05 | 0.523 | NA |
| pooled_cancer_spots | RAI_thyroid_z | 0.401 | 0.466 | 0 |
| slide_median_effect | RAI_thyroid_z | 0.353 | 0.653 | NA |
| pooled_cancer_spots | DM1_like_score | -0.159 | -0.352 | 3.87e-18 |
| slide_median_effect | DM1_like_score | -0.273 | -0.667 | NA |
| pooled_cancer_spots | RAI_8_score | 0.159 | 0.352 | 3.87e-18 |
| slide_median_effect | RAI_8_score | 0.273 | 0.667 | NA |

## Arm-direction audit

| Arm | True-CNA expectation | delta high-low | matches? | d | p |
|---|---|---:|---:|---:|---:|
| 7p | high>low | 0.147 | True | 1.15 | 0 |
| 7q | high>low | 0.114 | True | 1.38 | 0 |
| 12q | high>low | 0.167 | True | 1.41 | 0 |
| 16p | high>low | 0.145 | True | 1.46 | 0 |
| 16q | high>low | 0.156 | True | 1.32 | 0 |
| 2p | high<low | 0.119 | False | 1.17 | 0 |
| 2q | high<low | 0.111 | False | 1.15 | 0 |

## What this means for the paper

Use the territory analysis as **spatial expression-CNV evidence**, not as definitive copy-number clone calling. The title can safely say CNV-state or CNV-like spatial ecosystem only if the bulk TCGA/cBio ARMDRIVER result is the anchor and the spatial layer is explicitly framed as expression-CNV/territory proxy. The strongest next upgrade is allele-specific CalicoST or direct DNA/IHC/IF validation.

## Output files

- `/home/seungho/personal/THCA_data_analysis/project/results/high_impact_topic_pilots_2026_05_08/spatial_cnv_territory_gate_2026_05_08/spatial_cnv_territory_arm_direction_audit.tsv`
- `/home/seungho/personal/THCA_data_analysis/project/results/high_impact_topic_pilots_2026_05_08/spatial_cnv_territory_gate_2026_05_08/spatial_cnv_territory_coherence.tsv`
- `/home/seungho/personal/THCA_data_analysis/project/results/high_impact_topic_pilots_2026_05_08/spatial_cnv_territory_gate_2026_05_08/spatial_cnv_territory_ecosystem_enrichment.tsv`
- `/home/seungho/personal/THCA_data_analysis/project/results/high_impact_topic_pilots_2026_05_08/spatial_cnv_territory_gate_2026_05_08/spatial_cnv_territory_per_spot.tsv.gz`
- `/home/seungho/personal/THCA_data_analysis/project/results/high_impact_topic_pilots_2026_05_08/spatial_cnv_territory_gate_2026_05_08/spatial_cnv_territory_stability.tsv`