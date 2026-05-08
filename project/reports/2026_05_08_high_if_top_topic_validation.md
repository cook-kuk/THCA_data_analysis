# High-IF top-topic validation - 2026-05-08

Scope: H1-H7 high-probability directions from the high-IF strategy memo. This run added a spatial expression-CNV bridge using downloaded GENCODE v48 and UCSC hg38 cytoband annotations.

## New data added

| Data | Local file | Use |
|---|---|---|
| GENCODE v48 gene annotation | `/home/seungho/personal/THCA_data_analysis/project/data/external_annotations/gencode.v48.annotation.gtf.gz` | map Visium genes to chromosomes and arms |
| UCSC hg38 cytoband | `/home/seungho/personal/THCA_data_analysis/project/data/external_annotations/hg38_cytoBand.txt.gz` | assign p/q arm labels |
| Existing public thyroid spatial corpus | `project/data/processed/GSE250521/*.h5ad` | infer spot-level expression-CNV and ecosystems |
| Existing UNI morphology embeddings | `phase1_gse250521/uni_embeddings_size224.npz` | evidence-grounded multimodal retrieval pilot |
| Existing TCGA/cBio/GATCI/GeoMx evidence | `high_impact_topic_pilots_2026_05_08/*` | orthogonal validation layers |

## Topic verdict

| Rank | Topic | Grade | Key validation | Decision |
|---:|---|---:|---|---|
| 1 | H1 CNV-residual spatial ecosystem | A- | spatial expression-CNV stage effect max |d|=0.28; normal-q95 high rate max=23.31%; median Moran I=0.107 | Flagship; now has a spatial-CNV bridge, still needs orthogonal wet/IF validation for A/A+. |
| 2 | H2 spatial phylogeography of driver-negative THCA | B+ | spatial-CNV territories are measurable; signature Moran I median=0.107 | Promising but still expression-CNV proxy, not allele-specific CalicoST-grade clone phylogeny. |
| 3 | H3 multimodal visibility boundary | B+ | UNI retrieval best signature uni_only rho=0.044/AUC=0.511; TROP2 text_only rho=0.494/AUC=0.800 | Companion method boundary; use as evidence-grounded retrieval, not as main claim. |
| 4 | H4 thyroid spatial ecotypes | B+ | k-means ecotypes recovered n=6 ecosystem states with labels: CAF+TAM, TAM+CAF, TROP2+TAM, TAM+TROP2, TLS+TAM, TAM+TLS | Good integrated figure, but not enough as standalone without external ecotype cohort. |
| 5 | H5 CAF/ECM lineage suppression | B+ | spatial-CNV vs CAF_ECM pooled rho=0.084; UNI retrieval uni_text_a08 rho=0.348/AUC=0.748 | Mechanism layer inside H1; keep stromal abundance caveat. |
| 6 | H6 TROP2 spatial niche | A- | spatial-CNV vs TACSTD2 pooled rho=-0.010; UNI retrieval text_only rho=0.494/AUC=0.800 | Strong sub-aim; do not lead as ADC story. |
| 7 | H7 immune/TLS protective thyroid ecotype | B+ | spatial-CNV vs TLS_B rho=0.124; retrieval uni_text_a06 rho=0.245/AUC=0.565; TCGA Cox retained separately. | Useful ecosystem axis; still needs Hashimoto-vs-antitumor separation. |

## H1/H2 spatial expression-CNV bridge

| Comparison | mean condition | mean normal | Cohen d | q95 high rate | OR | p |
|---|---:|---:|---:|---:|---:|---:|
| PTC_vs_N | -0.00349 | 4.62e-17 | -0.107 | 0.06 | 1.21 | 0.000122 |
| LPTC_vs_N | -0.00581 | 4.62e-17 | -0.166 | 0.0836 | 1.73 | 1.69e-31 |
| ATC_vs_N | 0.0117 | 4.62e-17 | 0.285 | 0.233 | 5.77 | 0 |

Top arm-level spatial-expression shifts:

| Arm | Condition | Direction | Delta vs normal | d | p |
|---|---|---|---:|---:|---:|
| 12q | ATC | gain_like | 0.0951 | 0.853 | 0 |
| 7p | ATC | gain_like | 0.092 | 0.848 | 0 |
| 7p | PTC | gain_like | 0.0516 | 0.578 | 0 |
| 7p | LPTC | gain_like | 0.0689 | 0.562 | 0 |
| 2p | ATC | loss_like_negative_expected | 0.0492 | 0.521 | 1.67e-304 |
| 16p | PTC | gain_like | 0.0414 | 0.442 | 0 |
| 2q | PTC | loss_like_negative_expected | 0.0307 | 0.417 | 0 |
| 16q | ATC | gain_like | 0.0463 | 0.415 | 6.81e-263 |
| 2q | ATC | loss_like_negative_expected | 0.0298 | 0.385 | 1.64e-259 |
| 2p | PTC | loss_like_negative_expected | 0.0267 | 0.342 | 1.43e-280 |
| 2q | LPTC | loss_like_negative_expected | 0.0312 | 0.321 | 1.06e-178 |
| 12q | PTC | gain_like | 0.0228 | 0.266 | 1.64e-217 |

Spatial organization:

| Condition | n slides | median Moran I | mean Moran I |
|---|---:|---:|---:|
| ATC | 4 | 0.087 | 0.143 |
| LPTC | 4 | 0.161 | 0.174 |
| N | 4 | 0.0776 | 0.085 |
| PTC | 4 | 0.0849 | 0.101 |

## Ecosystem coupling

| Scope | Target | Spearman r | p/n |
|---|---|---:|---:|
| pooled_cancer_spots | TACSTD2_z | -0.00957 | p=0.0468 |
| slide_median_r | TACSTD2_z | 0.0212 | n_slides=12 |
| pooled_cancer_spots | CAF_ECM_z | 0.0836 | p=8.15e-68 |
| slide_median_r | CAF_ECM_z | 0.126 | n_slides=12 |
| pooled_cancer_spots | TLS_B_z | 0.124 | p=1.93e-147 |
| slide_median_r | TLS_B_z | 0.114 | n_slides=12 |
| pooled_cancer_spots | Tcell_cytotoxic_z | 0.0921 | p=5.24e-82 |
| slide_median_r | Tcell_cytotoxic_z | 0.0585 | n_slides=12 |
| pooled_cancer_spots | Macrophage_TAM_z | 0.179 | p=1.26e-306 |
| slide_median_r | Macrophage_TAM_z | 0.146 | n_slides=12 |
| pooled_cancer_spots | RAI_thyroid_z | -0.0241 | p=5.77e-07 |
| slide_median_r | RAI_thyroid_z | 0.0272 | n_slides=12 |
| pooled_cancer_spots | DM1_like_score | 0.000916 | p=0.964 |
| slide_median_r | DM1_like_score | -0.00755 | n_slides=12 |
| pooled_cancer_spots | RAI_8_score | -0.000916 | p=0.964 |
| slide_median_r | RAI_8_score | 0.00755 | n_slides=12 |

## Spatial ecotypes

| Ecotype | Label | n spots | n slides | CNV | TROP2 | CAF | TLS | Tcell | TAM | RAI | high-q95 rate |
|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0 | CAF+TAM | 6575 | 11 | 0.00445 | 0.548 | 2.43 | 0.657 | 0.163 | 1.14 | -1.72 | 0.136 |
| 1 | TAM+CAF | 8387 | 12 | -0.0188 | -0.334 | 0.207 | 0.184 | 0.0173 | 0.258 | -1.93 | 0.0281 |
| 2 | TROP2+TAM | 17954 | 11 | -0.00519 | 3.75 | 0.0451 | 0.137 | -0.0445 | 0.343 | -1.01 | 0.0477 |
| 3 | TAM+TROP2 | 6962 | 12 | 0.0317 | 0.397 | 0.39 | 0.306 | 0.514 | 1.9 | -2.74 | 0.38 |
| 4 | TLS+TAM | 1552 | 12 | -0.00569 | 0.533 | 1.25 | 3.95 | 0.545 | 1.25 | -1.79 | 0.0928 |
| 5 | TAM+TLS | 1750 | 11 | 0.000373 | -0.239 | 0.086 | 0.882 | 6 | 1.84 | -2.99 | 0.102 |

## H3 UNI morphology + metadata retrieval

| Target | Best mode | Spearman r | AUROC q75 | n |
|---|---|---:|---:|---:|
| CAF_ECM_z | uni_text_a08 | 0.348 | 0.748 | 3200 |
| DM1_like_score | uni_only | 0.164 | 0.637 | 3200 |
| RAI_8_score | uni_only | 0.164 | 0.543 | 3200 |
| RAI_thyroid_z | uni_text_a06 | 0.802 | 0.922 | 3200 |
| TACSTD2_z | text_only | 0.494 | 0.8 | 3200 |
| TLS_B_z | uni_text_a06 | 0.245 | 0.565 | 3200 |
| t00_spatial_cnv_signature | uni_only | 0.044 | 0.511 | 3200 |

## External evidence snapshot

| Layer | Metric | Effect | p |
|---|---|---|---:|
| TCGA_cBio_ARMDRIVER | cbio_ARMDRIVER_signature_any | DM2=0.354, DM1=0.0132, OR=41.1 | 1.61e-07 |
| TCGA_cBio_ARMDRIVER | cbio_ARMDRIVER_12q | DM2=0.25, DM1=0, OR=inf | 4.37e-06 |
| TCGA_cBio_ARMDRIVER | cbio_ARMDRIVER_7q | DM2=0.292, DM1=0.0132, OR=30.9 | 4.73e-06 |
| TCGA_cBio_ARMDRIVER | cbio_ARMDRIVER_16p | DM2=0.229, DM1=0, OR=inf | 1.34e-05 |
| TCGA_cBio_ARMDRIVER | ARM_SCNA_CLUSTER=Many SCNA | DM2=0.271, DM1=0.0132, OR=27.9 | 1.39e-05 |
| TCGA_cBio_ARMDRIVER | cbio_ARMDRIVER_7p | DM2=0.25, DM1=0.0132, OR=25 | 4.00e-05 |
| TCGA_cBio_ARMDRIVER | cbio_ARMDRIVER_2p | DM2=0.208, DM1=0, OR=inf | 4.01e-05 |
| TCGA_cBio_ARMDRIVER | cbio_ARMDRIVER_16q | DM2=0.188, DM1=0, OR=inf | 0.000118 |
| TCGA_stratified_robustness | full_Class6_DM2_vs_DM1 | DM2=0.354, DM1=0.0132, OR=41.12903225806452 | 1.61e-07 |
| TCGA_stratified_robustness | follicular_pct<median(70.000) | DM2=0.571, DM1=0, OR=inf | 0.000165 |
| TCGA_stratified_robustness | stage=Stage I | DM2=0.348, DM1=0.02, OR=26.133333333333333 | 0.000261 |
| TCGA_stratified_robustness | purity>=median(0.685) | DM2=0.4, DM1=0, OR=inf | 0.000613 |
| TCGA_stratified_robustness | histology_subtype=cPTC | DM2=0.278, DM1=0.0182, OR=20.76923076923077 | 0.00288 |
| TCGA_stratified_robustness | BRAFV600E_RAS=Ras-like | DM2=0.39, DM1=0, OR=inf | 0.00571 |
| TCGA_stratified_robustness | purity<median(0.685) | DM2=0.364, DM1=0.025, OR=22.285714285714285 | 0.00582 |
| TCGA_stratified_robustness | histology_subtype=FVPTC | DM2=0.407, DM1=0, OR=inf | 0.00714 |
| TCGA_covariate_AUC | mean_cross_validated_delta | covariates=0.935, covariates+CNV=0.947, delta_AUC=0.0116 | NA |
| GSE301163_GeoMx | TROP2 PanCK+ vs VIM+ (microPTC) | mean_a=2.13, mean_b=2.37 | 0.714 |
| GSE301163_GeoMx | TROP2 PanCK+ vs VIM+ (PDTC) | mean_a=1.9, mean_b=1.71 | 0.251 |
| GSE301163_GeoMx | TROP2 PanCK+ vs VIM+ (Normal thyroid) | mean_a=2.13, mean_b=2.12 | 0.748 |
| GSE301163_GeoMx | TROP2 PanCK+ vs VIM+ (Non-neoplastic nodular area) | mean_a=2.59, mean_b=2.06 | 0.0511 |
| GSE301163_GeoMx | TROP2 microPTC PanCK+ vs Normal PanCK+ | mean_a=2.13, mean_b=2.13 | 0.714 |
| GSE301163_GeoMx | DM1_axis microPTC PanCK+ vs Normal PanCK+ | mean_a=4.39, mean_b=4.35 | 1 |
| GSE301163_GeoMx | RAI_8 microPTC PanCK+ vs Normal PanCK+ | mean_a=4.83, mean_b=4.69 | 0.247 |
| GSE301163_GeoMx | Eight_gene_DM1 microPTC PanCK+ vs Normal PanCK+ | mean_a=2.23, mean_b=2.55 | 0.0173 |
| GSE301163_GeoMx | HLA_II microPTC PanCK+ vs Normal PanCK+ | mean_a=2.84, mean_b=2.5 | 0.0303 |
| GSE301163_GeoMx | Cell_cycle microPTC PanCK+ vs Normal PanCK+ | mean_a=2.2, mean_b=2.27 | 0.429 |
| TCGA_DSS_Cox_univariate | CAF | HR=1.37 | 0.268 |

## Bottom line

H1 remains the main high-IF bet. The new run adds the missing bridge: a chromosome-arm spatial expression-CNV layer that can be overlaid with TROP2, CAF/ECM, TLS/immune, RAI, and morphology retrieval. H2 is now a real next experiment, but it remains a proxy until allele-specific CalicoST-grade clone inference is run. H3-H7 should be kept as integrated layers, not separate headline papers.

## Output files

- `/home/seungho/personal/THCA_data_analysis/project/results/high_impact_topic_pilots_2026_05_08/top_topic_validation_2026_05_08/external_evidence_summary.tsv`
- `/home/seungho/personal/THCA_data_analysis/project/results/high_impact_topic_pilots_2026_05_08/top_topic_validation_2026_05_08/gencode_v48_gene_arm_map.tsv`
- `/home/seungho/personal/THCA_data_analysis/project/results/high_impact_topic_pilots_2026_05_08/top_topic_validation_2026_05_08/spatial_cnv_arm_condition_tests.tsv`
- `/home/seungho/personal/THCA_data_analysis/project/results/high_impact_topic_pilots_2026_05_08/top_topic_validation_2026_05_08/spatial_cnv_arm_gene_counts.tsv`
- `/home/seungho/personal/THCA_data_analysis/project/results/high_impact_topic_pilots_2026_05_08/top_topic_validation_2026_05_08/spatial_cnv_axis_correlations.tsv`
- `/home/seungho/personal/THCA_data_analysis/project/results/high_impact_topic_pilots_2026_05_08/top_topic_validation_2026_05_08/spatial_cnv_condition_summary.tsv`
- `/home/seungho/personal/THCA_data_analysis/project/results/high_impact_topic_pilots_2026_05_08/top_topic_validation_2026_05_08/spatial_cnv_condition_tests.tsv`
- `/home/seungho/personal/THCA_data_analysis/project/results/high_impact_topic_pilots_2026_05_08/top_topic_validation_2026_05_08/spatial_cnv_gene_arm_map_used.tsv`
- `/home/seungho/personal/THCA_data_analysis/project/results/high_impact_topic_pilots_2026_05_08/top_topic_validation_2026_05_08/spatial_cnv_morans.tsv`
- `/home/seungho/personal/THCA_data_analysis/project/results/high_impact_topic_pilots_2026_05_08/top_topic_validation_2026_05_08/spatial_ecotype_per_spot.tsv.gz`
- `/home/seungho/personal/THCA_data_analysis/project/results/high_impact_topic_pilots_2026_05_08/top_topic_validation_2026_05_08/spatial_ecotype_silhouette.tsv`
- `/home/seungho/personal/THCA_data_analysis/project/results/high_impact_topic_pilots_2026_05_08/top_topic_validation_2026_05_08/spatial_ecotype_summary.tsv`
- `/home/seungho/personal/THCA_data_analysis/project/results/high_impact_topic_pilots_2026_05_08/top_topic_validation_2026_05_08/spatial_expression_cnv_per_spot.tsv.gz`
- `/home/seungho/personal/THCA_data_analysis/project/results/high_impact_topic_pilots_2026_05_08/top_topic_validation_2026_05_08/spatial_geneset_gene_counts.tsv`
- `/home/seungho/personal/THCA_data_analysis/project/results/high_impact_topic_pilots_2026_05_08/top_topic_validation_2026_05_08/top_topic_validation_scorecard.tsv`
- `/home/seungho/personal/THCA_data_analysis/project/results/high_impact_topic_pilots_2026_05_08/top_topic_validation_2026_05_08/uni_retrieval_metrics.tsv`
- `/home/seungho/personal/THCA_data_analysis/project/results/high_impact_topic_pilots_2026_05_08/top_topic_validation_2026_05_08/uni_retrieval_predictions.tsv.gz`