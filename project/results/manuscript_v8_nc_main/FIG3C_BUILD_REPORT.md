# FIG3C BUILD REPORT

## Status: SUCCESS

## Data sources
- **Methylation (HM450 beta):** `/home/seungho/personal/THCA_data_analysis/project/results/audit_2026_04_30/round5/r5_2_sample_methylation_8gene.tsv`  (503 rows, 8 gene columns)
- **RNA expression:** `/home/seungho/personal/THCA_data_analysis/project/results/ncomm_push_2026_05_08/cbioportal_sweep/panel_expression_thpa_tcga_gdc.tsv`  (513 rows, log2-normalised)
- **DM labels:** `/home/seungho/personal/THCA_data_analysis/project/results/manuscript_v8_nc_main/master_tcga.tsv`  (500 rows; DM1/DM2 calls)

## Sample overlap
- Expression–methylation matched: **509** samples
  - DM1: 138
  - DM2: 366
  - DM label missing: 5

## Spearman correlations (beta vs RNA expression)
| Gene | r | p | n | Direction |
|------|-----|-----|---|-----------|
| TPO | -0.680 | 3.11e-70 | 509 | NEGATIVE (expected) |
| DIO1 | -0.332 | 1.56e-14 | 509 | NEGATIVE (expected) |
| TSHR | -0.359 | 6.11e-17 | 509 | NEGATIVE (expected) |
| TG | -0.464 | 1.66e-28 | 509 | NEGATIVE (expected) |

## Notes
- All 4 genes show negative Spearman correlation (higher beta → lower expression), consistent with epigenetic silencing.
- The panel replaces synthetic `(1 - beta) + rng.normal(0, 0.08, n)` noise used in fig3_epigenetic.py panel C.
- Sample ID matching: expression `sampleId` (e.g. TCGA-BJ-A0YZ-01) → first 12 chars = patient ID matching methylation `sample_short`.

## Output files
- `/home/seungho/personal/THCA_data_analysis/project/results/manuscript_v8_nc_main/Fig3C_beta_expression_scatter.png`
- `/home/seungho/personal/THCA_data_analysis/project/results/manuscript_v8_nc_main/Fig3C_beta_expression_scatter.pdf`
- `/home/seungho/personal/THCA_data_analysis/project/results/manuscript_v8_nc_main/fig3c_beta_expression_scatter.py` (this script)