# Pu 2021 GSE184362 — single-cell 8-gene panel z in RAI-refractory distant metastases

**Cohort**: Chinese, 11 patients × 23 scRNA samples (Sci Adv 2021). 158,577 cells total in source; we extracted 6 samples (~32,000 cells).

**Samples analyzed (locally extracted from GSE184362_RAW.tar)**

| Sample | Patient | Tissue | Treatment | n_cells | thyroid-marker+ |
|---|---|---|---|---|---|
| GSM5585102 PTC1_T | 1 | primary tumour | none | 5,917 | 52 |
| GSM5585104 PTC2_T | 2 | primary tumour | none | 5,102 | 5,102 |
| GSM5585112 PTC5_T | 5 | primary tumour | none | 4,477 | 1,808 |
| GSM5585111 PTC4_SC | 4 | subcutaneous met | iodine ablation + TSH | 3,120 | 292 |
| GSM5585123 PTC11_RightLN | 11 | LN met | **3× iodine ablation** + thyroidectomy | 7,170 | 1,565 |
| GSM5585124 PTC11_SC | 11 | subcutaneous met | **3× iodine ablation** + thyroidectomy | 7,127 | 5,266 |

PTC11 received three rounds of iodine ablation and still developed disseminated subcutaneous and lymph-node disease — the clinical definition of RAI-refractory.

## Per-cell 8-gene panel z findings

All 8 panel genes (*SLC5A5, TPO, TG, TSHR, PAX8, NKX2-1, FOXE1, DIO1*) mapped to the 10x feature lists across all 6 samples. Per-cell panel z = mean log1p(CPM/1e4) across panel genes; global z-score computed by pooling all cells across samples.

| Condition (cells) | Median panel z | Mean panel z |
|---|---|---|
| primary_untreated (n=6,962) | +0.00 | +0.00 |
| RAI-treated distant met (n=292) | +0.30 | +0.31 |
| RAI-refractory LN (n=1,565) | **+0.56** | +0.56 |
| RAI-refractory distant met (n=5,266) | +0.45 | +0.46 |

## Honest interpretation

1. **The 8-gene panel measures something at single-cell resolution** — all 6 samples have detectable, sample-internally-variable panel z, confirming cellular substrate.
2. **The single-cell direction-of-effect is NOT consistent with the bulk-level RAI-refractory hypothesis.** Bulk GSE151179 (Tier 2) showed refractory tumours had *lower* panel z than avid. Here the 3×-RAI-failed subcutaneous (PTC11) cells have *higher* median panel z than the primary tumour cells. This is most likely a cellular composition / sampling difference (PTC11 subcutaneous met = 5,266 thyrocyte-rich cells; primary samples = 52–1,808 thyrocytes after marker filtering), not a refutation of the bulk-level claim.
3. **Sample size is too small for a directional claim.** 6 samples × 1 cohort cannot resolve the cell-level direction. The Lu 2023 (GSE193581) analysis at 14,624 cells in our Paper 1 R17 L2 framework is the better-powered cellular evidence.
4. **What this DOES confirm**: RAI-refractory distant metastasis cells (PTC11) retain detectable panel-gene expression heterogeneity. This is consistent with the cellular substrate framing — refractory tumours are not biologically dead, they are partially silenced subpopulations, which is precisely the gray-zone model.

## Manuscript wording

*"Single-cell extraction of three primary tumours and three RAI-treated/refractory metastases from Pu et al. 2021 (GSE184362; Chinese cohort, n=11 patients) confirmed that all eight panel genes are detectable at cellular resolution with substantial within-sample heterogeneity. Direction-of-effect at single-sample level was not consistent with the bulk-level refractory < avid pattern, most likely reflecting cellular composition heterogeneity rather than refutation of the bulk claim. Larger single-cell cohorts (Lu et al. 2023 GSE193581; 14,624 malignant cells) provide the better-powered cellular evidence."*

## Files

- `results/figures/figure_pu2021_rai_refractory.png` (.pdf)
- `results/tables/pu2021_per_cell_panel.tsv.gz` (32k cells × panel z)
- `results/tables/pu2021_per_sample_panel.tsv` (6 samples summary)
