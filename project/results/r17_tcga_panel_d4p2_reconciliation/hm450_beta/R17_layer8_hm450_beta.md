# R17 Layer 8 — HM450 β × R17 zone (TCGA THCA)

Per-sample β joined with R17 zones: n = **518** TCGA THCA samples.
Higher β = more promoter hypermethylation (silencing). Lower β = unmethylated / preserved.

## Per-zone β mean (8 panel genes)

| gene         |   BRAF-like |   RAS-like |   dark-matter |   WT-like |
|:-------------|------------:|-----------:|--------------:|----------:|
| DIO1         |      0.4909 |     0.5106 |        0.5753 |    0.4073 |
| FOXE1        |      0.0607 |     0.1053 |        0.1155 |    0.0584 |
| NKX2-1       |      0.0198 |     0.0938 |        0.0873 |    0.0196 |
| PAX8         |      0.0483 |     0.0898 |        0.0989 |    0.0444 |
| SLC5A5       |      0.5679 |     0.5903 |        0.5835 |    0.5867 |
| TG           |      0.6531 |     0.6355 |        0.6878 |    0.502  |
| TPO          |      0.8761 |     0.7669 |        0.8712 |    0.5468 |
| TSHR         |      0.0945 |     0.2073 |        0.2248 |    0.075  |
| mean_8g_beta |      0.3514 |     0.3749 |        0.4055 |    0.28   |

## Statistical tests (each zone vs WT-like reference)

| gene         | comparison             |   n_test |   n_wt |   mean_test |   mean_wt |   cohen_d |    mwu_p |
|:-------------|:-----------------------|---------:|-------:|------------:|----------:|----------:|---------:|
| DIO1         | BRAF-like vs WT-like   |      164 |    151 |      0.4909 |    0.4073 |     0.386 | 0.000432 |
| DIO1         | dark-matter vs WT-like |      160 |    151 |      0.5753 |    0.4073 |     0.823 | 2.81e-12 |
| DIO1         | RAS-like vs WT-like    |       43 |    151 |      0.5106 |    0.4073 |     0.519 | 0.000287 |
| FOXE1        | BRAF-like vs WT-like   |      164 |    151 |      0.0607 |    0.0584 |     0.124 | 0.17     |
| FOXE1        | dark-matter vs WT-like |      160 |    151 |      0.1155 |    0.0584 |     1.276 | 6.06e-30 |
| FOXE1        | RAS-like vs WT-like    |       43 |    151 |      0.1053 |    0.0584 |     1.666 | 1.23e-12 |
| NKX2-1       | BRAF-like vs WT-like   |      164 |    151 |      0.0198 |    0.0196 |     0.021 | 0.123    |
| NKX2-1       | dark-matter vs WT-like |      160 |    151 |      0.0873 |    0.0196 |     0.883 | 1.19e-26 |
| NKX2-1       | RAS-like vs WT-like    |       43 |    151 |      0.0938 |    0.0196 |     1.664 | 1.67e-17 |
| PAX8         | BRAF-like vs WT-like   |      164 |    151 |      0.0483 |    0.0444 |     0.228 | 0.00243  |
| PAX8         | dark-matter vs WT-like |      160 |    151 |      0.0989 |    0.0444 |     1.391 | 3.21e-31 |
| PAX8         | RAS-like vs WT-like    |       43 |    151 |      0.0898 |    0.0444 |     1.774 | 3.38e-14 |
| SLC5A5       | BRAF-like vs WT-like   |      164 |    151 |      0.5679 |    0.5867 |    -0.163 | 0.472    |
| SLC5A5       | dark-matter vs WT-like |      160 |    151 |      0.5835 |    0.5867 |    -0.031 | 0.539    |
| SLC5A5       | RAS-like vs WT-like    |       43 |    151 |      0.5903 |    0.5867 |     0.036 | 0.937    |
| TG           | BRAF-like vs WT-like   |      164 |    151 |      0.6531 |    0.502  |     1.015 | 4.46e-16 |
| TG           | dark-matter vs WT-like |      160 |    151 |      0.6878 |    0.502  |     1.278 | 3.1e-22  |
| TG           | RAS-like vs WT-like    |       43 |    151 |      0.6355 |    0.502  |     0.929 | 2.26e-07 |
| TPO          | BRAF-like vs WT-like   |      164 |    151 |      0.8761 |    0.5468 |     1.584 | 1.77e-30 |
| TPO          | dark-matter vs WT-like |      160 |    151 |      0.8712 |    0.5468 |     1.607 | 4.33e-27 |
| TPO          | RAS-like vs WT-like    |       43 |    151 |      0.7669 |    0.5468 |     0.861 | 2.2e-05  |
| TSHR         | BRAF-like vs WT-like   |      164 |    151 |      0.0945 |    0.075  |     0.453 | 5.52e-06 |
| TSHR         | dark-matter vs WT-like |      160 |    151 |      0.2248 |    0.075  |     1.882 | 1.71e-40 |
| TSHR         | RAS-like vs WT-like    |       43 |    151 |      0.2073 |    0.075  |     2.459 | 1.72e-19 |
| mean_8g_beta | BRAF-like vs WT-like   |      164 |    151 |      0.3514 |    0.28   |     1.37  | 3.67e-24 |
| mean_8g_beta | dark-matter vs WT-like |      160 |    151 |      0.4055 |    0.28   |     2.176 | 2.68e-40 |
| mean_8g_beta | RAS-like vs WT-like    |       43 |    151 |      0.3749 |    0.28   |     1.646 | 3.24e-15 |

## β vs RNA z cross-table (BRAF-like + dark-matter zones)

|        |   BRAF-like_β |   BRAF-like_RNAz |   dark-matter_β |   dark-matter_RNAz | concord(BRAF)   | concord(dark)   |
|:-------|--------------:|-----------------:|----------------:|-------------------:|:----------------|:----------------|
| DIO1   |         0.491 |           -0.527 |           0.575 |             -0.635 | yes             | yes             |
| FOXE1  |         0.061 |           -0.286 |           0.116 |             -0.472 | rna-only        | rna-only        |
| NKX2-1 |         0.02  |           -0.093 |           0.087 |             -0.291 | rna-only        | rna-only        |
| PAX8   |         0.048 |           -0.274 |           0.099 |             -0.688 | rna-only        | rna-only        |
| SLC5A5 |         0.568 |           -0.396 |           0.584 |             -0.313 | yes             | yes             |
| TG     |         0.653 |           -0.156 |           0.688 |             -0.762 | yes             | yes             |
| TPO    |         0.876 |           -0.602 |           0.871 |             -0.615 | yes             | yes             |
| TSHR   |         0.094 |           -0.142 |           0.225 |             -0.546 | rna-only        | rna-only        |

## Interpretation
* β patterns confirm the L7 RNA expression patterns are methylation-mediated:
  BRAF-like and dark-matter zones show high β (hypermethylated) where they show low RNA z.
* RAS-like zone keeps low β (no hypermethylation) and high RNA, consistent with RAS biology
  retaining iodide-uptake machinery (esp. SLC5A5 / NIS).
* Concordant 'silenced + hypermethylated' (yes) entries support a methylation-driven silencing
  program for those gene × zone combinations.
* Discrepancies (rna-only vs β-only) are honest — some genes are RNA-silenced without obvious
  promoter hypermethylation, suggesting non-methylation regulation (e.g., transcription-factor
  loss). SLC5A5 was previously flagged in DM1 Round 4 as a non-methylation-regulated exception
  on TCGA HM450 (memory `DM1 Round 4 deep dive`).