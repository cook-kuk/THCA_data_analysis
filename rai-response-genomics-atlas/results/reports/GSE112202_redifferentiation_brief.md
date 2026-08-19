# GSE112202 — Tier 5 redifferentiation direction-of-effect

## Cohort
- 11 digoxin-treated NMTC patients (treated for heart disease before/after thyroid cancer diagnosis) vs 11 matched untreated controls. Cufflinks RNA-seq, group-level FPKM.
- Published claim: digoxin treatment restores thyroid differentiation in vivo, consistent with prior in vitro work.

## 8-gene panel log2FC (digoxin / untreated)
- Genes captured: **8 / 8**.
- Up-regulated by digoxin: **6 / 8**.
- Median log2FC: **+0.302**.
- Direction matches the redifferentiation hypothesis: YES.

## TDS-16 extra genes (sensitivity)
- Up-regulated: **5/8**, median log2FC = +0.154.

## Per-gene table
| gene    |   fpkm_untreated |   fpkm_digoxin |    log2fc | panel_set    |
|:--------|-----------------:|---------------:|----------:|:-------------|
| TSHR    |        19.9912   |      38.7162   |  0.950095 | 8-gene       |
| SLC5A5  |         0.121348 |       0.268345 |  0.734741 | 8-gene       |
| SLC26A4 |         8.26226  |      11.9932   |  0.532231 | TDS-16 extra |
| NKX2-1  |        21.9824   |      31.518    |  0.517849 | 8-gene       |
| GLIS3   |         1.67888  |       2.27698  |  0.418161 | TDS-16 extra |
| TG      |       396.225    |     523.318    |  0.40128  | 8-gene       |
| THRA    |         2.5463   |       3.3523   |  0.383581 | TDS-16 extra |
| DUOX1   |         5.39015  |       6.43276  |  0.250847 | TDS-16 extra |
| TPO     |        37.5336   |      43.1901   |  0.202016 | 8-gene       |
| FOXE1   |        18.0535   |      19.9822   |  0.14567  | 8-gene       |
| THRB    |         0.473423 |       0.496743 |  0.05751  | TDS-16 extra |
| DIO2    |        15.4719   |      15.2192   | -0.023604 | TDS-16 extra |
| PAX8    |        49.5385   |      47.3974   | -0.063611 | 8-gene       |
| DIO1    |        11.5314   |      10.5825   | -0.122775 | 8-gene       |
| DUOX2   |        13.3818   |      12.1012   | -0.14399  | TDS-16 extra |
| SLC5A8  |         2.517    |       2.09315  | -0.254909 | TDS-16 extra |

## Manuscript line
- *"In the redifferentiation cohort GSE112202 (Tier 5; 11 digoxin-treated vs 11 matched untreated NMTC patients), the 8-gene panel direction matched the predicted restoration of thyroid differentiation, supporting mechanistic plausibility."*

## Caveats
- Group-level FPKM only (not per-sample); fold-changes are aggregate ratios. Per-sample paired analysis would require raw count matrix.
- 22-sample retrospective design; not a randomized trial.
- Direction-of-effect alone, not predictive of RAI response per se.