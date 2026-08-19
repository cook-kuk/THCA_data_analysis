# R17 Layer 4 — Mun 2025 proteome (n=336) two-axis projection

Per-sample protein modules, n = **336** Mun 2025 thyroid tumors.

Axes: x = −thyroid_differentiation (protein panel silencing), y = mean(HLA_class_II + IFNG_T_cell_inflamed) at protein level.

## Per-group zone fractions (%)

| group   |   BRAF-like |   RAS-like |   dark-matter |   WT-like |
|:--------|------------:|-----------:|--------------:|----------:|
| ATC     |        36.3 |        3.5 |          58.4 |       1.8 |
| PDTC    |        37   |        8.7 |          23.9 |      30.4 |
| PTC     |        13   |       26.6 |          14.1 |      46.3 |

Fisher (ATC dark-matter vs PTC dark-matter, one-sided): OR = **8.54**, p = **2.72e-15**.

## Interpretation
* PROTEIN layer replicates the R17 two-axis frame: ATC samples concentrate in silenced zones
  (BRAF-like + dark-matter), PTC samples in preserved zones. The dark-matter zone exists
  cross-pillar (RNA + sc + Korean + now protein).
* This is the **4th independent replication layer** of R17, at the proteomic level (n=336).
  The two-axis biology is not RNA-platform specific.

## Caveats
* Protein modules are pre-z-scored; thresholds at 0 inherit from the original Mun analysis.
* PDTC group bridges PTC and ATC and shows mixed-zone occupancy as expected.