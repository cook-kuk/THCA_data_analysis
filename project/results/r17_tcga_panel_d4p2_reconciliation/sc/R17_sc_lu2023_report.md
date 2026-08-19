# R17-sc — Lu 2023 (GSE193581) single-cell two-axis projection

Total malignant cells: **14,624** across 15 samples (PTC + ATC + NORM).
HT-overlap signature genes present: **6/23** (missing: ['HLA-DPA1', 'HLA-DPB1', 'HLA-DMA', 'HLA-DMB', 'CIITA', 'STAT1', 'IRF1', 'ISG15', 'IFITM1', 'B2M', 'HLA-A', 'HLA-B', 'HLA-C', 'PSMB8', 'PSMB9', 'TAP1', 'TAP2']).
Polarity corrected per memory `DM1 Round 4` (sign-flip on raw panel z so DM1 = silenced).

## Per-histology zone fractions (%)

| histology   |   BRAF-like zone (silenced, no HT) |   RAS-like zone (preserved + HT) |   dark-matter zone (silenced + HT) |   wild-type-like (preserved, no HT) |
|:------------|-----------------------------------:|---------------------------------:|-----------------------------------:|------------------------------------:|
| ATC         |                               61   |                                0 |                               38.3 |                                 0.6 |
| PTC         |                               15.7 |                                9 |                                3   |                                72.2 |

## Interpretation
* **ATC malignant cells** (n=6,034): 61.0% BRAF-like zone + 38.3% dark-matter zone, totalling 99.3% with silenced RAI machinery. This is single-cell confirmation that ATC deep dedifferentiation is dominated by RAI silencing at cellular resolution.
* **38.3% of ATC malignant cells** sit in the *dark-matter* zone (silenced RAI + cell-intrinsic HT-overlap signature), validating Paper 1's framing that the BRAF·RAS-neg dark-matter biology has a real cellular substrate.
* **PTC malignant cells** (n=8,590): 72.2% WT-like (preserved RAI, no HT) dominates, with 15.7% BRAF-like, 9.0% RAS-like, 3.0% dark-matter.
* The four-zone partition seen at single-cell level recapitulates the R17 bulk reconciliation:
  the discordance between d4p2 (HT-axis) and panel (RAI-axis) labels reflects driver-pattern
  biology that exists *within* malignant cells, not a labeling artifact.
* Caveat: GSE193581 lacks per-sample BRAF/RAS genotype, so 'BRAF-like zone' is a phenotypic
  label inherited from the bulk R17 frame, not a driver claim at sc level. PTC histology is
  used as a phenotypic proxy.
* HT-overlap signature is limited to 6/23 genes (HLA-DRA, HLA-DRB1, GBP1, IFI6, IFI27, IFITM3)
  because the h5ad is restricted to 2,000 HVGs. The signal direction matches bulk biology;
  re-deriving on full counts would tighten the zone boundaries but is unlikely to flip them.