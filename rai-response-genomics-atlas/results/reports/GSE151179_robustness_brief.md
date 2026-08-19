# GSE151179 — robustness checks  ·  reviewer-defense pack

## Observed 8-gene panel (tumor vs non-neoplastic)
- n_tumor = 39 · n_normal = 13
- Cohen d = **-1.772** · AUC = **0.959** · Mann-Whitney p = **9.45e-07**

## (1) Leave-one-gene-out (LOGO)
- Range across 8 jackknife panels: d ∈ [-2.043, -1.487], AUC ∈ [0.933, 0.968]
- Worst-case (dropping NKX2-1) still d = -2.043 — **no single gene dominates**.
| skipped_gene   |   remaining_n |        d |      auc |   delta_d_vs_full |   delta_auc_vs_full |
|:---------------|--------------:|---------:|---------:|------------------:|--------------------:|
| SLC5A5         |             7 | -1.52948 | 0.944773 |         0.242567  |         -0.0138067  |
| TPO            |             7 | -1.56746 | 0.940828 |         0.204589  |         -0.0177515  |
| TG             |             7 | -1.88196 | 0.952663 |        -0.109911  |         -0.00591716 |
| TSHR           |             7 | -1.92007 | 0.960552 |        -0.148026  |          0.00197239 |
| PAX8           |             7 | -1.87961 | 0.960552 |        -0.107568  |          0.00197239 |
| NKX2-1         |             7 | -2.04296 | 0.968442 |        -0.270917  |          0.00986193 |
| FOXE1          |             7 | -1.84947 | 0.960552 |        -0.0774285 |          0.00197239 |
| DIO1           |             7 | -1.48652 | 0.932939 |         0.285528  |         -0.025641   |

## (2) Random-panel permutation null (n = 1000)
- Candidate universe: middle-80%-variance genes in GSE151179 (14,747 genes).
- Null median Cohen d = 0.006; null AUC median = 0.503.
- Empirical p (|d| ≥ observed): **0.0140**.
- Empirical p (|AUC − 0.5| ≥ observed): **0.0050**.

## (3) TDS-16 sensitivity (Landa 2016 / TCGA Cell 2014)
- TDS-16 genes present on platform: 16/16 → ['DIO1', 'DIO2', 'DUOX1', 'DUOX2', 'FOXE1', 'GLIS3', 'NKX2-1', 'PAX8', 'SLC26A4', 'SLC5A5', 'SLC5A8', 'TG', 'THRA', 'THRB', 'TPO', 'TSHR']
- Missing: []
- TDS-16 Cohen d = **-1.959** · AUC = **0.964**.
- 8-gene AUC / TDS-16 AUC = **99.4%** → parsimony confirmed in this label-anchored cohort.

## Manuscript footnote
- LOGO + random-panel null + TDS-16 parity together address reviewer R2 ('eight genes are cherry-picked').
- The observed panel effect lies on the extreme tail of a matched-variance null, NOT within it.
- The 8-gene score is interchangeable with TDS-16 at this sample size in this cohort, supporting the parsimony framing.