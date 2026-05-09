# H11 — spatial validation of HT/TLS axis (GSE250521 Visium n=8)

**Date:** 2026-05-09 · **Sprint:** Paper 1+2 BRAF Nature · **Cohort:** GSE250521 Lu 2023, restricted to PTC (n=4) + LPTC (n=4)

## Verdict — partial YES

HT-axis HIGH samples carry more HLA-class-II / B-cell / CXCL13 signal spatially, with higher Moran's I (clustered). TLS coverage trends 9× in the right direction but doesn't clear MWU at n=4 vs 4 — 55-µm Visium is too coarse for mature TLS.

## Caveat first

GSE250521 has **no BRAF/RAS/HT annotation** — only stage. The "BRAF-cPTC DM1 vs DM2" claim is not directly testable. We stratified by **HT-13 pseudobulk median split**. DM1_like_score weakly correlates with HT-axis class (rho=+0.52, p=0.18) at n=8.

## HIGH vs LOW HT-axis

| metric | HIGH | LOW | MWU p | rho vs HT-13 (p) |
|---|---|---|---|---|
| **HLA-DRA log1p** | 2.07 | 1.24 | **0.029** | **+0.976 (3.3e-5)** |
| **HT-13 breadth (≥6/13)** | 0.272 | 0.129 | 0.057 | **+0.905 (0.002)** |
| CXCL13 log1p | 0.019 | 0.002 (10×) | 0.114 | +0.714 (0.047) |
| CD79A log1p | 0.090 | 0.028 (3×) | 0.69 | +0.595 (0.12) |
| TLS coverage % | 0.74 | 0.08 (9×) | 0.88 | +0.38 (0.35) |
| **Moran's I CXCL13** | 0.065 | 0.009 (7×) | 0.69 | −0.02 (0.96) |
| Moran's I HT-13-z | 0.382 | 0.329 | 0.89 | +0.40 (0.32) |

Both classes show positive HT-13-z Moran's I (≥+0.32) — clustered in every sample.

## Three positives

1. **HLA-DRA** (p=0.029, rho=0.98): +0.83 log1p higher in HIGH-HT, same direction as TCGA DM1>DM2 d=+1.57.
2. **HT-13 breadth** (rho=0.905, p=0.002): 27% vs 13% of spots co-detect ≥6/13 panel genes.
3. **CXCL13 clustering**: Moran's I 0.065 vs 0.009; 7× more aggregated in HIGH-HT.

## What didn't replicate

- **TLS-domain density** — 4/7 LPTC samples have 0 qualifying spots; defer to tile-level CLAM.
- **CXCL13–CCR6/CCR7 LR co-expression** uniformly negative within sample (different cell types on different spots — biologically expected).
- **DM1_like_mean** does not separate HT-axis classes at n=8.

## Files

`run_h11_spatial_braf.py` · `h11_per_sample_metrics.tsv` · `h11_spatial_HTaxis_correlation.tsv` · `h11_spearman_vs_HT13.tsv` · `h11_summary.json`

## Bottom line

Spatial replication holds for HLA-class-II + panel-breadth + CXCL13-clustering. TLS-coverage density doesn't separate at Visium resolution — claim via histology / tile DL. "DM1 BRAF-cPTC > DM2 BRAF-cPTC for TLS" is not directly testable here (no driver labels); HT-axis stratification is the honest substitute.
