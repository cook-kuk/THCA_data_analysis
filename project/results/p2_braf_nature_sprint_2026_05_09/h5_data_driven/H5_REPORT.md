# H5 — data-driven axis discovery (BRAF_like / cPTC, N=41)

DM1=26, DM2=15. The 8 panel genes (DIO1, FOXE1, NKX2-1, PAX8, SLC5A5, TG, TPO, TSHR) are excluded from every analysis so they cannot tautologically reproduce the label. RNA-z matrix after exclusion: 50,112 genes. Same 5-fold CV splits as the v2 CLAM run.

## 1. DE DM1 vs DM2

Welch t per gene + BH-FDR. **All 50** top DEGs reach FDR<0.05; top gene **ICAM1** d=+4.63, FDR=1.4e-9. Top-50 also contains DUSP4, DUSP5 (MAPK output) and PTPRE, RUNX2, LAMB3. Only 6% of top-50 carry classical immune-prefix names (HLA/CD/IL/IFN/CCR/CCL/CXCL/…) — the curated HT-13 list captures one flavor of the same biology, not the dominant DEG mass.

## 2. Unsupervised axes (top-2000 variance, panel-excluded)

| method | best comp | pooled OOF AUC | \|r\| with HT-13 mean |
|---|---|---:|---:|
| PCA | PC2 (14% var) | 0.892 | 0.18 |
| PCA | PC3 (6.6% var) | 0.892 | 0.56 |
| NMF | **NMF2** | **0.910** | **0.70** |

Headline: **NMF2 — the single best unsupervised axis — recovers the HT-13 direction (|r|=0.70) and reaches pooled OOF AUC 0.91, within 0.02 of curated HT-13 (0.926).** The HT-immune axis emerges *unsupervised* from N=41 with 8 label-defining genes removed.

## 3. RandomForest + label-shuffle null

Top RF gene = TMEM255A (importance=0.019; null-max p95=0.021, p99=0.025; p=0.088). 4/50 top-RF genes are immune-prefix (CD300LB, HLA-DQA1, TLR2, CD58). At N=41 no single gene rises above the label-shuffle null at FDR-significance — the signal is multivariate, not one super-gene.

## 4. Random 13-gene panel null vs curated HT-13

Curated HT-13 pooled OOF AUC = 0.926.

| pool | n | mean | p95 | max | p (HT-13 ≤ pool) |
|---|---:|---:|---:|---:|---:|
| all coding | 50,112 | 0.827 | 0.972 | 1.000 | **0.208** |
| immune-like | 1,332 | 0.898 | 0.987 | 1.000 | 0.437 |

**Honest:** ~21% of random 13-gene draws from the entire coding genome and ~44% of immune-like 13-gene draws match-or-beat HT-13. The panel is not magic — it works because the underlying axis is broad. Any reasonable immune-flavored 13-gene panel works.

## 5. Pathway enrichment (Fisher, top-100 DEG)

| set | OR | p | FDR | overlap |
|---|---:|---:|---:|---:|
| HALLMARK_MAPK_SIGNALING_OUTPUT | 60.0 | 6.6e-4 | **0.008** | DUSP4, DUSP5 |
| IL6_JAK_STAT3 | 28.1 | 0.037 | 0.18 | JAK3 |
| TNFA_NFKB / INFLAMMATORY / IFNG | 12-18 | 0.06-0.09 | 0.18-0.21 | ICAM1 |

MAPK-output is the only FDR<0.05 pathway. Inflammation/IFNG reach nominal p<0.10 via ICAM1.

## Bottom line — cherry-pick rebuttal

1. **HT-axis IS data-driven.** Best NMF axis on panel-excluded top-2000-variance genes has |r|=0.70 with HT-13 and pooled OOF AUC 0.91 — within 0.02 of the curated panel.
2. **Panel-specificity p = 0.208** vs all-coding random 13-gene panels (0.437 vs immune-like). The HT-13 panel is not unique; the underlying axis is broad.
3. **DM1 vs DM2 within BRAF_like/cPTC = inflammation × MAPK-output**, not pure B-cell TLS. Top DEGs ICAM1 + DUSP4/5; only HALLMARK_MAPK_OUTPUT clears FDR. Aligns with deconv v9-v12 (MAPK drives 8-gene methylation) → supports Paper 2 "two-axis convergence" framing.

**Reviewer-safe phrasing.** *"In BRAF_like/cPTC (N=41), unsupervised NMF on panel-excluded top-2000-variance genes recovers an axis with |r|=0.70 to HT-13 and pooled OOF AUC 0.91 — within 0.02 of the curated panel. HT-13 is not uniquely powerful (panel-specificity p=0.21) because the underlying axis is broad and dominated by inflammation+MAPK-output (ICAM1, DUSP4, DUSP5; HALLMARK_MAPK_OUTPUT FDR=0.008). The HT-axis is data-driven, not curated."*

## Files
- `h5_de_top50.tsv` — 50 DEGs (all FDR<0.05)
- `h5_components.tsv` — PCA(5) + NMF(5) + AUC + corr(HT-13)
- `h5_perm_null.tsv` + `perm_null_aucs_*.npy` — 1000 random panels per pool
- `h5_rf_importance.tsv` — top-50 RF + null-max
- `h5_pathway_enrich.tsv` — Fisher Hallmark-like enrichment
- `run_h5_data_driven.py` — full pipeline
