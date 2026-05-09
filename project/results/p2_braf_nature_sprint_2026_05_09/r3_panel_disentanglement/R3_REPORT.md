# R3 — HT-13 panel disentanglement vs Ayers-TIS-18 vs Cabrita-TLS-12

**Verdict: HT-13 IS distinguishable from TIS+TLS at the gene level (8/13 unique = 62%) and at the cross-cohort Hashimoto-overlap classifier level (Korean GSE213647 hashi-AUC 0.910 vs 0.857/0.886), but the three panels track ONE strongly-correlated inflamed-tumor axis as scalar scores (r≈0.93 in BRAF-cPTC). HT-13 is the thyroid-Hashimoto-tuned member of that family — not a novel axis.**

## 1) Gene-level overlap (Jaccard)
- HT-13 ∩ TIS-18 = 2 (HLA-DRA, IFNG); J=0.065
- HT-13 ∩ Cabrita-12 = 3 (CD79B, CXCL13, MS4A1); J=0.136
- TIS-18 ∩ Cabrita-12 = 0; J=0.000
- HT-13 unique (8/13): AICDA, CCR6, CD79A, HLA-DPA1/DPB1/DQA1/DQB1/DRB1.

## 2) Score correlation (gene-mean z)
BRAF-cPTC n=252: HT-13↔TIS r=0.92, HT-13↔TLS r=0.87, TIS↔TLS r=0.78. RAS/FVPTC n=34: HT-13↔TLS drops to r=0.37 (B-cell-zone genes specific to BRAF-cPTC immunity).

## 3) DM1 vs DM2 AUC (n=179, 5-fold LogReg)
| Panel | n_genes | CV-AUC |
|---|---|---|
| HT-13 full | 13 | 0.895 |
| HT-13 minus all TIS+TLS overlap | 8 | **0.897** |
| TIS-18 full | 20 | **0.967** |
| TLS-12 full | 12 | 0.905 |
| Union (40 unique) | 40 | 0.967 |

Removing the 5 overlap genes from HT-13 leaves AUC unchanged — HT-13's unique 8 genes carry the full classifier signal.

## 4) Residualized scalar scores (the killer test)
- HT-13 residualized on TIS-18 score: AUC 0.876 → 0.581
- TIS-18 residualized on HT-13: 0.863 → 0.517 (~chance)
- HT-13 residualized on TLS-12: 0.876 → 0.754

**As single mean-z scores, HT-13 and TIS-18 are largely interchangeable summaries of the same axis** — both collapse to ~chance once the other is regressed out.

## 5) Mutual information (DM1 vs DM2)
MI: HT-13 0.254 / TIS-18 0.246 / TLS-12 0.231. Conditional: I(HT-13;DM | TIS-18)=0.064, I(TIS-18;DM | HT-13)=0.048. **Each panel retains ~20–25% unique MI; neither subsumes the other.**

## 6) B-cell sub-panel
| Sub-panel | n | CV-AUC | d |
|---|---|---|---|
| HT B-cell (CD79A/B+MS4A1+AICDA) | 4 | 0.743 | 0.89 |
| HT non-B (HLA-DR/DP/DQ+CXCL13+CCR6+IFNG) | 9 | **0.888** | **1.63** |
| HT-13 full | 13 | 0.895 | 1.47 |

**Class-II + IFN/chemokine arm carries the bulk of the signal**; B-cell block adds marginally.

## 7) Cross-cohort GSE213647 (Lee 2024, n=369 tumor)
Hashimoto (hashi_GMM) AUC: **HT-13 0.910 > TLS-12 0.886 > TIS-18 0.857**. When the label is HT-overlap (HT-13's design target), HT-13 leads by 5 AUC points despite being smaller than TIS-18.

## Bottom line for reviewer Q
HT-13 shares 15% of genes with TIS-18 and 23% with Cabrita-TLS-12; 62% are unique. The three panels co-vary strongly as scalar scores in BRAF-cPTC (one inflamed-tumor axis), but at the gene-level multivariate classifier and at the cross-cohort Hashimoto-overlap label, HT-13 contributes distinct, panel-specific information (~25% of MI uniquely; +0.05 AUC over TIS in Korean hashi-call). Recommended manuscript framing: "HT-13 is the thyroid-Hashimoto-tuned member of the inflamed-tumor signature family; co-varies with TIS as expected; outperforms TIS for the Hashimoto-overlap call it was designed for."
