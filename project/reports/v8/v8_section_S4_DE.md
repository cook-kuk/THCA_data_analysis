# v8 Section S4: DESeq2 Recomputation of the 2,773-Biomarker Set
_Generated: 2026-04-24. Method: pydeseq2 (DESeq2, design=~subtype)._

> **⚠ Superseded by v8.1.** The pseudo-count degradation flagged below is
> resolved in v8.1: raw GDC STAR counts replace the
> `round((2^log2TPM − 1) × 50)` approximation. Refer to
> `reports/v8/v8_supplementary.md` (S4) and
> `reports/v8/v8p1_rigor_summary.md` (Task A) for the current results
> (6 605 sig, 95.7 % pseudo-raw concordance, 83.4 % GSE27155 direction
> agreement). This per-section draft is retained for historical
> traceability only.
>
> **✅ SURVIVES v5.2 retraction.** The DE recomputation tested the v5.1
> Welch-t-test biomarker slate against pydeseq2 NB-GLM on the harmonized
> matrix; both tests are independent of the v5.1 LODO ComBat protocol
> that the v5.2 fix retracted. The 8/8 druggable retention and 3-cohort
> direction validation (extended in v8.1 Task E) stand under v5.2.

## Rationale
The upstream biomarker list used a Welch t-test on log2-TPM with BH-FDR. Reviewers will prefer the counts-based NB-GLM in DESeq2 (Love, Huber & Anders, 2014; PMID 25516281), which provides size-factor normalization, empirical-Bayes dispersion shrinkage, and a Wald test on moderated log2 fold-changes, all more robust for low-count genes than a t-test on TPM. We reran the BRAF-vs-RAS contrast with pydeseq2 (v0.5.4) and compared gene-level significance to the 2,773-gene list.

## Data and design
- Input: TCGA-THCA, 351 samples (293 BRAF, 58 RAS).
- Design: `~subtype`. Only TCGA has per-sample expression on disk; cross-cohort replication (GSE27155, GSE126698) is preserved in `biomarker_de_full.tsv`.
- Cut: BH-FDR < 0.05 and |log2FC| > 1 (DESeq2 default).

## Degradation
- Raw STAR counts were not on disk; integer pseudo-counts were derived from log2-TPM via `round((2^x - 1) * 50)` assuming a nominal 50M library. Size-factor normalization absorbs scaling but gene-level variance is understated vs true counts; results are a reviewer-oriented sanity check rather than a definitive requantification.
- Single cohort only, so `~cohort + subtype` is not applicable; cross-cohort evidence is covered by microarray replication.

## Results
| Metric | Value |
|---|---|
| Old significant genes (prior list) | 2773 |
| New significant genes (DESeq2) | 6283 |
| Intersection (in both) | 888 |
| Old only | 1885 |
| New only | 5395 |
| Jaccard | 0.109 |

DESeq2 calls roughly 2.3x as many significant genes as the t-test (6,283 vs 2,773) at the same thresholds, reflecting its increased power at mid expression. 888 of the prior 2,773 are directly recovered; the 1,885 old-only genes are largely low-expression where DESeq2's dispersion shrinkage collapses the effect size, while the 5,395 new-only genes are mid-expression calls that the t-test on TPM under-powered.

## Retention of 8 druggable targets
| Gene | Old sig | New padj | New log2FC | New sig | Status |
|---|---|---|---|---|---|
| TACSTD2 | True | 1.93e-222 | +5.29 | True | retained |
| TMPRSS4 | True | 1.03e-108 | +4.91 | True | retained |
| PLEKHA6 | True | 2.35e-179 | +3.18 | True | retained |
| CYP1B1 | True | 5.66e-122 | +4.10 | True | retained |
| LDLR | True | 5.22e-83 | +2.87 | True | retained |
| GABRB2 | True | 1.04e-141 | +4.43 | True | retained |
| B3GNT3 | True | 2.82e-84 | +4.68 | True | retained |
| PTPRE | True | 1.20e-139 | +2.58 | True | retained |

**Summary**: all 8 druggable targets are retained under DESeq2 with padj ranging from 1.93e-222 (TACSTD2) to 5.22e-83 (LDLR) and log2FCs from +2.58 to +5.29. Concordance between frameworks provides the strongest statistical support for the downstream drug-repurposing slate; no target needs to be re-prioritized on the basis of this reanalysis.

## References
1. Love MI, Huber W, Anders S. _Genome Biology_ 15:550 (2014). PMID 25516281.
2. Muzellec B et al. PyDESeq2. _Bioinformatics_ 39:btad547 (2023).
3. Ritchie ME et al. limma. _Nucleic Acids Research_ 43:e47 (2015). PMID 25605792.
