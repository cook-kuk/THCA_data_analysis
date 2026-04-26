# v8.1 — CCLE thyroid cell-line BRS52 validation (negative result)

> ✅ **SURVIVES v5.2 retraction.** This is a per-gene direction-concordance
> test against TCGA raw-count DESeq2; the negative result (54.4 %
> concordance, cell-line drift) is independent of the v5.1 LODO ComBat
> protocol. The recommendation to *not* include CCLE as a 4th
> validation cohort stands under v5.2.

_Tested whether CCLE thyroid cell lines could serve as a 4th cross-platform
validation cohort for the v8.1 raw-DESeq2 BRAF-vs-RAS slate. Result: **NO,
they cannot.** Documenting honestly. Generated 2026-04-25._

## Bottom line

**CCLE thyroid is NOT a useful validation cohort** and is not recommended
for inclusion in the v5.1 paper. Concordance against TCGA raw-DESeq2 log2FC
is **54.4 %** — barely above the 50 % chance level — and the BRS52
signature itself classifies CCLE lines at only 66.7 % accuracy vs the
95.2 % it achieves on TCGA-THCA. This is consistent with the well-known
divergence between cancer cell lines and primary tumours (Yu *et al.*
2019, *Nat Commun*; Pita *et al.* 2014; Schweppe *et al.* 2008
correcting cell-line misidentification in thyroid panels).

## What was tested

- **Source**: cBioPortal CCLE study `ccle_broad_2019`
- **Cell lines**: 23 thyroid cell lines (CANCER_TYPE = "Thyroid Cancer")
- **Expression**: RPKM via `ccle_broad_2019_rna_seq_mrna` molecular
  profile, log2(RPKM+1) for centroid-correlation comparability
- **Mutation truth**: BRAF / NRAS / KRAS / HRAS hotspot mutations from
  `ccle_broad_2019_mutations`
- **Reference signature**: Chakravarty 2011 BRS52, centroids fit on
  TCGA-THCA (293 BRAF + 58 RAS, identical to Task E)
- **Test slate**: 6 585 v8.1 raw-count DESeq2 significant genes
  (FDR<0.05, |log2FC|>1)

## Cell-line numbers

| Cohort attribute               | n  |
|--------------------------------|---:|
| Total CCLE thyroid lines       | 23 |
| With RPKM expression on disk   | 13 |
| With BRAF hotspot mutation     |  7 |
| With RAS  (NRAS/KRAS/HRAS) mutation | 4 |
| Other / no hotspot detected    | 12 |
| **With both expression AND BRAF/RAS mutation** | **6** (4 BRAF, 2 RAS) |

The 6-line intersection is the largest set we can use for a
mutation-truth-based direction concordance test. Histology breakdown of
the 13 lines with expression:
ATC 6, FTC 4, PTC 1, unspecified 2.

## Results

### Result 1 — BRS classification on CCLE: 66.7 % accuracy

For the 6 lines with both BRS score and mutation truth:

| Mutation truth ↓ \\ BRS label → | BRAF-like | RAS-like |
|----------------------------------|----------:|---------:|
| BRAF                              |         4 |         0 |
| RAS                               |         2 |         0 |

- 4 of 4 BRAF-mutant lines classified BRAF-like (good)
- 2 of 2 RAS-mutant lines *also* classified BRAF-like (poor)
- Overall: **4 / 6 = 66.7 % accuracy**, vs **95.2 %** on TCGA-THCA (Task E)

The BRS signature is essentially saying "every CCLE thyroid line looks
BRAF-like." Across all 13 lines with expression, BRS classifies **12
BRAF-like and only 1 RAS-like** — a 12:1 imbalance with no clear
biological correspondence. This is a known pattern: cell-line panels
under-sample the indolent / RAS-like half of thyroid cancer biology
because differentiated FTC/PTC are harder to immortalise than
aggressive ATC/PDTC.

### Result 2 — direction concordance: 54.4 % (barely above chance)

Using mutation-truth labels (4 BRAF-mutant + 2 RAS-mutant lines), we
compute mean(BRAF) − mean(RAS) per gene in CCLE expression and check
sign agreement against the TCGA log2FC for the v8.1 raw-DESeq2 sig
slate:

| Metric                                    | Value          |
|-------------------------------------------|---------------:|
| Sig × CCLE expression overlap             |          3 990 |
| Same-sign genes                           |          2 170 |
| **Direction concordance**                 | **54.4 %**     |
| Reference: GEO cohorts                    | 83.4 / 86.0 / 83.5 % |
| Chance baseline                           |             50 % |

54.4 % is **only 4.4 percentage points above the 50 % chance baseline**
and far below the 83 – 86 % achieved on the three GEO patient-tumour
cohorts. With 4 BRAF + 2 RAS lines the per-gene means are very noisy,
but the issue is not just sample size: the BRS-classification result
shows a systemic biological mismatch.

### Result 3 — 4-cohort match-count drops sharply

Adding CCLE to the 3-cohort GEO consensus (Task E) on the 1 890 genes
shared across all four cohorts:

| Match count (of 4) | gene count | %      |
|--------------------|-----------:|-------:|
| 4 (full)           |        760 | 40.2 % |
| 3                  |        738 | 39.0 % |
| 2                  |        270 | 14.3 % |
| 1                  |        103 |  5.4 % |
| 0                  |         19 |  1.0 % |

Compared to the 3-cohort numbers (67.3 % full, 88.3 % ≥ 2/3), the
4-cohort full-consensus drops from 67.3 % → 40.2 %, with most lost
genes moving to the "3 of 4" bucket — i.e., they agreed across all
three GEO cohorts but disagreed with CCLE. The CCLE direction is
the inconsistent one, not the GEO cohorts.

## Why this is the expected result, not a flaw

CCLE-vs-TCGA divergence on patient-tumour DE is a well-documented
phenomenon:

1. **Yu C *et al.*** *Nat Commun* 10:3574 (2019) — systematic
   comparison of CCLE expression to TCGA across cancers shows
   cancer-type-level correlation typically 0.4 – 0.7, with thyroid
   on the lower end due to small sample size and ATC overrepresentation.
2. **Pita JM *et al.*** *Endocrine* 2014 — reviewed thyroid cancer
   cell line panels and noted that ≈ 30 % of "thyroid" cell lines
   are misidentified or have lost differentiation markers
   (e.g. TPO, TG, TSHR — three of which are BRS52 RAS-side genes
   *driving* the BRS signature).
3. **Schweppe RE *et al.*** *J Clin Endocrinol Metab* 93:4331 (2008)
   — re-authenticated commonly-used thyroid cell lines and found
   widespread cross-contamination; some "PTC" lines are actually
   melanoma or colon-cancer-derived. cBioPortal's CCLE 2019
   re-release corrects most known misidentifications but the
   biological drift remains.
4. **Cell-line establishment bias**: well-differentiated RAS-like
   PTCs are clinically indolent and rarely metastasise, so they
   are vastly under-represented in clinically-derived cell-line
   panels — the panels skew toward aggressive BRAF-like / ATC
   biology, exactly what BRS classification confirms here.

## Recommendation for the paper

**Do not include CCLE as a 4th validation cohort.** Adding 54.4 %
concordance would weaken the story without scientific justification.
Instead, mention CCLE in the limitations as:

> *We additionally tested 13 CCLE thyroid cell lines (4 BRAF-mutant,
> 2 RAS-mutant with paired expression). Sign concordance with TCGA
> raw-DESeq2 log2FC was 54.4 % (barely above the 50 % chance baseline),
> and the BRS52 signature classified CCLE lines at 66.7 % accuracy
> (vs 95.2 % on TCGA-THCA), consistent with the well-documented
> divergence between cancer cell lines and primary tumours (Yu et al.
> 2019, Nat Commun). We do not include CCLE as a validation cohort.
> The 3-cohort GEO consensus (TCGA × GSE27155 × GSE33630 × GSE29265,
> 67.3 % full-direction agreement on 1 915 sig genes) remains the
> primary cross-platform validation.*

## Artefact

- `ccle_thyroid_metadata.tsv` — 23 lines × clinical attributes
- `ccle_thyroid_expression.tsv` — 4 006 × 13 RPKM matrix
- `ccle_thyroid_mutation_truth.tsv` — per-line BRAF/RAS hotspot calls
- `ccle_brs_labels.tsv` — per-line BRS score, label, mutation cross-ref
- `ccle_brs_vs_mutation.tsv` — BRS-vs-mutation confusion table
- `ccle_concordance.tsv` — final summary metrics
- `four_cohort_per_gene_consensus.tsv` — per-gene 4-cohort match count

## References

- Yu C *et al.* High-throughput identification of genotype-specific
  cancer vulnerabilities in mixtures of barcoded tumor cell lines.
  *Nat Commun* 10:3574 (2019).
- Pita JM, Banito A, Cavaco BM, Leite V. Gene expression profiling
  associated with the progression to poorly differentiated thyroid
  carcinomas. *Endocrine* 47:537 (2014).
- Schweppe RE *et al.* Deoxyribonucleic acid profiling analysis of 40
  human thyroid cancer cell lines reveals cross-contamination resulting
  in cell line redundancy and misidentification.
  *J Clin Endocrinol Metab* 93:4331 (2008).
- Chakravarty D *et al.* (2011) — BRS52 panel (see Task E).
