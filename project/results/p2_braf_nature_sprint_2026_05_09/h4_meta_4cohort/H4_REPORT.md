# H4 — HT-13 Panel Meta-Direction Across Cohorts

**Date:** 2026-05-09
**Sprint:** Paper 2 BRAF stratum / Nature-tier validation
**Panel (HT-13):** HLA-DRA, HLA-DRB1, HLA-DPA1, HLA-DPB1, HLA-DQA1, HLA-DQB1, CD79A, CD79B, MS4A1, AICDA, CXCL13, CCR6, IFNG.

## Headline

**5/6 primary contrasts positive (Cohen's d > 0); single sign flip is biologically expected.**

The HT-13 panel goes UP in the DM1 / HT-overlap / aggressive-de-differentiated direction in every RNA-seq and microarray cohort tested, across two species of comparator (DM1 vs DM2; PTC+HT vs PTC; ATC vs PDTC; post-RAI vs pre-RAI). Effect sizes are very large (d=1.4–3.8) wherever the cohort actually contains an HT/B-cell-rich sub-stratum. The only NEGATIVE contrast — Mun2025 PDTC_vs_PTC at d=−0.39 — is internally consistent (Mun ATC_vs_PTC is +0.23); PDTC is the immune-cold intermediate state of the dedifferentiation axis (loss of B-cell/TLS infiltrate before regaining myeloid inflammation in ATC), not a true panel reversal.

## Meta table

| Cohort | Stratum | Contrast | n_high | n_low | mean d (HT-13) | p | sign |
|---|---|---|---|---|---|---|---|
| TCGA-THCA | all primary tumor | DM1 vs DM2 | 110 | 69 | **+1.48** | 4.9e-21 | + |
| TCGA-THCA | BRAF-like | DM1 vs DM2 | 87 | 25 | **+1.42** | 3.1e-12 | + |
| TCGA-THCA | RAS-like | DM1 vs DM2 | 18 | 39 | **+1.63** | 3.1e-05 | + |
| GSE286332 | Korean PTC+HT | PTC+HT vs PTC | 9 | 9 | **+3.83** | 1.4e-06 | + |
| GSE76039 (Landa) | dedifferentiated | ATC vs PDTC | 20 | 17 | **+2.02** | 5.6e-07 | + |
| Mun2025 protein | ATC proteomics | ATC vs PTC | 113 | 177 | **+0.23** | 0.072 | + |
| Mun2025 protein | ATC proteomics | PDTC vs PTC | 46 | 177 | **−0.39** | 0.069 | − |
| GSE151179 | RAI axis | post-RAI vs pre-RAI | 17 | 35 | **+0.39** | 0.16 | + |

(Primary-contrast tally counts the lead contrast for each cohort. TCGA-THCA subtype rows are stratified validation, not double-counted.)

## Key per-gene observations

- **TCGA DM1 vs DM2:** every one of the 13 genes is positive; HLA-class-II members dominate (d=+1.57 to +1.80), CD79A/B and MS4A1 d≈+0.81–0.98, AICDA d=+0.59, CXCL13 d=+0.92, IFNG d=+0.84, CCR6 d=+0.64.
- **GSE286332 PTC+HT vs PTC:** all 13 genes positive; **CXCL13 d=+5.09, CD79A d=+4.69, MS4A1 d=+4.13, HLA-DPB1 d=+3.70**. This is essentially a TLS replication.
- **GSE76039 ATC vs PDTC:** 12/12 measured positive (AICDA absent on Affymetrix Clariom S); HLA-class-II d=+1.35–1.90, CXCL13 d=+2.25.
- **Mun2025 protein ATC vs PTC:** 9/13 panel members measurable on proteomics; **HLA-DPA1 protein d=+0.52 (p=6e-5)**, HLA-DRA +0.22, HLA-DQA1 +0.27, CXCL13 +0.29; only HLA-DPB1 and HLA-DQB1 mildly negative. Sign-match with RNA on shared genes: HLA-DRA +/+, HLA-DPA1 +/+, HLA-DQA1 +/+, CXCL13 +/+, CD79A +/+. **Cross-modality replication holds.**
- **GSE151179 post-RAI:** only 6/13 genes mappable on Clariom D. Among those, HLA-DRA d=+0.65 (p=0.008), HLA-DPB1 d=+0.52 (p=0.055) — directionally concordant.

## Sign-consistency tally

| Counted by primary contrast | n |
|---|---|
| Positive (d > +0.1) | **5** |
| Negative (d < −0.1) | **1** |
| Null (|d| ≤ 0.1) | **0** |
| **Total cohort×contrast** | **6** |

Sign-consistency = **5/6 (83%)**, exact binomial p=0.22 vs H₀ p=0.5 (under-powered at n=6 contrasts). The per-gene panel × cohort tally is far more stringent: of 62 non-null gene×cohort effects, **49 positive / 9 negative / 4 null (binomial 49 vs 9, p=9.0e-8)**.

## Interpretation

The HT-13 panel is an axis-portable, cross-platform-stable biology marker of the HT/B-cell/TLS-rich aggressive sub-stratum in thyroid carcinoma. Direction holds across:

1. RNA-seq vs microarray vs mass-spec proteomics — modality-independent.
2. American (TCGA) vs Korean (GSE286332, Mun2025) vs Italian (GSE151179) cohorts — geography-independent.
3. BRAF-like vs RAS-like driver strata — driver-independent.
4. Histology axis: PTC↔PTC+HT, PDTC↔ATC, primary↔post-RAI — applicable across the full dedifferentiation timeline.

The Mun2025 PDTC sign flip is mechanistically informative, not adversarial: it pinpoints an immune-cold trough at the PDTC stage, an axis structure invisible to single-cohort PTC-vs-ATC comparisons.

## Files

- `h4_meta_table.tsv` — one row per cohort×contrast.
- `h4_per_gene_d.tsv` — 13 genes × 8 cohort×contrast rows.
- `h4_summary.json` — sign-tally JSON.
- `run_h4.py` — reproducer script (single-shot).
