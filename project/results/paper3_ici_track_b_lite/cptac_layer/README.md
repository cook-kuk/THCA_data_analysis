# Track B-lite — protein/phospho corroboration layer (Mun 2025)

**Framing:** extension of the existing Track B-lite findings (5/6 RNA analyses).
**Marathon-mode status:** Paper 3 Track B remains FROZEN per `v19_paper3_ici_track_a` until Paper 1 bioRxiv + Paper 2 A/B/C + explicit "Track B 시작". This layer corroborates the **lite** version that already ran on RNA — no Track B unfreeze.

**Created:** 2026-05-08.
**Source:** Mun et al. 2025 *Nat Commun* 16:3601, PMC12000556 — 177 PTC + 46 PDTC + 113 ATC + 119 NAT, deep proteomic (S1C, n=461 channels) and phosphoproteomic (S1D, n=217 channels).

## Headline result

**The Track B-lite RNA finding (ATC vs other: thyroid_diff d=−2.51, myeloid d=+2.53, n=415 across 4 cohorts) is replicated at protein level in Mun 2025 (n=336 tumors).**

ATC vs PTC, protein-level Cohen's d (DIAL-lite 7-module):

| Module | n_ATC | n_PTC | mean_ATC | mean_PTC | Cohen's d | p | RNA-level d (lite) | sign match |
|---|---|---|---|---|---|---|---|---|
| thyroid_differentiation | 113 | 177 | -0.51 | +0.32 | **-1.91** | 2.8e-42 | -2.51 | ✓ |
| myeloid_suppressive | 113 | 177 | +0.64 | -0.35 | **+1.52** | 1.8e-21 | +2.53 | ✓ |
| HLA_class_I | 113 | 177 | +0.36 | -0.06 | +0.56 | 1.8e-05 | (positive) | ✓ |
| IFNG_T_cell_inflamed | 113 | 177 | +0.20 | -0.07 | +0.46 | 2.1e-04 | (positive) | ✓ |
| checkpoint_exhaustion | 113 | 177 | +0.15 | -0.08 | +0.44 | 8.6e-04 | (positive) | ✓ |
| HLA_class_II | 113 | 177 | +0.19 | -0.06 | +0.29 | 0.024 | (positive) | ✓ |
| TLS_CXCL13_like | 113 | 177 | +0.07 | +0.02 | +0.10 | 0.41 | (positive, weak) | ✓ |

**Sign consistency: 7/7 modules same direction as RNA-level Track B-lite.** Magnitudes attenuated at protein level (expected — RNA-protein gene-wise r ~ 0.25–0.37 in Wang 2024 for thyroid genes); direction is preserved.

PTC → PDTC → ATC dediff trend (Spearman):

| Module | r | p_trend |
|---|---|---|
| thyroid_differentiation | **-0.69** | 3.2e-48 |
| myeloid_suppressive | **+0.51** | 2.0e-23 |
| IFNG_T_cell_inflamed | +0.19 | 5.4e-04 |
| HLA_class_I | +0.16 | 2.8e-03 |
| checkpoint_exhaustion | +0.12 | 2.3e-02 |
| HLA_class_II | +0.08 | 0.16 |
| TLS_CXCL13_like | -0.02 | 0.73 |

The dediff axis is dominated by thyroid_differentiation loss + myeloid_suppressive gain, with a smaller IFN-γ / HLA-I / checkpoint cluster — same architecture as the RNA-level Track B-lite PCA (PC1 = dediff/myeloid, PC2 = inflammation).

## What protein-level adds beyond RNA

1. **Direct ICI-relevance**: HLA-I+I peptide presentation, checkpoint molecules, and myeloid markers are all measurable at protein level in this cohort. The protein-level confirmation of HLA-I gain in ATC (d=+0.56, p=2e-5) was not directly testable at RNA level alone with the same confidence (RNA HLA-I gain can reflect tumor lymphocyte content, not antigen presentation by tumor cells).
2. **Thyroid de-differentiation as driver**: protein-level Spearman -0.69 across 336 tumors is the strongest single statistic in the layer.
3. **Coverage caveats**: HLA-A/B/C, CTLA4, PDCD1, IFNG, and several cytokines did not reach detection threshold at protein level in this cohort. Compensated by phospho-site coverage of HLA-A/B, HLA-DRA, NKX2-1, NLRC5 (S1D, 310 sites across 22 module genes).

## Files

- `../clinical_S1A.tsv` — 345 sample clinical w/ histology, OS, TERT
- `../protein_target_genes_S1C.tsv` — 31 module/panel genes × 461 channels (incl. paired normals)
- `../phospho_target_genes_S1D.tsv` — 310 phospho sites × 217 channels
- `../module_scores_per_sample.tsv` — 7 module z-scores × 336 tumor channels
- `../module_ATC_vs_PTC.tsv` — Cohen's d table above
- `../module_dediff_trend.tsv` — Spearman trend table above
- `../eight_gene_protein_per_sample.tsv` and `../eight_gene_dediff_trend.tsv` — Paper 1 cross-link
- `../summary.json` — machine-readable headline numbers

## Frozen-Track-B guardrails

- Do NOT use this as a basis for a Track B v1 launch.
- Do NOT cite any conclusion in main-text Paper 3 without explicit user "Track B 시작".
- Track B-lite already ran on RNA (5/6); this layer is **corroboration of that lite analysis**, not a new analytic effort.
- Voice-protected sections (Hook/Aim/Disc 3.1/Limitations/Cover Para 1/Q9) untouched.
