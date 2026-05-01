# P3 — GSE286332 PTC vs PTC+HT (CRITICAL gate result)

**Date:** 2026-04-30 PM
**Status:** ✅ STRONG GO — 시나리오 2 enabled. Hashimoto-overlap PTC IS molecularly distinct, extremely so.

---

## ★ TL;DR (one line)

GSE286332 (n=9 PTC vs 9 PTC+HT, Korean RNA-seq) shows **10,380 DEGs (padj < 0.05)**, 8-gene RAI score Cohen d = **−1.60** (p=0.008), HLA-II module Cohen d = **+3.65** (p=0.0004), and Hallmark IFN-γ Response / KEGG Type I Diabetes Mellitus / Reactome TCR signaling all top-FDR up — a textbook autoimmune-driven dedifferentiation phenotype. Reviewer Q "is the autoimmune-PTC sub-axis real?" answered YES with high effect size despite n=18.

---

## 1. Differential expression (PyDESeq2 Wald, Benjamini–Hochberg)

| | n |
|---|---|
| Total genes tested | 29,672 |
| **DEGs (padj < 0.05)** | **10,380** |
| Up in PTC+HT | 6,004 |
| Down in PTC+HT | 4,376 |

### Top up in PTC+HT (B-cell + tertiary lymphoid signature)

| gene | log2FC | padj |
|---|---|---|
| IGHV3-66 | 7.25 | 1.4e-35 |
| BLK | 5.97 | 6.9e-33 |
| IGHV3-13 | 7.10 | 3.2e-32 |
| IGHV3-16 | 7.86 | 8.5e-32 |
| IGHV2-70 | 8.04 | 5.1e-28 |
| IGKV1-8 | 7.94 | 1.4e-27 |
| EOMES | 4.85 | 9.4e-25 |
| HLA-DOB | 5.43 | 7.3e-24 |
| IGLV6-57 | 6.34 | 1.7e-23 |
| (multiple Ig V/J/C chains) | | |

→ B-cell receptor + TLS (tertiary lymphoid structure) signature, classic Hashimoto.

### Top down in PTC+HT
RIMS4, CHST1, KCNC3, SNCA, BFSP1 (neuronal/secretory loss); KLK2, HSD17B3 (hormone metabolism loss); PON2, TUSC1.

---

## 2. GSEA pre-ranked

### MSigDB Hallmark — UP top 10 (immune cascade)
| Term | NES | FDR |
|---|---|---|
| Allograft Rejection | +2.12 | 0 |
| E2F Targets | +1.98 | 0 |
| G2-M Checkpoint | +1.96 | 0 |
| **Interferon Gamma Response** | **+1.80** | **1.9e-4** |
| Inflammatory Response | +1.75 | 1.3e-4 |
| IL-6/JAK/STAT3 | +1.75 | 1.5e-4 |
| Complement | +1.73 | 2.2e-4 |
| TNF-α/NF-κB | +1.65 | 1.2e-3 |
| Mitotic Spindle | +1.63 | 1.7e-3 |
| IL-2/STAT5 | +1.59 | 2.8e-3 |

### MSigDB Hallmark — DN top (metabolic dedifferentiation)
| Term | NES | FDR |
|---|---|---|
| Fatty Acid Metabolism | −1.85 | 0 |
| Adipogenesis | −1.66 | 0.015 |
| Oxidative Phosphorylation | −1.34 | 0.12 |

### KEGG — UP top 10 (autoimmune cluster)
| Term | NES | FDR |
|---|---|---|
| **Type I diabetes mellitus** | +1.92 | **0** |
| Staphylococcus aureus infection | +1.91 | 0 |
| Epstein-Barr virus infection | +1.89 | 0 |
| Osteoclast differentiation | +1.88 | 0 |
| **B cell receptor signaling** | +1.88 | 0 |
| **NF-κB signaling** | +1.88 | 0 |
| Asthma | +1.88 | 0 |
| Viral myocarditis | +1.88 | 0 |
| Toxoplasmosis | +1.86 | 0 |
| Hematopoietic cell lineage | +2.02 | 0 |

### Reactome — UP top (TCR + MHC-II)
- Phosphorylation Of CD3 And TCR Zeta
- ZAP-70 Translocation To Immunological Synapse
- Generation Of Second Messenger Molecules (TCR signal)
- Immunoregulatory Interactions Lymphoid–non-Lymphoid
- IL-10 Signaling

---

## 3. 8-gene RAI panel score

PTC+HT shows **strong dedifferentiation** in the 8-gene transcriptional differentiation axis.

| metric | PTC (n=9) | PTC+HT (n=9) |
|---|---|---|
| 8-gene Z-mean | +0.50 | −0.50 |
| Cohen's d | | **−1.60** |
| MW p | | **0.008** |
| Welch's t p | | 0.005 |

### Per-gene (PTC+HT vs PTC, log2 FPKM)
| gene | mean_PTC | mean_PTC+HT | Cohen d | MW p |
|---|---|---|---|---|
| **PAX8** | 8.68 | 7.73 | **−2.32** | 7.9e-4 |
| **NKX2-1** | 7.02 | 6.11 | **−1.92** | 3.6e-3 |
| **FOXE1** | 7.09 | 6.16 | **−1.75** | 6.2e-3 |
| **TG** | 12.17 | 11.29 | **−1.64** | 4.7e-3 |
| **TSHR** | 7.26 | 6.60 | **−1.60** | 6.2e-3 |
| TPO | 8.66 | 7.61 | −1.11 | 0.064 |
| DIO1 | 5.25 | 4.66 | −0.56 | 0.22 |
| SLC5A5 | 2.94 | 3.31 | +0.25 | 0.54 |

→ Transcription factor backbone (PAX8/NKX2-1/FOXE1) collapses; iodide transporter (SLC5A5) actually preserved → **dedifferentiation is TF-driven, not transporter-loss-driven**.

---

## 4. DM1/DM2 prediction (TCGA-trained centered-profile classifier)

| group | DM call | n |
|---|---|---|
| PTC | DM2 | 9/9 |
| PTC+HT | DM2 | 9/9 |

But **P(DM1) drifts strongly within DM2**:

| group | mean P(DM1) |
|---|---|
| PTC | 0.18 |
| PTC+HT | 0.05 |

MW p = **0.0036** → PTC+HT pushes deeper into DM2 territory. Consistent with **"Hashimoto-overlap PTC = extreme DM2 sub-cluster"**.

---

## 5. HLA-I / HLA-II module scores

| module | mean PTC | mean PTC+HT | Cohen d | MW p |
|---|---|---|---|---|
| HLA-I (HLA-A/B/C, B2M, TAP1/2, PSMB8/9, NLRC5) | −0.74 | +0.74 | **+2.34** | 1.5e-3 |
| **HLA-II (HLA-DRA/B1, DPA1/B1, DQA1/B1, DMA/B, CIITA, DOB)** | −0.83 | +0.83 | **+3.65** | 4.1e-4 |

→ HLA-II up-regulation is **near-perfect group separation** (d=+3.65 is exceptionally large for n=18). Direct mechanistic link to autoimmune CD4+ T-cell activation.

---

## 6. Concordance with our K2/Lee Q12 finding (17% Hashimoto-like in TCGA-THCA)

GSE286332 reproduces the **same direction** as our Q12 TCGA Hashimoto-like sub-cluster:
- HLA-II up (already noted in K2 + Lee2024 cohorts)
- IFN-γ + B-cell + TCR (matches Q12 immune subtype)
- 8-gene RAI score down (consistent with DM2 ⊃ Hashimoto-like)

3-cohort meta is now feasible: **TCGA-THCA Hashimoto-like (n≈80)** + **K2 Hashimoto-like portion** + **GSE286332 PTC+HT (n=9)** → all in same direction.

---

## 7. Paper venue impact (per Yu meeting calibration)

**Pre-P3:** Sci Rep (IF 4) baseline. Cell Rep Med (IF 14) reach.
**Post-P3:** **Cell Rep Med / JCI Insight reach is realistic.** Three pillars now in hand:
1. K2+Lee HLA imputation (n=890) — Pan-Asian cohort
2. **GSE286332 PTC+HT molecular dissection (10,380 DEGs, IFN-γ + HLA-II)** ← P3
3. P1 driver-mRNA neutrality (BRAF d=−0.04; orthogonal to mutation) — answers reviewer Q4

→ Recommend escalating venue plan **before** Bundang Graves' arrival. If Bundang Graves' lands, target Nat Commun.

---

## 8. Day 2 PM action items

- ✅ P3 done
- Next: P2 (power table + Plan B map for Bundang outreach branches)
- Next: P4 (pan-genome MAD top-5000 vs TIERA67) — reviewer Q3 robustness sanity
- Next: P5 (Finding 2 8-gene vs HLA-II autocorrelation; residualize HLA-II out and re-test 8-gene)

---

## 9. Outputs

- `results/p3_gse286332/deg_ptcht_vs_ptc.tsv` — full DEG table (29,672 genes)
- `results/p3_gse286332/gsea_MSigDB_Hallmark_2020.tsv` — 50 terms
- `results/p3_gse286332/gsea_KEGG_2021_Human.tsv` — 311 terms
- `results/p3_gse286332/gsea_Reactome_2022.tsv` — 1,403 terms
- `results/p3_gse286332/8gene_panel_per_sample.tsv` — sample × gene panel + RAI score
- `results/p3_gse286332/8gene_per_gene_compare.tsv` — per-gene Cohen d table
- `results/p3_gse286332/dm12_predictions.tsv` — 18-sample DM1/DM2 calls
- `results/p3_gse286332/hla_module_scores.tsv` — HLA-I/II per sample
- `results/p3_gse286332/P3_summary.json` — machine-readable summary
