# Supplementary Tables — Manuscript v8 Draft Structure

**Date:** 2026-05-03 (marathon mode prep)
**Format:** Cell Press supplementary table specs (Excel `.xlsx` / TSV `.tsv`, 1 sheet per table)

---

## Suppl Table 1 — Cohort assembly

| Cohort | Source / Study | n | Modality | Population | DM cluster avail | HLA avail | Use case |
|---|---|---|---|---|---|---|---|
| TCGA-THCA | TCGA-Cancer Network 2014 | 500 | RNA-seq Illumina | EUR-dominant | DM1=140, DM2=360 | 4-digit imputed (HLA-LA) | Discovery, classifier training |
| K2 (PRJEB11591) | Yoo SK 2016 SNU-GMI | 260 (235 valid HLA) | RNA-seq Illumina | Korean | DM1=14, DM2=246 | arcasHLA 4-digit (5/2) | Korean PTC validation |
| Lee 2024 (GSE213647) | Lee SE et al. 2024 Macrogen | 632 (630 valid HLA) | RNA-seq Illumina | Korean | DM proxy via 8-gene | arcasHLA 4-digit | Korean PTC replication |
| GSE286332 | Lim DW 2025 Dongguk Univ | 18 (9 PTC + 9 PTC+HT) | RNA-seq NovaSeq X | Korean | DM2 100% (all 18) | arcasHLA 4-digit (5/2) | PTC vs PTC+HT discovery |
| Chu 2018 J Med Genet | Chu X et al. 2018 | 2,958 (1,468 GD / 1,490 ctrl) | SNP2HLA Pan-Asian panel | Han Chinese | NA | summary stats published | Pan-Asian forest replication |

**Source file:** `project/results/p2_pillar1_forest/cohort_assembly.tsv` (built 2026-05-04, 5 cohorts × 10 columns)

---

## Suppl Table 2 — TIERA67 7-category 67-gene candidate pool

| Category | Genes (n) | Members |
|---|---|---|
| TDS_core | 16 | DIO1, DIO2, DUOX1, DUOX2, FOXE1, GLIS3, NKX2-1, PAX8, SLC26A4, SLC5A5, SLC5A8, TG, THRA, THRB, TPO, TSHR |
| MAPK_output_ERK | 10 | DUSP4, DUSP5, DUSP6, SPRY1, SPRY2, SPRY4, ETV4, ETV5, PHLDA1, FOSL1 |
| **Driver_anchor** | 12 | BRAF, NRAS, HRAS, KRAS, RET, NTRK1, NTRK3, ALK, PAX8, PPARG, TERT, EIF1AX |
| Aggressive_marker | 10 | TP53, CDKN2A, CDKN2B, PIK3CA, AKT1, PTEN, ATM, CTNNB1, APC, MSH2 |
| Dediff_invasion | 10 | VIM, ZEB1, ZEB2, SNAI1, SNAI2, TWIST1, CDH1, CDH2, MMP9, LOX |
| Immune_stromal_light | 5 | CD274, CD8A, FOXP3, IDO1, HLA-DRA |
| Thyroid_lineage_extra | 4 | IYD, THADA, MET, KLK10 |

**8-gene panel** (subset of TDS_core): SLC5A5, TPO, TG, TSHR, PAX8, NKX2-1, FOXE1, DIO1

**Source file:** `project/metadata/tierA67_genes.txt`

---

## Suppl Table 3 — GSE286332 top DEGs (top 200 up + 100 dn)

Columns: `gene | log2FC | log2FC_SE | stat | pvalue | padj | rank_up | rank_dn | category | note`

Top 25 up:
| gene | log2FC | padj | category |
|---|---|---|---|
| IGHV3-66 | 7.25 | 1.4e-35 | B-cell receptor |
| BLK | 5.97 | 6.9e-33 | B-cell signaling |
| IGHV3-13 | 7.10 | 3.2e-32 | B-cell receptor |
| IGHV3-16 | 7.86 | 8.5e-32 | B-cell receptor |
| ... | ... | ... | ... |
| HLA-DOB | 5.43 | 7.3e-24 | HLA-II |
| EOMES | 4.85 | 9.4e-25 | T-cell effector |

(full table 200+100 in supplementary `.xlsx`)

**Source files:**
- `project/results/p3_gse286332/deg_ptcht_vs_ptc.tsv` (29,672 genes, TSV)
- `project/results/p3_gse286332/SuppTable_S3_GSE286332_DEGs.xlsx` (built 2026-05-04, single-sheet, sorted by padj asc)

---

## Suppl Table 4 — Pan-Asian HLA per-allele full forest

| locus | allele | Chu GD freq | Chu ctrl freq | Chu OR | Korean PTC freq | OR Kr vs ctrl | p Kr vs ctrl | Pooled OR | I²% |
|---|---|---|---|---|---|---|---|---|---|
| A | A*02:07 | 9.7% | 4.9% | 2.10 [1.70, 2.59] | 8.1% | 1.72 [1.20, 2.46] | 0.002 | 1.99 [1.66, 2.37] | 0% |
| B | B*46:01 | 14.1% | 6.5% | 2.38 [1.99, 2.86] | 10.3% | 1.65 [1.22, 2.24] | 1e-3 | 2.02 [1.41, 2.89] | 76% |
| C | C*01:02 | 18.4% | 10.9% | 1.83 [1.57, 2.12] | 24.1% | 2.61 [2.05, 3.32] | <1e-30 | 2.16 [1.53, 3.06] | 85% |
| DPA1 | DPA1*02:02 | 59.5% | 44.8% | 1.90 [1.70, 2.12] | NA (locus n/a) | NA | NA | NA | NA |
| **DPB1** | **DPB1*05:01** | **44.0%** | **31.3%** | **1.90 [1.69, 2.14]** | **53.2%** | **2.50 [2.07, 3.02]** | **4e-26** | **2.16 [1.65, 2.83]** | **85%** |
| DQA1 | DQA1*02:01 | 7.0% | 15.2% | 0.43 [0.36, 0.51] | NA (locus n/a) | NA | NA | NA | NA |
| DQB1 | DQB1*02:01 | 10.9% | 17.8% | 0.57 [0.49, 0.66] | 0.0% | 0.003 [0.0, 0.04] | 3e-5 | 0.57 [0.49, 0.66] | 0% |
| DRB1 | DRB1*07:01 | 7.1% | 15.3% | 0.43 [0.36, 0.51] | 11.4% | 0.72 [0.55, 0.93] | 9e-3 | 0.55 [0.33, 0.91] | 91% |

(Korean sub-cohort breakdown — K2/Lee/GSE286332-PTC — in same table)

**Source files:**
- `project/results/p2_pillar1_forest/forest_meta_results.tsv`
- `project/results/p2_pillar1_forest/random_effects_pooled.tsv`
- `project/results/p2_pillar1_forest/korean_PTC_pool_per_subcohort.tsv`

---

## Suppl Table 5 — Mediation analysis full (P_DM1 ~ HLA-II + 8-gene + immune)

| Mediator | a (treat→med) | b (med→outcome\|treat) | c_total | c_direct | indirect | %mediated | Boot 95% CI | p_emp |
|---|---|---|---|---|---|---|---|---|
| **HLA-II** | +1.666 | −0.110 | −0.131 | +0.052 | **−0.183** | **140%** | [−0.313, −0.032] | **0.023** |
| g8_RAI | −1.003 | +0.082 | −0.131 | −0.048 | −0.082 | 63% | [−0.152, −0.030] | 0.002 |
| immune | +1.615 | −0.071 | −0.131 | −0.016 | −0.115 | 88% | [−0.264, +0.039] | 0.120 |
| HLA-I | +1.474 | −0.077 | −0.131 | −0.017 | −0.114 | 87% | [−0.203, −0.006] | 0.042 |

OLS decomposition (z-standardized):
| Predictor | β | SE | t | p |
|---|---|---|---|---|
| const | 0.117 | 0.013 | 9.08 | <0.001 |
| HLA-II | −0.130 | 0.057 | −2.28 | 0.039 |
| g8_RAI | +0.054 | 0.027 | +2.00 | 0.065 |
| immune | +0.097 | 0.058 | +1.65 | 0.120 |

R² = 0.756 (R²_adj = 0.704), F-stat = 14.50, p = 1.4e-4, n = 18

**Source file:** `project/results/d3p5_pdm1_gradient/mediation_results.json`

---

## Suppl Table 6 — TCGA Hashimoto-like × DM cluster (4 thresholds)

| Method | Hashi+ n / 500 | DM1 hashi+% | DM2 hashi+% | OR | Fisher p |
|---|---|---|---|---|---|
| GMM 2-component | 90 (18.0%) | 5.7% | 22.8% | 0.205 | 2.1e-6 |
| Otsu threshold | 98 (19.6%) | 7.1% | 24.4% | 0.238 | 4.5e-6 |
| Top 10% | 50 (10.0%) | 4.3% | 12.2% | 0.322 | 7.4e-3 |
| Top 20% | 100 (20.0%) | 7.9% | 24.7% | 0.260 | 1e-5 |
| **Top 30%** | **150 (30.0%)** | **10.7%** | **37.5%** | **0.20** | **6.4e-10** |
| Resid Otsu (Stromal+immune residualized) | 219 (43.8%) | 23.6% | 51.7% | 0.289 | 8.0e-9 |

HLA-II Cohen d residualization:
| Subset | HLA-II d (DM1 vs DM2) |
|---|---|
| Full TCGA | −1.41 |
| Excluding Hashimoto+ (Otsu) | −1.60 |
| Within Hashimoto+ only | +0.14 (NS) |

**Source file:** `project/results/d4p2_tcga_hashimoto_signature/D4P2_summary.json`

---

## Suppl Table 7 — Korean GSE213647 (n=632) + GSE286332 sub-B replication

| Cohort | n | Hashimoto-like GMM% | Hashimoto-like Otsu% | Sub-B-like GMM% | Sub-B-like Otsu% |
|---|---|---|---|---|---|
| TCGA-THCA | 500 | 18.0% | 19.6% | NA | NA |
| **GSE213647 (Korean Lee 2024)** | **632** | **22.8%** | **28.2%** | **47.2%** | **52.5%** |
| GSE286332 PTC+HT only | 9 | 100% | 100% | NA (DM2 100%) | NA |

**Source files:**
- `project/results/d8b_korean_replication/D8B_summary.json`
- `project/results/d8c_dm1_subB_x_K2_NBNR/D8C_summary.json`

---

## Suppl Table 8 — DM1 sub-A vs sub-B mutation × signature

| Group | n | BRAF+ | RAS+ | mut-neg | Hashi+ Otsu% | 8-gene RAI median |
|---|---|---|---|---|---|---|
| **sub-A** | 84 | 1 (1%) | 51 (61%) | 32 (38%) | 3 (3.6%) | +0.03 |
| **sub-B** | 56 | 1 (2%) | 2 (4%) | 53 (94%) | 7 (12.5%) | +0.25 |

Top 15 sub-B vs sub-A DEGs (full 8,935 sig at padj<0.05 from 51,711 tested):
populated from `project/results/d6p7_dm1_subcluster/dm1_subBvA_deg.tsv` (TSV) and packaged as Sheet S8d in the multi-sheet XLSX.

**Source files:**
- `project/results/d6p7_dm1_subcluster/subcluster_score_profile.tsv` → S8a
- `project/results/d6p7_dm1_subcluster/subcluster_clinical.tsv` → S8b
- (S8c hardcoded from SUPP T8 reconciliation; mutation × Hashimoto cross-tab)
- `project/results/d6p7_dm1_subcluster/dm1_subBvA_deg.tsv` → S8d
- `project/results/d6p7_dm1_subcluster/SuppTable_S8_DM1_subcluster.xlsx` (built 2026-05-04, 4 sheets)

---

## Notes for v8 manuscript

- All numerical values cross-verified with `results/*/summary.json` files
- Wilson 95% CI used for all proportions
- DerSimonian-Laird random-effects for forest meta
- Cohen d = pooled SD method
- BH-FDR for multiple testing throughout
- All TSV files committable with submission (under `project/submission_data/`)
