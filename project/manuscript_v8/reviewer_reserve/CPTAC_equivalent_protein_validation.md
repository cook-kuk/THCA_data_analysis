# Paper 1 reviewer-reserve — multi-modality corroboration of 8-gene panel

**Status:** reserve (not in main manuscript). Use only for reviewer Q on RNA-only validity.
**Created:** 2026-05-08 (marathon mode); methylation pillar pre-existed (audit_2026_04_30 round5).
**Modalities corroborated:** RNA (main) + protein + phospho + DNA methylation. Metabolomics (Wang MTBLS3339) and 2nd proteogenomic (Cell Rep Med 2026 advanced DTC) are blocked from programmatic access — see Limits.

| Pillar | Cohort | n | Status | File pointer |
|---|---|---|---|---|
| RNA | TCGA-THCA + Korean K2 + GSE286332 | 504+260+18 | main manuscript | — |
| **Protein** | Wang 2024 (PTC) + Mun 2025 (full dediff) | 102+336 tumors | **added 2026-05-08** | `project/results/proteogenomic_v1/` |
| **Phospho** | Wang 2024 + Mun 2025 (S1D 47k sites) | 217 channels | **added 2026-05-08** | `project/results/proteogenomic_v1/paper3_mun2025_dediff_layer/phospho_*` |
| **Methylation (HM450)** | TCGA-THCA | 503 with Beta | pre-existing (r5_2) | `project/results/audit_2026_04_30/round5/r5_2_*` |
| Metabolomics | Wang MTBLS3339 | 102 PTC | BLOCKED — token req'd | — |
| 2nd proteogenomic | Cell Rep Med 2026 advanced DTC | 113 | **added 2026-05-08** (PDF table extraction) | `project/results/proteogenomic_v1/processed/cellrepmed2026_*` |

**CPTAC pan-cancer caveat:** the formal CPTAC consortium does not include thyroid; the two proteogenomic studies above are the closest equivalent. Mun 2025 has the larger and more directly dediff-relevant sample frame.

## What the panel looks like at protein level

### Wang 2024 (PTC tumor vs paired normal)

8-gene panel: SLC5A5, TPO, TG, TSHR, PAX8, NKX2-1, FOXE1, DIO1.

| Gene | protein log2FC (T vs N) | FDR | mRNA-protein gene-wise corr |
|---|---|---|---|
| TPO   | -3.18 | 1.2e-06 | 0.37 |
| TG    | -2.68 | 3.0e-05 | 0.33 |
| TSHR  | -0.31 | 8.6e-03 | 0.34 |
| FOXE1 | -0.23 | 2.0e-04 | 0.25 |
| SLC5A5, PAX8, NKX2-1, DIO1 | not in filtered protein DEG list | — | not in filtered list |

Source: Fig S2d (protein DEG, Wilcox paired) and Fig S6b (gene-wise correlation), Source Data MOESM7.

Direction at protein level matches the RNA-based panel (loss of differentiation in tumor). 4/8 reach FDR significance in this filtered table; the others may be below detection or below filter — full protein matrix locked behind controlled-access supplementary files.

### Mun 2025 (PTC → PDTC → ATC dediff axis)

7/8 panel genes detected at protein level (NKX2-1 missing). All 7 are below in ATC vs PTC; **5 of 7 show monotonic dediff gradient** PTC > PDTC > ATC.

| Gene | mean PTC | mean PDTC | mean ATC | Spearman r (axis) | p |
|---|---|---|---|---|---|
| TG    | 1.29 | 0.52 | 0.48 | -0.74 | 7e-59 |
| TSHR  | 1.29 | 1.02 | 0.52 | -0.71 | 2e-52 |
| FOXE1 | 1.30 | 1.09 | 0.75 | -0.66 | 2e-43 |
| TPO   | 1.05 | 0.91 | 0.66 | -0.61 | 4e-36 |
| PAX8  | 1.07 | 0.98 | 0.76 | -0.52 | 2e-24 |
| SLC5A5 | 0.65 | 0.76 | 0.41 | -0.21 | 1.5e-04 |
| DIO1  | 0.07 | 0.24 | 0.06 | -0.01 | 0.87 |
| NKX2-1 | not detected at protein | — | — | — | — |

Source: Mun MOESM3 Tables S1A/S1C; analysis at `project/results/proteogenomic_v1/paper3_mun2025_dediff_layer/eight_gene_dediff_trend.tsv`.

NKX2-1 IS detected at phosphoprotein level (Mun S1D), so it is expressed but escaped peptide identification at protein level — not "not present in tumors."

### TCGA-THCA HM450 methylation (DM1 hyper vs DM2 hypo)

503 TCGA-THCA samples with Illumina HM450 Beta values, stratified by DM1 vs DM2 cluster (the same RNA-defined clusters in the manuscript). 7 of 8 panel genes are hypermethylated in DM1 (the differentiation-loss cluster), 6 reach Bonferroni FDR.

| Gene | n_DM1 | n_DM2 | mean β DM1 | mean β DM2 | Δβ | Cohen's d | MW p |
|---|---|---|---|---|---|---|---|
| TPO    | 90 | 54 | 0.83 | 0.41 | **+0.42** | **+2.30** | 1.9e-18 |
| DIO1   | 90 | 54 | 0.50 | 0.27 | +0.23 | +1.24 | 6.5e-11 |
| TSHR   | 90 | 54 | 0.20 | 0.08 | +0.12 | +1.20 | 9.8e-12 |
| PAX8   | 90 | 54 | 0.09 | 0.05 | +0.045 | +0.97 | 4.5e-08 |
| TG     | 90 | 54 | 0.65 | 0.52 | +0.13 | +0.86 | 2.3e-06 |
| FOXE1  | 90 | 54 | (data on file) | — | +0.042 | +0.83 | 1.0e-05 |
| NKX2-1 | 90 | 54 | 0.08 | 0.03 | +0.053 | +0.63 | 8.9e-07 |
| SLC5A5 | 90 | 54 | 0.59 | 0.56 | +0.026 | +0.23 | 0.42 (NS) |

Source: `project/results/audit_2026_04_30/round5/r5_2_methylation_DM.json` and `r5_2_per_gene_methylation_DM.tsv` (already on disk; pre-existing analysis from audit round 5).

**Triangulation:** RNA panel scores low → protein abundance low (Mun) → phosphosite detection sparse (Mun) → CpG hypermethylated (TCGA HM450). All four orthogonal layers point to *biological* loss of differentiation in DM1, not to RNA technical artifact. SLC5A5 is the gene with weakest methylation differential — consistent with its weakest protein-level dediff Spearman in Mun (-0.21); honest to disclose.

### Forward reference: pan-cancer DM1 axis (Paper 11 in preparation)

The DM1 axis defined here for THCA generalizes across cancer lineages (Paper 11, target Nat Commun). Lineage-portable DM1 scoring (cancer-specific TF panel + 21-gene Module 4 architecture: lineage_TF + JAK/STAT + SFK + Epigenetic + Metabolic) reproduces in 12 cancer types × 6,216 TCGA samples, with **independent prognostic effect in 4 cancers** (skin melanoma HR=0.78 p=4e-7 protective, lower-grade glioma HR=1.19 p=2e-3 risk, pancreatic adeno HR=1.18 p=5e-3 risk, lung adeno HR=1.11 p=0.05 risk). DepMap CRISPR essentiality across 1,141 cell lines identifies **MYC** (d=−0.50, p=1e-11) and **NAMPT** (d=−0.45, p=1.5e-9) as druggable dependencies enriched in DM1-high lineages. See `project/results/paper11_pancancer/PAPER11_OUTLINE.md`.

### Cell Rep Med 2026 advanced DTC (Zhang et al., 113 patients, FUSCC cohort)

Independent proteogenomic study published March 2026, n=113 advanced DTC. Authors performed proteomic consensus clustering and identified three subtypes (CC1 canonical, CC2 stromal, CC3 immunogenic). Their differentiation read-out is a **TDS (Thyroid Differentiation Score)** computed by AddModuleScore over 20 genes (Table S5 of their supplements):

> DIO1, DIO2, DUOX1, DUOX2, **FOXE1**, GLIS3, **NKX2-1**, **PAX8**, SLC26A4, **SLC5A5**, SLC5A8, **TG**, THRA, THRB, **TPO**, **TSHR**, GLIPR1, SLIT2, TIMP3, ITGB5

**8/8 of our panel is contained in TDS** (bold above). Our 8-gene panel is therefore a strict subset of an independently-derived 20-gene proteomic differentiation score — corroborating that the panel genes are the right targets for capturing this axis.

CC subtype × RAI sensitivity (extracted from Cell Rep Med 2026 Table S1):

| CC subtype | n | RAI Avid | RAI Refractory | Refractory % |
|---|---|---|---|---|
| CC1 (canonical, high TDS)   | 43 | 34 |  9 | 21 % |
| CC2 (stromal)               | 40 | 17 | 23 | 58 % |
| CC3 (immunogenic, low TDS)  | 30 |  4 | 26 | **87 %** |

CC subtype × histology is mixed (cPTC dominates all three CC subtypes; FTC and fvPTC distribute across all CC), confirming CC is a **molecular** subtype not a histology proxy.

**Cross-cohort axis alignment:**

| Mun 2025 PTC→PDTC→ATC | CRM 2026 advanced DTC | Our DM1/DM2 (TCGA) |
|---|---|---|
| ATC = +myeloid / -thyroid_diff (d ±1.5–1.9) | CC3 = high immune / low TDS / 87 % RAI-refractory | DM1 = low panel / hypermethylated TDS genes |
| PTC = thyroid_diff intact | CC1 = high TDS / 79 % RAI-Avid | DM2 = panel-positive |

Three independent cohorts (TCGA + Mun + CRM) converge on the same dediff axis with consistent direction at every modality (RNA, protein, phospho, methylation, clinical RAI response).

## Reviewer-Q-ready talking points

1. **Panel is RNA-defined but biologically grounded at protein level** — TG/TSHR/FOXE1/TPO/PAX8 each loss-of-differentiation Spearman r < -0.5 (FDR ≪ 1e-20) along PTC→PDTC→ATC at protein level in Mun 2025 (n=336 tumors).
2. **Modest mRNA-protein correlation** (~0.25–0.37 in Wang) is the well-known biological reality and explains why a multi-gene mean (panel) is more robust than any single gene at either RNA or protein level.
3. **DIO1's flat protein trajectory** in Mun is honest to flag — could indicate post-translational regulation or detection coverage; the panel does not depend on any single gene.
4. **NKX2-1 protein-level non-detection** is a coverage caveat (still detected at phospho), not a biological null.

## Limits

- Wang full protein abundance matrix not publicly accessible (MOESM3-6 access-blocked from Springer); only Source Data MOESM7 + PDF.
- Mun WES/RNA in GSA-Human under controlled access; we only used the public protein/phospho matrices (S1C/S1D).
- Wang **metabolomics** (MetaboLights MTBLS3339) — programmatic fetch blocked, public API requires token in 2026. User can browser-download `s_*.txt` + `m_*.tsv` if metabolomic pillar is wanted.
- Cell Rep Med 2026 advanced DTC (113 patients, CC1/CC2/CC3) — supplements (mmc1.pdf, mmc2.pdf) acquired via user browser-download 2026-05-08. PDF table extraction gives Table S1 (clinical + CC subtype × 113), Table S5 (TDS gene set). xlsx-form per-sample protein abundance matrix is NOT included in the PDF supplements — would require separate Data S1-S7 download for full proteomic-level integration.
- TCGA miRNA-seq vs 8-gene targeting — feasible from GDC but deferred (sprint-class new analysis under marathon mode).
- Lu 2023 ATC scRNA per-cell-type module decomposition — feasible from raw tar at `/data/thca/v17_lu2023_GSE193581/` but deferred (sprint-class).
- Two proteogenomic cohorts are Chinese; methylation is TCGA (mixed); transferability to Korean K2 is by analogy not direct reanalysis (K2 has no public protein/phospho).
- This memo is reviewer-reserve only — does NOT belong in main text unless explicitly asked.
