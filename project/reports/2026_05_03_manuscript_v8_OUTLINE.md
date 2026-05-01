# Manuscript v8 Outline — Phase 0 Cancer Paper (5-Pillar)

**Date:** 2026-05-03 (marathon mode prep, week 0)
**Target venue:** Cell Rep Med (IF 14) primary / JCI Insight (IF 8) backup
**Length target:** Cell Rep Med ~5,000 words (excl. methods); 7 main figures + 6 suppl figures
**Estimated draft timeline:** W1-W2 outline → W3-W5 first draft → W6 revision → W7-8 submission

---

## Working Title (3 candidates, 본인 voice 결정)

1. **"A transcriptional differentiation axis stratifies papillary thyroid cancer by autoimmune-driven dedifferentiation, independent of canonical driver mutations"**
2. **"DM1/DM2 — an eight-gene transcriptional axis reveals autoimmune-PTC overlap as a distinct molecular subtype"**
3. **"Beyond BRAF/RAS/TERT: a Pan-Asian transcriptional and HLA framework for thyroid cancer subtyping"**

---

## Abstract (250 words target, structured)

**Background:** Differentiated thyroid carcinoma (DTC) shows 15-35% recurrence yet anatomic risk stratification offers no transcriptional axis to predict molecular trajectory.

**Methods:** TCGA-THCA n=500 dynamic-risk-stratified into DM1/DM2 by 8-gene panel (SLC5A5/TPO/TG/TSHR/PAX8/NKX2-1/FOXE1/DIO1, anchored in canonical RAI biology). Integrated with K2 (PRJEB11591 n=260) + Lee 2024 (GSE213647 n=632) + GSE286332 (n=18) Korean RNA-seq cohorts and Chu et al. 2018 Han Chinese GD HLA summary statistics (n=2,958).

**Results:** (1) **Pan-Asian HLA cohort n=874** with DPB1\*05:01 53.2% replicates Korean Graves' risk continuum (vs Chu Han Chinese GD 44%, ctrl 31%). (2) **GSE286332 PTC+HT** (n=18) shows 10,380 DEGs, 8-gene RAI Cohen d=−1.60, HLA-II d=+3.65, GSEA Hallmark IFN-γ FDR=2e-4 + KEGG Type I diabetes FDR=0. (3) **Driver mRNA neutrality**: BRAF transcript Cohen d=−0.04 vs WT (n=273 vs 182, p=0.57); driver_anchor cluster ARI=0. (4) **Pan-genome cluster robustness**: TIERA67 ARI=0.90 ≈ unrestricted top-5000 MAD ARI=0.92. (5) **★ Autoimmune-PTC mechanism**: TCGA Hashimoto-like signature DM2-enriched (OR up to 5×, p=6e-10); HLA-II 140% mediation of PTC+HT → P(DM1) ↓ pathway; antigen-driven BCR clonal expansion (TLS d=+1.96, AICDA up); DM1 sub-B = 96% mutation-negative NBNR cluster.

**Conclusion:** A transcriptional axis (DM1/DM2) anchored in canonical thyroid differentiation biology stratifies PTC orthogonally to driver mutation status, with autoimmune-overlap defining a mechanistically distinct DM2 sub-population. Independent Korean replication (GSE213647 n=632, Hashimoto-like 22.8%; sub-B-like 47-53%) confirms generalizability.

**Keywords:** thyroid cancer, papillary thyroid carcinoma, differentiation axis, HLA, autoimmune, Hashimoto's thyroiditis, transcriptomics, BRAF, TERT, NBNR

---

## Introduction (3-4 paragraphs, ~600 words)

### ¶1 Hook (Discussion §3.1 framing — 본인 voice)
Thyroid cancer epidemiology + ATA 2015 risk stratification gap (3-5% to 50-75% recurrence range, no molecular axis) + Landa 2016 PDTC/ATC "molecular endpoints" framing.
- 본인 직접: "no mechanistic compass" vs "without molecular guidance" 결정

### ¶2 Existing molecular landscape
- TCGA Cancer Network 2014 BRAF/RAS dichotomy (BRS score)
- Yoo 2016 16-gene transcriptional differentiation panel
- Krishnamoorthy / Landa 2016 PDTC/ATC genomic hallmarks
- Pan/Ge 2025 proteogenomic ATC subtypes (RRP9, C5AR1)
- **Gap:** None of these address autoimmune-PTC overlap as a distinct mechanism

### ¶3 Korean / Asian-specific HLA architecture
- DPB1\*05:01 + B\*46:01 Asian Graves' alleles (Chu 2018)
- Korean PTC + Hashimoto's overlap clinical observation (Yu professor cohort)
- cookHLA Nat Commun pipeline (rheumatoid + T1D + Crohn's first-author work)
- **Bridge:** Cross-disease HLA imputation pipeline applied to thyroid cancer

### ¶4 Study aim + 5-Pillar overview
"Here we establish a transcriptional axis..." [본인 voice 영역, last paragraph]

---

## Methods (1,500-2,000 words, dense)

### M1 — Study design + cohort assembly
- TCGA-THCA n=500 (R1A DM1/DM2 cluster definition, 5/1 P-series)
- K2 (PRJEB11591 n=260, Yoo 2016 SNU-GMI; n=235 valid arcasHLA after QC)
- Lee 2024 (GSE213647 n=632, n=630 valid arcasHLA)
- GSE286332 (Lim 2025 Dongguk Univ n=18: 9 PTC + 9 PTC+HT)
- Chu et al. 2018 J Med Genet 55(10):685-692 (1,468 GD vs 1,490 Han Chinese ctrl, public summary stats)

### M2 — TIERA67 7-category candidate pool
- 67 genes from 7 thyroid-relevant categories (TDS_core 16, MAPK_output 10, Driver_anchor 12 [BRAF/NRAS/HRAS/KRAS/RET/NTRK1/3/ALK/PAX8/PPARG/TERT/EIF1AX], Aggressive_marker 10, Dediff_invasion 10, Immune_stromal_light 5, Thyroid_lineage_extra 4)
- 8-gene panel: subset of TDS_core based on canonical RAI biology

### M3 — Differential expression + GSEA (GSE286332)
- PyDESeq2 v0.5.4 Wald test, BH-FDR
- gseapy 1.1.13 prerank with MSigDB Hallmark / KEGG 2021 / Reactome 2022 (1,000 permutations, seed=42)
- Score = −log10(p) × sign(log2FC)

### M4 — DM1/DM2 classifier (centered profile)
- Within-sample-centered 8-gene profile
- LogisticRegression on TCGA-THCA, 5-fold CV AUC = 0.962
- Applied to GSE286332 + Lee for cross-cohort prediction

### M5 — TCGA Hashimoto-like signature transfer
- Top 150 up + 50 down DEGs (padj<0.01, |LFC|>1) Z-mean scored on TCGA n=500
- Bimodality coefficient + GMM/Otsu/quartile thresholds
- Stromal + generic immune-proxy residualization (Pillar 5 confounder check)

### M6 — Mediation analysis (Baron-Kenny + bootstrap)
- 5,000-iteration bootstrap for indirect effect
- HLA-II / 8-gene RAI / immune-proxy as candidate mediators

### M7 — BCR repertoire + clonality (gene-level proxy)
- Per-sample IGH/IGK/IGL V/J/C diversity (Shannon entropy)
- TLS 12-gene signature (Cabrita 2020)
- AICDA somatic hypermutation enzyme expression

### M8 — DM1 sub-cluster (KMeans k=2 on TIERA67)
- Re-clustered DM1 (n=140) into sub-A (n=84) + sub-B (n=56)
- DESeq2-equivalent Welch t-test + BH-FDR per gene

### M9 — Pan-Asian HLA forest meta
- arcasHLA RNA-seq imputation (4-digit, IPD-IMGT/HLA v3.44.0)
- Korean PTC pool n=874 + Chu 2018 Han Chinese GD published OR
- DerSimonian-Laird random-effects + I² + Cochran's Q
- 4-scenario sensitivity (full / excl GSE286332 / Lee only / K2 only)

### M10 — Pan-genome MAD selection
- Top 5000 MAD on TCGA-THCA × KMeans k=2
- ARI vs original DM1/DM2; hypergeometric enrichment p

### M11 — Statistics + reproducibility
- All Cohen's d (pooled SD), MW two-sided
- Wilson 95% CI for proportions
- Code: github.com/seunghocook/thyca-paper-2026 (on submission)

---

## Results (5 sections, paragraph-per-pillar)

### R1 — Pillar 1: Korean Pan-Asian HLA cohort n=874

- 3-arm forest meta (Korean PTC vs Chu ctrl vs Chu GD)
- DPB1\*05:01 53.2% Korean / 31.3% Chu ctrl / 44.0% Chu GD; **OR Kr vs ctrl=2.50, p=4e-26**
- 6 alleles direction-consistent; pooled OR table
- Korean sub-cohort heterogeneity (DPB1\*05:01 I²=0% perfect)
- Sensitivity 4-scenario stable
- → **Figure 1** (Pan-Asian forest plot) + Suppl Table 1 (per-allele full)

### R2 — Pillar 2: GSE286332 PTC vs PTC+HT molecular dissection

- 10,380 DEGs (padj<0.05): 6,004 up + 4,376 down
- Top up: IGHV/IGKV/IGLV multiple, BLK, EOMES, HLA-DOB
- GSEA: Hallmark IFN-γ Response (NES=+1.80, FDR=2e-4); KEGG Type I diabetes (NES=+1.92, FDR=0)
- 8-gene RAI Cohen d=−1.60 (p=0.008); per-gene PAX8 d=−2.32, NKX2-1 d=−1.92
- HLA-II Cohen d=+3.65 (p=4e-4)
- → **Figure 2** (multi-panel: volcano + GSEA + 8-gene boxplot + HLA-II heatmap) + Suppl Table 2 (top 200 DEGs)

### R3 — Pillar 3: Driver mRNA neutrality

- BRAF transcript V600E vs WT: Cohen d=−0.04 (p=0.57, n=273 vs 182)
- HRAS/NRAS/KRAS transcripts: all |d| < 0.4
- Driver single-feature AUC for DM: BRAF 0.602, TERT 0.578, others < 0.55
- TIERA67 univariate Cohen d ranking: 8-gene at #4-50 (top tier); driver_anchor at #52-67 (bottom)
- Driver_anchor cluster ARI=−0.007 (random)
- → **Figure 3** (driver mRNA × mutation + AUC + TIERA67 ranking)

### R4 — Pillar 4: Pan-genome cluster robustness

- 8-gene panel ARI=0.49 (modest, clinically interpretable)
- TIERA67 ARI=0.90, pan-genome top-5000 MAD ARI=0.92 (cluster definition robust)
- Driver_anchor only ARI=−0.007 (drivers cannot define DM)
- TIERA67 hypergeometric enrichment in pan-genome top 100: p=3e-4
- → **Figure 4** (ARI bar + top-N coverage ladder) + Suppl Fig 1 (full pan-genome ranking)

### R5 — ★ Pillar 5: Autoimmune-PTC mechanism layer

#### 5a. P_DM1 mediation (GSE286332 n=18)
- Spearman ρ(P_DM1, HLA-II) = −0.81 (p=5e-5)
- Linear regression R²=0.756
- Baron-Kenny: HLA-II 140% mediation (boot p=0.023), 8-gene 63% (boot p=0.002)

#### 5b. TCGA Hashimoto-like generalization
- 18-30% TCGA samples Hashimoto-like (GMM 18%, Otsu 19.6%, top30% = 30%)
- DM2-enriched: top30% threshold OR=0.20, p=6e-10 (3-5× higher rate than DM1)
- Confounder-residualized OR=0.29 (still significant)
- Korean GSE213647 replication: Hashimoto-like 22.8% (Otsu 28.2%) ≈ TCGA 18-20%

#### 5c. BCR clonal + TLS (GSE286332)
- TLS Cohen d=+1.96 (p<0.001)
- IGHV clonality d>0.5; AICDA up-regulated
- ρ(IGHV clonality, 8-gene RAI) = −0.67 (p=0.002)

#### 5d. DM1 sub-B = NBNR cluster
- KMeans k=2 on DM1 (n=140) → sub-A (n=84, 69% RAS+) + sub-B (n=56, 96% mut-neg)
- sub-B Hashimoto-like 12.5% vs sub-A 3.6%
- Korean GSE213647 sub-B-like rate 47-53%

→ **Figure 5** (4-panel: P_DM1 mediation + TCGA Hashimoto-DM2 + BCR/TLS + sub-B mutation) + Suppl Fig 2 (mediation table) + Suppl Fig 3 (cross-cohort generalization)

---

## Discussion (4 sections, ~1,500 words)

### D1 — Five pillars synthesis (¶1)
- 본인 voice: "Five orthogonal lines of evidence converge..." [original framing]
- DM1/DM2 axis as transcriptional + autoimmune complement to anatomic risk

### D2 — ATA 2015 + 8-gene panel complementarity (¶2)
- "ATA 2015 risk stratifies recurrence by anatomy; DM1/DM2 stratifies dedifferentiation by transcription"
- BRAF V600E + TERT promoter "marginal benefit" gap → transcriptional axis fills this
- 본인 voice + ATA cheatsheet phrasing

### D3 — Pan-Asian autoimmune-thyroid susceptibility continuum (¶3)
- DPB1\*05:01 53% Korean / 44% Chu GD / 31% ctrl
- cookHLA Nat Commun cross-disease application
- Korean PTC ≠ Chinese GD ≠ EUR Hashimoto: population-specific HLA architecture

### D4 — Limitations + future work (¶4)
- Korean PTC pool size imbalance vs Chu (874 vs 2,958)
- Phenotype heterogeneity (PTC vs GD distinct)
- Allele-level only (haplotype future)
- BCR-seq for definitive clonal expansion (future)
- Bundang prospective Graves' cohort access pending

---

## Figures (7 main + 6 suppl)

### Main figures
- **F1** Pan-Asian HLA forest meta (Pillar 1) ★ NEW
- **F2** GSE286332 PTC+HT multi-panel (Pillar 2) ★ NEW (built today)
- **F3** Driver mRNA neutrality (Pillar 3)
- **F4** Pan-genome ARI + 8-gene clinical-interpretability (Pillar 4)
- **F5** Autoimmune-PTC mechanism 4-panel (Pillar 5) ★ NEW (built today, paper-shaping)
- **F6** TCGA Hashimoto-like distribution + cross-cohort (Pillar 5 expansion)
- **F7** DM1 sub-A vs sub-B + K2 NBNR clinical phenotype

### Suppl figures
- **SF1** TIERA67 pan-genome rank ladder
- **SF2** Mediation Baron-Kenny detail
- **SF3** Korean GSE213647 replication detail
- **SF4** GSE286332 BCR repertoire detail
- **SF5** ATA 2015 risk tier × DM cluster cross-tab
- **SF6** Sensitivity analyses (4-scenario forest)

---

## Author contributions (placeholder)

- **Seungho Cook (first/corresponding author):** Conceptualization, methodology, data curation, formal analysis, software, writing – original draft, writing – review & editing
- **Yu professor:** Supervision, project administration, resources, writing – review & editing

---

## Acknowledgments
- ARIA project allocation
- TCGA Research Network for thyroid cancer dataset
- Yoo SK / SNU-GMI (PRJEB11591), Lee Y / Macrogen (GSE213647), Lim DW / Dongguk (GSE286332)
- cookHLA / arcasHLA tool maintainers
- Computing: Azure Korea Central burst pattern (~$5.80 total)

---

## Data + code availability

- All processed data: `project/results/` (GitHub commit on submission)
- arcasHLA results: `project/results/v17_korean/arcasHLA*/`
- Source code: github.com/seunghocook/thyca-paper-2026
- TCGA-THCA: GDC portal
- GSE213647, GSE286332: NCBI GEO
- PRJEB11591: ENA (Yoo SK 2016 SNU-GMI)

---

## v8 timeline (5/4 → 6/13, 6 weeks)

| Week | Focus | Deliverable |
|---|---|---|
| W1 (5/4-5/10) | Title + Abstract + Outline confirm | Methods §M1-M11 first pass |
| W2 (5/11-5/17) | Introduction + Methods first draft | Result R1-R5 first pass |
| W3 (5/18-5/24) | Results 5-pillar narrative | Figure F1-F5 polish |
| W4 (5/25-5/31) | Discussion D1-D4 | Suppl figures + tables |
| W5 (6/1-6/7) | Self-revision + Yu professor 1차 review | v8 → v9 |
| W6 (6/8-6/13) | Final revision + bioRxiv submission | bioRxiv DOI |

→ **6/13 bioRxiv preprint submission target.**
→ Cell Rep Med formal submission Q3 2026.
