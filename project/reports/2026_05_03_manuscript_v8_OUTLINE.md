# Manuscript v8 Outline — Phase 0 Cancer Paper (5-Pillar)

**Date:** 2026-05-03 (marathon mode prep, week 0)
**Target venue:** Cell Rep Med (IF 14) primary / JCI Insight (IF 8) backup
**Length target:** Cell Rep Med ~5,000 words (excl. methods); 7 main figures + 6 suppl figures
**Estimated draft timeline:** W1-W2 outline → W3-W5 first draft → W6 revision → W7-8 submission

---

> **Scope status (2026-05-04 update).** This v8 outline was authored 2026-05-03 as a 5-pillar mega-paper.
> Per `paper_numbering_2026_05_04` and `2026_05_04_8gene_curated_vs_denovo_final_strategy.md`, the mega-paper
> has since been carved into Paper 1 (DM1 molecular dark matter; this outline's primary scope going forward),
> Paper 2 (Hashimoto-overlap PTC; `v18_paper2_HT_isolated`), and Paper 4 (Korean GD HLA Pan-Asian; backlog).
> Sections below are now tagged `[PAPER 1]` / `[PAPER 2 — moved]` / `[PAPER 4 — moved]`. Paper 1 marathon
> (5/4–6/13 → bioRxiv 6/13) operates only on `[PAPER 1]` sections. Existing body prose for moved sections is
> retained for Paper 2 / Paper 4 future use; do not delete.

---

## Working Title (4 candidates, 본인 voice 결정)

1. ★ **T1 (primary, recommended):** "A transcriptional differentiation axis stratifies thyroid cancer orthogonally to canonical driver mutations and identifies a tumor-population TROP2 vulnerability"
2. **T2 (backup, BRAF/RAS-beyond hook):** "Beyond BRAF and RAS: a thyroid-lineage silencing axis defines a TROP2-targetable papillary thyroid cancer subtype"
3. **T3 (mechanism-led):** "Lineage transcription-factor collapse defines a dedifferentiation subtype of papillary thyroid cancer with TROP2 vulnerability"
4. **T4 (safe/bland):** "Driver-mutation-independent dedifferentiation in thyroid cancer: an RNA-defined axis with prognostic and therapeutic implications"

(Old candidates T5 "8-gene-headline" / T6 "dark-matter-headline" / T7 "spatial-led" are explicitly deprecated — see `2026_05_04_8gene_curated_vs_denovo_final_strategy.md` §10.)

---

## Abstract (250 words target, structured)

**Background:** Differentiated thyroid carcinoma (DTC) shows 15-35% recurrence yet anatomic risk stratification offers no transcriptional axis to predict molecular trajectory.

> **Note (2026-05-04 [PAPER 1] re-scope).** The Methods / Results / Conclusion below are the Paper 1-scoped
> abstract skeleton. The pre-2026-05-04 5-pillar version (HLA forest + GSE286332 + autoimmune-PTC mechanism)
> has been moved to Paper 2 / Paper 4 abstracts; see `2026_05_04_8gene_curated_vs_denovo_final_strategy.md` §12.2.
> The pre-2026-05-04 5-pillar Abstract is preserved below the keywords block as
> `[PAPER 2/4 — moved] legacy abstract` for cross-paper reference.

**Methods.** TCGA-THCA (n=500) was clustered on a 67-gene thyroid-relevant candidate pool (TIERA67) with KMeans k=2, yielding DM1/DM2 labels. Driver mutation orthogonality was tested by univariate Cohen's d, single-feature AUC, and driver-only k-means. Cluster robustness was tested by pan-genome top-N MAD selection (top-100 / -200 / -1000 / -5000). A compact eight-gene (SLC5A5, TPO, TG, TSHR, PAX8, NKX2-1, FOXE1, DIO1) within-sample-centered profile classifier was trained for DM1_like scoring. Mechanism (transcription-factor activity, Hallmark GSEA, TF network) and TROP2 / TACSTD2 differential expression were assessed in the n=261 expression-complete subset. Outcome was modeled by Kaplan–Meier and multivariate Cox (BRAF / RAS / age / stage adjusted) on n=461. Spatial validation used 28 Visium ST slides across three external cohorts (Moran's I, bivariate Moran with a zero-overlap thyroid-lineage gene set, autoimmune-without-cancer negative control). Methods M2 panel-selection rationale clarifies that driver genes were retained in the candidate pool but ranked at the bottom of univariate Cohen's d (not "excluded by design").

**Results.** (1) BRAF / RAS / TERT mutation status does not propagate to discriminative driver-transcript signal (\|d\| < 0.4, AUC < 0.7); driver-only clustering yields ARI = −0.007. (2) DM1/DM2 cluster definition is robust at the unrestricted transcriptome level (pan-genome top-5000 MAD ARI = 0.918 ≈ TIERA67 ARI = 0.903). (3) The eight-gene compact readout retains ≈99% of the canonical 16-gene differentiation panel's discriminative power (ΔAUC = 0.013, NS) with drop-one-out r ∈ [−0.987, −0.981]. (4) Lineage-TF collapse (FOXE1, NKX2-1, PAX8 down) co-occurs with DNMT1, STAT3, FOSL1 up-regulation, IL-6 / JAK / STAT3 Hallmark enrichment (NES = +8.0), and OXPHOS down (NES = −5.5). (5) DM1_like supports an independent prognostic association in multivariate Cox (PFI HR = 2.04 [1.15–3.61], p = 0.015; DFI HR = 1.41, p = 0.025). (6) TROP2 / TACSTD2 is a tumor-population vulnerability (Δmean = +4.33, FDR = 1.2 × 10⁻³²; tumor 1.4× normal). (7) Spatial validation across 28 ST slides supports the axis: bulk r = −0.885 with a zero-overlap lineage gene set (n = 561); Moran's I 26/28 significant; sample-mean ρ ≤ −0.85 in 12 external slides; HT / GD inflammation alone does not raise DM1_like above control.

**Conclusion.** A transcriptional differentiation axis, robust at the unrestricted transcriptome level and orthogonal to canonical driver mutations, stratifies papillary thyroid cancer along a thyroid-lineage / RAI-uptake silencing dimension. The axis is operationalisable as a compact eight-gene readout, mechanistically anchored in lineage-TF collapse with downstream DNMT1 / STAT3 activation, supports an independent prognostic association after BRAF / RAS / age / stage adjustment, and identifies a tumor-population TROP2 vulnerability that motivates evaluation of TROP2-directed antibody-drug-conjugate strategies.

**Keywords:** thyroid cancer, papillary thyroid carcinoma, differentiation axis, dedifferentiation, BRAF, RAS, TROP2, TACSTD2, lineage transcription factors, spatial transcriptomics

---

### `[PAPER 2/4 — moved]` legacy abstract (pre-2026-05-04 5-pillar; retained for Paper 2 / Paper 4 cross-reference, NOT for Paper 1 submission)

**Methods (legacy):** TCGA-THCA n=500 dynamic-risk-stratified into DM1/DM2 by 8-gene panel (SLC5A5/TPO/TG/TSHR/PAX8/NKX2-1/FOXE1/DIO1, anchored in canonical RAI biology). Integrated with K2 (PRJEB11591 n=260) + Lee 2024 (GSE213647 n=632) + GSE286332 (n=18) Korean RNA-seq cohorts and Chu et al. 2018 Han Chinese GD HLA summary statistics (n=2,958).

**Results (legacy):** (1) **Pan-Asian HLA cohort n=874** with DPB1\*05:01 53.2% replicates Korean Graves' risk continuum (vs Chu Han Chinese GD 44%, ctrl 31%). (2) **GSE286332 PTC+HT** (n=18) shows 10,380 DEGs, 8-gene RAI Cohen d=−1.60, HLA-II d=+3.65, GSEA Hallmark IFN-γ FDR=2e-4 + KEGG Type I diabetes FDR=0. (3) **Driver mRNA neutrality**: BRAF transcript Cohen d=−0.044 vs WT (n=273 vs 182, p=0.567); driver_anchor cluster ARI=0. (4) **Pan-genome cluster robustness**: TIERA67 ARI=0.903 ≈ unrestricted top-5000 MAD ARI=0.918. (5) **Autoimmune-PTC mechanism**: TCGA Hashimoto-like signature DM2-enriched (OR up to 5×, p=6.4×10⁻¹⁰); HLA-II 140% mediation of PTC+HT → P(DM1) ↓ pathway; antigen-driven BCR clonal expansion (TLS d=+1.96, AICDA up); DM1 sub-B = 53/56 (94.6%) mutation-negative NBNR cluster.

**Conclusion (legacy):** A transcriptional axis (DM1/DM2) anchored in canonical thyroid differentiation biology stratifies PTC orthogonally to driver mutation status, with autoimmune-overlap defining a mechanistically distinct DM2 sub-population. Independent Korean replication (GSE213647 n=632, Hashimoto-like 22.8%; sub-B-like 47-53%) confirms generalizability.

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

## Results (Paper 1 = R3 / R4 + new R6 / R7 / R8 / R9 scaffold; R1 / R2 / R5 moved)

> **[2026-05-04 re-scope]** Per `2026_05_04_8gene_curated_vs_denovo_final_strategy.md` §5, Paper 1 results
> are: driver landscape orthogonality (R3), axis robustness + compact readout (R4), mechanism (new R6),
> outcome + TROP2 (new R7), spatial validation (new R8), honest negatives (new R9). R1 / R2 / R5 below are
> moved to Paper 2 / Paper 4; their body prose is retained in `2026_05_03_results_R1_R5_prose.md` and
> `2026_05_03_figure_captions_all.md` for those papers' future use, NOT for Paper 1 marathon submission.

### R1 [PAPER 4 — moved] — Pillar 1: Korean Pan-Asian HLA cohort n=874

- 3-arm forest meta (Korean PTC vs Chu ctrl vs Chu GD)
- DPB1\*05:01 53.2% Korean / 31.3% Chu ctrl / 44.0% Chu GD; **OR Kr vs ctrl=2.50, p=4e-26**
- 6 alleles direction-consistent; pooled OR table
- Korean sub-cohort heterogeneity (DPB1\*05:01 I²=0% perfect)
- Sensitivity 4-scenario stable
- → **Figure 1** (Pan-Asian forest plot) + Suppl Table 1 (per-allele full)

### R2 [PAPER 2 — moved] — Pillar 2: GSE286332 PTC vs PTC+HT molecular dissection

- 10,380 DEGs (padj<0.05): 6,004 up + 4,376 down
- Top up: IGHV/IGKV/IGLV multiple, BLK, EOMES, HLA-DOB
- GSEA: Hallmark IFN-γ Response (NES=+1.80, FDR=2e-4); KEGG Type I diabetes (NES=+1.92, FDR=0)
- 8-gene RAI Cohen d=−1.60 (p=0.008); per-gene PAX8 d=−2.32, NKX2-1 d=−1.92
- HLA-II Cohen d=+3.65 (p=4e-4)
- → **Figure 2** (multi-panel: volcano + GSEA + 8-gene boxplot + HLA-II heatmap) + Suppl Table 2 (top 200 DEGs)

### R3 [PAPER 1] — Driver landscape orthogonality (BRAF / RAS / TERT mutation × transcript)
> Re-labeled 2026-05-04 (was: "Pillar 3: Driver mRNA neutrality"). Strategy memo §5 Layer 1: drivers cannot
> define the axis (Driver_anchor ARI = −0.007), establishing orthogonality of DM1/DM2 to canonical mutation
> status. Body prose in `2026_05_03_results_R1_R5_prose.md` §R3 stays valid; only the heading framing is
> updated to "orthogonality" rather than "neutrality" per strategy memo §11.

- BRAF transcript V600E vs WT: Cohen d=−0.044 (p=0.567, n=273 vs 182)
- HRAS/NRAS/KRAS transcripts: all |d| < 0.4
- Driver single-feature AUC for DM: BRAF 0.602, TERT 0.578, others < 0.55
- TIERA67 univariate Cohen d ranking: 8-gene at #4-50 (top tier); driver_anchor at #52-67 (bottom)
- Driver_anchor cluster ARI=−0.007 (random)
- → **Figure 3** (driver mRNA × mutation + AUC + TIERA67 ranking)

### R4 [PAPER 1] — Axis robustness — pan-genome cluster + compact eight-gene readout
> Re-labeled 2026-05-04 (was: "Pillar 4: Pan-genome cluster robustness"). Strategy memo §5 Layers 2–3:
> the DM1/DM2 cluster is robust at the unrestricted transcriptome level (pan-genome top-5000 MAD
> ARI = 0.918 ≈ TIERA67 ARI = 0.903); the eight-gene panel is the *compact RAI-lineage readout* of that
> axis, not a discovery panel (8/8 ⊂ Yoo 2016 TDS-16; ΔAUC vs k=16 = 0.013, NS; drop-one-out r ∈ [−0.987,
> −0.981]). Body prose in `2026_05_03_results_R1_R5_prose.md` §R4 stays valid; framing emphasises that 8-gene
> is a clinical-interpretability trade-off rather than maximum-power optimisation.

- 8-gene panel ARI=0.49 (modest, clinically interpretable)
- TIERA67 ARI=0.903, pan-genome top-5000 MAD ARI=0.918 (cluster definition robust)
- Driver_anchor only ARI=−0.007 (drivers cannot define DM)
- TIERA67 hypergeometric enrichment in pan-genome top 100: p=3e-4
- → **Figure 4** (ARI bar + top-N coverage ladder) + Suppl Fig 1 (full pan-genome ranking)

### R5 [PAPER 2 — moved] — ★ Pillar 5: Autoimmune-PTC mechanism layer

#### 5a. P_DM1 mediation (GSE286332 n=18)
- Spearman ρ(P_DM1, HLA-II) = −0.81 (p=5e-5)
- Linear regression R²=0.756
- Baron-Kenny: HLA-II 140% mediation (boot p=0.023), 8-gene 63% (boot p=0.002)

#### 5b. TCGA Hashimoto-like generalization
- 18-30% TCGA samples Hashimoto-like (GMM 18%, Otsu 19.6%, top30% = 30%)
- DM2-enriched: top30% threshold OR=0.20, p=6.4×10⁻¹⁰ (3-5× higher rate than DM1)
- Confounder-residualized OR=0.29 (still significant)
- Korean GSE213647 replication: Hashimoto-like 22.8% (Otsu 28.2%) ≈ TCGA 18-20%

#### 5c. BCR clonal + TLS (GSE286332)
- TLS Cohen d=+1.96 (p<0.001)
- IGHV clonality d>0.5; AICDA up-regulated
- ρ(IGHV clonality, 8-gene RAI) = −0.67 (p=0.002)

#### 5d. DM1 sub-B = NBNR cluster
- KMeans k=2 on DM1 (n=140) → sub-A (n=84; 51/74 mut-tested RAS+ = 69%, 51/84 = 61% of total) + sub-B (n=56; 53/56 mut-neg = 94.6%)
- sub-B Hashimoto-like 12.5% vs sub-A 3.6%
- Korean GSE213647 sub-B-like rate 47-53%

→ **Figure 5** (4-panel: P_DM1 mediation + TCGA Hashimoto-DM2 + BCR/TLS + sub-B mutation) + Suppl Fig 2 (mediation table) + Suppl Fig 3 (cross-cohort generalization)

---

### R6 [PAPER 1 — NEW SCAFFOLD ROW] — Mechanism: lineage-TF collapse → DNMT / STAT3 / FOSL1 → TROP2 re-expression

> **Source.** `molecular_only_lock` §2 claims 4 (TF activity collapse: FOXE1 −1.21 / NKX2-1 −0.65 / STAT3 +1.55 / FOSL1 +0.98 / DNMT1 +0.74), 5 (TROP2 / TACSTD2 Δmean +4.33, FDR = 1.2 × 10⁻³²), 7 (Hallmark GSEA — IL6_JAK_STAT3 +8.0, OXPHOS −5.5, EMT / IFN-γ / TNF-α up), 8 (TF network — FOXE1-PAX8 r = +0.55, PAX8-DIO1 r = +0.64, HHEX-PAX8 r = +0.54).
> **Frame (per strategy memo §5 Layer 4 + §11).** "Consistent with a 4-step framework (TF collapse → DNMT silencing → STAT3 / AP-1 activation → TROP2 re-expression)" — NOT "demonstrates." Methylation absent; functional perturbation absent. Use `consistent with` / `supports`.
> **Body prose deferred** to author / Yu professor pass — not generated in this scaffolding update.
> → **Figure 4** (mechanism + TROP2; see Figure plan §below).

### R7 [PAPER 1 — NEW SCAFFOLD ROW] — Outcome and tumor-population TROP2 vulnerability

> **Source.** `molecular_only_lock` §2 claims 3 (PFI dichotomized HR = 2.04 [1.15–3.61], p = 0.015; DFI multivariate HR = 1.41, p = 0.025, BRAF / RAS / age / stage adjusted, n = 461), 5 (TROP2 DE), 6 (tumor vs normal — DM1 p = 1.8 × 10⁻²⁰; TROP2 1.4× tumor > normal, n = 517).
> **Frame.** "Supports an independent prognostic association" — NOT "proves." TROP2 = "tumor-level / tumor-population vulnerability" (Q3 spot-level ρ = −0.016 negative — DO NOT claim spot-level co-localisation). Sacituzumab / TROP2-ADC framing is hypothesis-generating, not "actionable."
> **Body prose deferred** to author / Yu professor pass.
> → **Figure 5** (outcome panel; see Figure plan §below).

### R8 [PAPER 1 — NEW SCAFFOLD ROW] — Spatial validation across 28 ST slides

> **Source.** `molecular_only_lock` §2 claims 1 (TCGA bulk DM1 vs THYROID_NONOVERLAP r = −0.885, n = 561), 9 (bootstrap 95% CI 12-slide [−0.997, −0.916]), 11 (Moran's I per slide mean 0.36, 26/28 perm p < 0.05), 12 (bivariate Moran DM1 × NONOVERLAP all 12 negative), 13 (margin gradient 23/28 negative ρ), 16 (Harmony 28-slide UMAP). Plus `paper1_spatial_supplement_freeze_2026_05_03.md` Verdict C (non-overlap lineage cross-val ρ = −0.85 to −1.00 in 12 external slides; autoimmune-no-cancer negative control GSE248205, HT / GD inflammation alone does not raise DM1_like above ctrl).
> **Frame.** "Supportive spatial validation" — NOT "spatial validation alone resolves." Cancer progression in GSE250521 is *direction-consistent but underpowered* (sample-mean Spearman ρ = −0.35, p = 0.18 in 16-slide cohort) — keep that hedging.
> **Body prose deferred** to author / Yu professor pass.
> → **Figure 6 (Supp candidate)** spatial 4-panel; see Figure plan §below.

### R9 [PAPER 1 — NEW SCAFFOLD ROW] — Honest negatives and pre-tested no-go

> **Source.** `molecular_only_lock` §2 claim 15 (N1 BRAF-/RAS- subset, n = 156, PFI HR = 1.20, p = 0.80 NS — score does not improve risk-stratification within the BRAF-/RAS- subset). Plus `2026_05_04_image_dm1_final_nogo_decision.md` (closure battery: H&E → DM1 r = 0.022, AUROC = 0.511, real ≤ random max — pre-tested and disclaimed). Plus methylation absent / functional validation absent (per `molecular_only_lock` §3).
> **Frame.** This row exists to make the honest-negative envelope visible to reviewers up front: "BRAF/RAS-negative dark matter" is NOT a Paper 1 headline; H&E-based image inference was pre-tested and will be explicitly disclaimed in Limitations; functional / methylation work is stated as future. Forbidden language list (per `molecular_only_lock` §5) reaffirmed.
> **Body prose deferred** to author / Yu professor pass — including the Limitations text which is voice-protected.
> → **Figure 5C** (N1 honest-negative inset, or move to Supp); `[Internal-only]` H&E / RunPod content explicitly excluded from submission package.

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

## Figures (Paper 1 = 5 main + Supp; legacy 7+6 plan moved)

### Main figures — Paper 1 scope (per `2026_05_04_8gene_curated_vs_denovo_final_strategy.md` §9)

- **F1 [PAPER 1]** Driver landscape orthogonality. A: driver mRNA × mutation status (BRAF / HRAS / NRAS / KRAS); B: single-feature driver AUC for DM cluster; C: TIERA67 univariate Cohen's d ranking with Driver_anchor highlighted at #52–67; D: Driver_anchor-only cluster ARI = −0.007.
- **F2 [PAPER 1]** Axis robustness. A: pan-genome top-N MAD ARI ladder (top-100 / -200 / -1000 / -5000) vs TIERA67 vs 8-gene vs Driver_anchor; B: TIERA67 enrichment in pan-genome top-100 (hypergeometric p = 3 × 10⁻⁴).
- **F3 [PAPER 1]** Axis anchor and non-overlap validation. A: TCGA bulk DM1_like vs THYROID_NONOVERLAP scatter (r = −0.885, n = 561) — non-overlap, 8-gene-free; B: TCGA DM1 cluster vs DM1_like score boxplot.
- **F4 [PAPER 1]** Mechanism + TROP2. A: TF activity volcano (FOXE1 / NKX2-1 / PAX8 down; STAT3 / FOSL1 / DNMT1 up); B: Hallmark GSEA bars (IL6_JAK_STAT3 / OXPHOS / EMT / IFN-γ / TNF-α); C: TROP2 / TACSTD2 DE volcano + drug-target overlay; D: tumor vs normal panel.
- **F5 [PAPER 1]** Outcome. A: KM PFI by DM1_like dichotomy; B: multivariate Cox forest (BRAF / RAS / age / stage / DM1_like); C (or Supp): N1 BRAF-/RAS- subset HR = 1.20 NS — honest-negative inset.

### Supp candidate — Paper 1 spatial layer

- **F6 [PAPER 1 — Supp]** Spatial. Moran's I per slide; bivariate Moran (DM1 × NONOVERLAP); 8-gene drop-one-out per slide; autoimmune-no-cancer negative control (GSE248205).

### F8 [PAPER 1 — NEW, 2026-05-04] External advanced-disease replication (Landa 2016 GSE76039, n=37)

> Source: `2026_05_04_gse76039_first_pass_report.md` (Steps 5–7); script `project/notebooks_or_scripts/p_landa_2016_first_pass.py`; raw artifacts in `project/results/p_landa_2016/`. Cohort: 17 PDTC + 20 ATC (no PTC controls in series); Affymetrix HG-U133 Plus 2.0 microarray, gcRMA log2 (Landa lab pre-processed); within-cohort z normalisation; **TCGA-trained absolute classifier explicitly NOT applied across platform** per `v17_korean_K2_calibration`. **Status:** Main figure candidate (3-panel composite); awaiting author/Yu approval for inclusion in submission package.

- **F8A** ATC-vs-PDTC box+strip across four axis variants (RAI_8 / DM1_like / THYROID_NONOVERLAP / TDS_like). Cohen's d range −3.32 to −3.57 (axis); MW p ≤ 6.3 × 10⁻⁷.
- **F8B** Scatter: DM1_like vs THYROID_NONOVERLAP across all 37 samples — **Spearman ρ = −0.925, p = 3.1 × 10⁻¹⁶** (zero-overlap lineage cross-validation; stronger than the TCGA-bulk reference r = −0.885 on n=561).
- **F8C** Mechanism heatmap (10 genes: FOXE1, NKX2-1, PAX8, HHEX, DNMT1, DNMT3B, STAT3, FOSL1, JUNB, TACSTD2; within-cohort z; samples ordered by histology then DM1_like). Direction-consistent with TCGA TF activity volcano (lock-claim 4) and TROP2 DE (lock-claim 5).

**Honest framing constraints (per `molecular_only_lock` §4 + first-pass report §6/§9):**
- Frame as "external **advanced-disease** replication" — NOT as "PTC → PDTC → ATC monotonic gradient" (no PTC controls in series).
- TROP2 d = +1.13 is the smallest effect of the seven axes; report honestly as "tumor-population vulnerability replication," NOT as "DM1-spot-level co-localisation" (Q3 NEG remains the spatial finding).
- Cross-platform note: within-cohort z is the only operation done on this cohort; the TCGA-trained absolute-form classifier is NOT applied here.
- No survival / mutation claims from this cohort (metadata absent in GEO).

### Supplementary figures — Paper 1 scope

- 28-slide Harmony UMAP (claim 16) · 8-gene drop-one-out per slide referee-Q (claim 10) · pan-cancer THCA outlier (claim 14) · TF network FOXE1-PAX8 / PAX8-DIO1 (claim 8) · bootstrap 95% CI 12-slide (claim 9) · margin gradient (claim 13) · N1 dark-matter subset KM + forest (claim 15, if not in F5C) · GSE250521 cancer-progression PT→ATC underpowered (spatial freeze §a) · GSE230424 PTC+HT non-overlap lineage cross-val (spatial freeze §a) · GSE248205 autoimmune-no-cancer negative control (spatial freeze §b) · closure battery negative-feasibility 1-page Methods supp.

### Internal-only — must NOT enter submission package

- §9 / §11 of `PAPER1_DM1_FULL_2026_05_04.md` (RunPod G1 / G2 / G3 — superseded by closure NO-GO)
- venue-probability framing
- TROP2 spot-level co-localisation (Q3 NEG)
- H&E-DM1 / WSI-projection / `H&E-inferable` / `morphology-derived DM1` in any active-claim form
- Fusion-map figure (no fusion calls available)

### `[PAPER 2/4 — moved]` superseded figure list (legacy 5/3 mega-paper plan; retained for cross-paper reference, NOT for Paper 1 submission)

- (legacy) **F1** Pan-Asian HLA forest meta (Pillar 1) → **Paper 4**
- (legacy) **F2** GSE286332 PTC+HT multi-panel (Pillar 2) → **Paper 2**
- (legacy) **F3** Driver mRNA neutrality (Pillar 3) → **Paper 1 (kept; absorbed into new F1 above)**
- (legacy) **F4** Pan-genome ARI + 8-gene clinical-interpretability (Pillar 4) → **Paper 1 (kept; absorbed into new F2 above)**
- (legacy) **F5** Autoimmune-PTC mechanism 4-panel (Pillar 5) → **Paper 2**
- (legacy) **F6** TCGA Hashimoto-like distribution + cross-cohort (Pillar 5 expansion) → **Paper 2**
- (legacy) **F7** DM1 sub-A vs sub-B + K2 NBNR clinical phenotype → **Paper 2**
- (legacy) **SF1** TIERA67 pan-genome rank ladder → **Paper 1 (Suppl)**
- (legacy) **SF2** Mediation Baron-Kenny detail → **Paper 2**
- (legacy) **SF3** Korean GSE213647 replication detail → **Paper 2**
- (legacy) **SF4** GSE286332 BCR repertoire detail → **Paper 2**
- (legacy) **SF5** ATA 2015 risk tier × DM cluster cross-tab → **Paper 1 (Suppl candidate, not yet built)**
- (legacy) **SF6** Sensitivity analyses (4-scenario forest) → **Paper 4**

---

## Author contributions (placeholder)

- **Seungho Cook (first/corresponding author):** Conceptualization, methodology, data curation, formal analysis, software, writing – original draft, writing – review & editing
- **Yu professor:** Supervision, project administration, resources, writing – review & editing

---

## Acknowledgments
- ARIA project allocation
- TCGA Research Network for thyroid cancer dataset
- Yoo SK / SNU-GMI (PRJEB11591), Lee SE / Macrogen (GSE213647), Lim DW / Dongguk (GSE286332)
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
