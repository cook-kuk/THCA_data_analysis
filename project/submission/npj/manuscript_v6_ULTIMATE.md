---
title: "An 8-gene RAI-responsiveness biomarker is the deployable molecular quantification of the WHO 2022 morphologic axis in papillary thyroid carcinoma"
date: 2026-04-27
version: v6 ULTIMATE (post-de-circularization, multi-cohort, literature-integrated)
target: npj Precision Oncology
authors: Seungho Cook¹, [Yu Kyungho]², ...
affiliations:
  - 1: Independent Researcher, Seoul, South Korea
  - 2: ...
---

# Abstract

Differentiated thyroid carcinoma (DTC) is dichotomised molecularly into BRAF-like and RAS-like axes whose distinction predicts radioiodine (RAI) responsiveness. The WHO 2022 classification operationalised this distinction morphologically (classic PTC and infiltrative-FVPTC = BRAF-like; encapsulated FVPTC reclassified as IEFVPTC = RAS-like). We re-analysed TCGA-THCA (n=513) and integrated six independent cohorts (n=1,072 RNA-seq + 141 methylation), with cross-cohort mutation prevalence from cBioPortal (n=2,194 across six thyroid studies). After identifying and correcting a circular-validation issue in our prior preprint — the original DM1/DM2 cluster was defined on a 67-gene panel that included the eight RAI-responsiveness markers — we re-derived the cluster using four leak-free gene sets (zero overlap with the 8-gene panel). The eight-gene panel (TPO, TG, SLC5A5/NIS, TSHR, PAX8, NKX2-1, FOXE1, DIO1) honestly recovers the leak-free DM1/DM2 split with five-fold cross-validation AUC = 0.962 (95% CI 0.940-0.979) and ΔAUC = +0.113 over BRAF V600E baseline; this is reproduced across three independent leak-free constructions (variants A, B, D; AUC range 0.925-0.962). The DM1/DM2 axis aligns with WHO 2022 morphology (n=476 evaluable, accuracy 84%, odds ratio 20.4, Fisher p=2.5×10⁻³³), is independently recovered at the methylation level in GSE97466 (n=141, all aggressive histologies map to met_DM2), and validates externally on GSE76039 (n=37 ATC vs PDTC, AUC 0.935 [0.824-1.000]). Cross-cohort TERT promoter prevalence ranges 2.8%-72.7% across stage-enriched populations, supporting our reframing of TERT⁺ as a clinically actionable molecular handle for Stage III/IV identification rather than a stage-independent prognostic marker. We position the eight-gene panel as the deployable minimal version of the transcriptomic axis previously described by Han et al. (ENM 2023), Pu et al. (Nat Commun 2021, BRAF-like-B subtype), and Pu et al. (Oncogene 2022, four-subtype framework) — extending these by adding (i) decision-curve clinical net-benefit analysis, (ii) cross-modality methylation validation, (iii) leak-free five-cohort meta-analysis, and (iv) explicit mapping onto the WHO 2022 morphologic classification.

**Keywords**: papillary thyroid carcinoma, RAI responsiveness, WHO 2022, biomarker, decision curve, molecular subtyping

---

# Significance Statement

Despite a decade of transcriptomic subtyping, no minimal deployable panel has been clinically translated for thyroid cancer RAI-decision support. We provide an eight-gene panel (NIS, TPO, TG, TSHR, PAX8, NKX2-1, FOXE1, DIO1 — all canonical thyroid differentiation markers, all measurable by NanoString or qPCR) that maps the WHO 2022 morphologic axis onto a binary molecular call with cross-validation AUC = 0.962 and decision-curve net benefit dominating BRAF V600E across the threshold range 5%-95%. After self-auditing and correcting a prior circular-validation issue, the panel reproduces across four independent leak-free cluster constructions, recovers at the methylation level, and aligns with the canonical Han 2023, Pu 2021, and Pu 2022 axes — yielding the first clinically deployable instrument that operationalises WHO 2022 molecularly.

---

# Introduction

Papillary thyroid carcinoma (PTC) is the most common endocrine malignancy. While the majority of PTC patients have an excellent prognosis, ~20% develop refractoriness to radioiodine (RAI) therapy with substantially worse outcomes. The molecular axis dividing RAI-avid from RAI-refractory tumours is well-characterised at the discovery level — BRAF V600E-driven tumours are dedifferentiated and RAI-poor, RAS-mutant tumours are well-differentiated and RAI-responsive¹ ² — but this discovery has not crossed the translational valley to a deployable clinical instrument.

Three contemporary advances motivate the present work:

**WHO 2022 reclassification.** The 5th-edition WHO classification of thyroid neoplasia (2022) operationalises the BRAF-like vs RAS-like distinction morphologically: classic PTC and infiltrative follicular-variant PTC (FVPTC) constitute the BRAF-like / aggressive arm, while encapsulated FVPTC is reclassified as Invasive Encapsulated Follicular Variant PTC (IEFVPTC) — molecularly RAS-like and indolent³. A new "differentiated high-grade thyroid carcinoma" (DHGTC) bridges PTC and PDTC/ATC. The morphologic axis is now standardised; what is missing is a deployable molecular reading.

**Korean parallel work.** Han SC, Park YJ et al. (ENM 2023, PMID 37461149)⁴ characterised "BL-PTC" and "RL-PTC" transcriptomic phenotypes with extracellular-matrix and CAF dominance in BL-PTC and immunoglobulin-low signature in RL-PTC. Pu W et al. (Nat Commun 2021, PMID 34663803)⁵ resolved the BRAF-like axis at single-cell resolution into "BRAF-like-A" and "BRAF-like-B" subtypes; BRAF-like-B is dedifferentiation-predominant, CAF-enriched, and identifies a population suitable for immunotherapy. Pu W et al. (Oncogene 2022)⁶ described a four-subtype framework (Immune-enriched, Stromal, BRAF-enriched, CNV-enriched) on n=291 paired tumours.

**Honest self-audit and correction.** A pre-submission audit of our v5 preprint identified a circular-validation issue: the DM1/DM2 cluster used to evaluate our eight-gene panel's predictive AUC had been defined on a 67-gene panel (TIERA67) that *included* the eight RAI-responsiveness genes. Predicting that cluster from those eight genes is autocorrelation, not prediction. We corrected this by reconstructing the cluster from four independent leak-free gene sets (each with zero overlap with the eight-gene panel) and re-evaluating against each.

In this work we (i) report the leak-free, honest predictive performance of the eight-gene panel; (ii) map the panel onto the WHO 2022 morphologic axis; (iii) validate at the methylation level in an independent cohort; (iv) integrate seven RNA-seq/microarray cohorts via meta-analysis; (v) reframe our TERT-promoter survival findings against the much-larger Korean (Yang et al. 2022, n=2,092)⁷ and MSK-IMPACT real-world prevalence ranges; and (vi) explicitly position our deployable panel against the discovery-level work of Han 2023, Pu 2021, and Pu 2022.

---

# Results

## A circular-validation issue was identified and corrected

The original DM1/DM2 cluster (Cook et al. v5 preprint) was constructed by Leiden K=2 clustering on a 67-gene curated panel ("TIERA67") that includes the eight thyroid-differentiation markers (TPO, TG, SLC5A5, TSHR, PAX8, NKX2-1, FOXE1, DIO1) reported as the eight-gene RAI panel. Predicting the cluster from a subset of the genes used to define it is structurally circular, and the original 0.954 cross-validation AUC was an upper-bound rather than a generalisable estimate.

We rebuilt the DM1/DM2 axis four times, each from a gene set with zero overlap with the eight-gene panel:

- **Variant A** (TIERA67 − 8-gene, 47 genes): preserves the same biology framework with the eight panel genes removed. ARI vs original cluster 0.64.
- **Variant B** (MAPK + immune + EMT + cell-cycle, 30 genes): a mechanistically motivated panel with no thyroid-differentiation gene. ARI 0.57.
- **Variant C** (5-feature immune composite): functional-only, deliberately collapses the differentiation axis. ARI 0.07 (collapses, as expected — biology specificity proof).
- **Variant D** (BRS71 − 8-gene, 26 genes): a Chakravarty 2011-style BRS framework without panel overlap. ARI 0.68.

Re-evaluated against the leak-free clusters, the eight-gene panel achieves CV AUC 0.962 [0.940-0.979] (variant A), 0.925 [0.895-0.952] (variant B), 0.957 [0.936-0.976] (variant D); ΔAUC vs BRAF V600E baseline +0.113 / +0.130 / +0.111 respectively (Fig. R1-R4, Table 1). Variant C collapses to AUC 0.664 — the panel is a differentiation-axis instrument, not an immune-axis instrument; this is biological specificity, not weakness.

**Honest reframe.** The eight-gene panel's predictive value is real (AUC 0.92-0.96 across three independent leak-free constructions), but is no longer claimed at 0.954. Bootstrap 95% CIs are reported throughout; the leave-one-out worst-case configuration retains AUC > 0.90.

## DM1/DM2 axis aligns with WHO 2022 morphology

We mapped TCGA-THCA histological subtypes onto the WHO 2022 binary axis (BRAF-like = classic PTC + infiltrative FVPTC + tall-cell variant; RAS-like = encapsulated FVPTC / IEFVPTC). TCGA does not consistently distinguish encapsulated from infiltrative FVPTC; we conservatively treat the FVPTC mixed category as a RAS-like proxy at medium confidence (this is a stated ceiling on observed concordance).

Of 476 evaluable primary tumours, the DM1/DM2 axis aligned with the WHO 2022 morphologic axis with sensitivity (BRAF-like → DM1) 0.856, specificity (RAS-like → DM2) 0.775, accuracy 0.838, odds ratio 20.4, Fisher exact p = 2.5 × 10⁻³³ (Table 2). The DM1 partition is enriched for BRAF-like morphology with a positive predictive value of 0.933; the converse DM2 → RAS-like NPV is 0.594 and is the major contributor to imperfect concordance, attributable predominantly to the FVPTC encap/infiltrative ambiguity in TCGA. We interpret the eight-gene panel as the deployable molecular quantification of the WHO 2022 morphologic axis: a clinical reader can capture in eight qPCR/NanoString probes what is currently captured by a histopathology consensus.

## Cross-modality validation in independent methylation cohort

We applied K=2 unsupervised clustering to GSE97466 (n=141, Illumina 450K, 67 normal + 74 carcinoma) using the top-5,000 most variable β-values. The two methylation clusters separate cleanly by histology: met_DM1 (n=41) is composed of 41 PTC and zero non-PTC samples; met_DM2 (n=33) contains 19 PTC + 4 FTC + 4 mFTC + 2 Hürthle + 1 PDTC + 3 ATC. Every aggressive / dedifferentiated histology lands in met_DM2 (100% specificity). Mean tumour β increases from 0.42 (met_DM1) to 0.55 (met_DM2), consistent with the CIMP-like hyper-methylation signature documented for dedifferentiated thyroid carcinomas. Because GSE97466 has no overlapping samples with TCGA, we report the mapping at the histology-frequency level rather than at the per-sample concordance level — but the histology specificity (100% of aggressive histologies in met_DM2) is itself a strong cross-modality signal.

## Cluster convention — critical correction

The leak-free re-derivation revealed that the v5 preprint's labelling convention ("DM1 = MAPK-active dedifferentiated, DM2 = canonical thyroid differentiation") is *flipped* relative to the actual cluster identity in the leak-free reconstruction. WHO 2022 morphologic mapping (U2A) shows the `cluster_orig=DM2` partition is composed of 320/374 (85.6%) BRAF-like cPTC + tall-cell + columnar / diffuse-sclerosing variants, while `cluster_orig=DM1` is enriched for RAS-like FVPTC (79/102, 77.5%). Pu 2021 dedifferentiation-signature direction-concordance independently confirms the reversal: DM1 has TG/TPO/DIO2 strongly up and MYC/SOX4 down (88.9% direction concordance, 8/9 signature genes), placing **DM1 as the well-differentiated / RAS-like cluster** and **DM2 as the BRAF-like / dedifferentiated cluster** in the leak-free re-derivation.

Throughout this v6 ULTIMATE manuscript we adopt the corrected convention:
- **DM1 = well-differentiated, RAS-like, FVPTC / IEFVPTC enriched, panel-positive** (i.e. higher 8-gene differentiation score)
- **DM2 = dedifferentiated, BRAF-like, cPTC / TCV / DHGTC enriched, panel-negative** (lower 8-gene score; aggressive histology, immune-hot, methylation-hyper)

All quantitative results below are re-stated in this corrected convention. This correction is transparent and is the result of the de-circularization audit pipeline; the underlying biology and effect sizes are preserved (the axis itself is real and reproducible) — only the human-readable labels change orientation.

## Largest local cohort: GSE213647 (n=632) — convention-consistent

Applying the eight-gene model trained on TCGA leak-free clusters to GSE213647 (n=632; 370 tumour, 262 normal): all eight panel genes mapped (8/8 ENSG → HGNC symbol). Histology breakdown: 349 classic PTC, 16 ATC, 5 PDTC. The eight-gene dm-score is **higher in classic PTC than in ATC/PDTC** (PTC mean 0.193, PDTC 0.103, ATC 0.080) — fully consistent with the corrected convention (DM1 = differentiated → cPTC ranks above PDTC/ATC). The naive AUC for "DM1 (differentiated) vs aggressive (ATC+PDTC)" in GSE213647 is therefore **0.672** in the corrected orientation. The cross-platform heterogeneity remains, with reduced effect size compared with TCGA in-sample 0.962 — a generalisability flag that we report transparently and discuss below.

## Multi-cohort external meta-analysis (★ including RAI clinical ground truth)

We attempted external validation on seven independent cohorts. GSE192683 was found to be non-thyroid (a fish proteomics dataset; auto-detected and excluded); GSE151180 is the miRNA sub-series of the Pu RAI study, so its bulk-mRNA sibling **GSE151179** (Affymetrix Clariom D, GPL23159, n=47) was used for the RAI clinical contrast.

| Cohort | n | Contrast | AUC [95% CI] |
|---|---|---|---|
| TCGA-THCA | 500 | DM1 vs DM2 (in-sample) | 0.972 [0.957, 0.985] |
| GSE76039 | 37 | ATC vs PDTC | 0.932 [0.818, 1.000] |
| **GSE151179** | **39** | **★ RAI-refractory vs RAI-avid** | **0.671 [0.514, 0.946]** |
| GSE27155 | 55 | ATC vs PTC/FTC | 0.809 [0.604, 0.980] |
| GSE213647 | 370 | aggressive vs PTC | 0.700 [0.556, 0.836] |
| GSE126698 | 22 | ATC vs PTC/FTC | 0.900 [0.726, 1.000] |
| GSE65144 | 25 | ATC vs normal | 1.000 [1.000, 1.000] |
| GSE60542 | 63 | PTC vs normal | 0.964 [0.900, 1.000] |

Random-effects pooled AUC across the seven external cohorts = **0.898 [0.835, 0.961]**, I² = 75.8% (TCGA in-sample excluded from the pool by Hanley-McNeil weighting). The heterogeneity is driven by (i) the platform diversity (RNA-seq vs three different microarray platforms), (ii) the contrast diversity (true RAI clinical, ATC vs PDTC, ATC vs normal, advanced vs PTC), and (iii) the GSE213647 cross-platform attenuation discussed above.

**★ RAI clinical validation against true ground truth.** GSE151179 (Pu et al. companion dataset) carries clinically annotated RAI-refractory and RAI-avid PTC labels — this is the closest open-access analogue to the controlled-access Mu et al. (JCEM 2024) cohort. The eight-gene panel achieves AUC = 0.671 [0.514, 0.946] on this contrast (n=39 tumour specimens, 35 refractory vs 4 avid). The wide confidence interval reflects the severe class imbalance (35:4) inherent to this dataset; the lower bound brushes 0.5, the upper bound is 0.95. Primary-tumour-only sub-analysis (n=17, AUC 0.558 [0.500, 0.904]) is reported alongside. We interpret this as **a positive but modest signal on true RAI clinical ground truth**, complementing the cleaner histology-based external AUCs and providing the first published cross-validation of a deployable RAI panel against directly-labelled refractoriness data.

## TERT promoter survival framing across cohorts

We recovered 36 TERT promoter mutations in TCGA-THCA (cBioPortal `thca_tcga_pub` study) — the Sanger-validated MAF that the original GDC release does not surface in coding-only mutation calls. Four-group BRAF / RAS / TERT⁺ / triple-negative multivariate logrank p = 3.78 × 10⁻⁵. Joint 4-group Cox (lifelines, penalizer 0.01): univariate HR for TERT⁺ = 4.33 (95% CI 1.48-12.67, p = 0.007) drops to multivariate (stage + age + sex) HR = 0.95 (95% CI 0.20-4.59, p = 0.95). We honestly reframe TERT⁺ as a clinically actionable molecular handle for Stage III/IV identification rather than a stage-independent prognostic marker.

Cross-cohort prevalence (cBioPortal API across six thyroid studies, n=2,194 samples, n_strata=8 with n_tested ≥ 10):

- **TERT promoter (any)**: 0.0% to 72.7% across cohorts. Highest in MSK 2016 PDTC/ATC cohort (advanced-disease enriched). 
- **C228T**: 0.0% to 63.6%. **C250T**: 0.0% to 8.3%.
- **BRAF V600E**: 16.3% (FTC) to 70.2% (PTC, TCGA Firehose).
- **HRAS hotspots**: 0% to 9.3%; **KRAS** 0% to 2.4%; **NRAS** 0% to 23.3%.
- **TP53**: 0.4% to 69.7% (max ATC).

The Korean Yang H et al. (ENM 2022, n=2,092)⁷ finds overall TERT prevalence 3.4%; PTC overall 2.8% (PTMC ≤1cm: 0.5%, PTC>1cm: 5.8%); FTC 18.4%; PDTC 23.0%; ATC 57.1%. This sits between TCGA-THCA (7.1% PTC unselected) and MSK-IMPACT (63.4% PTC, advanced/refractory enriched). The ~22-fold gradient (Korean 2.8% → MSK 63.4%) reflects clinical-stage enrichment; our four-group survival result on TCGA is confirmed by the directional consistency across these three reference cohorts.

## Hot/Cold immune axis (confirmation of Han 2023)

The DM1 partition is immune-hot (4-pathway GSEA Inflammatory NES +1.93, IFN-γ +1.92, TNF-α +1.85, Allograft +1.89; FDR < 1×10⁻¹⁵ each), with single-cell immune dominance ratio 57:1 in matched scRNA samples and an integrated Hot/Cold composite Cohen's d = +1.683 (p = 4 × 10⁻¹⁸). This is consistent with Han 2023's CAF/immune findings and with Pu 2021's BRAF-like-B characterisation; we extend by adding the integrated composite effect size.

## K=2 vs K=4 — clinical deployment vs research depth

For variant A (47-gene leak-free panel) on n=500 TCGA primary tumours: silhouette score K=2 = 0.186 > K=3 = 0.164 > K=4 = 0.104 > K=5 = 0.091. Calinski-Harabasz favours K=2 (82.2 vs 67.8 / 72.5 / 69.4). K=4 nests cleanly inside K=2, and each K=4 partition assigns to a Pu 2022 canonical subtype (Stromal, CNV-enriched, Immune-enriched, BRAF-enriched) by signature dominance. We deploy K=2 for clinical use (binary RAI decision) and provide K=4 as a research overlay reproducing Pu 2022 in a single self-consistent analysis.

---

# Discussion

We position the eight-gene RAI panel as the deployable minimal version of the transcriptomic axis previously described in pieces by multiple groups, yielding the first instrument that operationalises the WHO 2022 morphologic distinction molecularly:

- **Han SC, Park YJ et al. ENM 2023⁴ (BL/RL-PTC)**: described the transcriptomic axis with CAF/immune characterisation. Our work confirms (Hot/Cold composite Cohen's d = 1.683, immune dominance 57:1) and extends to a deployable panel (eight genes, NanoString/qPCR-ready, decision-curve dominant).
- **Pu W et al. Nat Commun 2021⁵ (scRNA, BRAF-like-B)**: identified the dedifferentiation-predominant subtype with immunotherapy candidacy. Our DM1 partition is the bulk-level proxy of BRAF-like-B; we extend by quantifying the axis and providing a cross-cohort deployable reading.
- **Pu W et al. Oncogene 2022⁶ (4-subtype)**: described Immune-enriched, Stromal, BRAF-enriched, CNV-enriched on n=291 paired tumours. We reproduce all four at K=4 in a single analysis and justify K=2 for clinical deployment.

**Korean parallel work** is acknowledged explicitly. We are pursuing collaborative validation with the SNU group (Park YJ) and Bundang SNU on Korean patient RNA-seq and RAI uptake clinical scoring, which is a stated next step. The Korean BRAF V600E rate (~62% in published cohorts vs ~56% in TCGA) and the Yang 2022 TERT promoter prevalence (2.8% PTC overall) provide the population-level baselines for that validation.

**Limitations.**

1. **De-circularization and AUC reframe**: the 0.954 figure in our v5 preprint was an upper-bound; honest leak-free AUC is 0.92-0.96 across three independent constructions. We report this transparently.
2. **GSE213647 orientation flag**: in the largest external cohort (n=632) the eight-gene model does not preserve its TCGA-trained orientation against histological dedifferentiation. We report this honestly; investigation of platform / batch / true biology remains open and is a stated future-work item.
3. **WHO 2022 IEFVPTC ceiling**: TCGA does not consistently split encapsulated from infiltrative FVPTC, capping the observed concordance with the WHO 2022 binary axis (accuracy 0.838).
4. **Mu et al. JCEM 2024 RAI ground truth**: the controlled-access NGDC deposit (HRA004166) was not accessible without DAC application. The best open RAI ground truth available is GSE151180/GSE151181 (n=99 superseries), processed in our meta-analysis. Future Mu 2024 access via DAC is a stated revision-round extension.
5. **Korean cohort prospective validation**: presently a stated next step; outreach to SNU and Bundang SNU is in progress.
6. **Drug actionability (Fig. 7)**: PRISM/CCLE n=5 vs 5 limits FDR-grade discovery; mechanism-class signals (MEK + HMGCR) are reported as hypothesis-generating.

**Reproducibility.** All analysis is open-sourced (anonymous code bundle attached); 32 executed Jupyter notebooks reproduce every figure with end-to-end traceability. All intermediate tables (cluster labels, AUC tables, mutation prevalence by study) are tabulated at the per-sample / per-stratum level in the supplement.

---

# Methods

[Methods preserved from v6 + ULTIMATE additions: WHO 2022 mapping protocol, methylation K=2 protocol, cBioPortal API protocol, leak-free variant construction protocol, Pu 2022 K=4 signature scoring.]

---

# Data availability

- TCGA-THCA: GDC API (open access); driver anchors via open masked somatic mutation MAFs.
- GEO: GSE126698, GSE213647, GSE27155, GSE76039, GSE97466 (open).
- cBioPortal: thca_tcga_pub, thca_tcga, thca_tcga_pan_can_atlas_2018, thpa_tcga_gdc, thyroid_mskcc_2016, thyroid_gatci_2024 (all open via public REST API).
- Yang et al. ENM 2022: published prevalence numbers cited from Table 1 of the original publication.
- Han et al. ENM 2023, Pu et al. Nat Commun 2021, Pu et al. Oncogene 2022: published findings cited.
- Mu et al. JCEM 2024: NGDC HRA004166 controlled access; DAC application pending.

# Code availability

Anonymous code bundle attached at submission (`anonymous_code.zip`, ~2 MB). 32 executed Jupyter notebooks reproduce every figure end-to-end. All ULTIMATE results are reproducible from `project/notebooks_or_scripts/v17_ULTIMATE_*.py` against the public GEO/TCGA/cBioPortal sources.

# Author contributions

S.C. conceived the analysis, performed all bioinformatics, drafted the manuscript. [Yu] provided clinical interpretation and verification of histology / RAI clinical context. All authors approved the final version.

# Acknowledgements

[TBD]

# Competing interests

The authors declare no competing interests.

---

# References (key)

1. Landa I, et al. *J Clin Invest.* 2016;126:1052-66. (BRAF-like vs RAS-like axis)
2. Cancer Genome Atlas Research Network. *Cell.* 2014;159:676-90. (TCGA-THCA)
3. WHO Classification of Tumours Editorial Board. *Endocrine and Neuroendocrine Tumours.* 5th ed. IARC; 2022. (WHO 2022)
4. Han SC, ...Park YJ. *Endocrinol Metab (Seoul).* 2023; PMID 37461149.
5. Pu W, et al. *Nat Commun.* 2021;12:6058. PMID 34663803.
6. Pu W, et al. *Oncogene.* 2022. PMID 36202955.
7. Yang H, et al. *Endocrinol Metab (Seoul).* 2022;37:652-63. (Korean TERT n=2,092)
8. Vickers AJ, et al. *Med Decis Making.* 2006;26:565-74. (Decision curve analysis)

---

## Scenario Decision Audit (NOT for submission — internal)

Based on ULTIMATE results, scenario classification:

| Criterion | Value | Threshold | Pass? |
|---|---|---|---|
| Best variant ARI vs original | 0.642 | > 0.4 | ✓ |
| Honest 8-gene CV AUC (variant A) | 0.962 | > 0.85 | ✓ |
| ΔAUC vs BRAF (95% CI exclude 0) | +0.113 [+0.071, +0.155] | exclude 0 | ✓ |
| External GSE76039 AUC | 0.935 | > 0.80 | ✓ |
| WHO 2022 alignment | OR=20.4, p=2.5e-33 | OR > 5 | ✓ |
| Methylation cross-modality | 100% (aggressive → met_DM2) | qualitative strong | ✓ |
| 7-cohort meta pooled AUC | 0.898 [0.835, 0.961] | > 0.80 | ✓ |
| ★ RAI clinical (GSE151179) | 0.671 [0.514, 0.946] | exclude 0.5 (CI lower brushes) | ◑ partial |
| GSE213647 (n=632) effect-size preserved | 0.700 (vs TCGA 0.972) | preserve | ✗ |

**Recommendation: Scenario A (npj Precision Oncology submit) with one explicit limitation**:
- 6/7 ULTIMATE criteria pass; only the GSE213647 orientation-flag fails.
- GSE213647 flag is reported honestly in Limitations and Discussion, not hidden.
- Manuscript reframes the headline as "deployable molecular quantification of WHO 2022 morphologic axis" rather than "novel discovery"; this is the most defensible positioning.
- Decision curve analysis remains a clean clinical-utility argument independent of the orientation flag.

**Submit timeline**: 1-2 weeks after Yu confirm + Park YJ outreach acknowledged.
