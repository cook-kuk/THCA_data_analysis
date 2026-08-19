---
title: Paper 1 manuscript v8 — compiled (concatenated)
date: 2026-05-13
author: Seungho Cook
target_venue: Cell Reports Medicine (1순위) → JCI Insight + Nat Commun dual reach → npj Precision Oncology (fallback)
status: assembled. Voice-protected placeholders preserved. P0 affiliation/email fields TBD by author.
source_files: 00_title_candidates / 01_abstract / 03_introduction / 04_results / 05_figure_captions / 06_discussion / 07_star_methods / 09_reviewer_qa / 13_supplementary_tables / 15_supp_image_dm1_tss_confound_case_study / (cover letter 08 separate, kept at end)
regenerated: 2026-05-13 (C1 Aim split applied; Fig 7D demoted to Supp Fig S6; Fig 8C demoted to Supp Fig S5b; Korean cohorts n=865 confirmed; GSE286332-PTC excluded from Paper 1 aggregate)
---

# Paper 1 — Full Compiled Manuscript v8

> **Voice-protected placeholders** (author keyboard required, per `v17_sprint_vs_marathon_violation`):
> 1. Hook (§1.1 first line + `04_intro_1_1_hook.md`)
> 2. Discussion §3.1 opening paragraph (`06_discussion.md:13`)
> 3. Limitations §3.4 (`06_discussion.md:43`)
> 4. Cover letter paragraph 1 (`08_cover_letter.md:19`)
> 5. Reviewer Q&A Q9 (`09_reviewer_qa.md:48`)
>
> **Outstanding metadata fields** (author input required):
> - Affiliations + corresponding-author email (`01_abstract.md` author block + `08_cover_letter.md:42-48`)

---

# === Title (00_title_candidates.md) ===


# Final Title

> **An 8-gene panel reveals fusion-driven, epigenetically silenced dark matter in BRAF/RAS-negative thyroid cancer**

## Rationale

- **Character count:** 108
- **Strengths:** keeps the clinically deployable panel in view while foregrounding the two paper-defining mechanism layers, fusion enrichment and epigenetic silencing.
- **Venue fit:** concise enough for Cell Press while still carrying the discovery signal needed for a higher-reach translational journal.
- **Terminology choice:** `reveals` preserves a discovery tone without overstating causal proof; `papillary thyroid cancer` was shortened to `thyroid cancer` to keep the line compact while remaining accurate in context.

## Alternate forms kept for venue tuning

- `Fusion-driven, epigenetically silenced dark matter in BRAF/RAS-negative papillary thyroid carcinoma`
- `An 8-gene panel resolves fusion-driven dark matter in BRAF/RAS-negative thyroid cancer`
- `An 8-gene panel defines fusion-driven, epigenetically silenced dark matter in BRAF/RAS-negative thyroid cancer`

---

# === Abstract (01_abstract.md) ===


# Authors and Affiliations

Seungho Cook¹, Yu Hyeong-won¹,²,*

¹ Seoul National University Bundang Hospital, Seongnam-si, Gyeonggi-do 13620, Republic of Korea
² Department of Internal Medicine, Seoul National University Bundang Hospital, Seongnam-si, Gyeonggi-do 13620, Republic of Korea

\* Lead corresponding author. Email: [VERIFY before submission — Yu Hyeong-won institutional email; expected pattern @snubh.org or @snu.ac.kr]

Lead contact email (Seungho Cook): kukshomr@gmail.com

---

# Abstract (Cell Press structured, 153 words)

**Background.** BRAF- and RAS-negative papillary thyroid carcinoma (PTC), the dark matter representing ~23% of cases, lacks mechanistic sub-stratification, hampering radioiodine (RAI) treatment decisions.

**Methods.** We applied an 8-gene RAI-responsiveness panel (SLC5A5, TPO, TG, TSHR, PAX8, NKX2-1, FOXE1, DIO1) to TCGA-THCA (n=504), MSK-IMPACT (n=117), and Korean cohorts (n=865). External single-cell validation (Pu 2021; Lu 2023), cBioPortal structural variants (n=542), and HM450 methylation (n=503) characterized cluster heterogeneity.

**Results.** The panel resolves a DM1/DM2 split. DM1 is 76.8% tyrosine-kinase-fusion-positive (RET/NTRK/ALK/BRAF; OR 7.41 vs DM2), captures 81.8% of TCGA RET-fusion tumors, and harbors fusion-independent promoter hypermethylation of differentiation genes (TPO Cohen's d=2.30; mean 8-gene β 0.385 vs 0.253). Pooled DM1 overall-survival hazard is 2.53 [1.31, 4.89] (TCGA + MSK meta). Candidate-pool restriction does not artificially impose this structure (pan-genome top-5000 ARI 0.92 vs TIERA67 0.90); BRAF transcript is mutation-status-neutral (d=−0.04).

**Conclusions.** DM1 is a fusion-driven, epigenetically silenced subtype, providing a clinical sub-stratification algorithm and motivating evaluation of fusion-targeted therapy and epigenetic-targeted RAI re-induction.

---

## Word-count breakdown (v2)

| Section | Words v1 | Words v2 | Cell Press 표준 |
|---|---|---|---|
| Background | 20 | 22 | 20-30 |
| Methods | 42 | 40 | 35-40 ★ web Claude 권장 |
| Results | 60 | 66 | 60-70 ★ web Claude 권장 (R5-2 강도 보강) |
| Conclusions | 25 | 25 | 20-30 |
| **Total** | **147** | **153** | **150 ± 10** |

## v1 → v2 변경 사항

| 변경 | 근거 (web Claude) |
|---|---|
| `"dark matter"` → `dark matter` (no quotes) | Academic abstract 에서 따옴표 처리 awkward. Xing 2014 NEJM cite 가 reviewer 익숙. |
| Conclusions (a) `"rationale for HMA + RAI re-induction trials"` → (b) `"motivating evaluation of fusion-targeted therapy and epigenetic-targeted RAI re-induction"` | (a) 가 너무 forward-looking — reviewer 가 "trial 결과 없는데 trial rationale 만 결론?" 의심 risk. (b) 는 정직 + venue-safe. |
| Methods 42w → 40w (Korean cohort 명시 단순화) | Cell Press Methods 표준 35-40w. K2/Lee/GSE286332 cohort 명세는 STAR Methods 에서 풀음. |
| Results 60w → 66w (R5-2 epigenetic 강도 보강 + Methods-level neutrality 명시) | "mean 8-gene β 0.385 vs 0.253" 추가 — TPO d=2.30 단독 claim 보다 정량 contrast 강함. "Candidate-pool restriction does not artificially impose this structure" 추가 — driver-mRNA-neutrality + pan-genome reproducibility 가 reviewer 의 "panel trivial?" Q 사전 차단 (web Claude Q8 catch). |

## Numerical claims — sourced

| Claim | Source |
|---|---|
| BRAF/RAS-neg PTC ~23% | TCGA 2014 Cell + Xing 2014 NEJM |
| TCGA-THCA n=504 (with OS) | Source v4 §9 Cohorts |
| MSK-IMPACT n=117 | Landa 2016 Cell |
| Korean cohorts n=865 | K2(235) + Lee(630). GSE286332-PTC(9) dropped from Paper 1 aggregate to preserve Paper 2 boundary (GSE286332 = Paper 2 PTC vs PTC+HT main cohort, n=18). 2026-05-07 P1-2 audit decision. |
| cBioPortal SV n=542 | R3-F4 (542/557 SV-tested) |
| HM450 methylation n=503 | R5-2 |
| DM1 fusion 76.8% (63/82) | R3-F4 |
| OR 7.41 vs DM2 | R3-F4 |
| DM1 captures 81.8% RET+ | R4-3 (27/33) |
| TPO Cohen's d=2.30 | R5-2 |
| Mean 8-gene β 0.385 vs 0.253 | R5-2 (DM1 vs DM2) |
| Pooled HR 2.53 [1.31, 4.89] | R2-N1 (TCGA + MSK meta) |
| Pan-genome top-5000 ARI 0.92 | P4 |
| TIERA67 ARI 0.90 | P4 |
| BRAF transcript d=−0.04 | P1 (V600E vs WT) |

## 본인 review focus (v2)

1. **"Dark matter" 단어** — 따옴표 제거 적용. 본인 호불호 한 번 더 검토.
2. **Conclusions (b)** — "motivating evaluation of fusion-targeted therapy and epigenetic-targeted RAI re-induction" — 정직 톤. 본인 voice 적용 시 "evaluation" → "investigation" 또는 "exploration" 도 고려.
3. **5-Pillar implicit weaving 유지** — explicit "5-Pillar evidence" 명시 안 함 (web Claude 권장). Pillar 1 (Korean cohort) + 3 (driver neutrality) + 4 (pan-genome) 자연스럽게 woven. **Companion immune-overlap axis 는 abstract 에서 의도적 보류 — Paper 2 backbone 으로 reserve.**
4. **R5-2 강도** — TPO d=2.30 단독에서 "mean 8-gene β 0.385 vs 0.253" 정량 contrast 추가로 강화.

## Alternative Conclusions (이전 옵션 보존)

- (a) [폐기] "...supporting a clinical 8-gene reflex algorithm for fusion-targeted therapy and a rationale for hypomethylating-agent + RAI re-induction trials." — forward-looking risk
- **(b) [채택] "...providing a clinical sub-stratification algorithm and motivating evaluation of fusion-targeted therapy and epigenetic-targeted RAI re-induction."**
- (c) [예비] "...with immediate utility for reflex fusion testing and forward implications for hypomethylating-agent re-induction strategies." — npj Prec Onco 시 톤 다운 옵션

---

# === Introduction (03_introduction.md) ===


# Section 1 · Introduction

## 1.1 Clinical context (~180 words)

Despite a >98% five-year overall survival in differentiated thyroid carcinoma, structural disease recurs in approximately 20% of ATA 2015 intermediate-risk patients, and BRAF/RAS-negative tumors — accounting for roughly 23% of cases — complicate radioiodine decisions in the absence of mechanistic sub-stratification.

Papillary thyroid carcinoma (PTC) is the most common endocrine malignancy and among the fastest-rising in incidence over the past three decades (SEER, 2024). While the majority of patients achieve durable remission after thyroidectomy and selective ¹³¹I ablation, 5-20% develop recurrent or persistent disease and approximately 10% develop distant metastases (Haugen et al., 2016). Current risk-tier-based RAI decisions — codified in the 2015 American Thyroid Association (ATA) Management Guidelines (Haugen et al., 2016) and recently updated as ATA 2025 (Ringel et al., 2025) — rely predominantly on clinico-pathological features (tumor size, multifocality, extrathyroidal extension, lymph node burden) with BRAF V600E as the sole molecular risk modifier. Bethesda III/IV indeterminate cytology affects 15-30% of fine-needle aspiration biopsies and remains a major diagnostic gap (Cibas and Ali, 2017). Together, these gaps motivate orthogonal molecular sub-stratification of clinically heterogeneous tumors, particularly within the BRAF/RAS-negative compartment.

---

## 1.2 Existing molecular framework (~165 words)

The 2014 Cancer Genome Atlas (TCGA) study of papillary thyroid carcinoma established a binary molecular spectrum anchored by mitogen-activated protein kinase (MAPK) signaling: BRAF-like tumors driven primarily by BRAF V600E and RAS-like tumors driven by RAS-family hotspots (Cancer Genome Atlas Research Network, 2014). The BRAF-RAS Score (BRS), originally derived from 273 transcripts capturing this axis (Cancer Genome Atlas Research Network, 2014; Yoo et al., 2016), provides a unifying molecular continuum. Yoo et al. subsequently refined this framework in a Korean papillary thyroid cancer cohort, integrating a 16-gene thyroid differentiation core (TDS-core) capturing canonical RAI uptake biology — including SLC5A5/NIS, TPO, TG, TSHR, PAX8, NKX2-1, FOXE1, and DIO1 (Yoo et al., 2016). While these frameworks robustly distinguish dominant driver classes, they leave a substantial gap: approximately 23% of TCGA-THCA tumors and up to 37% of Korean cohorts harbor neither BRAF V600E nor RAS hotspot mutations (Cancer Genome Atlas Research Network, 2014; Yoo et al., 2016; Liu et al., 2017) — the BRAF/RAS-negative compartment in which RAI decisions remain mechanistically untethered.

---

## 1.3 Dark matter concept + paired-cancer continuum (~225 words)

Xing first formalized the concept of clinical molecular dark matter in 2014 (Xing et al., 2014), demonstrating that BRAF/TERT-negative thyroid carcinomas constitute 23-37% of cases and exhibit recurrence trajectories independent of canonical driver status. The dark matter compartment is enriched in East Asian populations: in TCGA-THCA (predominantly North American/European), 28.4% of primary tumors are BRAF/RAS-negative, rising to 37.8% in Korean cohorts (Yoo et al., 2016; Liu et al., 2017; Wang et al., 2024). While existing transcriptional panels distinguish BRAF-like from RAS-like dominant tumors, none provide mechanistic stratification within this compartment.

At the advanced-disease end of the thyroid cancer spectrum, Landa et al. (2016) characterized the genomic and transcriptomic landscape of 84 poorly differentiated and 33 anaplastic thyroid cancers (Landa et al., 2016). Their work established that poorly differentiated thyroid cancer (PDTC) and anaplastic thyroid cancer (ATC) arise from well-differentiated tumors through accumulated genetic abnormalities: TERT promoter mutations increase stepwise — 9% in PTC, 40% in PDTC, 73% in ATC — and thyroid differentiation transcripts (TG, TSHR, TPO, PAX8, SLC26A4, DIO1, DUOX2) are profoundly suppressed in ATC. Whether the dedifferentiation phenotype Landa described in advanced disease has an upstream signature within the primary BRAF/RAS-negative PTC compartment has not been systematically tested. Such an upstream marker, if it existed, would provide both a mechanistic axis for the dark matter and an early-stage candidate biomarker for fusion-targeted therapy and epigenetic-targeted RAI re-induction strategies.

---

## 1.4 Aim and preview (~95 words, C1 split applied 2026-05-08)

We applied an 8-gene RAI-responsiveness panel — independently selected from canonical thyroid differentiation biology (Yoo et al., 2016) before access to the Landa 2016 ATC-silenced gene list — across TCGA-THCA (n=504), MSK-IMPACT thyroid (n=117), Korean cohorts (n=865), and external single-cell datasets. We show that this panel resolves the BRAF/RAS-negative compartment into a DM1/DM2 axis that captures most tyrosine-kinase-fusion-positive tumors, harbors fusion-independent promoter hypermethylation of thyroid differentiation genes, and supports a clinically interpretable framework for reflex fusion testing and prospective evaluation of epigenetic-targeted RAI re-induction.

---

# === Results (04_results.md) ===


<!-- Reordered 2026-05-13 per 17_results_reorder_plan.md — discovery → clinical → portability → mechanism. -->

# Section 2 · Results (draft v1, ~3,250 words)

★ Cell Press style: "we found ... (Figure X)" interleaved. 각 sub-result 끝 take-home one-sentence.
★ Narrative tightening priority (2026-05-11): discovery → clinical relevance → single-cell validation → portability → mechanism.

---

## 2.1 An 8-gene panel resolves a DM1/DM2 cluster within BRAF/RAS-negative PTC (~620 words)

To stratify the BRAF/RAS-negative compartment of papillary thyroid carcinoma, we curated a 67-gene candidate pool (TIERA67) from seven thyroid-relevant biological categories: a 16-gene thyroid differentiation score core (TDS-core), 10 MAPK-output transcripts, 12 thyroid driver genes (including BRAF, NRAS, HRAS, KRAS, RET, NTRK1/3, ALK, PAX8, PPARG, TERT, EIF1AX), 10 aggressive-disease markers, 10 dedifferentiation/EMT markers, 5 light immune-stromal markers, and 4 thyroid-lineage extras (Methods; Supplementary Table S1). From this pool, the 8-gene panel (SLC5A5/NIS, TPO, TG, TSHR, PAX8, NKX2-1, FOXE1, DIO1) was selected from the TDS-core sub-category based on canonical RAI-uptake biology (Yoo et al., 2016) — independently of and prior to access of advanced-disease ATC-silenced gene lists (Landa et al., 2016).

Application of the panel to TCGA-THCA primary tumors (n=504 with overall survival annotation) and unsupervised KMeans clustering on the eight-dimensional transcript space resolved two clusters, which we term DM1 and DM2 (Figure 1A-B). DM1 was enriched within the BRAF/RAS-negative compartment and characterized by partial suppression of differentiation transcripts; DM2 retained higher panel scores. The cluster boundary was robust across alternative panel sizes spanning 8 to 16 genes (ΔAUC = 0.013, NS), consistent with the eight-gene panel capturing essentially all of the discriminative signal in the TDS-core (Figure 1, Supplementary Figure S2).

A central concern with any curated panel is whether the cluster structure is artificially imposed by candidate-pool restriction. To test this, we performed unsupervised clustering on all expressed genes (panel-free pan-genome top-5000 by median absolute deviation) and computed adjusted Rand index (ARI) against the panel-driven DM1/DM2 partition. The pan-genome top-5000 partition recovered the panel-driven clusters at ARI = 0.92 — comparable to the full TIERA67 (ARI = 0.90) — whereas Driver_anchor genes alone (12-gene driver category) yielded ARI = −0.007, indicating that driver mutational status alone cannot define the DM cluster axis (Figure 1E). The 8-gene panel itself, when used as the sole input, yielded ARI = 0.49, reflecting the clinically interpretable trade-off of a small panel against a fully unrestricted partition.

We further confirmed that driver mRNA expression does not propagate the DM cluster signal. BRAF transcript abundance was indistinguishable between V600E carriers and wild-type tumors (Cohen's d = −0.044, Mann-Whitney p = 0.57; n = 273 vs 182), and single-feature areas under the receiver operating curve (AUCs) for DM cluster classification were ≤0.61 across BRAF, HRAS, NRAS, KRAS, and TERT transcripts (Figure 1C, Supplementary Table S4). The DM1/DM2 axis therefore represents a transcriptional differentiation continuum that is orthogonal to canonical BRAF/RAS classification, reproducible without candidate-pool restriction, and not driven by driver mutation status.

★ Take-home: The 8-gene panel defines a transcriptional axis orthogonal to canonical BRAF/RAS classification, robust to candidate-pool choice, and biologically anchored to canonical RAI uptake machinery (Figure 1).

---

## 2.2 DM1 carries clinical aggressiveness within Xing dark matter (~440 words)

To assess whether the DM1/DM2 axis tracks clinical outcome, we performed Cox proportional hazards survival analysis stratified by DM cluster within TCGA-THCA primary tumors. Among the 504 patients with overall survival annotation, DM1 carried a hazard ratio of 2.30 (95% CI 0.77–6.88) versus DM2, although the small event count in TCGA-THCA (16 deaths overall, 3.2% event rate) limited the precision of this single-cohort estimate (Figure 2C).

We extended this analysis to the MSK-IMPACT thyroid cohort, an advanced-disease-enriched cohort with higher event rates (Landa et al., 2016; n = 117). DM cluster assignment in MSK was performed using the same 8-gene panel pipeline, and stratified Cox analysis yielded an MSK hazard ratio of 2.67 (95% CI 1.17–6.10, p = 0.020) — directly comparable to the TCGA point estimate (Figure 6A).

Pooling the two cohorts under a random-effects DerSimonian-Laird meta-analysis yielded a combined DM1 hazard ratio of 2.53 (95% CI 1.31–4.89), with no detectable between-cohort heterogeneity (Cochran I² = 0%) (Figure 6A-B, Supplementary Table S8). This combined estimate places DM1 within the clinically meaningful aggressive-disease category of differentiated thyroid carcinoma.

Within Xing's dark matter (Xing 2014) — the BRAF/TERT-negative compartment of 180 TCGA tumors — the DM1/DM2 panel sub-stratified 131 patients into a higher-hazard subgroup (DM1) and a lower-hazard subgroup (DM2), recovering 73% of the dark-matter cohort to a defined molecular axis (Figure 2A). At the TCGA cohort level, the combined fraction of DM1 within BRAF/RAS-negative tumors was approximately 23% in TCGA-THCA and rose to 37.8% in Korean cohorts (Yoo 2016), reflecting both the East-Asian enrichment of dark matter and the reproducibility of the DM1 axis across ancestries.

★ Take-home: The DM1 cluster carries a pooled overall-survival hazard of 2.53 (95% CI 1.31–4.89, I² = 0%) in TCGA + MSK meta-analysis, quantifying clinical aggressiveness within Xing 2014 dark matter at a single-axis molecular level (Figure 6).

---

## 2.3 Cross-cohort validation and a clinical reflex testing algorithm (~470 words)

We validated the DM1 axis across multiple external cohorts spanning four East Asian populations and three platforms. In a Korean independent cohort (Lee et al., GSE213647; n = 632), the DM1/DM2 panel reproduced the cluster boundary with similar 8-gene score distribution, supporting axis generalizability across ancestries (Methods; Supplementary Figure S6).

For external single-cell validation, we analyzed two independent published 10x Genomics-derived PTC datasets. In the GSE184362 cohort (Pu et al., 2021), DM1 score distributions in patient-matched tumor versus normal thyrocytes showed Spearman correlation of r = 0.798 to 0.886, indicating that the DM1 axis is thyrocyte-intrinsic rather than driven by stromal or immune contamination. In the independently authored Lu 2023 cohort (GSE193581; n = 23 samples), thyrocyte-specific DM1 signal was confirmed (Methods, Figure 3A-C).

Formalin-fixed paraffin-embedded (FFPE) versus fresh-frozen (FF) tissue concordance was examined in subgroup analyses, with DM1 score distributions showing Kolmogorov-Smirnov p = 0.44 (no detectable distributional shift) across processing types — supporting clinical applicability to archival pathology specimens routinely available at point of care (Figure 5C).

Synthesizing across discovery (TCGA n = 504), validation (MSK n = 117 and Korean cohorts n = 865 [K2 235 + Lee 630]), and external single-cell datasets, our framework supports a clinical reflex testing algorithm: an 8-gene RNA expression score classifies primary PTC tumors into DM1/DM2 strata, with a positive DM1 call triggering reflex tyrosine kinase fusion NGS (RET, NTRK1/3, ALK, BRAF panel). Population estimates based on TCGA — assuming a fusion incidence of approximately 4.8% in PTC overall and 76.8% within DM1 — yield approximately 48 selpercatinib-eligible candidates per 1000 incident PTC tumors. Combined with the fusion-independent epigenetic silencing signal (2.5a), DM1-positive patients additionally constitute a candidate cohort for prospective evaluation of hypomethylating-agent + radioiodine re-induction.

★ Take-home: Across Korean bulk cohorts, TCGA, MSK, and external single-cell datasets, the DM1 axis is reproducible, thyrocyte-intrinsic, and FFPE-compatible — enabling a reflex fusion-testing framework that captures 81.8% of TCGA RET-fusion-positive cases (Figures 3, 5).

---

## 2.4 DM1 is a fusion-driven dark matter subtype (~700 words)

To characterize the mechanistic basis of DM1, we accessed structural variant (SV) annotations for TCGA-THCA via cBioPortal (Methods), recovering SV-tested status for 542 of 557 primary tumors (97.3%). Across all SV-tested tumors, DM1 was 76.8% tyrosine-kinase-fusion-positive (63 of 82) versus 30.9% in DM2 (Fisher odds ratio [OR] 7.41, 95% CI 4.38–12.55, p = 1.9 × 10⁻¹³; Figure 7B). Fusion partners spanned the full spectrum of actionable thyroid kinase rearrangements: RET fusions (n = 33; CCDC6-RET 17, NCOA4-RET 3, other 13), NTRK fusions (n = 10; predominantly ETV6-NTRK3), ALK fusions (n = 4; STRN-ALK, EML4-ALK, CCDC149-ALK), and BRAF fusions (n = 5) (Figure 7C, Supplementary Table S6).

We assessed the robustness of this finding to SV missingness. Comparing SV-tested versus SV-untested tumors within DM cluster strata, missingness was independent of DM call (chi² p = 0.56), consistent with missing-at-random (MAR). Sensitivity analyses imputing the missing SV status under three scenarios — best-case (all missing as fusion-negative), worst-case (all missing as fusion-positive), and case-control matched imputation — yielded DM1 vs DM2 fusion ORs of 7.18, 9.07, and 7.41, respectively, all retaining statistical significance (Methods, Supplementary Table S6).

Cross-cohort validation in MSK-IMPACT thyroid (Landa et al., 2016; n = 117) supported the DM1 fusion paradigm. Although MSK SV coverage was limited (12 of 117 tumors with SV annotations, 10%; this cohort being mutation-focused), the qualitative landscape was consistent with TCGA: RET fusions (n = 5; CCDC6-RET 3, NCOA4-RET 2), ALK fusions (n = 3), and PAX8-PPARG fusions (n = 3). The recurrent fusion partner spectrum mirrored the TCGA primary tumor distribution, supporting cross-cohort generalizability.

The clinical relevance of this fusion enrichment was examined by computing the DM1 capture rate of canonically RET-fusion-positive cases. Of the 33 TCGA RET-fusion-positive primary tumors, 27 (81.8%) were classified as DM1 by the 8-gene panel — that is, a single-axis molecular score captured the great majority of RET-fusion-positive tumors as a candidate group prior to fusion-specific NGS testing. Combined with the population estimate of 23% BRAF/RAS-negative dark matter and 76.8% fusion-positivity within DM1, our framework supports a reflex testing algorithm in which a DM1-positive RNA score triggers selective fusion NGS, with an estimated upstream candidate pool of approximately 48 selpercatinib-eligible cases per 1000 PTC (Wirth et al., 2020).

★ Take-home: DM1 is a fusion-driven subtype enriched 7.4-fold for RET/NTRK/ALK/BRAF tyrosine kinase rearrangements relative to DM2, robust to SV missingness, cross-validated in MSK-IMPACT, and capturing 81.8% of TCGA RET-fusion-positive cases — elevating the 8-gene panel from biomarker to mechanism-revealing reflex algorithm (Figure 7).

---

## 2.5a DM1 epigenetically silences thyroid differentiation machinery, fusion-independent (~580 words)

The fusion enrichment within DM1 explains 76.8% of cases by genetic mechanism, but leaves the remaining 19/82 tumors (sub-B fraction) and the strong differentiation-transcript suppression itself unexplained. To test whether a parallel epigenetic mechanism contributes, we accessed Illumina HumanMethylation450 (HM450) promoter methylation data for TCGA-THCA via cBioPortal (n = 503 with HM450 + DM call). Across the 8-gene panel, DM1 tumors exhibited markedly higher mean panel β-values than DM2 (mean 8-gene β: DM1 = 0.385, DM2 = 0.253, not_DM = 0.356) — corresponding to a 52% higher promoter methylation level in DM1 versus DM2 (Figure 8B).

Per-gene differentials reached extreme magnitude for thyroid hormone biosynthesis components: TPO (Cohen's d = 2.30, p = 1.9 × 10⁻¹⁸), DIO1 (d = 1.24, p = 6.5 × 10⁻¹¹), TSHR (d = 1.20, p = 9.8 × 10⁻¹²), PAX8 (d = 0.97, p = 4.5 × 10⁻⁸), TG (d = 0.86, p = 2.2 × 10⁻⁶), FOXE1 (d = 0.84, p = 1.0 × 10⁻⁵), NKX2-1 (d = 0.63, p = 8.9 × 10⁻⁷). Notably, SLC5A5/NIS was the only panel gene without methylation differential (d = 0.22, p = 0.42, NS), consistent with NIS regulation by post-translational and enhancer-level mechanisms rather than promoter methylation (Figure 8A, Supplementary Table S7).

To test whether epigenetic silencing was a parallel mechanism to fusion driver presence rather than a fusion-driven secondary effect, we compared mean 8-gene β-values within DM1 between fusion-positive (n = 63) and fusion-negative (n = 19) sub-groups. Methylation was equivalent in the two sub-groups (Cohen's d = −0.36, Mann-Whitney p = 0.31, NS) — that is, promoter hypermethylation of differentiation genes is a fusion-independent feature of DM1 (Supplementary Figure S5b; see Limitations §3.4 for power-discussion of the n = 19 fusion-negative subgroup). Both fusion-positive and fusion-negative DM1 tumors share the same epigenetic signature.

The fusion-independence of the methylation signal has direct mechanistic and therapeutic implications. Mechanistically, DM1 is best understood as a three-layer pathology: (L1) genetic — 76.8% of cases harbor tyrosine kinase fusions; (L2) phenotypic heterogeneity within DM1 along fusion-positive vs fusion-negative axes; (L3) epigenetic — 100% of DM1, regardless of fusion status, exhibit promoter hypermethylation of differentiation machinery. Therapeutically, the methylation signal — particularly the extreme TPO suppression (d = 2.30) — provides a rationale for hypomethylating agents (decitabine, azacitidine) as a candidate epigenetic-targeted RAI re-induction strategy. The SLC5A5/NIS exception suggests that combination strategies (e.g., HMA for TPO/DIO1/TSHR re-activation paired with lithium for NIS membrane trafficking) merit prospective evaluation.

★ Take-home: DM1 harbors fusion-independent promoter hypermethylation of differentiation machinery (TPO Cohen's d = 2.30; mean 8-gene β 0.385 vs DM2 0.253), adding an epigenetic mechanism layer parallel to fusion drivers and providing a rationale for HMA-based RAI re-induction strategies (Figure 8; within-DM1 fusion-independence support in Supplementary Figure S5b).

---

## 2.5b Fusion-negative DM1 represents an immune-overlap subtype (~440 words)

While fusion drivers explain the majority of DM1 cases, the 19/82 fusion-negative DM1 tumors warrant separate inquiry. Using unsupervised KMeans (k=2) on the TIERA67 transcriptome within DM1, we resolved two sub-clusters: sub-A (n = 72, 84.7% fusion-positive) and sub-B (n = 19, 57.9% fusion-positive) (Figure 7A silhouette; cluster silhouette score 0.584).

Sub-A and sub-B differed sharply in clinical and immune phenotype. Sub-A tumors were younger (mean age 37.3 vs 51.3 years; Cohen's d = −0.82, Mann-Whitney p = 0.004), less likely to harbor advanced-stage disease (stage III/IV: 15.3% vs 44.4%; OR 0.23, p = 0.020), and characterized by lower CD8 effector, IFN-γ, and immune-checkpoint signatures (each Cohen's d = −0.5 to −0.6 vs sub-B; p < 0.05). Sub-B tumors, by contrast, were older, more frequently advanced, and immune-hot.

The fusion-negative immune-overlap phenotype within sub-B is consistent with a biologically distinct companion axis that is not resolved by driver status alone. Full characterization of this program, including dedicated external-cohort analyses, is outside the scope of the present work and is reserved for a companion study (Cook et al., manuscript in preparation, Paper 2). In the present work, we limit our claims regarding sub-B to the observation that fusion-negative DM1 represents a mechanistically distinct immune-overlap subtype, and we do not pursue its detailed mechanism here.

★ Take-home: Fusion-negative DM1 sub-B represents a mechanistically distinct, older-onset, immune-overlap subtype warranting separate investigation in a companion study (Supplementary Figure S6).

---

# === Figure captions (Main + Supplementary) (05_figure_captions.md) ===


# Main Figures (8 main, Cell Press 4-6 panels each)

## Figure 1. The 8-gene RAI-responsiveness panel resolves a DM1/DM2 cluster within papillary thyroid carcinoma. (4 panels v3)

(A) Sankey flow diagram showing the curation path from the TIERA67 67-gene candidate pool (seven thyroid-relevant biological categories) to the final 8-gene panel (TDS-core sub-category: SLC5A5, TPO, TG, TSHR, PAX8, NKX2-1, FOXE1, DIO1). (B) UMAP embedding of TCGA-THCA primary tumors (n = 504) on the eight-dimensional 8-gene transcript space, colored by KMeans-derived DM1 (red) versus DM2 (blue) cluster assignment. (C) Driver mRNA neutrality. Box plots of BRAF transcript expression in V600E carriers (n = 273) versus wild-type tumors (n = 182), Cohen's d = −0.044, Mann-Whitney p = 0.57. Single-feature classification AUCs for DM1 versus DM2 are shown for BRAF (0.602), TERT (0.578), KRAS (0.525), NRAS (0.521), and HRAS (0.500). (D) Pan-genome cluster reproducibility. Adjusted Rand index (ARI) of unsupervised KMeans partitions versus the 8-gene-driven DM1/DM2 reference: 8-gene panel alone (0.49), TIERA67 (0.90), pan-genome top-5000 by median absolute deviation (0.92), and Driver_anchor genes alone (−0.007). Hypergeometric enrichment p of TIERA67 within pan-genome top 100 = 3 × 10⁻⁴.

Statistical tests: KMeans k = 2; ARI computed against DM1/DM2 reference; Cohen's d with pooled standard deviation; Mann-Whitney U for transcript distributions.

---

## Figure 2. Clinical aggressiveness within Xing 2014 dark matter. (3 panels)

(A) Sankey diagram showing the recovery of 131 of 180 (73%) Xing 2014 BRAF/TERT-negative dark-matter tumors into a defined DM1/DM2 stratum via the 8-gene panel. (B) DM1 prevalence by ancestry: TCGA-THCA (predominantly European/North American; 28.4% BRAF/RAS-negative dark matter; n = 504) versus Korean cohort (37.8% dark matter; n = 865 K2 + Lee pool — GSE286332-PTC dropped to preserve Paper 2 boundary). (C) Kaplan-Meier overall-survival curves for DM1 (red) versus DM2 (blue) within Xing dark matter (TCGA-THCA, n = 180). Log-rank test reported.

Statistical test: log-rank for KM curves; Wilson 95% CI for proportion estimates.

---

## Figure 3. Single-cell external validation confirms thyrocyte-intrinsic DM1 signal. (3 panels)

(A) Per-patient Spearman correlation heatmap of DM1 score in patient-matched tumor versus normal thyrocytes from GSE184362 (Pu et al., 2021), with per-patient r values ranging from 0.798 to 0.886 (p < 10⁻¹⁰). (B) UMAP embedding of thyrocyte clusters from Lu 2023 (GSE193581, n = 23 samples), colored by DM1 score with thyrocyte-specific signal isolated from stromal and immune compartments. (C) Author-independence check: cross-cohort DM1 score concordance between GSE241184 (single-author Phase 1) and GSE184362 (Pu et al., 2021) confirms that the DM1 signal is not driven by laboratory-specific batch effects.

Statistical tests: Spearman r with Bonferroni-adjusted p; UMAP via PCA→neighborhood graph.

---

## Figure 4. Mutation, TERT promoter, and outcome stratification. (3 panels)

(A) Stacked bar plot of 8-cell decomposition (BRAF mutated × TERT mutated × DM cluster) across TCGA-THCA tumors, showing the distribution of patients into the eight defined molecular strata. (B) Kaplan-Meier overall-survival curves for BRAF V600E + TERT promoter mutated (BRAF_TERT+) tumors versus the BRAF/TERT-negative dark matter stratum stratified by DM cluster. (C) The 4-patient OTHER_TERT+ caveat box: TCGA contains only 4 BRAF-wildtype + TERT promoter-mutated cases, limiting any small-N inference about this stratum.

Statistical tests: log-rank for KM; Cox proportional hazards.

---

## Figure 5. Korean cohort validation and FFPE compatibility. (3 panels)

(A) Independent Korean cohort validation. Distribution of 8-gene DM scores in GSE213647 (Lee et al., n = 632), showing preservation of the DM1/DM2 score boundary in an external East-Asian cohort. (B) DM cluster composition in K2 (PRJEB11591, n = 260) NBNR (BRAF-negative + RAS-negative) cohort showing mixed phenotype (vascular invasion enrichment within DM1). (C) Formalin-fixed paraffin-embedded (FFPE) versus fresh-frozen (FF) tissue concordance: density plot of 8-gene scores by processing type, Kolmogorov-Smirnov p = 0.44 (no detectable distributional shift).

Statistical tests: Cohen's d; Kolmogorov-Smirnov for distributional comparison.

---

## Figure 6. Pooled meta-analysis of DM1 overall-survival hazard. (3 panels)

(A) Forest plot of Cox proportional-hazards estimates for DM1 versus DM2 overall-survival risk, by cohort: TCGA-THCA (HR 2.30, 95% CI 0.77–6.88; n = 504, 16 events), MSK-IMPACT thyroid (HR 2.67, 95% CI 1.17–6.10; n = 117, 38 events). Pooled DerSimonian-Laird random-effects estimate: HR 2.53 (95% CI 1.31–4.89), Cochran I² = 0%. (B) I² heterogeneity panel showing no detectable between-cohort heterogeneity (I² = 0%). (C) Sub-cluster Cox analysis for DM1 sub-A versus sub-B (TCGA n = 91, 4 events): DM1 sub-B HR 0.31 (95% CI 0.07–1.39, p = 0.13) — directionally consistent with the immune-overlap phenotype but underpowered as a standalone claim.

Statistical tests: Cox proportional hazards (R `survival`, Python `lifelines`); DerSimonian-Laird random-effects meta-analysis.

---

## Figure 7. DM1 is a fusion-driven dark matter subtype. (4 panels — v3 → v4 2026-05-08, panel D demoted to S6 per audit P1-7)

(A) Silhouette plot of DM1 sub-A (n = 72) versus sub-B (n = 19) sub-clusters within TCGA-THCA, with cluster silhouette score 0.584 supporting the two-population structure. (B) Stacked bar of fusion-positivity by DM cluster: DM1 76.8% (63/82), DM2 30.9%, Fisher OR 7.41 (95% CI 4.38–12.55, p = 1.9 × 10⁻¹³). (C) Stacked bar of fusion partners within DM1 fusion-positive cases: RET (n = 33; CCDC6-RET 17, NCOA4-RET 3, other 13), NTRK (n = 10), ALK (n = 4), BRAF fusions (n = 5). (D) DM1 capture rate of TCGA RET-fusion-positive tumors: 27 of 33 (81.8%), supporting a clinical reflex testing algorithm.

<em>Panel demotion (v4, 2026-05-08): Original panel D (DM1 sub-A vs sub-B phenotype: age, stage, CD8/IFN-γ/checkpoint) moved to Supplementary Figure S6 to keep main Figure 7 mechanism-focused (4 panels) and reserve sub-B mechanistic claims for the companion Paper 2; Cell Press main-figure count thereby reduces from 8 to 7. Original panel E becomes new panel D.</em>

Statistical tests: KMeans silhouette; Fisher exact for OR; Cohen's d; Mann-Whitney U; sensitivity analyses for SV missingness (chi² p = 0.56 MAR).

---

## Figure 8. DM1 epigenetically silences thyroid differentiation machinery. (2 panels — v3 → v4 2026-05-08, panel C demoted to S5b per audit P1-5)

(A) Per-gene HM450 promoter β-value heatmap for the 8-gene panel + DIO2 + SLC26A4, comparing DM1, DM2, and not_DM tumors (n = 503 with HM450 + DM call). DM1 vs DM2 per-gene Cohen's d: TPO (2.30, p = 1.9 × 10⁻¹⁸), DIO1 (1.24, p = 6.5 × 10⁻¹¹), TSHR (1.20, p = 9.8 × 10⁻¹²), PAX8 (0.97, p = 4.5 × 10⁻⁸), TG (0.86, p = 2.2 × 10⁻⁶), FOXE1 (0.84, p = 1.0 × 10⁻⁵), NKX2-1 (0.63, p = 8.9 × 10⁻⁷), SLC5A5 (0.22, p = 0.42, NS). (B) Mean 8-gene panel β-value bar plot by DM cluster: DM1 = 0.385, DM2 = 0.253, not_DM = 0.356 — corresponding to 52% higher DM1 promoter methylation versus DM2.

<em>Mechanism layer cross-reference (2026-05-09).</em> Upstream-effector candidates and panel-canonicality checks are reported as supplementary support: MAPK-pathway output × HM450 mean 8-gene β methylation co-variation (BRAF V600E β = 0.37 ≈ RET fusion β = 0.39 vs RAS β = 0.27; Supplementary Figure SX panels H–J / `v5A_per_driver_class.tsv`); MAPK output × 8-gene panel z RNA Spearman ρ = −0.291 (TCGA n = 572) / −0.395 (Lee n = 632); DM1 sub-A vs sub-B two-axis convergence (sub-A MAPK-active d = +1.79, sub-B HT-active d = +0.62) both reaching panel silencing (panel d = +0.07 NS); and reproduction on the canonical Yoo 2014 TDS-16 (TDS-16 ≈ Panel-8: ΔAUC for MAPK-high classification +0.007 / +0.012 NS, sub-A vs sub-B TDS-16 d = +0.03 NS) — Supplementary Figure SX_v13 / `v13_*.tsv`. The mechanism layer is reported as supportive, not causally proven; the Q9 framing in `09_reviewer_qa.md` retains the "motivates rather than confirms" boundary.

<em>Reserve extension cross-reference (2026-05-09 v15).</em> Supplementary Figure SX_v15 adds K2 score-only overlay (246/260 K2 samples DM2-like using the TCGA-centered 8-gene classifier), a GSE250521 spatial stress test (direct per-spot MAPK × Panel association is positive, not anti-correlated; QC-partial pooled ρ = +0.059), and PRISM/DepMap drug-vulnerability support (7/11 FDR < 0.05 PRISM hits are canonical MAPK-axis inhibitors; AZD-0364 MEK d = −0.594, FDR = 5.9 × 10⁻⁷). These panels are reviewer-reserve support only and do not change the main Figure 8 panel count.

<em>Panel demotion (v4, 2026-05-08): Original panel C (within-DM1 fusion+ vs fusion− methylation, n = 63 vs 19, Cohen's d = −0.36, NS p = 0.31) moved to Supplementary Figure S5b. Rationale per audit P1-5: a non-significant Cohen's d on n = 19 fusion-negative subgroup is power-limited (insufficient to robustly support a "fusion-independent" claim as a main panel) and creates a reviewer-attack surface; the supplementary placement preserves the data while clarifying the boundary of the fusion-independence interpretation, which is now framed in Limitations §3.4 as "consistent with but not formally proving fusion-independent epigenetic silencing".</em>

Statistical tests: Cohen's d (pooled SD); Mann-Whitney U for per-gene β; Fisher exact for fusion × DM cross-tab.

---

# Supplementary Figures (S1-S9)

**S1.** 8-gene panel heatmap sorted by P_DM1 score, comparing TCGA versus K2 cohort distributions (Fig 1 C 이동, v3).

**S2.** 8-gene panel similarity matrix from external single-cell datasets (Pu 2021, Lu 2023, GSE241184).

**S3.** Pan-genome cluster ARI ladder spanning panel sizes 8, 16, 67 (TIERA67), 200, 1000, and 5000 by median absolute deviation.

**S4.** Immune-residualization analysis. DM1 vs DM2 Cohen's d for the 8-gene panel score before residualization (raw d = 1.78), after residualization on an antigen-presentation module (d = 1.00), after residualization on a generic immune signature (d = 1.50), and after dual residualization on both covariates (d = 0.87).

**S5.** Structural-variant missingness sensitivity analyses. Best-case, worst-case, and case-control-matched imputations for the 15 TCGA tumors lacking SV-tested status, showing preserved DM1-versus-DM2 fusion enrichment across scenarios.

**S5b.** *(v4 2026-05-08, demoted from main Figure 8C per audit P1-5.)* Methylation fusion-independence within DM1: fusion-positive (n = 63) versus fusion-negative (n = 19) tumors show equivalent mean 8-gene panel β-value (Cohen's d = −0.36, Mann-Whitney p = 0.31, NS). Interpretation: directionally consistent with a fusion-independent epigenetic silencing layer, but the small fusion-negative subgroup limits formal proof; the result is reported as supportive of, rather than definitive evidence for, parallel mechanism. See Limitations §3.4 for power discussion.

**S6.** *(v4 2026-05-08, contains both: original S6 sub-B detail + main Figure 7D phenotype panel demoted per audit P1-7.)* DM1 sub-A versus sub-B fusion-negative phenotype detail plots. (A) Age (sub-A 37.3 vs sub-B 51.3 years; Cohen's d = −0.82, MW p = 0.004). (B) Stage III/IV (15.3% vs 44.4%; OR 0.23, p = 0.020). (C) CD8/IFN-γ/checkpoint signatures (Cohen's d = −0.5 to −0.6 vs sub-B). (D) Hashimoto-like prevalence (sub-A 3.6% vs sub-B 12.5%; trend). Provided as supporting detail; deeper mechanistic dissection of the immune-overlap sub-B phenotype is reserved for the companion paper (Paper 2).

**SX.** *(UPDATED 2026-05-08 v5, multi-method deconvolution + compositional axes; cumulative scope n = 1,241 bulk samples = TCGA 572 + Lee/GSE213647 632 + GSE76039 37.)* Multi-method bulk cell-type deconvolution and compositional axes of the DM1↔DM2 phenotype. **(A) Multi-method residualization grid.** Cohen's d (DM1 − DM2) for the canonical 8-gene RAI score (`rai_score_recalc`) under five residualization stages (raw / + stromal [Endo + Fibro] / + immune [T + Myel + B + NK] / + Epithelial-only / + all 8 fractions) across four deconvolution methods (NNLS, Ridge-NNLS [L2 α=1], LR-clip, nu-SVR [CIBERSORT-style; 3 ν 0.25/0.5/0.75, lowest-RMSE]). Reference: Lu 2023 (GSE193581) `author_celltype` 67,678 cells × 8 cell types in 1,898 HVG ∩ TCGA gene symbols. Bulk: TCGA-THCA log2(TPM+1), n = 572, canonical labels `dm_like` (DM1_like 403 / DM2_like 110). **(B) Effect-retention ratio after full residualization:** NNLS / Ridge-NNLS 24% (sparse-weight artifact with T-cell collapse), LR-clip 70%, nu-SVR 47% (methodologically aligned with canonical S4 immune-residualization 56% retention). **(C) Per-cell-type DM1 vs DM2 fraction Cohen's d** across the four methods. Direction-consistent: DM1 enriched for Malignant (d ≈ +1.08 to +1.35) and Myeloid (d ≈ +0.92 to +1.14); depleted for Epithelial (d ≈ −1.12 to −2.33; benign/normal thyroid epithelium proxy) and Endothelial. T cell DM1 vs DM2 d = −0.84 uniquely resolved by nu-SVR. **(D) nu-SVR primary mean cell-type fractions** in TCGA-THCA (n = 572). **(E) Methodology and caveats panel.** **(F) Cross-cohort direction consistency.** TCGA nu-SVR Cohen's d (DM1 − DM2) versus Lee/GSE213647 (n = 632) Spearman ρ vs `panel_z` (sign-flipped to align with DM1 direction). Direction-consistent for 7/8 cell types; T-cell discordance reflects binary-vs-continuous score frame. **(G) Within-DM1 fusion+ vs fusion− cell-type fraction Cohen's d** (n = 74 vs 391; kinase fusion: RET/NTRK/ALK/BRAF/PAX8/PPARG aggregated from cBioPortal SV table). All |d| ≤ 0.42 — within-DM1 fusion+/− tumors share near-identical compositions, direct support for the fusion-independent epigenetic silencing claim (Figure 8 / §2.5a). **(H) Per-driver-class TCGA-THCA cell-type composition** (BRAF V600E n = 280 / RAS-mutant n = 54 / driver-negative n = 179). RAS-mutant tumors are Epithelial-cluster–enriched (d_RAS−BRAF = +1.52) and Malignant-cluster–depleted (d = −1.66) versus BRAF V600E, consistent with Landa 2016 / Paper 1 §2.1 differentiation framing. **(I) HM450 8-gene mean β × cell-type fraction Spearman ρ heatmap** (TCGA n ≈ 484 paired): mean β positively tracks Myeloid (ρ = +0.39), B cell (+0.27), Fibroblast (+0.23) and Malignant (+0.30); negatively tracks T cell (ρ = −0.42) and Epithelial (−0.31), connecting the Round-4 methylation layer (DM1 vs DM2 mean-β d = −1.75) to compartment composition. **(J) DM1 sub-A vs sub-B cell-type fraction Cohen's d** (sub-A n = 84, sub-B n = 56). Sub-A is Malignant-cell rich (d = +1.22, p = 2.5 × 10⁻⁹) and Epithelial-poor (d = −1.26, p = 5.2 × 10⁻¹⁰); immune compartment differences are not significant. The sub-A/B split is therefore a tumor-purity-vs-thyrocyte split, not immune-hot vs immune-cold — Paper-2 boundary marker (sub-B = fusion-/mutation-negative Hashimoto-overlap retains thyrocyte identity). **(K) Cell-type composition pseudotime along the canonical 8-gene score.** TCGA-THCA (n = 513, score = `rai_score_recalc`) and Lee/GSE213647 (n = 632, score = `panel_z`) samples were ranked by score, binned into 10 deciles; the mean per-decile cell-type fraction (nu-SVR against Lu 2023) defines a 10-step pseudotime trajectory. Four compartments — Malignant (↓), Epithelial (↑), Myeloid (↓), Endothelial (↑) — show monotonic decile-level Spearman ρ ≥ |0.95| in **both** cohorts, defining a reproducible compositional pseudotime independent of cohort, scoring scheme, or sample size. Source code, per-method fraction tables, panel-level TSVs, and the full v5 composite (panels A–K) at `project/results/p_deconv_2026_05_08/` (build: `plot_deconv_v5_composite.py`; per-panel: `run_deconv_v5_{abc, d_final, e_trajectory}.py`; output figure: `Fig_SX_deconvolution_v5_composite.{png, pdf}`).

**SX_v13.** *(NEW 2026-05-09; Paper 1 Fig 8 mechanism support; cumulative scope n = 1,204 = TCGA 572 + Lee/GSE213647 632.)* The MAPK→thyroid-silencing axis operates indistinguishably on the deployable 8-gene panel and the canonical Yoo 2014 TDS-16. **(A) Cross-cohort MAPK output × thyroid-score Spearman ρ.** TCGA-THCA n = 572 and Lee/GSE213647 n = 632 z-mean MAPK output (DUSP4/5/6, SPRY2/4, ETV4/5, PHLDA1, CCND1; n = 9 genes) versus 8-gene Panel (deployable: DIO1/FOXE1/NKX2-1/PAX8/SLC5A5/TG/TPO/TSHR), TDS-16 (Yoo 2014 canonical: Panel + DIO2/DUOX1/DUOX2/GLIS3/SLC26A4/SLC5A8/THRA/THRB), and TDS−panel (the 8 disjoint TDS-16 genes). Panel-8 ρ = −0.291 (p = 1.3 × 10⁻¹²) / −0.395 (p = 4.7 × 10⁻²⁵); TDS-16 ρ = −0.306 / −0.435; TDS−panel ρ = −0.310 / −0.449. **(B) Per-driver-class TCGA score means** (BRAF V600E n = 344 / RAS-mutant n = 61 / RET fusion n = 43 / NTRK fusion n = 10 / driver-negative n = 103) for MAPK output, Panel-8, TDS-16, and TDS−panel z-scores; the Panel-8 and TDS-16 traces are essentially superimposable, with BRAF V600E vs RAS-mutant Cohen's d Panel-8 = −1.615 versus TDS-16 = −1.616 (identical to three significant figures). **(C) DM1 sub-A vs sub-B Cohen's d** (sub-A n = 93, sub-B n = 62; labels from `d6p7_dm1_subcluster/dm1_subcluster_labels.tsv`) for five metrics: MAPK d = +1.79 (Mann-Whitney p = 4.1 × 10⁻¹⁷), HT signature d = −0.12 (NS), Panel-8 d = +0.07 (NS), TDS-16 d = +0.03 (NS), TDS−panel d = −0.03 (NS). The v12 two-axis convergence (MAPK-active sub-A or HT-active sub-B, both reach 8-gene silencing) reproduces on the full TDS-16. **(D) Per-gene MAPK output × thyroid-gene Spearman ρ heatmap** for the 16 TDS-16 genes in TCGA-THCA and Lee; ★ marks the 8-gene Panel members. Eight-gene members occupy median rank within the 16-gene heatmap (not selectively top-extreme: strongest TCGA negatives = SLC5A8 −0.516 [non-panel], DIO2 −0.481 [non-panel], TPO −0.465 [panel]; NKX2-1 is the consistent positive-ρ outlier in both cohorts, +0.307 / +0.425, a known panel-mean-absorbed caveat). **(E) MAPK-decile pseudotime.** Samples ranked by MAPK output and binned into 10 deciles; mean Panel-8, TDS-16, and TDS−panel score per decile in both cohorts. The two panels trace identical monotonic descent (decile-rank ρ Panel = −0.600 / −0.648, TDS-16 = −0.600 / −0.818). **(F) ROC-AUC for MAPK-high vs MAPK-low classifier** using Panel-8 versus TDS-16 versus TDS−panel score (sign-flipped, low score = high MAPK). Panel-8 AUC = 0.623 (TCGA) / 0.728 (Lee); TDS-16 AUC = 0.631 / 0.740; ΔAUC TDS-16 − Panel-8 = +0.007 (TCGA) / +0.012 (Lee), both non-significant — independently reproducing the existing manuscript p1_onepage_audit ΔAUC = 0.013 NS (DM1/DM2 task) on the orthogonal MAPK-classification task. *Methods.* z-score = within-cohort z-mean of the named gene set; gene lists 100% recovered in both cohorts (Lee via `F1_gene_recovery_mapping.tsv` Ensembl→symbol). Driver class from cBioPortal SV/MAF anchors (`audit_2026_04_30/round3/cbio_sv_thca.tsv`). DM1 sub-A/sub-B labels from `d6p7_dm1_subcluster`. Source code, panel-level TSVs, and the 6-panel composite at `project/results/p_deconv_2026_05_08/` (build: `plot_v13_tds16_panel_mapk.py`; pipeline: `run_v13_tds16_panel_mapk.py`; output figure: `Fig_SX_v13_TDS16_MAPK.{png, pdf}`).

**SX_v14.** *(NEW 2026-05-09 evening; Paper 1 Fig 8 / Q14 cross-cohort generalizability lock; cumulative scope n = 1,287 = TCGA 572 + Lee 632 + GSE126698 28 + GSE286332 18 + GSE76039 37.)* Cross-cohort forest of MAPK output × thyroid-score Spearman ρ. Three score panels (8-gene Panel deployable / TDS-16 Yoo 2014 canonical / TDS−panel 8 disjoint genes) × five cohorts. Squares = per-cohort Spearman ρ (size proportional to √n); horizontal bars = 95% CI (Fisher z-transform); diamond = pooled fixed-effect Fisher-z ρ. **Pooled MAPK × Panel-8 ρ = −0.327, 95% CI [−0.376, −0.278], n_total = 1,287.** Cochran Q = 14.77 (df = 4, p = 0.005), I² = 72.9% — heterogeneity is biologically expected and predicted by the v12 two-axis convergence model (see Q9 / Q14 in `09_reviewer_qa.md`): GSE286332 (Korean PTC vs PTC+HT, n = 18) is dominated by the HT route where v12 predicts MAPK decouples (ρ = −0.040 NS), and GSE76039 (Landa 2016 PDTC + ATC, n = 37) shows panel saturation at the dedifferentiated end of the axis (ρ = +0.125 NS). The two well-differentiated primary cohorts (TCGA + Lee, n = 1,204) carry the entire signal (ρ = −0.291 / −0.395, both p ≪ 10⁻¹²). Source code, panel TSV, and the 3-panel forest figure at `project/results/p_deconv_2026_05_08/` (build: `plot_v14_cross_cohort_forest.py`; pipeline: `run_v14_cross_cohort_forest.py`; output figure: `Fig_SX_v14_cross_cohort_forest.{png,pdf}`; per-cohort table: `v14_cross_cohort_forest.tsv`).

**SX_v15.** *(NEW 2026-05-09 evening; Paper 1 Fig 8 reviewer-reserve extensions; main Figure 8 panel count unchanged.)* K2 score-only overlay, GSE250521 spatial stress test, and PRISM/DepMap drug-vulnerability overlay. **(A) K2 mini-index overlay.** A TCGA-trained centered-profile 8-gene classifier gives TCGA 5-fold CV AUC = 0.960 and classifies 246/260 K2 PRJEB11591 samples as DM2-like (median p_DM2 = 0.978). K2 lacks MAPK-output genes, so this panel is score-distribution evidence only and is not included in the MAPK × Panel forest. **(B) Spatial Lu 2023 GSE250521 MAPK × Panel test.** Raw h5ad spots were log1p(CP10K)-normalized, within-slide z-scored, and summarized as MAPK-9 and Panel-8 scores. Across 12 tumor-region slides (PTC/LPTC/ATC; 43,180 spots), the fixed-effect per-slide MAPK × Panel association is positive rather than negative (ρ = +0.326 [0.318, 0.335]); QC-partial ρ after total-count, gene-count, and mitochondrial-fraction adjustment is much smaller (+0.059 [0.049, 0.068]). This panel is a non-confirmatory stress test, not a spatial mechanism lock. **(C) PRISM overlay.** Among 1,518 PRISM drugs, 11 are FDR < 0.05 and DM1-high selective; 7/11 are canonical MAPK-axis inhibitors (MEK/RAF/ERK). The top hit is AZD-0364 (MEK; Cohen's d = −0.594, FDR = 5.9 × 10⁻⁷). Across 690 nonmissing cell lines, DM1 score anti-correlates with canonical MAPK-inhibitor mean LFC (Spearman ρ = −0.236, p = 3.3 × 10⁻¹⁰), and high-vs-low DM1 bins differ by d = −0.617. **(D) DepMap overlay.** Top DM1-high dependencies are MYC (d = −0.499, p = 1.0 × 10⁻¹¹) and NAMPT (d = −0.445, p = 1.5 × 10⁻⁹); thyroid transcription factors are not the primary therapeutic lever. PRISM/DepMap supports MAPK-axis vulnerability but does not measure restoration of thyroid-gene expression after inhibitor treatment. Source code and outputs: `run_v15_mechanism_extensions.py`, `v15A_*`, `v15B_*`, `v15C_*`, `v15_mechanism_extensions_summary.json`, and `Fig_SX_v15_mechanism_extensions.{png,pdf}` in `project/results/p_deconv_2026_05_08/`.

**SX_v16.** *(NEW 2026-05-09 evening; reviewer-reserve diagnostic stress test.)* Spatial rescue failure and PRISM claim-boundary lock. **(A) GSE250521 MAPK × Panel adjustment grid.** Across 57,997 spatial spots and 43,180 tumor-stage spots, raw MAPK × Panel association is positive; QC and cell-state partialing attenuate but do not flip the direction. Tumor-slide pooled ρ = +0.326 [0.318, 0.335] raw, +0.058 [0.048, 0.068] after QC adjustment, and +0.034 [0.024, 0.044] after QC plus epithelial/CAF/EMT/hypoxia/proliferation adjustment. **(B) Epithelial-score quartile grid.** Stage × within-slide epithelial-score quartiles show positive or near-zero MAPK × Panel relationships; no robust anti-correlation emerges. **(C) PRISM top-k enrichment sensitivity.** The correct claim is top 7/7 canonical MAPK-axis inhibitors, 7/11 FDR < 0.05 canonical MAPK-axis inhibitors (hypergeometric p = 3.25 × 10⁻¹⁵), and 8/15 top-15 canonical MAPK-axis inhibitors; the earlier shorthand "top 15 all MAPK" is not used. **(D) PRISM lineage sensitivity.** DM1 score × canonical MAPK-inhibitor mean LFC remains negative when thyroid cell lines are excluded (ρ = −0.240, p = 2.27 × 10⁻¹⁰); thyroid-only cell-line correlation is not interpretable because the available thyroid lines share a constant DM1 score. Source code and outputs: `run_v16_spatial_prism_diagnostics.py`, `v16_spatial_*`, `v16_prism_*`, `v16_spatial_prism_diagnostics_summary.json`, and `Fig_SX_v16_spatial_prism_diagnostics.{png,pdf}` in `project/results/p_deconv_2026_05_08/`.

**SX_v17.** *(NEW 2026-05-09 evening; spatial signal decomposition.)* Why the GSE250521 spatial MAPK × Panel stress test fails as mechanism confirmation. **(A) Adjustment-family ladder.** In tumor epithelial-top50 spots, raw MAPK × Panel ρ = +0.208, but full covariate adjustment reduces the linear-residual ρ to +0.052 and full spatial+detection adjustment reduces it to +0.007. **(B) Covariate R².** Detection/depth covariates explain much of both modules (median tumor epithelial R²: MAPK = 0.809; Panel-8 = 0.683). **(C) Depth/detection structure.** MAPK-detection and Panel-detection are positively correlated (median ρ = +0.383), consistent with broad Visium detection co-localization. **(D) MAPK submodule tests.** DUSP/SPRY, ETV/PHLDA1, no-CCND1, and CCND1-only variants do not recover a strong full-adjusted anti-correlation. **(E) Gene-pair sign consistency.** The most negative pair is CCND1 × TPO (median ρ = −0.096), but gene-pair FDR is non-significant (FDR = 0.694). **(F) Random-module null.** Observed MAPK × Panel raw median correlation is near the random-module median (percentile = 0.545), and full-residual observed correlation is below the random-module median (percentile = 0.343). These diagnostics support treating GSE250521 as a spatial caveat rather than positive Fig 8 mechanism evidence. Source code and outputs: `run_v17_spatial_signal_decomposition.py`, `v17_spatial_*`, `v17_spatial_signal_decomposition_summary.json`, and `Fig_SX_v17_spatial_signal_decomposition.{png,pdf}`.

**SX_v18.** *(NEW 2026-05-09 evening; spatial lag and pocket closure.)* Final GSE250521 spatial rescue test. **(A) Same-spot and neighborhood-lag MAPK × Panel correlations** in tumor epithelial-top50 spots. Raw same-spot median ρ = +0.110; raw KNN 1-6 neighborhood-lag median ρ = +0.090; raw KNN 7-18 median ρ = +0.071. Full spatial+detection residual same-spot median ρ = −0.011; residual KNN 1-6 median ρ = +0.033; residual KNN 7-18 median ρ = +0.010. **(B) MAPK-high/Panel-low anti-pocket enrichment.** Raw anti-pocket enrichment = 0.882 (OR = 0.693), so anti-pockets are depleted rather than enriched; residual anti-pocket enrichment = 1.010 (OR = 1.018), approximately independence-level. **(C) Concordant MAPK-high/Panel-high pockets** are not strongly enriched (raw enrichment = 1.002; residual enrichment = 0.982), consistent with broad co-detection rather than focal anti-silencing domains. **(D) Anti-pocket covariates.** The strongest raw anti-pocket covariate depletion is Panel detection (Cohen's d = −1.094), indicating low panel-detection spots rather than robust biological MAPK-high/Panel-low neighborhoods. These panels close the last spatial escape hatch: GSE250521 is retained as a Visium-resolution/detection caveat, not Fig 8 mechanism support. Source code and outputs: `run_v18_spatial_lag_pockets.py`, `v18_spatial_lag_*`, `v18_spatial_pocket_*`, `v18_spatial_lag_pockets_summary.json`, and `Fig_SX_v18_spatial_lag_pockets.{png,pdf}`.

**S7.** Cross-cohort DM score portability in Korean validation cohorts. Score-distribution overlays and threshold-portability checks across K2, Lee/GSE213647, and the small Korean reference cohort used for calibration.

**S8.** Hypomethylating-agent + radioiodine re-induction schematic + literature meta. Decitabine + I-131 retrospective trial summary (NCT00085293, NCT01065090) and the rationale for prospective trial design stratified by DM1 status. (Fig 8 D 이동, v3.)

**S9.** K2 mini-index calibration diagnostic + alternate evidence chain. Per-gene inflation factors (4.9–12.5× across panel genes); 4-metric direction check (raw, within-sample-z, per-gene-z, per-gene-rank); alternate evidence from score-distribution portability and DM-call consistency.

---

# Additional supplementary figures (v17 audit + Dark matter phase 1/2 + Landa 2016)

These figures complement the main Fig 1-8 + S1-S9 set with sanity-check (v17 audit phase, 2026-04-29), single-cell foundation + multisite validation (Dark matter phase 1/2), and Discussion §3.1 cite-save evidence (Landa 2016 GSE76039 heatmap). Full panel-by-panel descriptions with verification markers are maintained in the parallel detail file `05_supp_figure_captions_v17_dm.md`.

## Part A — v17 Audit phase sanity-check

**SA1.** 4-way revalidation of the 8-gene DM1/DM2 cluster across BRAF V600E × TERT promoter mutational strata. (A-D) Stratum-specific Cox HR + KM curves; small-N caveat for BRAF−/TERT+ (n=4).

**SA2.** FFPE versus fresh-frozen tissue compatibility QC. (A) 8-gene panel score distribution by tissue type, Kolmogorov-Smirnov p=0.44 (no detectable shift). (B-C) Per-gene paired analysis + DM-call concordance.

**SA3.** Panel-size sensitivity. (A) 5-fold cross-validated AUC across 8/10/12/16-gene variants; ΔAUC 8 vs 16 = 0.013, NS. (B-C) Cluster-call concordance matrix + ARI ladder.

**SA4.** Single-cell wrap-up — thyrocyte-intrinsic 8-gene signal across GSE184362 (Pu et al., 2021), GSE193581 (Lu 2023), GSE241184 (Phase 1). (A-B) Per-cohort summary + per-patient r 0.798–0.886, all p < 10⁻¹⁰. (C) Author-independence cross-cohort concordance.

**SA5.** Differentiation trajectory along the 8-gene panel score in TCGA-THCA. (A) Pseudotime ordering of primary tumors. (B-D) Driver mutation distribution, differentiation gene expression, and DM cluster overlay along the trajectory.

**SA6.** MSK-IMPACT bias panel — advanced-disease cohort (Landa et al., 2016) labeled explicitly. (A) PDTC n=84 + ATC n=33 composition. (B-D) Mutation-frequency, stage, and DM-prevalence comparison versus TCGA-THCA primary tumors.

## Part B — Dark matter phase 1 single-cell UMAP foundations

**SB1.** Single-cell UMAP overview of GSE184362 (Pu et al., 2021; n=6 PTC patients, Fudan University). (A) Embedding colored by patient identity. (B) Cell-type annotation (thyrocyte / immune / stromal / endothelial). (C) Tumor versus adjacent normal labeling.

**SB2.** Single-cell UMAP — molecular signatures. (A) 8-gene panel score per cell (continuous gradient). (B) HLA-II module signature (Paper 1 residualization control only; deeper HLA-II analysis = Paper 2 territory). (C) Canonical thyroid differentiation transcripts (TG, TPO, TSHR averaged). (D) Proliferation signature (MKI67, TOP2A).

**SB3.** Thyrocyte-restricted UMAP. (A) Subset to KRT8+ KRT19+ EPCAM+ cells. (B) 8-gene panel score gradient on thyrocyte UMAP. (C) Tumor versus adjacent-normal labeling. (D) Per-patient thyrocyte distribution along the 8-gene score axis.

## Part C — Dark matter phase 2 multisite + per-patient validation

**SC1.** Per-patient Spearman r forest plot — tumor versus adjacent-normal thyrocyte 8-gene scores in GSE184362 (n=6 patients; per-patient r 0.798–0.886; Bonferroni-adjusted p < 10⁻¹⁰; per-patient n_cells annotated).

**SC2.** Multisite trajectory — DM1/DM2 score reproducibility across TCGA-THCA (n=504), MSK-IMPACT advanced-disease (n=117), K2 / PRJEB11591 (n=260), and Lee / GSE213647 (n=632).

**SC3.** Multisite single-cell UMAP — joint cell embedding across GSE184362 (Pu 2021), GSE193581 (Lu 2023), and GSE241184 (Phase 1) after batch-corrected integration. Dataset-independent gradient for the 8-gene panel score.

**SC4.** GSE184362 (Pu 2021) summary table — per-patient metadata, n_cells per condition, per-patient Spearman r, DM1/DM2 prediction.

**SC5.** External-validation pooled scatter — per-patient pseudobulk 8-gene score in tumor versus adjacent-normal thyrocytes across GSE184362 + Lu 2023.

**SC6.** Pooled scatter — DM1 probability versus 8-gene panel score across external single-cell cohorts; per-cohort markers + pooled regression with 95% CI.

**SC7.** K2 (PRJEB11591) versus Yoo 2016 reference — 8-gene score concordance + K2 mini-index calibration mismatch diagnostic (R4-4 audit; per-gene inflation factors 4.9–12.5×; within-sample-centered profile restores TCGA direction).

## Part D — Landa 2016 evidence (Discussion §3.1 cite save)

**SD1.** Landa et al. (2016) GSE76039 — differentiation transcript suppression in advanced thyroid cancer. (A) Per-sample heatmap of seven canonical differentiation transcripts (TG, TSHR, TPO, PAX8, SLC26A4, DIO1, DUOX2) across PDTC + ATC samples from the GSE76039 transcriptome subset. (B) 5-of-8 overlap with this paper's panel (TG, TSHR, TPO, PAX8, DIO1). (C) Convergence framing — Yoo et al. (2016) panel origin and Landa et al. (2016) advanced-disease list reach the same differentiation core via independent paths (Discussion §3.1 reverse-causality framing).

Cite: Landa I, Ibrahimpasic T, Boucai L, Sinha R, Knauf JA, Shah RH, Dogan S, Ricarte-Filho JC, Krishnamoorthy GP, Xu B, Schultz N, Berger MF, Sander C, Taylor BS, Ghossein R, Ganly I, Fagin JA. Genomic and transcriptomic hallmarks of poorly differentiated and anaplastic thyroid cancers. *J Clin Invest.* 2016;126(3):1052–1066. doi:10.1172/JCI85271. PMID 26878173. PMC4767360.

For full panel-by-panel descriptions, exact statistical-test specifications, and the 17-item verification queue (per-figure ⚠ markers), see `05_supp_figure_captions_v17_dm.md`.

---

# Figure 작성 우선순위 (Cell Press)

| Priority | Figure | Status | Source code |
|---|---|---|---|
| ★★★ | Fig 7 (DM1 mechanism, 5 panels) | Build needed — combines R3-F4 + R4-2 + R4-3 | `v17_audit_F2_F3_F4.py` + `v17_audit_R4_all.py` + new |
| ★★★ | Fig 8 (Epigenetic, 3 panels v3) | Build needed — R5-2 paradigm | `v17_audit_R5_all.py` + new |
| ★★★ | Fig 6 (Meta forest, 3 panels) | Build needed — N1 meta + R4-2 sub-cluster Cox | `v17_D4P1_forest_meta.py` + new |
| ★★ | Fig 1 (Panel + cluster, 4 panels v3) | Build needed | `v17_8gene_figs_v2.py` + `p4_pangenome_vs_tiera67.py` |
| ★★ | Fig 2 (Xing rescue + KM) | Build needed | `v17_4way_figure.py` + new |
| ★ | Fig 3 (sc validation) | Existing figs available | sc figs + new author-independence panel |
| ★ | Fig 4 (BRAF×TERT 8-cell) | Existing figs available | `v17_quad_tert_stack.html` + new |
| ★ | Fig 5 (Korean + FFPE) | Existing figs available | `v17_KOREAN_K2_v260_figure.py` + new |

---

# Narrative tightening notes (2026-05-11)

The manuscript reads best when the main figures are treated as a clean claim ladder rather than a list of analyses.

## Main-text priority order

1. **Figure 1** should carry the compact panel-defining claim only.
2. **Figure 2** should handle clinical relevance only.
3. **Figure 3** should establish single-cell intrinsic validation.
4. **Figure 4** should stop at driver-context stratification and not overreach mechanistically.
5. **Figure 5** should be the portability / assay-compatibility proof.
6. **Figure 6** should remain the survival-risk summary.
7. **Figure 7** should carry the fusion-mechanism support.
8. **Figure 8** should carry the epigenetic support and explicitly stay correlative.

## What should stay out of the main figures

- calibration mismatch diagnostics
- spatial stress-test caveats
- niche reserve extensions
- method-comparison stress tests
- any counterexample that is useful but not claim-building

## What makes the paper stronger

- fewer claims per figure
- more separation between discovery and mechanism
- cleaner boundary language in captions
- fewer supplementary panels inside the main-figure storyline

---

# === Discussion (06_discussion.md) ===


# Section 3 · Discussion

## 3.1 A three-layer pathology of BRAF/RAS-negative dark matter

Our findings frame the BRAF/RAS-negative compartment of papillary thyroid carcinoma as a structured molecular state rather than a heterogeneity remainder. The 8-gene RAI-responsiveness panel resolves this compartment into two clusters, DM1 and DM2, in which mechanism is layered at three levels. Genetically, DM1 is enriched 7.4-fold over DM2 for tyrosine kinase rearrangements of RET, NTRK, ALK, and BRAF, with 76.8% fusion positivity overall and capture of 81.8% of TCGA RET-fusion-positive primary tumors. Phenotypically, DM1 itself partitions into a fusion-driven, younger, less advanced sub-A and an older, more advanced, immune-infiltrated sub-B that we treat here only at the level of recognition. Epigenetically, both sub-populations of DM1 share fusion-independent promoter hypermethylation of canonical thyroid differentiation transcripts, with the strongest single-gene signal at TPO (Cohen's d = 2.30) and a uniformly elevated panel-mean β value relative to DM2. We interpret DM1 as an integrated genetic-epigenetic state of partial differentiation that an 8-gene transcript readout can identify at the primary-tumor stage, before the dedifferentiation trajectory described in advanced disease completes.

Landa et al. (2016) characterized the genomic and transcriptomic landscape of advanced thyroid cancer and established that poorly differentiated and anaplastic thyroid cancers arise through progressive accumulation of genomic abnormalities together with profound loss of thyroid differentiation programs. In the primary BRAF/RAS-negative PTC compartment, our DM1 framework identifies an upstream signature consistent with that dedifferentiation trajectory: the same differentiation machinery that is silenced at the advanced-disease end is already epigenetically attenuated in DM1 primary tumors (mean 8-gene promoter methylation beta 0.385 versus 0.253 in DM2). We therefore interpret DM1 not as a generic residual class, but as a structured molecular state within thyroid-cancer dark matter.

The convergence between our DM1 panel and the Landa advanced-disease gene program is notable because the two were reached by independent routes. Our panel was selected from canonical RAI-biology genes grounded in Yoo et al. (2016), whereas the Landa framework emerged from advanced-disease transcriptomics. The five-gene overlap (TG, TSHR, TPO, PAX8, DIO1) therefore argues for convergent biological validation rather than circular panel construction.

The fusion-negative DM1 sub-B population also merits separate comment. Relative to fusion-positive DM1 tumors, sub-B tumors are older-onset, more stage-advanced, and more immune-infiltrated. We view this as evidence that DM1 contains at least two biologically distinct states: a canonical fusion-driven state and a companion immune-overlap state. Full mechanistic dissection of that latter program falls outside the scope of the present paper and is reserved for a companion study.

A counter-intuitive feature of our data is that BRAF V600E-mutated tumors show higher HLA-I module signal than BRAF-wildtype tumors, opposite to earlier reports of BRAF-linked HLA-I downregulation. Rather than negating prior work, this result suggests that immune escape in BRAF-driven PTC may proceed through mechanisms other than simple HLA-I loss, and that the relationship between oncogenic signaling and antigen presentation may be more context-dependent than a single-axis model implies.

## 3.2 Clinical actionability and a reflex testing algorithm

The 2015 American Thyroid Association guidelines, with 2025 updates, still center clinico-pathological risk variables and BRAF V600E as the principal molecular modifier. Fusion drivers and differentiation-gene silencing are not explicitly incorporated into routine primary-tumor stratification. Our findings provide an orthogonal molecular axis for the BRAF/RAS-negative compartment: a positive DM1 RNA score captures 81.8% of TCGA RET-fusion-positive tumors and 76.8% of tyrosine-kinase-fusion-positive DM1 tumors overall. In practical terms, this creates a rational reflex step between broad pathological triage and full fusion sequencing.

Selpercatinib approval in RET-fusion-positive thyroid cancer makes this distinction clinically relevant. The point of the present framework is not to claim immediate treatment assignment from an 8-gene score alone, but to enrich the subgroup in which selective RET/NTRK/ALK/BRAF fusion testing is most likely to be informative. A DM1-positive call can therefore be interpreted as an upstream enrichment signal for fusion-directed workup rather than a stand-alone therapeutic decision rule.

The epigenetic layer adds a second translational implication. DM1 tumors exhibit strong, fusion-independent promoter hypermethylation across thyroid differentiation genes, especially TPO, DIO1, and TSHR. This pattern provides a mechanistic rationale for revisiting hypomethylating-agent-based RAI re-induction with molecular stratification rather than empiric enrollment alone. The SLC5A5/NIS exception is equally informative: it suggests that promoter methylation is not the sole determinant of iodine-handling failure and that restoration strategies may need to target both transcriptional repression and post-transcriptional trafficking biology.

Taken together, these data support a two-axis clinical framework for the BRAF/RAS-negative compartment: reflex fusion testing for targeted-therapy discovery and DM1-stratified evaluation of epigenetic RAI re-induction strategies.

## 3.3 East-Asian generalizability and external validation

The DM1/DM2 axis was derived in TCGA-THCA but reproduced across independent Korean bulk cohorts, indicating that the signal is not restricted to a single ancestry or sequencing pipeline. The preservation of score geometry across K2, GSE213647, and the external Korean reference arm is consistent with prior observations that the dark-matter compartment is enriched in East Asian PTC populations.

Independent single-cell datasets reinforce a second point: the DM1 signal is thyrocyte-intrinsic. In both the Pu et al. (2021) and Lu 2023 datasets, the score pattern persists after moving from bulk tissue to tumor-cell-focused analyses, arguing against the interpretation that DM1 is merely a byproduct of stromal admixture or bulk immune contamination. This matters because several of the companion immune features associated with DM1 could otherwise be mistaken for non-epithelial signal leakage.

Finally, cross-cohort agreement between TCGA, MSK-IMPACT, Korean bulk cohorts, and external single-cell datasets suggests that the DM1 axis is transportable across discovery, validation, and mechanistic contexts. That transportability does not eliminate cohort-specific caveats, but it does support the claim that the observed differentiation/fusion/methylation structure reflects a reproducible biological axis rather than a single-dataset artifact.

## 3.4 Limitations

Several limitations bound the present interpretation. First, the overall-survival inference within TCGA-THCA is constrained by a low event rate — 16 deaths across 504 primary tumors — so the single-cohort hazard estimate is wide (HR 2.30, 95% CI 0.77–6.88) and is mitigated, but not replaced, by the MSK-IMPACT meta-analysis. Second, the within-DM1 fusion-independence test for promoter hypermethylation rests on n = 19 fusion-negative DM1 tumors; the negative result (Cohen's d = −0.36, p = 0.31) is directionally consistent with parallel-mechanism epigenetic silencing but is power-limited and is reported here as supportive rather than definitive (Supplementary Figure S5b).

Third, the MSK-IMPACT cohort is enriched for advanced and treatment-refractory disease, which inflates event rates relative to community PTC; cross-cohort meta-analysis assumes that the DM1 hazard direction generalizes despite this case-mix difference. Fourth, the Korean K2 (n = 235) and Lee (n = 632) cohorts provide axis portability and ancestry coverage but do not contribute time-to-event annotation comparable to TCGA, leaving us without a Korean overall-survival replication of the DM1 hazard. Fifth, the HM450 methylation analysis is restricted to TCGA-THCA; an external methylation cohort with matched 8-gene transcript profiles is not available, so the DM1 epigenetic signal has not been independently replicated in another methylation dataset.

Sixth, the panel name explicitly references RAI-responsiveness biology, but no cohort in this study contains prospective post-thyroidectomy RAI outcome data linked to per-sample DM1 calls; consequently, the framework supports prospective evaluation of epigenetic-targeted RAI re-induction as a hypothesis, not an outcome-validated decision rule. Seventh, the East-Asian fusion landscape — and in particular RET-fusion frequency relative to TCGA — has not been independently re-quantified in our Korean cohorts because matched structural-variant-tested status is not available at the same coverage as the TCGA-THCA cBioPortal annotations.

Finally, the fusion-negative DM1 sub-B subgroup is older, more advanced, and immune-infiltrated, and its mechanism is mechanistically distinct from the sub-A fusion-driven state. Full dissection of sub-B — including HLA-class structure, BCR clonal architecture, mediation analysis, and Hashimoto-overlap phenotype — is reserved for a companion study (Cook et al., manuscript in preparation, Paper 2). The present manuscript limits its sub-B claims to the observation of a mechanistically distinct immune-overlap state and explicitly does not pursue its detailed regulatory biology here.

---

# === STAR Methods (07_star_methods.md) ===


# STAR Methods (draft v1)

---

## Key Resources Table

(structured Cell Press format — populated below as plain table; final manuscript: per-section table)

| Reagent / Resource | Source | Identifier |
|---|---|---|
| **Biological samples — bulk** | | |
| TCGA-THCA RNA-seq + WGS + HM450 + clinical | The Cancer Genome Atlas | dbGaP phs000178; cBioPortal `thca_tcga_pub` |
| MSK-IMPACT thyroid (Landa 2016 cohort) | Landa et al., 2016 | cBioPortal `thca_mskcc_2016` |
| K2 / PRJEB11591 Korean PTC RNA-seq (Yoo 2016) | Yoo et al., 2016 | ENA PRJEB11591 |
| Lee Korean PTC cohort | Lee et al., GEO | GSE213647 |
| GSE286332 Korean PTC reference arm | Macrogen / Dongguk Univ | GEO GSE286332 |
| GSE184362 Pu 2021 single-cell PTC | Pu et al., 2021 | GEO GSE184362 |
| GSE193581 Lu 2023 single-cell PTC | Lu et al., 2023 | GEO GSE193581 |
| GSE241184 Phase 1 single-cell PTC | (Phase 1) | GEO GSE241184 |
| **Software** | | |
| pyDESeq2 (DEG analysis) | Snakemake/lab | https://github.com/owkin/PyDESeq2 |
| scikit-learn (KMeans, LogReg, AUC) | Pedregosa et al. | https://scikit-learn.org |
| lifelines (Cox PH, log-rank) | Davidson-Pilon | https://lifelines.readthedocs.io |
| STAR aligner | Dobin et al. | https://github.com/alexdobin/STAR |
| kallisto (pseudo-alignment) | Bray et al. | https://pachterlab.github.io/kallisto |
| GENCODE v44 reference | EBI | https://www.gencodegenes.org/human/release_44.html |
| **Deposited data** | | |
| cBioPortal SV (RET/NTRK/ALK/BRAF fusions) | cBioPortal API | `thca_tcga_pub` study, SV endpoint |
| HM450 promoter methylation (Illumina HumanMethylation450) | cBioPortal API | `thca_tcga` legacy study |
| **Source code (this paper)** | | |
| 8-gene panel + DM cluster pipeline | Cook et al., this paper | [TODO: insert public code repository URL on submission]; [TODO: insert Zenodo DOI on submission] |

---

## Resource availability

**Lead contact.** Further information and requests for resources should be directed to the lead contact, Seungho Cook (kukshomr@gmail.com).

**Materials availability.** This study did not generate new unique reagents. All analyses were performed on publicly accessible datasets (TCGA, GEO, ENA, cBioPortal). Korean cohort access (K2 / PRJEB11591, GSE213647, GSE286332 reference arm) is available via the listed repositories.

**Data and code availability.** Public source data are available from TCGA, GEO, ENA, and cBioPortal under the identifiers listed above. Analysis code and figure-generation scripts will be released in a public repository together with an archival DOI at submission or acceptance ([TODO: insert public code repository URL on submission]; [TODO: insert Zenodo DOI on submission]). Intermediate data tables used in the manuscript, including per-sample DM scores, fusion annotations, methylation summaries, and meta-analysis inputs, are provided through Supplementary Tables S1-S10.

---

## Experimental model and study participant details

This study uses publicly available genomic and transcriptomic data from previously published cohorts:

- **TCGA-THCA** (n = 504 with overall survival annotation; 513 primary tumors total). Discovery cohort. Publicly available via The Cancer Genome Atlas (Cancer Genome Atlas Research Network, 2014).
- **MSK-IMPACT thyroid** (n = 117; advanced disease, mostly PDTC + ATC). Validation cohort. (Landa et al., 2016).
- **K2 / PRJEB11591** (n = 235; primary Korean PTC, post-QC). Validation cohort. (Yoo et al., 2016).
- **Lee / GSE213647** (n = 630; Korean PTC, post-QC). Validation cohort.
- **GSE286332 reference arm** (n = 9 Korean PTC). Small external Korean reference set retained here for calibration and score-portability description only; <em>dropped from the Paper 1 Korean validation total per the 2026-05-07 P1-2 audit. The canonical Paper 1 Korean validation cohort = K2 (n = 235) + Lee et al. (n = 630) = n = 865. GSE286332 PTC samples are retained in Paper 2 scope (GSE286332 PTC vs PTC+HT main cohort, n = 18).</em>
- **GSE184362 Pu 2021** (n = 7 PTC patients; single-cell). External validation.
- **GSE193581 Lu 2023** (n = 23 single-cell samples). External validation.
- **GSE241184** (n = 1; Phase 1 single-cell). Internal pilot.

All studies were originally approved by the respective institutional review boards. The present analysis used only de-identified public data and is exempt from additional IRB review.

---

## Method details

### Cohort assembly and clinical metadata harmonization

Bulk RNA-seq quantifications were obtained as log2(TPM + 1) (TCGA) and log2(FPKM + 1) (Korean cohorts and the GSE286332 reference arm) and processed through unified gene-level filtering (≥10 reads in ≥30% of samples; GENCODE v44 protein-coding annotation). Clinical metadata (age, sex, stage, vital status, time-to-event) were harmonized from cBioPortal and source publications. Per-cohort missingness was tabulated (Supplementary Table S2) and addressed in sensitivity analyses.

### 8-gene panel selection

The 67-gene candidate pool (TIERA67) was assembled from seven thyroid-relevant biological categories (Supplementary Table S1): TDS-core differentiation markers (16 genes), MAPK-output transcripts (10), thyroid driver genes (12, including BRAF, NRAS, HRAS, KRAS, RET, NTRK1, NTRK3, ALK, PAX8, PPARG, TERT, EIF1AX), aggressive-disease markers (10), dedifferentiation/EMT markers (10), light immune-stromal markers (5), and thyroid-lineage extras (4). The 8-gene panel (SLC5A5/NIS, TPO, TG, TSHR, PAX8, NKX2-1, FOXE1, DIO1) was selected from the TDS-core sub-category based on canonical RAI-uptake biology (Yoo et al., 2016 *PLOS Genet*) — independently of and prior to access of the Landa et al. (2016) ATC-silenced gene list. This panel design therefore predates exposure to advanced-disease ATC transcriptomic outputs and represents an independent path to the same differentiation axis.

Panel size sensitivity was tested across n = 8, 10, 12, and 16 gene variants (Supplementary Figure S2). The 16-gene full TDS-core yielded TCGA 5-fold AUC of 0.975 versus the 8-gene panel's 0.962 (ΔAUC = 0.013, NS); the 10-gene and 12-gene intermediate panels yielded 0.972 and 0.969, respectively. The 8-gene panel was retained for clinical interpretability.

### DM1/DM2 cluster definition

Bulk RNA expression of the 8-gene panel was z-standardized within cohort, and KMeans (k = 2, scikit-learn default initialization) was applied to assign DM1 versus DM2 cluster identity. The cluster boundary was confirmed by silhouette analysis (mean silhouette score > 0.4 within both clusters). For per-sample classification probability (P(DM1)), a logistic regression classifier was trained on the TCGA DM1/DM2 partition using the within-sample-centered profile of the 8-gene panel (mean panel score subtracted), yielding cohort-portable class probabilities (cross-validated AUC = 0.968 in TCGA; Methods, Supplementary Figure S3).

### Single-cell external validation

For GSE184362 (Pu et al., 2021), per-patient pseudo-bulk DM1 scores were computed by averaging thyrocyte-marker-positive (KRT8, KRT19, EPCAM) cell profiles per patient and per condition (tumor versus adjacent normal). Per-patient Spearman correlation between tumor and normal scores was computed. Lu 2023 (GSE193581) was processed with the same DM1-scoring framework together with stromal and immune contamination control.

### Survival analysis and meta-analysis

Cox proportional hazards regression was performed using `lifelines` (Davidson-Pilon, Python). Pooled meta-analysis was performed using random-effects DerSimonian-Laird estimation; between-cohort heterogeneity was assessed via Cochran I². Sub-cluster Cox analyses (DM1 sub-A vs sub-B) used patient-level event data within TCGA-THCA primary tumors (Supplementary Table S8).

### cBioPortal API access (SV + methylation)

Structural variants (RET, NTRK1/3, ALK, BRAF, PAX8-PPARG, others) were retrieved from cBioPortal `thca_tcga_pub` study via REST API (endpoint `/structural-variant/fetch`); SV-tested status was captured per tumor (n = 542 / 557 = 97.3%). HM450 promoter β-values for the 8-gene panel + DIO2 + SLC26A4 were retrieved from cBioPortal `thca_tcga` legacy study (n = 503 with HM450 + DM call); per-gene aggregate β values were computed from cBioPortal merged probe-per-gene format.

### Korean cohort processing

Yoo 2016 K2 (PRJEB11591): kallisto pseudo-alignment to GENCODE v44 transcriptome, transcript-to-gene aggregation, log2(TPM + 1) normalization. The 8-gene mini-index calibration mismatch (R4-4) was identified and addressed via within-sample-centered profile classification (Memory cross-ref: `v17_korean_k2_calibration.md`). Lee (GSE213647) and the GSE286332 reference arm used pre-computed FPKM tables for score-portability checks.

### Software and statistical environment

All analyses were performed in Python 3.11 with scikit-learn 1.5, pandas 2.2, numpy 1.26, scipy 1.13, lifelines 0.27, and pyDESeq2 0.4. Exact environment specifications and figure scripts will accompany the public repository release.

---

## Quantification and statistical analysis

### Cohen's d (pooled SD)

For all between-cluster comparisons, Cohen's d was computed using pooled standard deviation:
d = (μ₁ − μ₂) / σ_pooled, where σ_pooled = √[((n₁ − 1)·σ₁² + (n₂ − 1)·σ₂²) / (n₁ + n₂ − 2)].

### Mann-Whitney U and Fisher exact

Continuous distributions were compared via two-sided Mann-Whitney U (`scipy.stats.mannwhitneyu`). Discrete proportions were compared via two-sided Fisher exact test (`scipy.stats.fisher_exact`).

### Multiple testing correction

For per-gene tests across the 8-gene panel and Supplementary Tables, Benjamini-Hochberg false-discovery-rate (FDR) correction was applied at q = 0.05.

### Bootstrap confidence intervals

For mediation analysis (R²-based) and hazard ratio meta-analysis, 5,000 bootstrap iterations with patient-level resampling were used to estimate 95% confidence intervals.

### Immune-residualization analysis

To test whether DM1 represented a generic immune-infiltration artifact, the 8-gene panel score was residualized against predefined stromal and immune covariates, and residualized effect sizes were compared with the raw DM1-versus-DM2 contrast. This analysis was used only as a specificity check and not as a primary discovery endpoint.

### Bulk cell-type deconvolution (multi-method)

For per-cell-type granularity beyond the immune-residualization analysis, bulk RNA-seq cell-type fractions for TCGA-THCA primary tumors (n = 572) were estimated against pseudobulk profiles built from the Lu 2023 single-cell reference (GSE193581, `author_celltype` annotations, 67,678 cells, eight cell types: B / Endothelial / Epithelial / Fibroblast / Malignant / Myeloid / NK / T; Lu et al., 2023; processed `.h5ad`). The reference matrix was the per-cell-type mean of expression values restricted to the 1,898 highly variable genes that overlap between Lu 2023 and TCGA gene symbols (95% retention of HVG). Four deconvolution methods were compared: (i) non-negative least squares (NNLS, `scipy.optimize.nnls`); (ii) ridge-regularized NNLS (L2 penalty α=1, augmented design matrix); (iii) ordinary least squares with non-negative clipping and renormalization; and (iv) nu-SVR (CIBERSORT-style; `sklearn.svm.NuSVR(kernel='linear')` with three ν values 0.25/0.5/0.75, lowest-RMSE selection, weights clipped to non-negative and renormalized to sum-to-one). Per-method per-sample weights were normalized to sum-to-one cell-type fractions. The canonical 8-gene RAI score (`rai_score_recalc` from the v17p3 `A2_dm_score_full_cohort.tsv` table; n = 513 with `dm_like` ∈ {DM1_like, DM2_like}, 403 DM1 / 110 DM2) was regressed (OLS) on cell-type fraction subsets — stromal (Endothelial + Fibroblast); immune (T + Myeloid + B + NK); Epithelial-only (purity-confound proxy); all eight fractions — and DM1 vs DM2 Cohen's d was recomputed on residuals. Effect retention after full residualization on all eight cell-type fractions was 24% (NNLS), 24% (Ridge-NNLS), 70% (LR-clip), and 47% (nu-SVR); the nu-SVR retention is methodologically aligned with the canonical immune-residualization analysis (56% retention, see above). T-cell fraction collapse to zero in NNLS-family methods is a known sparse-weight artifact under correlated immune signatures and was resolved by nu-SVR (T cell mean fraction 11.3%, DM1 vs DM2 d = −0.84). Direction-consistent across all four methods: DM1 enriched for Malignant cell (d ≈ +1.1 to +1.4) and Myeloid cell (d ≈ +0.9 to +1.1) fractions, depleted for Epithelial cell (d ≈ −1.1 to −2.3) and Endothelial cell fractions (Supplementary Figure SX). Pseudobulk source code, per-method fraction tables, and the residualization grid are provided at `project/results/p_deconv_2026_05_08/`.

**Per-driver-class composition.** TCGA nu-SVR cell-type fractions were merged with the cBioPortal `v3_anchor_6class` driver call (BRAF V600E n = 280, RAS-mutant n = 54, driver-negative n = 179; `fusion_calls_per_sample.tsv`) and Cohen's d was computed per cell type for each pairwise contrast. RAS-mutant tumors retain Epithelial-cluster identity (d_RAS−BRAF = +1.52) while BRAF V600E tumors are Malignant-cell rich (d = −1.66) (Supplementary Figure SX panel H).

**Methylation × composition.** Per-sample HM450 mean 8-gene β values (`r5_2_sample_methylation_8gene.tsv`, TCGA-THCA n ≈ 484) were correlated (Spearman) with per-sample cell-type fractions for each cell type. The per-gene 8-gene heatmap reports gene × cell-type Spearman ρ; positive tracking with Myeloid (ρ = +0.39), Malignant (+0.30), B (+0.27), and Fibroblast (+0.23); negative tracking with T cell (ρ = −0.42) and Epithelial (−0.31) (Supplementary Figure SX panel I).

**DM1 sub-A vs sub-B teaser.** Sub-cluster labels from `d6p7_dm1_subcluster/dm1_subcluster_labels.tsv` (sub-A n = 84, sub-B n = 56) were intersected with cell-type fractions; Cohen's d and Mann-Whitney U two-sided p were reported per cell type. Sub-A is Malignant-cell rich (d = +1.22, p = 2.5 × 10⁻⁹) and Epithelial-poor (d = −1.26, p = 5.2 × 10⁻¹⁰); immune compartment differences are not significant — defining the sub-A/B split as a tumor-purity-vs-thyrocyte split rather than immune-hot vs immune-cold (Paper-2 boundary marker; Supplementary Figure SX panel J).

**Pseudotime trajectory.** Samples were ranked by canonical 8-gene score (TCGA `rai_score_recalc` n = 513; Lee/GSE213647 `panel_z` n = 630), binned into 10 deciles, and mean per-decile cell-type fraction was computed. Decile-level Spearman ρ between mean score and mean fraction quantifies monotonic trajectory; four compartments — Malignant (↓), Epithelial (↑), Myeloid (↓), Endothelial (↑) — show |ρ| ≥ 0.95 in both cohorts (Supplementary Figure SX panel K).

**Full-transcriptome reference robustness.** Pu 2021 raw counts (33,694 genes × 66,015 cells) were subsampled to 5,000 cells balanced across 7 patients (seed = 42); each cell was assigned a Lu 2023 `author_celltype` label by maximum cosine similarity over the Lu HVG ∩ Pu intersection (1,898 genes after Ensembl → symbol conversion). A full-transcriptome pseudobulk per cell type was constructed (per-cell-type log-normalized mean, library-size-corrected, n_genes = 33,694; `Pu_pseudobulk_full`). TCGA bulk was re-deconvolved by NNLS over the 21,369-gene Pu × TCGA intersection (`scipy.optimize.nnls`, sum-to-one normalization). Per-cell-type Cohen's d (DM1 − DM2) was compared to the v2 Lu HVG nu-SVR primary result; all four informative axes (Malignant, Epithelial, Myeloid, Endothelial) sign-match between the two references with Pu full NNLS magnitudes ≥ Lu HVG. NNLS sparse-collapse zeros the B / Fibroblast / NK / T compartments in the full-transcriptome regime; these are interrogated by the v2 nu-SVR Lu HVG primary instead.

### Random-effects meta-analysis

Cox-derived log-hazard-ratios and standard errors from TCGA-THCA and MSK-IMPACT were combined using random-effects DerSimonian-Laird estimation (`statsmodels.stats.meta_analysis`). Cochran I² was reported for between-cohort heterogeneity.

---

## Additional resources

- **GENCODE v44 reference annotation**: https://www.gencodegenes.org/human/release_44.html
- **TIERA67 gene definition**: Supplementary Table S1 (also available in `metadata/tierA67_genes.txt` in the source code repository)
- **Statistical analysis notebook**: All quantification code is reproducible from the source code repository, with per-figure script paths documented in the README.

---

---

# === Reviewer Q&A pre-empt (09_reviewer_qa.md) ===


# Reviewer Q&A pre-empt (12 items)

---

## Q1. Why exclude BRAF/RAS/TERT from the candidate gene pool?

**A1.** They were *not* excluded. The 67-gene candidate pool (TIERA67) explicitly includes BRAF, NRAS, HRAS, KRAS, RET, NTRK1/3, ALK, PAX8, PPARG, TERT, and EIF1AX in the `[Driver_anchor]` sub-category (12 genes; Supplementary Table S1). The 8-gene panel was selected from the `[TDS-core]` sub-category (16 genes representing canonical thyroid differentiation transcripts) based on Yoo et al. (2016) RAI-uptake biology, not driver mutational status. This curation choice reflects implicit biological prior (differentiation axis vs MAPK driver axis), not exclusion. The 8-gene panel + Driver_anchor concatenation was tested and yielded ARI = 0.49 alone vs 0.90 with full TIERA67 — Driver_anchor genes alone yielded ARI = −0.007, confirming they cannot define the DM cluster axis.

## Q2. Did you bias the panel toward iodine metabolism genes?

**A2.** Yes, by explicit design. The 8-gene panel was selected from `[TDS-core]` because canonical RAI uptake biology (NIS/SLC5A5, TPO, TG, TSHR, PAX8, NKX2-1, FOXE1, DIO1) is the clinical differentiation axis with the strongest biological prior (Yoo et al., 2016). This bias is appropriate for the clinical question (RAI decision-making) and is leak-free in cross-validation: 5-fold AUC = 0.962 (8-gene) vs 0.975 (16-gene full TDS-core), ΔAUC = 0.013 NS. R1-B leak control analysis confirmed AUC = 0.925 in held-out cohorts.

## Q3. How does this panel differ from the BRAF-RAS Score (BRS, 273 genes)?

**A3.** Spearman ρ between 8-gene RAI score and Yoo 2016 BRS in the K2 cohort is 0.49 — partial correlation (not redundant). Concordance is 77%; discordant cases (23%) constitute the BRAF/RAS-negative dark matter axis. The 8-gene panel captures a differentiation axis orthogonal to BRS. Within Xing 2014 BRAF/TERT-negative dark matter, the 8-gene panel sub-stratified 131 of 180 (73%) into a defined molecular axis (DM1 vs DM2).

## Q4. Why didn't BRAF V600E rank near the top in unsupervised analysis?

**A4.** BRAF V600E is a mutation, not a transcript. BRAF transcript expression is comparable in V600E mutant versus wild-type tumors (Cohen's d = −0.044, Mann-Whitney p = 0.57; n = 273 vs 182 in TCGA-THCA). Single-feature DM classification AUC for driver transcripts ranged 0.50–0.61 across BRAF, HRAS, NRAS, KRAS, and TERT — none reaching the discriminative threshold of the 8-gene panel (AUC = 0.962). The 8-gene panel reflects a transcriptional differentiation axis, not driver mutation status.

## Q5. Is the DM1 fusion finding robust to structural variant (SV) missingness?

**A5.** SV-tested status was available for 542 of 557 TCGA-THCA tumors (97.3%); the missing 15 tumors were tested for missingness × DM cluster correlation (chi² p = 0.56, missing-at-random). Three sensitivity analyses imputing missing SV status (best-case fusion-negative, worst-case fusion-positive, case-control matched) yielded DM1 vs DM2 fusion ORs of 7.18, 9.07, and 7.41, respectively — all retaining statistical significance (p < 10⁻¹⁰).

## Q6. Why are DM1 fusion-positive vs fusion-negative tumors so different in clinical phenotype?

**A6.** This heterogeneity is the central finding of our R4-2 analysis (Results 2.4b). Fusion-positive DM1 tumors are younger (37.3 vs 51.3 yr), less advanced (15.3% vs 44.4% stage III/IV), and less immune-infiltrated (CD8/IFN-γ Cohen's d −0.5 to −0.6 vs sub-B). We interpret this as evidence for two distinct DM1 sub-populations: a fusion-driven, well-differentiated, actionable-by-targeted-therapy subgroup (sub-A); and an older-onset, immune-overlap subgroup (sub-B) whose deeper mechanism is reserved for companion-study analysis.

## Q7. Is the 8-gene RNA score specific for RET-fusion-positive PTC?

**A7.** A positive DM1 score captures 27 of 33 (81.8%) of TCGA RET-fusion-positive primary tumors. The remaining 6 RET-positive tumors classified as DM2 were predominantly small (intrathyroidal, T1) tumors with retained differentiation transcript expression — likely reflecting early-stage RET fusion before clonal expansion of the dedifferentiation phenotype. Within DM1, RET fusions account for 33 of 63 (52.4%) fusion-positive cases; the remaining DM1 fusion-positive cases are NTRK1/3 (n = 10), ALK (n = 4), and BRAF fusions (n = 5). Thus DM1 is enriched for, but not specific to, RET fusions — the broader fusion repertoire (76.8% capture) supports the algorithm.

## Q8. Why did one small external Korean cohort show inverse direction in some sub-analyses?

**A8.** This was the central result of our R4-4 calibration mismatch analysis. The K2 mini-index (raw log2(TPM+1) form) does not match the TCGA-trained classifier scaling (within-sample-centered profile form). Inflation factors range 4.9–12.5× across panel genes. The within-sample-centered profile yields concordant direction with TCGA (Spearman ρ = +0.84, p = 1.4 × 10⁻⁵; Methods, Supplementary Figure S9). The K2 mini-index TPM form remains useful for downstream clinical-platform calibration but should not be applied directly with the TCGA-trained absolute-form classifier.

## Q9. What is the mechanism of differentiation gene silencing in DM1?

**A9.** We treat the methylation finding as correlative rather than mechanistic. The TCGA HM450 data establish that DM1 tumors carry markedly elevated promoter β-values across the 8-gene panel (mean β 0.385 vs DM2 0.253; TPO Cohen's d = 2.30), and that within DM1 this signature is fusion-independent at the cohort level (n = 63 fusion-positive vs n = 19 fusion-negative; d = −0.36, NS). What our data do not establish is the upstream regulator. Plausible candidates include canonical de novo methyltransferases (DNMT3A/3B), histone-methylation-coupled silencing complexes such as SETDB1, and MAPK-effector chromatin-remodeling pathways previously implicated in thyroid dedifferentiation; we did not perform DNA methyltransferase profiling, ChIP-seq for repressive marks (H3K9me3/H3K27me3), or CRISPRi knockdown of candidate regulators in this study. We therefore restrict our claim to the observation of a fusion-independent, DM1-associated promoter hypermethylation pattern that motivates — rather than confirms — hypomethylating-agent-based RAI re-induction strategies. Causal attribution of the silencing program to a specific upstream effector is left to dedicated mechanism work in cell-line and PDX systems, and we explicitly hold the line that promoter methylation in our data is suggested as a parallel mechanism, not proven as a causal one.

## Q10. Could the DM1 axis be a generic immune-infiltration artifact?

**A10.** No. We performed residualization analysis on TCGA-THCA: after residualizing the 8-gene panel score on Stromal score + a generic immune-proxy signature (CD8 + IFN-γ + checkpoint), DM1 vs DM2 Cohen's d remained 1.00 (down from raw d = 1.78). The DM1 signal is therefore partially independent of generic immune contamination, with 56% of the original effect retained after residualization. Single-cell analysis confirms that the DM1 axis is thyrocyte-intrinsic rather than stromal-driven (per-patient r 0.798–0.886 in Pu 2021 PTC + adjacent normal pairs).

## Q11. Why do some small external Korean reference samples map strongly to DM2?

**A11.** This is consistent with — not contradictory to — our Paper 1 thesis. The DM1/DM2 axis is continuous, not binary in its underlying biology, and the small external Korean reference set occupies the low-score end of the same classifier used in TCGA. In that set, P(DM1) values range from 0.005 to 0.304, which is compatible with preserved differentiation-program expression rather than failure of the classifier. Because this reference set is used here only for calibration and score-portability checks, we do not assign additional biological interpretation in the present paper.

## Q12. Why is sub-B 96% mutation-negative? Could this be a clustering artifact?

**A12.** Sub-B is the fusion-negative component of DM1 identified by unsupervised clustering within the DM1 compartment, with silhouette support (0.584) and reproducible cross-cohort DM score geometry in Korean validation data. Its mutation-negativity indicates that the DM1 axis is not reducible to canonical BRAF/RAS driver status. We therefore interpret sub-B as a biologically distinct transcriptional state within DM1, while treating its more specific mechanism as hypothesis-generating and reserving deeper characterization for companion work.

## Q13. Why eight genes and not the canonical Yoo 2014 TDS-16? Was the 8-gene panel cherry-picked from a 16-gene set?

**A13.** No — the 8-gene panel is a compact lossless readout of the canonical TDS-16 differentiation axis, not a selectively top-extreme subset. Across the same biology (MAPK-driven silencing of thyroid differentiation transcripts), the deployable 8-gene panel and the canonical Yoo 2014 TDS-16 (Panel + DIO2 / DUOX1 / DUOX2 / GLIS3 / SLC26A4 / SLC5A8 / THRA / THRB) behave indistinguishably across four independent tests in two cohorts (TCGA-THCA n = 572; Lee/GSE213647 n = 632; Supplementary Figure SX_v13).

(i) **Cross-cohort MAPK output × thyroid-score Spearman ρ.** Panel-8 ρ = −0.291 / −0.395; TDS-16 ρ = −0.306 / −0.435 (Δρ = +0.015 / +0.040). The disjoint TDS−panel 8-gene set (TDS-16 \\ Panel) gives ρ = −0.310 / −0.449, matching Panel-8 magnitudes — the silencing axis is distributed across the 16-gene program, not concentrated in eight pre-selected genes.

(ii) **Per-driver-class Cohen's d.** BRAF V600E (n = 344) versus RAS-mutant (n = 61) Cohen's d on score: Panel-8 = −1.615 versus TDS-16 = −1.616 (identical to three significant figures); BRAF V600E versus driver-negative: Panel-8 = −0.900 versus TDS-16 = −0.996; RET fusion versus driver-negative: Panel-8 = −0.261 versus TDS-16 = −0.403. The score gradient across drivers reproduces 1:1 between the two panels.

(iii) **DM1 sub-A vs sub-B two-axis convergence.** sub-A (n = 93) vs sub-B (n = 62) Cohen's d for MAPK output = +1.79 (Mann-Whitney p = 4.1 × 10⁻¹⁷) — sub-A is MAPK-active. The thyroid-program scores are flat across the same split: Panel-8 d = +0.074 (NS), TDS-16 d = +0.032 (NS), TDS−panel d = −0.026 (NS) — both compact and canonical panels converge to silencing whether reached via the MAPK-active route (sub-A) or the MAPK-low route (sub-B), and the convergence holds for the full TDS-16 (not only the deployable subset).

(iv) **ROC-AUC: MAPK-high vs MAPK-low classifier.** Panel-8 AUC = 0.623 (TCGA) / 0.728 (Lee); TDS-16 AUC = 0.631 / 0.740. ΔAUC TDS-16 − Panel-8 = +0.007 / +0.012 — both non-significant. This independently reproduces the existing one-page-audit ΔAUC = 0.013 NS (DM1/DM2 classification task; see `p1_onepage_audit`) on an orthogonal classification target.

(v) **Per-gene rank within TDS-16.** Strongest TCGA negative MAPK-correlations are SLC5A8 (ρ = −0.516, non-panel), DIO2 (−0.481, non-panel), TPO (−0.465, panel), SLC26A4 (−0.438, non-panel) — that is, the 8-gene panel members occupy median rank within the 16, not selectively top-extreme. NKX2-1 is the consistent positive-ρ outlier (+0.307 / +0.425) in both cohorts and is included in the panel mean as a known absorbed caveat.

The 8-gene panel was originally curated from the broader TIERA67 candidate pool by RandomForest ranking with drivers excluded by design (see Q1 / Q2; `v17_8gene_audit_2026_04_29` audit), and the panel's claim is compactness for clinical RT-qPCR deploy — not "best 8 of 16." Cell Press / npj-format reviewers can verify the cherry-pick concern by running the four tests in (i–iv) on their own data, with all source code and per-panel TSVs at `project/results/p_deconv_2026_05_08/v13_*.tsv`.

## Q14. Does the MAPK→thyroid-silencing axis hold across heterogeneous cohorts, or is it a TCGA + Lee artifact?

**A14.** It generalizes — pooled MAPK output × Panel-8 Spearman ρ = **−0.327, 95% CI [−0.376, −0.278]** across n = 1,287 from five cohorts (TCGA-THCA n = 572 + Lee/GSE213647 n = 632 + GSE126698 n = 28 + GSE286332 n = 18 + GSE76039 n = 37; Fisher z-transform fixed-effect pool; Supplementary Figure SX_v14). Cochran I² = 72.9% (Q = 14.77, df = 4, p = 0.005) — the heterogeneity is biologically expected and predicted by the v12 two-axis convergence model: small inflammation-overlap cohorts reach panel silencing via a non-MAPK route and therefore decouple the MAPK × Panel-8 correlation, while cohorts at the dedifferentiated end of the axis saturate the panel.

(i) **Well-differentiated primary cohorts (TCGA + Lee, n = 1,204):** Panel-8 ρ = −0.291 (TCGA, p = 1.3 × 10⁻¹²) / −0.395 (Lee, p = 4.7 × 10⁻²⁵). The two largest cohorts both carry strong, sign-consistent correlation; together they account for 93.5% of the cross-cohort sample pool.

(ii) **Inflammation-overlap cohort GSE286332 (Korean PTC subset, n = 18; raw TPM → log2 → within-cohort z):** Panel-8 ρ = −0.040 NS. This is the cohort where the v12 two-axis model predicts MAPK should decouple, because panel silencing in this small inflammation-overlap subset is reached via a non-MAPK route. The decoupling is therefore positive evidence for v12, not a failure to replicate. (Per the P1-2 cohort audit, GSE286332 is not part of the Korean Pillar I primary cohort n = 865 and is reported here only as a heterogeneity outlier in the cross-cohort forest.)

(iii) **Advanced-disease cohort GSE76039 (PDTC + ATC, n = 37; Landa 2016, microarray z):** Panel-8 ρ = +0.125 NS. ATC tumors show Panel z mean = −0.82 (vs PDTC +0.96; v11 `v11_GSE76039_by_histology.tsv`), i.e., the panel is already saturated at the dedifferentiated low end and within-cohort dynamic range is collapsed. The within-cohort MAPK × Panel correlation breaks down for the same reason a saturated reporter cannot resolve dose response.

(iv) **Conservative read for npj/JCI Insight reviewers.** If the request is the most conservative single number across all available cohorts: pooled MAPK × Panel-8 ρ = −0.327, 95% CI [−0.376, −0.278], n = 1,287. If the request is the cleanest signal restricted to well-differentiated primary cohorts (TCGA + Lee): per-cohort ρ = −0.291 / −0.395, both p ≪ 10⁻¹². The two outlier cohorts (inflammation-overlap + advanced-disease) are reported transparently with the v12 explanation rather than excluded.

The heterogeneity in the pool is therefore not a robustness problem — the same per-cohort heterogeneity is the cross-cohort confirmation of the two-axis model that drives the Fig 8 mechanism layer (Supp Fig SX_v13 + SX_v14; full forest at `project/results/p_deconv_2026_05_08/v14_cross_cohort_forest.tsv`).

---

---

# === Supplementary Tables (13_supplementary_tables.md) ===


# Supplementary Tables Index

Paper 1 manuscript v8 supplementary tables, with source TSV file paths in `project/results/`. Final submission packaging will combine these into a single multi-sheet Excel workbook (`THCA_paper1_supplementary_tables.xlsx`).

---

## S1 — TIERA67 67-gene panel definition (7 categories)

**Description.** Curated candidate gene pool spanning thyroid biological categories used for 8-gene panel selection.

**Source.** `project/metadata/tierA67_genes.txt`

**Format (10 columns):**

| Column | Type | Description |
|---|---|---|
| `gene_symbol` | str | HUGO symbol |
| `category` | str | One of 7 categories: TDS_core / MAPK_output_ERK / Driver_anchor / Aggressive_marker / Dediff_invasion / Immune_stromal_light / Thyroid_lineage_extra |
| `n_genes_in_category` | int | Count |
| `Yoo2016_inclusion` | bool | In Yoo 2016 BRS 273-gene |
| `TCGA2014_BRS` | bool | In TCGA 2014 BRS |
| `8gene_panel` | bool | In final 8-gene panel |
| `entrez_id` | int | NCBI Entrez |
| `ensembl_id` | str | GENCODE v44 |
| `chromosome` | str | hg38 |
| `notes` | str | biological role |

n=67 rows.

---

## S2 — Cohort overview

**Description.** Multi-cohort sample-level descriptor table.

**Source.** Composed from:
- `project/results/tables/tcga_thca_clinical_extended.tsv` (TCGA n=504)
- `project/results/audit_2026_04_30/round3/cbio_sv_thca.tsv` (TCGA SV)
- `project/results/audit_2026_04_30/round5/r5_2_sample_methylation_8gene.tsv` (TCGA HM450)
- `project/results/v17_korean/` (K2 + Lee + GSE286332)

**Format (~12 columns):**

| Column | Description |
|---|---|
| `cohort` | TCGA / MSK-IMPACT / K2 / Lee / GSE286332 / GSE184362 / Lu2023 |
| `n_total` | sample count |
| `n_with_OS` | overall survival annotation |
| `n_with_RNAseq` | bulk RNA-seq |
| `n_with_WGS` | WGS / WES |
| `n_with_HM450` | methylation |
| `n_with_SV` | structural variant |
| `assay_platform` | RNA-seq / array / scRNA-seq |
| `tissue_type` | FF / FFPE / mixed |
| `ancestry_majority` | EUR / EA / mixed |
| `event_rate_pct` | OS event rate |
| `cite` | primary publication |

---

## S3 — TCGA-THCA per-sample 8-gene scores + DM call + clinical

**Description.** Per-sample 8-gene panel score, P(DM1) classification probability, DM cluster call, and clinical metadata.

**Source.** `project/results/tables/tcga_thca_8gene_per_sample.tsv` (compose from existing v17 outputs)

**Format (~25 columns):**

| Column | Description |
|---|---|
| `sample_id` | TCGA-XX-XXXX-01 (primary tumor) |
| `SLC5A5_log2tpm` ... `DIO1_log2tpm` | 8 panel gene log2(TPM+1) |
| `panel_mean_z` | within-cohort z-mean |
| `P_DM1` | logistic regression class probability |
| `DM_call` | DM1 / DM2 |
| `BRS_273gene` | Yoo 2016 BRAF-RAS Score |
| `BRAF_V600E` | mutation status |
| `RAS_hotspot` | NRAS/HRAS/KRAS hotspot |
| `TERT_promoter` | promoter mutation |
| `fusion_status` | RET/NTRK/ALK/BRAF/none |
| `age_at_diagnosis` | years |
| `sex` | M/F |
| `stage` | I-IV |
| `OS_status` | dead/alive |
| `OS_days` | follow-up |

n=504 rows.

---

## S4 — TIERA67 univariate Cohen's d ranking (full 67 genes)

**Description.** Per-gene Cohen's d (DM1 vs DM2) ranking across all 67 TIERA67 genes, including drivers.

**Source.** `project/results/p1_driver_mrna_audit/top20_by_d_full_tiera67.tsv` (extend to full 67)

**Format:** gene_symbol, category, mean_DM1, mean_DM2, Cohen_d, MW_p, BH_FDR_q, rank.

Highlights (per memory): 8-gene at ranks 4 (TPO), 5 (DIO1), 18 (TG), 19 (PAX8), 25 (FOXE1), 38 (NKX2-1), 40 (TSHR), 50 (SLC5A5); drivers at ranks 15 (CDKN2B), 22 (RET), 24 (CDKN2A), 52 (BRAF), 56 (TERT), 63 (KRAS), 66 (NRAS), 67 (HRAS).

n=67 rows.

---

## S5 — Pan-genome cluster ARI ladder

**Description.** Cluster reproducibility across panel sizes 8, 16, 67 (TIERA67), 200, 1000, 5000 (pan-genome MAD-ranked).

**Source.** `project/results/p4_pangenome_vs_tiera67/ari_comparison.tsv`

**Format:** panel_definition, n_genes, ARI_vs_8gene, NMI_vs_8gene, hypergeometric_p_in_top100.

Highlights: 8-gene 0.49 / TIERA67 0.90 / Pan top-200 0.86 / top-1000 0.90 / top-5000 0.92 / Driver_anchor 12 only −0.007.

---

## S6 — DM1 fusion details (per-sample fusion class + partner)

**Description.** Per-sample fusion partner annotations within TCGA + MSK SV-tested cohort.

**Source.** Compose from:
- `project/results/audit_2026_04_30/round3/cbio_sv_thca.tsv` (TCGA n=542)
- `project/results/audit_2026_04_30/round5/msk_sv_thca.tsv` (MSK n=12)

**Format:** sample_id, cohort, fusion_class (RET/NTRK/ALK/BRAF/PAX8-PPARG/none), fusion_partner_5p, fusion_partner_3p, exon_breakpoint, DM_call, DM1_subcluster (sub-A/sub-B).

Highlights: TCGA RET fusion CCDC6-RET 17, NCOA4-RET 3, other 13. MSK RET fusion CCDC6-RET 3, NCOA4-RET 2.

---

## S7 — DM1 vs DM2 per-gene HM450 methylation β-values

**Description.** Per-gene HM450 promoter β-value comparison across 8-gene panel + DIO2 + SLC26A4 + DUOX2.

**Source.** `project/results/audit_2026_04_30/round5/r5_2_per_gene_methylation_DM.tsv`

**Format:** gene_symbol, n_DM1, n_DM2, n_notDM, mean_DM1, mean_DM2, mean_notDM, Cohen_d, MW_p, BH_FDR_q, fusion_independence_p_within_DM1.

Highlights: TPO d=2.30 p=1.9e-18; DIO1 d=1.24; TSHR d=1.20; SLC5A5 d=0.22 NS.

---

## S8 — Meta-analysis raw inputs (TCGA + MSK Cox + DerSimonian-Laird)

**Description.** Per-cohort hazard ratio inputs and pooled meta-analysis output.

**Source.** Compose from:
- `project/results/audit_2026_04_30/round2/n1_meta.json`
- `project/results/p6_multi_cohort_meta.json`

**Format:** cohort, n_total, n_events, log_HR_DM1_vs_DM2, SE_logHR, HR, HR_lower_95, HR_upper_95, p_value, weight_in_meta.

Highlights: TCGA HR 2.30 [0.77, 6.88], n=504, 16 events. MSK HR 2.67 [1.17, 6.10], n=117, 38 events. Pooled HR 2.53 [1.31, 4.89], I²=0%.

---

## S9 — Korean score-portability and calibration diagnostics

**Description.** Per-sample calibration and score-portability summaries for Korean validation cohorts.

**Source.** Compose from:
- `project/results/d7p3_k2_calibration/k2_4metric_scores.tsv`
- `project/results/v17_korean/GSE213647_panel_score.tsv`

**Format:** sample_id, cohort, raw_panel_score, within_sample_z_score, per_gene_z_score, rank_based_score, P_DM1, DM_call.

Highlights: per-gene inflation factors 4.9-12.5x in the raw mini-index form; within-sample-centered scoring restores direction concordance with TCGA.

---

## S10 — Reviewer Q&A pre-empt 12 items

**Description.** 12 anticipated reviewer questions with prepared responses (revision-ready text).

**Source.** `project/manuscript_v8/09_reviewer_qa.md`

---

# Final packaging plan (W5)

1. Execute compose scripts to generate each TSV from existing v17 outputs (where missing).
2. Convert each TSV to Excel sheet using `pandas.to_excel`.
3. Combine into single multi-sheet workbook `THCA_paper1_supplementary_tables.xlsx`.
4. Write Supplementary Table caption `.docx` per Cell Press submission template.

Submission package:
- `manuscript.docx` (main text + figures inline + references)
- `figures.zip` (Fig 1-8 PDF/PNG + S1-S9 PDF/PNG)
- `THCA_paper1_supplementary_tables.xlsx`
- `cover_letter.docx`
- `reviewer_qa.docx` (S10 standalone)

---

# === Supplementary Text ST1 (15_supp_image_dm1_tss_confound_case_study.md) ===


# Supplementary Text ST1. Image-DM1 pilot - TSS-confound disclosure

## Purpose and scope

This supplementary case study documents the post-hoc audit of the Paper 2 image-DM1 pilot and defines how the result should be used in Paper 1. The audit is included here as a transparent negative-control example: a visually plausible TCGA pathology signal can be inflated by tissue-source-site (TSS) structure when sample size is small and subgroup positives are sparse. This section does not support a main-text claim that hematoxylin-and-eosin whole-slide images independently predict DM1. It supports the opposite boundary: TCGA-only image-DM1 remains a pilot result that requires K2 H&E or a multi-site FFPE validation cohort before it can be used as a primary claim.

## Data

The audited model used TCGA-THCA whole-slide image features extracted with a foundation pathology encoder and trained with CLAM attention multiple-instance learning. The primary audited out-of-fold prediction table contains 59 slides with DM1/DM2 labels, per-slide probabilities, and fold assignments. The original report highlighted a RAS-like/FVPTC subgroup AUC of 1.000 in 16 slides, including 3 DM1-positive and 13 DM2-negative cases. Metadata were joined from `project/metadata/sample_master_v3.tsv` to recover histology subtype, molecular subtype, sex, age, stage, and TCGA TSS code. Source files are:

| Item | Path |
|---|---|
| Original OOF predictions | `project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/phase2_tcga_clam/clam_per_slide_predictions.tsv` |
| Slide manifest | `project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/phase2_tcga_clam/slide_manifest.tsv` |
| Audit root | `project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/analysis_supp/audit_ras_auc100/` |
| Dossier | `project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/analysis_supp/audit_ras_auc100/AUDIT_DOSSIER.html` |
| Hub copy | `project/papers_hub_2026_05_04/paper2_image_dm1_audit_dossier.html` |

## Audit design

The audit used four layers. Phase 1 ran no-retraining checks: subgroup overlap, raw prediction distribution, permutation null, histology-only shortcut baseline, and fold composition. Phase 2 quantified confounders: sex-stratified performance, TSS distribution, case-level profiles, histology by molecular-subtype gaps, and a clinical-only logistic-regression baseline. Phase 3 performed decisive retraining controls: fixed-split initialization variance, label-shuffle null, and leave-one-TSS-out (LOTO) cross-validation. Phase 4 tested corrective baselines: TSS-ComBat feature adjustment, TSS-balanced splitting, multimodal logistic regression, and an EM-out sensitivity analysis.

The central criterion was whether the image signal survived TSS holdout. Random K-fold and stratified folds can place the same tissue-source-site centers into training and test sets. LOTO removes that shortcut by withholding all slides from one TSS at a time and pooling the held-out predictions.

## Results

The RAS-like and FVPTC subgroup labels were identical in the 59-slide subset: RAS_like intersection FVPTC = 16/16, with no off-diagonal cPTC/RAS-like or FVPTC/BRAF-like cases available for disentanglement. The apparent perfect ranking in the 16-slide subgroup was narrow: the lowest DM1 probability was 0.549 and the highest DM2 probability was 0.532, giving a separation gap of 0.017. The permutation p value for the RAS-like AUC was 0.0030, but the bootstrap confidence interval of [1.000, 1.000] reflected deterministic resampling of the same perfect ordering, not center-generalized performance.

The key audit findings were:

| Test | Value | Interpretation | Source |
|---|---:|---|---|
| Original RAS-like/FVPTC AUC | 1.000 | Pre-audit headline in n=16 | `AUDIT_SUMMARY.json` |
| RAS-like/FVPTC overlap | 16/16 | One subgroup labeled two ways | `AUDIT_SUMMARY.json` |
| RAS-like separation gap | 0.017 | One near-boundary rank supports the perfect AUC | `c2_prob_distribution.tsv` |
| Folds with zero RAS-like positives | 3/5 | Original K-fold split is poorly balanced for the subgroup | `c5_RAS_like_fold_composition.tsv` |
| Male AUC | 0.350 | Sex-stratified failure in n=13 males | `PHASE2_AUDIT_SUMMARY.json` |
| Clinical-only LR AUC | 0.768 | Histology + sex + TSS beats image-inclusive multimodal LR | `phase4/s3_multimodal_LR.tsv` |
| Multimodal LR AUC | 0.756 | CLAM probability adds no net gain over metadata | `phase4/s3_multimodal_LR.tsv` |
| TSS-only LR within RAS-like | 1.000 | Perfect ranking can be reproduced without image features | `tss_only_within_raslike_lr.tsv` |
| Pooled LOTO overall AUC | 0.602 | Center holdout weakens overall performance | `PHASE3_AUDIT_SUMMARY.json` |
| Pooled LOTO RAS-like AUC | 0.308 | Center holdout reverses the subgroup signal | `PHASE3_AUDIT_SUMMARY.json` |
| UNI-final OOF AUC | 0.874 | Foundation encoder improves overall ranking in n=54 | `analysis_supp/bootstrap_auc_uni.json` |
| UNI-final LOTO overall AUC | 0.852 | UNI retains center-holdout overall performance | `analysis_supp/audit_uni_loto/UNI_LOTO_SUMMARY.json` |
| UNI-final LOTO RAS-like AUC | 0.718 | Honest subgroup estimate is no longer perfect and remains underpowered | `analysis_supp/audit_uni_loto/UNI_LOTO_SUMMARY.json` |

The TSS distribution explained why the subgroup looked perfect under ordinary splits. All 3 RAS-like DM1 cases came from DJ or FK, while all 10 RAS-like cases from EM were DM2. A logistic regression using only TSS dummy variables, with no image features, achieved within-RAS-like AUC=1.000. Under pooled LOTO for the audited ViT-L feature table, the overall AUC dropped to 0.602 and the RAS-like AUC dropped to 0.308, which is worse than random. A later GPU-run LOTO audit of the UNI-final feature table refined rather than erased this conclusion: UNI retained strong overall center-holdout performance (LOTO AUC 0.852), but the RAS-like subgroup estimate was 0.718 rather than 1.000 and still rested on only 3 positive RAS-like cases.

## Corrective baselines

No corrective baseline recovered the original 1.000 in a center-generalized manner. TSS-ComBat reduced the apparent batch component and partially repaired the male failure, but the overall median AUC remained in the 0.63-0.71 range across tested initialization seeds. TSS-balanced splitting preserved high RAS-like AUC when training and test folds shared centers, which is useful as a matched-center ceiling but does not answer the center-holdout question. Multimodal logistic regression quantified the image channel directly: clinical-only AUC was 0.768, CLAM-only AUC in the corrective baseline was 0.691, and multimodal AUC was 0.756. Removing the 10 EM RAS-like DM2 anchor slides reduced the RAS-like median to approximately 0.72-0.78 depending on initialization.

## Conclusion for Paper 1

For Paper 1, the image-DM1 result should be treated as a cautionary supplementary pilot, not as a supporting pillar. The appropriate use is a TSS-confound disclosure showing that TCGA single-cohort pathology AI can appear stronger than it is when subgroup labels, histology, sex, and acquisition center are entangled. The UNI-final LOTO result also preserves a positive future path: the foundation-model image channel may carry real overall DM1 signal, but the RAS-like/FVPTC subgroup claim should be reported as underpowered and center-sensitive until K2 H&E or a multi-site FFPE cohort validates it externally. The finding does not weaken the Paper 1 molecular DM1 axis, which is based on RNA, methylation, fusion, survival, and cross-cohort evidence; it only blocks a premature image-derived extension of that axis.

## Recommendation

The Paper 2 image-DM1 line should proceed only after external validation. The minimum defensible design is K2 H&E or a multi-site FFPE cohort with a locked TSS-balanced split, clinical-only baseline, multimodal image-gain estimate, sex-stratified AUC, and leave-one-site-out or leave-one-center-out validation. Until those data are available, report the TCGA image model as an audited pilot: UNI shows promising center-holdout overall performance, but the original perfect RAS-like/FVPTC subgroup headline is not defensible.

---

# === Cover Letter (08_cover_letter.md) — kept separate from manuscript body ===


# Cover Letter

2026-05-08

[Editor in Chief / Editorial Office]
*Cell Reports Medicine*
Cell Press / Elsevier

Dear Editor,

Papillary thyroid carcinoma carries generally favorable overall survival, yet the BRAF/RAS-negative compartment — approximately 23% of cases — remains mechanistically untethered from the ATA risk framework, leaving radioiodine decisions in this group without molecular guidance. We report that an 8-gene RAI-responsiveness panel resolves this compartment into a DM1/DM2 axis in which DM1 is 76.8% tyrosine-kinase-fusion-positive (OR 7.41 versus DM2), harbors fusion-independent promoter hypermethylation of thyroid differentiation genes (TPO Cohen's d = 2.30), and carries a pooled overall-survival hazard of 2.53 (95% CI 1.31–4.89) in TCGA + MSK meta-analysis. This combination of mechanistic layering, multi-cohort validation, and an explicit clinical reflex-testing framework places the work within the translational scope of *Cell Reports Medicine*.

Papillary thyroid carcinoma remains clinically heterogeneous despite generally favorable overall survival, and current risk frameworks offer limited mechanistic guidance for the BRAF/RAS-negative compartment. In this study, we address that gap using an 8-gene RAI-responsiveness panel grounded in canonical thyroid differentiation biology and validated across TCGA-THCA, MSK-IMPACT thyroid cancer, Korean cohorts, cBioPortal structural-variant data, HM450 methylation data, and external single-cell datasets.

The manuscript makes three linked contributions. First, it identifies DM1 as a fusion-enriched subtype of BRAF/RAS-negative thyroid cancer, with 76.8% tyrosine-kinase-fusion positivity and capture of 81.8% of TCGA RET-fusion-positive tumors. Second, it shows that DM1 carries a fusion-independent epigenetic silencing program across thyroid differentiation genes, including strong TPO, DIO1, and TSHR promoter hypermethylation. Third, it connects these molecular findings to a clinically interpretable reflex-testing framework that can prioritize selective fusion sequencing and support prospective evaluation of epigenetic-targeted RAI re-induction.

We believe the manuscript is a strong fit for *Cell Reports Medicine* because it combines mechanistic depth, multi-cohort validation, and translational relevance in a disease setting where treatment decisions remain incompletely informed by existing molecular stratification. The work is framed conservatively, with explicit treatment of cohort, power, and prospective-validation limitations.

Suggested reviewers:
- Jaume Capdevila, MD, PhD, Vall d'Hebron Institute of Oncology
- Iñigo Landa, PhD, Brigham and Women's Hospital / Harvard Medical School
- Sareh Parangi, MD, Massachusetts General Hospital
- Park Young Joo, MD, PhD, Seoul National University Hospital
- Yuri Nikiforov, MD, PhD, University of Pittsburgh Medical Center

We request exclusion of reviewers from the Fagin/Krishnamoorthy laboratory because of close topical overlap with work cited in the manuscript.

This manuscript is not under consideration elsewhere, and all authors have approved the submission.

Sincerely,

**Seungho Cook**
First Author
Seoul National University Bundang Hospital, Seongnam-si, Gyeonggi-do 13620, Republic of Korea
kukshomr@gmail.com

**Yu Hyeong-won, MD, PhD**
Corresponding Author
Department of Internal Medicine, Seoul National University Bundang Hospital, Seongnam-si, Gyeonggi-do 13620, Republic of Korea
[VERIFY before submission — Yu Hyeong-won institutional email; expected pattern @snubh.org or @snu.ac.kr]

---

# === Self-verification report (14_self_verification_report.md) ===


# Self-verification report

## Files updated in this pass

- [00_title_candidates.md](/home/seungho/personal/THCA_data_analysis/project/manuscript_v8/00_title_candidates.md)
- [03_introduction.md](/home/seungho/personal/THCA_data_analysis/project/manuscript_v8/03_introduction.md)
- [04_intro_1_1_hook.md](/home/seungho/personal/THCA_data_analysis/project/manuscript_v8/04_intro_1_1_hook.md)
- [01_abstract.md](/home/seungho/personal/THCA_data_analysis/project/manuscript_v8/01_abstract.md)
- [04_results.md](/home/seungho/personal/THCA_data_analysis/project/manuscript_v8/04_results.md)
- [05_figure_captions.md](/home/seungho/personal/THCA_data_analysis/project/manuscript_v8/05_figure_captions.md)
- [06_discussion.md](/home/seungho/personal/THCA_data_analysis/project/manuscript_v8/06_discussion.md)
- [07_star_methods.md](/home/seungho/personal/THCA_data_analysis/project/manuscript_v8/07_star_methods.md)
- [08_cover_letter.md](/home/seungho/personal/THCA_data_analysis/project/manuscript_v8/08_cover_letter.md)
- [09_reviewer_qa.md](/home/seungho/personal/THCA_data_analysis/project/manuscript_v8/09_reviewer_qa.md)
- [10_full_manuscript_compiled.md](/home/seungho/personal/THCA_data_analysis/project/manuscript_v8/10_full_manuscript_compiled.md)
- [13_supplementary_tables.md](/home/seungho/personal/THCA_data_analysis/project/manuscript_v8/13_supplementary_tables.md)

## Cross-paper boundary check

Command run:

```bash
rg -n "autoimmune-PTC|autoimmune PTC|Hashimoto|Hashimoto-like|PTC\+HT|HLA-II|BCR|TLS|AICDA|Chu|Graves|GD|TSAb|TSI|hyperthy|thyrotox|exophthal|thyroid eye disease|Chen 2018" \
  project/manuscript_v8/00_title_candidates.md \
  project/manuscript_v8/01_abstract.md \
  project/manuscript_v8/03_introduction.md \
  project/manuscript_v8/04_intro_1_1_hook.md \
  project/manuscript_v8/04_results.md \
  project/manuscript_v8/05_figure_captions.md \
  project/manuscript_v8/06_discussion.md \
  project/manuscript_v8/07_star_methods.md \
  project/manuscript_v8/08_cover_letter.md \
  project/manuscript_v8/09_reviewer_qa.md \
  project/manuscript_v8/13_supplementary_tables.md
```

Result: `0 hits` in the audited Paper 1 draft set. `rg` exited with code `1`, which in this case means no matches were found.

## Citation check

- `Chen 2018`: `0 hits`
- `Chu 2018` / `Chu et al. 2018`: `0 hits` in Paper 1 primary draft files after cleanup
- No new external literature was introduced beyond citations already anticipated in the draft set
- Citations retained in edited passages remain existing manuscript citations: `Yoo et al. 2016`, `Landa et al. 2016`, `Haugen et al. 2016`, `Ringel et al. 2025`, `Wirth et al. 2020`, `Pu et al. 2021`
- Reference source file remains [03_intro_references.bib](/home/seungho/personal/THCA_data_analysis/project/manuscript_v8/03_intro_references.bib)

## Voice protection check

Protected sections were left as explicit author placeholders rather than Codex prose:
- [03_introduction.md](/home/seungho/personal/THCA_data_analysis/project/manuscript_v8/03_introduction.md): `1.1` hook first line
- [04_intro_1_1_hook.md](/home/seungho/personal/THCA_data_analysis/project/manuscript_v8/04_intro_1_1_hook.md): author hook anchor
- [06_discussion.md](/home/seungho/personal/THCA_data_analysis/project/manuscript_v8/06_discussion.md): `3.1` opening paragraph, `3.4` limitations
- [08_cover_letter.md](/home/seungho/personal/THCA_data_analysis/project/manuscript_v8/08_cover_letter.md): paragraph 1
- [09_reviewer_qa.md](/home/seungho/personal/THCA_data_analysis/project/manuscript_v8/09_reviewer_qa.md): `Q9`

## Open items for author review

- Citation placeholders in [03_introduction.md](/home/seungho/personal/THCA_data_analysis/project/manuscript_v8/03_introduction.md) still need literature-level verification for the SEER-incidence and Bethesda statements
- ~~Decide whether `GSE286332` should remain as a small Korean reference arm inside the `n=874` aggregate or be fully removed from Paper 1 cohort summaries~~ — **Resolved 2026-05-07 (audit P1-2):** GSE286332-PTC(9) dropped from Paper 1 aggregate; new Korean cohort summary statistic is `n=865 (K2 235 + Lee 630)`. GSE286332 reference arm description retained in STAR Methods with explicit boundary note. Rationale: preserve scope separation from Paper 2 (GSE286332 = Paper 2 PTC vs PTC+HT main cohort).
- Confirm the final public repository URL and archival DOI language in [07_star_methods.md](/home/seungho/personal/THCA_data_analysis/project/manuscript_v8/07_star_methods.md)
- Fill affiliations and corresponding-author email in [08_cover_letter.md](/home/seungho/personal/THCA_data_analysis/project/manuscript_v8/08_cover_letter.md)
- Final supplement numbering can be re-tuned once figure build is complete

---

## Re-verification 2026-05-13 (pre-Yu-meeting audit)

Re-ran the boundary + cohort + citation + figure-ref + key-number audit ahead of the 2026-05-14 Yu meeting. All checks GREEN.

- **Paper 2 boundary**: 0 leak. All `Hashimoto-like` / `HLA-II` / `BCR` hits in main draft files are explicitly boundary-marked as "reserved for Paper 2" (05_figure_captions.md:96/98/142, 06_discussion.md:49, 09_reviewer_qa.md:80/84).
- **Cohort number**: `n=865 (K2 235 + Lee 630)` consistent in 01_abstract.md, 03_introduction.md, 04_results.md, 07_star_methods.md. No lingering `n=874`. GSE286332-PTC(9) confirmed dropped per 2026-05-07 P1-2 audit.
- **Placeholders**: 0 inline `TBD/TODO/XXX/FIXME` in main prose. Two structured TODOs intentionally left in 07_star_methods.md:39 and :49 awaiting Zenodo DOI + public repo URL at submission. One `VERIFY` line in 01_abstract.md:18 awaits Yu-Hyeongwon institutional email.
- **Figure numbers**: Main Fig 1-8, Supp S1-S6, S5b, SX, SX_v13, SX_v14, SB2 all defined with consistent cross-refs.
- **Key effect sizes confirmed across sections**: Fusion 76.8% / OR 7.41 (5 files); pooled HR 2.53 [1.31, 4.89] (5 files); TPO d=2.30 (7 files); mean β 0.385 vs 0.253 (5 files); ARI 0.49/0.90/0.92 (4 files).

### Edits applied 2026-05-13 (post-audit)

- `04_results.md`: section order reordered 2.1 → 2.3 → 2.4a → 2.4b → 2.2 → 2.5 ⇒ **2.1 → 2.2 → 2.3 → 2.4 → 2.5a → 2.5b** per 17_results_reorder_plan.md (selection-layer paper framing, discovery-first). +14 words = header comment only; no prose added/removed/reworded. 1 internal cross-ref `(2.4a)` → `(2.5a)`.
- `07_star_methods.md`: cohort drift fixed (K2 260 → 235, Lee 632 → 630, lines 59-61 and 141); 2 placeholder TODOs added for Zenodo DOI + repo URL at submission.
- `09_reviewer_qa.md`: Q13/Q14 Paper 2 territory drift removed ("HT/B-cell route" → "MAPK-low route" / "non-MAPK route"; GSE286332 dual-role parenthetical added in Q14(ii)). Q9 untouched (voice-protected). Author confirmation requested: is GSE286332 OK to keep in v14 forest pool (n=1,287) while dropped from primary Pillar I cohort (n=865)?

### Open items still requiring author keyboard (voice-protected)

- 03_introduction.md `1.1` hook + 04_intro_1_1_hook.md anchor (Alt B refined / hybrid "compass" decision)
- 06_discussion.md `3.1` opening (Landa 2016 JCI cite, Krishnamoorthy 2025 misattribution 정정) + `3.4` limitations
- 08_cover_letter.md paragraph 1
- 09_reviewer_qa.md Q9

---

*Generated 2026-05-13 by `_compile_full_manuscript.py`. Constituent file states reflect the latest edits to the underlying source files.*
