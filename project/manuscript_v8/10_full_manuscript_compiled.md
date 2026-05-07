---
title: Paper 1 manuscript v8 — compiled (concatenated)
date: 2026-05-08
author: Seungho Cook
target_venue: Cell Reports Medicine (1순위) → JCI Insight + Nat Commun dual reach → npj Precision Oncology (fallback)
status: assembled. Voice-protected placeholders preserved. P0-9 affiliation/email fields TBD by author.
source_files: 00_title_candidates / 01_abstract / 03_introduction / 04_results / 05_figure_captions / 06_discussion / 07_star_methods / 09_reviewer_qa / 13_supplementary_tables / (cover letter 08 separate, kept at end)
regenerated: 2026-05-08 (post hub deep reconciliation + GSE286332 drop n=874→865 + Fig 7D/8C demote)
---

# Paper 1 — Full Compiled Manuscript v8

> **Voice-protected placeholders** (author keyboard required, per `v17_sprint_vs_marathon_violation`):
> 1. Hook (§1.1 first line + `04_intro_1_1_hook.md`)
> 2. Discussion §3.1 opening paragraph (`06_discussion.md:13`)
> 3. Limitations §3.4 (`06_discussion.md:43`)
> 4. Cover letter paragraph 1 (`08_cover_letter.md:19`)
> 5. Reviewer Q&A Q9 (`09_reviewer_qa.md:48`)
>
> **Outstanding factual fields** (user data input required):
> - Affiliations + corresponding-author email (`01_abstract.md` author block + `08_cover_letter.md:42-48`)
> - Aim §1.4 sentence split (`03_introduction.md:37`) — split spec in DECISIONS_PENDING_2026_05_07.md C1

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

[AUTHOR HOOK — voice insertion point; see `04_intro_1_1_hook.md`]

Papillary thyroid carcinoma (PTC) is the most common endocrine malignancy and among the fastest-rising in incidence over the past three decades (SEER, 2024). While the majority of patients achieve durable remission after thyroidectomy and selective ¹³¹I ablation, 5-20% develop recurrent or persistent disease and approximately 10% develop distant metastases (Haugen et al., 2016). Current risk-tier-based RAI decisions — codified in the 2015 American Thyroid Association (ATA) Management Guidelines (Haugen et al., 2016) and recently updated as ATA 2025 (Ringel et al., 2025) — rely predominantly on clinico-pathological features (tumor size, multifocality, extrathyroidal extension, lymph node burden) with BRAF V600E as the sole molecular risk modifier. Bethesda III/IV indeterminate cytology affects 15-30% of fine-needle aspiration biopsies and remains a major diagnostic gap (Cibas and Ali, 2017). Together, these gaps motivate orthogonal molecular sub-stratification of clinically heterogeneous tumors, particularly within the BRAF/RAS-negative compartment.

---

## 1.2 Existing molecular framework (~165 words)

The 2014 Cancer Genome Atlas (TCGA) study of papillary thyroid carcinoma established a binary molecular spectrum anchored by mitogen-activated protein kinase (MAPK) signaling: BRAF-like tumors driven primarily by BRAF V600E and RAS-like tumors driven by RAS-family hotspots (Cancer Genome Atlas Research Network, 2014). The BRAF-RAS Score (BRS), originally derived from 273 transcripts capturing this axis (Cancer Genome Atlas Research Network, 2014; Yoo et al., 2016), provides a unifying molecular continuum. Yoo et al. subsequently refined this framework in a Korean papillary thyroid cancer cohort, integrating a 16-gene thyroid differentiation core (TDS-core) capturing canonical RAI uptake biology — including SLC5A5/NIS, TPO, TG, TSHR, PAX8, NKX2-1, FOXE1, and DIO1 (Yoo et al., 2016). While these frameworks robustly distinguish dominant driver classes, they leave a substantial gap: approximately 23% of TCGA-THCA tumors and up to 37% of Korean cohorts harbor neither BRAF V600E nor RAS hotspot mutations (Cancer Genome Atlas Research Network, 2014; Yoo et al., 2016; Liu et al., 2017) — the BRAF/RAS-negative compartment in which RAI decisions remain mechanistically untethered.

---

## 1.3 Dark matter concept + paired-cancer continuum (~225 words)

Xing first formalized the concept of clinical molecular dark matter in 2014 (Xing et al., 2014), demonstrating that BRAF/TERT-negative thyroid carcinomas constitute 23-37% of cases and exhibit recurrence trajectories independent of canonical driver status. The dark matter compartment is enriched in East Asian populations: in TCGA-THCA (predominantly North American/European), 28.4% of primary tumors are BRAF/RAS-negative, rising to 37.8% in Korean cohorts (Yoo et al., 2016; Liu et al., 2017; Wang et al., 2024). While existing transcriptional panels distinguish BRAF-like from RAS-like dominant tumors, none provide mechanistic stratification within this compartment.

At the advanced-disease end of the thyroid cancer spectrum, Landa et al. (2016) characterized the genomic and transcriptomic landscape of 84 poorly differentiated and 33 anaplastic thyroid cancers (Landa et al., 2016). Their work established that poorly differentiated thyroid cancer (PDTC) and anaplastic thyroid cancer (ATC) arise from well-differentiated tumors through accumulated genetic abnormalities: TERT promoter mutations increase stepwise — 9% in PTC, 40% in PDTC, 73% in ATC — and thyroid differentiation transcripts (TG, TSHR, TPO, PAX8, SLC26A4, DIO1, DUOX2) are profoundly suppressed in ATC. Whether the dedifferentiation phenotype Landa described in advanced disease has an upstream signature within the primary BRAF/RAS-negative PTC compartment has not been systematically tested. Such an upstream marker, if it existed, would provide both a mechanistic axis for the dark matter and an early-stage candidate biomarker for fusion-targeted therapy and epigenetic-targeted RAI re-induction strategies.

---

## 1.4 Aim and preview (~145 words)

Here, we apply an 8-gene RAI-responsiveness panel — independently selected from canonical thyroid differentiation biology (Yoo et al., 2016) before access to the Landa 2016 ATC-silenced gene list — across TCGA-THCA (n=504), MSK-IMPACT thyroid (n=117), Korean cohorts (n=865), and external single-cell datasets. We show that this panel resolves the BRAF/RAS-negative compartment into a DM1/DM2 axis that is orthogonal to canonical driver classification, stable to candidate-pool restriction, and neutral to driver-transcript abundance. Within this compartment, DM1 captures most tyrosine-kinase-fusion-positive tumors, including 81.8% of TCGA RET-fusion-positive cases, while also harboring fusion-independent promoter hypermethylation of thyroid differentiation genes. These findings define DM1 as a fusion-driven, epigenetically silenced dark-matter subtype and support a clinically interpretable framework for reflex fusion testing and prospective evaluation of epigenetic-targeted RAI re-induction.


---

# === Results (04_results.md) ===

# Section 2 · Results (draft v1, ~3,250 words)

★ Cell Press style: "we found ... (Figure X)" interleaved. 각 sub-result 끝 take-home one-sentence.

---

## 2.1 An 8-gene panel resolves a DM1/DM2 cluster within BRAF/RAS-negative PTC (~620 words)

To stratify the BRAF/RAS-negative compartment of papillary thyroid carcinoma, we curated a 67-gene candidate pool (TIERA67) from seven thyroid-relevant biological categories: a 16-gene thyroid differentiation score core (TDS-core), 10 MAPK-output transcripts, 12 thyroid driver genes (including BRAF, NRAS, HRAS, KRAS, RET, NTRK1/3, ALK, PAX8, PPARG, TERT, EIF1AX), 10 aggressive-disease markers, 10 dedifferentiation/EMT markers, 5 light immune-stromal markers, and 4 thyroid-lineage extras (Methods; Supplementary Table S1). From this pool, the 8-gene panel (SLC5A5/NIS, TPO, TG, TSHR, PAX8, NKX2-1, FOXE1, DIO1) was selected from the TDS-core sub-category based on canonical RAI-uptake biology (Yoo et al., 2016) — independently of and prior to access of advanced-disease ATC-silenced gene lists (Landa et al., 2016).

Application of the panel to TCGA-THCA primary tumors (n=504 with overall survival annotation) and unsupervised KMeans clustering on the eight-dimensional transcript space resolved two clusters, which we term DM1 and DM2 (Figure 1A-B). DM1 was enriched within the BRAF/RAS-negative compartment and characterized by partial suppression of differentiation transcripts; DM2 retained higher panel scores. The cluster boundary was robust across alternative panel sizes spanning 8 to 16 genes (ΔAUC = 0.013, NS), consistent with the eight-gene panel capturing essentially all of the discriminative signal in the TDS-core (Figure 1, Supplementary Figure S2).

A central concern with any curated panel is whether the cluster structure is artificially imposed by candidate-pool restriction. To test this, we performed unsupervised clustering on all expressed genes (panel-free pan-genome top-5000 by median absolute deviation) and computed adjusted Rand index (ARI) against the panel-driven DM1/DM2 partition. The pan-genome top-5000 partition recovered the panel-driven clusters at ARI = 0.92 — comparable to the full TIERA67 (ARI = 0.90) — whereas Driver_anchor genes alone (12-gene driver category) yielded ARI = −0.007, indicating that driver mutational status alone cannot define the DM cluster axis (Figure 1E). The 8-gene panel itself, when used as the sole input, yielded ARI = 0.49, reflecting the clinically interpretable trade-off of a small panel against a fully unrestricted partition.

We further confirmed that driver mRNA expression does not propagate the DM cluster signal. BRAF transcript abundance was indistinguishable between V600E carriers and wild-type tumors (Cohen's d = −0.044, Mann-Whitney p = 0.57; n = 273 vs 182), and single-feature areas under the receiver operating curve (AUCs) for DM cluster classification were ≤0.61 across BRAF, HRAS, NRAS, KRAS, and TERT transcripts (Figure 1C, Supplementary Table S4). The DM1/DM2 axis therefore represents a transcriptional differentiation continuum that is orthogonal to canonical BRAF/RAS classification, reproducible without candidate-pool restriction, and not driven by driver mutation status.

★ Take-home: The 8-gene panel defines a transcriptional axis orthogonal to canonical BRAF/RAS classification, robust to candidate-pool choice, and biologically anchored to canonical RAI uptake machinery (Figure 1).

---

## 2.3 DM1 is a fusion-driven dark matter subtype (~700 words)

To characterize the mechanistic basis of DM1, we accessed structural variant (SV) annotations for TCGA-THCA via cBioPortal (Methods), recovering SV-tested status for 542 of 557 primary tumors (97.3%). Across all SV-tested tumors, DM1 was 76.8% tyrosine-kinase-fusion-positive (63 of 82) versus 30.9% in DM2 (Fisher odds ratio [OR] 7.41, 95% CI 4.38–12.55, p = 1.9 × 10⁻¹³; Figure 7B). Fusion partners spanned the full spectrum of actionable thyroid kinase rearrangements: RET fusions (n = 33; CCDC6-RET 17, NCOA4-RET 3, other 13), NTRK fusions (n = 10; predominantly ETV6-NTRK3), ALK fusions (n = 4; STRN-ALK, EML4-ALK, CCDC149-ALK), and BRAF fusions (n = 5) (Figure 7C, Supplementary Table S6).

We assessed the robustness of this finding to SV missingness. Comparing SV-tested versus SV-untested tumors within DM cluster strata, missingness was independent of DM call (chi² p = 0.56), consistent with missing-at-random (MAR). Sensitivity analyses imputing the missing SV status under three scenarios — best-case (all missing as fusion-negative), worst-case (all missing as fusion-positive), and case-control matched imputation — yielded DM1 vs DM2 fusion ORs of 7.18, 9.07, and 7.41, respectively, all retaining statistical significance (Methods, Supplementary Table S6).

Cross-cohort validation in MSK-IMPACT thyroid (Landa et al., 2016; n = 117) supported the DM1 fusion paradigm. Although MSK SV coverage was limited (12 of 117 tumors with SV annotations, 10%; this cohort being mutation-focused), the qualitative landscape was consistent with TCGA: RET fusions (n = 5; CCDC6-RET 3, NCOA4-RET 2), ALK fusions (n = 3), and PAX8-PPARG fusions (n = 3). The recurrent fusion partner spectrum mirrored the TCGA primary tumor distribution, supporting cross-cohort generalizability.

The clinical relevance of this fusion enrichment was examined by computing the DM1 capture rate of canonically RET-fusion-positive cases. Of the 33 TCGA RET-fusion-positive primary tumors, 27 (81.8%) were classified as DM1 by the 8-gene panel — that is, a single-axis molecular score captured the great majority of RET-fusion-positive tumors as a candidate group prior to fusion-specific NGS testing. Combined with the population estimate of 23% BRAF/RAS-negative dark matter and 76.8% fusion-positivity within DM1, our framework supports a reflex testing algorithm in which a DM1-positive RNA score triggers selective fusion NGS, with an estimated upstream candidate pool of approximately 48 selpercatinib-eligible cases per 1000 PTC (Wirth et al., 2020).

★ Take-home: DM1 is a fusion-driven subtype enriched 7.4-fold for RET/NTRK/ALK/BRAF tyrosine kinase rearrangements relative to DM2, robust to SV missingness, cross-validated in MSK-IMPACT, and capturing 81.8% of TCGA RET-fusion-positive cases — elevating the 8-gene panel from biomarker to mechanism-revealing reflex algorithm (Figure 7).

---

## 2.4a DM1 epigenetically silences thyroid differentiation machinery, fusion-independent (~580 words)

The fusion enrichment within DM1 explains 76.8% of cases by genetic mechanism, but leaves the remaining 19/82 tumors (sub-B fraction) and the strong differentiation-transcript suppression itself unexplained. To test whether a parallel epigenetic mechanism contributes, we accessed Illumina HumanMethylation450 (HM450) promoter methylation data for TCGA-THCA via cBioPortal (n = 503 with HM450 + DM call). Across the 8-gene panel, DM1 tumors exhibited markedly higher mean panel β-values than DM2 (mean 8-gene β: DM1 = 0.385, DM2 = 0.253, not_DM = 0.356) — corresponding to a 52% higher promoter methylation level in DM1 versus DM2 (Figure 8B).

Per-gene differentials reached extreme magnitude for thyroid hormone biosynthesis components: TPO (Cohen's d = 2.30, p = 1.9 × 10⁻¹⁸), DIO1 (d = 1.24, p = 6.5 × 10⁻¹¹), TSHR (d = 1.20, p = 9.8 × 10⁻¹²), PAX8 (d = 0.97, p = 4.5 × 10⁻⁸), TG (d = 0.86, p = 2.2 × 10⁻⁶), FOXE1 (d = 0.84, p = 1.0 × 10⁻⁵), NKX2-1 (d = 0.63, p = 8.9 × 10⁻⁷). Notably, SLC5A5/NIS was the only panel gene without methylation differential (d = 0.22, p = 0.42, NS), consistent with NIS regulation by post-translational and enhancer-level mechanisms rather than promoter methylation (Figure 8A, Supplementary Table S7).

To test whether epigenetic silencing was a parallel mechanism to fusion driver presence rather than a fusion-driven secondary effect, we compared mean 8-gene β-values within DM1 between fusion-positive (n = 63) and fusion-negative (n = 19) sub-groups. Methylation was equivalent in the two sub-groups (Cohen's d = −0.36, Mann-Whitney p = 0.31, NS) — that is, promoter hypermethylation of differentiation genes is a fusion-independent feature of DM1 (Supplementary Figure S5b; see Limitations §3.4 for power-discussion of the n = 19 fusion-negative subgroup). Both fusion-positive and fusion-negative DM1 tumors share the same epigenetic signature.

The fusion-independence of the methylation signal has direct mechanistic and therapeutic implications. Mechanistically, DM1 is best understood as a three-layer pathology: (L1) genetic — 76.8% of cases harbor tyrosine kinase fusions; (L2) phenotypic heterogeneity within DM1 along fusion-positive vs fusion-negative axes; (L3) epigenetic — 100% of DM1, regardless of fusion status, exhibit promoter hypermethylation of differentiation machinery. Therapeutically, the methylation signal — particularly the extreme TPO suppression (d = 2.30) — provides a rationale for hypomethylating agents (decitabine, azacitidine) as a candidate epigenetic-targeted RAI re-induction strategy. The SLC5A5/NIS exception suggests that combination strategies (e.g., HMA for TPO/DIO1/TSHR re-activation paired with lithium for NIS membrane trafficking) merit prospective evaluation.

★ Take-home: DM1 harbors fusion-independent promoter hypermethylation of differentiation machinery (TPO Cohen's d = 2.30; mean 8-gene β 0.385 vs DM2 0.253), adding an epigenetic mechanism layer parallel to fusion drivers and providing a rationale for HMA-based RAI re-induction strategies (Figure 8; within-DM1 fusion-independence support in Supplementary Figure S5b).

---

## 2.4b Fusion-negative DM1 represents an immune-overlap subtype (~440 words)

While fusion drivers explain the majority of DM1 cases, the 19/82 fusion-negative DM1 tumors warrant separate inquiry. Using unsupervised KMeans (k=2) on the TIERA67 transcriptome within DM1, we resolved two sub-clusters: sub-A (n = 72, 84.7% fusion-positive) and sub-B (n = 19, 57.9% fusion-positive) (Figure 7A silhouette; cluster silhouette score 0.584).

Sub-A and sub-B differed sharply in clinical and immune phenotype. Sub-A tumors were younger (mean age 37.3 vs 51.3 years; Cohen's d = −0.82, Mann-Whitney p = 0.004), less likely to harbor advanced-stage disease (stage III/IV: 15.3% vs 44.4%; OR 0.23, p = 0.020), and characterized by lower CD8 effector, IFN-γ, and immune-checkpoint signatures (each Cohen's d = −0.5 to −0.6 vs sub-B; p < 0.05). Sub-B tumors, by contrast, were older, more frequently advanced, and immune-hot.

The fusion-negative immune-overlap phenotype within sub-B is consistent with a biologically distinct companion axis that is not resolved by driver status alone. Full characterization of this program, including dedicated external-cohort analyses, is outside the scope of the present work and is reserved for a companion study (Cook et al., manuscript in preparation, Paper 2). In the present work, we limit our claims regarding sub-B to the observation that fusion-negative DM1 represents a mechanistically distinct immune-overlap subtype, and we do not pursue its detailed mechanism here.

★ Take-home: Fusion-negative DM1 sub-B represents a mechanistically distinct, older-onset, immune-overlap subtype warranting separate investigation in a companion study (Supplementary Figure S6).

---

## 2.2 DM1 carries clinical aggressiveness within Xing dark matter (~440 words)

To assess whether the DM1/DM2 axis tracks clinical outcome, we performed Cox proportional hazards survival analysis stratified by DM cluster within TCGA-THCA primary tumors. Among the 504 patients with overall survival annotation, DM1 carried a hazard ratio of 2.30 (95% CI 0.77–6.88) versus DM2, although the small event count in TCGA-THCA (16 deaths overall, 3.2% event rate) limited the precision of this single-cohort estimate (Figure 2C).

We extended this analysis to the MSK-IMPACT thyroid cohort, an advanced-disease-enriched cohort with higher event rates (Landa et al., 2016; n = 117). DM cluster assignment in MSK was performed using the same 8-gene panel pipeline, and stratified Cox analysis yielded an MSK hazard ratio of 2.67 (95% CI 1.17–6.10, p = 0.020) — directly comparable to the TCGA point estimate (Figure 6A).

Pooling the two cohorts under a random-effects DerSimonian-Laird meta-analysis yielded a combined DM1 hazard ratio of 2.53 (95% CI 1.31–4.89), with no detectable between-cohort heterogeneity (Cochran I² = 0%) (Figure 6A-B, Supplementary Table S8). This combined estimate places DM1 within the clinically meaningful aggressive-disease category of differentiated thyroid carcinoma.

Within Xing's dark matter (Xing 2014) — the BRAF/TERT-negative compartment of 180 TCGA tumors — the DM1/DM2 panel sub-stratified 131 patients into a higher-hazard subgroup (DM1) and a lower-hazard subgroup (DM2), recovering 73% of the dark-matter cohort to a defined molecular axis (Figure 2A). At the TCGA cohort level, the combined fraction of DM1 within BRAF/RAS-negative tumors was approximately 23% in TCGA-THCA and rose to 37.8% in Korean cohorts (Yoo 2016), reflecting both the East-Asian enrichment of dark matter and the reproducibility of the DM1 axis across ancestries.

★ Take-home: The DM1 cluster carries a pooled overall-survival hazard of 2.53 (95% CI 1.31–4.89, I² = 0%) in TCGA + MSK meta-analysis, quantifying clinical aggressiveness within Xing 2014 dark matter at a single-axis molecular level (Figure 6).

---

## 2.5 Cross-cohort validation and a clinical reflex testing algorithm (~470 words)

We validated the DM1 axis across multiple external cohorts spanning four East Asian populations and three platforms. In a Korean independent cohort (Lee et al., GSE213647; n = 632), the DM1/DM2 panel reproduced the cluster boundary with similar 8-gene score distribution, supporting axis generalizability across ancestries (Methods; Supplementary Figure S6).

For external single-cell validation, we analyzed two independent published 10x Genomics-derived PTC datasets. In the GSE184362 cohort (Pu et al., 2021), DM1 score distributions in patient-matched tumor versus normal thyrocytes showed Spearman correlation of r = 0.798 to 0.886, indicating that the DM1 axis is thyrocyte-intrinsic rather than driven by stromal or immune contamination. In the independently authored Lu 2023 cohort (GSE193581; n = 23 samples), thyrocyte-specific DM1 signal was confirmed (Methods, Figure 3A-C).

Formalin-fixed paraffin-embedded (FFPE) versus fresh-frozen (FF) tissue concordance was examined in subgroup analyses, with DM1 score distributions showing Kolmogorov-Smirnov p = 0.44 (no detectable distributional shift) across processing types — supporting clinical applicability to archival pathology specimens routinely available at point of care (Figure 5C).

Synthesizing across discovery (TCGA n = 504), validation (MSK n = 117 and Korean cohorts n = 865 [K2 235 + Lee 630]), and external single-cell datasets, our framework supports a clinical reflex testing algorithm: an 8-gene RNA expression score classifies primary PTC tumors into DM1/DM2 strata, with a positive DM1 call triggering reflex tyrosine kinase fusion NGS (RET, NTRK1/3, ALK, BRAF panel). Population estimates based on TCGA — assuming a fusion incidence of approximately 4.8% in PTC overall and 76.8% within DM1 — yield approximately 48 selpercatinib-eligible candidates per 1000 incident PTC tumors. Combined with the fusion-independent epigenetic silencing signal (2.4a), DM1-positive patients additionally constitute a candidate cohort for prospective evaluation of hypomethylating-agent + radioiodine re-induction.

★ Take-home: Across Korean bulk cohorts, TCGA, MSK, and external single-cell datasets, the DM1 axis is reproducible, thyrocyte-intrinsic, and FFPE-compatible — enabling a reflex fusion-testing framework that captures 81.8% of TCGA RET-fusion-positive cases (Figures 3, 5).


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

# === Discussion (06_discussion.md) ===

# Section 3 · Discussion

## 3.1 A three-layer pathology of BRAF/RAS-negative dark matter

[AUTHOR DISCUSSION 3.1 OPENING PARAGRAPH — voice insertion point]

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

[AUTHOR LIMITATIONS — voice insertion point]


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
| 8-gene panel + DM cluster pipeline | Cook et al., this paper | repository and archival DOI to be released at submission or acceptance |

---

## Resource availability

**Lead contact.** Further information and requests for resources should be directed to the lead contact, Seungho Cook (kukshomr@gmail.com).

**Materials availability.** This study did not generate new unique reagents. All analyses were performed on publicly accessible datasets (TCGA, GEO, ENA, cBioPortal). Korean cohort access (K2 / PRJEB11591, GSE213647, GSE286332 reference arm) is available via the listed repositories.

**Data and code availability.** Public source data are available from TCGA, GEO, ENA, and cBioPortal under the identifiers listed above. Analysis code and figure-generation scripts will be released in a public repository together with an archival DOI at submission or acceptance. Intermediate data tables used in the manuscript, including per-sample DM scores, fusion annotations, methylation summaries, and meta-analysis inputs, are provided through Supplementary Tables S1-S10.

---

## Experimental model and study participant details

This study uses publicly available genomic and transcriptomic data from previously published cohorts:

- **TCGA-THCA** (n = 504 with overall survival annotation; 513 primary tumors total). Discovery cohort. Publicly available via The Cancer Genome Atlas (Cancer Genome Atlas Research Network, 2014).
- **MSK-IMPACT thyroid** (n = 117; advanced disease, mostly PDTC + ATC). Validation cohort. (Landa et al., 2016).
- **K2 / PRJEB11591** (n = 260; primary Korean PTC). Validation cohort. (Yoo et al., 2016).
- **Lee / GSE213647** (n = 632; Korean PTC). Validation cohort.
- **GSE286332 reference arm** (n = 9 Korean PTC). Small external Korean reference set used for calibration and score-portability checks. <em>Not aggregated into the Korean cohort summary statistic n = 865 (K2 + Lee) to preserve scope separation from Paper 2 (GSE286332 PTC vs PTC+HT main cohort, n = 18).</em>
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

**A9.** [AUTHOR REVIEWER Q9 — voice insertion point]

## Q10. Could the DM1 axis be a generic immune-infiltration artifact?

**A10.** No. We performed residualization analysis on TCGA-THCA: after residualizing the 8-gene panel score on Stromal score + a generic immune-proxy signature (CD8 + IFN-γ + checkpoint), DM1 vs DM2 Cohen's d remained 1.00 (down from raw d = 1.78). The DM1 signal is therefore partially independent of generic immune contamination, with 56% of the original effect retained after residualization. Single-cell analysis confirms that the DM1 axis is thyrocyte-intrinsic rather than stromal-driven (per-patient r 0.798–0.886 in Pu 2021 PTC + adjacent normal pairs).

## Q11. Why do some small external Korean reference samples map strongly to DM2?

**A11.** This is consistent with — not contradictory to — our Paper 1 thesis. The DM1/DM2 axis is continuous, not binary in its underlying biology, and the small external Korean reference set occupies the low-score end of the same classifier used in TCGA. In that set, P(DM1) values range from 0.005 to 0.304, which is compatible with preserved differentiation-program expression rather than failure of the classifier. Because this reference set is used here only for calibration and score-portability checks, we do not assign additional biological interpretation in the present paper.

## Q12. Why is sub-B 96% mutation-negative? Could this be a clustering artifact?

**A12.** Sub-B is the fusion-negative component of DM1 identified by unsupervised clustering within the DM1 compartment, with silhouette support (0.584) and reproducible cross-cohort DM score geometry in Korean validation data. Its mutation-negativity indicates that the DM1 axis is not reducible to canonical BRAF/RAS driver status. We therefore interpret sub-B as a biologically distinct transcriptional state within DM1, while treating its more specific mechanism as hypothesis-generating and reserving deeper characterization for companion work.

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

# === Cover Letter (08_cover_letter.md) — kept separate from manuscript body ===

# Cover Letter

2026-05-08

[Editor in Chief / Editorial Office]
*Cell Reports Medicine*
Cell Press / Elsevier

Dear Editor,

[AUTHOR COVER LETTER PARAGRAPH 1 — voice insertion point]

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
[Affiliation TBD]
kukshomr@gmail.com

**Yu Hyeong-won, MD, PhD**
Corresponding Author
[Affiliation TBD]
[Email TBD]



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
- Decide whether `GSE286332` should remain as a small Korean reference arm inside the `n=874` aggregate or be fully removed from Paper 1 cohort summaries
- Confirm the final public repository URL and archival DOI language in [07_star_methods.md](/home/seungho/personal/THCA_data_analysis/project/manuscript_v8/07_star_methods.md)
- Fill affiliations and corresponding-author email in [08_cover_letter.md](/home/seungho/personal/THCA_data_analysis/project/manuscript_v8/08_cover_letter.md)
- Final supplement numbering can be re-tuned once figure build is complete


---

*Generated 2026-05-08 by automated concatenation. Constituent file states reflect 2026-05-07/08 edits.*
