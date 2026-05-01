---
title: "Paper 1 manuscript v8 — Section 1 Introduction (draft v1)"
date: 2026-05-01
author: Seungho Cook
target_venue: Cell Reports Medicine (1순위) → JCI Insight + Nat Commun dual stretch → npj Precision Oncology (fallback)
target_words: 600-900
draft_words: 715
status: v1 draft — 본인 voice 영역 명시. 수치/cite placeholder 5/4 본인 read 후 verify.
voice_protected_areas: 1.1 첫 줄 (Hook), 1.4 (Aim + preview)
---

# Section 1 · Introduction (draft v1, ~715 words)

## 1.1 Clinical context (~180 words)

> ★ **Hook 첫 줄 — 본인 voice 영역.** 아래 Alt B refined 가 default. Hybrid "compass" 옵션 본인 결정.

[**Alt B refined / default**] Despite a >98% 5-year overall survival in differentiated thyroid carcinoma, structural disease recurrence reaches ~20% in ATA 2015 intermediate-risk patients (Haugen et al., 2016), and BRAF/RAS-negative tumors — accounting for ~23% of cases — complicate radioiodine (RAI) treatment decisions in the absence of mechanistic sub-stratification.

[**Hybrid "compass" 옵션**] Despite a >98% 5-year overall survival, structural disease recurrence reaches ~20% in ATA 2015 intermediate-risk papillary thyroid carcinoma (Haugen et al., 2016), and clinicians continue to make radioiodine decisions on heterogeneous BRAF/RAS-negative tumors — ~23% of cases — without a mechanistic compass.

Papillary thyroid carcinoma (PTC) is the most common endocrine malignancy and among the fastest-rising in incidence over the past three decades [SEER cite — verify]. While the majority of patients achieve durable remission after thyroidectomy and selective ¹³¹I ablation, 5-20% develop recurrent or persistent disease and approximately 10% develop distant metastases (Haugen et al., 2016). Current risk-tier-based RAI decisions — codified in the 2015 American Thyroid Association (ATA) Management Guidelines (Haugen et al., 2016) and recently updated as ATA 2025 (Ringel et al., 2025) — rely predominantly on clinico-pathological features (tumor size, multifocality, extrathyroidal extension, lymph node burden) with BRAF V600E as the sole molecular risk modifier. Bethesda III/IV indeterminate cytology affects 15-30% of fine-needle aspiration biopsies and remains a major diagnostic gap [Cibas 2017 — verify]. Together, these gaps motivate orthogonal molecular sub-stratification of clinically heterogeneous tumors, particularly within the BRAF/RAS-negative compartment.

---

## 1.2 Existing molecular framework (~165 words)

The 2014 Cancer Genome Atlas (TCGA) study of papillary thyroid carcinoma established a binary molecular spectrum anchored by mitogen-activated protein kinase (MAPK) signaling: BRAF-like tumors driven primarily by BRAF V600E and RAS-like tumors driven by RAS-family hotspots (Cancer Genome Atlas Research Network, 2014). The BRAF-RAS Score (BRS), originally derived from 273 transcripts capturing this axis (Cancer Genome Atlas Research Network, 2014; Yoo et al., 2016), provides a unifying molecular continuum. Yoo et al. subsequently refined this framework in a Korean papillary thyroid cancer cohort, integrating a 16-gene thyroid differentiation core (TDS-core) capturing canonical RAI uptake biology — including SLC5A5/NIS, TPO, TG, TSHR, PAX8, NKX2-1, FOXE1, and DIO1 (Yoo et al., 2016). While these frameworks robustly distinguish dominant driver classes, they leave a substantial gap: approximately 23% of TCGA-THCA tumors and up to 37% of Korean cohorts harbor neither BRAF V600E nor RAS hotspot mutations (Cancer Genome Atlas Research Network, 2014; Yoo et al., 2016; Liu et al., 2017) — the BRAF/RAS-negative compartment in which RAI decisions remain mechanistically untethered.

---

## 1.3 Dark matter concept + paired-cancer continuum (~225 words)

Xing first formalized the concept of clinical molecular dark matter in 2014 (Xing et al., 2014), demonstrating that BRAF/TERT-negative thyroid carcinomas constitute 23-37% of cases and exhibit recurrence trajectories independent of canonical driver status. The dark matter compartment is enriched in East Asian populations: in TCGA-THCA (predominantly North American/European), 28.4% of primary tumors are BRAF/RAS-negative, rising to 37.8% in Korean cohorts (Yoo et al., 2016; Liu et al., 2017; Wang et al., 2024). While existing transcriptional panels distinguish BRAF-like from RAS-like dominant tumors, none provide mechanistic stratification within this compartment.

At the advanced-disease end of the thyroid cancer spectrum, Landa et al. (2016) characterized the genomic and transcriptomic landscape of 84 poorly differentiated and 33 anaplastic thyroid cancers (Landa et al., 2016). Their work established that poorly differentiated thyroid cancer (PDTC) and anaplastic thyroid cancer (ATC) arise from well-differentiated tumors through accumulated genetic abnormalities: TERT promoter mutations increase stepwise — 9% in PTC, 40% in PDTC, 73% in ATC — and thyroid differentiation transcripts (TG, TSHR, TPO, PAX8, SLC26A4, DIO1, DUOX2) are profoundly suppressed in ATC. Whether the dedifferentiation phenotype Landa described in advanced disease has an upstream signature within the primary BRAF/RAS-negative PTC compartment has not been systematically tested. Such an upstream marker, if it existed, would provide both a mechanistic axis for the dark matter and an early-stage candidate biomarker for fusion-targeted therapy and epigenetic-targeted RAI re-induction strategies.

---

## 1.4 Aim and preview (~145 words)

> ★ **본인 voice 영역.** Paper 의 identity 를 정의하는 paragraph — 본인이 voice 다시 적기 권장.

Here, we introduce an 8-gene RAI-responsiveness panel — independently selected from canonical RAI uptake biology (Yoo et al., 2016) prior to access of the Landa 2016 ATC-silenced gene list — and apply it to TCGA-THCA (n=504), MSK-IMPACT thyroid (n=117), three Korean cohorts (K2, Lee, GSE286332; n=874), and external single-cell datasets (Pu et al., 2021; Lu et al., 2023). The panel resolves a DM1/DM2 split that is orthogonal to canonical BRAF/RAS classification, stable to candidate-pool restriction (pan-genome top-5000 ARI 0.92 vs panel-driven TIERA67 0.90), and driver-mRNA-neutral (BRAF transcript Cohen's d=−0.04 in V600E carriers vs wild-type). Within the BRAF/RAS-negative compartment, DM1 captures 76.8% of tyrosine-kinase-fusion-positive tumors (RET/NTRK/ALK/BRAF), 81.8% of TCGA RET-fusion-positive cases, and harbors fusion-independent promoter hypermethylation of differentiation genes — together defining a fusion-driven, epigenetically silenced subtype amenable to reflex fusion testing and rationale-supported epigenetic-targeted RAI re-induction.

---

## Word-count breakdown

| Section | Target | Draft v1 |
|---|---|---|
| 1.1 Clinical context | 150-200 | 180 |
| 1.2 Existing molecular framework | 150-200 | 165 |
| 1.3 Dark matter + paired continuum | 200-250 | 225 |
| 1.4 Aim + preview | 100-150 | 145 |
| **Total** | **600-900** | **715** ✅ |

## Cite placeholder verify list (5/4 본인 read 후)

| Cite | 위치 | 검증 필요 |
|---|---|---|
| SEER cancer statistics (PTC incidence rise 30 years) | 1.1 | ⚠️ verify — alternate cite: Davies & Welch 2014 *JAMA*, Lim et al. 2017 *JAMA Oncol* |
| Cibas & Ali 2017 *Thyroid* (Bethesda III/IV 15-30%) | 1.1 | ⚠️ verify — *Bethesda system 2017 update* |
| Haugen 2016 ATA 2015 (5-20% recurrence, 10% distant) | 1.1 | ✅ confirmed (PMID 26462967) — 본인 PDF Tables 11-15 verbatim 검증 |
| Ringel 2025 ATA 2025 | 1.1 | ✅ confirmed (PMID 40844370) |
| Cancer Genome Atlas Research Network 2014 *Cell* | 1.2 | ✅ confirmed (PMID 25417114) — 본인 PMC4243044 read |
| Yoo 2016 *PLOS Genet* | 1.2, 1.4 | ✅ confirmed (PMID 27494611) — 본인 read |
| Liu et al. 2017 (East Asian) | 1.2, 1.3 | ⚠️ verify — Liu Z et al. 2017 cohort 정확 cite |
| Wang et al. 2024 (Shanghai) | 1.3 | ⚠️ verify — Wang YL Endocr Connect 13(11):e240301 |
| Xing et al. 2014 (dark matter) | 1.3 | ⚠️ verify — Xing M et al. *NEJM* 2014 또는 *JCO* 2015 정확 cite |
| Landa et al. 2016 *JCI* | 1.3 | ✅ confirmed (PMID 26878173) — 본인 read |
| Pu et al. 2021 *Nat Commun* | 1.4 | ✅ confirmed (PMID 34663816) |
| Lu et al. 2023 | 1.4 | ⚠️ verify — Lu et al. 2023 *Cancer Cell* 또는 *Nat Cancer* 정확 cite |

## 본인 voice 적용 영역 (W1 5/5-5/10)

- [ ] **1.1 첫 줄 hook**: Alt B refined vs hybrid "compass" 1개 선택, 본인 voice 미세 조정
- [ ] **1.4 aim paragraph**: Paper identity — "fusion-driven, epigenetically silenced" 표현 본인 voice 검증, "amenable to" 표현 톤 (vs "supports", "motivates", "enables") 본인 결정
- [ ] **1.3 마지막 줄**: "Such an upstream marker..." → 본 paper 의 thesis statement 직전 line — Cell Press style 검증
- [ ] **Section 간 transitions** — 본인 voice 자연스러운 흐름 검증

## 다음 step

1. 본인 5/4 read (Landa, Yoo, ATA, TCGA, Pu) → cite verify
2. 본인 voice 적용 (1.1 hook + 1.4 aim)
3. v2 draft → Web Claude review → v3 final
4. Prompt 3 (Results 2.1-2.5) 진입 (W2-W3)
