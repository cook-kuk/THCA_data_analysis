---
title: "Paper 1 manuscript v8 — Section 1 Introduction (draft v1)"
date: 2026-05-04
author: Seungho Cook
target_venue: Cell Reports Medicine (1순위) → JCI Insight + Nat Commun dual stretch → npj Precision Oncology (fallback)
target_words: 600-900
draft_words: 715
status: clean draft
---

# Section 1 · Introduction

## 1.1 Clinical context (~180 words)

[AUTHOR HOOK — voice insertion point; see `04_intro_1_1_hook.md`]

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

Here, we apply an 8-gene RAI-responsiveness panel — independently selected from canonical thyroid differentiation biology (Yoo et al., 2016) before access to the Landa 2016 ATC-silenced gene list — across TCGA-THCA (n=504), MSK-IMPACT thyroid (n=117), Korean cohorts (n=874), and external single-cell datasets. We show that this panel resolves the BRAF/RAS-negative compartment into a DM1/DM2 axis that is orthogonal to canonical driver classification, stable to candidate-pool restriction, and neutral to driver-transcript abundance. Within this compartment, DM1 captures most tyrosine-kinase-fusion-positive tumors, including 81.8% of TCGA RET-fusion-positive cases, while also harboring fusion-independent promoter hypermethylation of thyroid differentiation genes. These findings define DM1 as a fusion-driven, epigenetically silenced dark-matter subtype and support a clinically interpretable framework for reflex fusion testing and prospective evaluation of epigenetic-targeted RAI re-induction.
