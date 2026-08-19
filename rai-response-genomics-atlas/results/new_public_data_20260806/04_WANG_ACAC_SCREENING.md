# 04 — Wang acetoacetate metabolomics: screening result

**Grade D · NOT_USABLE, and it carries a warning that matters more than the dataset itself.**

## Identification

Wang J, Xu Q, Xuan Z, Mao Y, Tang X, Yang K, Song F, Zhu X.
*Metabolomics reveals the implication of acetoacetate and ketogenic diet therapy in
radioiodine-refractory differentiated thyroid carcinoma.*
**The Oncologist** 2024, doi 10.1093/oncolo/oyae075 · PMID 38760956 · PMC11379656.

Design: untargeted LC-MS metabolomics of **tissue** from **24 RAIR-DTC versus 18 non-RAIR
DTC** patients, drawn from the tissue bank of Zhejiang Cancer Hospital and reconfirmed by two
pathologists. RAIR classification follows the 2015 ATA criteria. Downstream work is
cell-line and xenograft: γ-counter radioiodine uptake, NIS and TSHR western blot, CCK8,
colony formation, scratch/transwell, Annexin V/PI, and a ketogenic-diet xenograft arm.

## Why it is not usable

**No per-patient data of any kind is published.** Data availability states only:
*"The data underlying this article will be shared on reasonable request to the corresponding
author."* The two supplementary tables were checked: **Table S1 is the mouse ketogenic-diet
composition**, and **Table S2 is the clinical-characteristics comparison of the two groups**
(sex, tumour size, extrathyroidal extension, lymph node metastasis, TNM stage — all
non-significant between groups). Neither contains a metabolite matrix or a patient-level row.

**It is metabolomics, not transcriptomics.** As with Zheng, the eight-gene panel cannot be
tested here even if the matrix were obtained. The γ-counter iodine-uptake measurements —
the one tier-1-like quantity in the paper — are **cell-line measurements, not patient
measurements**, so they do not pair with any patient molecular profile.

## The finding that actually matters: Zheng and Wang are the same laboratory

This was the specific screening question — possible patient overlap with the Zheng cohort.
The answer is that overlap is not merely possible, it should be assumed.

| | Zheng 2024 | Wang 2024 |
|---|---|---|
| Journal | *Sci Rep* 2024, PMC11079026 | *The Oncologist* 2024, PMC11379656 |
| Cohort | **20 RAIR / 14 non-RAIR** | **24 RAIR / 18 non-RAIR** |
| Specimen | serum, collected before first surgery | tissue, from the institutional tissue bank |
| Institution | **Zhejiang Cancer Hospital, Hangzhou** | **Zhejiang Cancer Hospital, Hangzhou** |
| Lab | Key Laboratory of Head and Neck Cancer Translational Research of Zhejiang Province | **same laboratory** |
| Shared authors | **Tang Xi**, **Zhu Xin** (corresponding) | **Tang Xi**, **Zhu Xin** (corresponding) |
| Ethics approval | IRB-2020-438 Ke | Zhejiang Cancer Hospital ethics committee |
| RAIR definition | 2015 ATA criteria | 2015 ATA criteria |

Same institution, same laboratory, same senior author, same period, same refractory
definition, similar group sizes. A third paper from the same group —
Wang J et al., *Preliminary study on ketone body metabolism in anaplastic thyroid cancer*,
PMC13130877 (2026) — extends the same programme.

**Consequence: Zheng and Wang must never be entered as two independent cohorts in any
pooled analysis or evidence count.** Doing so would double-count patients and shrink a
meta-analytic confidence interval on an overlap that the papers do not disclose. Whether the
serum and tissue series are the same individuals is not stated in either paper and would have
to be asked.

This is logged as a screening exclusion rather than a candidate.

## What to do with it

Nothing analytically. Two uses remain:

1. **Cite as convergent metabolic evidence, once, with the overlap disclosed.** If the
   acetoacetate axis is mentioned in a discussion of RAI-refractory metabolism, cite Zheng and
   Wang together as one group's programme, not as two replications.
2. **If the Zheng Baidu matrix is ever obtained**, ask the same corresponding author
   (Zhu Xin) for the Wang tissue matrix in the same message — it is one request, not two.

## Screening verdict

`WANG_ACAC_2024 · EXCLUDE · no patient-level data published; metabolomics cannot test the
panel; non-independent from ZHENG_2024 (same laboratory and senior author).`
