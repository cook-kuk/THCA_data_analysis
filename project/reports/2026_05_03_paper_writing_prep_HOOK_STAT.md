# Hook Stat Fact-Check — "5-15% recur/metastasize"

**Date:** 2026-05-03 (5/4 prep, Claude Code 30 min fact-check)
**Source:** ATA 2015 Guidelines (Haugen et al. 2016, Thyroid 26:1-133, PMC4739132) + JCEM 2019 review (Tuttle/Alzahrani, oup.com/jcem/104/9/4087)

---

## ★ Critical finding — Hook stat needs revision

**본인 brief 의 "5-15% recur/metastasize" 수치는 부정확.** 정확 stat:

### Overall PTC recurrence (ATA 2015 Haugen 2016)

> "Overall recurrence rates for PTC are **15–35%**, with tumor recurrence typically occurring in the tumor bed, cervical lymph nodes, or (rarely) distant sites."

→ **15-35% NOT 5-15%.**

### Microcarcinoma (PTMC, ≤1 cm) — possibly the source of "5-15%" misremembering

| Outcome | Microcarcinoma rate |
|---|---|
| Disease-specific mortality | <1% |
| Locoregional recurrence | 2-6% |
| Distant recurrence | 1-2% |

→ Microcarcinoma combined recurrence could be ~3-8%, but this is for ≤1 cm only — not all PTC.

---

## Recurrence rates by ATA 2015 risk tier (JCEM 2019 review)

Note: ATA 2019 review reports **response-to-therapy outcomes** (not traditional recurrence) since dynamic risk stratification became standard:

### Low-Risk
- No evidence of disease (NED) at follow-up: **80-90%**
- Biochemical incomplete response: ~15%
- **Structural incomplete response: 3-5%**

### Intermediate-Risk
- Excellent response: ~60%
- Biochemical incomplete: 15-20%
- **Structural incomplete: ~20%**

### High-Risk
- Excellent response: <30%
- Biochemical incomplete: 10-15%
- **Structural incomplete: 50-75%**

→ Traditional ATA 2015 estimates approximately mapped: **low ~3%, intermediate ~20%, high ~50-75% structural recurrence.**

---

## Hook revision options for paper

### Option 1 — Use ATA 2015 overall PTC rate (most defensible)

> "Although 5-year overall survival exceeds 98%, **15-35% of patients with papillary thyroid carcinoma develop recurrent or metastatic disease** (American Thyroid Association 2015 Guidelines¹), and current risk stratification offers no mechanistic compass for predicting which differentiated tumors will follow this trajectory."

### Option 2 — Use stratified rate (more specific, more defensible)

> "Despite excellent overall survival, **structural disease recurrence ranges from 3-5% in ATA 2015 low-risk to 50-75% in high-risk patients**¹, and current risk stratification systems do not clarify the molecular axis underlying this heterogeneity."

### Option 3 — Stay with "5-15%" but reframe as range estimation

> "Approximately **5-30% of differentiated thyroid carcinomas recur or metastasize**¹·², an order-of-magnitude variability that current anatomic risk stratification only partially explains."

→ **권고: Option 1 + 2 hybrid.** "15-35% recurrence" 가 ATA published exact stat; if want lower-bound for less-aggressive cohort, "3-5% (low-risk) to 50-75% (high-risk)" 사용. **5-15% 단일 수치는 source 불명확하므로 제거 권장.**

---

## Citation

```bibtex
@article{Haugen2016ATA,
  author = {Haugen, Bryan R and Alexander, Erik K and Bible, Keith C and others},
  title = {2015 American Thyroid Association Management Guidelines for Adult Patients with Thyroid Nodules and Differentiated Thyroid Cancer},
  journal = {Thyroid},
  volume = {26},
  number = {1},
  pages = {1-133},
  year = {2016},
  doi = {10.1089/thy.2015.0020},
  pmid = {26462967},
  pmcid = {PMC4739132}
}

@article{Tuttle2019,
  author = {Tuttle, R Michael and Alzahrani, Ali S},
  title = {Risk Stratification in Differentiated Thyroid Cancer: From Detection to Final Follow-Up},
  journal = {Journal of Clinical Endocrinology and Metabolism},
  volume = {104},
  number = {9},
  pages = {4087-4100},
  year = {2019},
  doi = {10.1210/jc.2019-00177}
}
```

---

## 본인 voice 영역 (Hook 첫 단락 표현 — Claude Code 안 함)

| Phrase | Tone |
|---|---|
| "no mechanistic compass" | Dramatic, distinctive |
| "without a mechanistic basis" | Safe, slightly dry |
| "without molecular guidance" | Middle ground |
| "lacks a transcriptional axis to predict" | Specific to our claim |

→ 본인이 manuscript v2.5 reviewer-tone preference 따라 결정.

---

## Quality check ✅
- ATA 2015 source: Haugen et al. 2016 Thyroid PMC4739132 verified
- Stat 15-35% PTC recurrence: confirmed via PMC fetch
- Microcarcinoma 2-6% locoregional / 1-2% distant: confirmed
- ATA 2015 modification (low risk includes <0.2 cm LN, <4 vascular foci, BRAFV600E microcarcinoma): confirmed
- 5-15% original brief stat: **NOT directly attributable to ATA 2015**, source unclear
