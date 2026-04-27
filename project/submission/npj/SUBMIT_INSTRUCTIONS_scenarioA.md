# Scenario A · npj Precision Oncology submit instructions

**Status**: Submit-ready (post-ULTIMATE sprint)
**Date**: 2026-04-27 KST
**Decision**: Scenario A · npj Precision Oncology 1차
**ETA from now to click**: 1 week (Yu confirm 후)

---

## Pre-flight

- [ ] **Yu (공동저자) confirm** — 5 결정 항목 (`reports/v17_ultimate/yu_decision_5_items.md`)
- [ ] **Park YJ outreach** 발송 (timing per Yu decision)
- [ ] **분당서울대 outreach v2** 발송
- [ ] **Manuscript v6 ULTIMATE final read-through** (~2-3시간)
- [ ] **Figure 정합성** — 모든 inline 통계 + figure caption 일치
- [ ] **References list** 완성도 (Han 2023 PMID 37461149, Pu 2021 PMID 34663803, Pu 2022 PMID 36202955, Yang 2022 PMID lookup, Mu 2024 JCEM citation)
- [ ] **Author affiliations + ORCID** — 본인 + Yu

---

## Submit URL + account

- **Journal**: npj Precision Oncology
- **Submit portal**: https://mts-npjprecisiononcology.nature.com/
- **Account**: 본인 ORCID 등록 (없으면 https://orcid.org/register)
- **Article type**: Article (full primary research)
- **Word limit**: 4,000-5,000 (current v6 ULTIMATE = ~3,800)
- **Figure limit**: 8 main + unlimited supp (current = 8 main + 15 supp + 8 R-figures = 23 supp 가능 충분)
- **Reference limit**: ~50 (current ≈ 40 with Korean parallel work additions)

---

## Document checklist (5 PDFs to upload)

| # | File | Description | Status |
|---|---|---|---|
| 1 | `manuscript_v6_ULTIMATE.pdf` | Main manuscript (compile from .md → docx → pdf via pandoc) | ☐ pending compile |
| 2 | `cover_letter_v5.pdf` | Cover letter (3 시나리오 통합 honest reframe) | ☐ to write |
| 3 | `figures.zip` | 8 main + 23 supp PNG (300 DPI) + PDF (vector) | ✓ ready |
| 4 | `Supplementary_Tables.xlsx` | 7 sheets (R2 honest AUC + R1 concordance + R3 RAI + WHO 2022 + cBioPortal S20 + ...) | ☐ to compile |
| 5 | `anonymous_code.zip` | Anonymous code bundle (~2MB) — full reproducibility | ✓ ready |

추가 (요청 시):
- `competing_interests_statement.pdf` — 빈 (no COI)
- `data_availability_statement.pdf` — 모두 open access cohorts; Mu 2024 controlled

---

## Suggested reviewers (5)

| # | Reviewer | Affiliation | Rationale |
|---|---|---|---|
| 1 | **Park YJ** | SNU (Korea) | Han SC 2023 senior author · Korean parallel work · 본인 outreach 진행 중 |
| 2 | **Pu W** | Sichuan Univ (China) | scRNA BRAF-like-B subtype (Nat Commun 2021) + 4-subtype framework (Oncogene 2022) |
| 3 | **Landa I** | Memorial Sloan Kettering | BRAF-like vs RAS-like axis 원 정의 (JCI 2016) |
| 4 | **Fagin JA** | Memorial Sloan Kettering | RAI refractoriness clinical mechanism · NIS biology |
| 5 | **Xing M** | Mayo Clinic | TERT promoter prognosis 갑상선 (BRAF/TERT framework) |

**자제 (excluded)**: 동일 institution · 본인 또는 Yu 직간접 협력 관계 있는 lab

---

## Cover letter v5 — 핵심 단락

**작성 priorities**:

1. **Pre-submission self-audit + correction** explicit
   > "Following pre-submission self-audit, we identified a circular validation issue in our v5 preprint and corrected it through four independent leak-free cluster reconstructions. The reported AUC values (0.92-0.96 range across three leak-free constructions) now reflect honest predictive performance; the v5 preprint's 0.954 figure is no longer claimed."

2. **Cluster orientation correction** transparent
   > "The de-circularization pipeline also surfaced a labelling correction: the v5 preprint's 'DM1 = aggressive' convention is reversed in the leak-free re-derivation, where DM1 corresponds to the well-differentiated / RAS-like cluster and DM2 to the dedifferentiated / BRAF-like cluster. Four independent lines of evidence confirm this corrected orientation (WHO 2022 mapping, Pu 2021 direction-concordance, methylation cross-modality, GSE213647 effect direction); the underlying biology and effect sizes are preserved."

3. **★ RAI clinical validation** highlight
   > "We provide the first cross-validation of a deployable RAI panel against directly-labelled refractoriness data (GSE151179, n=39, AUC 0.671 [0.514, 0.946]). This complements the histology-based external validations (GSE76039 ATC vs PDTC AUC 0.935; 7-cohort pooled AUC 0.898 [0.835, 0.961])."

4. **WHO 2022 alignment** as positioning
   > "The eight-gene panel operationalises the WHO 2022 morphologic axis molecularly (concordance OR = 20.4, p = 2.5×10⁻³³, n=476 evaluable). This is the deployable minimal version of the transcriptomic axis previously described by Han et al. (ENM 2023), Pu et al. (Nat Commun 2021, BRAF-like-B subtype), and Pu et al. (Oncogene 2022, four-subtype framework) — extending these by adding decision-curve clinical net benefit, cross-modality methylation validation, and seven-cohort meta-analysis with true RAI clinical ground truth."

5. **Korean parallel work outreach** acknowledged
   > "Collaborative validation with the Han / Park YJ group (SNU) and Bundang SNU on Korean patient RNA-seq and RAI uptake clinical scoring is in progress (revision-round commitment). The Korean BRAF V600E rate (~62%) and Yang et al. 2022 (n=2,092) TERT prevalence (2.8% PTC) provide the population-level baselines for that validation."

6. **Limitations** — short, honest paragraph

---

## After submit — what happens

- **Review timeline**: npj Precision Oncology median review 6-8 weeks
- **Likely outcome**:
  - **Strong accept** (P 25%): minor revisions only
  - **Major revision** (P 60%): expected — Park YJ Korean cohort validation + RAI ground truth strengthening + cluster convention reframe
  - **Reject + transfer** (P 15%): BMC Cancer 또는 Frontiers Oncology가 fall-back
- **Revision plan**: GSE213647 추가 분석 + Korean cohort validation (Park YJ outreach success 시) + Mu 2024 DAC 신청 결과
- **Do not** add new claims during revision; only deepen existing

---

## Alternative venues (if reject)

| Venue | Fit | ETA submit-after-reject | Word limit |
|---|---|---|---|
| BMC Cancer | High (broader scope OK) | 1 week | 4,500 |
| Frontiers Oncology | Medium-high (open access) | 1 week | 6,000 |
| Bioinformatics | Methodology pivot (DIAL framework) | 2 weeks (rewrite) | 4,000 |
| Cancer Research | Aim higher (long shot) | not recommended | 5,000 |

---

## SUBMIT 직전 마지막 체크 (1 hour before click)

- [ ] PDF 모든 figure 300 DPI (실측)
- [ ] manuscript word count 5,000 미만
- [ ] 모든 inline 통계 supplementary table에 매칭됨
- [ ] References 모두 PMID 또는 DOI 표기 정확
- [ ] Code repo URL 작동 (또는 anonymous bundle 첨부)
- [ ] Cover letter Yu confirm 받은 최종 버전
- [ ] Submission portal 5 PDF 모두 업로드 + suggested reviewers 입력
- [ ] **호흡 한 번 + click submit**
