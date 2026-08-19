---
title: Author Voice Approval Record — DM1 Paper 1
date: 2026-07-30
author: Seungho Cook
status: ALL 6 SLOTS APPROVED — 2026-07-30 — inserted into manuscript
rule: Do not merge any slot into the manuscript until the author-approved? column is YES and approval_date is filled.
---

# Author Voice Approval Record

아래 6개 섹션은 Seungho Cook의 저자 목소리가 필요한 voice-protected 구역이다.
현재 문서는 기존 연구 방향, 미팅 결론, audit-locked 결과를 바탕으로 **저자 검토용 초안**을 채운 상태다.

**중요:** 아래 영문은 아직 승인본이 아니다. 저자 확인 전 원고에 병합하지 않는다.

---

## V1 — Introduction Opening Hook

**파일**: `NATURE_COMMUNICATIONS_FULL_DRAFT_v2_2026_07_08.md` line ~44

### 승호님 raw input

가장 억울한 지점은 RAI가 실제로 듣지 않을 가능성이 높은 환자에게도, 그 사실을 미리 알 수 없다는 이유로 고용량 방사성요오드를 반복 투여한다는 것이다. 갑상선암 환자는 대체로 오래 생존하기 때문에, 치료 효과가 불확실한 상태에서 누적되는 독성은 단순한 부작용이 아니라 장기간 남는 치료 부담이다.

첫 세 문장의 중심 톤은 **예방 가능한 피해와 진단적 불확실성**이다. 과도한 긴급성을 만들기보다, 현재의 임상 의사결정이 분자적 근거 없이 경험적으로 이뤄지는 문제를 보여주고 싶다.

독성은 타액선 손상, 골수 억제, 백혈구 감소, 장기적으로는 이차 악성종양 위험을 언급할 수 있다. 다만 모든 환자에게 흔하거나 확정적인 결과처럼 쓰지 않고, 반복적 또는 고용량 노출에서 고려되는 누적 위험으로 제한해야 한다.

### Claude polished draft A — Restrained NC style

Most patients with differentiated thyroid cancer survive for many years, yet a clinically important subset receives repeated high-dose radioiodine without a reliable means of determining, before treatment, whether their tumours retain the molecular machinery required for iodine handling. For these patients, the absence of a pretreatment response biomarker can lead to additional exposure despite limited expected benefit, with cumulative risks that include salivary-gland injury, bone-marrow toxicity and secondary malignancy. A molecular state that identifies loss of thyroid-lineage and iodine-handling capacity could therefore help distinguish patients for whom repeated radioiodine is biologically plausible from those in whom further empiric treatment may add harm without commensurate benefit.

### Claude polished draft B — Direct authorial style

The central clinical problem is not simply that some thyroid cancers fail to respond to radioiodine, but that this failure is often recognized only after patients have already received repeated treatment. Because many patients with differentiated thyroid cancer live for decades, avoidable exposure to ineffective high-dose radioiodine can impose lasting salivary, haematologic and secondary-cancer risks. We therefore asked whether the molecular loss of thyroid-lineage function could be detected before repeated treatment and used to identify tumours in which further radioiodine is unlikely to be biologically productive.

### Claude-added content note

- "Clinically important subset"은 유병률을 직접 제시하지 않고 임상 중요성만 표현함.
- "Biologically plausible"과 "may add harm"으로 직접 치료 효용 입증처럼 과장하지 않음.
- 직접적인 pretreatment RAI response 데이터가 없으므로 "predict non-response" 대신 "identify loss of capacity"로 제한함.

| 항목 | 내용 |
|---|---|
| author-approved? | **YES — 2026-07-30** |
| preferred draft | **A** |
| approval_date | 2026-07-30 |
| inserted_to | `NATURE_COMMUNICATIONS_FULL_DRAFT_v2_2026_07_08.md` §Introduction ¶1 |

---

## V2 — Introduction Final Aim Paragraph

**파일**: `NATURE_COMMUNICATIONS_FULL_DRAFT_v2_2026_07_08.md` line ~52

### 승호님 raw input

이 논문의 1차 기여는 **driver와 독립적인 thyroid-lineage state를 발견하고 외부 코호트에서 재현한 것**이다. 임상 결과 층화와 IHC 3-plex는 중요한 번역적 확장이지만, single-cohort interaction과 in-silico proxy이므로 논문의 중심 발견보다 앞에 두면 과장이 된다.

"In this study, we…" 이후에는 `define and validate a thyroid-lineage state axis`가 가장 정확하다. 이어서 이 state가 methylation, fusion/MAPK biology, RAI-refractory tumour expression과 어떻게 정렬되는지 보여주고, BRAF V600E-positive 환자에서의 PFI interaction은 hypothesis-generating clinical extension으로 배치한다.

### Claude polished draft A — Restrained NC style

In this study, we defined and validated an eight-gene thyroid-lineage state axis across transcriptomic, epigenomic, proteomic and single-cell datasets. Using TCGA-THCA as the discovery cohort, we identified iodine-handling-low and iodine-handling-high states that were robust to genome-wide feature selection and largely orthogonal to individual oncogenic drivers. We then examined their association with kinase fusions, promoter methylation, external thyroid-cancer cohorts and post-radioiodine refractory tumours. Finally, we assessed whether this lineage state stratified clinical outcomes within molecularly defined subgroups and explored a three-marker expression proxy as a potential path toward tissue-based validation. Together, these analyses establish a reproducible lineage-state framework and generate a testable hypothesis for avoiding repeated radioiodine in tumours with diminished iodine-handling biology.

### Claude polished draft B — Direct authorial style

Here, we define a reproducible thyroid-lineage state that is not reducible to BRAF, RAS or kinase-fusion status. Starting from TCGA-THCA, we identify an eight-gene axis that separates iodine-handling-low and iodine-handling-high tumours, trace this state across methylation, bulk RNA, single-cell and proteomic datasets, and show concordance with tumours that became refractory after radioiodine. We further show that the state stratifies progression-related outcomes specifically in BRAF V600E-positive disease and can be approximated in silico using three routinely assessed thyroid-lineage proteins. These findings provide a biological framework for prospective testing of whether lineage-state assessment can reduce empiric repetition of radioiodine in patients unlikely to benefit.

### Claude-added content note

- 1차 novelty를 "state discovery and validation"으로 고정함.
- BRAF+ PFI는 "specifically" / "within subgroup"으로 기술하되 predictive biomarker 확정 안 함.
- IHC는 "approximated in silico"로 한정함.
- 마지막 문장은 임상 적용이 아닌 prospective testing 필요성 명시.

| 항목 | 내용 |
|---|---|
| author-approved? | **YES — 2026-07-30** |
| preferred draft | **A** |
| approval_date | 2026-07-30 |
| inserted_to | `NATURE_COMMUNICATIONS_FULL_DRAFT_v2_2026_07_08.md` §Introduction 마지막 문단 |

---

## V3 — Discussion §3.1 Mechanism Interpretation

**파일**: `NATURE_COMMUNICATIONS_FULL_DRAFT_v2_2026_07_08.md` line ~134

### 승호님 raw input

가장 타당한 설명은 강한 MAPK signalling이 유지되는 종양에서 thyroid-lineage programme이 억제되고, 그 상태가 promoter methylation과 함께 안정화된다는 모델이다. BRAF와 RET fusion에서 평균 methylation이 높고 RAS에서 낮다는 점은 이 해석과 일치한다. 그러나 현재 자료만으로 MAPK activation → methylation → gene silencing이라는 방향성과 인과 순서를 증명할 수는 없다.

다른 plausible mechanism으로는 lineage transcription-factor network의 붕괴, enhancer 접근성 변화, repressive chromatin regulator의 동원, tumour evolution 과정에서의 state selection이 있다. methylation은 원인이라기보다 이미 형성된 탈분화 상태를 유지하거나 반영하는 현상일 수도 있다.

이 모델을 포기하거나 크게 수정했을 결과: (1) driver 보정 후 DM1/DM2 차이가 사라지는 경우, (2) promoter methylation과 RNA 발현이 같은 방향으로 움직이거나 무관한 경우, (3) 독립 코호트에서 lineage axis가 재현되지 않는 경우, (4) RAI-refractory 종양이 DM1이 아니라 DM2에 정렬되는 경우, (5) fusion-negative DM1이 존재하지 않고 state가 fusion status로 완전히 설명되는 경우.

### Claude polished draft A — Restrained NC style

The DM1 state is consistent with a model in which sustained MAPK-pathway activity accompanies loss of the thyroid-lineage programme and its stabilization through epigenetic repression. The higher mean methylation observed in BRAF V600E- and RET-fusion-positive tumours than in RAS-mutant tumours supports this interpretation, as does the inverse relationship between promoter methylation and expression for several iodine-handling genes. These findings extend the dedifferentiation programme described by Landa and colleagues in poorly differentiated and anaplastic thyroid cancers, with five of the eight genes in our panel overlapping the lineage genes suppressed in that study. However, our data do not establish a causal sequence from oncogenic signalling to DNA methylation and transcriptional silencing. Methylation may contribute to maintenance of the state, arise downstream of lineage-factor loss, or mark a broader repressive chromatin configuration. The persistence of DM1 across driver classes and among fusion-negative tumours further argues that the state is not equivalent to a single genomic alteration. Thus, our results motivate, rather than confirm, a mechanistic model linking MAPK signalling, epigenetic stabilization and loss of iodine-handling capacity.

### Claude polished draft B — Direct authorial style

We interpret DM1 as a stable loss of thyroid-lineage identity rather than as a transcriptional proxy for any single driver. Its enrichment in BRAF V600E- and RET-fusion-positive tumours, together with higher promoter methylation in these groups than in RAS-mutant disease, is compatible with a model in which strong MAPK signalling favours lineage collapse and epigenetic stabilization. The overlap with the dedifferentiation programme reported by Landa et al. in poorly differentiated and anaplastic thyroid cancer places DM1 on a biologically credible continuum of thyroid-cell identity loss. Nevertheless, the present data cannot determine whether methylation initiates silencing, maintains an already repressed state, or merely records it. Alternative mechanisms — including disruption of the PAX8–NKX2-1–FOXE1 transcriptional network, enhancer remodelling and selection of pre-existing low-lineage states — remain plausible. We therefore regard methylation as an epigenetic correlate that makes the state experimentally actionable, not as a demonstrated causal mechanism.

### Claude-added content note

- **Landa 2016 JCI 126(3):1052-1066** 정확히 인용. Krishnamoorthy 2025 오인용 제거.
- DNMT3A/3B, SETDB1은 직접 분석 없으므로 본문에서 제거하고 reviewer response 또는 future work에만 허용.
- HMA-based re-induction은 "experimentally testable" 수준으로만 언급.

| 항목 | 내용 |
|---|---|
| author-approved? | **YES — 2026-07-30** |
| preferred draft | **B** |
| approval_date | 2026-07-30 |
| inserted_to | `NATURE_COMMUNICATIONS_FULL_DRAFT_v2_2026_07_08.md` §Discussion 3.1 |

---

## V4 — Limitations

**파일**: `NATURE_COMMUNICATIONS_FULL_DRAFT_v2_2026_07_08.md` line ~154

### 승호님 raw input

임상 해석을 가장 크게 제한하는 것은 **치료 전 RAI 흡수 또는 반응과 연결된 독립 환자 코호트가 없다는 점**이다. DM1이 iodine-handling-low biology와 post-RAI refractory state에 정렬된다는 것은 보였지만, 이 패널이 개별 환자의 RAI 반응을 치료 전에 예측한다고 결론 내릴 수 없다.

두 번째 핵심 한계는 BRAF+ PFI interaction이 TCGA 단일 코호트 후향 결과라는 점이다. predictive biomarker 규정을 위해서는 독립적인 treatment-linked cohort에서 재현이 필요하다.

첫 번째 현실적 검증 연구: 다기관 또는 SNUBH 후향 FFPE 코호트에서 임상 진단용 TG, PAX8, NKX2-1 IHC를 측정하고, 투여 RAI activity, post-therapy scan, stimulated Tg, structural response 및 반복 RAI 여부와 연결. 단순 power simulation이 아니라 실제 endpoint와 cut-off를 사전 정의. 가능하다면 독립 RNA 기반 코호트에서 8-gene score 병행 검증.

### Claude polished draft A — Restrained NC style

This study has several limitations. First, all clinical analyses were retrospective, and no independent cohort contained pretreatment tumour expression together with quantitative radioiodine uptake or prospectively defined treatment response. Accordingly, the observed concordance with post-radioiodine refractory tumours supports biological relevance but does not establish that DM1 predicts response to radioiodine in individual patients. Second, overall-survival analyses were constrained by the low number of events in TCGA-THCA, whereas the MSK cohort was enriched for advanced disease and may overestimate risk relative to an unselected population. The BRAF-specific progression-free interval interaction was identified in a single cohort and requires external replication before it can be interpreted as a treatment-selection biomarker. Third, promoter methylation was associated with, but not shown to cause, repression of lineage genes; analyses of the fusion-negative DM1 subset were also underpowered because of its small sample size. Fourth, cross-platform application required within-sample calibration in the K2 dataset, and the spatial transcriptomic analysis was non-confirmatory. Finally, no internal RNA-sequencing or pathology validation was completed. The three-marker model is an in-silico expression proxy rather than a validated immunohistochemical assay, and TSO500 cannot quantify the RNA expression state examined here. A next step will be a prespecified FFPE study linking routinely available thyroid-lineage markers to administered radioiodine activity, post-treatment uptake and structural response, followed by prospective validation of the full eight-gene score.

### Claude polished draft B — Direct authorial style

The principal limitation is that we did not have a pretreatment cohort in which the DM1 state could be tested directly against measured radioiodine uptake and clinical response. We therefore show that DM1 represents diminished iodine-handling biology and resembles tumours sampled after development of radioiodine refractoriness, but we do not show that it prospectively identifies non-responders. The clinical interaction in BRAF V600E-positive disease is likewise based on a single retrospective cohort and should be treated as hypothesis-generating. Survival estimates are further limited by only 16 overall-survival events in TCGA and by enrichment for advanced disease in the MSK cohort. Mechanistically, methylation remains correlative; the fusion-negative DM1 comparison was underpowered; cross-platform calibration was imperfect in K2; and Visium analysis did not provide confirmatory spatial evidence. We also did not complete internal RNA or tissue validation. In particular, the TG–PAX8–NKX2-1 model is an in-silico proxy and not a clinically validated IHC assay, while TSO500 is a DNA-focused platform that cannot recover the required expression state. The most direct next study is therefore a prespecified retrospective FFPE validation linked to administered radioiodine dose, post-therapy uptake and response, followed by independent prospective testing.

### Claude-added content note

- "전향적 SNUBH IHC validation pending"은 실험 중단 확정과 충돌 → "future validation required"로 교체.
- 13개 한계를 임상·통계·기전·플랫폼·검증 5개 군으로 압축.
- "TSO500 p=0.30"은 본문보다 Supplementary에 두는 편이 자연스러움.

| 항목 | 내용 |
|---|---|
| author-approved? | **YES — 2026-07-30** |
| preferred draft | **B** |
| approval_date | 2026-07-30 |
| inserted_to | `NATURE_COMMUNICATIONS_FULL_DRAFT_v2_2026_07_08.md` §Limitations |

---

## V5 — Cover Letter Opening Paragraph

**파일**: `08_cover_letter.md`

### 승호님 raw input

이 논문이 지금 존재해야 하는 이유는, 갑상선암에서 RAI 치료 여부와 반복 여부가 여전히 tumour의 실제 iodine-handling state를 직접 반영하지 못한 채 결정되는 경우가 많기 때문이다. 우리는 driver mutation과 별개로 유지되는 lineage state를 정의했고, 이 state가 여러 플랫폼과 코호트에서 재현되며 RAI-refractory biology와 정렬된다는 근거를 제시한다.

Nature Communications에 보내는 이유: 단순한 thyroid biomarker 논문이 아니라 cell-lineage state, epigenetic stabilization, multi-platform validation, clinically testable treatment-avoidance framework를 연결했기 때문이다. 이 접근은 lineage dependency가 치료 반응을 규정하는 다른 암종에도 개념적으로 적용될 수 있다. 다만 pan-cancer 적용을 직접 증명한 것은 아니므로 broad conceptual relevance로만 표현한다.

### Claude polished draft A — Restrained NC style

We submit this work because a central problem in differentiated thyroid cancer remains unresolved: radioiodine is administered on the basis of clinicopathological risk, yet the tumour's retained capacity for thyroid-lineage and iodine-handling function is rarely measured directly before repeated treatment. Our study identifies a driver-orthogonal lineage state that is reproducible across transcriptomic, epigenomic, proteomic and single-cell datasets and that aligns with radioiodine-refractory tumour biology. By linking oncogenic context, lineage loss, epigenetic correlates and a feasible path to tissue-based testing, the work provides a broadly relevant example of how cell state can complement genotype in treatment-oriented cancer stratification.

### Claude polished draft B — Direct authorial style

Patients should not have to undergo repeated radioiodine before clinicians learn that their tumours have already lost the molecular programme required to handle iodine. In this study, we define a thyroid-lineage state that is largely independent of canonical driver status, reproduce it across multiple cohorts and molecular platforms, and connect it to radioiodine-refractory biology and outcome stratification in BRAF V600E-positive disease. We believe this work is suited to Nature Communications because it addresses a general problem in precision oncology: genotype alone may not capture the cell state that ultimately determines whether a lineage-dependent therapy remains biologically plausible.

### Claude-added content note

- Cover letter ¶1은 abstract 수치 반복 없이 "why now / why NC" 집중.
- Draft B ("Patients should not…")는 강한 목소리 — 저자 성향에 따라 A가 더 안전.
- "Treatment-oriented stratification"은 즉시 임상 적용이 아닌 표현.

| 항목 | 내용 |
|---|---|
| author-approved? | **YES — 2026-07-30** |
| preferred draft | **A** |
| approval_date | 2026-07-30 |
| inserted_to | `08_cover_letter.md` ¶1 |

---

## V6 — Reviewer Q9 (Mechanism of Differentiation Gene Silencing)

**파일**: `09_reviewer_qa.md` line 46
**현재 상태**: 기존 Q9 답변 전문이 이미 작성되어 있음. 방향 수정 필요.

### 승호님 raw input

기존 Q9의 핵심 입장은 현재 연구 해석과 대체로 일치한다. methylation을 기전적으로 증명된 원인으로 보지 않고 epigenetic correlate로 보는 것이 맞다. DNMT3A/3B, SETDB1 또는 MAPK-effector-dependent chromatin remodelling은 가능한 후보지만, 우리 데이터에서 직접 검증한 upstream regulator가 아니므로 본문에서 강조하면 과장이 된다.

HMA 기반 RAI re-induction 역시 현재 결과가 임상 효능을 뒷받침하는 것이 아니라, methylation과 lineage loss의 연관성이 실험적으로 검증 가능한 가설을 만든다는 의미다.

"Fusion-independent methylation"은 더 조심스럽게 써야 한다. fusion-negative DM1이 존재한다는 사실은 state가 fusion 하나로 완전히 설명되지 않음을 지지하지만, n=19 비교에서 유의하지 않았으므로 independence를 입증하지 못한다. → "observed outside fusion-positive disease, although the subgroup comparison was underpowered" 정도가 적절.

### Claude polished draft A — Restrained reviewer response

We agree that our data do not establish the upstream mechanism responsible for differentiation-gene silencing. We have therefore revised the manuscript to describe promoter methylation as an epigenetic correlate rather than a causal mediator. Although MAPK-dependent recruitment of DNA-methylation or repressive chromatin machinery is biologically plausible, we did not perform perturbation, chromatin-immunoprecipitation or accessibility experiments that would identify a specific regulator. We now state explicitly that the observed pattern motivates mechanistic testing of DNMT- and chromatin-dependent lineage repression rather than confirming such a pathway.

We have also tempered our description of fusion independence. DM1 and elevated methylation were observed outside fusion-positive tumours, indicating that the state is not restricted to kinase-fusion disease. However, the fusion-negative DM1 subgroup was small (n=19), and the direct methylation comparison was underpowered and non-significant. We therefore no longer present this analysis as proof of fusion-independent methylation. Instead, we describe it as evidence that the lineage state can occur in the absence of a detected fusion, while emphasizing the need for larger cohorts.

Finally, we clarify that hypomethylating-agent-based redifferentiation is a testable implication of the observed association, not a therapeutic conclusion supported by the present study.

### Claude polished draft B — Direct authorial reviewer response

We agree with the reviewer that methylation should not be interpreted as the demonstrated cause of lineage loss. Our data establish co-occurrence between promoter hypermethylation and reduced expression of several thyroid-lineage genes, but they do not resolve directionality or identify an upstream chromatin regulator. We have therefore removed mechanistic wording that implied a causal MAPK–methylation–silencing sequence and now describe this model as one of several experimentally testable explanations.

We have also revised the term "fusion-independent methylation." The presence of DM1 among fusion-negative tumours shows that DM1 is not synonymous with fusion status, but the methylation comparison in this subgroup was limited by only 19 fusion-negative DM1 samples and was not statistically significant. The revised text acknowledges this limitation directly and avoids claiming independence. Any reference to hypomethylating-agent-mediated re-induction is now framed solely as a motivation for future perturbation studies.

### Claude-added content note

- 기존 Q9 방향은 대체로 옳으나 "fusion-independent methylation" 표현은 반드시 약화.
- 특정 enzyme 후보(DNMT3A/3B, SETDB1)는 직접 분석 없으므로 reviewer가 요구할 때만 제한적 언급.
- 최종 삽입 시 manuscript revision 위치와 변경 문구 line reference 추가 필요.

| 항목 | 내용 |
|---|---|
| 현재 Q9 전체 승인 여부 | **YES — wording revised per Draft A** |
| author-approved? | **YES — 2026-07-30** |
| preferred draft | **A** |
| approval_date | 2026-07-30 |
| inserted_to | `09_reviewer_qa.md` §Q9 |

---

# 저자 승인 체크

승호님이 각 슬롯을 확인한 뒤 아래 표를 직접 갱신.
선택: **A** / **B** / **수정안** (원하는 수정 내용 간략히)

| Slot | 선택 | 수정 필요 사항 | author-approved? | approval_date |
|---|---|---|---|---|
| V1 Introduction hook | **A** | — | **YES** | 2026-07-30 |
| V2 Introduction aim | **A** | — | **YES** | 2026-07-30 |
| V3 Mechanism discussion | **B** | — | **YES** | 2026-07-30 |
| V4 Limitations | **B** | — | **YES** | 2026-07-30 |
| V5 Cover letter opening | **A** | — | **YES** | 2026-07-30 |
| V6 Reviewer Q9 | **A** | — | **YES** | 2026-07-30 |

---

# 권장 선택 (Claude 제안)

| Slot | 권장 | 이유 |
|---|---|---|
| V1 | **A** (Restrained) | 직접 pretreatment RAI 데이터 없으므로 톤 억제 |
| V2 | **A** (Restrained) | 1차 기여 = state discovery 명확히 |
| V3 | **B** (Direct) | Discussion에서 저자 판단 분명히 드러내는 편이 강함 |
| V4 | **B** (Direct) | Limitations는 솔직할수록 reviewer 신뢰 높아짐 |
| V5 | **A** (Restrained) | Cover letter는 과장 없이 NC fit 설명 |
| V6 | **A** (Restrained) | Reviewer response는 방어적보다 동의+수정 구조가 안전 |
