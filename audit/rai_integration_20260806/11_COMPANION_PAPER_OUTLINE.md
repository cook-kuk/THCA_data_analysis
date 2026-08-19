---
title: "11 — Companion paper outline: radioiodine-response public-data reappraisal"
date: 2026-08-06
author: Claude (managing-editor 역할), 검토자: Seungho Cook
status: 초안 아웃라인 — 저자 승인 및 author-voice 섹션 별도 작성 필요
depends_on: |
  07_MANUSCRIPT_CHANGE_PLAN.md (Strategy B 채택 전제).
  독립적으로 생성된 2차 감사 00_EXECUTIVE_DECISION.md 항목 10("Does a separate radioiodine
  paper stand up? Yes, on today's data alone.")이 이 outline의 전제를 재확인함. 항목 11이
  §13 scope (ii)의 Zhang iProX/GSA 접근을 "가장 가치 있는 다음 데이터"로 독립적으로 동일하게
  지목함.
scope_note: |
  전략 문서 07의 권고에 따라 이 논문은 DM1 8-gene 패널을 검증하는 논문이 아니다.
  "driver mutation class와 발현/분화 축 중 어느 쪽이 갑상선암 방사성옥소(RAI) 불응성을
  더 잘 설명하는가, 그리고 오늘 공개 데이터로 이 질문에 얼마나 답할 수 있는가"에 대한
  체계적 재평가(systematic reappraisal) + 자원(resource) 논문이다. 새 통계 분석 없이
  2026-08-06 자 기존 산출물(그림/표/보고서)만 재구성한다.
---

# 11 — Companion paper outline

## 0. 이 논문이 하는 일과 하지 않는 일

**한다:** (1) TCGA-THCA GDC BCR Biotab에 묻혀 있던 RAI course-level 임상 데이터를 커뮤니티가
쓸 수 있는 자원으로 처음 노출한다. (2) driver mutation class(BRAF/RAS/negative)가 RAI
불응성을 예측하는지 3개 독립 코호트(n=294)에서 재검증하고 null임을 보인다. (3) 8-gene
분화축을 RAI 관련 endpoint에 대해 직접 검증한 오늘까지의 4개 시도(TCGA 1차 반응, GSE151179
uptake, GSE138042, MERAIODE 재분화 시험)를 한 자리에 모아 왜 전부 null 또는 검증력 부족인지
보인다. (4) Zhang 2026의 단백체 기반 분자 서브타입이 같은 코호트·같은 endpoint에서 driver
class보다 훨씬 강하게 불응성과 연관됨을 재현하고, 이것이 "발현/분화 축이 옳은 축"이라는
가설과 어떻게 연결되는지 논한다. (5) 16개 후보 코호트를 스크리닝한 로그를 투명하게 공개해
다음 연구자가 같은 탐색을 반복하지 않도록 한다.

**하지 않는다:** DM1 8-gene 패널 자체를 이 논문의 주인공으로 쓰지 않는다(그 패널의 OS/PFI/
methylation/fusion/single-cell 검증은 전부 Paper 1 소유). TCGA 구조적 재발 prognostic
결과(Firth OR)는 Paper 1의 1차 결과이므로 이 논문에서는 1회 인용만 하고 재도출하지 않는다.

---

## 1. 제목 후보 (영문, 5개 이상)

| # | 제목 | 프레이밍 |
|---|---|---|
| T1 | *Driver mutation class does not predict radioiodine refractoriness in differentiated thyroid cancer: a systematic reappraisal of public data* | 부정적 결론을 정면에 내세움 — 가장 직설적 |
| T2 | *The state of public evidence for molecular prediction of radioiodine response in thyroid cancer: a systematic data audit and resource* | 감사(audit) + 자원(resource) 이중 정체성을 제목에서부터 명시 |
| T3 | *An unexploited TCGA-THCA radioiodine treatment-response resource, and what public data can and cannot yet answer about radioiodine refractoriness* | 자원 우선 프레이밍 — Scientific Data류 저널에 적합 |
| T4 | *Expression state, not driver mutation, tracks radioiodine refractoriness across three independent thyroid cancer cohorts* | 양성 결론(Zhang) 중심 프레이밍 — 가장 낙관적이나 패널 미검증이라는 한계와 다소 긴장 |
| T5 | *Why radioiodine-response prediction has stalled: a data-landscape audit of thyroid cancer public repositories* | 필드 전체를 향한 진단형 제목 — 리뷰어 흥미 유도, 다소 과감 |
| T6 | *A public-data resource and reappraisal of molecular predictors of radioiodine response in differentiated thyroid cancer* | T2의 보수적 변형 — 가장 무난 |

**권고: T2 1차, T6 대안.** 이유: T1과 T4는 각각 부정/긍정 한쪽으로 치우쳐 논문의 실제
내용(양쪽 다 있음: driver 축 null + 발현 축 양성 + 직접검증 4개 null)을 왜곡한다. T2는
"감사"와 "자원"이라는 이 논문의 두 실제 기여를 제목에서부터 정확히 반영하며, T5의 도발적
어조 없이도 필드-차원 진단이라는 임팩트를 유지한다. T3은 자원 요소만 강조된 저널(Scientific
Data 등)에 낼 때의 대안 제목으로 보관한다.

---

## 2. Abstract (영문, ~250 단어)

> Radioiodine (RAI) response prediction in differentiated thyroid cancer has relied heavily on
> driver-mutation class (BRAF V600E versus RAS versus driver-negative), yet no large public
> cohort has directly tested this assumption against adjudicated RAI outcome. We performed a
> systematic search of public repositories for datasets pairing molecular data with a
> radioiodine-specific outcome, screening 16 candidate sources and retaining seven with a usable
> variable (candidate ledger provided as a resource). We first show that driver mutation class
> does not consistently discriminate RAI-refractory disease across three independent cohorts
> (Siraj 2022, n=158; Zhang 2026, n=113; Boucai 2023, n=23; pooled random-effects odds ratio
> 0.75, 95% CI 0.30-1.87, P=0.53, I^2=63%). We then report a previously unexploited resource: the
> TCGA-THCA GDC BCR Biotab clinical-radiation file contains course-level RAI treatment records
> (293 courses, 274 patients) that are absent from every curated cBioPortal view of the same
> study. Using this resource, a thyroid-differentiation transcriptional score was unassociated
> with best response to the first RAI course (n=167, Cohen's d=-0.03, 95% CI -0.50 to 0.46,
> P=0.73). Three further independent attempts to link differentiation-axis expression to a
> direct RAI-effect endpoint (RAI uptake at the metastatic site, GSE151179; a matched
> radioresistant-versus-radiosensitive contrast, GSE138042; response to BRAF/MEK-inhibitor RAI
> redifferentiation therapy, E-MTAB-12837/12900) were each directional but non-significant and
> underpowered. In contrast, within the same advanced-disease cohort and endpoint, a
> proteome-wide molecular subtype identified RAI-refractory disease far more strongly than driver
> class (Zhang 2026, chi-square P=1.4x10^-7, CC1-versus-CC3 odds ratio 24.6), indicating that an
> expression-based axis, not driver identity, is the more promising target. No public dataset
> yet combines panel-compatible expression, tumour purity and adjudicated RAI response at
> sufficient scale (an estimated 226 patients) to test a compact differentiation panel directly.
> We provide the TCGA course-level resource, the pooled driver-class analysis, and the full
> screening ledger to accelerate the next attempt.

(약 255 단어. 저자 검토 후 정확한 단어 수 확인 필요.)

---

## 3. Introduction — 4문단 아웃라인

이 섹션은 실제 영문 산문이 아니라 각 문단이 담아야 할 내용과 논증 순서를 한국어로 정리한
설계도다. 최종 영문 집필은 저자 검토 후 진행한다(Abstract·제목과 달리 이 섹션은 voice-protected
지정 대상은 아니지만, Paper 1의 Introduction hook과 겹치지 않도록 별도로 새로 써야 한다).

**문단 1 — 임상 배경과 미해결 질문.** 분화갑상선암에서 RAI는 여전히 표준 보조치료이지만
불응성 예측은 오랫동안 driver mutation class(특히 BRAF V600E vs RAS)에 의존해 왔다. 이
가정 자체가 대규모 공개 코호트에서 직접 검증된 적이 없다는 공백을 짚는다. Paper 1(DM1
원고)과 겹치지 않도록, 여기서는 "8-gene 패널"이 아니라 "driver class 가정 자체"를 질문
대상으로 설정한다.

**문단 2 — 이 논문이 취한 접근.** 체계적 검색(16개 후보 → 7개 채택, `16_SCREENED_EXCLUSIONS.tsv`)을
수행했고, 그 과정에서 TCGA-THCA GDC biotab에 지금까지 어떤 cBioPortal 큐레이션에도 노출되지
않은 course-level RAI 임상 데이터가 존재함을 발견했다는 것을 이 논문의 방법론적 기여로 제시한다.

**문단 3 — 무엇을 발견했는지 미리 요약(질문 3개).** (a) driver class가 3개 코호트에서 불응성을
일관되게 구분하는가 — 아니다. (b) 발현/분화 축을 RAI 관련 endpoint에 직접 연결한 시도들이
지금까지 무엇을 보여주었는가 — 방향은 일관되나 전부 null 또는 검증력 부족. (c) 그렇다면 어떤
축이 유망한가 — driver가 아니라 발현/단백체 기반 서브타입(Zhang)이다, 그러나 이것이 압축된
패널로 재현되는지는 아직 미검증.

**문단 4 — 이 논문의 기여와 경계를 명시.** 이 논문은 특정 8-gene 패널을 검증하는 논문이
아니라 필드 전체에 대한 데이터 감사이며, 부수적으로 TCGA course-level RAI 자원, driver-축
메타분석, 다음 연구가 필요로 하는 표본크기 목표(≈226명)를 제공한다는 점을 명확히 한다. Paper
1과의 관계(자매 논문, 결과 비중복)를 이 문단 마지막에 1문장으로 명시한다.

---

## 4. Results — 섹션 지도 (6개 주요 결과 섹션 = Fig 1–6과 1:1 대응)

| Results 소제목 | 대응 Figure | 핵심 수치 |
|---|---|---|
| R1. 공개 데이터에서 RAI-특이 표현형을 찾는 일은 예상보다 어렵다 | Fig 1 | 16개 후보 스크리닝, 7개 채택, EXC-01~11 배제 사유 |
| R2. TCGA-THCA는 노출되지 않은 course-level RAI 임상 자원을 갖고 있다 | Fig 2 | 293 courses / 274 patients / 235명 패널-매칭 가능; 최우량 반응 CR 167·PR 21·SD 6·PD 7 |
| R3. Driver mutation class는 3개 독립 코호트에서 RAI 불응성을 일관되게 구분하지 못한다 | Fig 3 | pooled OR 0.75 [0.30,1.87], P=0.53, I²=63%, n=294 |
| R4. 분화/패널 축을 RAI 관련 endpoint에 직접 연결한 4개 시도는 모두 null 또는 검증력 부족이다 | Fig 4 | TCGA 1차반응 d=-0.03; GSE151179 uptake d=+0.37(purity로 소거); GSE138042 d=-0.42; MERAIODE d=+0.34(5/8 유전자만 측정) |
| R5. 그러나 같은 축이 RAI-치료군 내부의 예후 신호와는 연관된다(참조: 자매 논문) | Fig 5 | TCGA Firth OR 0.20 신규 종양사건, 0.22 지속성 질환 — **Paper 1 1차 소유, 여기서는 1회 인용만** |
| R6. Driver가 아니라 발현/단백체 기반 서브타입이 실제로 RAI 불응성을 구분한다 | Fig 6 | Zhang 2026 CC1→CC2→CC3 21%→57%→87%, χ² P=1.4×10⁻⁷, OR=24.6 |

---

## 5. Main Figures — 5–6개, 기존 산출물 매핑

모두 `rai-response-genomics-atlas/results/figures/`에 **이미 2026-08-06 자로 생성되어 있다**.
신규 그림 제작이 필요 없고, 패널 재배치·캡션 작성만 남았다.

| Fig | 제목(안) | 사용할 기존 파일 | 비고 |
|---|---|---|---|
| **Fig 1** | 데이터 랜드스케이프: 발견된 모든 후보 코호트, 크기 vs RAI-효과 직접성, 접근성 | `figure_synthesis_landscape_2026_08_06.{png,pdf}` | 이미 "크기 × 직접성 × access status" 3축으로 설계되어 있음(스크립트 `synthesis_figures_2026_08_06.py` Figure C) — 그대로 Fig 1로 승격 가능 |
| **Fig 2** | TCGA-THCA course-level RAI 임상 자원: 코호트 구성 + 용량/반응 분포 | `figure_tcga_rai_best_response_2026_08_06.{png,pdf}` (핵심 패널) + `03_TCGA_COHORT_CONSTRUCTION_AUDIT.md`의 코호트 원장(293/274/249/237/235)을 표 패널로 추가 | 자원 소개용이므로 CONSORT류 흐름도 패널 1개 추가 권장(신규 제작 필요, 유일한 미존재 패널) |
| **Fig 3** | Driver class는 RAI 불응성을 예측하지 않는다 (3-코호트 pooled) | `figure_pooled_driver_vs_rai_2026_08_06.{png,pdf}` | 이미 forest plot으로 존재(`pooled_driver_vs_rai_2026_08_06.py`) |
| **Fig 4** | 분화축 vs RAI endpoint: 4번의 직접 검증, 모두 null/underpowered | `figure_synthesis_effects_2026_08_06.{png,pdf}` | "모든 RAI-endpoint 검증을 하나의 effect-size 축에" 정확히 이 그림의 설계 목적(스크립트 Figure A) |
| **Fig 5** | (자매 논문 인용 패널) RAI-치료군 내 구조적 재발과의 연관 — 예후적, RAI-특이적 아님 | `figure_tcga_rai_structural_2026_08_06.{png,pdf}` | **주의:** Paper 1이 이 결과의 1차 소유자. 이 논문에서는 "reproduced from [Paper 1, in preparation]"로 명시하고 자체 신규 해석을 추가하지 않는다 (중복출판 방지, §07 참조) |
| **Fig 6** | 발현/단백체 서브타입이 불응성을 구분한다 (Zhang 2026) | `figure_zhang2026_driver_rai_2026_08_06.{png,pdf}` | 이미 존재 |

Fig 5를 제외하면 5개 그림 모두 오늘 이미 생성된 파일 그대로 승격 가능하다. Fig 5를 포함할지
여부(자매 논문과의 인용 경계가 편집적으로 부담스러우면 제외 가능)에 따라 **5개 또는 6개**
main figure 구성 중 선택할 수 있다 — 과제 요구사항(5–6개)과 정확히 일치.

---

## 6. Supplementary Figures

| Supp Fig | 내용 | 소스 |
|---|---|---|
| S1 | 검증력 문제: 각 코호트가 검출 가능했던 최소 효과크기 vs 관찰된 효과크기가 요구하는 표본크기(≈226명 목표) | `figure_synthesis_power_2026_08_06.{png,pdf}` (이미 존재) |
| S2 | TCGA index-course 민감도 분석(첫 course/마지막 course/최고용량/최선반응/최악반응 5가지 규칙, 전부 동일 null) | `03_TCGA_COHORT_CONSTRUCTION_AUDIT.md` 표 → 그림화 |
| S3 | TCGA rare-event 재분석: Firth vs 일반 ML 로지스틱 비교, LOO, bootstrap 안정성 | `04_STATISTICAL_REANALYSIS.md` 표 → 그림화 |
| S4 | GSE151179 상세: 음성대조군(비종양 갑상선) 포함 전체 6개 대비, purity 보정 모델 계수 | `figure_gse151179_uptake_at_met_site_2026_08_06.{png,pdf}` |
| S5 | GSE138042 상세: per-gene 표(TG/TSHR/SLC5A5 등), **estimand 정정**(최초 보고된 "composition-adjusted OR=0.288, P=0.020"은 13 vs 10의 올바른 head-to-head가 아니라 95개 전 라이브러리 중 "RAIR-* 명명 여부"라는 다른 질문에 적합된 모델이었음 — 올바른 대비에 동일 공변량 적용 시 OR=0.44, P=0.31로 비유의), **배치 교란**(refractory 7/13이 다른 군과 다른 조직학 라벨·라이브러리 명명을 사용 — 군 소속이 주석 배치와 완전히 얽힘) | `figure_gse138042_rair_panel_2026_08_06.{png,pdf}` + `04c_GSE138042_ESTIMAND.md` |
| S6 | MERAIODE/재분화 시험 assay 한계: HTG EdgeSeq 패널이 8개 유전자 중 5개만 측정 가능함을 보이는 유전자 커버리지 표 | `figure_meraiode_redifferentiation_2026_08_06.{png,pdf}` |
| S7 | Boucai 2023 코호트 흐름 해소(n=23 vs 24 불일치, ERAI_3 배제, A/B pair 독립성 검증) | `05_EXTERNAL_COHORT_AUDIT.md` §1 → 그림/표화 |
| S8 | 전체 스크리닝 원장(16개 후보, 포함/배제 사유 전체) | `16_SCREENED_EXCLUSIONS.tsv` 그대로 Supplementary Table로 |
| S9 | Siraj 2022 단독 상세(refractory vs avid 내부 타당성 점검: 누적 용량, Tg, PFS) + PFS 순환성(circularity) 경고 | `figure_siraj2022_driver_rai_2026_08_06.{png,pdf}` |

---

## 7. Methods 아웃라인

1. **체계적 검색 및 스크리닝 프로토콜.** 검색 전략(키워드 + phenotype-기반 GEO 재탐색,
   `GSE173248` 사례처럼 "radioiodine"이라는 단어 없이도 관련 phenotype이 존재하는 경우 포함),
   포함/배제 기준, `16_SCREENED_EXCLUSIONS.tsv` 형식(candidate_id, decision,
   exclusion_reason 등) 그대로 방법 절 표준으로 기술.
2. **TCGA-THCA GDC BCR Biotab 자원화.** 파일 경로, `radiation_adjuvant_units`(mCi/Gy/cGy
   구분), `treatment_best_response` 필드가 course-specific이지 patient-level이 아니라는
   발견(코호트 구성 감사 §1), index-course 선택 규칙과 5가지 민감도 분석.
3. **외부 코호트 획득 절차의 투명한 기록.** GSE138042 임상주석이 PMC의 JS 게이트를 우회하는
   EuropePMC supplementary-files 엔드포인트로만 얻어졌다는 것, Zhang 2026 표가 gated
   accession(iProX/GSA) 대신 오픈 access supplemental PDF에서 `pdftotext -layout` +
   fixed-width 파싱으로 추출되었다는 것, Siraj 2022가 CC-BY 오픈 표라는 것 — 각각 재현
   가능하도록 정확히 기술.
4. **통계 방법.** Cohen's d + bootstrap CI(5,000 draws), Mann-Whitney, Firth-penalized
   로지스틱(희귀 event 보정), random-effects 메타분석(driver-class pooling, DerSimonian-Laird
   또는 사용된 방법 명시), leave-one-out/bootstrap 안정성 진단, 검정력 계산 공식.
5. **패널 점수 정의.** Paper 1과 동일한 8-gene z-score 평균 정의를 재사용한다는 점을 명시하고
   Paper 1을 인용(패널 정의 자체를 여기서 재도출하지 않음).
6. **재현성.** seed(20260806), bootstrap 반복수, 소프트웨어 버전(pandas/scipy/statsmodels/
   lifelines), 각 분석의 원본 스크립트 파일명을 표로 제공.

---

## 8. Discussion 포인트

1. **주요 발견.** Driver class는 재현성 있게 RAI 불응성을 구분하지 못한다(3개 코호트, n=294,
   I²=63%로 이질적이지만 방향은 null 쪽). 발현/단백체 축은 같은 코호트·같은 endpoint에서
   훨씬 강하게 구분한다(Zhang). 그러나 압축된 8-gene 전사체 패널로 이를 직접 재현하려는 4번의
   시도는 모두 null 또는 검증력 부족이었다.
2. **선행 연구와의 관계.** ESTIMABL2(NEJM 2022)와 IoN(Lancet 2025)이 저위험군에서 RAI
   생략이 non-inferior임을 보인 맥락에서, TCGA-THCA처럼 86% 완전관해를 보이는 인구 기반
   코호트는 애초에 RAI-response 판별력이 거의 없는 집단이라는 점(이것이 TCGA null의 가장
   경제적인 설명). Boucai 2023/Mu 2024/Siraj 2022/Laschinsky 2023 등 양성 결과들은 모두
   전이성/진행성 질환에서 나왔다는 패턴을 명시.
3. **그럴듯한 해석.** 발현/분화 상태 축이 RAI 처리 생물학과 관련된 진짜 축일 가능성이 높지만,
   현재 공개 데이터와 현재 형태의 8-gene 패널로는 이를 증명할 수 없다 — 데이터 공백이지
   가설 반증이 아니다.
4. **경쟁 설명.** 8-gene 패널(effector 5 + TF 3) 자체가 RAI-response 예측에는 최적이 아닌
   feature set일 가능성 — Zhang의 단백체 consensus subtype은 전체 프로테옴 기반이라 훨씬
   넓은 정보를 담고 있어, 압축 패널이 포착하지 못하는 신호를 포함할 수 있음을 정직하게 논한다.
5. **추론의 경계.** 전부 후향적, 코호트 간 불응성 정의 상이(Siraj 7-기준 vs Zhang 영상 기반
   vs Boucai RECIST), Zhang 수치는 미보정(disease burden/전이부위/병기/세포충실도 등 통제
   불가), MERAIODE는 애초에 패널의 3개 유전자를 측정할 수 없는 assay라는 구조적 한계.
6. **구체적 다음 검증 단계.** (a) Zhang iProX/GSA 데이터 접근 신청(교신저자 1명이 데이터
   접근 위원회 담당자와 동일 — 이메일 1통으로 해결 가능하다고 보고서에 명시됨) → 실제 패널을
   n=113 불응성 라벨에 직접 검증. (b) Mu 2024(HRA004166, n=214) DAC 신청. (c) 설계 목표
   ≈226명(GSE151179 power 계산에서 도출)을 충족하는 전향적 코호트 확보.

---

## 9. Limitations

- 후향적, 코호트 간 이질적 불응성 정의(7-기준/영상기반/RECIST) — pooled 메타분석 I²=63%로
  강한 주장 불가.
- TCGA course-level 주석은 GDC biotab 계층에만 존재하며 반응이 평가 가능한 course는
  176/249건뿐 — 완전한 자원이 아니라 부분 자원.
- 구조적 outcome event 수가 매우 작음(9, 13건) — 이 결과 자체는 Paper 1 소유이며 여기서는
  1회 인용에 그치므로 이 논문의 결론(driver 축 null, 발현 축 promising)에는 이 작은 표본이
  결정적으로 기여하지 않는다는 점을 분명히 한다.
- MERAIODE는 구조적으로 8개 유전자 중 3개(NIS·TG·DIO1)를 측정할 수 없는 assay라서, 이 논문의
  "null"이라는 라벨은 "증거 없음"이 아니라 "측정 불가능"에 가깝다는 것을 명시.
- Zhang 수치는 gated 데이터가 아니라 오픈 supplemental PDF에서 추출한 값 — 원 저자 승인
  하의 재현이 아니므로 방법 절에서 QC(파싱된 행 수가 보고된 코호트 크기와 일치하는지 등)를
  투명하게 밝힌다.
- Siraj 2022의 PFS 대비 refractoriness 정의는 순환적(refractoriness 정의 자체에 구조적
  진행이 포함됨) — 이 논문은 이 조합을 아예 사용하지 않는다는 것을 명시(이미 `05_EXTERNAL_
  COHORT_AUDIT.md`가 지정한 금지 문구 계승).
- 8-gene 패널 자체를 이 논문의 어떤 분석에서도 직접 RAI outcome에 대해 검증하지 못했다 —
  이것이 이 논문 전체의 핵심 한계이자, 동시에 존재 이유(그 검증이 왜 아직 불가능한지 설명하는
  논문이기 때문).

---

## 10. Data / Code Availability

- TCGA-THCA GDC BCR Biotab 방사선 임상 파일: 공개, 경로
  `nationwidechildrens.org_clinical_radiation_thca.txt`(GDC 공개 다운로드) — 파생 course-level
  TSV(`tcga_rai_response_main_2026_08_06.tsv` 등)를 Zenodo에 자원으로 별도 기탁.
- GSE138042, GSE151179, E-MTAB-12837/E-MTAB-12900: 모두 공개 accession, 그대로 인용.
- Siraj 2022 Supplementary Tables S1/S3: CC-BY, 파생 `siraj2022_patient_level_2026_08_06.tsv`
  (158행, 재사용 가능)를 Zenodo에 기탁.
- Zhang 2026: 1차 데이터(iProX IPX0011848000, GSA HRA011340)는 **제한 접근**이며 이 논문은
  이를 재배포하지 않는다. 여기서 쓰인 표는 오픈 supplemental PDF(`mmc1.pdf`)에서 저자가 직접
  전사한 것이며, 이 파생 표(`zhang2026_patient_table_2026_08_06.tsv`,
  `zhang2026_mutations_2026_08_06.tsv`)만 기탁한다 — 원 저자 및 원 논문에 대한 명시적 귀속을
  Methods와 Data Availability 양쪽에 기재.
- 16개 후보 스크리닝 전체 원장(`16_SCREENED_EXCLUSIONS.tsv`)을 Supplementary Table로 공개 —
  이 자체가 다음 연구자를 위한 자원.
- 전체 분석 스크립트(`scripts/*2026_08_06*.py`, `scripts/audit/*.py`)를 공개 저장소에 기탁,
  seed 및 소프트웨어 버전 명시.

---

## 11. 목표 저널과 근거

이 논문의 실제 증거 구성(3개 코호트의 이질적 null 메타분석 + 새로 노출된 자원 데이터셋 +
타인 데이터(Zhang)의 재현 + 자매 논문 결과 1회 인용)은 **NC급을 채우지 못한다** — 사용자의
출판 기준(NC 이상 우선, npj/JCI Insight를 fallback으로 허용)에 비추어 정직하게 이 등급으로
분류한다.

| 순위 | 저널 | 근거 |
|---|---|---|
| **1차 권고** | **npj Precision Oncology** | 재평가/음성결과/자원 성격의 논문을 명시적으로 환영하는 포맷 보유, 오픈 데이터 재현성 강조와 논문의 방법론적 기여(스크리닝 로그, TCGA biotab 노출)가 잘 맞음. Cell Press/NPG 계열이라 Paper 1(NC 목표)과 게재 시 상호 인용 및 "companion paper" 언급이 자연스러움 |
| **대안 (임상 중심)** | **Thyroid** 또는 **JCEM** | 임상 갑상선학 독자에게 "driver class 가정을 재검증했다"는 메시지가 정확히 도달하는 저널. 통계적으로 화려하지 않은 재평가 논문도 임상 실천 함의가 있으면 받아들이는 전통이 있음 |
| **자원 강조 시 대안** | **Scientific Data** | TCGA course-level RAI 자원 자체를 data descriptor 형식으로 낸다면 이 저널이 최적. 단, driver-축 메타분석/논증 서사를 크게 축소해야 하므로 이 논문 전체를 옮기기보다는 자원 부분만 별도 data descriptor로 분리하는 4번째 옵션으로 고려 |

**결정적 조언:** NC급을 목표로 이 논문의 스코프를 인위적으로 부풀리지 말 것. 이 논문의 정직한
가치는 "필드가 가정해온 것(driver class)이 틀렸고, 유망한 축(발현/단백체)이 아직 검증 안
됐다는 것을 처음으로 체계적으로 보였다"는 데 있으며, 이는 npj Precision Oncology 급에서
충분히 인용될 가치가 있는 기여다. §13의 scope (ii)/(iii)가 실현되면 저널 목표를 재상향한다.

---

## 12. Reviewer 공격 지도

| 예상 공격 | 준비된 답변 |
|---|---|
| "결국 null 결과 모음 아닌가" | 같은 가설에 대한 4번의 독립적·투명하게 문서화된 시도가 전부 null인 것은 그 자체로 필드에 유용한 정보다. 특히 이전에는 "그런 데이터가 아예 없다"고 가정되었던 상황(우리 자신의 Paper 1 초안도 그렇게 썼다)에서, 데이터가 실제로 존재하고 직접 검증했더니 null이라는 것은 침묵보다 강한 정보다. 더불어 새로 노출한 TCGA 자원과 driver-축 null 재현(n=294)이라는 양성 기여가 있다. |
| "TCGA는 애초에 잘못된 인구집단이라 null이 아무것도 말해주지 않는다" | 동의하며 본문에서 명시한다(86% 완전관해, ESTIMABL2/IoN 적격 인구와 일치). TCGA를 RAI-response 질문에 대한 "적절히 검정력 있는 음성대조군"으로 프레이밍하고, 진짜 표적 인구(진행성/전이성 질환)는 Zhang/Mu/Siraj/Boucai임을 명시한다. |
| "Zhang 서브타입 분석에 교란변수 보정이 없다" | `05_EXTERNAL_COHORT_AUDIT.md`가 이미 명시한 한계(질병 부담, 전이 부위, 병기, 세포충실도, 배치 효과 통제 불가)를 그대로 채택해 "unadjusted association"으로만 보고하고, "confounding으로 설명 안 됨" 같은 금지 문구를 쓰지 않는다. |
| "3-코호트 driver 메타분석 I²=63%는 신뢰할 수 없다" | 동의. `05_EXTERNAL_COHORT_AUDIT.md`가 지정한 완화된 문구("did not consistently discriminate", 세 코호트 모두 일관적이라는 말 금지)를 그대로 사용하고, CI [0.30, 1.87]의 폭을 본문에서 숨기지 않는다. |
| "이게 왜 DM1(Paper 1) 원고와 별개인가 — 사실상 같은 프로젝트의 중복 출판 아닌가" | Methods 마지막 문단에 명시적 non-overlap 선언 포함: Paper 1은 8-gene 패널이 TCGA/MSK/5개 외부 코호트에서 driver-orthogonal 상태·생존과 갖는 관계만 다루며 이 논문의 어떤 수치도 재사용하지 않는다. 이 논문은 8-gene 패널을 검증하지 않으며 Paper 1의 결과를 오직 R5/Fig 5에서 1회, 출처를 명시하고 인용한다. 두 원고 모두 상호 각주로 자매 논문임과 결과 비중복을 명시. |
| "TCGA 구조적 재발 결과(Firth OR)를 여기 실었는데 event가 9/13건이라 노이즈 아닌가" | 이 결과는 Paper 1의 1차 소유이며 이 논문에서는 해석을 추가하지 않고 "다른 논문에서 보고된, RAI-특이성이 확립되지 않은 예후 신호"로만 인용한다 — 이 논문 자체의 결론(driver 축 null, 발현 축 promising)은 이 작은 표본에 의존하지 않는다. |
| "Zhang 표를 PDF에서 pdftotext로 뽑았다는데 신뢰할 수 있나" | Methods에 QC 절차를 명시: 파싱된 113개 행이 논문이 보고한 코호트 크기와 정확히 일치, 14개 필드 모두 렌더링됨, 175개 변이 행도 별도로 파싱해 교차 검증. 완전한 검증은 향후 iProX/GSA 접근 신청(§8, discussion 6번) 이후 가능함을 명시. |
| "MERAIODE 섹션은 null이고 패널도 못 재는데 왜 넣었나" | 이 자체가 이 논문의 논지: 공개된 유일한 재분화-시험 전사체가 구조적으로 압축 패널을 검증할 수 없는 assay(HTG EdgeSeq, NIS/TG/DIO1 결측)를 썼다는 사실은, 다음 시험을 설계하는 사람에게 직접적으로 유용한 정보다 — 결과가 아니라 공백(gap) 진술로 명시적으로 위치시킨다. |
| "GSE138042의 composition-adjusted 결과(OR=0.288, P=0.020)를 인용했던데, 이것도 유의하다는 것 아닌가" | 아니다 — 이 값은 다른 질문(13 vs 10 head-to-head가 아니라 95개 전체 라이브러리 중 RAIR-* 명명 여부)에 적합된 모델이며, 본 논문은 이 오염된 추정치를 인용하지 않는다. 올바른 대비에 동일 공변량을 적용하면 OR=0.44, P=0.31로 비유의함을 Methods와 S5에서 명시하고, 나아가 refractory군이 조직학 라벨·라이브러리 명명 배치와 완전히 교란되어 있다는 사실을 추가로 공개해 이 코호트를 "지지도 반박도 못 하는" 코호트로 명확히 재분류한다. |

---

## 13. 스코프 3버전

### (i) 최소 제출가능 — 오늘 공개 데이터만

- §5의 5–6개 main figure 전부 이미 존재하는 파일로 구성, §7 Methods 전부 오늘 자 스크립트로
  재현 가능.
- Zhang은 오픈 PDF 추출 표만 사용(§10에 명시된 한계 그대로 보고).
- 목표 저널: npj Precision Oncology(1차) / Thyroid·JCEM(대안).
- 실행 시간: 신규 분석 없이 그림 조합 + 원고 집필만 필요 — 수 주 내 제출가능 스코프.

### (ii) 확장 — Zhang 113명 프로테오믹스(iProX/GSA) 접근 승인 시

- 추가 Fig 7: **실제 8-gene 패널**(전사체 또는 단백체 매칭 가능한 부분집합)을 Zhang의
  n=113 불응성 라벨에 직접 검증 — 이는 이 논문 전체에서 유일하게 "패널 자체의 RAI-outcome
  직접검증"이 되는 분석이며, §9 Limitations의 핵심 항목("패널을 어떤 분석에서도 직접
  검증하지 못했다")을 해소한다.
- 이 경우 논문의 정체성이 "재평가"에서 "재평가 + 최초의 진짜 패널 검증"으로 상승하며,
  Zhang 원 논문이 이미 Cell Reports Medicine(2026)에 실렸다는 점을 고려하면 목표 저널을
  **Cell Reports Medicine** 또는 동급으로 상향 검토할 수 있다.
- 소요: 교신저자 1인(Xiao Shi)에게 이메일 1통으로 양쪽 접근 신청 가능하다고 보고서에 이미
  명시됨 — 승인 소요 기간이 유일한 병목.

### (iii) 결정적 — 한국 전향적 코호트 확보 시

- SNUBH 등에서 전향적으로 모집한 코호트에 대해 사전등록된 8-gene 패널(RNA 또는 IHC) 점수와
  ATA 기준에 따라 판정된 RAI 반응을 직접 짝지어 측정. §6-S1에서 도출한 설계 목표(≈226명,
  또는 관찰된 효과크기에 맞춘 현실적 부분표본)를 그대로 검정력 계산의 기준으로 사용.
- 이 버전이 실현되면 이 companion paper가 최초로 "treatment-predictive"라는 단어를 정당하게
  쓸 수 있는 유일한 데이터가 되며, 이 시점에는 Paper 1과 재통합해 하나의 상향된 원고로
  묶는 것을 재검토할 가치가 있다(현재는 이르다 — §07 Strategy 결정은 오늘 시점 데이터 기준).
- 목표 저널: 이 시점에는 **Nature Communications 또는 Nature Medicine** 재상향이 정당화됨.

---

## 14. Paper 1과의 상호 참조 문구 (양쪽 원고에 삽입 권장, 저자 확정 필요)

Paper 1 Methods/Discussion 말미 권장 문구(초안):
> "A companion analysis (Cook et al., in preparation) systematically re-evaluates driver
> mutation class and expression-based molecular subtypes as predictors of radioiodine
> refractoriness across independent public cohorts; the present manuscript does not repeat
> that analysis and restricts its own radioiodine-related claims to direct tests of the
> eight-gene panel described here."

Companion paper Methods 말미 권장 문구(초안):
> "The eight-gene thyroid-differentiation panel used for one cross-reference in this study
> (Fig. 5) is defined and independently validated for survival and driver-orthogonality in a
> companion manuscript (Cook et al., in preparation); the present study does not re-derive or
> re-validate that panel and confines its own contribution to driver-class and public-cohort
> radioiodine-response evidence."

이 문구는 저자가 두 원고의 실제 제출 순서와 "in preparation" 상태를 확인한 뒤 확정해야 한다.
