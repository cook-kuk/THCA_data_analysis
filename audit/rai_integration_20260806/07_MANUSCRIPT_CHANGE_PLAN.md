---
title: "07 — RAI 통합 전략: Paper 1(DM1) 반영 방식 결정 + 제목 재평가"
date: 2026-08-06
author: Claude (managing-editor 역할), 검토자: Seungho Cook
status: 결정 권고안 — 저자 최종 승인 필요
inputs:
  - rai-response-genomics-atlas/results/reports/tcga_rai_best_response_brief_2026_08_06.md
  - rai-response-genomics-atlas/results/reports/tcga_rai_structural_disease_brief_2026_08_06.md
  - rai-response-genomics-atlas/results/reports/gse138042_rair_panel_brief_2026_08_06.md
  - rai-response-genomics-atlas/results/reports/gse151179_uptake_at_met_site_brief_2026_08_06.md
  - rai-response-genomics-atlas/results/reports/meraiode_redifferentiation_brief_2026_08_06.md
  - rai-response-genomics-atlas/results/reports/siraj2022_driver_vs_rai_brief_2026_08_06.md
  - rai-response-genomics-atlas/results/reports/zhang2026_driver_vs_rai_brief_2026_08_06.md
  - audit/rai_integration_20260806/03_TCGA_COHORT_CONSTRUCTION_AUDIT.md
  - audit/rai_integration_20260806/04_STATISTICAL_REANALYSIS.md
  - audit/rai_integration_20260806/05_EXTERNAL_COHORT_AUDIT.md
  - audit/rai_integration_20260806/16_SCREENED_EXCLUSIONS.tsv
  - project/manuscript_v8/NATURE_COMMUNICATIONS_FULL_DRAFT_v2_2026_07_08.md (v2, 2026-07-08/07-30)
  - project/manuscript_v8/SUBMISSION_CONTROL/*.md
  - audit/rai_integration_20260806/00_EXECUTIVE_DECISION.md (독립 2차 감사, 이 문서 작성 중 동시 생성됨)
  - audit/rai_integration_20260806/04b_GSE151179_INDEPENDENCE.md
  - audit/rai_integration_20260806/04c_GSE138042_ESTIMAND.md
  - audit/rai_integration_20260806/06_CLAIMS_SUMMARY.md + 06_CLAIMS_REGISTRY.csv (원고 전체 101건 claim 감사)
scope_note: |
  이 문서는 전략 판단 문서다. 새 통계 분석을 수행하지 않았다. 모든 수치는 위 소스 파일에서 직접
  인용했다. 작성 도중 동일 디렉터리에 독립적인 2차 감사 패스(00/04b/04c/06 파일)가 생성된 것을
  확인했다 — 이 문서의 Strategy A/B 권고 및 제목 판단(§3–§4)과 그 2차 감사의 결론(00_EXECUTIVE_
  DECISION.md 항목 8–10)은 서로 독립적으로 도출되었음에도 일치한다. 이 문서는 2차 감사에서 발견된
  추가 세부사항(GSE138042 estimand mismatch, 원고 전체 claim registry)을 반영해 갱신되었다.
---

# 07 — RAI 통합 전략 결정 문서

## 0. 한 줄 결론 (Executive summary)

**Strategy B(분리 논문)를 권고한다.** DM1 원고(Paper 1)에는 오늘(2026-08-06) 새로 발굴한 TCGA
RAI course-level 임상 데이터를 이용해 **패널 자체를 직접 검증한 결과만** — 즉 (1) 1차 RAI 반응은
null이고 (2) RAI-치료군 내에서 구조적 재발과는 연관이 있으나 이는 예후적(prognostic)일 뿐 치료
반응 예측적(predictive)이 아니라는 내용만 — 정직하게 편입한다. 그 외의 훨씬 방대한 오늘 자 작업
(TCGA biotab 자원 자체, 3개 코호트 driver-축 메타분석 n=294, Zhang 2026 단백체 서브타입,
GSE138042/GSE151179/MERAIODE의 null들, 16개 후보 코호트 스크리닝 로그)은 **DM1 패널을 직접
검증하지 않는 별도의 질문**("무엇이 RAI 불응성을 예측하는가에 대한 공개 데이터 현황 감사")이므로
독립된 companion paper로 분리한다.

이 판단의 핵심 근거는 통계가 아니라 **논문의 정체성**이다. Paper 1은 "8-gene 패널이 driver 축과
독립적인 갑상선 계통 상태를 정의한다"는 논문이지 "RAI 반응을 예측하는 진단검사를 개발했다"는
논문이 아니다(SUBMISSION_CONTROL/DM1_PROJECT_STATE.md가 이미 이 경계를 명시). 오늘 자 데이터는
이 경계를 두 가지 방향에서 재확인한다: RAI 반응 예측이라는 주장은 직접 검증에서 지지되지 않고,
반면 companion paper가 다루는 "driver 축이 아니라 발현/분화 축이 맞는 축이다"라는 더 큰 주장은
오히려 강해졌다(Siraj+Zhang+Boucai 드라이버 null, Zhang 단백체 서브타입 강한 양성). 이 두 주장을
한 논문에 넣으면 서로 다른 증거 기준(패널-특이적 vs 축-일반적)이 뒤섞여 리뷰어에게 "결과를
선택적으로 조합했다"는 공격 지점을 만든다.

---

## 1. 오늘 자 증거의 정직한 요약

| # | 분석 | 코호트 | 결과 | 해석 |
|---|---|---|---|---|
| 1 | 1차 RAI course 반응(non-CR vs CR) | TCGA-THCA, n=167 (5가지 index-course 규칙 전부 동일 결론) | first-course(가장 올바른 index): d=0.00, 95% CI [−0.45, +0.46], P=0.91; best/last/highest-dose 규칙: d=−0.03, P=0.73 | **깨끗한 null, 규칙 선택과 무관.** 2차 감사(00_EXECUTIVE_DECISION.md §1)가 "detects d≥0.63 at 80% power이므로 임상적 효과를 배제한다"는 최초 표현이 **동등성 검정과 최소검출효과 계산을 혼동**한 것이라며 철회함 — 정확한 표현은 "신뢰구간이 자신의 경계보다 큰 효과를 배제하며, 연관성이 검출되지 않았다"임. 권고 원고 문장: "The thyroid-lineage differentiation score showed no detectable association with the recorded response to the initial radioiodine course in this predominantly complete-response cohort (Cohen's d = 0.00, 95% CI −0.45 to +0.46, P = 0.91; 141 of 167 evaluable patients achieved complete response)." |
| 2 | RAI-치료군 내 구조적 재발(new tumor event) | TCGA-THCA, n=145, event=9 | Firth OR 0.20–0.24, P=0.017–0.024; LOO 방향 100% 안정, bootstrap 98.8%가 OR<1 | **실재하는 신호.** 그러나 non-RAI군 event가 2건뿐이라 treatment×score interaction 추정 불가 → **prognostic-within-treated이지 RAI-predictive가 아님** |
| 3 | RAI-치료군 내 지속성 질환(persistent disease, 3개월) | TCGA-THCA, n=72, event=13 | Firth OR 0.22–0.24, P=0.019–0.031 | 동일 해석. 단 이 endpoint는 수술 후 3개월 창이 RAI 시점과 겹쳐 순수 post-RAI 지표가 아님 |
| 4 | RAI 흡수(met-site uptake) | GSE151179, n=32 patients | d=+0.37, P=0.53; purity 보정 후 β=+0.035, P=0.88; 음성대조군(비종양 갑상선)이 오히려 d=+0.65로 더 큼 | **null, 그리고 purity 교란으로 완전히 설명됨.** 방향은 가설과 일치하나 해석 불가 |
| 5 | RAIR vs 방사성옥소 민감 (올바른 head-to-head) | GSE138042, n=23 (13 vs 10) | d=−0.42, 95% CI [−1.44, +0.42], P=0.34; power=0.16 | **null, underpowered, 게다가 배치(batch) 교란까지 확인됨.** 기존에 보고된 큰 효과(d=−0.94~−1.03)는 병기 교란(넓은 비교군에 양성/조기 종양 포함)이었음. **2차 감사(04c_GSE138042_ESTIMAND.md)가 추가로 발견:** 최초 보고된 "composition-adjusted OR=0.288, P=0.020"은 실제로는 **다른 질문**(95개 전체 라이브러리 중 "RAIR-* 이름을 가진 라이브러리인가")에 적합된 모델이며, 올바른 head-to-head(13 refractory vs 10 sensitive)에 동일 공변량을 적용하면 OR=0.44 (0.091–2.14), P=0.31로 **유의하지 않음.** 더 심각하게, refractory군 7/13는 다른 군이 전혀 쓰지 않는 조직학 라벨("Papillary **thyroid** cancer")을 쓰고 라이브러리 명명 규칙(`RAIR-*` vs `TC-*`)도 완전히 갈려 — **군 소속이 주석 배치(annotation batch)와 완전히 교란**되어 있음. 이 코호트는 "지지도 반박도 못 하는, 정보가 없는" 코호트로 재분류해야 함 |
| 6 | 재분화 시도(BRAF/MEK+RAI) 반응 | E-MTAB-12837/12900, n=21 (9 vs 12, RECIST) | d=+0.34, P=0.80; 5/8 유전자만 측정 가능(NIS·TG·DIO1 결측) | **null, 검증력 부족 + 자산(assay) 자체가 패널을 못 잼.** 결과가 아니라 공백(gap) 진술로만 사용 가능 |
| 7 | Driver class → RAI 불응성 | Siraj 2022 (n=158) + Zhang 2026 (n=113) + Boucai 2023 (n=23), pooled n=294 | pooled OR 0.75, 95% CI [0.30, 1.87], P=0.53, I²=63% | **driver 축은 불응성을 일관되게 구분하지 못함.** 이질성 커서 강하게 주장 불가 |
| 8 | 발현/단백체 서브타입 → RAI 불응성 | Zhang 2026, n=113, 동일 코호트·동일 endpoint | χ² P=1.4×10⁻⁷; CC1 vs CC3 OR=24.6 | **강한 양성.** 그러나 8-gene 패널이 아니라 전체 단백체 consensus subtype — 패널 검증이 아님 |

**핵심 긴장 관계 — Fig 6D 재조정 필요(양쪽 전략 공통, 반드시 선행):** 현재 원고 Fig 6D /
Abstract는 "GSE151179 post-RAI refractory tumours shifted toward DM1 (d≈−1.0, P≈10⁻⁴)"를
핵심 "생물학적 anchor"로 쓴다(`05_figure_captions_NC.md:125`, `NATURE_COMMUNICATIONS_FULL_DRAFT_v2:124`).
오늘 GSE151179를 **RAI 효과에 가장 직접적인 필드**(`rai uptake at the metastatic site`)로
재검증하니 null·purity-교란이었다(#4). 두 결과는 서로 다른 GEO 필드/비교(및 표본짝짓기)를
쓰고 있어 통계적으로 모순은 아니지만, **하나의 논문 안에 "d≈−1.0 강한 정렬"과 "d=+0.37 null·교란"을
같은 데이터셋에서 동시에 제시하면 리뷰어가 즉시 cherry-pick으로 읽는다.** 이 재조정은 Strategy
A/B 어느 쪽을 택하든 착수 전 반드시 코드 레벨에서 어느 비교가 1차 GSE151179 RAI-anchor
주장인지 확정해야 한다 (evidence-lock 규칙: 출처 상충 시 두 값 보존 + primary output에서 해결,
더 매력적인 값 선택 금지).

---

## 2. Strategy A vs Strategy B — 정면 비교

| 항목 | **Strategy A — DM1 원고에 보수적으로 통합** | **Strategy B — 별도 RAI 재평가 논문으로 분리** |
|---|---|---|
| **중심 주장 (한 문장)** | "8-gene 갑상선 계통 상태는 driver 축과 독립적인 예후 인자이며, RAI-치료군 내에서 구조적 재발과 연관되나 이는 예후적 신호일 뿐 RAI 반응을 예측한다는 증거는 아니다." | "Driver mutation class는 RAI 불응성을 재현성 있게 구분하지 못하는 반면(3개 코호트, n=294), 발현/단백체 기반 분자 서브타입은 구분한다(Zhang, P=1.4×10⁻⁷) — 그러나 이를 8-gene 패널로 직접 검증한 시도는 오늘까지 4번 모두 null 또는 검증력 부족이었다." |
| **본문(main text)에 들어가는 것** | ① Limitation 6 재작성(아래 §3) ② TCGA course-level RAI 코호트를 새 Results 문단으로 1개 추가 — 1차 반응 null + 구조적 재발 prognostic-within-treated (Firth OR, event 수, interaction 불가 명시) ③ Fig 6D 재조정(위 긴장 해소 후 확정) ④ 제목/Abstract에서 "predicts/marks radioiodine-refractory(-ness)" 어휘 제거(§4) | ① 데이터 랜드스케이프 감사(Fig 1) ② TCGA biotab RAI 자원 자체를 자원논문 형태로 소개(Fig 2) ③ driver-축 3코호트 메타분석 null(Fig 3) ④ 패널/축에 대한 4개 직접검증 compendium — 전부 null/underpowered(Fig 4) ⑤ TCGA 구조적 재발 prognostic 신호(Fig 5, Paper 1과 동일 분석을 인용만 하고 재도출하지 않음 — 아래 §5 참조) ⑥ Zhang 단백체 서브타입 양성(Fig 6) |
| **Supplement로 가는 것** | GSE138042, MERAIODE, driver-메타분석(1문장 인용), Zhang(1문장 인용) — 모두 "우리 패널이 아니라 문헌상 다른 접근"이라는 각주 수준. index-course sensitivity 표. | index-course sensitivity, Firth vs ML 비교/LOO/bootstrap, GSE151179 음성대조군 상세, GSE138042 per-gene 표, MERAIODE HTG 패널 유전자 커버리지, Siraj/Boucai 코호트 상세(circularity flag 포함), 16_SCREENED_EXCLUSIONS 전체 표, power/design-target 그림 |
| **제외되는 것** | Zhang 단백체 서브타입 전체 스토리, driver-축 메타분석 forest, 데이터 랜드스케이프 감사, MERAIODE의 assay-gap 진술(1문장 초과 분량), 16개 후보 스크리닝 로그 전체 | TCGA/MSK OS·PFI, methylation, fusion 농축, single-cell 검증 등 Paper 1 고유 결과 전부(중복 금지, §5) |
| **주요 신규 그림 수 (main)** | **0개.** 기존 Fig 5/6에 패널 1–2개 추가 또는 교체로 흡수 가능(신규 Figure 번호 불필요) | **6개** — 모두 오늘 자로 이미 생성된 파일 재사용 가능(§Deliverable 2 참조), 신규 조합 작업만 필요 |
| **오늘 데이터로 제출 가능한가** | **예 — 오히려 이전보다 제출 가능성이 높아짐.** 이전에는 "RAI 결과 데이터가 있는 코호트가 없다"는 (이제는 틀린) 문장으로 취약점을 방어했다. 오늘 자로 그 데이터가 실제 존재하고 직접 검증했다는 사실 자체가 리뷰어의 가장 뻔한 공격("GDC를 봤는데 왜 무시했나")을 선제 차단한다. 필요한 작업은 전부 편집·수치 교체이며 신규 분석이 필요 없다. | **예, 최소 제출가능 버전(scope i)으로 오늘 데이터만으로 가능.** 단 아래 §5 목표 저널 논의에서 보듯 NC급은 아니고 npj Precision Oncology/Thyroid/JCEM급 현실적. |
| **리스크** | (a) 제목·Abstract의 "가장 강한 clinical hook"을 스스로 약화시켜야 함 — 저자 설득 필요. (b) Fig 6D 재조정 미해결 시 그대로 제출하면 리뷰어가 GSE151179 내부 모순을 발견할 위험. (c) event 9/13의 Firth 결과를 새 Results 문단으로 추가하면 "그래서 표본이 너무 작지 않냐"는 새 공격면이 열림 — 그러나 이는 정직한 신호이므로 리스크가 아니라 정상적인 peer review 절차로 흡수 가능. | (a) 이 논문 자체가 "null 결과 모음"으로 읽혀 저널 임팩트가 제한됨(§Deliverable 2 목표 저널 참고). (b) Zhang 수치를 오픈 PDF에서 pdftotext로 추출했다는 방법론적 취약점 — 리뷰어가 신뢰성 문제 제기 가능(방법 절에서 QC 명시 필요). (c) 3-코호트 메타분석 I²=63% — 강하게 주장하면 안 됨(이미 05_EXTERNAL_COHORT_AUDIT.md가 금지 문구를 명시). |
| **중복출판 방지 방법** | Paper 1은 TCGA course-level 데이터에 대한 **패널 자체의 직접 검증**(1차 반응 null + 구조적 재발 prognostic)만 소유. Driver-축 메타분석, Zhang, GSE138042/MERAIODE의 세부 수치는 인용 1문장으로만 언급하고 "companion paper (in preparation)"로 참조, 수치 재제시 금지. | Companion paper는 Paper 1의 고유 결과(OS/PFI HR, ARI, methylation, fusion enrichment, single-cell)를 재도출하지 않고 인용만 한다. **TCGA 구조적 재발 prognostic 결과는 Paper 1이 1차 소유권을 가지며, companion paper는 이를 "다른 논문에서 보고된 결과"로 1회 인용**해 "그래서 4번의 직접 검증 중 유일하게 양성이었던 신호도 RAI-특이적이지 않다"는 문장에 활용한다. 두 논문 모두 상호 각주로 non-overlap을 명시("이 두 편은 동일 프로젝트의 자매 논문이며 결과 중복이 없음을 저자가 확인함"). |

---

## 3. 최종 권고와 그 이유

**권고: Strategy B (분리) + Strategy A의 최소 필수 편입.** 즉 이분법이 아니라, Paper 1에는
"패널이 검증되었는가"에 대한 정직한 답변만(§2 Strategy A 열) 넣고, "RAI 불응성 예측이라는 더 큰
질문에 대해 공개 데이터가 지금 무엇을 말해주는가"라는 훨씬 큰 별도 서사는 companion paper로
분리한다.

이유:
1. **범위 불일치.** Strategy A 열의 내용은 8-gene 패널과 TCGA 코호트 하나에 관한 것이다.
   Strategy B 열의 내용은 3개의 완전히 다른 코호트, 8-gene 패널을 전혀 쓰지 않는 분석(Zhang
   단백체, driver-메타분석), 그리고 데이터 랜드스케이프 감사 자체를 포함한다. 이걸 전부 DM1
   원고에 넣으면 원고의 정체성이 "driver-orthogonal state 발견 논문"에서 "RAI 예측 바이오마커
   문헌 리뷰 논문"으로 흔들린다.
2. **저널 적합성.** Paper 1은 이미 Nature Communications를 목표로 8개월 가까이 다듬어진 원고다
   (SUBMISSION_CONTROL 참조). companion paper의 실제 콘텐츠(§Deliverable 2에서 상세)는 이
   급의 저널에 맞지 않는 "재평가/자원" 성격이 강하다 — 억지로 합치면 Paper 1의 임팩트까지
   같이 희석된다.
3. **오늘 데이터가 실제로 하는 일이 다르다.** #1의 null은 Paper 1의 "RAI-refractory" 언어를
   **약화**시키는 정보다. #7/#8(driver 축 null, 발현 축 양성)은 Paper 1의 근본 전제("driver
   축보다 분화/발현 축이 옳다")를 **강화**하는 정보다. 하나는 축소해야 할 주장이고 다른 하나는
   별도로 키울 가치가 있는 주장이다 — 같은 논문 안에 넣으면 이 방향 차이가 리뷰어에게 안
   보인다.
4. **분리하면 두 논문 모두 더 방어 가능해진다.** Paper 1은 이제 "RAI 반응을 예측한다"는
   말을 하지 않으므로 그 축의 공격에 노출되지 않는다. Companion paper는 "이건 8-gene 패널
   검증 논문이 아니라 공개 데이터 현황 감사"라고 처음부터 선언하므로 "패널을 안 썼다"는
   비판이 성립하지 않는다.

---

## 4. 제목 재평가

### 4-1. 현재 제목의 문제

리포지토리에는 두 개의 유효 제목 문자열이 있다:

- `SUBMISSION_CONTROL/DM1_PROJECT_STATE.md:13` 및 `STATUS_INSILICO_FINALIZATION_2026_07_30.md:21`:
  **"A thyroid-lineage state predicts radioiodine-refractoriness in BRAF V600E-mutant papillary
  thyroid cancer"**
- 최신 실제 초고 `NATURE_COMMUNICATIONS_FULL_DRAFT_v2_2026_07_08.md:29`:
  **"A thyroid-lineage state marks radioiodine-refractory biology in BRAF V600E-mutant papillary
  thyroid cancer"**

두 제목 모두 "predicts radioiodine-refractoriness" 또는 그 완곡화("marks ... refractory
biology")를 포함한다. `SUBMISSION_CONTROL/CLAIM_LANGUAGE_MATRIX.md`는 이미 2026-07-30에 이
문구를 "vulnerable"로 지정했지만, 그 시점의 근거는 **"직접 측정 데이터가 없다"**는 부재
논증이었다.

**오늘 자 데이터는 이 논증을 근본적으로 바꾼다.** 이제는 "데이터가 없어서 증명 못 한다"가 아니라
**"데이터가 있고, 직접 검증했고, null이다"**다:

- 1차 RAI 반응(가장 직접적인 RAI-response endpoint): d=−0.03, 95% CI가 ±0.5 효과를 배제(§1 #1).
- Treatment×panel interaction: non-RAI군 event 2건뿐이라 **추정 자체가 불가능**(§1 #2 각주).
- GSE151179의 진짜 RAI-effect 필드(met-site uptake): null, purity로 완전히 설명됨(§1 #4).
- 두 개의 추가 코호트(GSE138042, MERAIODE)에서도 동일 방향이나 null(§1 #5, #6).

### 4-2. 판단

**"predicts radioiodine-refractoriness" 및 이와 동의어인 "marks radioiodine-refractory
biology"는 오늘 자 증거로 지지되지 않는다.** 이는 소극적 결론("아직 증명 안 됨")이 아니라
적극적 결론("직접 검증에서 null")이다. 이 어휘를 제목에 유지한 채 제출하면, 리뷰어가 GDC
biotab(공개 데이터, 특별한 접근 권한 불필요)을 확인하고 정확히 이 null을 재현할 위험이 실재하며,
이는 major revision이 아니라 신뢰성(scientific integrity) 문제로 읽힐 수 있다.

원고가 실제로 지지하는 것: (a) driver 축과 독립적인 갑상선 계통/분화 상태의 존재(ARI=0.92,
driver-only=−0.007), (b) 이 상태의 예후적 의미(pooled OS HR=2.53), (c) BRAF+ 환자군 내 PFI와의
연관성(단일 코호트, hypothesis-generating으로 이미 명시됨), (d) 이제 추가로 — RAI-치료군 내
구조적 재발과의 연관(Firth-보정, robust하나 RAI-특이적이지 않음). 제목은 이 네 가지 중
검증된 것에 근거해야지, 검증에 실패한 다섯 번째 주장(RAI 반응 예측)에 근거해서는 안 된다.

### 4-3. 제목 후보 3계열, 각 6개 이상

각 후보에 대해 **주장(commits to)**, **오늘 데이터가 뒷받침하는가**, **잔여 리스크**를 표기한다.

#### 계열 A — 분화 상태(differentiation-state) 중심

| # | 제목 | 주장 | 지지 여부 | 잔여 리스크 |
|---|---|---|---|---|
| A1 | *A driver-independent thyroid-differentiation state defines an aggressive subtype of BRAF V600E-mutant papillary thyroid cancer* | 분화 상태 존재 + driver 독립성 + (BRAF+ 내) 공격적 임상 경과 | **지지됨** — ARI=0.92/−0.007, PFI HR=0.66 p=0.013(BRAF+ 한정) | "aggressive"가 임상적 중증도 언어이므로 Abstract에서 반드시 단일 코호트·hypothesis-generating 단서를 동반해야 함. 낮음 |
| A2 | *An eight-gene thyroid-lineage state resolves driver-orthogonal risk in papillary thyroid cancer* | 패널이 driver-무관 위험(risk)을 구분함 | **지지됨** — OS, BRAF+ PFI 모두 근거 있음 | "risk"가 다소 모호 — Abstract에서 즉시 OS/PFI로 구체화 필요. 낮음. **1차 권고안** |
| A3 | *Loss of a thyroid-lineage transcriptional state marks aggressive, driver-independent papillary thyroid cancer* | 상태 손실 = 공격적 경과, driver 무관 | 지지됨, A1과 유사 | "marks"라는 동사 자체는 문제 없음(원인이 아니라 상태를 지칭) — 다만 A1과 동일하게 aggressive 단서 필요. 낮음 |
| A4 | *A compact thyroid-lineage panel identifies a driver-orthogonal prognostic state in papillary thyroid cancer* | 패널 → 예후 상태 식별, driver 직교성 | 지지됨 | "identifies"는 "predicts"보다 약한 동사라 안전. 낮음 |

#### 계열 B — 구조적 위험(structural-risk) 중심

| # | 제목 | 주장 | 지지 여부 | 잔여 리스크 |
|---|---|---|---|---|
| B1 | *A thyroid-lineage state predicts structural recurrence risk after radioiodine in BRAF V600E-mutant papillary thyroid cancer* | RAI 이후 구조적 재발을 "예측" | **부분 지지, 동사가 위험함** — treatment×score interaction 추정 불가이므로 "predicts...after radioiodine"는 RAI-특이성을 암시해 오늘 데이터가 명시적으로 금지하는 주장과 같은 함정에 다시 빠짐 | **높음 — 비권고.** §1 #2/#3의 핵심 caveat("RAI-specific이 아님")를 제목 차원에서 위반 |
| B2 | *A thyroid-lineage state is associated with structural disease recurrence in radioiodine-treated papillary thyroid cancer* | RAI-치료군 내 구조적 재발과 "연관"(동사를 association으로 낮춤) | **지지됨** — Firth OR 0.19–0.24, LOO/bootstrap 안정 | 중간 — event 9/13으로 작음, Abstract에서 즉시 명시 필요. "radioiodine-treated"가 population 기술이지 RAI 효과 주장이 아님을 Methods에서 분명히 해야 함 |
| B3 | *A thyroid-lineage state stratifies post-treatment structural risk in papillary thyroid cancer, independent of driver mutation* | 치료 후 구조적 위험 계층화 + driver 독립성 | 지지됨, B2와 유사하나 driver-직교성을 제목에 명시해 더 안전 | 중간 — "post-treatment"가 여전히 RAI를 암시할 수 있어 Methods 명확화 필요 |

> **계열 B에 대한 종합 판단:** B1은 §4-2의 판단을 정면으로 위반하므로 제외 권고. B2/B3는
> 사용 가능하지만, 이 신호(event 9/13)를 제목의 헤드라인 주장으로 올리는 것은 원고의 무게중심을
> 지나치게 작은 표본으로 옮기는 편집적 판단이 필요하다 — 1차 권고안으로는 채택하지 않는다
> (§4-4).

#### 계열 C — driver-독립적 생물학(driver-independent-biology) 중심

| # | 제목 | 주장 | 지지 여부 | 잔여 리스크 |
|---|---|---|---|---|
| C1 | *An eight-gene circuit reveals fusion-driven, epigenetically silenced iodine-handling loss in BRAF/RAS-negative thyroid cancer* | 순수 생물학: fusion 농축 + methylation 침묵 + 요오드 처리 소실, RAI 반응 주장 전혀 없음 | **가장 강하게 지지됨** — FUS-01/02/03, M-01~M-12 모두 audit-locked | **가장 낮음.** 다만 임상적 훅(clinical hook)이 사라져 NC 리뷰어의 흥미를 낮출 위험 — 이는 과학적 리스크가 아니라 임팩트 포지셔닝 리스크 |
| C2 | *Fusion-driven epigenetic silencing defines a driver-orthogonal iodine-handling-low state in papillary thyroid cancer* | C1과 동일 취지, "circuit" 대신 "state" | 지지됨 | 낮음. C1보다 약간 더 임상적으로 읽힘 |
| C3 | *A driver-orthogonal, epigenetically silenced iodine-handling state stratifies survival in papillary thyroid cancer* | 기전(driver-orthogonal, epigenetic silencing) + 생존 계층화(OS 결과만 사용, BRAF-interaction 논쟁 회피) | **지지됨** — OS pooled HR=2.53는 RAI 논쟁과 완전히 독립적인 결과 | 낮음. **A2와 함께 공동 1위 후보** |

### 4-4. 최종 권고

1차 권고: **A2** (*An eight-gene thyroid-lineage state resolves driver-orthogonal risk in
papillary thyroid cancer*). 이유: (i) 오늘 자 데이터가 요구하는 것 — RAI 반응 예측 주장의
완전한 제거 — 를 만족시키면서도, (ii) 원고가 실제로 갖고 있는 가장 강한 두 결과(OS, BRAF+ PFI
interaction)를 여전히 제목이 가리키는 범위 안에 남겨 두어 Abstract·Discussion에서 clinical
hook을 유지할 수 있다. 2차 대안: **C3**(생존만으로 제목을 완결시켜 BRAF-interaction 논쟁 자체를
제목 레벨에서 회피, 가장 방어적). 3차 대안(만약 저자가 RAI 언어를 제목에 남기고 싶다면): **B2**,
단 이 경우 Abstract 첫 문장에서 "association, not RAI-specific prediction; interaction not
estimable"을 명시하는 것을 조건으로 한다.

**비권고: 현재 두 제목 문자열(둘 다) 및 B1.** "predicts radioiodine-refractoriness"와
"marks radioiodine-refractory biology"는 오늘 자 직접 검증 결과와 정면으로 배치되므로 그대로
제출하면 안 된다.

제목은 저자가 최종 결정한다. 이 문서는 변경을 실행하지 않았다.

---

## 5. Paper 1 편입 시 구체적 실행 체크리스트 (Strategy A 최소 세트)

편집 실행 전 반드시 완료해야 하는 선행 항목:

1. **[BLOCKING]** Fig 6D / GSE151179 anchor 재조정. 현재 원고가 쓰는 d≈−1.0 비교와 오늘
   재검증한 d=+0.37(met-site uptake, null) 비교가 서로 다른 필드/짝짓기를 쓰는지 코드
   레벨에서 확인하고, 어느 것이 1차 anchor 주장인지 확정한다. 확정 전에는 이 문장을 어떤
   형태로도 손대지 않는다(evidence-lock 규칙).
2. Limitation 6 재작성. 현재 문장(`06_discussion.md:47`, `10_full_manuscript_compiled.md:511`)
   — *"no cohort in this study contains prospective post-thyroidectomy RAI outcome data linked
   to per-sample DM1 calls"* — 는 이제 사실이 아니다. TCGA GDC BCR Biotab에 course-level RAI
   임상 데이터가 존재하며 오늘 그것을 이용해 직접 검증했다. 새 문장은 "검증했고 null이었다"는
   내용과 그 power(§1 #1)를 담아야 한다.
3. 새 Results 문단 1개 추가(신규 Figure 불필요, 기존 Fig 5/6에 패널 흡수 또는 텍스트만):
   1차 반응 null(§1 #1) + 구조적 재발 prognostic-within-treated(§1 #2, #3, Firth 수치,
   event 수, interaction 미추정 명시)를 함께 제시해 "예후적이지 예측적이 아니다"라는 대비를
   본문에서 스스로 드러낸다 — 리뷰어가 지적하기 전에.
4. 제목/Abstract에서 "predicts/marks radioiodine-refractory(-ness)" 계열 어휘 제거,
   §4-4 권고안 중 하나로 교체(저자 승인 필요).
5. `CLAIM_LANGUAGE_MATRIX.md`의 "RAI link" 행 갱신 — "predicts RAI uptake" 금지 문구는
   유지하되, 새로 "RAI-treated 코호트 내 구조적 재발과 연관, RAI-특이성 미확정"이라는 허용
   문구를 추가.
6. GSE138042·MERAIODE·driver-메타분석·Zhang은 각 1문장 이하로만 인용(§2 표의 "제외" 열),
   본문 수치 전체를 재제시하지 않는다 — companion paper 출간 전까지는 "unpublished
   companion analysis" 또는 사전 인쇄본 인용으로 처리.

이 항목들은 모두 **문서 편집**이며 신규 분석을 요구하지 않는다. 1번(Fig 6D 재조정)만 primary
output 확인이 필요하므로 가장 먼저 처리해야 한다.

### 5-1. 2차 감사(claims registry)가 확인한 추가 구체 항목

`06_CLAIMS_SUMMARY.md` / `06_CLAIMS_REGISTRY.csv`가 원고 전체(NC v2 본문, 구 compiled본,
cover letter, reviewer QA, figure caption NC본)를 대상으로 101건의 RAI/predictive/prognostic
관련 문장을 file:line 단위로 감사했다. 이 중 **P0(fatal) 15건**이 위 1–4번 항목과 직접
겹치거나 이를 구체화한다. Strategy A 실행 시 다음도 함께 처리해야 한다(§1–4에 없던 신규 항목만
발췌):

7. **패널의 이름 자체가 문제.** "8-gene **RAI-responsiveness panel**"이라는 명칭이
   `10_full_manuscript_compiled.md:67,158,261`, `08_cover_letter.md:21`(CL-004, 현재
   cover letter에 남아있는 stale boilerplate — 2026-07-30 개정 시 누락됨),
   `05_figure_captions_NC.md:21`(FIGNC-001)에 남아 있다. NC v2 본문은 이미 "thyroid
   differentiation and iodine-handling axis/panel"로 고쳐져 있으나 나머지 문서가 따라가지
   못했다. 전부 동일 용어로 통일. CL-004의 "validated"는
   `CLAIM_LANGUAGE_MATRIX.md` 규칙대로 "reproduced"로 교체.
8. **"predictive interaction" / "treatment-selection biomarker" 프레이밍 전체 재라벨링
   (원고에서 가장 강한 단일 overclaim, NCv2-037).** 이것은 DM1-score × **BRAF genotype**
   interaction(p=0.022, 원고 자체 결과, 유효함)을 DM1-score × **radioiodine treatment**
   interaction(오늘 자 감사로 추정 불가로 확정됨, §1 #2)과 뒤섞어 부르는 문제다.
   `NATURE_COMMUNICATIONS_FULL_DRAFT_v2_2026_07_08.md`의 §3 표제(line 74), line 76, 78,
   80, 142, Abstract(line 40), Fig 7 제목(line 260)이 전부 해당. 수정 방향: "genotype-
   context-dependent prognostic" 표현으로 통일하고, RAI 치료-선택 역할은 treatment×score
   interaction이 실제로 검정되었으며 추정 불가였다는 문장을 명시적으로 동반해야만 언급
   가능.
9. **Fig 6D 캡션의 인과·정체성 초과 서술.** `05_figure_captions_NC.md:125`
   (FIGNC-005) — "...providing an external clinical anchor that the DM1 axis **names the
   same biology that defines RAI-refractory progression**"는 본문(line 124)이 이미 쓰고
   있는 절제된 표현("we interpret it as a biological anchor rather than as a validated
   response-prediction model")보다 더 나간다. 이는 §1의 Fig 6D 수치 재조정(BLOCKING 항목)과
   **별개의, 언어 수준의 문제**이며 수치가 어떻게 정리되든 캡션은 본문 수준의 절제된 어휘로
   맞춰야 한다.

이 세 항목은 §2 실행 체크리스트 3–4번과 함께 한 번의 편집 패스로 처리 가능하다(전부 어휘 교체,
신규 분석 없음). 전체 101건 목록과 정확한 file:line은 `06_CLAIMS_REGISTRY.csv`를 참조한다.
