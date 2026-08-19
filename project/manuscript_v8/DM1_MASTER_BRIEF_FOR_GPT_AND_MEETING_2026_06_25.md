---
title: "DM1 / DM2 종합 브리프 — GPT 컨텍스트 + 오늘 미팅 준비"
date: 2026-06-25
audience: GPT 브리핑용 · 오늘 미팅 준비용
purpose: 연구 핵심, 현재 분석 결과, 강민수 선생님 피드백, 미팅 액션 아이템, 모든 deploy 상태를 한 파일로 종합
companion_docs:
  - MANUSCRIPT_REORG_6FIG_PLAN_2026_06_05.md  (6-figure proposal)
  - EDITORIAL_REVIEW_RUTHLESS_CUT_2026_06_05.md  (NC editor 관점 5-figure cut)
  - 04_results.md  (current results draft)
  - 05_figure_captions_NC.md  (current captions)
---

# DM1 / DM2 종합 브리프 — 2026-06-25

> **한 줄 요약.** 우리는 갑상선암에서 RAI 가 잘 들을 사람과 안 들을 사람을, 갑상선 분화 / iodine-handling 관련 **8 개 유전자** 로
> 미리 나누는 분자 axis (DM1 / DM2) 를 정의했고, 여러 공개 코호트에서 재현성을 확인했다.
> **임상 메시지는 단순히 "더 빨리 항암으로 보낸다" 가 아니라, "RAI 가 안 들을 환자에게 불필요한 고용량 RAI 를 반복하지 않게 한다" 가 더 강하다.**

---

## 0. 이 문서의 사용법

본 문서는 두 가지 용도다.

1. **GPT 에게 전체 컨텍스트를 한 번에 브리핑** — §13 의 GPT prompt 를 그대로 복사해서 GPT 에 붙여넣으면 됨.
2. **오늘 미팅 (2026-06-25 오전 9시) 준비 노트** — §1~§12 를 미팅 전 5분에 훑으면 흐름 잡힘.

---

## 1. 연구 핵심 아이디어 (clinical message)

현재 갑상선암 환자는 수술 후 RAI 치료를 해보고, 반응이 안 좋으면 그제야 고위험 / refractory 로 판단된다.
**6 개월 ~ 1 년 이상 지연**될 수 있고, 특히 RAI 가 흡수되지 않는 환자에게 고용량 RAI 를 반복 투여하면
**골수억제, 백혈구 감소, 장기 부작용** 위험이 생긴다.

따라서 본 연구의 임상적 메시지는 두 줄이다.

1. **★ 주 메시지 — 강민수 선생님 강조.**
   > **"RAI 가 안 들을 환자에게 불필요한 고용량 RAI 를 반복하지 않게 한다."**
   > 갑상선암은 장기 생존 환자가 많아서, "불필요한 고용량 RAI 독성 회피" 가 환자에게 훨씬 직접적이고 임상적으로 설득력 있음.

2. **부 메시지.**
   > "RAI 가 안 들을 환자에게서는 systemic / targeted therapy 로의 전환이 지연되지 않도록 한다."

논문 framing 에서는 1 번을 더 앞세운다.

---

## 2. 8-Gene Signature — 공식 symbol 통일

### 공식 HGNC symbol (분석 / 코드 / 테이블 기준)

| # | Official symbol | Alias / clinical name | 기능 범주 |
|---|---|---|---|
| 1 | **TG**       |  thyroglobulin                          | thyroid hormone precursor · storage |
| 2 | **TPO**      |  thyroperoxidase                        | iodide oxidation · organification |
| 3 | **TSHR**     |  TSH receptor                           | TSH signaling |
| 4 | **SLC5A5**   |  **NIS**  (sodium/iodide symporter)     | iodide uptake (핵심) |
| 5 | **DIO1**     |  type 1 deiodinase                      | thyroid hormone metabolism |
| 6 | **PAX8**     |  —                                       | thyroid lineage transcription factor |
| 7 | **NKX2-1**   |  **TTF-1**                              | thyroid lineage transcription factor |
| 8 | **FOXE1**    |  **TTF-2**                              | thyroid lineage transcription factor |

### 논문 / 발표용 병기 표기
> TG, TPO, TSHR, **SLC5A5 / NIS**, DIO1, PAX8, **NKX2-1 / TTF-1**, **FOXE1 / TTF-2**

### 패널 구조 (개념)
- **Effector 5 개** — TG · TPO · TSHR · SLC5A5 · DIO1   (RAI 흡수 / iodine handling 의 실제 도구)
- **Transcription factor 3 개 (TF₃)** — PAX8 · NKX2-1 · FOXE1   (갑상선 분화 마스터)

### Reviewer 가 가장 먼저 공격할 부분
> "왜 하필 이 8 개 gene 인가?"

→ 답변용 단계 (manuscript 에 diagram 필요):

1. 문헌 기반 후보 풀 수집 (RAI uptake · iodine metabolism · thyroid differentiation 관련)
2. 기능별 분류 (effector / TF / metabolism)
3. 공개 cohort 에서 재현성 평가
4. Redundancy / collinearity / model parsimony 기준 축소
5. 최종 8-gene 선정
6. 16-gene 기존 signature / random 8-gene / 4-gene subset 등과 비교
7. 최종 signature 가 가장 robust 라는 근거 제시

→ 메시지: <em>"ML 이 아무 gene 이나 고른 게 아니라, thyroid biology 에 기반한 후보군에서 출발했고, 공개 데이터로 검증해서 compact 8-gene 으로 정리했다."</em>

---

## 3. DM1 / DM2 정의

DM1, DM2 는 **기존 임상 분류가 아니라 이번 연구팀이 만든 molecular state classification** 이다.

| 군 | 분자 상태 | 임상 해석 (조심스럽게) |
|---|---|---|
| **DM2** | thyroid differentiation 유지 · iodine-handling on | RAI-responsive-like |
| **DM1** | dedifferentiated · iodine-handling off | RAI-refractory-like · aggressive 가능성 |

### 이름 재검토 (논문용)

"DM1 / DM2" 는 내부 코드명처럼 느껴짐.  논문에서는 다음으로 바꿔도 좋음:

| 후보 | 의미 |
|---|---|
| RAI-sensitive-like / RAI-refractory-like | 직관적 |
| Iodine-avid / Iodine-poor | 임상 표현 |
| Thyroid-differentiated / De-differentiated | 분자 중립 |
| **Iodine-handling-high / Iodine-handling-low** | **★ 가장 안전 (논문용)** |

→ 권장: **Iodine-handling-high / Iodine-handling-low** (논문 본문에서 "previously termed DM2 / DM1" 으로 cross-reference).

---

## 4. 현재 분석 결과 — manuscript audit-locked 핵심 수치

### Discovery cohort
- **TCGA-THCA n = 504** primary PTC · 8-gene KMeans (k=2) → **DM1 28.4 % / DM2 71.6 %**
- Pan-genome top-5000 MAD clustering 으로도 **ARI = 0.92** 재현 → 패널 artifact 아님

### Driver-orthogonality (핵심 방어)
- BRAF · TERT · KRAS · NRAS · HRAS 단일 feature AUC ≈ 0.5 (모두 chance level)
- BRAF V600E 양성 종양의 DM1 prevalence = **0.7 %**
- RAS-mutant 종양의 DM1 prevalence = **96.4 %**
- BRAF / RAS-음성 dark matter = **DM1 49.1 % / DM2 50.9 %** ← <b>분층화 가치 가장 큰 환자군</b>

### Epigenetic mechanism (TCGA HM450, n = 503)
- TPO promoter β Cohen's d = **2.30**, p = 1.9 × 10⁻¹⁸
- DM1 mean 8-gene β = **0.385** vs DM2 0.253  (+52 %)
- Per-driver-class: BRAF (0.37) ≈ RET fusion (0.39) ≫ RAS (0.27)  → MAPK 활성에 비례, driver 정체성 아님

### 외부 검증 — 19 cohort cross-validation
- **Master cross-cohort forest**: 14 entries × 11 distinct cohorts, **평균 Cohen's d = 2.81, median = 2.37, 모두 ≥ 1.56**
- **Per-gene × cohort × contrast matrix**: 8 genes × 10 contrast × 4 cohort = **80 / 80 cell direction-consistent**
- Lee 2024 (GSE213647, n = 632 Korean PTC FFPE) within-cohort KMeans Cohen's d = **5.93**
- K2 (PRJEB11591, n = 260) within-cohort d = 1.94 (calibration mismatch 해결 후)
- GPL570 4-cohort (n = 205): 4 / 4 ρ ≤ −0.84 (Western microarray platform)
- **Mun 2025 proteogenomic (n = 336)**: thyroid_diff d = −1.91, 7/7 panel 유전자 sign-consistent → 단백질 수준 검증
- Landa 2016 PDTC + ATC 의 silenced gene list 와 **5 / 8 유전자 overlap** (TG · TSHR · TPO · PAX8 · DIO1) → reverse-causality lock

### 임상 의미
- TCGA + MSK-IMPACT pooled OS HR = **2.53 [1.31, 4.89]**, I² = 0 %  (후향 데이터)
- ATA 2015 / 2025 중간 위험군에 DM1 과대표현  (RAI 결정 불확실 zone)
- GSE151179 post-RAI refractory 종양이 transcriptionally DM1 state 와 일치 (d ≈ −1.0, MW p ≈ 10⁻⁴)
- FFPE vs FF score 분포: Kolmogorov-Smirnov p = **0.44** → FFPE 호환 가능

### 단일세포
- Lu 2023 thyrocyte UMAP (n = 14,624) — KRT8 ∩ KRT19 ∩ EPCAM 필터 후에도 DM score gradient 명확 → **stromal / immune confound 아닌 thyrocyte 본질적 신호**
- Pu 2021 paired 6 환자 모두 per-patient r = 0.798–0.886 (Bonferroni p < 10⁻¹⁰)

---

## 5. 강민수 선생님 피드백 4 가지

### A. 임상 메시지 수정
- 처음: "RAI refractory 를 빨리 알아내서 대체 치료를 빨리 한다"
- **★ 수정: "RAI 가 안 들을 환자에게 불필요한 고용량 RAI 를 하지 않게 해준다"**
- 이유: 갑상선암은 오래 사는 환자가 많아서 부작용 회피가 더 직접적

### B. 8 개 gene 선정 과정이 논리적으로 보여야 함
- "왜 하필 8 개?" reviewer 공격 차단 위해 §2 의 7-step 흐름 figure 필요
- 16-gene 기존 signature, random 8-gene, 4-gene subset 과 비교 필수

### C. 실제 임상 NGS panel 에 8 개 gene 이 들어있는지 확인 필요
- 분당병원 NGS panel = whole genome / exome 이 아닌 제한 panel
- 시나리오:
  1. **8 개 모두 포함** → 바로 임상 NGS 기반 validation 가능
  2. **일부만 포함** → 4-gene / 5-gene / 6-gene reduced model 만들어야 함 → "임상 panel 에서도 근사 예측 가능" 으로 스토리 조정
  3. **거의 안 포함** → 별도 RNA panel · targeted assay · IHC validation 으로 방향 전환

### D. IHC validation 이 매우 중요
- DNA / RNA 분석은 논문에는 좋지만, 실제 임상 적용성은 IHC 가 훨씬 강함
- 수술 조직 슬라이드 많음, IHC 는 빠르고 싸게 가능
- 장기 비전:  **8-gene RNA signature → protein expression / IHC score → RAI response prediction**
- IHC 로 RAI 반응성 예측 가능하면 **진단법 특허**까지 고려 가능

---

## 6. 공동연구 역할 분담

### 강민수 선생님
1. 8 개 gene 이 병원 NGS panel 에 들어있는지 확인
2. 갑상선암 환자 중 NGS 검사 받은 환자 수 확인
3. IRB 준비
4. 약 30 명 규모 internal validation 가능성 확인
5. RAI response / refractory 여부 + 임상 outcome 정리

### 승호님 (나)
1. 8-gene signature 선정 근거 정리
2. Gene alias / symbol 통일 (§2)
3. 공개 cohort 분석 결과 정리 (§4)
4. DM1 / DM2 classification 결과 정리
5. 기존 signature 와 비교 (16-gene 등)
6. 4-gene / 6-gene reduced model 가능성 테스트
7. NGS panel subset model 준비
8. 논문 figure / supplementary figure 정리

### 유 교수님 / PI
1. 전체 manuscript 방향 설정
2. IRB 문서 초안 준비
3. 임상 narrative 정리
4. 공동저자 구조 정리
5. 6 월 25 일 오전 9 시 회의 준비

---

## 7. 논문 스토리라인 (현 잠정)

### Title 후보

1. **A Compact Thyroid Differentiation Signature Predicts Radioiodine Refractoriness in Thyroid Cancer Across Multi-Platform Public Cohorts**
2. **A Clinically Translatable 8-Gene Signature for Avoiding Ineffective High-Dose Radioiodine Therapy in Thyroid Cancer** (★ 강민수 메시지 반영)
3. A driver-orthogonal differentiation axis defines an aggressive subtype of papillary thyroid cancer
4. Hidden in plain sight — a driver-independent thyroid differentiation axis underlying radioiodine failure

### Main claim (안전한 표현)

> A compact 8-gene thyroid differentiation and iodine-handling signature stratifies thyroid cancer patients into RAI-avid and RAI-refractory-like molecular states across multiple public cohorts and **may help avoid ineffective high-dose RAI therapy.**

### Introduction 메시지
- 대부분 갑상선암은 예후 좋음 · 일부는 aggressive
- RAI 는 중요 치료지만 일부는 primary resistance / poor uptake
- 현재는 RAI 치료 후 반응 보고 판단 → 지연 + 부작용
- 특히 RAI 안 듣는 환자에게 고용량 RAI 반복은 장기 독성
- 따라서 치료 전 / 초기에 RAI non-responder 예측할 compact, clinically translatable biomarker 필요

### Results 구조 (8-step)
1. Literature-guided iodine-handling / thyroid differentiation gene set 구성
2. 8-gene RAI signature 도출
3. DM1 / DM2 두 molecular state 식별
4. Cross-cohort validation (GEO + TCGA + 문헌 cohort)
5. 기존 16-gene 및 광역 signature 와 비교
6. RNA-seq / microarray platform 간 robust
7. Reduced-panel 분석 (clinical NGS 호환성)
8. Aggressive histology · ATC/PTC · 생존 · RAI response 와의 association
9. (Future) institutional cohort + IHC validation

---

## 8. 오늘 미팅 (2026-06-25 오전 9시) 체크리스트

### 우선순위 액션
1. **8 개 gene 최종명 통일** — HGNC symbol + alias + Ensembl ID + Entrez ID  (§2 표 그대로 사용)
2. **임상 NGS panel 포함 여부 확인** — 강민수 선생님에 즉시 질문
3. **8-gene 선정 flowchart 도식 만들기** — literature → candidate pool → filtering → model selection → validation
4. **DM1 / DM2 이름 재검토** — Iodine-handling-high / -low 권장
5. **"불필요한 고용량 RAI 회피" 임상 메시지 강화** — abstract / introduction / discussion 에 반영
6. **IRB 초안 작성** — retrospective chart review · NGS data · RAI dose/response · recurrence/metastasis · 혈액수치 (백혈구 감소 등)
7. **30 명 validation endpoint 정의** (§9)

### 오늘 미팅에서 말할 핵심 멘트 (그대로 사용 가능)

> 제가 정리해보니, 이 연구의 임상적 메시지는 두 가지로 잡는 게 좋을 것 같습니다. 첫째는 RAI 가 잘 듣지 않을 환자에게 불필요한 고용량 RAI 를 반복하지 않도록 해서 골수억제나 장기 부작용을 줄이는 것이고, 둘째는 예후가 나쁜 aggressive thyroid cancer 에서 RAI 를 오래 기다리지 않고 더 빠르게 systemic therapy 로 전환할 수 있게 하는 것입니다.
>
> 다만 공개 데이터셋에서는 실제 RAI response label 이 명확하지 않은 경우도 있기 때문에, 현재 claim 은 "RAI response 를 직접 확정 예측한다" 보다는 "thyroid differentiation 과 iodine-handling 상태를 반영하는 molecular state 를 조기 분류한다" 정도로 잡는 게 안전할 것 같습니다.
>
> 8 개 gene 은 서로 다른 set 이 아니라 명명법 차이로 보입니다. 공식 symbol 은 TG, TPO, TSHR, SLC5A5, DIO1, PAX8, NKX2-1, FOXE1 이고, 설명용으로는 SLC5A5/NIS, NKX2-1/TTF-1, FOXE1/TTF-2 처럼 병기하면 될 것 같습니다.
>
> 가장 중요한 부분은 왜 하필 8 개 gene 인지입니다. 이게 갑자기 나온 것처럼 보이면 리뷰어가 공격할 수 있기 때문에, 문헌 기반으로 RAI uptake · iodine metabolism · thyroid differentiation 관련 후보군을 만들고, 그중에서 biological relevance · cohort 간 재현성 · model parsimony 를 기준으로 최종 8 개를 선정했다는 흐름을 figure 로 보여주는 게 필요할 것 같습니다.
>
> 오늘 확인해야 할 핵심은 분당서울대병원 임상 NGS panel 에 이 8 개 gene 이 얼마나 포함되어 있는지, NGS 시행 갑상선암 환자가 몇 명인지, 그리고 RAI 치료 이력과 response label 을 어느 정도까지 retrospective 하게 정리할 수 있는지입니다.

---

## 9. 오늘 미팅에서 던질 질문 리스트

### 강민수 교수님께

1. 분당서울대병원 NGS panel 에 8 개 gene 이 각각 들어가 있나요?
   - TG · TPO · TSHR · SLC5A5 · DIO1 · PAX8 · NKX2-1 · FOXE1
2. 갑상선암 환자 중 NGS 시행 환자가 몇 명 정도 있나요?
3. 그중 RAI 치료 이력, cumulative RAI dose, post-therapy scan, Tg response, recurrence / persistence 정보를 추적할 수 있나요?
4. RAI 불응성 정의를 어떤 기준으로 잡는 게 제일 임상적으로 맞을까요?
5. 30 명 validation cohort 에서 가장 현실적으로 볼 수 있는 endpoint 는 무엇인가요?
   - RAI uptake · stimulated Tg · post-therapy scan · recurrence · distant metastasis
   - cumulative RAI dose · structural incomplete response
   - leukopenia / marrow suppression

### 병리 협업 관련

1. 수술 후 보관 FFPE block 에서 8 개 marker IHC 가 가능한가요?
2. 8 개 중 실제 IHC antibody 가 안정적인 marker 는 무엇인가요?
3. 병리 scoring 은 H-score, intensity 0/1/2/3, positive cell % 중 어떤 방식이 적절한가요?
4. Cancer cell-specific staining 과 stromal staining 을 분리해서 볼 수 있나요?

---

## 10. IRB 변수 (chart review 항목)

```
[기본 인구학]            나이 · 성별
[병리]                   histology (PTC / FTC / PDTC / ATC) · tumor size · LN met · distant met · TNM stage
[분자]                   BRAF / RAS / TERT 등 driver mutation · NGS panel 결과
[치료]                   수술일 · RAI 시행 여부 · RAI dose · cumulative RAI dose · post-therapy WBS uptake
[생화학]                 stimulated Tg · TgAb
[결과]                   recurrence · persistent disease · progression · systemic therapy · PFS
[독성  (강민수 강조)]    WBC · neutrophil · platelet 변화 · leukopenia · marrow suppression
```

---

## 11. 위험한 표현 vs 안전한 표현

| 위험 (쓰면 안 됨) | 안전 (이렇게 써야 함) |
|---|---|
| "이 8 개 gene 으로 RAI response 를 완벽하게 예측한다" | "This 8-gene thyroid differentiation and iodine-handling signature identifies molecular states associated with RAI avidity and RAI-refractory-like biology." |
| "RAI response predictor" | "RAI-refractory-like molecular state classifier" 또는 "iodine-handling / thyroid differentiation state classifier" |
| "treatment-selection biomarker" | "risk stratification axis" |
| 한국어: "이 8 개로 RAI 반응성을 예측한다" | 한국어: "이 8 개 유전자 패널은 RAI 반응성을 직접 확정 예측한다기보다, RAI 섭취 및 갑상선 분화 상태와 관련된 분자적 상태를 분류하는 도구로 보는 것이 타당하다." |

---

## 12. 현재 deploy / 작성물 inventory

### Manuscript 파일들 (`project/manuscript_v8/`)

| 파일 | 내용 |
|---|---|
| `MANUSCRIPT_REORG_6FIG_PLAN_2026_06_05.md` | 6-figure 재구성 strategic plan (PI 첫 라운드) |
| `EDITORIAL_REVIEW_RUTHLESS_CUT_2026_06_05.md` | NC editor 관점 5-figure cut (더 ruthless) |
| `04_results.md` | 현재 Results draft |
| `05_figure_captions.md` | Cell Press 8-fig v3 captions (CRM-anchor) |
| `05_figure_captions_NC.md` | NC 6-fig captions (NC-reach) |
| `06_discussion.md` | Discussion draft |
| `07_star_methods.md` | STAR Methods |
| `09_reviewer_qa.md` | 예상 reviewer Q&A |
| `10_full_manuscript_compiled.md` | 통합 draft |

### bioRxiv 한글 / 영문 preprint (`project/manuscript_biorxiv_2026_05_20/`)

| 파일 | 페이지 | 크기 |
|---|---|---|
| `paper1_biorxiv.pdf` (영문 master) | 37 | 7.8 MB |
| `paper1_biorxiv_kr.pdf` (한글 v1) | 37 | 8.4 MB |
| `paper1_biorxiv_kr_v2_expanded.pdf` (v2 외부검증 + R17 종합 추가) | 50 | 12.5 MB |
| `paper1_biorxiv_kr_v3_annotated.pdf` (v3 cv2 강조판) | 59 | 18.1 MB |

### 웹 배포 (외부 IP `40.82.129.113`)

| 페이지 | URL | 용도 |
|---|---|---|
| **DM1 한글 스토리 (메인)** | http://40.82.129.113/papers/dm1_story_web/ | **5-7 분 임상의 친화적 한글 다이제스트** (Vite + React + TS) |
| DM1 한글 스토리 (dev) | http://40.82.129.113:5173/ | HMR 가능한 개발 서버 |
| 한글 v3 PDF | http://40.82.129.113/papers/paper1_biorxiv_2026_05_20/paper1_biorxiv_kr_v3_annotated.pdf | 59-페이지 강조판 |
| Paper 1 NC 6-fig dossier | http://40.82.129.113/papers/paper1_nc_6fig_2026_05_27/ | 6 main + 15 ED + EV gallery dossier |

### Figure assets (`project/papers_hub_2026_05_04/assets/paper1_nc/`)

| 카테고리 | 개수 | 내용 |
|---|---|---|
| Main figures (annotated) | 6 | Fig 1-6 cv2 강조판 |
| Extended Data | 15 | ED1 ~ ED15 |
| External Validation gallery | 14 | EV-1 ~ EV-14 + Atlas |
| Panel thumbnails | 42 | Fig 1-6 의 panel-level crops |

---

## 13. GPT 에 붙여넣을 prompt (그대로 복사)

```
You are helping me prepare for a clinical research meeting today (2026-06-25) on a thyroid cancer
RAI-refractoriness 8-gene signature project. I will give you the full background. Read it carefully
and then answer my questions.

═══════════════════════════════════════════════════════════════════
PROJECT CONTEXT
═══════════════════════════════════════════════════════════════════

CORE QUESTION
We want to identify which thyroid cancer patients will fail radioactive iodine (RAI) therapy BEFORE
giving them high-dose RAI repeatedly. The clinical message we are converging on is:

  "Avoid unnecessary high-dose RAI in patients whose tumors cannot take it up — to prevent marrow
   suppression, leukopenia, and long-term toxicity in a disease where most patients live a long time."

This framing (avoid harm) is more clinically defensible than "switch to systemic therapy faster"
(escalation), per our clinical co-investigator Dr Kang Min-su.

THE 8-GENE PANEL
Official HGNC symbols:
  TG, TPO, TSHR, SLC5A5 (NIS), DIO1, PAX8, NKX2-1 (TTF-1), FOXE1 (TTF-2)

Decomposition:
  - 5 effectors (iodine handling / hormone synthesis): TG, TPO, TSHR, SLC5A5, DIO1
  - 3 transcription factors (thyroid lineage): PAX8, NKX2-1, FOXE1

The panel was distilled from a curated 67-gene candidate pool (TIERA67) built on RAI biology priors
(Yoo et al. 2016 PLOS Genet), NOT from outcome-driven gene selection. This is a critical
reverse-causality lock.

DM1 / DM2 STATES
The 8-gene panel partitions thyroid tumors into two molecular states (internal codes DM1, DM2):
  - DM2 = differentiated state, iodine-handling machinery on  → RAI-responsive-like
  - DM1 = dedifferentiated state, iodine-handling machinery off → RAI-refractory-like, aggressive

We may rename DM1/DM2 → "Iodine-handling-low / Iodine-handling-high" for the manuscript to avoid
internal-code feel. DM1/DM2 are NOT existing clinical entities — we created them.

KEY RESULTS (audit-locked numbers, verify with PI before external sharing)

Discovery (TCGA-THCA n = 504 primary PTC):
  - DM1 prevalence 28.4 %, DM2 71.6 % (8-gene KMeans k = 2)
  - Pan-genome top-5000 MAD clustering recovers same partition (ARI = 0.92) → not a panel artifact
  - All 5 canonical drivers (BRAF, TERT, KRAS, NRAS, HRAS) single-feature AUC ~ 0.5 (chance) → driver-orthogonal

Driver-prevalence asymmetry (the dark-matter compartment) :
  - BRAF V600E+ tumors: DM1 prevalence 0.7 %
  - RAS+ tumors:        DM1 prevalence 96.4 %
  - BRAF/RAS-negative:  DM1 prevalence 49.1 % (the patient subgroup where this axis adds maximum substratification value)

Epigenetic mechanism (TCGA HM450 n = 503):
  - TPO promoter β Cohen's d = 2.30, p = 1.9 × 10⁻¹⁸
  - Mean 8-gene β: DM1 0.385 vs DM2 0.253 (+52 %)
  - Per-driver β: BRAF V600E (0.37) ≈ RET fusion (0.39) ≫ RAS (0.27)
    → methylation tracks MAPK activity, NOT driver identity per se

External validation (19 external cohorts spanning 7 modalities, 4 ethnicities):
  - Master cross-cohort forest: 14 entries × 11 cohorts, mean Cohen's d = 2.81, median 2.37,
    every entry d ≥ 1.56
  - Per-gene × cohort × contrast matrix: 80 / 80 cells direction-consistent
  - Lee 2024 (GSE213647, n = 632 Korean PTC FFPE) within-cohort recovery d = 5.93
  - K2 (PRJEB11591, n = 260 Korean FF) within-cohort d = 1.94
  - GPL570 4-cohort microarray (n = 205): 4 / 4 ρ ≤ −0.84
  - Mun 2025 proteogenomic (n = 336): protein-level thyroid_diff d = −1.91, 7/7 sign-consistent
  - Landa 2016 PDTC + ATC: 5 / 8 panel genes overlap their silenced gene list
    (TG · TSHR · TPO · PAX8 · DIO1) — convergent biology, independent design

Clinical association:
  - Pooled OS HR (TCGA + MSK-IMPACT) = 2.53 [1.31, 4.89], I² = 0 % (RETROSPECTIVE)
  - DM1 over-represents the ATA 2015/2025 intermediate-risk tier (RAI decision uncertainty zone)
  - GSE151179 post-RAI refractory tumors transcriptionally align with DM1 state (d ≈ −1.0, p ≈ 10⁻⁴)

Deployability:
  - FFPE (Lee n = 632) vs FF (TCGA n = 504) KS test p = 0.44 → FFPE-compatible

Single-cell:
  - Lu 2023 thyrocyte UMAP (n = 14,624 cells, KRT8∩KRT19∩EPCAM filter): DM score gradient is
    thyrocyte-intrinsic, not stromal/immune confound
  - Pu 2021 paired 6 patients: per-patient r 0.798–0.886, all p < 10⁻¹⁰ Bonferroni

PI / CLINICAL FEEDBACK (Dr Kang Min-su)
1. Reframe primary clinical message → "avoid unnecessary high-dose RAI" (not "escalate to systemic faster")
2. Show 8-gene selection logic explicitly (literature pool → filter → parsimony → comparison
   with 16-gene / random 8-gene / 4-gene subset) — reviewer will attack "why 8?"
3. Check whether the 8 genes are in the hospital's clinical NGS panel:
   - All 8 → direct clinical NGS validation
   - Subset → 4/5/6-gene reduced model
   - Almost none → pivot to RNA panel / targeted assay / IHC
4. IHC validation will be crucial for true clinical translatability; consider diagnostic patent path.

LIMITATIONS WE ARE EXPLICITLY HONEST ABOUT
- Retrospective only; no prospective validation cohort yet (Bundang prospective = 0 % enrolled)
- Pooled OS HR driven heavily by MSK-IMPACT advanced cohort; primary PTC component underpowered alone
- DM1 is a RISK STRATIFICATION axis, NOT a treatment-selection biomarker
- Drug perturbation (DepMap/PRISM) is prioritization-only, not validated mechanism
- Spatial GSE250521 does not survive QC adjustment (caveat-only)
- "Why 8 genes" is the single biggest reviewer attack surface

WORK SPLIT
- Me (Cook): finalize 8-gene symbol/alias, organize cohort analyses, build DM1/DM2 classification,
  compare with existing signatures, build 4/6-gene reduced models, prepare NGS subset model,
  organize figures.
- Dr Kang Min-su: check 8 genes in hospital NGS panel, count thyroid patients with NGS,
  prepare IRB, scope ~30-patient internal validation, define RAI response endpoints.
- Prof Yu: manuscript direction, IRB draft, clinical narrative, authorship, today's 9 AM meeting.

TITLE CANDIDATES (working)
- "A Compact Thyroid Differentiation Signature Predicts Radioiodine Refractoriness in Thyroid
   Cancer Across Multi-Platform Public Cohorts"
- "A Clinically Translatable 8-Gene Signature for Avoiding Ineffective High-Dose Radioiodine
   Therapy in Thyroid Cancer"
- "A driver-orthogonal differentiation axis defines an aggressive subtype of papillary thyroid
   cancer"

═══════════════════════════════════════════════════════════════════
WHAT I NEED FROM YOU
═══════════════════════════════════════════════════════════════════

For today's meeting (and follow-up):

1. Sanity-check my clinical framing — is "avoid unnecessary high-dose RAI" the right primary
   message, or should we present both messages equally?

2. Draft the figure for "why these 8 genes" (literature pool → filter → parsimony → comparison).
   What is the clearest visual?

3. Help me articulate the safest possible main claim that still feels strong for a high-impact
   journal (Nature Communications or equivalent).

4. Anticipate the top 5 reviewer attacks against this manuscript and pre-empt them.

5. What should the IRB document emphasize given that we want to track RAI-induced marrow toxicity
   as a secondary endpoint?

6. Reduced-panel feasibility — if only 4 or 6 of our 8 genes are on the hospital NGS panel, which
   genes should we keep and what is the expected performance drop? (Use the biology: 3 TFs vs
   5 effectors.)

7. Recommend a name change: should we keep DM1/DM2 or switch to
   "Iodine-handling-low / Iodine-handling-high" for the manuscript?

Be rigorous, clinically grounded, and avoid hype. I would rather have a defensible weaker claim
than an attackable stronger one.
```

---

## 14. 한 줄 메모리 (오늘 미팅 1 분 정리)

> 우리는 RAI response 를 과장해서 맞힌다고 주장하기보다, **thyroid differentiation / iodine-handling state 를 robust 하게 분류** 하고,
> 이를 통해 **불필요한 고용량 RAI 를 줄이는 임상 의사결정 도구로 발전** 시키자는 방향.
>
> 8 genes = lens · DM1 = discovery · 메시지 = "harm avoidance, not escalation."
