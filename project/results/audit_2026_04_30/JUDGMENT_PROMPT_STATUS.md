---
title: 판단(Decision) Prompt 상황 정리 — 2026-04-29 ~ 05-01 audit cycle
author: Seungho Cook
date: 2026-05-01
status: Awaiting user judgment after R1–R7 completion
---

# 판단(Decision) Prompt 상황 정리

> **목적** 4월 29일 유형원 교수 미팅 10 concerns로 시작된 8-gene RAI-responsiveness 논문 (npj Precision Oncology) 자체검수 cycle이 R7까지 끝난 시점. 지금 **사용자(=저자)가 의사결정을 내려줘야 할 7개의 prompt 상황**을 한 곳에 모았습니다. 각 항목마다 **현 상태 / 옵션 / 권장 / 비용·리스크**가 표시되어 있어 그대로 답만 주시면 다음 round로 넘어갑니다.

---

## 0. 현재까지 진행된 audit round 요약 (배경)

| Round | 날짜 | Prompt 수 | 핵심 결과 | Status |
|---|---|---|---|---|
| Initial A–J | 04-29 | 10 | dashboard 67-chart + HIGH_IMPACT_SUMMARY (47KB) | ✅ done |
| **R1** P1–P7 | 04-30 | 7 | TERT⁺ Cox HR=4.33, multi-cohort meta | ✅ done |
| **R2** N1–N7 | 04-30 | 7 | MSK MAF 100% match 복구, autocorrelation OK | ✅ done |
| **R3** F1–F4 | 04-30 | 4 | **DM1 76.8% fusion+ paradigm shift**, GSE286332 framework=Scenario B (별도 paper) | ✅ done |
| **R4** | 04-30 | 4 | Missingness MAR, OR 7–9 robust across scenarios | ✅ done |
| **R5** | 04-30 | 4 | TPO promoter methylation Cohen's d=2.30 | ✅ done |
| **R6** | 05-01 | 4 | meth × expr ρ=−0.68 TPO causal link | ✅ done |
| **R7** | 05-01 | 5 | DM1 fusion+/TERT+ 0 overlap, composite AUC 0.831, **DM1 stemness HIGHER (counter-intuitive)** | ✅ done |

전체 결과물 → `FINAL_COMPREHENSIVE_SUMMARY_v6.md` (14.5KB, submission readiness 220%).

---

## 🟡 판단 #1 — R8를 자동으로 돌릴까, 멈출까?

**상황**
- v6 까지 cycle이 매끄럽게 닫혔음 (9-layer DM1 mechanism, 39 honest disclosures, 12 reviewer Q&As).
- 더 진행하면 **diminishing return 구간** 진입.
- 동시에 R7-3 (arm-level CNV) 404로 deferred 상태 → 이건 새 endpoint 시도하면 회수 가능.

**옵션**
| 옵션 | 추가 분석 | 시간 | 임팩트 |
|---|---|---|---|
| **A. STOP** | R8 안 함, v6 그대로 npj 제출 준비 | 0h | 0 |
| **B. R8-mini** | R7-3 GISTIC2 broad text 재시도 + R8-2 RET 파트너별 phenotype | 2–3h | 중간 (review Q4 보강) |
| **C. R8-full** | mini + R8-3 RAI clinical outcome + R8-4 cross-cohort composite validation | 6–8h | 큼 (translational layer 강화) |

**권장** B (R7-3 미해결 1개 + RET 파트너 1개만). 이미 paper 설득력 충분, 하지만 fusion partner 분리는 reviewer 자주 묻는 항목.

**리스크** C 가면 paper scope creep. 현재 frame "8-gene RAI-responsiveness mini-index"에서 점점 "DM1 multi-omic atlas" 쪽으로 끌려가고 있음.

---

## 🟡 판단 #2 — DM1 stemness 역설 (R7-5) 을 paper에 어떻게 frame할 것인가?

**상황 (R7-5 결과 요약)**
- 통념: DM1 = differentiated (RAI-responsive 잘 됨), DM2 = de-differentiated (de-novo dark matter)
- 측정: stemness 10-gene score → **DM1이 DM2보다 stemness HIGHER** (combined Cohen's d 0.67, p=8.1e-6)
- Drivers: BMI1 d=0.72, MYC d=0.72, CD44 d=0.63, KLF4 d=0.56, SOX2 d=0.30
- Counter: ALDH1A1 d=−0.99 (DM2가 더 높음), OCT4 d=−0.17

**옵션**
| 옵션 | 서술 전략 | 위험도 |
|---|---|---|
| **A. Hide** | Supplementary에만 묻고 main에 안 씀 | 낮음, but reviewer가 ssGSEA 돌려보면 들킴 |
| **B. Reframe** | "DM1=BRAF-driven proliferative stem-like, DM2=ALDH1A1-driven dormant" 로 양자 다 stemness 인정하되 다른 axis 강조 | 중간 (story 살아남) |
| **C. Pivot** | DM1=고스템 differentiated paradox 자체를 finding으로 격상 (Lan et al. 2020 BRAF-driven stemness 인용) | 큼 (impact↑ but Reviewer 1 critique target↑) |

**권장** **B**. ALDH1A1이 DM2 dormant marker, MYC/BMI1이 DM1 active proliferative stem — 두 개 다른 stemness program으로 구분하면 "RAI-responsiveness ≠ stemness 부재"라는 nuanced message 가능. Lan 2020 + Buishand 2018 인용으로 보강.

**Action 필요** 사용자 본인이 임상 의사로서 "DM1 환자가 RAI 잘 듣는데 stemness도 높다"가 본인 보는 환자에서 말이 되는지 직접 판단해 주셔야 함 — 영상의학/병리에서 이 모순이 보이는 경우 있음?

---

## 🟡 판단 #3 — GSE286332 (Korean Hashimoto+THCA) framework 결정 재확인

**상황 (R3-F1, R5-4 결과)**
- TCGA Hashimoto: DM1 risk **감소** (OR 0.31)
- GSE213647 Korean Hashimoto: DM1 risk **증가** (OR 2.4)
- → **방향이 정반대.** Single paper로 묶으면 reviewer Q "왜 cohort 따라 부호 다른가"
- 결정: **Scenario B (별도 paper)** 로 잠정 결정됨 (R3-F1)

**옵션 재확인**
| 옵션 | 처리 | 장점 | 단점 |
|---|---|---|---|
| **A. Scenario A** (현 paper에 통합) | DM1×Hashimoto interaction term + cohort effect 모델 | One paper, more complete | Reviewer가 cohort confounding 강하게 칠 가능성 |
| **B. Scenario B** (별도 paper) | 현 paper는 TCGA/MSK/Lu/Pozdeyev, GSE286332+한국군은 Graves' pivot paper에 합류 | Story clean, npj 빨리 나감 | GSE286332 데이터 일부 손해 |
| **C. Scenario C** (Discussion only) | Main 분석 안 함, Discussion 1 paragraph + Supplementary에 inverse direction 도표만 | 중간 안전 | reviewer가 "왜 분석 안 했냐" 물을 수 있음 |

**권장** **B 유지** — 현재 v17 Graves' pivot이 main paper trajectory로 결정됨 (memory: v17_graves_pivot.md). Hashimoto 데이터를 거기 합치는 게 자연스러움.

**Action 필요** 사용자 한 마디 confirm: "B 유지 / C로 후퇴 / A로 욕심" 중 하나.

---

## 🟡 판단 #4 — composite RNA-only score (R7-2, AUC 0.831) 를 paper Fig으로 올릴까?

**상황**
- 5-feature: 8-gene RAI score, HLA-I, HLA-II, mean 8-gene methylation β, age
- DM1 예측 AUC 0.831, fusion 예측 AUC 0.726
- 모두 **NanoString-feasible** (DNA-based methylation은 EPIC v2 nCounter pilot 필요)
- 현재 Fig 8개 / Table 4개 / Suppl. 12개

**옵션**
| 옵션 | 결과 | 임팩트 |
|---|---|---|
| **A. New Fig 9** | Composite 점수 ROC + nomogram (Korean K2 적용 결과) | translational angle 큼 |
| **B. Suppl. Fig S13** | 데이터만 추가 | 작음 |
| **C. 다음 paper용 hold** | translational paper 따로 | 0 |

**권장** **A**. AUC 0.831 + clinical NanoString feasibility는 reviewer가 정확히 좋아하는 "actionable" 메시지.

**비용** Fig 9 디자인 + caption 작성 약 1.5h. K2 cohort에서 score 계산해 적용한 뒤 KM 그려야 ROC 외에 "outcome impact" 보이려면 추가 0.5h.

---

## 🟡 판단 #5 — submission timing / outreach send

**상황 (memory v17_npj_ship_status.md)**
- npj submission folder 2026-04-27 ship-ready
- **outreach 4 drafts 작성됨 — but NOT 발송**
- 사용자 본인이 send + 클릭 submit 해야 함
- 4월 30일 ~ 5월 1일에 audit cycle 진행하면서 **추가로 보강된 결과물 (R3 fusion paradigm, R5 methylation, R6 causal, R7 composite)** 이 v6 MD에 정리됨

**옵션**
| 옵션 | 행동 | 시점 |
|---|---|---|
| **A. 즉시 제출** | v6 결과 반영 안 하고 ship-ready 그대로 제출 | 5/1 ~ 5/2 |
| **B. R8 mini 후 제출** | 1–2일 추가, R8 + v6 → manuscript 본문 update → 제출 | 5/3 ~ 5/4 |
| **C. R8 full 후 제출** | 일주일 더, full revision | 5/8 ~ 5/10 |
| **D. Graves' pivot 우선** | npj 8-gene 제출 backburner, Graves' paper에 force allocation | 5월 후순 |

**권장** **B**. 4월 30일~5월 1일 audit이 너무 강력해서 (DM1 76.8% fusion+, 9-layer mechanism) 안 반영하고 제출하면 손해. 그러나 D도 사용자가 이미 결정한 경로 (memory). **두 트랙 동시 진행** 가능: (i) npj 8-gene paper는 B로 5/4 제출, (ii) Graves' pivot은 그대로 main trajectory.

**Action 필요** "B + Graves' 병행" / "D (8-gene 완전 backburner)" / "A (그냥 빨리 제출)" 중 결정.

---

## 🟡 판단 #6 — outreach 4 drafts: 누구한테 언제 보낼지

**상황 (memory)**
- 4명 outreach drafts 준비됨 (구체 대상은 메모리에 저장 안 됨 — 사용자가 압니다)
- npj 제출 전/후 어느 시점에 보낼지 결정 안 됨

**옵션**
| 옵션 | 시점 | 효과 |
|---|---|---|
| **A. 제출 전 사전 노티** | 제출 전 1주 | reviewer pool 확보, 그러나 priority claim 위험 |
| **B. 제출 직후** | 제출 후 24h | scoop 위험 ↓, dialog 시작 |
| **C. accept 후** | 발표 직전 | 가장 안전, but 영향력 ↓ |

**권장** **B**. 보통 npj 단계에서는 제출 직후가 표준.

**Action 필요** 4명이 누구이고 어떤 메시지인지 사용자만 앎. 제 쪽에서 보낼 수 없음 — 직접 click send.

---

## 🟡 판단 #7 — Korean K2 mini-index TPM inflation 후속

**상황 (memory v17_korean_k2_calibration.md)**
- 8-gene mini-index TPM이 K2에서 ~10–100× inflated (kallisto idx 문제)
- TCGA-trained absolute-form LogReg는 "wrong-direction" 결과
- **현재 해결책: within-sample-centered profile** (TCGA AUC 0.968)
- 이게 paper에 반영됐는지 ?? (memory 만으로는 불확실, 코드 확인 필요)

**옵션**
| 옵션 | 행동 |
|---|---|
| **A. 현 within-sample 방식 그대로 paper에 명시** | Methods에 "TPM normalization sensitivity" 1 paragraph 추가 |
| **B. 전체 transcriptome index로 K2 재정량** | 추가 2일, 결과 안정성 검증 |
| **C. K2 결과 main에서 빼고 Suppl.로** | 빠른 해결, 그러나 multi-cohort 메시지 약화 |

**권장** **A** + **R8에서 GSE213647로 cross-validate** (R7-4 KS p=0.063 — Korean cohort끼리는 distribution 비슷, 즉 within-sample-z 통하면 둘 다 align 가능)

**Action 필요** B (재정량) 시간 투자 의향.

---

## 한 페이지 의사결정 요약 (사용자가 답해야 할 7개)

```
판단 #1  R8 진행?  →  [ A 멈춤  /  B mini  /  C full ]
판단 #2  Stemness frame?  →  [ A hide  /  B reframe  /  C pivot ]
판단 #3  GSE286332 framework?  →  [ A 통합  /  B 별도 paper(현)  /  C Discussion만 ]
판단 #4  Composite Fig 9?  →  [ A new Fig  /  B suppl  /  C 다음 paper ]
판단 #5  npj 제출 timing?  →  [ A 즉시  /  B R8 mini 후  /  C full  /  D backburner ]
판단 #6  Outreach 4명 시점?  →  [ A 사전  /  B 직후  /  C accept 후 ]
판단 #7  K2 TPM 처리?  →  [ A 현 within-sample 명시  /  B 재정량  /  C suppl 로 강등 ]
```

답 형식 예: "1B / 2B / 3B / 4A / 5B / 6B / 7A" 처럼 한 줄로 주시면 즉시 다음 round로 진입합니다.

---

## 빠른 권장 (제 한 줄짜리)

`1B / 2B / 3B / 4A / 5B (D 병행) / 6B / 7A`

= R8 mini로 RET partner + GISTIC arm CNV 회수 → composite Fig 9 추가 → manuscript R5–R7 결과 반영해 5/4 npj 제출 → 24h 내 outreach send → Graves' pivot paper로 main trajectory 이동.

---

## 참고 위치
- 본 audit 자료: `project/results/audit_2026_04_30/` (round2 ~ round7 폴더)
- 최신 종합 MD: `FINAL_COMPREHENSIVE_SUMMARY_v6.md`
- 미팅 발단 MD: `audit_2026_04_29/HIGH_IMPACT_SUMMARY_for_claude_web.md`
- npj 제출 폴더: `project/submission/npj/`
- 본 판단 MD: `project/results/audit_2026_04_30/JUDGMENT_PROMPT_STATUS.md`
