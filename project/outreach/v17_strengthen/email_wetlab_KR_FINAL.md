# Korean wet-lab outreach 5-template pack — FINAL DRAFT

본인 + 유 교수님 명의 한국어 outreach. 8-gene panel qPCR 또는 NanoString 검증 협력.

## 우선순위 (본인 결정)

| # | 후보 lab | 강점 | 비고 |
|---|----------|------|------|
| 1 | **SNUH 갑상선·내분비외과** (서울대학교병원) | 한국 갑상선암 환자수 1위, 외과 lead | 유 교수님 personal connection 우선 |
| 2 | **SNUH 병리학과 (분자병리)** | FFPE 보관량 풍부, qPCR/IHC validation 강점 | |
| 3 | **삼성서울병원 갑상선 lab** | RNA-seq 자체 platform 보유, 분자 진단 lead | |
| 4 | **서울아산병원 두경부암 lab** | ATC/PDTC 케이스 다수, NanoString 인프라 | |
| 5 | **가천길병원 갑상선·내분비외과** | 한국형 cohort registry, 협력 의지 높음 | |

본인 추천: **#1 SNUH 외과** 부터 1-2주 간격으로 1곳씩 발송. 한 번에 5곳에 보내면 conflict 우려.

---

## Template (5개 lab 공용 — [LAB] 부분만 변경)

To: [LAB] [PI 성함] 교수님 <[PI 이메일]>
From: 국승호 (Seungho Cook) <kukshomr@gmail.com>
CC: 유 교수님 (Yu Kyungho)
Subject: 갑상선 유두암 8-gene RAI 반응성 panel — [LAB] qPCR 또는 NanoString 협력 검증 제안

---

존경하는 [PI 성함] 교수님,

안녕하십니까. 저는 독립 컴퓨터과학 연구자 **국승호**이며, [Yu Kyungho 교수님 성함] 교수님과 공동 1저자·교신저자 관계로 갑상선암 유전체 분석 연구를 진행 중입니다.

최근 TCGA-THCA 513명 갑상선 유두암 환자 재분석을 통해 **8-gene RAI 반응성 panel** (NIS/SLC5A5, TPO, TG, TSHR, PAX8, NKX2-1, FOXE1, DIO1) 을 발견하였고, 5-fold CV AUC 0.954, BRAF V600E 단독 (AUC 0.822) 대비 ΔAUC = **+0.132** 의 정량적 우월성을 입증했습니다. 4-cohort meta-analysis 에서 pooled AUC = 0.980 (I² = 0%) 의 robustness 도 확인했으며, Decision Curve Analysis 에서 임상 threshold 0.05–0.95 의 91/91 지점에서 BRAF-only 대비 dominant 한 net benefit 을 보입니다.

본 연구는 **npj Precision Oncology** 제출 직전이며, **한국인 환자 sample 에서의 wet-lab 검증** (qPCR 또는 NanoString) 을 통해 임상 적용 가능성을 입증하고자 [LAB] 와의 협력을 제안드립니다.

## 1. 검증 방식 옵션

### Option 1 — qPCR (8-gene primer set)
- **시료**: PTC FFPE 또는 fresh frozen tissue, n = 20–50
- **소요 비용**: 본인 측 또는 유 교수님 grant 활용 가능 (논의 가능)
- **Timeline**: 1–2개월 (primer 검증 + 실험 + 분석)
- **결과물**: 8-gene log-Ct 값 → 본 panel 적용 → DM1/DM2 cluster 예측

### Option 2 — NanoString (8-gene panel)
- **시료**: FFPE n = 30–100
- **소요 비용**: NanoString 한국 distributor 와 협력 가능 (CRO 옵션 존재)
- **Timeline**: 2–3개월
- **결과물**: 정량적 expression count → 본 panel 적용

### Option 3 — 기존 RNA-seq 데이터 재분석 (있다면, 가장 빠름)
- Wet-lab 작업 없이 in-silico only
- **Timeline**: 2주 이내
- **결과물**: 8-gene 발현값 → DM1/DM2 prediction → 분당서울대/SNUH cohort 와 cross-comparison

## 2. 협력 framework

- **본 연구 1저자**: 국승호 (분석 lead)
- **공동 교신저자**: [Yu Kyungho 교수님] + [LAB PI 성함] (co-corresponding)
- **Wet-lab 검증 part 별도 figure/section** 으로 명시 ("Korean wet-lab validation" 단락 + 별도 figure)
- **저자 형식**:
  - npj revision round 에 통합되면 → 모든 contributor 가 published author
  - 별도 후속 paper 가 더 적합하다면 → [LAB PI] 가 1저자 또는 senior author 인 후속 paper

## 3. 임상적 의의

- **8-gene panel = 임상 즉시 적용 가능**: qPCR 수십 분 안 결과 산출 → RAI 치료 결정에 사용 가능.
- **BRAF V600E 단독 대비 정확**: ΔAUC +0.132 = 임상적으로 의미 있는 차이.
- **한국인 환자 적용 가능성 검증** = critical 한 last-mile validation. 한국 PTC 의 BRAF prevalence (~62%) 가 TCGA (~56%) 와 다르므로, 한국인 cohort 에서의 transferability 입증이 paper 의 마지막 약점 보완.

## 4. 첨부

- `manuscript_v4_summary_kr.pdf` — 한국어 요약
- `manuscript_v4.pdf` — 영문 manuscript draft
- `figures/figure6E_dca.png` — Decision Curve Analysis
- `figures/figure6F_multi_cohort_forest.png` — 4-cohort meta forest
- 인터랙티브 calculator: http://40.82.129.113:8765/rai_calculator.html

## 5. 미팅 가능 시점

본인 직접 [LAB] 방문 가능합니다. 짧은 미팅 (30분–1시간) 으로 가능성 논의 부탁드립니다.

협력 가능성에 대한 한 줄 답변만 받아도 충분합니다. 본인이 그 다음 step 을 정합니다.

감사합니다.

**국승호 (Seungho Cook)**
독립 연구자 + 박사과정 학생 (서울대 융합과학기술대학원)
이메일: kukshomr@gmail.com

**공동 교신저자**: [Yu Kyungho 교수님 성함], [소속]

---

## Lab-별 변형 사항

### #1 SNUH 갑상선·내분비외과
- **추가 문구**: "SNUH 갑상선암 cohort 의 양과 질이 한국 1위인 만큼, 본 panel 의 한국인 transferability 검증에 가장 적합한 site 라 판단합니다."
- **소요 sample**: n = 50–100 가능 (적극 협력 시)

### #2 SNUH 병리학과
- **변형**: Option 1 (qPCR) + 분자병리 IHC (TPO, NIS) 병행 가능 시 강점
- **추가 문구**: "FFPE archive 에서 IHC + qPCR 동시 검증을 통해 platform-portable biomarker 임을 입증하면 임상 적용 가능성이 더 높아집니다."

### #3 삼성서울병원
- **변형**: Option 3 (RNA-seq 재분석) 가능성 우선 검토
- **추가 문구**: "삼성서울병원의 자체 RNA-seq platform 이 있다면, in-silico re-analysis 만으로도 2주 안 결과 가능합니다."

### #4 서울아산병원 두경부암
- **변형**: PDTC/ATC 케이스 활용 강조
- **추가 문구**: "본 panel 의 differentiation continuum 가설은 PTC → PDTC → ATC trajectory 를 통과하므로, 아산의 dedifferentiated case 다수 활용이 strong validation 입니다."

### #5 가천길병원
- **변형**: 한국형 registry 활용 강조
- **추가 문구**: "한국 갑상선암 registry 데이터와의 통합이 한국인 epidemiologic 특성 반영의 가장 효과적인 경로입니다."

---

## 본인 액션

- [ ] 5개 lab 우선순위 결정 (현 추천 SNUH 외과 부터)
- [ ] 첫 1곳 PI 성함/이메일 확인
- [ ] 변형 부분 fill in
- [ ] 유 교수님께 endorsement 받기
- [ ] 발송 (1주에 1곳 cadence 권장)

## Follow-up

- 답 받은 곳: 미팅 일정 → 데이터 협상 → 실제 wet-lab 진행 (3–6개월)
- 답 없는 곳: 1주 후 정중한 reminder, 2주 후 다음 lab 으로 이동
