# Bundang SNUH Korean outreach email — FINAL DRAFT

To: 유형원 교수님 (분당서울대학교병원 외과·내분비외과)
From: 국승호 (Seungho Cook) <kukshomr@gmail.com>
CC: 유 교수님 (Yu Kyungho — 성함 + 이메일 fill)
Subject: 갑상선 유두암 RAI 반응성 8-gene biomarker 한국인 cohort 검증 협력 제안 — npj Precision Oncology 제출 직전

---

존경하는 유형원 교수님,

안녕하십니까. 저는 서울에서 활동 중인 독립 컴퓨터과학 연구자 **국승호**이며, 의료 AI 분야에서 약 10년의 경력을 가지고 있습니다 (Portrai 의료영상 CAIO 2.5년, Cornerstone Partners AI 솔루션 아키텍트). 현재 [Yu Kyungho — 교수님 성함] 교수님과 **공동 1저자/공동 교신저자 관계**로 갑상선암 유전체 분석 연구를 진행 중이며, 박사과정은 [서울대 융합과학기술대학원] 파트타임으로 병행하고 있습니다.

최근 TCGA-THCA 513명 갑상선 유두암 환자의 transcriptomic 재분석을 통해 **임상적으로 즉시 적용 가능한 RAI 반응성 biomarker**를 도출하였고, 본 연구를 **npj Precision Oncology** 에 제출하기 직전입니다. 한국인 cohort 에서의 cross-validation 을 통해 paper 의 임상적 신뢰도를 한 단계 높이고자, 분당서울대병원 갑상선 환자 데이터 활용에 대한 협력을 정중히 제안드립니다.

## 1. 연구 핵심 결과

| 항목 | 결과 |
|------|------|
| **8-gene RAI 반응성 패널** (NIS, TPO, TG, TSHR, PAX8, NKX2-1, FOXE1, DIO1) | 5-fold CV AUC = **0.954** |
| 대비 BRAF V600E 단독 baseline | AUC = 0.822 (**ΔAUC = +0.132**) |
| 외부 4-cohort meta-analysis (n = 290) | Pooled AUC **0.980 (95% CI 0.869–0.997, I² = 0%)** |
| Decision Curve Analysis (Vickers–Elkin) | 0.05–0.95 임계값 91/91 모두에서 8-gene 우월 |
| 9-strata subgroup forest | 7/9 strata AUC ≥ 0.85; **Stage III/IV 군 AUC 0.996** |
| 외부 검증 (GSE76039 PDTC+ATC, n = 37) | AUC 0.974 (correct-direction) |
| Hot/Cold immune composite (DM1 vs DM2) | Cohen's d = +1.683, MW p = 4.0×10⁻¹⁸ |
| TERT 4-그룹 생존 분석 (Liu–Xing 4-genotype) | Multivariate logrank p = 3.78×10⁻⁵; stage 보정 후 HR 1.88 (정직 보고) |

## 2. 분당서울대 협력 제안

### (a) 데이터 modality 우선순위
1. **1순위 — bulk RNA-seq** (full transcriptome, n ≥ 30 이상도 의미 있음)
2. **2순위 — target panel sequencing** (8-gene 중 적어도 6개 이상 포함)
3. **3순위 — DNA panel only** (TERT 프로모터, BRAF V600E mutation 정보만이라도 — TCGA prevalence 7.1% 와 비교 가능)

이상적 sample size 는 n = 50–100 수준이며, matched normal 동반 시 더 강력합니다. 임상 정보로는 RAI 반응성 (가능한 범위), 재발/생존, stage, age, histology subtype 정보가 도움됩니다.

### (b) 협력 framework 제안
- **본 연구 1저자**: 국승호 (분석 lead)
- **공동 교신저자**: [Yu Kyungho 교수님] + 유형원 교수님 (co-corresponding)
- **분당서울대 cohort 분석 결과**:
  - Option A — Paper revision round 에 추가 figure/section 으로 통합. npj 최종 published 시 한국 cohort 분석 contributor 들이 author 로 등재.
  - Option B — 별도 후속 paper 로 진행 (이 경우 분당서울대 PI 가 1저자 또는 senior author).
- **IRB**: 분당서울대 측 lead, 분석은 본인 협력. 데이터는 분당서울대 IRB 가 허용하는 범위에서만 활용.

### (c) Timeline
1. **1차 미팅** (2주 안): zoom 또는 분당서울대 직접 방문. 30분–1시간. 가능성 논의만.
2. **데이터 access 가능성 확인**: 1차 미팅 결과에 따라.
3. **IRB 진행** (필요 시): 분당서울대 lead.
4. **분석 완료**: 데이터 받은 후 3–6개월.

본인이 **분당서울대 직접 방문 가능합니다**. 짧은 미팅으로도 협력 가능성 논의 부탁드립니다.

## 3. 이 연구의 임상적 의미

- 현재 RAI 치료 결정은 **BRAF V600E 단독 + 임상의 판단** 에 의존하며, 정량적 사전 예측 biomarker 가 부재합니다.
- 본 8-gene panel 은 **RNA-seq, microarray, NanoString, qPCR** 어느 platform 에서도 deploy 가능하며, BRAF V600E 단독 대비 ΔAUC = +0.132 의 정량적 우월성을 입증했습니다.
- **두 개의 ongoing thyroid TROP2-ADC trial** (NCT06235216 SETHY, NCT07521670 STRAP) 이 분자 sub-stratification 없이 모집 중이며, 본 panel 을 correlative biomarker overlay 로 제안 가능합니다.

분당서울대의 한국인 환자 cohort 에서 이 panel 의 transferability 가 입증되면, 한국인 갑상선암 환자 임상에 즉시 적용 가능한 도구가 됩니다.

## 4. 첨부 파일

- `manuscript_v4_summary_kr.pdf` — 한국어 1-page 요약 (BOOST sprint 결과 포함)
- `manuscript_v4.pdf` — 영문 manuscript draft (~3,900 words, npj Article format)
- `cover_letter_v3.pdf` — Editor 용 cover letter (참고용)
- `key_findings_1pager_kr.pdf` — 한국어 1-page 핵심 결과 요약
- 인터랙티브 RAI calculator: http://40.82.129.113:8765/rai_calculator.html (browser-based, 서버/설치 불필요)
- 한국어 검토 대시보드: http://40.82.129.113:8765/yu_dashboard_KR.html

(상기 localhost 는 본인 server 이며, 외부 공유용으로는 해당 HTML 파일들을 직접 첨부 또는 이메일에 link 형태로 보내드립니다.)

## 5. 답변 부탁드립니다

협력 가능성에 대한 짧은 한 줄 의견 (긍정/부정/추가 정보 필요) 만 받아도 충분합니다. 본인이 그 다음 step 을 정합니다 (미팅 일정, IRB 절차, etc).

본 연구는 **한국 의료 AI 분야의 정직한 외부 발신**이라는 의미도 있으며, 한국인 cohort 활용은 학술적·임상적으로 모두 의미가 큰 기여라고 생각합니다.

감사합니다.

**국승호 (Seungho Cook)**
독립 연구자 + 박사과정 학생 (서울대 융합과학기술대학원, part-time)
이메일: kukshomr@gmail.com
연락처: [본인 휴대전화 번호]

**공동 교신저자**: [Yu Kyungho — 교수님 성함], [소속]

---

## 본인 액션 (메일 보내기 전 체크)

- [ ] 분당서울대 PI 정확한 성함 + 이메일 확인 (유 교수님께 여쭤보거나 본인 network 활용)
- [ ] [Yu Kyungho — 교수님 성함] 부분 fill in
- [ ] [소속] fill in
- [ ] [본인 휴대전화 번호] fill in
- [ ] 첨부 파일 4개 + 2개 link 준비
- [ ] 유 교수님께 cc 또는 한 줄 endorsement 받기 (가능하면)
- [ ] 보내기 (오늘 또는 내일)

## Follow-up 계획

- Day 0: 메일 보냄
- Day 3: 답 없으면 정중한 reminder
- Day 7: 답 없으면 다른 분당서울대 갑상선과 contact 시도
- Day 14: 답 없으면 그대로 npj 제출 진행
