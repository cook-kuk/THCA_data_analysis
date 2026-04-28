# Park YJ Group Outreach (한국어)

To: 박영주 교수님 (서울대학교, 갑상선 분야)
From: 국승호 (Seungho Cook) <kukshomr@gmail.com>
CC: 유형원 교수님 (Yu Hyeong Won, 분당서울대 외과·내분비외과)
Subject: Korean parallel work 협력 제안 — TCGA-THCA 8-gene RAI decision panel npj 제출 직전, Park YJ group 연구와의 alignment

---

존경하는 박영주 교수님,

안녕하십니까. 저는 서울에서 활동 중인 독립 컴퓨터과학 연구자 **국승호**이며, **유형원 교수님 (분당서울대 외과·내분비외과)**과 공동 1저자 / 공동 교신저자 관계로 갑상선암 transcriptomic 분석 연구를 진행 중입니다. 박사과정은 서울대학교 융합과학기술대학원에서 part-time 으로 진행 중이며, 의료 AI 분야 약 10년 경력 (Portrai 의료영상 CAIO 2.5년, Cornerstone Partners AI Solution Architect)을 보유하고 있습니다.

## 교수님 group 의 연구를 깊게 참고했습니다

본인 paper 의 framework 에 다음 work 들이 직접 연결되어 있어 정중하게 cite 와 협력 가능성을 함께 여쭙고자 합니다:

- **Han SC et al. Endocrinol Metab 2023** — BL-PTC vs RL-PTC progression framework
- **Yoo SK et al. PLOS Genet 2016 (PRJEB11591)** — 본인 paper 에 Korean primary cohort validation 으로 이미 통합
- **Yoo SK et al. Nat Commun 2019 (EGAD00001004845)** — EGA DAR 신청 진행 중 (revision-round 통합 예정)

## 본인 paper 의 핵심 결과 (TCGA-THCA n=513 + PRJEB11591 Korean validation)

| 항목 | 결과 |
|------|------|
| 8-gene RAI panel (NIS, TPO, TG, TSHR, PAX8, NKX2-1, FOXE1, DIO1) | TCGA **5-fold CV AUC** (held-out fold, leak-free) = **0.962** (95% CI 0.940–0.979) |
| vs BRAF V600E baseline | 5-fold CV AUC 0.849 → **ΔAUC = +0.113** (sustained across 4 leak-free cluster definitions) |
| GSE76039 **external held-out** validation (ATC vs PDTC, histology ground truth) | AUC = **0.935** (95% CI 0.824–1.000) |
| **PRJEB11591 Korean validation (Yoo 2016)** | n = 9 pilot subset 처리 완료 → 9/9 DM2 (mean p_DM2 = **0.898**, range 0.619–0.998). ground-truth-labeled AUC 는 전체 262-run 처리 + Yoo SK metadata sharing 후 revision round 측정 예정 |
| Hot/Cold immune composite | Cohen's d = **+1.683**, p = 4.0×10⁻¹⁸ (Han 2023 CAF/immune findings 와 alignment) |
| TERT 4-group survival (Liu-Xing framework) | multivariate logrank p = 3.78×10⁻⁵ |
| Cross-platform transfer success | 97% performance retention (4 transfer barriers) |

## 본인 paper 의 sales pitch

> "BL/RL framework (Han 2023, Park YJ group) 의 **clinical deployable minimal version** + immune integration + leak-free honest validation 입니다. 새 subtype 발견은 아니며, 기존 framework 위에 임상 적용 가능한 8-gene 검사 + 면역 axis 통합 + 정직한 검증 protocol 를 추가한 sub-stratification + integration layer 입니다."

본 paper 는 **npj Precision Oncology** 제출 직전이며, 본인 dashboard (한국어) 에서 모든 분석 결과 + Q&A 확인 가능합니다:
http://40.82.129.113:8765/

## 협력 제안 (3가지)

### (a) Discussion 의 정확한 cite + Korean parallel work positioning

본인 manuscript v6 Discussion 에 다음 단락 추가 검토 중이며, 교수님 의견을 받아 final 결정하겠습니다:

> "Our 8-gene panel is the deployable minimal version of the transcriptomic axis described by Han et al. (Endocrinol Metab 2023) and Yoo et al. (PLoS Genet 2016, Nat Commun 2019), with additional Decision Curve Analysis demonstrating clinical net benefit, PRJEB11591 Korean cross-validation, and a four-tier honest underperformance interpretive framework."

이 framing 이 교수님 group 의 work 를 정확하게 represent 하는지 한 줄 의견 부탁드립니다.

### (b) 추가 Korean cohort validation (revision-round 무기)

만약 교수님 연구실에 본인 8-gene panel 을 직접 검증 가능한 추가 Korean PTC RNA-seq cohort (n ≥ 30) 가 있으시다면:

- **Option 1 — revision round 에 통합**: n = 50–100 정도 cohort 면 본인 paper 의 revision-round 강화에 의미. 본 paper 에 contributor / co-author 로 추가.
- **Option 2 — 별도 후속 paper**: n > 100 또는 더 큰 cohort 면 "Korean-specific 8-gene panel deployment" paper 별도 진행. 박영주 교수님 또는 group 내 senior 가 1저자 또는 senior co-corresponding 가능.

### (c) 30분 미팅 (zoom 또는 SNU 직접 방문)

본인이 SNU 직접 방문 가능합니다. 짧은 미팅으로:
- 본인 진행 결과 공유 (manuscript v6 PDF, 한국어 dashboard)
- Korean parallel work positioning 합의
- 협력 framework 논의 (revision round vs 후속 paper)
- Han SC 박사님 또는 송영신 박사님 cc 또는 introduction 가능 시 그것도 환영

## 첨부

- **manuscript_v6.pdf** — 영문 manuscript draft (npj 제출 직전, ~12,000 words, 20+ figures)
- **한국어 review packet** — 1–2 페이지 결과 요약
- **본인 background** — CV (Portrai CAIO 2.5y, Cornerstone Partners AI Solution Architect, SNU 융합과기원 PhD candidate)
- **본인 분석 dashboard (한국어)**: http://40.82.129.113:8765/easy_explainer.html

## 답변 부탁드립니다

협력 가능성에 대한 한 줄 의견 (긍정/부정/추가 정보 필요) 만 받아도 충분합니다. 본인이 그 다음 step 을 정합니다.

본 paper 는 약 1–2주 안 npj submit 진행 예정이며, 교수님 의견이 본인 paper 의 진짜 strength 가 됩니다. 빠른 시일 내 회신 부탁드립니다.

감사합니다.

**국승호 (Seungho Cook)**
독립 연구자 + 박사과정 학생 (서울대 융합과학기술대학원)
이메일: kukshomr@gmail.com
연락처: [본인 휴대전화 번호]

**공동 교신저자**: **유형원 (Yu Hyeong Won) 교수님**
Department of Surgery, Seoul National University Bundang Hospital, Seongnam, Republic of Korea

---

P.S. **Han SC 박사님** 또는 **송영신 박사님** 께도 cc 또는 introduction 가능하시면 본인이 직접 contact 도 가능합니다.

---

## 본인 액션 (메일 보내기 전)

- [ ] 박영주 교수님 정확한 이메일 확인 (서울대 갑상선 외과 또는 학회 directory)
- [ ] [본인 휴대전화 번호] fill
- [x] PRJEB11591 n=9 pilot 측정값 반영 (mean p_DM2 = 0.898, 9/9 DM2). 전체 262-run + ground-truth AUC 는 revision round.
- [ ] manuscript_v6.pdf + 한국어 review packet 첨부 준비
- [ ] 유형원 교수님께 cc 또는 한 줄 endorsement 받기
- [ ] 발송

## Follow-up

- Day 0: 메일 발송
- Day 7: 답 없으면 정중한 reminder
- Day 14: 답 받으면 → 미팅 일정; 답 없으면 → 다른 contact 시도
