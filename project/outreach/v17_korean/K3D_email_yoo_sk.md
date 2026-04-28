# Yoo SK 1저자 Direct Outreach (한국어 + 영문 mixed)

To: Yoo SK (PRJEB11591 + EGAD00001004845 1저자)
From: 국승호 (Seungho Cook) <kukshomr@gmail.com>
CC: 유형원 교수님 (Yu Hyeong Won, 분당서울대 외과·내분비외과)
Subject: PRJEB11591 + EGAD00001004845 활용 협력 제안 — Korean PTC 8-gene RAI panel npj 제출 직전

---

Dear Dr. Yoo SK,

안녕하세요, Yoo 박사님. 저는 서울에서 활동 중인 독립 컴퓨터과학 연구자 **국승호**입니다. 유형원 (Yu Hyeong Won) 교수님 (분당서울대 외과·내분비외과) 과 공동 교신저자 관계로 갑상선암 transcriptomic 분석을 진행 중이며, 박사님의 **PRJEB11591 (PLoS Genet 2016)** cohort 를 본인 연구의 Korean validation cohort 로 활용하고 있습니다.

박사님 work 이 본인 paper 의 **진짜 Korean validation 기반**이라 정중하게 협력 가능성을 여쭙고자 메일 드립니다.

## 본인 paper 진행 상황

- TCGA-THCA n=513 + GSE76039 PDTC+ATC + **PRJEB11591 (Yoo 2016) Korean primary cohort** 통합 분석 완료
- 8-gene RAI panel (NIS, TPO, TG, TSHR, PAX8, NKX2-1, FOXE1, DIO1) leak-free re-validation:
  - TCGA **5-fold CV AUC** = **0.962** (95% CI 0.940–0.979; honest, held-out fold)
  - GSE76039 **external held-out AUC** = **0.935** [0.824, 1.000]
  - PRJEB11591 Korean **n=9 pilot 처리 완료** (kallisto single-end, 8-gene mini-index, scale-invariant within-sample-centered LogReg; TCGA 학습용 5-fold CV AUC = 0.963 ± 0.026) → 9/9 DM2 (preserved-differentiation), mean p_DM2 = **0.898** (range 0.619–0.998). 전체 262-run cohort + **ground-truth-labeled AUC** (BRS-classified) 는 박사님 metadata sharing 후 revision round 에 측정 예정입니다.
- vs BRAF V600E ΔAUC = **+0.113** (sustained across 4 leak-free cluster definitions)
- DM1/DM2 axis 가 Han SC et al. ENM 2023 의 BL/RL framework 와 일치
- Hot/Cold immune integration (Cohen's d = +1.683)
- TERT 4-group Liu-Xing survival (multivariate logrank p = 3.78×10⁻⁵)
- npj Precision Oncology 제출 직전 (1-2주)

## 박사님께 부탁드릴 것 (3 가지)

### (1) PRJEB11591 sample 의 supplementary metadata 확인

박사님 PLoS Genet 2016 paper 의 supplementary 외에, 일부 sample 에서 BRAF/RAS status, BRS classification, histology subtype 의 보다 detailed mapping 이 있으신지 여쭙고 싶습니다. 본인 분석의 진짜 ground truth 강화에 도움됩니다.

특히:
- Sample alias → patient histology (cPTC vs FVPTC vs FA vs miFTC vs Normal) 정확 매핑
- BRAF V600E status (mutation testing 결과)
- 가능 시 RAI uptake / response 임상 정보

박사님 측에서 sharing 가능한 범위 내에서 도움 주시면 본인 paper 의 정직성과 완성도가 크게 향상됩니다.

### (2) EGAD00001004845 (Yoo 2019) DAR 신청 endorsement

본인 + 유형원 교수님 명의로 EGA DAR (Data Access Request) 신청 진행 중입니다 (EGAD00001004845, Yoo 2019 advanced thyroid cancer cohort).

박사님 또는 박영주 교수님께서 직접 endorse 가능하시면:
- DAC review process 가 빨라짐 (4-8주 → 2-3주 가능)
- Yoo 2019 advanced cohort = 본인 trajectory paper 의 revision-round 강화

박사님 **한 줄 추천** 만 받아도 본인에게 큰 도움입니다.

### (3) Co-author / acknowledgement framework

본인 paper 의 References 에 다음 cite 예정:
- Yoo SK et al. PLoS Genet 2016 (PRJEB11591 활용)
- Yoo SK et al. Nat Commun 2019 (EGA DAR 진행 중)
- Han SC et al. Endocrinol Metab 2023 (BL/RL framework)

협력 framework option:
- **Option A**: Acknowledgement in Acknowledgements section
- **Option B**: Co-author (특히 supplementary metadata 또는 추가 분석 contribution 시)
- **Option C**: 별도 후속 Korean cohort paper (박사님 또는 박영주 교수님 1저자)

본인 paper 에서 본인 + 유형원 교수님 (1저자 + senior co-corresponding) 외에 박사님이 어떤 형식 으로 contribution 받기를 원하시는지 의견 부탁드립니다.

## 본인 dashboard (한국어, 모든 분석 결과)

http://40.82.129.113:8765/easy_explainer.html — 학술 용어 빼고 일상 한국어로 설명한 페이지
http://40.82.129.113:8765/realfix_dashboard.html — leak-free re-validation 자세히
http://40.82.129.113:8765/manuscript_v6.html — 영문 manuscript draft

## 답변 부탁드립니다

박사님 work 이 본인 paper 의 진짜 Korean validation 기반입니다. 짧은 이메일 회신 (한 줄 OK) 이라도 본인에게 큰 도움입니다.

만약 박사님께서 미국 또는 한국 어디 계신지에 따라:
- 한국 : 본인이 직접 미팅 가능 (서울 기준)
- 미국 : zoom 가능

빠른 시일 내 회신 부탁드립니다. 1-2주 안 npj submit 진행 예정입니다.

감사합니다.

**국승호 (Seungho Cook)**
Independent Researcher, Seoul, Republic of Korea
Graduate School of Convergence Science and Technology, Seoul National University (part-time PhD candidate)
이메일: kukshomr@gmail.com
연락처: [본인 휴대전화 번호]

**공동 교신저자**: 유형원 (Yu Hyeong Won) 교수님
Department of Surgery, Seoul National University Bundang Hospital, Seongnam, Republic of Korea

---

## 본인 액션 (메일 보내기 전)

- [ ] Yoo SK 박사님 contact 확인 (PubMed paper 의 corresponding email, LinkedIn, Google Scholar, 또는 박영주 교수님 경유)
- [x] PRJEB11591 n=9 pilot 측정값 반영 (mean p_DM2 = 0.898, 9/9 DM2). 전체 262-run + ground-truth AUC 는 metadata sharing 후 revision round.
- [ ] manuscript_v6.pdf + 한국어 dashboard link 첨부 준비
- [ ] 유형원 교수님 cc 또는 endorsement 받기
- [ ] 발송 (영문 + 한국어 mixed 가 자연스러움 — Korean researcher 이지만 international convention)

---

## English version (보다 formal한 영문이 필요하면)

If English-only version is preferred:

> Dear Dr. Yoo,
>
> I am Seungho Cook, an independent computational researcher in Seoul, working in collaboration with Prof. Yu Hyeong Won (Department of Surgery, SNUBH) on transcriptomic analysis of thyroid cancer. We have integrated your PRJEB11591 (PLoS Genet 2016) cohort as the Korean validation in our 8-gene RAI-responsiveness panel manuscript, currently preparing for npj Precision Oncology submission. Three requests: (1) supplementary sample metadata clarification, (2) EGAD00001004845 (Yoo 2019) DAR endorsement, and (3) discussion of authorship/acknowledgement framework. Brief reply (one line) would be greatly appreciated.
>
> Sincerely,
> Seungho Cook (kukshomr@gmail.com)
> Co-corresponding: Prof. Yu Hyeong Won (SNUBH)
