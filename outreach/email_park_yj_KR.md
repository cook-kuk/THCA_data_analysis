# Park YJ (서울대 교수) 한국어 outreach 메일

**Status**: 사용자 본인이 직접 발송 — 아래 본문 + 첨부 리스트 그대로 활용
**작성일**: 2026-04-27 KST
**보낸 사람**: 국승호 (Seungho Cook) · 독립 연구자 · kukshomr@gmail.com
**받는 사람**: yjparkmd@snu.ac.kr (확인 필요 — SNU 갑상선 외과 또는 분자내분비대사 lab)
**참조**: [유 교수님 성함 + 이메일]

---

## Subject (제목)

> Korean parallel work 협력 제안 — TCGA-THCA transcriptomic axis의 deployable 8-gene panel npj 제출 직전

---

## 본문

존경하는 박영주 교수님께,

안녕하십니까. 저는 서울에서 활동하는 독립 컴퓨터과학 연구자 **국승호**이며, [유 교수님 성함] 교수님과 파트너 관계로 갑상선암 유전체 재분석 연구를 6개월 진행 중입니다. (회신 시 [유 교수님 성함] 교수님 함께 참조 부탁드립니다.)

교수님 연구실 **Han SC, Park YJ et al. ENM 2023** ("Different Molecular Phenotypes of Progression in BRAF- and RAS-Like PTC", PMID 37461149)을 매우 깊게 참고하여 본인 분석을 진행했고, 본인 결과가 교수님 finding을 직접 confirm 하면서 일부 extend 한다 판단되어 짧게 brief 드리고 협력 가능성을 여쭙고자 메일 드립니다.

### 본인 finding 요약 (TCGA-THCA n=513 재분석 + 6 GEO + GSE97466 methylation)

1. **DM1/DM2 axis** — 교수님께서 정의하신 BL-PTC / RL-PTC와 거의 일치하는 transcriptomic axis. WHO 2022 morphologic 분류 (cPTC + 침습형 FVPTC = BRAF-like, IEFVPTC = RAS-like)와도 정렬됩니다.
2. **8-gene RAI panel** — BL/RL distinction의 deployable minimal version (TPO, TG, SLC5A5/NIS, TSHR, PAX8, NKX2-1, FOXE1, DIO1). De-circularized validation 후 진짜 5-fold CV AUC = **0.962** (95% CI 0.940–0.979), ΔAUC vs BRAF V600E baseline = **+0.113**.
3. **Hot/Cold immune composite** — Cohen's d = **+1.683** (p = 4 × 10⁻¹⁸). 교수님께서 ENM 2023에서 보고하신 BL-PTC의 CAF + immune-evasion finding을 confirm 하면서, scRNA scale에서 immune dominance ratio 57:1 로 extend 합니다.
4. **Methylation cross-modality 검증** — GSE97466 (n=141, Illumina 450K) 독립 cohort에서 K=2 clustering 시 모든 aggressive histology (ATC, PDTC, FTC, Hürthle, mFTC)가 met_DM2 군으로 100% 정렬되었습니다. RNA-seq 기반 axis가 methylation 수준에서도 robust 함을 시사합니다.
5. **TERT 4-group 생존 분석** — 36 TCGA-THCA TERT promoter mutation을 cBioPortal `thca_tcga_pub`에서 회수, multivariate logrank p = 3.78 × 10⁻⁵. Honest reframe: univariate Cox HR 6.31 → multivariate (stage + age + sex 보정) HR 1.88 (p = 0.29). TERT⁺을 stage-independent prognostic이 아닌 **Stage III/IV 식별의 분자적 handle**로 framing 합니다.
6. **외부 검증** — GSE76039 (n=37, ATC vs PDTC) AUC = **0.935** (95% CI 0.824–1.000); 7-cohort meta-analysis 진행 중.

### 본인 paper status

- **Target journal**: npj Precision Oncology (1차) — 1-2주 내 submit 예정
- **Manuscript v6** literature 통합 + de-circularization fix 완료
- **Pre-submission self-audit** 통과 (circular validation 발견 → 수정 → 정직 reframe)

### 협력 제안

a) **Han 2023 explicit cite** — 본인 paper Discussion에 다음 단락 명시:

> "Our 8-gene panel is the deployable minimal version of the transcriptomic BL/RL axis described by Han et al. (Endocrinol Metab 2023), with additional Decision Curve Analysis demonstrating clinical net benefit and cross-modality methylation validation in an independent GSE97466 cohort."

b) **Park YJ group의 한국 환자 RNA-seq 데이터** 본인 8-gene panel 검증 (revision round 또는 후속 paper)
   - n = 50–100 정도면 의미 있음
   - Korean BRAF V600E rate (~62%) vs TCGA (~56%) 차이 자체가 paper 한 단락의 가치
   - 본인 분석 + dashboard + reproducible code 100% 공유 가능

c) **Co-authorship 가능성** — 본인 paper의 Discussion 한 단락 자문 또는 후속 Korean cohort paper 공동 1저자 / 시니어 모두 논의 가능

### 미팅 가능 시

- SNU 직접 방문 (제가 가는 것이 맞다 판단됩니다)
- 또는 Zoom 30분
- 본인 분석 dashboard 한국어 버전 + manuscript v6 영문 draft 첨부 드리겠습니다

본인 + [유 교수님] 모두 npj submit 직전이라 timing이 다소 급합니다만, 교수님 일정 우선 맞춰드리겠습니다. 짧은 미팅으로 가능성 여부만이라도 논의 부탁드립니다.

감사합니다.

**국승호 (Seungho Cook)**
독립 연구자 · 파트타임 박사과정
이메일: kukshomr@gmail.com
GitHub: (project URL — submit 후 공개 예정)

공동연구자: [유 교수님 성함]
이메일: [...]

---

## 첨부 (보낼 때 함께)

1. `manuscript_v6.pdf` — 영문 draft (제출 직전)
2. `cover_letter_v5.pdf` — 영문
3. `Korean_dashboard.pdf` 또는 URL — 한국어 review dashboard
4. `interactive_report.html` 링크 — 인터랙티브 종합 리포트
5. `Han_2023_overlap_analysis.pdf` — 본인 finding ↔ Han 2023 직접 비교 (U2B 결과)

## Pre-send 체크리스트

- [ ] 박영주 교수님 정확한 이메일 주소 확인 (yjparkmd@snu.ac.kr 또는 SNU 홈페이지 분자내분비)
- [ ] 유 교수님 confirm + 함께 참조 추가
- [ ] manuscript_v6 read-through 1회
- [ ] Han 2023 PMID 37461149 정확히 cite 했는지 본문 확인
- [ ] U2B 결과 (Han marker overlap) 본문 숫자 fill-in
- [ ] 첨부 5개 모두 PDF 변환 + 압축

## 발송 후

- Park YJ group 회신 24-48h 대기
- 회신 없으면 1주 후 1회 polite follow-up
- 회신 시: 미팅 일정 잡고 + 분당서울대 outreach (U4B)와 시너지 가능성 논의
