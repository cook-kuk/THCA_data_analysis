# 분당서울대병원 한국어 outreach 메일 v2 (강화)

**Status**: 사용자 본인 직접 발송 — v1 이후 Yang 2024 + WHO 2022 + Han 2023 통합 강화
**작성일**: 2026-04-27 KST
**보낸 사람**: 국승호 (Seungho Cook) · kukshomr@gmail.com
**받는 사람 후보**:
  - 분당서울대병원 갑상선암센터장 / 외과 (1순위)
  - 분당서울대병원 병리과 갑상선 sub-specialty (2순위)
  - 분당서울대병원 핵의학과 RAI 임상 담당 (3순위 — RAI ground truth 협력 시)

---

## Subject

> 한국 갑상선암 transcriptomic axis + 8-gene RAI panel — 분당서울대 cohort 협력 제안 (npj Precision Oncology 제출 직전)

---

## 본문

존경하는 [수신 교수님 성함]교수님께,

안녕하십니까. 서울에서 활동 중인 독립 컴퓨터과학 연구자 **국승호**이며, [유 교수님 성함] 교수님과 6개월 파트너 관계로 갑상선암 유전체 재분석 연구를 진행하고 있습니다. (회신 시 [유 교수님] 함께 참조 부탁드립니다.)

분당서울대병원이 한국 갑상선암 임상 및 분자 검사 (NIS scintigraphy, RAI uptake clinical scoring, NanoString, qPCR 등)에서 갖는 위치 때문에, 본인 분석의 한국 환자 cohort 검증 단계에서 **분당서울대 협력이 가장 자연스럽다**고 판단되어 메일 드립니다.

### 본인 paper 핵심 finding (TCGA-THCA n=513 + 6 GEO + GSE97466 methylation)

1. **DM1/DM2 분자축** — 갑상선암을 분화축 (de-differentiated DM1 vs differentiated DM2)으로 binary 분류. 박영주 교수 (SNU) Han SC 2023의 BL-PTC/RL-PTC 정의와 일치하며, 본인은 deployable minimal version (8-gene)을 제안합니다.

2. **8-gene RAI panel** — TPO, TG, SLC5A5/NIS, TSHR, PAX8, NKX2-1, FOXE1, DIO1. 모두 canonical thyroid-differentiation marker. **De-circularized 5-fold CV AUC = 0.962** (95% CI 0.940–0.979), ΔAUC vs BRAF V600E = +0.113. NanoString 또는 qPCR 8 probe 1회 측정만으로 임상 deploy 가능합니다.

3. **WHO 2022 morphologic axis와 정렬** — IEFVPTC (encapsulated, RAS-like, indolent) → DM2; cPTC + 침습 FVPTC (BRAF-like, aggressive) → DM1. 본인 axis는 **WHO 2022의 분자적 quantification** 으로 framing 합니다.

4. **TERT promoter 4-group 생존 분석** — 36명 TCGA-THCA TERT⁺ 환자 cBioPortal `thca_tcga_pub`에서 회수, multivariate logrank p = 3.78 × 10⁻⁵. **Yang H et al. ENM 2024** (한국 2,092명 환자 TERT prevalence)와 본인 framing 직접 정렬 — TERT⁺을 stage-independent prognostic이 아닌 **Stage III/IV 식별 분자 handle**로 reframing.

5. **Methylation cross-modality** — GSE97466 (n=141, Illumina 450K) 독립 cohort에서 모든 aggressive histology (ATC, PDTC, FTC, Hürthle, mFTC)가 met_DM2 100% 정렬. RNA-seq 기반 axis가 methylation 수준에서도 robust.

6. **External validation** — GSE76039 (n=37, ATC vs PDTC) AUC = 0.935 (95% CI 0.824–1.000). 7-cohort meta-analysis 진행 중.

### 본인 paper status

- **Target**: npj Precision Oncology (1차) — 1-2주 내 submit
- **Manuscript v6** literature 통합 + de-circularization fix 완료
- **Pre-submission self-audit** 통과 — circular validation 발견 → 수정 → 정직 reframe (학문적 정직성 우선)

### 협력 제안 — 3가지 framework

**Framework A: 분당서울대 갑상선암 환자 RNA-seq 100명 (revision round)**
- 본인 8-gene panel + DM1/DM2 axis 검증
- 한국 환자 BRAF V600E rate (~62% in Korean cohorts vs ~56% TCGA) 차이 정량
- **본인이 제공**: full reproducible code, dashboard, 분석 자동화, manuscript co-authorship
- **분당이 제공**: 한국 환자 RNA-seq + clinical metadata (de-identified)

**Framework B: NanoString/qPCR 8-gene 임상 검증 (별도 paper)**
- 분당서울대 임상 환자 표본 (n = 50–100, FFPE 가능)
- 8 probe만 측정 → P(DM2) 분류 → RAI uptake 임상 결과와 비교
- **이 framework가 가장 임상 implementation에 가깝습니다**
- 별도 paper, co-first authorship 협상 가능

**Framework C: RAI ground truth (★ 가장 가치 큼)**
- Mu Z et al. JCEM 2024 ("RAI avidity 유전체 alteration") 데이터 access 시도 중이나 paywall 가능성
- **분당서울대 핵의학 RAI uptake clinical scoring + 본인 8-gene panel 직접 비교**가 best alternative
- n = 30–50 환자, RAI scintigraphy + 8-gene panel paired data → 진짜 RAI clinical AUC
- 이게 성공하면 본인 paper의 **결정적 game-changer**

### 미팅 제안

본인이 분당서울대 직접 방문 + 30분 brief 진행을 우선 제안 드립니다. 분당병원 갑상선암 다학제 회의 일정에 맞춰 timing 조율 가능합니다. 한국어 dashboard + manuscript v6 영문 draft 모두 미리 보내드리겠습니다.

본인 + [유 교수님] 모두 npj submit 직전이라 timing이 다소 급합니다만, 교수님 일정 우선 맞춰드리겠습니다.

감사합니다.

**국승호 (Seungho Cook)**
독립 연구자 · 파트타임 박사과정
이메일: kukshomr@gmail.com

공동연구자: [유 교수님 성함]
이메일: [...]

---

## 첨부

1. `manuscript_v6.pdf` — 영문 draft
2. `Korean_dashboard_v2.html` 또는 PDF — 한국어 review dashboard
3. `interactive_report.html` URL — 인터랙티브 종합 리포트
4. `RAI_clinical_validation_proposal.pdf` — Framework C 상세 (분량 1-2 page)
5. `cover_letter_v5.pdf`

## v1 vs v2 차이

| 항목 | v1 | v2 (강화) |
|---|---|---|
| Han 2023 cite | ✗ | ✓ (박영주 SNU group 명시) |
| WHO 2022 reframe | ✗ | ✓ (IEFVPTC, DHGTC) |
| Yang 2024 한국 TERT | ✗ | ✓ (2,092 환자 reference) |
| Methylation cross-modality | ✗ | ✓ (GSE97466 결과 추가) |
| De-circularization | 없음 | ✓ (정직 reframe 명시) |
| 협력 framework | 1개 | 3개 (A/B/C) |
| RAI ground truth (Mu 2024) | 없음 | ✓ (Framework C 별도) |

## Pre-send 체크리스트

- [ ] 분당서울대 갑상선암센터 정확한 담당자/이메일 확인
- [ ] [유 교수님] confirm + 공동 발송
- [ ] U2A WHO 2022 결과 본문 fill-in
- [ ] U3B Yang 2024 한국 TERT 정확한 prevalence 수치 확인 (현재 placeholder)
- [ ] Mu 2024 access 결과에 따라 Framework C 강조도 조정
- [ ] 첨부 5종 모두 PDF
