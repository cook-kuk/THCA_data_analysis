# 한국 wet-lab outreach 5 lab 한국어 메일 template

**Status**: 사용자 본인 직접 발송 — 각 lab에 customize 후
**작성일**: 2026-04-27 KST
**보낸 사람**: 국승호 (Seungho Cook) · kukshomr@gmail.com

---

## 5 Lab priority + framework

| 우선순위 | Lab | 강점 | 협력 framework | Likelihood |
|---|---|---|---|---|
| 1 | **SNUH 갑상선 외과** (박영주 group 별도) | n 큰 임상 cohort + Korean BRAF rate baseline | RNA-seq retrospective n=50–100 | 중 (Park YJ 본인이 별도 contact 시) |
| 2 | **삼성서울병원 갑상선** | NanoString + clinical trial 인프라 | NanoString 8-probe 임상 검증 | 중-상 |
| 3 | **서울아산병원 두경부암** | 분자병리 강점, ATC/PDTC cohort | 공격적 갑상선암 (ATC/PDTC) RNA-seq | 중 |
| 4 | **가천대 길병원** | TPO/TG IHC 자동화 + 갑상선 sub-specialty | qPCR 8-gene 또는 IHC validation | 낮-중 |
| 5 | **분당서울대** (별도 메일 v2 — Framework C) | RAI uptake clinical scoring | RAI ground truth paired data | 별도 발송 |

---

## ★ COMMON 본문 template (한국어, 각 lab customize)

### Subject

> 한국 갑상선암 8-gene 분자 panel 임상 검증 협력 제안 — npj Precision Oncology 제출 직전 / [LAB 이름] [방법]

### 본문

존경하는 [수신 교수님/연구실장님 성함]께,

안녕하십니까. 서울에서 활동 중인 독립 컴퓨터과학 연구자 **국승호**이며, [유 교수님 성함] 교수님과 6개월 파트너 관계로 갑상선암 transcriptomic 재분석 연구를 진행 중입니다. (회신 시 [유 교수님] 함께 참조 부탁드립니다.)

[LAB 이름]께서 [LAB 강점, 예: NanoString + 갑상선 clinical trial 인프라]에서 갖는 위치 때문에, 본인 분석의 **deployable 8-gene RAI panel** 한국 임상 검증 단계에서 [LAB 이름] 협력이 가장 자연스럽다 판단되어 메일 드립니다.

#### 본인 paper 핵심

- **DM1/DM2 axis** + **8-gene panel** (TPO, TG, SLC5A5/NIS, TSHR, PAX8, NKX2-1, FOXE1, DIO1)
- TCGA-THCA n=513 + 6 GEO + GSE97466 methylation cross-modality 검증
- De-circularized 5-fold CV AUC = **0.962** (95% CI 0.940–0.979) · ΔAUC vs BRAF V600E = **+0.113**
- WHO 2022 morphologic axis (cPTC/IEFVPTC/DHGTC)와 정렬
- Han SC, Park YJ et al. ENM 2023 (BL/RL-PTC) 의 **deployable minimal version**
- npj Precision Oncology 1-2주 내 submit 예정

#### [LAB 이름] 협력 framework

[Lab별로 다음 중 1-2개 선택 + customize]

**Option α (RNA-seq retrospective)**:
- [LAB] 보유 갑상선암 환자 RNA-seq 또는 FFPE 표본 n=50-100
- 본인 8-gene panel + DM1/DM2 axis 적용
- 한국 환자 cohort 검증 + Korean BRAF V600E rate (~62%) 차이 정량
- **본인이 제공**: full reproducible code, dashboard, 분석 자동화, manuscript co-authorship

**Option β (NanoString 8-probe)**:
- 8 probe + 2 housekeeping = 10 probe NanoString panel custom
- [LAB]의 NanoString 인프라 + 본인 분석 protocol
- 임상 deploy-ready validation (n=30-50 환자)
- 별도 paper, co-first authorship 가능

**Option γ (qPCR 8-gene)**:
- TaqMan 8-probe + GAPDH/ACTB
- FFPE 가능 (NanoString보다 cost low)
- pilot n=20-30, scale up to 50-100 후 paper

**Option δ (IHC 4-gene subset)** — TPO/TG/PAX8/NKX2-1만:
- 기존 [LAB] IHC 인프라 활용
- 본인 4-gene reduced panel AUC 0.91 (full 8-gene 0.96 대비 5점 손실)
- 진단 병리 deploy 가능성 가장 높음

[LAB이 어느 framework에 가장 fit인지 1-2개 선택해서 보냄]

#### 본인 + [유 교수님] 제공 가능

- ✓ Full reproducible Python pipeline + interactive dashboard + manuscript draft
- ✓ 데이터 분석 자동화 (의료기관 환경 설정 지원)
- ✓ 통계 컨설팅 + clinical trial protocol 자문
- ✓ Co-authorship — 본인 paper 또는 후속 한국 cohort paper, 협의 가능
- ✓ 발표 자료 (한국어 dashboard + 영문 manuscript draft) 즉시 공유

#### 미팅 제안

본인이 직접 [LAB 위치] 방문 + 30분 brief 진행을 우선 제안 드립니다. [LAB 이름] 다학제 회의 또는 lab meeting 일정 맞춰드릴 수 있습니다. 한국어 dashboard + manuscript v6 영문 draft 모두 미리 보내드리겠습니다.

본인 + [유 교수님] 모두 npj submit 직전이라 timing이 다소 급합니다만, 교수님 일정 우선 맞춰드리겠습니다.

감사합니다.

**국승호 (Seungho Cook)**
독립 연구자 · 파트타임 박사과정
이메일: kukshomr@gmail.com

공동연구자: [유 교수님 성함] · [...]

---

## Lab-specific 추가 단락 (각 메일에 차별화)

### 삼성서울병원

> 삼성서울병원의 NanoString 인프라 + 갑상선 clinical trial 운영 경험은 본인 분석의 **NanoString 8-probe 임상 검증** (Option β)에 가장 fit 한다 판단됩니다. NanoString panel 디자인 protocol + 본인 8-gene logistic regression coefficients + threshold tuning 모두 즉시 공유 가능합니다.

### 아산병원 두경부암

> 아산병원의 ATC/PDTC cohort + 분자병리 강점은 본인 8-gene panel을 **공격적 갑상선암 (ATC/PDTC) 식별** 측면에서 평가하기에 best fit 입니다. GSE76039 (n=37, ATC vs PDTC) 외부 검증 AUC = 0.935; 한국 ATC/PDTC cohort 추가 검증 시 본인 paper의 dedifferentiation framing 강력 reinforce 가능합니다.

### 가천대 길병원

> 길병원의 TPO/TG IHC 자동화 + 갑상선 sub-specialty 강점은 본인 분석의 **4-gene reduced panel IHC validation** (Option δ)에 best fit 입니다. 본인 4-gene 축소 panel (TPO/TG/PAX8/NKX2-1) AUC 0.91 — 진단 병리 deploy 가장 빠른 path가 IHC 입니다.

### SNUH 갑상선 외과

> 박영주 교수 group과의 별도 협력 채널 (별도 메일)이 진행 중이며, SNUH 갑상선 외과 본부와는 임상 cohort RNA-seq 협력 (Option α) 측면에서 별도 협의 부탁드립니다. Han SC 2023 framework 위에서 deployable version을 한국 환자 검증 — 자연스러운 next step 입니다.

### 분당서울대 (별도 메일 v2 — Framework C 우선)

> RAI uptake clinical scoring + 본인 8-gene panel paired data 가 paper의 결정적 game-changer 가 될 수 있습니다. 별도 상세 메일 (`email_bundang_KR_v2.md`) 발송 예정.

---

## Pre-send 체크리스트 (각 lab별)

- [ ] 정확한 수신자 이름 + 이메일 확인
- [ ] [유 교수님] confirm + 공동 발송
- [ ] Lab 강점 단락 customize
- [ ] 첨부 5종 PDF 변환 (manuscript, dashboard, cover letter, RAI proposal, 8-gene SOP)
- [ ] Subject 라인 lab 명 + 방법 specific
- [ ] 발송 timing: 분당서울대 (v2) 우선 발송 → 응답 봐서 SNUH/삼성/아산 → 마지막 길병원
