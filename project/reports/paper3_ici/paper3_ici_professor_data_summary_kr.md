# Paper 3 ICI — 교수님 보고용 데이터 요약 (한국어)

**작성일:** 2026-05-06
**작성자:** Seungho Cook
**한 줄 요약:** Paper 3 (ICI 면역유전체 vulnerability) 데이터 등록을 v1 → v2로 확장. 공개 thyroid TME 10건, 임상 evidence 7건, pan-cancer ICI 6건 등록. 분석은 마라톤 (Paper 1/2 출하) 이후로 게이팅 유지. 클레임은 ICI-readiness까지만, response prediction은 데이터상 불가.

---

## 1. 핵심 변경점

| 항목 | v1 (Track A 동결, 2026-05-04) | v2 (이번 sprint, 2026-05-06) |
|---|---|---|
| 마스터 registry | 8개 design 문서 안에 산재 | **단일 TSV** (`paper3_ici_data_registry_v2.tsv`) 23개 entry |
| 임상 evidence map | 없음 | **7개 thyroid ICI trial / paper** TSV 명시 |
| 공개 thyroid TME registry | design 문서 산재 | **단일 TSV** 14개 entry |
| 다운로드 manifest | 없음 | 우선순위 + 접근 가능성 명시 |
| 분석 가능성 매트릭스 | 모듈 단위 | **데이터셋 단위** + 모듈 매핑 |
| 강화의 정량 근거 | "design 동결" | GSE126698 n=365 추가, GSE193581 ATC scRNA 추가, RAI-refractory 축 추가 |

---

## 2. 새로 들어온 데이터 (v2 추가분)

### 2.1 분석에 직접 쓸 수 있는 공개 thyroid TME 데이터

1. **GSE126698 (n=365)** — 36 ATC + 18 PDTC + 132 PTC + 55 FTC + 124 normal. TCGA-THCA 외에 가장 큰 cross-subtype RNA-seq. Module A (ecotype 발견) 의 statistical power 핵심.
2. **GSE184362 (scRNA, 158,577 cells, 11 patients)** — localized + advanced + recurrent LN + RAI-refractory distant met 모두 포함. scRNA atlas의 1차 후보.
3. **GSE193581 (ATC transformation scRNA)** — 10 ATC tumors + 7 PTC + 6 normal + 9 ATC cell lines bulk. ATC scRNA의 희소 자원 — dedifferentiation trajectory 증거.
4. **GSE151179 (RAI-refractory vs RAI-avid PTC)** — aggressive 표현형 proxy. Module A ecotype × RAI status 교차 분석 가능.
5. **GSE126698 / GSE184362 / GSE193581** 는 v1 design 안에 없거나 우선순위가 낮았던 자원들 — **v2 sprint의 가장 큰 강화 포인트**.

### 2.2 임상 evidence (논문 인용용, 분석용 X)

| Trial / Paper | 약물 | 결과 |
|---|---|---|
| KEYNOTE-158 thyroid (Maio 2022) | pembrolizumab | 진행성 갑상선암 cohort, RECIST 기록 |
| KEYNOTE-028 thyroid (Mehnert 2019) | pembrolizumab | PD-L1+ n=22, ORR ~9% |
| NCT03246958 (Lorch 2020 ASCO, Sehgal/Chintakuntlawar) | nivo + ipi | DTC/ATC/MTC 분리 cohort |
| Dierks 2021 Thyroid | lenva + pembro | ATC/PDTC, WES + PD-L1 supplement |
| Spartalizumab ATC (Capdevila 2020 JCO) | spartalizumab | n=42 ATC, ORR ~19%, PD-L1 cutoff 분석 |
| Atezo + matched TT ATC (Cabanillas, NCT03181100) | atezolizumab + 표적 | ATC actionable mutation arm |
| Lenva + pembro RAI-refractory DTC | combo | 다수 trial |

→ Discussion + Introduction hook + Figure 1A 임상 맥락 cite. **분석 입력 아님.**

### 2.3 Pan-cancer ICI (DIAL audit용)

IMvigor210 (UC), Hugo (멜라노마, GSE78220), Riaz (GSE91061), Gide (PRJEB23709), Liu (dbGaP — 신청 필요), Kim GC (PRJEB25780). DIAL audit 의 direction-invariance anchor. **response modeling은 pan-cancer 안에서만 — thyroid에 직접 transfer 금지.**

---

## 3. 지금 분석할 수 있는 것 vs. 게이팅된 것

### 3.1 지금 분석할 수 있는 것 (그러나 마라톤 모드에서는 *실행 X*)

- 공개 GEO 데이터에서 immune module 스코어링 (registry §4 기준 8개 모듈).
- ecotype 후보 클러스터링 (NMF / consensus).
- pan-cancer ICI 코호트에서 DIAL audit.
- scRNA atlas 통합 (scVI/scANVI/Harmony).

→ 모두 **Track B 권한 (사용자 명시 "Track B 시작" + Paper 1 bioRxiv + Paper 2 Task A/B/C closure)** 필요. 이번 sprint는 등록 + 매니페스트 + plan만.

### 3.2 controlled access — 신청 필요

- TCGA-THCA paired BAM (LOHHLA) — dbGaP.
- Liu melanoma (phs000452) — dbGaP.
- Yoo 2019 Korean ATC — EGA 가능성 (verify 필요).

### 3.3 영원히 비공개일 가능성

- thyroid ICI-treated raw RNA-seq — KEYNOTE-158 / -028 / Spartalizumab / Sehgal / Dierks / Cabanillas 모두 raw 비공개. 이것이 **K1 (response predictor 금지) 의 단일 근본 원인.**

---

## 4. 어떤 데이터로 어떤 클레임이 정당화되는가

| 데이터 조합 | 정당화 가능 클레임 | 정당화 불가 클레임 |
|---|---|---|
| TCGA-THCA + GSE76039 + GSE126698 + GSE184362 + GSE193581 | "갑상선암에서 immune ecotype 후보가 dedifferentiation 축을 따라 분포한다" / "BRAF/RAS-negative dark matter 안에서 immune-rich 아군이 식별 가능하다" / **ICI-readiness prioritization hypothesis** | "ICI 치료시 반응 예측" 절대 불가 |
| 위 + IMvigor210 + Hugo + Riaz + Gide | "pan-cancer ICI 코호트에서 DIAL이 통과한 시그니처만 thyroid readiness score에 사용한다" / "TIDE/IMPRES 같은 composite는 thyroid에서 transfer 검증을 통과해야만 사용한다" | "thyroid 환자에서 ICI 반응을 예측한다" 불가 |
| 위 + 임상 evidence 7개 | "thyroid ICI 임상 시도에서 ORR ~9–20% 범위가 보고되었으며, ATC/dedifferentiated 가 PTC 보다 유망하다는 baseline 정황 증거가 존재한다" (인용) | "임상 의사 결정 지원" 불가 |

---

## 5. 강화의 한계 — 솔직한 진술

이번 sprint로 페이퍼는 **데이터의 폭과 깊이**가 모두 강화되었습니다. 그러나:

> **지금 방어할 수 있는 endpoint는 여전히 "ICI-readiness / immunogenomic vulnerability prioritization" 이며, "thyroid ICI response prediction" 은 아닙니다.**

이 한계는 단 하나의 이유로 결정됩니다: **thyroid ICI-treated raw RNA-seq 의 부재.** v2 sprint의 어떤 데이터 추가도 이 한계를 변경하지 않으며, 변경하기 위해서는 **연세 / 서울대 / 분당** 같은 기관에서의 비공개 코호트 협력이 유일한 경로입니다.

---

## 6. 다음 우선순위 (사용자 결정용 옵션)

### A. 마라톤 모드 보존 + Track A 강화 sprint 후속
1. 사전등록(Pre-registration) 초안 — DIAL verdict lock 시점, 가중치 선택지, 1차 vs 민감도 결과 정의
2. K3 (DIAL all-fail) negative-result 시나리오 사전 abstract 작성
3. claim guard forbidden words 리스트 강화

### B. controlled access 신청 (마라톤 비위반)
1. TCGA-THCA dbGaP paired-BAM 신청
2. Liu melanoma phs000452 신청
3. Yoo 2019 EGA 상태 확인 + 신청

### C. 비공개 thyroid ICI 코호트 협상 (Yu 교수 의사결정 필요)
- 연세 / 서울대 / 분당 / 삼성 어디든 thyroid ICI-treated 환자가 있는 기관 1곳 — raw RNA-seq 또는 처리된 expression이라도 가능하다면 K1 해제.
- 이것이 페이퍼 framing 을 "vulnerability" → "validated readiness" 로 격상시키는 *유일한* 경로.

### D. Track B 진입 여부 결정 (Paper 1 출하 ~6/13 이후 자연 트리거)
- G1–G6 6 게이트 closure 시점에 결정.
- 사용자 명시적 "Track B 시작" 명령 필요.

---

Track A 동결 (Paper 1/2 분석) 유지. 마라톤 비위반.
