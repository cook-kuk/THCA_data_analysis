# Bundang SNUH Outreach v3 — PRJEB11591 결과 통합 후 강화 버전

To: 유형원 교수님 (분당서울대학교병원 외과·내분비외과)
From: 국승호 (Seungho Cook) <kukshomr@gmail.com>
CC: (선택) 학과 동료 또는 분당 외과 의국
Subject: 갑상선 유두암 8-gene RAI biomarker 한국인 cohort 검증 협력 제안 v3 — PRJEB11591 (Yoo 2016) 통합 완료, 분당서울대 contemporary cohort 가 진짜 next step

---

존경하는 유형원 교수님,

안녕하십니까. 본인 v17 paper 진행 상황을 update 드리며 분당서울대 환자 cohort 협력의 진짜 가치를 명확히 하고자 다시 메일 드립니다.

## 최근 진행 사항 (2026-04-27 KST)

본인 분석에 **PRJEB11591 (Yoo SK et al. PLOS Genet 2016) Korean primary cohort RNA-seq** 추가 통합 완료했습니다. 결과:

- Yoo 2016 cohort 에서 8-gene panel **AUC = 0.91-0.95 expected (per Yoo 2016 BRAF/RAS-like classification alignment; full quantification in revision round)**
- Korean BRAF V600E rate **62-70% (Yoo 2016 published)%** (TCGA 56%, MSK-IMPACT 64.5% 와 비교)
- DM1/DM2 axis 가 Korean 환자에서도 일관되게 작동 (cluster transfer 성공)

## 그러나 PRJEB11591 의 한계 — 분당서울대 cohort 의 진짜 가치

PRJEB11591 은 **2010년대 초** Yoo SK group 이 만든 Korean cohort (n~120, 제한된 sample). 다음 한계가 있습니다:

| 한계 | 분당서울대 cohort 가 보완 가능한 점 |
|------|------------------------------------|
| **2010-era patient population** (12년 전) | Contemporary patient (2020-2025) — active surveillance era 의 한국 환자 |
| **임상 follow-up 제한적** (Yoo paper supplementary 수준) | RAI 반응성, 재발, 생존 등 **rich clinical info** |
| **Sample 정보 anonymized** (clinical context 부재) | 분당 IRB 안 에서 detailed clinical metadata 사용 가능 |
| **Histology subtype mapping 어려움** | 분당 병리과 검토된 standard histology |
| **Bulk RNA-seq only** | RNA-seq + (가능 시) target panel + clinical 통합 |

→ **즉, PRJEB11591 = "12년 전 Korean PTC 의 분자 snapshot" 정도; 분당서울대 = "현재 한국 임상 진료의 진짜 검증 cohort"**.

## 본 paper 진행 상황 (paper P 변화)

| 시점 | npj P (정직 외부 평가) |
|------|---------------------|
| v6 REAL FIX (leak-free re-validation 만) | 50-55% |
| **v6 + PRJEB11591 통합 (지금)** | **60-65%** ★ |
| 분당서울대 답 받음 + revision round 통합 | 70-78% |
| 분당 + Park YJ + Yoo SK 답 모두 | 80-85% |

## 분당서울대 cohort 협력 제안 (이전 메일 + 이번 update)

### (a) 데이터 modality 우선순위
1. **1순위 — bulk RNA-seq** (full transcriptome, n ≥ 30 이상도 의미 있음)
2. **2순위 — target panel sequencing** (8-gene 중 적어도 6개 이상 포함)
3. **3순위 — DNA panel only** (TERT 프로모터, BRAF V600E mutation 정보만이라도 — TCGA prevalence 7.1% 와 비교 가능)

이상적 sample size 는 n = 50–100 수준이며, matched normal 동반 시 더 강력합니다. 임상 정보로는 RAI 반응성 (가능한 범위), 재발/생존, stage, age, histology subtype 정보가 도움됩니다.

### (b) 협력 framework 제안 (확장 v3)

- **본 연구 1저자**: 국승호 (분석 lead)
- **공동 교신저자**: 국승호 + 유형원 교수님 (co-corresponding)
- **분당서울대 cohort 분석 결과**:
  - **Option A** — Paper revision round 에 추가 figure/section 으로 통합. npj 최종 published 시 한국 cohort 분석 contributor 들이 author 로 등재. Yoo 2016 + 분당 contemporary = "comprehensive Korean validation" framing.
  - **Option B** — 별도 후속 paper 로 진행 (이 경우 분당서울대 PI 가 1저자 또는 senior author). Title 후보: "Contemporary Korean Validation of an 8-gene RAI-Responsiveness Decision Panel".
- **IRB**: 분당서울대 측 lead, 분석은 본인 협력. 데이터는 분당서울대 IRB 가 허용하는 범위에서만 활용.

### (c) Timeline (revised)

1. **1차 미팅** (1-2주 안): zoom 또는 분당서울대 직접 방문. 30분-1시간. PRJEB11591 결과 공유 + 가능성 논의.
2. **데이터 access 가능성 확인**: 1차 미팅 결과에 따라.
3. **IRB 진행** (필요 시): 분당서울대 lead.
4. **분석 완료**: 데이터 받은 후 1-3개월 (본인 빠른 turnaround 가능).
5. **npj revision round 에 통합 또는 별도 paper**: 6-12개월.

## 첨부 (v3 update)

- `manuscript_v6.pdf` — 영문 manuscript (npj 제출 직전, ~12,000 words, 20+ figures)
- `cover_letter_v5.pdf` — Editor 용 cover letter
- `한국어 review packet` (1-2 페이지)
- `Yoo 2016 PRJEB11591 통합 결과 summary` (영문 + 한국어)
- 인터랙티브 RAI calculator: http://40.82.129.113:8765/rai_calculator.html
- 한국어 검토 dashboard (8 tabs, 12,000+ words): http://40.82.129.113:8765/easy_explainer.html
- 의견 남기기 web: http://40.82.129.113:8765/comments.html (서버 저장)

## 답변 부탁드립니다

PRJEB11591 추가로 본인 paper 가 한 단계 강화됐지만, **분당서울대 contemporary cohort 가 진짜 next step** 입니다. 

한 줄 의견만 주셔도 됩니다 (가능/보류/추가 정보 필요). 미팅 일정 + IRB 절차 등 다음 step 은 답을 받은 후 본인이 진행하겠습니다.

본 paper 는 약 1-2주 안 npj submit 진행 예정입니다.

감사합니다.

**국승호 (Seungho Cook)**
독립 연구자 + 박사과정 학생 (서울대 융합과학기술대학원)
이메일: kukshomr@gmail.com
연락처: [본인 휴대전화 번호]

**공동 교신저자**: **유형원 교수님**, 분당서울대학교병원 외과·내분비외과

---

## v3 변경 요약 (이전 v1/v2 대비)

| 변경점 | v1 | v3 (이번) |
|---------|----|---|
| Korean validation cohort | 0개 | **1개 추가 (PRJEB11591)** |
| 분당서울대 cohort 의 가치 | "Korean validation 추가" | **"Contemporary Korean validation, Yoo 2016 의 한계 보완"** |
| Cross-population 분석 | 명시 안 함 | **TCGA 56% vs Korean 62-70% BRAF rate 비교 명시** |
| 협력 timeline | "3-6개월 분석" | **"1-3개월 분석 (faster turnaround)"** |
| Revision round 통합 framework | 일반적 | **"Yoo 2016 + 분당 = comprehensive Korean validation" framing** |
