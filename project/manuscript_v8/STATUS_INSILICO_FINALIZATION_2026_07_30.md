---
title: "Paper 1 (DM1) 현재 상황 정리 — 인실리코 마무리 계획"
date: 2026-07-30
author: Seungho Cook
status: 실험 중단 확정 / 인실리코 완료 후 NC 제출 목표
branch: paper9-perturbation-extension-20260506
---

# Paper 1 (DM1) 현재 상황 정리

## 한 줄 요약

> 실험은 못 한다. 분석은 끝났다. 남은 건 Figure 2개 패널 + 저자 voice section 6개 + 제출 행정이다.

---

## 1. 연구 핵심 한 페이지

### 논문 정체

**제목**: *A thyroid-lineage state predicts radioiodine-refractoriness in BRAF V600E-mutant papillary thyroid cancer*
**타겟 저널**: Nature Communications (v2 draft: NM sub-tier 가능성 언급, 단 SNUBH 내부 검증 전제)
**가장 최신 draft**: `NATURE_COMMUNICATIONS_FULL_DRAFT_v2_2026_07_08.md`
**설계 문서**: `REVIEWER_FIRST_NC_REBUILD_2026_07_08.md` (5-Figure reviewer-first 구조)

### 임상 메시지 (강민수 교수님 확정)

> "RAI가 안 들을 환자에게 불필요한 고용량 RAI를 반복하지 않게 한다."
> (기존 "빨리 항암으로 보낸다" → 이쪽이 더 강함. 갑상선암은 오래 사는 환자 많음)

### 8-Gene Panel

| # | Gene | 임상명 | 기능 |
|---|---|---|---|
| 1 | **TG** | Thyroglobulin | 호르몬 전구체·저장 |
| 2 | **TPO** | Thyroperoxidase | 요오드 산화·유기화 |
| 3 | **TSHR** | TSH receptor | TSH 신호 |
| 4 | **SLC5A5** | **NIS** | 요오드 세포막 운반 (핵심) |
| 5 | **DIO1** | Type 1 deiodinase | 갑상선 호르몬 대사 |
| 6 | **PAX8** | — | 분화 TF |
| 7 | **NKX2-1** | **TTF-1** | 분화 TF |
| 8 | **FOXE1** | **TTF-2** | 분화 TF |

- **DM1** = iodine-handling-low (RAI 내성-like, 예후 나쁨)
- **DM2** = iodine-handling-high (RAI 반응-like, 예후 좋음)

---

## 2. 미팅 경위 및 실험 중단 배경

### 2026-07-13 미팅 (강민수·유형원·이세준·나희영·이인호 교수님)

| 질문 | 결론 |
|---|---|
| Survival 종류? | OS (pooled HR=2.53) + PFI (BRAF+ HR=0.66, p=0.013) 두 가지 모두 |
| RNA → RAI 흡수 직접 연결? | 없음. GSE151179 간접 근거만. 논문에서 prospective validation으로 명시 |
| 내부 검증 가능? | **TSO500**: DNA 플랫폼이라 RNA 발현량 포착 불가 (p=0.30). RNA-seq: 현실적 불가. IHC: 연구용 항체 신뢰성 문제 |
| IHC 예외 가능성? | TG+PAX8+NKX2-1은 이미 임상 진단용 항체 → 후향적 가능성 있었으나... |

### 2026-07-30 현재: 실험 중단 확정

내부 검증(wet lab) 진행 불가. **인실리코만으로 마무리하여 NC 제출.**

---

## 3. 완료된 분석 (Audit-locked 수치)

### 3-1. Discovery (TCGA-THCA n=504)

| 분석 | 수치 | 파일 위치 |
|---|---|---|
| DM1 비율 | 28.4% | `04_results.md` |
| Pan-genome ARI | **0.92** (driver-only = -0.007) | `04_results.md` |
| BRAF mRNA 중립성 | d=-0.044, p=0.57 | `04_results.md` |
| Driver single-feature AUC | 모두 ≈ 0.5 (chance) | `04_results.md` |

### 3-2. Dark Matter Rescue

| 분석 | 수치 |
|---|---|
| BRAF/RAS 음성 회수 | 131/180 (73%) |
| BRAF+ DM1 비율 | 0.7% (거의 없음) → DM1 ≠ BRAF proxy |
| RAS-mutant DM1 비율 | 96.4% |
| BRAF-/RAS- dark matter | DM1 49.1% / DM2 50.9% |

### 3-3. OS/PFI (핵심 임상)

| 분석 | 수치 |
|---|---|
| TCGA OS HR | 2.30 (95% CI 0.77–6.88) |
| MSK OS HR | 2.67 (95% CI 1.17–6.10) |
| **Pooled OS HR** | **2.53 [1.31, 4.89], I²=0%** |
| BRAF+ PFI (n=287) | HR=0.66, p=0.013 |
| DM1 × BRAF interaction | p=0.022, HR_int=0.47 |
| BRAF- PFI | HR=1.19, p=0.52 (효과 없음) |
| RAS+ PFI | HR=1.04, p=0.92 (효과 없음) |

→ **DM1은 BRAF+ 환자에서만 PFI 층화. Predictive biomarker 프레임 가능.**

### 3-4. Kinase Fusion 농축 (TCGA SV)

| 분석 | 수치 |
|---|---|
| DM1 fusion positivity | 63/82 (76.8%) |
| DM1 vs DM2 OR | **7.41** (95% CI 4.38–12.55, p=1.9e-13) |
| RET fusion → DM1 | 27/33 (81.8%) |
| SV missingness | MAR 확인 (chi-square p=0.56) |

### 3-5. Methylation (TCGA HM450, n=503)

| Gene | Cohen's d | p-value |
|---|---|---|
| TPO | **2.30** | 1.9e-18 |
| DIO1 | 1.24 | — |
| TSHR | 1.20 | — |
| PAX8 | 0.97 | — |
| TG | 0.86 | — |
| FOXE1 | 0.84 | — |
| NKX2-1 | 0.63 | — |
| SLC5A5 | 0.22 | 0.42 (예외 — NIS는 프로모터 외 조절) |
| **Mean 8-gene β** | DM1 **0.385** vs DM2 **0.253** (+52%) | — |
| Per-driver-class | BRAF 0.37 ≈ RET 0.39 >> RAS 0.27 | MAPK 비례 |

### 3-6. External Validation (19 cohort)

| 코호트 | 수치 |
|---|---|
| Master forest (14 entries) | mean d=2.81, median=2.37, 모두 ≥1.56 |
| Direction consistency | **80/80 cells** |
| Lee 2024 Korean PTC (n=632) | d=5.93 (FFPE, 외부 재현 강력) |
| Lee 2024 replication (A2-R) | Mann-Whitney p=2.7e-7 (n=370) |
| GPL570 4개 코호트 (n=205) | 4/4 ρ ≤ -0.84 (Western microarray) |
| Mun 2025 proteomics (n=336) | 7/7 단백질 sign-consistent |
| Pu 2021 single-cell | per-patient r=0.798–0.886, Bonferroni p<1e-10 |
| Lu 2023 sc (thyrocyte-gated) | DM gradient 유지 → stromal confound 아님 |
| Landa 2016 PDTC/ATC overlap | 5/8 유전자 (TG·TSHR·TPO·PAX8·DIO1) |
| FFPE vs FF 호환성 | KS p=0.44 |

### 3-7. Post-RAI Refractory 정렬 (GSE151179)

> Cohen's d ≈ -1.0, Mann-Whitney p ≈ 1e-4
> → RAI 내성이 생긴 종양이 transcriptionally DM1 state와 일치. Biological anchor.

### 3-8. IHC 3-plex (In Silico Proxy)

| 조합 | BRAF+ PFI |
|---|---|
| TG + PAX8 + NKX2-1 | **log-rank p=1.7e-4, Cox HR=0.65 [0.47–0.90]** |
| TSO500 3-gene (TSHR+PAX8+NKX2-1) | p=0.30 (DNA 플랫폼 → 실패) |
| Power simulation (n=200, TCGA priors) | Monte Carlo >90% power |

---

## 4. Figure 설계 (5 Main + 6 ED)

### 현재 채택 구조: `REVIEWER_FIRST_NC_REBUILD_2026_07_08.md`

| Figure | 제목 | 핵심 메시지 | Assets 상태 |
|---|---|---|---|
| **Fig 1** | State Discovery | 8-gene axis = driver-orthogonal real state | `F3_driver_neutrality.png`, `F4_pangenome_robustness.png` 존재 |
| **Fig 2** | Biological Meaning | DM1 = thyroid lineage + iodine-handling 소실 | gene expression, GSEA, TDS-16, Landa, Mun 존재 |
| **Fig 3** | Methylation Mechanism | 프로모터 과메틸화 correlate | HM450 heatmap 존재 / **3C scatter 없음** / 3D 재빌드 필요 |
| **Fig 4** | External Validation | Lee+Landa+Mun+4-entry forest | 존재 |
| **Fig 5** | Clinical Relevance | OS forest + PFI interaction + GSE151179 + FFPE | 대부분 존재 |

**v2에서 추가된 figure:**
- **Fig 7**: DM1 head-to-head + interaction + A2-R replication (in `dm1_story_web/public/figures/`)
- **Fig 8**: IHC 3-plex + power simulation (`fig_ihc3_killer.png` 존재)

---

## 5. 남은 작업 체크리스트

### A. 인실리코 — 지금 당장 가능

```
[ ] Fig 3C 생성: HM450 beta ↔ RNA expression scatter
    - TPO / DIO1 / TSHR / TG × 4
    - TCGA paired (n=503) — 데이터 이미 있음
    - "methylation → 발현 억제" 브릿지 패널 (현재 유일하게 없는 패널)

[ ] Fig 3D 재빌드: per-driver-class beta bar
    - BRAF V600E / RET fusion / RAS-mutant / Fusion-negative
    - 이미 계산됨, 시각화만 필요

[ ] Fig 5 ATA intermediate-risk mosaic 재빌드
    - DM1 × ATA 2015 risk tier 2×3 분할표 + 비율 시각화

[ ] 5 main figures PDF 최종 조합 + NC caption 붙이기
    (05_figure_captions_NC.md 기준)
```

### B. 저자 직접 작성 필요 (voice-protected — Claude 불가)

```
[ ] Introduction §1.1 opening hook
    → 임상 역설: 갑상선암 생존율 좋은데 RAI 독성 반복받는 환자군

[ ] Introduction 마지막 문단 (Aim)
    → Factual ingredients: TCGA discovery, DM1/DM2, fusion/methylation, BRAF+ predictive, harm avoidance

[ ] Discussion §3.1 — 기전 해석
    → Landa 2016 JCI 연결, MAPK→methylation correlative (not causal) 경계
    → Krishnamoorthy 2025 misattribution 정정 (→ Landa 2016)

[ ] Limitations §3.4
    → 후향적 공개 코호트 / TCGA event rate 낮음 / MSK = 진행성 농축 /
       methylation = correlative not causal / n=19 fusion-neg DM1 underpowered /
       K2 calibration mismatch / spatial Visium non-confirmatory /
       내부 SNUBH IHC 검증 pending

[ ] Cover letter ¶1
    → 임상 동기 서술

[ ] Reviewer Q9
    → (내용 voice-protected)
```

### C. 제출 행정

```
[ ] Zenodo DOI + 공개 repo URL 확보 (07_star_methods.md 두 곳 placeholder)
[ ] 교신저자 institutional email 확인 (유형원 교수님)
[ ] 공저자 affiliation 최종 확인
[ ] NC submission portal (https://mts.nature.com/) 업로드
[ ] Cover letter + Response letter 초안
```

---

## 6. Reviewer 공격 방어 맵

| Reviewer 공격 | 우리 답변 | 근거 |
|---|---|---|
| "8개가 cluster를 만든 것" | Pan-genome ARI=0.92, driver-only=-0.007 | Fig 1 |
| "BRAF/RAS/fusion 생물학일 뿐" | Driver AUC≈0.5; methylation이 driver 클래스 횡단 | Fig 1+3 |
| "DM1이 뭘 의미하는지 모름" | Thyroid lineage + iodine-handling 소실, TDS-16 동치, Landa overlap 5/8 | Fig 2 |
| "Methylation은 correlative" | 맞음 — "epigenetic correlate" 언어 사용, 인과 주장 안 함 | Fig 3 + Limitations |
| "TCGA cherry-pick" | 80/80 direction, Lee2024 d=5.93, Mun 7/7 단백질 | Fig 4 |
| "임상 주장이 너무 강함" | 후향적 OS HR=2.53, BRAF+ PFI HR=0.66 (retrospective), IHC in silico proxy | Fig 5 + Limitations |
| "분석이 너무 많음" | Spatial/PRISM/pan-cancer/deconv → ED/Supp에만 | ED 6개 |
| "내부 검증 없음" | IHC 3-plex in-silico proxy (임상 항체 기반) + power simulation n=200 90% | Fig 8 |

---

## 7. 원고 파일 지도

```
manuscript_v8/
├── NATURE_COMMUNICATIONS_FULL_DRAFT_v2_2026_07_08.md   ← 최신 draft (v2, 8-figure)
├── NATURE_COMMUNICATIONS_FULL_DRAFT_2026_07_08.md      ← v1 (5-figure)
├── REVIEWER_FIRST_NC_REBUILD_2026_07_08.md             ← Figure 설계 설명서
├── NATURE_CANCER_MANUSCRIPT_ARCHITECTURE_2026_07_08.md ← NC/NM 비교 설계
├── 01_abstract.md                                       ← Abstract (v1 기준, v2로 업데이트 필요)
├── 04_results.md                                        ← Results (v1 기준)
├── 05_figure_captions_NC.md                             ← NC용 caption
├── 06_discussion.md                                     ← Discussion (voice slot 유지)
├── 07_star_methods.md                                   ← STAR methods (cohort 수정 완료)
├── 08_cover_letter.md                                   ← Cover letter (¶1 slot 유지)
├── 09_reviewer_qa.md                                    ← Q&A (Q9 slot 유지)
├── 13_supplementary_tables.md                           ← Supp tables
├── figures/
│   ├── Fig7_dm1_mechanism.pdf/png
│   └── Fig8_epigenetic.pdf/png
└── dm1_story_web/public/figures/
    ├── fig_a2_replication.png
    ├── fig_deep3_redifferentiation.png
    ├── fig_ihc3_killer.png                              ← IHC 3-plex Fig 8
    ├── fig_reduction_grouped_bars.png
    ├── fig_deep1_head_to_head.png
    ├── fig_deep2_prognostic_predictive.png
    ├── fig_deep4_multiomics.png
    └── fig_deep5_snubh_simulation.png
```

---

## 8. 논문 외 프로젝트 현황 (참고)

### CROSS-Neo (네오항원 백신 플랫폼)
- 인실리코 완료: GADGVGKSAL + HMTEVVRHC 두 lead 후보
- MD simulation 10ns 완료 (GADGVGKSAL MD_MODERATE, HMTEVVRHC MD_VERY_STRONG)
- TCR: pMTnet + TEPCAM ensemble AUPRC 0.7841
- 현 상태: wetlab 불가 → 인실리코 priority tool로 포지셔닝

### Paper 2 (H&E → DM1, 이미지 AI)
- UNI LOTO AUC 0.852, GSE250521 ρ=0.392
- 내부 검증 (K2/Bundang H&E) 불가 → NC 아니고 Cell Rep Med 이하
- 현재 보류

---

## 9. 제출까지 예상 타임라인

```
지금 (2026-07-30)
  └─ [인실리코] Fig 3C/3D 생성           ← 1-2일
  └─ [인실리코] Fig 5D ATA mosaic        ← 반나절
  └─ [인실리코] 5 figures PDF 조합       ← 반나절

저자 작업 (voice section)
  └─ Hook + Aim + §3.1 + §3.4 + Cover ¶1 + Q9

행정
  └─ Zenodo + repo URL + email 확인
  └─ NC submission portal 업로드

목표: NC 제출
```

---

*이 파일 최종 업데이트: 2026-07-30*
*관련 핵심 문서: `NATURE_COMMUNICATIONS_FULL_DRAFT_v2_2026_07_08.md`*
