# Spatial Transcriptomics Brief — Yu 교수님께 (2026-05-07)

대상: Yu Hyeong-won 교수님
보고: Seungho Cook
관련: Paper 2 (HT-overlap PTC) Pillar I + Pillar II 보강
목적: 공간 전사체 데이터 finding 의 한 페이지 보고 + 다음 단계 결정 요청

---

## 1. Why TROP2?

- **이전 (bulk)**: TCGA / GPL570 4-cohort bulk validation 에서 5/5 contrast 중 1만 expected direction → bulk-level enrichment claim 폐기 (2026-05-04)
- **이번 (spatial, n=28 / 3 cohorts)**: GSE250521 PTC + LPTC **8/8 슬라이드** Moran's I 0.27–0.58 (post-fix hex 6-nbr 기준; 정량값은 audit 후 freeze). PT (n=4) + autoimmune-only (n=8) Moran's I < 0.13.
- **해석 가능 범위 (signature-level)**: TROP2 는 scattered tumor population 이 아니라 **공간적으로 조직된 niche** 로 존재. Bulk 가 평균화 시켜 못 본 신호.
- **해석 NOT 가능**: validated biomarker 가 아님. Patient selection / clinical prediction / drug target claim 아직 X.

## 2. Why spatial transcriptomics?

6가지 distinct reason — **bulk 로는 답할 수 없는 질문들**:
1. DM1 region heterogeneity 직접 관찰 (slide 내 영역 간 score 차이)
2. TROP2 scattered vs niche 직접 분리
3. **HT-overlap PTC TLS / B-cell / IGHV niche 직접 확인** (GSE230424 4/4 슬라이드 Moran's I 0.27–0.84) — Paper 2 Pillar II 의 TCGA-bulk-only 약점 보강
4. Driver-orthogonal axis 의 spot-level 분리
5. H&E → DM1 closure battery (Paper 2A NO-GO) 의 self-evidence: LOSO ResNet50 ρ ≈ 0.06
6. Autoimmune-only negative control (GSE248205 8/8 시그널 없음)

## 3. What is new (2026-05-07 추가) — exploratory only

| # | Finding | 수치 | 신뢰도 |
|---|---|---|---|
| New-1 | 분화 program 의 **공간 coherence 붕괴** PT→PTC→LPTC→ATC | 평균 Moran I: 0.46 / 0.57 / 0.61 / **0.12** | exploratory; n=4 ATC caveat |
| New-2 | 24-axis stage gradient | 떨어지는 axis (Thyroid_differen, Thyroid_TF), 오르는 axis (Cell_cycle, T_cell, Senescence_SASP, Macrophage_TAM, Tumor_suppressor) | exploratory |
| New-3 | ATC-2 (GSM7980873) **single-slide 8-axis super-organization** | Cellular_stress 0.84 / Hypoxia 0.80 / EMT 0.73 / Apoptosis_DDR 0.69 / Stemness 0.59 / Cell_cycle 0.79 / Immune_escape 0.82 / Glycolysis 0.85 | **case-like outlier**; n=1 |
| New-4 | LPTC-3 두 번째 outlier 후보 | Hypoxia 0.68 / EMT 0.72 / Stemness 0.59 / Cellular_stress 0.74 | n=1 |

→ 위 4 finding 모두 **signature-level / hypothesis-generating** 으로만 사용. Paper 2 main claim 으로 끌어올리지 않음.

## 4. What is safe to claim

| Allowed | Forbidden |
|---|---|
| spatially organized niche | biomarker |
| spatial coherence | clinical prediction |
| signature-level evidence | causal mechanism |
| exploratory axis | drug target |
| case-like outlier (n=1 caveat) | patient selection |
| hypothesis-generating | validated therapeutic vulnerability |

본 패키지는 **spatial niche 의 존재 + Pillar II 보강** 까지만 주장. ATC-2 같은 outlier 는 case study 형식으로만 supp 진입.

## 5. What should stay supplementary?

- v12 신규 finding 4개 (분화 공간붕괴, ATC-2, LPTC-3, 24-axis gradient) → **모두 supp**
- 8 cancer pathway 공간 maps (S_F81–S_F88) → web-only archive
- 8 deep gene-set maps (S_F73–S_F80) → web-only archive 또는 supp 1 summary
- Pre-fix 의심 figure (S_F31 niche cluster spacing) → web-only / 추후 재생성 옵션

curated 결과:
- **MAIN_FIGURE: 8** (Q1, Q2, Pillar I, Pillar II, closure NO-GO, autoimmune negative control, graphical abstract)
- **SUPPLEMENTARY: 28** (cohort breakdown · QC · v12 supp 4 · methods)
- **WEB_ONLY_ARCHIVE: 64** (redundant · exploratory · pathway atlas)

상세 표는 `project/reports/spatial_figure_curation_2026_05_07.tsv` 와 `PAPER2_SPATIAL_CURATION_DECISION_2026_05_07.md`.

## 6. Manuscript decision request

**Yu 교수님 confirm 필요한 4 항목**:

| # | 결정 | 추천 (Seungho) |
|---|---|---|
| D1 | Paper 2 base venue 선정 | **Scientific Reports** (signature-level + n=4 cohort 적정) |
| D2 | reach venue 도전 여부 | JCI Insight reach (independent validation cohort 1개 확보 시) |
| D3 | v12 finding (분화 공간붕괴 등) 처리 | **Paper 2 supp only** + 별도 preprint 옵션 보유 |
| D4 | Paper 1 outreach 4건 발송 + npj submit click | Marathon 1주차 (5/4–5/10) 마감 전 진행 |

## 7. Reviewer risk caveat (미리 확인 권장)

| Risk | Mitigation |
|---|---|
| n=4 PTC+HT TLS cohort | caption + Limitations 한 줄 |
| n=4 ATC differentiation collapse | supp only + exploratory label |
| n=1 ATC-2 outlier | case-like only; main claim 진입 X |
| v1–v6 Moran's I hex offsets bug | post-fix 재계산 후 main figure freeze; 별도 audit MD 보유 |
| 100-figure firehose | 8 main + 28 supp 만 forward; 나머지 web archive link |

---

## 8. 한 페이지 요약 한 문장

> 28 슬라이드 / 3 cohort spatial transcriptomics 가 (i) TROP2 의 tumor-specific 공간 niche 존재, (ii) HT-overlap PTC 의 TLS/B-cell niche 의 직접 spatial 증거, (iii) H&E → DM1 closure battery NO-GO 의 self-evidence 를 동시에 제공한다. v12 신규 finding (분화 공간붕괴, ATC-2 super-organizer 등) 은 흥미롭지만 **signature-level exploratory** 로만 supp 진입하고 main claim 에 들어가지 않는다.

---

## 9. 다음 step

1. (가장 먼저) Paper 1 outreach 발송 + npj click — Yu 교수님 확인
2. Paper 2 figure curation freeze (8 main + 28 supp)
3. Hex audit 의 NEEDS_REMAKE 7 figure 재생성 (≤ 1시간)
4. v12 supplement vs preprint 결정
5. Independent validation cohort plan (JCI Insight reach gate)

---

## 10. 첨부

- `project/reports/spatial_figure_curation_2026_05_07.tsv` — 100-figure curation table
- `project/reports/PAPER2_SPATIAL_CURATION_DECISION_2026_05_07.md` — 결정 메모
- `project/reports/SPATIAL_HEX_ADJACENCY_AUDIT_2026_05_07.md` — bug audit
- `project/spatial-full/index.html` — 100 figures 시각자료 (web)

— Seungho Cook · 2026-05-07
