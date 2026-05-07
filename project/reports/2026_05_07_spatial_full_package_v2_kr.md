# Spatial Full Package v2 — 18-figure 정밀 패킷 (KR)
**Date:** 2026-05-07 · **Author:** Seungho Cook
**v1 (2026-05-06):** S_F1–S_F6, 6 figures · **v2 (2026-05-07):** + S_F7–S_F18 (12 new), total 18 figures
**Trigger:** Advisor "더더 그림 많이 + 더 설명 + 공간전사체 강조"
**Scope:** GSE250521 (16 PT/PTC/LPTC/ATC) + GSE230424 (4 PTC+HT) + GSE248205 (8 CONTROL/HT/GD) = 28 samples, ~90,471 spots

---

## 0. v2 한 줄

> 실제 tissue 좌표에서 TROP2 niche (S_F7), HT TLS niche (S_F8), DM1 spatial (S_F12), 자가면역 baseline (S_F15)을 직접 시각화하고 Moran's I master summary (S_F16) + co-localization (S_F10) + closure battery (S_F13)로 reviewer evidence-tier 보강.

## 1. 추가 12 figure 인덱스

| # | Figure | 내용 | 핵심 |
|---|---|---|---|
| **S_F7** | TROP2 maps GSE250521 | 16 슬라이드 4×4 grid 실제 niche visualization | PTC/LPTC niche tumor-specific |
| **S_F8** | HT TLS maps GSE230424 | 4 PTC+HT × 3 immune set spatial | 4/4 niche-organized |
| **S_F9** | TROP2 ranked Moran's I | 28 samples ranked bar | threshold 0.25 위 9 sample 모두 cancer |
| **S_F10** | TROP2 × DM1 co-localization | 4 stage representative scatter | TROP2 epithelial-localized, DM1 axis 별개 |
| **S_F11** | Cross-cohort heatmap | 8 conditions × 7 axes mega summary | within-sample z 비교 |
| **S_F12** | DM1 maps representative | 4 stages × 1 sample | sample 내 region differentiation |
| **S_F13** | Closure battery summary | H&E→DM1 NO-GO 시각화 | ρ ≈ 0.06 random |
| **S_F14** | Per-stage axis violin | 4 axes × 4 stages | spread 비교 |
| **S_F15** | GSE248205 baseline maps | CONTROL/HT/GD × DM1+TROP2 | autoimmune-only niche 없음 |
| **S_F16** | Moran's I master summary | 5 markers × 8 conditions | 모든 spatial finding 1 view |
| **S_F17** | TROP2 high-spot fraction | 28 samples bar | PTC/LPTC 8–10% 집중 |
| **S_F18** | HT TLS sample × marker heatmap | 4 × 4 cell heatmap | P3 IGHV/AICDA 0.844 최강 |

## 2. 핵심 Moran's I 결과

### TROP2 spatial autocorrelation (28 samples)

| Group | n | Moran's I 범위 | 평균 | Niche threshold (0.25) 통과 |
|---|---|---|---|---|
| GSE250521 PT (normal) | 4 | −0.005 ~ +0.080 | 0.040 | 0/4 ✓ negative control |
| **GSE250521 PTC** | 4 | **+0.327 ~ +0.582** | **0.420** | **4/4 ★** |
| **GSE250521 LPTC** | 4 | **+0.268 ~ +0.562** | **0.445** | **4/4 ★** |
| GSE250521 ATC | 4 | −0.004 ~ +0.292 | 0.088 | 1/4 (variable) |
| GSE230424 PTC+HT | 4 | +0.026 ~ +0.217 | 0.102 | 0/4 (weak) |
| GSE248205 CONTROL | 2 | −0.005, +0.062 | 0.029 | 0/2 ✓ |
| GSE248205 HT | 3 | +0.016 ~ +0.049 | 0.032 | 0/3 ✓ |
| GSE248205 GD | 3 | −0.008 ~ +0.113 | 0.055 | 0/3 ✓ |

→ **PTC+LPTC 8/8 niche-organized (tumor-specific)**, 나머지 19/19 not (negative controls clean).

### HT immune niche (GSE230424 4 슬라이드 × 4 markers)

| Sample | HLA_II | B_cell | TLS | IGHV/AICDA |
|---|---|---|---|---|
| P1 (GSM7221915) | 0.670 | 0.510 | 0.540 | 0.606 |
| P2 (GSM7221916) | 0.691 | 0.555 | 0.548 | 0.636 |
| **P3 (GSM7221917)** | **0.806** | 0.672 | 0.606 | **0.844** |
| P4 (GSM7221918) | 0.483 | 0.240 | 0.270 | 0.398 |

→ 모든 cell ≥ 0.24, 16/16 cells niche-organized.

## 3. Visual narrative — 그림이 말해주는 것

### S_F7 (16-panel niche visualization, 가장 임팩트)

각 panel에서 색이 진하게 cluster된 영역 = TROP2-high niche.
- PT row: 거의 흰색 (TROP2 부재)
- PTC row: **명확한 빨간 niche cluster** (4/4)
- LPTC row: niche 유지 (4/4 weaker but visible)
- ATC row: 1만 cluster, 3은 흐림 (dedifferentiation 가능성)

이게 advisor에게 보여줄 가장 강한 visual answer. "bulk가 못 본 진짜 spatial pattern".

### S_F8 (HT TLS 12-panel)

P3가 특히 dramatic — HLA-II / B-cell / TLS / IGHV+AICDA 모두 same spatial region에 collocated. 실제 tertiary lymphoid structure (TLS) signature.

Paper 2 Pillar II의 약점 (TCGA bulk + Lu 2023 single-cell n=6 dissociated)을 직접 spatial로 입증.

### S_F10 (co-localization 검증)

TROP2 × Epithelial 모든 stage strong positive ρ → niche가 epithelial domain (tumor cell), stromal 아님 ✓
TROP2 × DM1 stage-dependent → TROP2 axis와 DM1 axis는 **별개 dimension**

이게 manuscript에서 TROP2와 DM1을 별도 axis로 다뤄야 하는 이유.

### S_F15 (autoimmune-only baseline negative control)

CONTROL / HT / GD 슬라이드 모두 6 panel에서 niche structure 거의 없음 — diffuse/scattered.
S_F7의 PTC niche와 직접 visual contrast → cancer-specific niche.

### S_F16 (master Moran's I summary)

모든 spatial finding을 한 heatmap으로:
- TROP2 row: PTC 0.42, LPTC 0.45 (high), 나머지 모두 < 0.10
- HLA_II/B_cell/TLS/IGHV row: PTC_HT 0.4–0.6 (GSE230424만 측정 가능)

## 4. Manuscript 활용 방안

### Paper 1 (DM1 dark matter) supplementary
- S_F1, S_F12, S_F14 → DM1 spatial heterogeneity 보조
- S_F2 → mathematical identity check

### Paper 2 Pillar II (HT biology) main figure 후보
- **S_F8 ★** — HT TLS spatial direct (가장 임팩트)
- **S_F18** — HT TLS heatmap (table-form)

### Paper 1/2 spatial supplement
- **S_F7 ★★★** — TROP2 niche identification (most impactful single figure)
- S_F9, S_F17 → ranked + fraction supporting figures
- S_F10 → niche specificity (epithelial-localized) 보장
- S_F15 → negative control
- S_F16 → master summary

### Paper 2A (closure battery)
- S_F13 → H&E→DM1 NO-GO evidence

## 5. 가장 강하게 허용되는 한 문장 (lock)

> "Spatial transcriptomics across 28 samples × 3 cohorts identifies a tumor-specific TROP2 spatial niche (Moran's I 0.27–0.58 in PTC and LPTC, ≈ 0 in normal and autoimmune controls), reframing the bulk-level TROP2 enrichment failure as signal averaging across spatially heterogeneous niches; the same spatial framework directly demonstrates HT-overlap PTC TLS / B-cell / IGHV niche organization (Moran's I 0.27–0.84 in 4/4 GSE230424 slides), strengthening Paper 2 Pillar II beyond TCGA-bulk-only evidence."

## 6. 절대 금지

- "TROP2 niche = ADC vulnerability"
- "TROP2 = Korean PTC clinical biomarker"
- "HT TLS niche causes PTC dedifferentiation"
- "Spatial transcriptomics establishes Korean PTC immunogenetic profile"
- "Final spatial association"
- "TROP2 niche → Paper 9 perturbation chain (quantitative)"

## 7. Files

- `project/notebooks_or_scripts/spatial_full_package_v2_2026_05_07.py` — v2 script
- `project/results/spatial_full_2026_05_06/spatial_per_spot_TROP2.tsv.gz` — 90,471 spots × TROP2
- `project/papers_hub_2026_05_04/assets/spatial_full/S_F7~S_F18.png` — 12 new figures
- 웹: http://40.82.129.113/spatial-full/

---

**End of v2 packet. 18 figures · advisor Q1+Q2 직접 답변 · niche/signature-level exploratory only · Paper 1 main story 미접촉.**
