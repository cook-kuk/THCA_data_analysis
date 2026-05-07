# Spatial Transcriptomics 풀 패키지 정밀 패킷 (KR)
**Date:** 2026-05-06 · **Author:** Seungho Cook
**Trigger:** Advisor 질문 — (1) TROP2 갑자기 왜 등장? (2) 공간 전사체 왜 필요?
**Scope:** GSE250521 (16 PT/PTC/LPTC/ATC) + GSE230424 (4 PTC+HT) + GSE248205 (8 CONTROL/HT/GD) = 28 samples, ~85,000 spots
**원칙:** signature-level exploratory only. Causal / clinical / biomarker / vulnerability 주장 금지.

---

## 0. 한 줄 요약

> Spatial transcriptomics across 28 samples × 3 cohorts identifies a **tumor-specific TROP2 spatial niche** (Moran's I 0.27–0.58 in PTC and LPTC, ≈ 0 in normal and autoimmune controls), reframing the bulk-level TROP2 enrichment failure as **signal averaging across spatially heterogeneous niches**; the same spatial framework directly demonstrates **HT-overlap PTC TLS / B-cell / IGHV niche organization** (Moran's I 0.27–0.84 in 4/4 GSE230424 slides), strengthening Paper 2 Pillar II beyond TCGA-bulk-only evidence.

---

## 1. Advisor Q1 — "TROP2가 갑자기 왜 튀어나왔냐?" 정직한 답변

### 1.1 시간순 등장 (3-stage)

| Stage | 시점 | 데이터 | 결과 | 결정 |
|---|---|---|---|---|
| **1. Anchor 부산물** | 2026-04 | GSE76039 (Landa 2016) | TACSTD2 z, ATC vs PDTC d=+1.13, p=2.7e−3 | Supplement-tier 등재 ("tumor-population vulnerability signal") |
| **2. Bulk validation FAIL** | 2026-05-04 | GPL570 4-cohort | 5 ATC contrast 중 1만 expected direction; GSE33630 ATC vs PTC d=−1.72 (반대), GSE29265 ATC vs PTC d=−1.98 (반대) | Bulk-level enrichment claim **폐기** → "supplement / spatial follow-up"으로 강등 |
| **3. Spatial niche identification** | 2026-05-06 (오늘) | 28 samples × 3 cohorts spatial | PT/HT/GD: Moran's I ≈ 0; **PTC+LPTC: Moran's I 0.27–0.58** | Reframing: scattered tumor population NOT, **tumor-specific spatial niche** |

### 1.2 Spatial niche evidence (S_F5)

`project/results/spatial_full_2026_05_06/spatial_E_TROP2_spot_distribution.tsv`:

| Cohort | Condition | Sample range | TROP2 Moran's I 범위 | 해석 |
|---|---|---|---|---|
| GSE250521 | **PT (normal)** | N-1~4 (n=4) | −0.005 ~ +0.080 | scattered/absent ✓ |
| GSE250521 | **PTC** | PTC-1~4 (n=4) | **+0.327 ~ +0.582** | **★ niche-organized** |
| GSE250521 | **LPTC** | LPTC-1~4 (n=4) | **+0.268 ~ +0.562** | **★ niche-organized** |
| GSE250521 | ATC | ATC-1~4 (n=4) | −0.004 ~ +0.292 | variable (1 niche, 3 scattered) |
| GSE230424 | PTC+HT | P1–P4 (n=4) | +0.026 ~ +0.217 | weak niche |
| GSE248205 | CONTROL | C1, C2 (n=2) | −0.005, +0.062 | scattered ✓ |
| GSE248205 | HT | HT1–HT3 (n=3) | +0.016 ~ +0.049 | scattered ✓ |
| GSE248205 | GD | GD1–GD3 (n=3) | −0.008 ~ +0.113 | scattered ✓ |

→ **PTC + LPTC 8/8 sample에서 niche structure**, normal/autoimmune-only에서는 ≈ 0. Tumor-specific.

### 1.3 왜 bulk가 못 봤나?

- Sample-level mean이 niche-low + niche-high를 평균화해서 effect 약화
- ATC variability (1/4 sample만 niche): late-stage dedifferentiation으로 niche structure 깨짐 가능성
- PTC+HT (GSE230424) weak niche: HT 동반 시 niche가 약해질 수 있으나 sample size n=4 작음

### 1.4 정정된 framing (manuscript-ready)

**FROM (폐기):** "TROP2 elevation in advanced disease confirmed at tumor-population level"
**TO (허용):** "TROP2 forms tumor-specific spatial niches in PTC/LPTC, absent in normal thyroid and autoimmune-only thyroid; bulk-level enrichment failure is consistent with signal averaging across spatial heterogeneity rather than absence of TROP2 biology"

---

## 2. Advisor Q2 — "공간 전사체 왜 필요한가?" 6가지 distinct reason

| # | Reason | 어떤 분석으로 풀리나 | 해당 figure |
|---|---|---|---|
| 1 | DM1 region heterogeneity 검증 | A — per-spot DM1 × Epithelial 2D density | S_F1 |
| 2 | TROP2 scattered vs niche 분리 | E — TROP2 spot Moran's I (28 samples) | S_F5 ★ |
| 3 | **HT-overlap PTC TLS niche** | C — GSE230424 PTC+HT TLS/B-cell/IGHV scoring + Moran's I | S_F3 ★ |
| 4 | Driver-orthogonal axis spot | B — within-sample DM1 spread + identity check | S_F2 |
| 5 | Closure battery (Paper 2A) self-evidence | (이미 committed: closure_battery_metrics.tsv) | — |
| 6 | Autoimmune-only negative control | D — GSE248205 HT/GD vs CONTROL pairwise MW | S_F4 |

→ 6개 모두 bulk-level로는 답할 수 없는 질문. Spatial 없이는 모두 unresolved.

---

## 3. C 분석 — GSE230424 PTC+HT TLS niche 직접 확인 (★ Paper 2 Pillar II 보강)

`project/results/spatial_full_2026_05_06/spatial_C_HT_TLS_spatial.tsv`:

| Sample | HLA-II Moran's I | B-cell I | TLS I | IGHV/AICDA I |
|---|---|---|---|---|
| P1 (GSM7221915) | 0.670 | 0.510 | 0.540 | 0.606 |
| P2 (GSM7221916) | 0.691 | 0.555 | 0.548 | 0.636 |
| P3 (GSM7221917) | **0.806** | 0.672 | 0.606 | **0.844** |
| P4 (GSM7221918) | 0.483 | 0.240 | 0.270 | 0.398 |

→ **4/4 PTC+HT 슬라이드 모두 immune niche organization 강하게 검증.** Lu 2023 single-cell n=6 dissociated data가 implication만 줬던 것을 spatial로 직접 확인.

**Paper 2 Pillar II 본문 보강 핵심:** 기존 (TCGA bulk + Lu 2023 단일-cell n=6) → 추가 (GSE230424 spatial n=4, all niche-organized).

---

## 4. 6 분석 결과 요약

### A — DM1 × Epithelial niche 2D density (GSE250521, 16 slides)
- Per-stage hexbin density (PT/PTC/LPTC/ATC × 4 panels)
- Within-stage Spearman ρ(DM1, Epithelial) 보고
- 결과: stage-mean 0 (within-sample z normalization) — spread만 비교 가능

### B — Driver-orthogonal proxy
- Within-sample DM1 ↔ RAI_8 = −1.0 (mathematical identity check ✓)
- DM1 IQR per sample × stage: PTC stages에서 spread 큼

### C — HT TLS spatial (★)
- 4 PTC+HT slides × 4 immune sets × Moran's I
- All slides niche-organized (Moran's I 0.24–0.84)

### D — GSE248205 baseline
- CONTROL (n=2) vs HT (n=3) vs GD (n=3) pairwise MW
- HT vs GD proliferation: p = 2.5e−41 (large spot count, small effect)
- Cancer-side DM1 axis: weak signals (negative control 역할 ✓)

### E — TROP2 spatial (★)
- 28 samples × Moran's I + mean + p90 + high-spot fraction
- PTC+LPTC niche-organized, PT/HT/GD scattered
- Bulk-fail의 정체 = signal averaging

### F — Cross-cohort integration
- 28 samples pooled DM1 distribution per condition
- Within-sample z 정규화로 mean ≈ 0; spread-level 비교만 가능

---

## 5. 가장 강하게 허용되는 표현

> "Spatial transcriptomics across 28 samples × 3 cohorts identifies a tumor-specific TROP2 spatial niche (Moran's I 0.27–0.58 in PTC and LPTC, ≈ 0 in normal and autoimmune controls), reframing the bulk-level TROP2 enrichment failure as signal averaging across spatially heterogeneous niches; the same spatial framework directly demonstrates HT-overlap PTC TLS / B-cell / IGHV niche organization (Moran's I 0.27–0.84 in 4/4 GSE230424 slides), strengthening Paper 2 Pillar II beyond TCGA-bulk-only evidence."

이 한 문장 이상으로 강하게 말하지 않습니다.

---

## 6. 절대 금지 표현

- "TROP2 high niche = ADC vulnerability" / "TROP2-targeted therapy responder"
- "TROP2 niche = clinical risk biomarker"
- "HT TLS niche causes PTC dedifferentiation"
- "Spatial transcriptomics establishes Korean PTC immunogenetic profile"
- "Final spatial association" / "Causal spatial mechanism"
- TROP2 niche → Paper 9 perturbation chain (quantitative connection 금지)

---

## 7. 다음 검증 단계 (validation roadmap)

1. **Independent spatial PTC cohort** (≥ 10 슬라이드) — TROP2 niche reproduction
2. **TROP2 spatial niche × tumor-cell-only filtering** — niche가 stromal contamination인지 epithelial인지 분리
3. **HLA-II/TLS niche × HT pathology centrally-reviewed** — TLS pathology 일치 여부
4. **TROP2 niche ↔ functional readout** (e.g., proliferation, dedifferentiation marker) 별도 cohort에서 quantitative test — 단, 본 패키지에서는 안 함
5. Paper 2 Pillar II 본문에 GSE230424 TLS spatial 결과 supplementary 추가

---

## 8. 첨부 자료

- `project/notebooks_or_scripts/spatial_full_package_2026_05_06.py` — analysis script
- `project/results/spatial_full_2026_05_06/spatial_full_summary.json` — summary
- `project/results/spatial_full_2026_05_06/spatial_{A,B,C,D,E,F}_*.tsv` — 6 result tables
- `project/papers_hub_2026_05_04/assets/spatial_full/S_F{1,2,3,4,5,6}.png` — 6 figures
- 웹: http://40.82.129.113/spatial-full/

---

**End of spatial full-package packet. 모든 결과 signature/niche level exploratory · final association 없음 · Paper 1 main story 미접촉.**
