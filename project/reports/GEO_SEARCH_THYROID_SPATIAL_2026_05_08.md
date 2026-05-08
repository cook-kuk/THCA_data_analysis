# GEO / BioStudies / ENA 검색 — Thyroid Spatial Transcriptomics 추가 cohort

작성: 2026-05-08
작업 시간: 1시간 (사용자 요청)
범위: 공개 repo (GEO · BioStudies/ArrayExpress · ENA · Zenodo · HRA/GSA · Europe PMC)
이미 보유: GSE250521 (16) · GSE230424 (4) · GSE248205 (8) = **28 슬라이드 / 3 cohort**

---

## 1. 한 줄 결론 (TL;DR)

> **공개 repo 에 human thyroid carcinoma Visium spot-level cohort 추가로 사용 가능한 것 = 0 개.**
> 우리는 이미 **공개된 human thyroid Visium 전체 corpus** 를 보유 중. JCI Insight reach 를 위해 기대했던 "1 public Visium cohort 추가" 가능성은 가용하지 않다.

추가 cohort 가 필요하면: **Bundang prospective FFPE outreach** 만 path 임 (Marathon 6주 안 도착 어려움).

---

## 2. 검색 결과 정량

| Repo | Query 종류 | unique hits | 진짜 thyroid Visium cohort |
|---|---|---|---|
| GEO (E-utilities `gds`) | 15 queries (thyroid + spatial / Visium / CosMx / Xenium / ST / Hashimoto / Graves / ATC) | 79 unique UIDs | **0 추가** |
| Europe PMC | 9 queries (papillary thyroid / Visium / Slide-seq / FFPE / 2024-2026) | 약 200 results | **0 추가** |
| BioStudies / ArrayExpress | 5 queries (thyroid + Visium / CytAssist / spatial transcriptomics) | 약 250 results | **0 추가** |
| ENA (read_run + study_title filter) | thyroid + Visium / spatial transcriptomics | **0 rows** | 0 |
| Zenodo | thyroid + Visium | 0 thyroid match | 0 |
| HRA / GSA (NGDC China) | thyroid + spatial | 0 indexed | 0 |

**Output files**:
- `project/results/geo_search_2026_05_08/geo_search_results.tsv` (79 UIDs, full metadata)
- `project/results/geo_search_2026_05_08/geo_candidates.tsv` (36 thyroid + spatial UIDs, mostly noise)
- `project/results/geo_search_2026_05_08/geo_candidates_drilldown.tsv` (9 manual drill-down)
- `project/results/geo_search_2026_05_08/GSE301163_samples.tsv` (78 GeoMx ROIs)
- `project/results/geo_search_2026_05_08/epmc_round2.tsv`
- `project/results/geo_search_2026_05_08/biostudies_round2.tsv`

---

## 3. 자세히 — 후보 별 verdict

### 3.1 ❌ GSE301163 — 78 samples, **GeoMx DSP** (NOT Visium)

- Title: *"Tracing the molecular route to progression in miRNA biogenesis-defective thyroid lesions [Spatial_transcriptomics]"*
- Public: 2025-12-10 · PubMed 41657311
- Organism: Homo sapiens
- Platform: GPL18573 (Illumina HiSeq 2500)
- **Modality**: Each "sample" = ROI (Region of Interest) × {PanCK+ epithelial OR VIM+ stromal} segmentation. ROI-level pooled gene expression. **NOT spot-level Visium.**
- Sample structure (Patient 1 example): Normal thyroid (3 ROI × 2 compartment = 6) → Non-neoplastic nodular (3 ROI × 2 = 6) → PDTC DICER1-mutated (6 ROI × 2 = 12) → Tumor Capsule (6 ROI × 1 VIM+ = 6)
- Disease coverage: Normal · Non-neoplastic nodular · PDTC (DICER1-mutated) · Tumor Capsule
- **Use case (cross-platform external validation, Plan B)**:
  - PanCK+ epithelial compartment 으로 TROP2 / DM1 axis ROI-level validation 가능
  - VIM+ stromal compartment 으로 stromal niche 분리 가능
  - DICER1 mutation 을 orthogonal genotype axis 로 활용 가능
  - **Moran I 직접 비교 불가** (ROI-level pooled, no spot grid)
- **권장**: Paper 2 supplementary 의 **adjacent platform external validation** 으로 검토 가능 (n=78 ROI from real human thyroid lesions). 단 main claim 진입은 제한적.

### 3.2 ❌ GSE231952 / GSE231953 — Mouse, NF-κB thyroid dysgenesis

- Mus musculus, n=4 / n=1 spatial slides
- 발달 학적 thyrocyte migration 연구
- **Use**: biology citation only (Paper 2 Discussion 에 인접 reference 가능). 데이터 통합 불가 (mouse).

### 3.3 ❌ GSE308145 — FFPE benchmark, n=1, off-topic

- "Systematic benchmarking of imaging spatial transcriptomics platforms in FFPE tissues - scRNA-seq Data"
- 8-plex 작은 slide, methodology paper
- Thyroid carcinoma 관련 X
- **Use**: 0

### 3.4 ❌ GSE305978 / GSE305186 — Fetal brain neurogenesis (CosMx), thyroid hormone biology, off-topic

- Maternal T4 / fetal brain
- Thyroid carcinoma 관련 X

### 3.5 ❌ GSE145958 / GSE273927 — Mouse, off-topic

- TSH receptor learning · uterine RTHa
- 데이터 무관

### 3.6 ❌ Europe PMC `thyroid Visium` 결과 (50+) — 모두 **brain / glioblastoma / HCC / breast / cardiac / liver / opisthosoma / spider** Visium

- 키워드 매칭 noise. 실제 thyroid Visium paper 0 건.

### 3.7 ❌ Cell line / organoid (S-CMO* in BioStudies) — Visium 아님, organoid bulk

### 3.8 ❌ Other GEO thyroid carcinoma — bulk RNA-seq / microarray / miRNA / ChIP / WES — Visium 아님

---

## 4. 왜 이런 결과인가

- **Visium 은 비교적 최근 platform** (10x Genomics 2019 출시, FFPE CytAssist 2022). 갑상선암은 fresh-frozen 이 어렵고 (small thyroidectomy specimen) FFPE 가 표준 → Visium 적용 늦음.
- **갑상선암 spatial 연구의 publication bottleneck** = 데이터가 deposit 안 됨 (preprint stage 또는 controlled-access).
- **GSE250521 (16)** 가 사실상 가장 큰 공개 thyroid Visium cohort. **GSE230424 (4)** = HT-overlap PTC 유일. **GSE248205 (8)** = autoimmune-only 유일.
- 우리 28-sample 패키지는 이미 **공개 corpus 의 100% 포함**.

---

## 5. 영향 — venue tier 재평가

전 분석에서 (`PAPER2_SPATIAL_CURATION_DECISION_2026_05_07.md`):
- Sci Rep base / JCI Insight reach (≥1 추가 cohort 또는 H&E pathology layer)
- Cell Rep Med = NOT recommended

오늘 검색 결과 → **"1 public Visium cohort 추가" path 가 막힘**.

| Path | 현 상태 |
|---|---|
| **Sci Rep base** | 그대로 가능 — 28-sample / 3 cohort 가 공개 corpus 전체 → 오히려 강점 framing 가능 ("we exhaust the publicly available human thyroid Visium corpus") |
| **JCI Insight reach** | **public Visium 추가는 불가**. 단 하나 남은 path = **GSE301163 GeoMx DSP cross-platform external validation + Bundang H&E pathology overlay (medical imaging 본인 강점)** |
| **Cell Rep Med +** | functional validation 없이 불가. 기존 진단 그대로 |

**중요 reframe 가능 표현 (Sci Rep submission 강점)**:
> *"We integrate 28 slides across 3 publicly available cohorts — to our knowledge the entire publicly indexed human thyroid Visium corpus as of 2026-05 — providing comprehensive spatial coverage of normal thyroid, papillary/anaplastic thyroid carcinoma, Hashimoto-overlap PTC, and autoimmune-only baseline (Hashimoto, Graves, control)."*

→ "We exhaust the public corpus" 라는 framing 자체가 **n 약점을 cohort completeness 강점으로 전환** 가능.

---

## 6. 다음 step 권장 (ranked)

| # | Action | Cost | Lift |
|---|---|---|---|
| 1 | **현 상태로 Paper 2 freeze + Sci Rep 정조준** + "exhausts public corpus" framing | 0 (이미 보유) | **+** Sci Rep submission 신뢰도 + reviewer pushback 약화 |
| 2 | **GSE301163 GeoMx DSP** ROI-level external validation 추가 — TROP2 / DM1 axis 가 PanCK+ epithelial compartment 에 specific 하게 응집하는지 확인 | 1–2일 분석 (78 ROI 작음) | **+** JCI Insight reach 살짝 회복 — cross-platform support |
| 3 | **Bundang prospective FFPE outreach** 계속 push (timing risk; Marathon 안 어려움) | weeks–months | reach venue 외부 |
| 4 | **TCGA-THCA H&E + DM1 closure 재시도** with stronger model (이미 NO-GO 였으나 재 audit 가능) | 1주 | low (이미 NO-GO 명확) |
| 5 | **Functional validation (organoid / IHC TROP2)** — Marathon 안 불가, Paper 2 v2 후속용 | months | Cell Rep Med 까지 reach |

**가장 현실적 path**: **Action 1 + Action 2** 조합. Public corpus exhaustion framing + GSE301163 GeoMx ROI-level external validation = JCI Insight 가능성 일부 복원. Marathon 안에 가능.

---

## 7. Discipline self-audit

| Rule | Compliance |
|---|---|
| Voice-protected sections 무손상 | ✓ (본 문서 = decision/research memo only) |
| Yu 분리벽 (Paper 2 HT vs Paper 4 GD) | ✓ (GSE301163 후보 도입 시도 시 Hashimoto 별도 데이터 무관) |
| Marathon mode = paper-blocking only | ✓ (reach venue 결정 직접 영향) |
| Signature-level only | ✓ |
| Paper 1 main story 무손상 | ✓ |
| 새 large analyses 시작 금지 | ✓ (이번은 search/decision; analysis 시작 X) |

---

## 8. 한 페이지 요약

- 1시간 search 결과: **추가 가능한 public Visium thyroid cohort = 0**
- GSE301163 (n=78) = GeoMx DSP, Visium 아님 → cross-platform validation 으로만 사용 가능
- Paper 2 reach venue 도달은 **public Visium 추가 path 가 닫혔다**
- 살아있는 path: (a) "exhaust public corpus" framing 으로 Sci Rep base 강화, (b) GSE301163 GeoMx ROI-level external validation, (c) Bundang prospective (Marathon 후)
- **추천**: Action 1 + 2 조합 — 1–2일 추가 작업으로 JCI Insight reach 가능성 일부 복원

— `GEO_SEARCH_THYROID_SPATIAL_2026_05_08.md`
