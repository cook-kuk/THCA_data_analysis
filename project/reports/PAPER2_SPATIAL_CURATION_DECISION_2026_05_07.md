# Paper 2 Spatial Curation Decision — 2026-05-07 (UPDATED 2026-05-08)

작성: 2026-05-07 / 갱신: 2026-05-08
저자: Seungho Cook (논문)
지도교수: Yu Hyeong-won
Mode: Marathon writing (2026-05-04 → 2026-06-13)
Boundary: Paper 2 = HT-overlap PTC ONLY · Paper 4 backlog = Korean GD HLA · Yu 2026-05-04 분리벽

---

## ★ 2026-05-08 UPDATE — GeoMx DSP cross-platform external validation 추가

(상세는 §10 + `GEO_SEARCH_THYROID_SPATIAL_2026_05_08.md`)

- **GEO 1시간 검색 결과**: 추가 가능한 public Visium PTC cohort = 0개. **공개 corpus 전체 보유** 확인 → "exhausts public Visium thyroid corpus" 강점 framing 가능.
- **GSE301163** (n=78 GeoMx DSP, 2 patients, microPTC 11 / PDTC 18 / Normal 18 / Nodular 18 / Capsule 6) 추가 분석 완료.
- **Cross-platform validated finding** (signature-level, MW p):
  - DM1_axis PDTC < microPTC PanCK+ **p=0.001** (3.98 vs 4.39) → v12 differentiation collapse 직접 재현
  - RAI_8 same direction p=0.001
  - HLA-II microPTC > Normal PanCK+ p=0.030
  - 8-gene DM1 microPTC vs Normal PanCK+ p=0.017
- **Honest negative**: TROP2 ROI-mean PanCK+ vs VIM+ ns (microPTC p=0.714). ROI-pooled 측정으로는 Visium spot-level niche 직접 검출 불가. Visium claim 을 refute 하지 않음 (다른 layer).
- **6 figures 추가** (S_F101–S_F106). MAIN 후보 = S_F106 (9th main).
- **Tier**: MAIN=**9** / SUPP=**33** / WEB=**64**. Total 106.

---

## 1. Executive decision (2026-05-08 갱신)

| 항목 | 결정 |
|---|---|
| **Paper 1** | submit-ready 유지. 본문/figure/supp 무수정. 4 outreach 발송 + npj submit 클릭은 user 의 키보드. |
| **Paper 2 base venue** | **Scientific Reports** — "exhausts public Visium thyroid corpus" framing + GeoMx cross-platform consistency 가 강점화. |
| **Paper 2 reach venue** | **JCI Insight** — 도달 가능성 회복 (GSE301163 GeoMx ROI validation 의 v12 finding 재현 + cross-platform 4-panel summary). |
| **Cell Reports Medicine** | **NOT 추천** (still). Functional / IHC 검증 부재. 단 Bundang prospective FFPE 도착 시 재평가. |
| **v12 findings (분화 공간붕괴 ⭐ cross-platform 재현됨)** | supp 진입 + cross-platform 강화. Main claim 에는 차분히 (n=4 ATC + 2 GeoMx patient caveat). |
| **Spatial-full 106 figures** | curated supplementary engine (9 main + 33 supp + 64 web-only archive). |

추천 route 한 줄:
> **Paper 1 click → Paper 2 base (Sci Rep) curation freeze → independent validation 1개 확보 → JCI Insight reach**

---

## 2. Current grade (2026-05-08 갱신)

| Subject | Grade | Comment |
|---|---|---|
| Paper 1 (DM1 dark matter, 8-gene + DM1 axis) | **A− / npj-ready** | submission package frozen; outreach 4건 발송 + click 대기 |
| Paper 2 base (HT-overlap PTC, Pillar I + Pillar II) | **B+ → A− 도달** (post GeoMx) | post-fix Moran + cross-platform GeoMx external = JCI Insight reach 도달 가능 |
| Spatial-full standalone | **A− as supplementary engine** | 106 figures 단독 publication 어려움; supplementary engine 으로는 매우 강 |
| v12 findings (분화 공간붕괴 cross-platform 재현 ★) | **EXPLORATORY → MODERATE_SUPPORT** (post-platform) | DM1/RAI_8 PDTC vs microPTC p=0.001; v12 finding 신뢰도 상승 |
| Cross-platform consistency (Visium ↔ GeoMx) | **STRONG (only on differentiation/immune axes)** | TROP2 spot-level niche 는 platform 차이로 GeoMx 측정 불가; honest caveat |
| Closure battery (H&E → DM1 NO-GO) | **STRONG (LOSO ρ ≈ 0.06)** | Paper 2A NO-GO self-evidence; 사용 안전 |

---

## 3. Main figure shortlist (9, 2026-05-08 갱신)

| # | Figure | Question answered | Cohort | Claim |
|---|---|---|---|---|
| **M1** | **S_F7** TROP2 maps GSE250521 (16-slide) | Why TROP2? (visual proof) | GSE250521 16 | TROP2 is spatially organized in PTC + LPTC |
| **M2** | **S_F11** Cross-cohort TROP2 heatmap | Why TROP2? (quantitative) | ALL 28 / 3 cohorts | TROP2 niche tumor-specific (PTC+LPTC 8/8 Moran's I 0.27–0.58; PT/autoimmune ~0) |
| **M3** | **S_F8** HT TLS maps GSE230424 (4-slide) | Why HT-overlap PTC special? | GSE230424 4 | HT TLS niche directly visible 4/4 slides |
| **M4** | **S_F18** HT TLS Moran I heatmap | Why HT-overlap PTC special? (quantitative) | GSE230424 4 | Moran's I 0.27–0.84 across HT immune axes |
| **M5** | **S_F19** DM1 axis maps GSE250521 (16-slide) | Why spatial transcriptomics? | GSE250521 16 | DM1 region heterogeneity exists at spot level |
| **M6** | **S_F30** Closure battery LOSO breakdown | Why H&E → DM1 closure failed? | GSE250521 + TCGA H&E | LOSO ResNet50 ρ ≈ 0.06 (NO-GO self-evidence) |
| **M7** | **S_F36** GSE248205 full maps | Why this is signature-level not overclaimed? | GSE248205 8 | Autoimmune-only negative control: no TROP2/TLS niche |
| **M8** | **S_F56** Four-key-finding collage | Graphical abstract | ALL 28 / 3 cohorts | Combined Q1+Q2+Pillar I+Pillar II in one panel |
| **M9** ★ | **S_F106** Cross-platform consistency summary | External validation? | Visium 28 + GeoMx 78 | Visium spot-level Moran ↔ GeoMx ROI-mean 일치 (differentiation/immune axes) |

> **선정 원칙**: 5개 question 1:1 대응 + cross-platform external validation (M9). ATC-2 / LPTC-3 outlier 는 main 에서 제외. v12 finding 은 main 에서 제외 (단 supp S_F91/S_F97/S_F98/S_F99 에 cross-platform GeoMx 재현 caveat 추가).

---

## 4. Supplementary figure shortlist (33 = 28 Visium + 4 GeoMx external + S_F31 promotion, 2026-05-08 갱신)

> **2026-05-08 추가 SUPP**: S_F101 (GeoMx 78×14 heatmap) · S_F103 (DM1 by histology, p=0.001 cross-platform) · S_F104 (8-gene + genotype) · S_F105 (3-axis differentiation collapse, p=0.001/0.001/0.007). 추가로 S_F31 (post-fix hex 6-nbr remake) 가 web → supp 로 promoted.


| # | Figure | Role | Note |
|---|---|---|---|
| S1 | S_F2 | Driver-orthogonal spot-level box | Pillar I support |
| S2 | S_F5 | TROP2 spot distribution per sample | Pillar I support |
| S3 | S_F9 | TROP2 ranked Moran I per sample | Quantitative companion to M1 |
| S4 | S_F10 | TROP2 × DM1 spatial co-localization | TROP2-DM1 link |
| S5 | S_F13 | Closure battery summary | Companion to M6 |
| S6 | S_F22 | Moran I per condition box | Cohort separation |
| S7 | S_F24 | HT immune marker co-localization | Pillar II support |
| S8 | S_F25 | TROP2 niche cluster size distribution (post-fix) | Hex 6-nbr; technical QC |
| S9 | S_F29 | TROP2 maps GSE230424 | Cross-cohort TROP2 |
| S10 | S_F38 | Per-sample QC | Methods |
| S11 | S_F40 | Moran axis correlation | Methods |
| S12 | S_F41 | TROP2 28 samples mega panel | Pillar I overview |
| S13 | S_F46 | Niche threshold sensitivity | Methods QC |
| S14 | S_F47 | Per-axis Moran I summary | Compact form of S_F64 |
| S15 | S_F48 | Depth confounding control | Methods QC |
| S16 | S_F49 | DM1 28 samples mega panel | Pillar II overview |
| S17 | S_F50 | Niche robustness (post-fix) | Methods QC |
| S18 | S_F53 | HT full-axis Moran summary | Pillar II quantitative |
| S19 | S_F61 | Closure target × model heatmap | Extends M6 |
| S20 | S_F63 | 28-sample summary table | Methods table |
| S21 | S_F64 | Master 28-sample × 11-axis Moran I matrix | Overview |
| S22 | S_F65 | TLS individual chemokine receptor maps (P3) | Pillar II supp; n=1 case caveat |
| S23 | S_F67 | Sample PCA in axis-Moran space | Overview |
| S24 | S_F71 | Within-slide cross-axis Spearman heatmap | Overview |
| S25 | S_F91 | Thyroid differentiation spatial coherence (v12) | Exploratory finding (n=4 ATC caveat) |
| S26 | S_F97 | 24-axis stage gradient lineplot (v12) | Exploratory finding |
| S27 | S_F98 | ATC-2 8-axis deep dive (v12) | Case-like outlier (n=1 caveat) |
| S28 | S_F99 | Cross-axis spatial co-occurrence Jaccard (v12) | Overview |

> **v12 finding 4개** (S_F91, S_F97, S_F98, S_F99) 는 supp 에 포함되지만 caption 에서 **exploratory · n caveat** 명시 필수.

---

## 5. Web-only / archive list (2026-05-08 갱신)

| Tier | Count | Rationale |
|---|---|---|
| MAIN_FIGURE | **9** | 8 Visium + S_F106 cross-platform (M9) |
| SUPPLEMENTARY_FIGURE | **33** | 28 Visium + S_F101/S_F103/S_F104/S_F105 GeoMx external + S_F31 promotion |
| WEB_ONLY_ARCHIVE | **64** | Redundancy / superseded / exploratory; S_F102 GeoMx TROP2 compartment (honest negative) |
| DROP_OR_REMAKE | 0 | (after curation + remake) |
| **Total** | **106** | 100 Visium + 6 GeoMx external (S_F101–S_F106) |

WEB_ONLY_ARCHIVE 64개 구성 (paragraph 요약):
- v1 redundant (S_F1, 3, 4, 6) — superseded by v2/v3/v5 main+supp
- v2 redundant (S_F12, 14, 15, 17) — superseded by S_F19/S_F11/S_F36
- v3 mega summary (S_F23, 26) — superseded by S_F56/S_F100
- v4 niche-cluster pre-fix-suspect (S_F31) — see Hex Adjacency Audit
- v4–v6 ATC collapse / RAI / epi / prol / Hypoxia/CAF/EMT 16-slide — biological context but redundant w/ S_F91 v12 supp
- v7 redundant fixes (S_F51, 52, 54, 55) — superseded
- v8 cluster IDs / GSE230424 RGB / TROP2 CDF — redundant
- v9 stage 3-axis spatial / niche compactness / ultimate summary (S_F69, 70, 72) — superseded
- v10 8 deep gene-set maps (S_F73–S_F79) — exploratory; S_F80 master also web-only
- v11 8 cancer pathways (S_F81–S_F88) — exploratory pathway maps; signature-level only
- v12 exploratory non-finding maps (S_F89, 90, 92–96, 100) — superseded by 4 supp v12 finding figures

---

## 6. Claim boundary

### 6.1 Allowed language
- *"spatially organized niche"*
- *"spatial coherence"*
- *"signature-level evidence"*
- *"exploratory axis"*
- *"case-like outlier"*
- *"hypothesis-generating"*
- *"region-specific signal"*
- *"tumor-specific spatial pattern"* (with Moran I cite)
- *"directly observable in space"* (vs bulk averaging)

### 6.2 Forbidden language (자동 grep 가드)
- ❌ biomarker
- ❌ clinical prediction
- ❌ predicts outcome
- ❌ patient stratification
- ❌ patient selection
- ❌ causal susceptibility / mechanism
- ❌ drug target
- ❌ therapeutic vulnerability (validated)
- ❌ companion diagnostic
- ❌ prognostic signature
- ❌ predictive marker

### 6.3 ATC-2 / LPTC-3 outlier 문구 제한
- 본문 main figure 에서 인용 X
- supp caption 에서만 *"case-like outlier (single slide; exploratory)"* 형식 사용
- "ATC-2 demonstrates" 같은 generalization 문장 금지

---

## 7. Reviewer risk map

| Risk | Severity | Mitigation |
|---|---|---|
| **R1. n=4 PTC+HT TLS** (M3 + M4) | HIGH | Caption 에 *"4-slide cohort, replication pending"* 필수. Limitations 본문 (voice-protected) 에서 별도 한 줄. |
| **R2. n=4 ATC differentiation collapse** (S_F91 v12) | HIGH | Supp only. *"exploratory; small ATC sub-cohort"* 명시. main claim 진입 금지. |
| **R3. n=1 ATC-2 outlier** (S_F98) | HIGH | Supp only. *"case-like single sample"* 명시. main claim 진입 금지. |
| **R4. S_F25 hex adjacency pre-fix figures** | MEDIUM | Hex Adjacency Audit (별도 문서) 따라 처리. SUPP/MAIN 진입 figure 는 모두 post-fix 만. |
| **R5. 100-figure firehose to advisor** | MEDIUM | 1-page brief + 8 main figure tier 만 forward. Supp tier 는 별도 ZIP. Web archive 는 link only. |
| **R6. Paper 1 timing risk** | MEDIUM | Paper 1 outreach 발송 + npj submit click 우선. Paper 2 spatial work 가 Paper 1 timing 을 추가 지연시키지 않도록. |
| **R7. v12 finding scope creep** | MEDIUM | "분화 공간붕괴" / "ATC-2 super-organizer" 를 Paper 2 main claim 으로 끌어올리지 않음. preprint 분리 옵션 보유. |
| **R8. Cross-paper inference (Paper 2 ↔ Paper 4 GD)** | MEDIUM | GSE248205 GD samples 는 negative control 로만. Paper 4 framing 인용 금지. |
| **R9. Bulk-level TROP2 enrichment failure 와 spatial niche claim 의 reframe** | LOW | "bulk averaging vs spatial niche" 표현은 reframe 으로 OK; "validated" 는 X. |

---

## 8. Next actions (ranked)

| # | Action | Owner | Gating |
|---|---|---|---|
| 1 | **Paper 1 outreach 4 발송 + npj submit click** | user | other actions blocked until completed |
| 2 | **Paper 2 figure curation freeze** — 8 main + 28 supp confirm | user (sign-off) | 본 문서 + TSV 검토 후 freeze |
| 3 | **S_F25 hex adjacency audit** — 별도 SPATIAL_HEX_ADJACENCY_AUDIT_2026_05_07.md 따라 figure 분류 | claude scaffolding | 2 와 병렬 가능 |
| 4 | **Advisor 1-page brief 전달** — ADVISOR_SPATIAL_BRIEF_2026_05_07.md | user | 2 sign-off 후 |
| 5 | **v12 supplement-vs-preprint 결정** | user | 2, 4 후 |
| 6 | **Independent validation 1 cohort 확보 plan** (JCI Insight reach gate) | user + Yu | 5 후 |
| 7 | **Paper 2 본문 figure 캡션 작성 (voice-protected = author keyboard only; non-voice = supp captions OK by claude)** | user (voice) / claude (non-voice) | 2 sign-off 후 |

> **Action 1 이 모든 다음 step 의 gate**. Marathon 1주차 (5/4–5/10) 마감 전 클릭 권장.

---

## 9. Discipline self-audit

| Rule | Compliance | 비고 |
|---|---|---|
| Voice-protected sections 무손상 | ✓ | 본 문서는 decision memo only |
| Yu 분리벽 (Paper 2 HT vs Paper 4 GD) | ✓ | GSE248205 GD samples 는 negative control 로만 사용 |
| Signature-level only | ✓ | claim boundary §6 명시 |
| Paper 1 main story 무손상 | ✓ | submission package 무수정 |
| Marathon mode = paper-blocking | ✓ | curation 은 Paper 2 진행 directly 보강 |
| n caveat 표시 | ✓ | n=4 PTC+HT, n=4 ATC, n=1 outlier 모두 §7 risk map 명시 |

— 본 문서 (`PAPER2_SPATIAL_CURATION_DECISION_2026_05_07.md`).

---

## 10. GeoMx GSE301163 cross-platform external validation (2026-05-08 추가)

### 10.1 Cohort

| 항목 | 값 |
|---|---|
| Accession | GSE301163 |
| Public | 2025-12-10 |
| PubMed | 41657311 |
| Platform | NanoString GeoMx Digital Spatial Profiler (GPL18573) |
| n ROIs | 78 |
| Patients | 2 (Patient 1 + Patient 2) |
| Histology | Normal thyroid 12 / Normal adjacent 6 / Non-neoplastic nodular 18 / microPTC 11 / PDTC 18 / Tumor Capsule 6 / immune ROI (macrophage 3 / lymphocyte 3 / high stromal 1) |
| Compartment | PanCK+ epithelial 32 / VIM+ stromal 46 |
| Genotype | DICER1-mutated 18 / DGCR8-E518K-mutated 40 / non-DICER1-mutated 20 |

### 10.2 Cross-platform stats (signature-level)

| Comparison | Direction | MW p | 의미 |
|---|---|---|---|
| DM1_axis PDTC PanCK+ vs microPTC PanCK+ | PDTC < microPTC (3.98 vs 4.39) | **0.001** | v12 differentiation collapse 직접 재현 ★ |
| RAI_8 PDTC PanCK+ vs microPTC PanCK+ | PDTC < microPTC (4.46 vs 4.83) | **0.001** | same direction ★ |
| Thyroid_TF PDTC PanCK+ vs microPTC PanCK+ | PDTC > microPTC (4.73 vs 4.28) | **0.007** | TF 보존, function loss (Landa 2016 일치) |
| HLA_II microPTC PanCK+ vs Normal PanCK+ | microPTC > Normal (2.84 vs 2.50) | 0.030 | microPTC immune activation |
| 8-gene_DM1 microPTC PanCK+ vs Normal PanCK+ | (mean 차이) | 0.017 | Paper 1 main signature cross-platform reactivity |
| TROP2 microPTC PanCK+ vs Normal PanCK+ | (mean 동일) | 0.714 (ns) | ROI-mean 으로 spot-level niche 검출 불가 (caveat) |
| TROP2 microPTC PanCK+ vs VIM+ | (mean 동일) | 0.714 (ns) | compartment specificity ROI-mean 차원에서 ns |
| TROP2 PDTC PanCK+ vs VIM+ | (mean 동일) | 0.251 (ns) | same |

### 10.3 Discipline 강조

- **TROP2 ns 결과는 Visium spot-level niche claim 을 refute 하지 않는다** — 측정 layer 가 다르다 (ROI-pooled vs spot grid + Moran's I). Limitations / Discussion 의 한 줄 caveat 권장.
- **Cross-platform 일치는 differentiation/immune axes 에서만** — TROP2 spatial niche 는 Visium 단독 evidence 로 남는다.
- **n=2 patients** caveat: 2 patient cohort 라는 점이 reach venue reviewer 의 추가 risk 항목.
- **DICER1/DGCR8 genotype** 은 orthogonal layer 로 사용 가능하나, Paper 2 main claim 인 HT-overlap PTC 와 cross-inference 금지 (genotype 별 PTC 분류는 Paper 4 backlog 와도 무관 — 별도 차원).

### 10.4 Output files

```
project/data/external_geomx_GSE301163/GSE301163_Normalized_ST_matrix.txt.gz   (input)
project/results/spatial_full_2026_05_06/GSE301163_meta_aligned.tsv             (78 ROI metadata)
project/results/spatial_full_2026_05_06/GSE301163_per_ROI_scores.tsv           (78 × 14 axes scores)
project/results/spatial_full_2026_05_06/GSE301163_crossplatform_stats.tsv      (cross-platform MW p)
project/papers_hub_2026_05_04/assets/spatial_full/S_F101_GSE301163_ROI_heatmap.png
project/papers_hub_2026_05_04/assets/spatial_full/S_F102_GSE301163_TROP2_compartment.png
project/papers_hub_2026_05_04/assets/spatial_full/S_F103_GSE301163_DM1_by_histology.png
project/papers_hub_2026_05_04/assets/spatial_full/S_F104_GSE301163_8gene_DM1.png
project/papers_hub_2026_05_04/assets/spatial_full/S_F105_GSE301163_differentiation.png
project/papers_hub_2026_05_04/assets/spatial_full/S_F106_cross_platform_summary.png
```

### 10.5 영향

- **Sci Rep base submission**: 가능성 강화 — "we triangulate spatial niche evidence with an adjacent platform (GeoMx DSP, n=78 ROI), confirming differentiation/immune axes consistency while acknowledging that ROI-pooled measurements cannot capture spot-level niche organization."
- **JCI Insight reach**: **현실적 가능성 회복**. 단 reviewer 는 n=2 patient cohort + ROI 가 작은 cohort 라는 점을 지적 가능 — 답변 가능 (PDTC 18 ROI / microPTC 11 ROI 자체는 적정).
- **Cell Reports Medicine**: 여전히 NOT (functional / IHC validation 부재).

— `PAPER2_SPATIAL_CURATION_DECISION_2026_05_07.md` (updated 2026-05-08)
