# Per-figure data audit — Paper 1 manuscript v8 (2026-05-07)

**Purpose:** 각 figure가 어떤 dataset을 썼는지 명시적 검증. caption 원문 (`05_figure_captions.md` + `05_supp_figure_captions_v17_dm.md`) 기준으로 추출.

**Notation:** ★ = 다운로드/처리 완료 / ⚠ = caption에 등장하지만 본 repo 내 가공된 데이터로 verify 필요 / ✗ = manuscript 등장하지만 사전 v3에 빠져있던 항목 (이번에 보강 완료) / ☐ = request-required (Section E).

---

## 1. Main Figures (Fig 1–8)

| Fig | 핵심 분석 | Dataset | n | Modality | Source | Sprint status |
|-----|-----------|---------|---|----------|--------|----|
| **1A** | TIERA67 → 8-gene curation Sankey | (no data; schematic) | — | — | analyst-curated | ★ |
| **1B** | UMAP DM1/DM2 KMeans | TCGA-THCA | 504 | RNA-seq | TCGA / GDC | ★ |
| **1C** | BRAF V600E vs WT mRNA neutrality + driver AUC | TCGA-THCA | 273 V600E + 182 WT (=455) | RNA-seq + mut | TCGA / cBioPortal | ★ |
| **1D** | Pan-genome ARI ladder (8 / TIERA67 / top-5000 / driver_anchor) | TCGA-THCA | 504 | RNA-seq | TCGA | ★ |
| **2A** | Xing 2014 dark-matter rescue Sankey (131/180) | TCGA-THCA Xing dark-matter subset | 180 → 131 rescued | RNA-seq + mut | TCGA + Xing 2014 list | ★ |
| **2B** | DM1 prevalence by ancestry (28.4% vs 37.8%) | TCGA-THCA + Korean pool (K2 + Lee/GSE213647 + GSE286332) | 504 + 874 | RNA-seq | TCGA / PRJEB11591 / GSE213647 / GSE286332 | ★ (✗ GSE213647 Lee 2024 v3 보강) |
| **2C** | KM OS DM1 vs DM2 within Xing dark matter | TCGA-THCA dark-matter | 180 | clinical + RNA-seq | TCGA | ★ |
| **3A** | Per-patient Spearman r tumor↔normal thyrocyte | **GSE184362 (Pu 2021)** | 6 patients / 158K cells | scRNA | GEO Pu et al. Nat Commun 2021 | ★ (✗ v3 보강) |
| **3B** | UMAP thyrocyte cluster DM1 score | **GSE193581 (Lu 2023)** | 23 samples | scRNA | GEO Lu 2023 | ★ (✗ v3 보강) |
| **3C** | Author-independence cross-cohort | **GSE241184 (Phase 1)** vs GSE184362 | 1 vs 6 | scRNA | GEO | ★ (✗ v3 보강) |
| **4A** | 8-cell stack (BRAF × TERT × DM) | TCGA-THCA + cBioPortal `thca_tcga_pub` (TERT⁺ recovery n=36) | 504 | RNA-seq + mut | TCGA + cBioPortal | ★ (✗ cBioPortal v3 보강) |
| **4B** | KM OS BRAF_TERT⁺ vs BRAF/TERT-neg dark matter × DM | TCGA-THCA | 504 | clinical + mut | TCGA | ★ |
| **4C** | OTHER_TERT⁺ N=4 caveat box | TCGA-THCA | 4 | mut + clinical | TCGA | ★ |
| **5A** | Korean cohort 8-gene DM-score distribution | **GSE213647 (Lee 2024)** | **632** | RNA-seq | GEO Lee et al. | ★ (✗ v3 보강) |
| **5B** | K2 NBNR DM cluster composition + vascular invasion | K2 / PRJEB11591 | 260 | RNA-seq | ENA Yoo 2016 | ★ |
| **5C** | FFPE vs FF concordance density (KS p=0.44) | TCGA-THCA + (FFPE source verify) | — | RNA-seq | TCGA + ⚠ FFPE source | ⚠ FFPE source verify |
| **6A** | Forest Cox HR DM1 vs DM2 | TCGA-THCA + **MSK-IMPACT** (Landa 2016) | 504 + 117 | clinical + mut | TCGA + cBioPortal `thca_mskcc_2016` | ★ (✗ MSK-IMPACT v3 보강) |
| **6B** | I² heterogeneity panel | derived from 6A | — | meta-stats | — | ★ |
| **6C** | DM1 sub-A vs sub-B Cox | TCGA-THCA DM1 sub-set | 91 (4 events) | clinical + RNA-seq | TCGA | ★ |
| **7A** | DM1 sub-A vs sub-B silhouette | TCGA-THCA DM1 | 72 + 19 = 91 | RNA-seq | TCGA | ★ |
| **7B** | Fusion-positivity DM1 vs DM2 (OR 7.41) | TCGA-THCA + cBioPortal SV table | 504 (SV-tested 489) | structural variants | cBioPortal | ★ |
| **7C** | Fusion partner stack within DM1 | TCGA-THCA fusion+ DM1 | 63 | SV | cBioPortal | ★ |
| **7D** | DM1 sub-A vs sub-B age/stage/CD8/IFNγ/checkpoint panel | TCGA-THCA DM1 | 72 + 19 | clinical + RNA-seq | TCGA | ★ |
| **7E** | DM1 captures 81.8% TCGA RET-fusion+ | TCGA-THCA RET-fusion subset | 27/33 | SV | cBioPortal | ★ |
| **8A** | HM450 promoter β heatmap (8-gene + DIO2 + SLC26A4) | **TCGA-THCA HM450** | **503** | methylation array | TCGA | ★ (✗ HM450 v3 보강) |
| **8B** | Mean panel β by DM cluster | TCGA-THCA HM450 | 503 | methylation | TCGA | ★ |
| **8C** | Methylation fusion-independence within DM1 | TCGA-THCA HM450 ∩ DM1 | 63 fusion+ vs 19 fusion− | methylation + SV | TCGA + cBioPortal | ★ |

---

## 2. Supplementary Figures S1–S9 (manuscript core)

| Fig | 핵심 분석 | Dataset | n | Source |
|-----|-----------|---------|---|--------|
| **S1** | 8-gene heatmap sorted by P_DM1; TCGA vs K2 distribution | TCGA-THCA + K2 / PRJEB11591 | 504 + 260 | TCGA + ENA |
| **S2** | 8-gene panel similarity matrix sc external | GSE184362 + GSE193581 + GSE241184 | 6 + 23 + 1 | GEO |
| **S3** | Pan-genome ARI ladder (8/16/67/200/1000/5000) | TCGA-THCA | 504 | TCGA |
| **S4** | Immune-residualization Δd (raw 1.78 → dual 0.87) | TCGA-THCA + immune signature module | 504 | TCGA + signature library |
| **S5** | SV missingness sensitivity (15 SV-untested) | TCGA-THCA | 489 + 15 | cBioPortal |
| **S6** | DM1 sub-B age/stage/immune detail | TCGA-THCA DM1 sub-set | 72 + 19 | TCGA |
| **S7** | Cross-cohort DM-score portability | K2 + Lee/GSE213647 + Korean reference (small) | 260 + 632 + (verify) | ENA + GEO + ⚠ small Korean ref |
| **S8** | Decitabine + I-131 schematic + literature meta | NCT00085293, NCT01065090 (literature only) | n/a | citation only |
| **S9** | K2 mini-index calibration diagnostic (R4-4) | K2 vs Yoo 2016 reference | 260 + Yoo 2016 panel | ENA + Yoo 2016 |

---

## 3. Additional supp — v17 audit (SA1–SA6) + Dark matter (SB1–SB3, SC1–SC7) + Landa (SD1)

### Part A — v17 Audit phase

| Fig | 핵심 분석 | Dataset | n |
|-----|-----------|---------|---|
| **SA1** | 4-way BRAF × TERT revalidation | TCGA-THCA + cBioPortal TERT⁺ | 504 (BRAF−/TERT⁺ n=4) |
| **SA2** | FFPE vs FF QC (KS p=0.44) | TCGA-THCA + ⚠ FFPE source | verify |
| **SA3** | Panel-size sensitivity (8/10/12/16) | TCGA-THCA | 504 |
| **SA4** | Single-cell wrap-up (3 sc cohorts) | GSE184362 + GSE193581 + GSE241184 | 6 + 23 + 1 |
| **SA5** | Pseudotime trajectory along 8-gene | TCGA-THCA | 504 |
| **SA6** | MSK-IMPACT bias panel (84 PDTC + 33 ATC) | **MSK-IMPACT (Landa 2016)** | **117** |

### Part B — Dark matter phase 1 (single-cell foundation)

| Fig | 핵심 분석 | Dataset | n |
|-----|-----------|---------|---|
| **SB1** | sc UMAP overview by patient + cell-type + tumor/normal | GSE184362 (Pu 2021) | 6 patients |
| **SB2** | sc 8-gene + HLA-II + thyroid TFs + proliferation | GSE184362 | 6 patients |
| **SB3** | Thyrocyte-restricted UMAP (KRT8+ KRT19+ EPCAM+) | GSE184362 | 6 patients |

### Part C — Dark matter phase 2 (multisite)

| Fig | 핵심 분석 | Dataset | n |
|-----|-----------|---------|---|
| **SC1** | Per-patient Spearman r forest tumor↔normal thyrocyte | GSE184362 | 6 |
| **SC2** | Multisite trajectory DM1/DM2 score | TCGA-THCA + MSK-IMPACT + K2 + Lee/GSE213647 | 504 + 117 + 260 + 632 |
| **SC3** | Joint sc UMAP (batch-corrected, ⚠ method verify) | GSE184362 + GSE193581 + GSE241184 | 6 + 23 + 1 |
| **SC4** | GSE184362 per-patient summary table | GSE184362 | 6 |
| **SC5** | External-validation pooled scatter (tumor vs normal) | GSE184362 + GSE193581 | 6 + 23 |
| **SC6** | Pooled scatter DM1 prob vs 8-gene score | GSE184362 + GSE193581 | 6 + 23 |
| **SC7** | K2 vs Yoo 2016 reference R4-4 audit | K2 / PRJEB11591 vs Yoo 2016 panel | 260 |

### Part D — Landa 2016 cite save

| Fig | 핵심 분석 | Dataset | n |
|-----|-----------|---------|---|
| **SD1** | GSE76039 differentiation transcript heatmap (TG, TSHR, TPO, PAX8, SLC26A4, DIO1, DUOX2) | **GSE76039 (Landa 2016)** | ~37 PDTC+ATC RNA-seq subset |

---

## 4. Cross-figure dataset usage matrix

| Dataset | Used in figures |
|---------|-----------------|
| **TCGA-THCA RNA-seq (n=504)** | Fig 1B, 1C, 1D, 2A, 2B (component), 2C, 4A, 4B, 4C, 6A, 6C, 7A, 7B, 7C, 7D, 7E, S1, S3, S4, S5, S6, SA1, SA3, SA5, SC2 |
| **TCGA-THCA HM450 (n=503)** ✗ | Fig 8A, 8B, 8C |
| **cBioPortal `thca_tcga_pub` (TERT⁺ n=36, SV)** ✗ | Fig 4A, 4B, 7B, 7C, 7E, 8C, S5, SA1 |
| **MSK-IMPACT thyroid / Landa 2016 (n=117)** ✗ | Fig 6A, SA6, SC2 |
| **K2 / PRJEB11591 (n=260)** | Fig 2B, 5B, S1, S7, S9, SC2, SC7 |
| **Lee 2024 / GSE213647 (n=632)** ✗ | Fig 2B (component), 5A, S7, SC2 |
| **GSE286332 (n=18)** | Fig 2B (component) |
| **GSE76039 (Landa 2016, ~37)** | SD1 |
| **GSE184362 (Pu 2021, n=6)** ✗ | Fig 3A, 3C, S2, SA4, SB1, SB2, SB3, SC1, SC3, SC4, SC5, SC6 |
| **GSE193581 (Lu 2023, n=23)** ✗ | Fig 3B, S2, SA4, SC3, SC5, SC6 |
| **GSE241184 (Phase 1, n=1)** ✗ | Fig 3C, S2, SA4, SC3 |
| **Yoo 2016 panel reference** | S9, SC7 |
| **NCT00085293, NCT01065090 (decitabine + I-131 lit)** | S8 (citation only) |
| **⚠ FFPE source** | Fig 5C, SA2 — verify within-repo source |
| **⚠ Small Korean reference (S7)** | S7 — verify (likely subset of K2 calibration arm) |

---

## 5. Per-figure verification flags

| Flag | Fig | Issue |
|------|-----|-------|
| ⚠ | 5C, SA2 | FFPE 시료 source (어떤 cohort의 FFPE인지) — repo 내 정확한 출처 verify 필요 |
| ⚠ | S7 | "small Korean reference cohort" 정체 — K2 calibration arm subset인지 별도 cohort인지 verify |
| ⚠ | SC3 | Joint sc UMAP integration method (Harmony / scVI / scANVI) — caption에 verify 표시 |
| ✗ → ★ | 3A/B/C, 4A, 5A, 6A, 8A/B/C, SA4, SA6, SB*, SC1–7, SD1 | v3 사전에 cohort 보강 완료 (Lee 2024, MSK-IMPACT, Pu 2021, Lu 2023, GSE241184, cBioPortal, HM450) |
| ☐ | (no main figure currently uses request-required data) | dbGaP / EGA / pharma 항목은 Paper 3 Track B 또는 Paper 4 backlog에서 활용 예정 — 본 audit 시점 기준 Paper 1 figure 어디에도 직접 등장 안 함 |

---

## 6. 결정 행렬 (Paper 1 figure 관점)

| 항목 | 영향 | 결정 대기 |
|------|------|-----------|
| Fig 5C / SA2 FFPE source verify | caption에 cohort 명시 필요 | repo 내 FFPE script grep |
| Fig S7 "small Korean reference" 정체 | caption 명시 필요 | K2 calibration arm 자료 verify |
| Fig SC3 integration method | caption에 ⚠ 그대로 → 명시로 변경 | Harmony / scVI 중 실제 사용 method confirm |
| MSK-IMPACT (n=117) public access path 명시 | Methods section 강화 | cBioPortal `thca_mskcc_2016` 재확인 |
| GSE213647 access path | Data Availability 강화 | GEO 직접 access 가능 (controlled X) confirm |

---

## 7. 산출물 / 다음 행보

- `project/reports/2026_05_07_per_figure_data_audit.md` (이 파일)
- 위 ⚠ 3건은 **manuscript caption** 수정 사항 — voice-protected 영역 아님 (Methods/caption verify는 마라톤 모드 허용).
- 사전 v3는 이미 모든 ✗ 항목 보강 완료 (47 페이지 재주입됨).
- 새 분석 명령 안 함. 새 데이터 다운로드 안 함. Commit 안 함.
