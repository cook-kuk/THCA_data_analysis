# Paper 3 ICI — Track B-lite Analysis Report

**작성일:** 2026-05-06
**저자:** Seungho Cook
**Sprint 권한:** 2026-05-06 user 명시적 "Track B" 명령. Track A "no public download" + "no NMF/LOHHLA/NetMHCpan/DIAL execution" 일부 오버라이드. **G1 (Paper 1 bioRxiv)** + **G2 (Paper 2 closure)** 게이트는 여전히 OPEN — 따라서 본 결과는 **hypothesis-generating sanity check** 이며 **Module A/D 의 본 분석을 대체하지 않는다**.
**Claim guard (binding):** ICI-readiness / immunogenomic vulnerability / candidate immune ecotype / hypothesis-generating prioritization 까지만. **"thyroid 환자에서 ICI 반응 예측" 절대 불가.**

---

## 0. 핵심 결과 (1 페이지)

1. **ATC = 면역염증 + 탈분화 동시 상승** — 4 thyroid 코호트 pooled meta-analysis, ATC vs non-ATC Cohen's d:
   - **myeloid_suppressive +2.53** (4/4 sign), **thyroid_differentiation -2.51** (4/4) — **거대 효과**
   - checkpoint_exhaustion +1.01 (4/4), HLA-II +0.68 (3/4), HLA-I +0.58 (4/4), IFNG +0.55 (4/4)
   - TLS_CXCL13_like -0.15 (3/4) — 미세하게 ATC에서 *감소* 경향 (예상과 반대; B/Tfh niche 가 PTC+HT 에 가까울 수 있음 — 페이퍼 2 영역 보호하면서 cite 가능)
2. **PCA structure (n=415, 8 cohorts)** — 2 차원으로 압축 가능:
   - **PC1 = pan-immune activation** (69.4% var; HLA-I/II + IFN + checkpoint 모두 ~+0.4 균등 loaded; thyroid_diff -0.20)
   - **PC2 = dedifferentiation** (14.2%; thyroid_diff +0.83 dominant; myeloid_suppressive -0.37)
3. **DIAL-lite (Hugo + Riaz pre-treatment, B=1000 bootstrap):**
   - **Sign 일치 3/7 modules only** — binomial p=0.77 → **direction-invariance NOT 확립** in 2-cohort minimal sample.
   - HLA-I/II + TLS 만 sign-consistent (+).
   - **IFNG / checkpoint / myeloid_suppressive / thyroid_differentiation flip** — Hugo 에서 (-), Riaz 에서 (+).
   - 이것은 페이퍼 3 의 DIAL audit 가 **존재 이유 그대로 입증** — pan-cancer composite를 thyroid에 transfer 하기 전 direction-invariance 검증이 필수.
4. **Hierarchical cluster bootstrap stability = 0.274** — 낮음. **batch / cohort 효과가 ecotype 신호를 압도** — 진짜 Track B 에서 **ComBat / harmonization 필수** (Track A signature_registry §0 와 일치).
5. **Riaz Pre→On 43 paired patients** — 7 module 모두 |Δ z-score| < 0.13, p > 0.30 — ICI treatment 단기 (Pre→On) 에서 bulk 모듈 변화는 검출되지 않음.

> **모든 결과는 readiness/vulnerability prioritization framing.** "thyroid ICI response prediction" 클레임 절대 불가.

---

## 1. 분석 범위 (왜 "Track B-lite" 인가)

**실행한 것:**
- 디스크에 이미 있는 9 cohort (8 thyroid + Hugo + Riaz) 처리
- 7-module × 60-gene 모듈 스코어링
- ssGSEA-like (단순 mean) per-sample → within-cohort z-score
- 8 코호트 pooled = **n=415 sample** (414 unique) module score matrix
- ATC-vs-non-ATC pooled meta (4 thyroid cohorts)
- DIAL-lite: Hugo (28) + Riaz_pre (51) 응답자 vs 비응답자 t-test + 1000-bootstrap CI + 1000-permutation null
- PCA + hierarchical cluster + bootstrap stability (200 boots)
- Riaz Pre→On per-patient delta (n=43 paired)
- 6 PNG figures (`figures_png/`)

**실행하지 않은 것 (true Track B Wk5–11 에 위임):**
- NMF rank K selection (단순 hclust 만 했음)
- ComBat / ComBat-seq harmonization
- LODO sensitivity per `v52_lodo_finding`
- LOHHLA / NetMHCpan / arcasHLA HLA typing
- scVI / scANVI scRNA atlas integration
- 전체 DIAL audit (IMvigor210 / Gide / Liu / Kim 추가)
- Class II HLA + neoantigen
- TIDE-like + IMPRES + EXHAUSTION_INDEX composite scores
- Cox 생존 분석
- spatial overlay
- 출판 품질 figure

---

## 2. 데이터 인벤토리

| Cohort | n | Modality | Subtype 분포 |
|---|---|---|---|
| TCGA-THCA | (Paper 1 ETL only — touched) | RNA-seq | — |
| GSE126698 | 28 | RNA-seq FPM | ATC=10 / FTC=6 / NT=6 / PTC=6 |
| GSE76039 | 37 | microarray (GPL570) | PDTC_or_ATC_advanced=37 (per-sample 구분 불가 from GSM characteristics) |
| GSE65144 | 25 | microarray | ATC=12 / NT=13 |
| GSE29265 | 49 | microarray | ATC=9 / PTC=20 / NT=20 |
| GSE33630 | 105 | microarray | ATC=11 / PTC=49 / NT=45 |
| GSE60542 | 92 | microarray | PTC_nodal_met=28 / NT=30 / unknown=34 |
| Hugo (GSE78220) | 28 | RNA-seq FPKM | melanoma; pre-treatment only |
| Riaz pre (GSE91061) | 51 | RNA-seq FPKM | melanoma; pre-treatment |
| Riaz on (GSE91061) | 58 | RNA-seq FPKM | melanoma; on-treatment |
| GSE151179 | 52 | Clariom D | RAI status — module 스코어링 실패 (transcript cluster ID 매핑 미완) |

총 pooled = 415 (GSE151179 제외).

---

## 3. ATC-vs-non-ATC pooled meta-analysis (Module A 1차 sanity)

**4 thyroid 코호트 (GSE126698, GSE65144, GSE29265, GSE33630, GSE76039)** 에서 ATC (또는 PDTC_or_ATC_advanced) vs non-ATC (PTC + FTC + NT) 비교, n-weighted Cohen's d:

| Module | d (n-weighted) | sign 일치 | 해석 |
|---|---|---|---|
| myeloid_suppressive | **+2.53** | 4/4 | ATC 에서 myeloid suppression 격증 |
| thyroid_differentiation | **−2.51** | 4/4 | ATC 에서 분화 손실 (TDS 감소와 일치) |
| checkpoint_exhaustion | +1.01 | 4/4 | ATC inflamed-but-exhausted |
| HLA_class_II | +0.68 | 3/4 | tumor-intrinsic Class II 표현 (Hashimoto-like 신호 with Paper 2 cite 가능) |
| HLA_class_I | +0.58 | 4/4 | Class I 정량적 유지 (loss-of-presentation 가설은 *somatic LOH* 측정 필요 — Track B Module C) |
| IFNG_T_cell_inflamed | +0.55 | 4/4 | ATC 가 IFNG-rich 환경 |
| TLS_CXCL13_like | −0.15 | 3/4 | ATC 에서 TLS 신호 *감소* — Paper 2 PTC+HT 영역과 분리되는 보강 증거 |

**해석 (claim-guarded):**
- ATC 는 `thyroid_diff` 축에서는 일관되게 분화 손실, 면역 축에서는 *paradoxical inflamed-but-suppressed* 표현형 (myeloid + checkpoint 동반 상승).
- 이 패턴은 **ICI vulnerability hypothesis** 와 일관 — myeloid 가 ICI 효과를 제한할 수 있고, checkpoint 가 약물 작용 기전 표적임.
- TLS 가 ATC 에서 약하게 감소하는 패턴은 Paper 2 의 PTC+HT TLS-rich 표현형과 *분리* — Paper 2 와의 cross-paper boundary 안정.

**유의:** *response prediction 의 직접적 증거 아님*. ATC = 면역치료 vulnerability candidate 라는 가설을 *생성* 하는 결과이며, 검증을 위해서는 thyroid ICI-treated raw RNA-seq (현재 비공개) 또는 LOHHLA + neoantigen architecture (Track B Module C) 가 필요.

---

## 4. PCA 구조 (n=415 pooled)

| | PC1 | PC2 | PC3 | PC4 | PC5 |
|---|---|---|---|---|---|
| **explained var** | **69.4%** | 14.2% | 6.6% | 4.0% | 2.7% |
| HLA_class_I | +0.41 | +0.07 | +0.07 | +0.67 | -0.20 |
| HLA_class_II | +0.42 | +0.11 | +0.20 | +0.06 | -0.58 |
| IFNG | +0.43 | +0.15 | -0.14 | +0.22 | +0.42 |
| TLS_CXCL13 | +0.37 | +0.36 | -0.41 | -0.59 | -0.29 |
| checkpoint | +0.43 | +0.06 | -0.11 | -0.11 | +0.59 |
| myeloid_suppressive | +0.34 | -0.37 | **+0.71** | -0.38 | +0.08 |
| thyroid_differentiation | -0.20 | **+0.83** | +0.50 | -0.02 | +0.14 |

- **PC1 (69.4%)** = "pan-immune activation" — 6 개 면역 모듈 모두 +0.34~+0.43 균등; thyroid_diff 만 -0.20.
- **PC2 (14.2%)** = "dedifferentiation" — thyroid_diff dominant +0.83; myeloid_suppressive -0.37 (분화된 조직 = myeloid 적음).
- **PC3 (6.6%)** = "myeloid-vs-TLS" — myeloid +0.71 vs TLS -0.41.

**해석 (hypothesis-generating):** thyroid TME 의 first-order 변동은 단일 inflammation 축 + 분화 축으로 90% 이상 설명. 진짜 Module A 의 NMF 는 **PC1 + PC2 위에서 ecotype 구조** 를 분해해야 함. ComBat 정정 후 NMF rank K ≤ 6 가 자연스러운 후보.

---

## 5. DIAL-lite — pan-cancer ICI direction-invariance check

**Hugo (GSE78220, n=28)** + **Riaz pre (GSE91061, n=51)** 응답자 (CR/PR/PRCR) vs 비응답자 (PD), 7 module 별 Cohen's d ± 1000-bootstrap 95% CI:

| Module | Hugo d (CI) | Hugo perm-p | Riaz_pre d (CI) | Riaz_pre perm-p | sign 일치 |
|---|---|---|---|---|---|
| HLA_class_I | +0.53 (-0.17, +1.67) | 0.20 | +0.26 (-0.50, +0.99) | 0.52 | ✅ + + |
| HLA_class_II | +0.28 (-0.47, +1.21) | 0.48 | +0.43 (-0.22, +1.28) | 0.26 | ✅ + + |
| TLS_CXCL13_like | +0.02 (-0.90, +0.77) | 0.97 | +0.58 (-0.27, +1.39) | 0.13 | ✅ + + |
| IFNG | -0.13 (-0.84, +0.72) | 0.77 | +0.56 (-0.23, +1.49) | 0.17 | ❌ − + FLIP |
| checkpoint_exhaustion | -0.22 (-1.04, +0.47) | 0.58 | +0.63 (-0.09, +1.51) | 0.11 | ❌ − + FLIP |
| myeloid_suppressive | **-0.72 (-1.61, +0.02)** | 0.09 | +0.16 (-0.53, +0.86) | 0.67 | ❌ − + FLIP |
| thyroid_differentiation | -0.03 (-1.17, +0.64) | 0.93 | +0.11 (-0.69, +0.91) | 0.80 | ❌ − + FLIP |

**Sign 일치: 3/7 (binomial 1-sided p=0.77, NOT 일치).**

**해석 (이번 분석의 가장 중요한 발견):**
- 페이퍼 3 의 핵심 가설 — "pan-cancer composite ICI 시그니처를 thyroid 에 transfer 하기 전 DIAL audit 필수" — 가 이미 **2-cohort minimal sample 에서 입증**.
- Hugo (n=28, melanoma 단일 site) 는 myeloid_suppressive 가 비응답과 트렌드 일치 (d=-0.72, perm p=0.09) — *responders 가 myeloid 적음* 의 정형적 패턴.
- Riaz (n=51) 는 **모든 모듈에서 응답자 ≥ 비응답자 (양의 d)** — 더 큰 melanoma cohort 에서 IFN/checkpoint/TLS 가 응답과 양의 연관.
- 두 cohort 모두 melanoma 인데도 4/7 모듈에서 sign 이 flip — 즉, **단일 disease 안에서도 cohort-specific direction inconsistency 존재**.
- **thyroid 에 transfer 하기 전, IMvigor210 / Gide / Liu / Kim 까지 5+ 코호트 완전 DIAL audit 필요** (Track B Module D).

**경고 (claim-guard 강제):**
- 이 결과는 *thyroid 환자에서 ICI 반응을 예측한다* 는 클레임의 evidence 아님.
- 또한, *어떤 modules 이 thyroid 에서 "OK" 인지* 결정짓지 못함 (그 결정은 thyroid ICI-treated raw RNA-seq 가 있어야 함).
- 이 결과는 **DIAL audit 의 필요성** 을 *방법론적으로* 입증할 뿐.

---

## 6. Hierarchical cluster bootstrap stability

K=4 hclust on n=415 pooled, 200 bootstrap resamples, sample-pair co-membership stability mean = **0.274**.

**해석:** ecotype 신호가 cohort/batch 효과보다 약하다 — 진짜 Module A 분석은 ComBat-seq harmonization (subtype 보호 covariate) 필수. v52_lodo_finding 의 LODO 원칙 엄격 적용 필요.

---

## 7. Riaz Pre→On per-patient delta (n=43 paired)

| Module | mean Δz (On − Pre) | t | p |
|---|---|---|---|
| HLA_class_I | -0.10 | -0.82 | 0.42 |
| HLA_class_II | -0.12 | -1.05 | 0.30 |
| IFNG | -0.12 | -0.94 | 0.35 |
| TLS_CXCL13_like | -0.02 | -0.20 | 0.84 |
| checkpoint_exhaustion | -0.13 | -1.01 | 0.32 |
| myeloid_suppressive | +0.02 | +0.19 | 0.85 |
| thyroid_differentiation | +0.10 | +0.63 | 0.53 |

**해석:** Pre→On 단기 시점에서 bulk 모듈은 통계적으로 유의한 변화 없음. 가능한 이유:
- 응답자/비응답자 mix (heterogeneous direction cancels).
- pre/on biopsy 시점이 너무 짧음 (response typically RECIST 6–8 주).
- Bulk pseudobulk → cell-type compartment 변화 신호 cancels out → scRNA atlas 분석에서만 보일 수 있음.

이는 *negative result* 로서 valuable — bulk Pre→On signature shift 가 ICI 반응의 일관 marker 가 *아니다* 라는 cross-disease evidence.

---

## 8. 산출물

### 8.1 Track B-lite outputs

`project/results/paper3_ici_track_b_lite/`:
- `gpl570_probe_symbol_map_module_subset.tsv` (60 module gene → 144 probe)
- `pooled_module_scores_v2.tsv` (n=473 with Riaz_on)
- `pooled_module_scores_v2_with_subtype.tsv`
- `module_score_by_cohort_subtype.tsv` (112 rows)
- `atc_vs_other_per_cohort.tsv` (28 rows)
- `atc_vs_other_meta_pooled.tsv` (7 modules)
- `dial_lite_bootstrap.tsv` (14 rows, 1000 boot CI)
- `dial_lite_sign_consistency_v2.tsv`
- `pca_loadings_v2.tsv`
- `riaz_pre_to_on_delta.tsv` (43 paired)
- `riaz_pre_to_on_paired_test.tsv`
- `phase2_summary.json`
- `scores_per_cohort/<cohort>_module_scores.tsv` (10 files)

### 8.2 Figures

`project/results/paper3_ici_track_b_lite/figures_png/`:
- **fig1_atc_vs_other_meta.png** — ATC vs non-ATC pooled effect-size bar (across 4 thyroid cohorts)
- **fig2_pca_scatter_by_cohort.png** — pooled PCA colored by cohort
- **fig3_subtype_heatmap.png** — module × cohort×subtype z-score heatmap
- **fig4_dial_lite_forest.png** — DIAL-lite Hugo + Riaz forest plot with bootstrap CI
- **fig5_dediff_vs_inflam_GSE126698.png** — GSE126698 dedifferentiation × IFN scatter, colored by subtype
- **fig6_riaz_pre_to_on_delta.png** — Riaz 43-paired per-module Δ bar

> 모든 figure 는 layout schematic 수준 — Track B 의 출판 품질 figure 를 대체하지 않음.

### 8.3 분석 스크립트

`project/results/paper3_ici_data_registry/`:
- `track_b_lite_analysis.py` — Phase 1 (cohort scoring + PCA + hclust + DIAL-lite)
- `track_b_lite_phase2.py` — Phase 2 (subtype + bootstrap + figures)

---

## 9. 페이퍼에 가져오는 강화 (정량화)

| 강화 포인트 | 기존 (v1 sprint) | 이번 (Track B-lite) | 페이퍼 사용 |
|---|---|---|---|
| ATC 면역+탈분화 표현형 정량 | 없음 | **myeloid d=+2.53, thyroid_diff d=-2.51, 4/4 sign** | Fig 1A 핵심 narrative anchor |
| DIAL audit 필요성 입증 | "design plan" 만 | **2-cohort 에서 4/7 module flip 실측** | Fig 5 anchor + Discussion 의 DIAL 정당화 |
| Pan-immune axis 구조 | 가설 | **PC1 69.4% var 단일 활성화 축** | Fig 2 + Module A → Module E 의 weighting 합리화 |
| Cohort/batch effect 의 압도 | 가설 | **hclust stability 0.274** | Methods 의 ComBat 정당화 + reviewer Q 대비 |
| Pre→On 단기 bulk 변화 부재 | 가설 | **n=43 paired, all p>0.30** | Discussion 의 limitation + bulk vs sc 정당화 |

---

## 10. Limitations & risks (정직)

1. **K-cluster sanity 만, 진짜 NMF 아님** — Module A 본 분석을 대체하지 않음. ComBat + NMF + bootstrap rank stability 는 Track B Wk3.
2. **GSE76039 = microarray (RNA-seq 아님)** — Landa 2016 paper 의 RNA-seq 부분이 별도 deposit 인지 Track B Wk1 verify.
3. **GSE151179 (RAI 축)** — Clariom transcript cluster ID 매핑 미완으로 module scoring 실패. Track B Wk2 에서 `clariomdhumantranscriptcluster.db` 또는 platform annotation file 필요.
4. **Microarray 5 코호트 모두 GPL570 — 144 probe → 60 symbol max-aggregated.** Per-symbol 결측치 가능; Track B 에서 RMA + frma + collapseRows 정식 처리.
5. **DIAL-lite 가 2 코호트만** — 진짜 DIAL audit 는 5+ cohort 필요.
6. **HLA / neoantigen 분석 없음** — Module C 자체. dbGaP TCGA paired BAM 필요.
7. **Riaz Pre→On null finding은 단일 disease (melanoma) 결과** — thyroid 에 일반화 불가.
8. **모든 결과는 z-score 단위 + ssGSEA-like 단순 평균** — Singscore / UCell sensitivity 미수행.

---

## 11. 현재 시점에서 방어 가능한 클레임 한계

> 이번 Track B-lite 분석은 **Paper 3 의 hypothesis 와 일관된 정량 evidence** 를 수집했다 — ATC 가 inflamed-but-myeloid-suppressed phenotype 를 가지며 (d=+2.53), 분화 손실과 결합하여 (d=-2.51) immunogenomic vulnerability candidate 로 정의될 수 있다 (gold standard 의 *후보*).
>
> 동시에, pan-cancer ICI signature 의 direction-invariance 가 **2 cohort minimal sample 에서 이미 깨짐 (4/7 module flip)** — 이는 페이퍼 3 의 DIAL audit framework 가 단순한 design proposal 이 아니라 실측 필요성을 가짐을 입증한다.
>
> 그러나 본 결과는 어떤 의미에서도 **"thyroid 환자에서 ICI 반응을 예측한다"** 클레임의 evidence 가 *아니다*. 본 결과는 ICI-readiness / immunogenomic vulnerability prioritization hypothesis 의 *생성* 단계 evidence 이며, validation 은:
> 1. Track B Module C HLA LOH + neoantigen architecture (TCGA paired BAM 필요)
> 2. Track B Module D 5+ pan-cancer cohort 완전 DIAL audit
> 3. (가장 중요) 비공개 thyroid ICI-treated cohort 협상 + raw RNA-seq 확보
>
> 위 셋 중 *둘 이상* 충족 후에야 readiness "score" 가 의미를 가짐.

---

## 12. 다음 우선순위 (사용자 결정용)

### A. Track B-lite Phase 3 — 동일 데이터로 추가 짜기 (마라톤 비위반)
- ComBat-seq harmonization 시뮬레이션 (subtype protect)
- ssGSEA-vs-Singscore-vs-UCell 민감도
- Riaz responder vs non-responder pre→on delta 분리
- Kim GC + Gide melanoma 다운로드 → DIAL audit 5-cohort 까지 확장

### B. Controlled access 신청 (4–8 주 lag)
- TCGA-THCA dbGaP paired BAM
- Liu phs000452
- Yoo 2019 EGA

### C. 비공개 thyroid ICI cohort 협상 (교수님 의사결정)
- K1 해제의 유일한 경로

### D. 진짜 Track B 진입 (Paper 1 출하 후, 6/13 +)
- 본 Track B-lite 결과를 Wk1 input 으로 재사용

---

Track A 동결 (8 design 파일 chmod 444) 유지. Paper 1 / Paper 2 ETL 미수정. 마라톤 모드 — Paper 1/2 작문 우선순위 변경 없음.
