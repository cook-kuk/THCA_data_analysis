# Paper 3 ICI — Scoring-Ready Matrix Plan

**Date:** 2026-05-06
**Author:** Seungho Cook
**Status:** **Plan only — no scoring execution in this sprint.** Track B (NMF / DIAL / scoring runs) 은 G1–G6 closure 후에만.
**Claim guard:** 모든 스코어링은 ICI-readiness / immunogenomic vulnerability 안에서만 해석. response prediction 의 입력으로 사용 금지.

---

## 1. 스코어링 모듈 (Phase 4 §1 binding)

| Module | gene_count | source registry |
|---|---|---|
| HLA_class_I | 9 | `paper3_ici_module_gene_list.tsv` |
| HLA_class_II | 8 | 동일 |
| IFNG_T_cell_inflamed | 10 | 동일 |
| TLS_CXCL13_like | 9 | 동일 |
| checkpoint_exhaustion | 8 | 동일 |
| myeloid_suppressive | 8 | 동일 |
| thyroid_differentiation | 8 | 동일 |

총 60 gene. Track A signature_registry §1–§7 의 더 큰 카탈로그와 호환 (subset).

---

## 2. 데이터셋별 스코어링 가능성 (gene coverage 검증 결과 기반)

### 2.1 RNA-seq / 직접 symbol-mappable

| Dataset | n_features | id_type | 검증된 module hit / total | 우선순위 |
|---|---|---|---|---|
| GSE126698 | 57,773 | symbol_or_other | IFN 10/10, checkpoint 8/8, TLS 7/9, myeloid 6/8, thyroid 7/8, HLA-I 2/9, HLA-II 2/8 | 1 (가장 큰 hit, sprint 검증 완료) |
| GSE91061 (Riaz) | 22,187 | numeric_likely_entrez | 0 hit (Entrez ID conversion 필요) | 1 (Track B Wk1 conversion) |
| GSE193581 cell_line | 60,675 | ensembl | 0 hit (ENSG → symbol conversion 필요) | 2 (Track B) |
| GSE78220 (Hugo) | (xlsx not yet parsed) | per_supp | 검증 필요 | 1 (Track B) |
| TCGA-THCA | (Paper 1 ETL) | symbol | 60/60 expected | 1 (read-only Paper 1) |

### 2.2 microarray (probe → symbol mapping 필요)

| Dataset | platform | probes | mapping 방법 |
|---|---|---|---|
| GSE76039 | GPL570 (Affy HG U133 Plus 2) | 54,675 | hgu133plus2.db Bioconductor 패키지 (Track B Wk2) |
| GSE65144 | GPL570 | 54,675 | 동일 |
| GSE29265 | GPL570 | 54,675 | 동일 |
| GSE33630 | GPL570 | 54,675 | 동일 |
| GSE60542 | GPL570 | 54,675 | 동일 |
| GSE151179 | GPL23159 (Clariom D Human) | 27,189 | clariomdhumantranscriptcluster.db (Track B Wk2) |

### 2.3 scRNA — UMI count → pseudobulk → module score

| Dataset | RAW.tar 크기 | sprint 수행 | Track B 수행 |
|---|---|---|---|
| GSE184362 | 970 MB | series_matrix만 | per-sample mtx pull + Seurat / scanpy QC + module UCell |
| GSE193581 | (RAW defer) | celltype_annotation + cell_line bulk | tumor scRNA 통합 + UCell |
| GSE232237 | (RAW defer) | series_matrix만 | per-sample pull + UCell |
| GSE191288 | (RAW defer) | series_matrix만 | per-sample pull + UCell |
| GSE148673 | (RAW defer) | series_matrix만 | ATC subset만 추출 |

---

## 3. 스코어링 알고리즘 (Track A signature_registry §0 binding)

- **Bulk:** ssGSEA (GSVA) primary; Singscore sensitivity.
- **Single-cell:** UCell primary; AUCell sensitivity.
- **정규화:** log2(TPM+1) for RNA-seq; RMA log-intensity for microarray; log-normalized counts for scRNA.
- **Score sign convention:** higher = more of feature.
- **Cohort z-score:** within-cohort z-score 후 cross-cohort 비교.
- **Direction-flip handling:** DIAL audit (Module D) 통과 후에만 thyroid 에 적용.

---

## 4. 스코어링 → ecotype 발견 → integrated readiness score 의 흐름

```
[모든 thyroid TME bulk + microarray + scRNA pseudobulk]
        ↓ ssGSEA (bulk) / UCell (sc)
[7-module per-sample score matrix]
        ↓ within-cohort z-score
[harmonized score matrix N_samples × 7_modules]
        ↓ NMF rank K=3..8 + bootstrap
[bulk ecotype labels E1..EK + per-sample assignment confidence]
        ↓ cross-tab vs subtype / driver / TDS / RAI status
[Module A 결과: Fig 2 patterns]

[parallel: HLA typing arcasHLA / OptiType + LOHHLA on TCGA-THCA paired BAM]
        ↓
[HLA Class I/II calls + LOH calls]

[parallel: MC3 MAF + RNA expression filter + NetMHCpan]
        ↓
[neoantigen Class I/II tables, expression-filtered, self-filtered]

[parallel: pan-cancer ICI cohorts ssGSEA → DIAL audit]
        ↓
[DIAL verdict per signature → DIAL-passing subset]

[Integrated readiness score]
ICI_VULN_DIAL = sum_S∈DIAL_PASS w_S · z_S(thyroid_sample)
              + neo_load_log · HLA_intactness − HLA_LOH penalty
              − M2_TAM/MDSC/Treg suppressive penalty
              + (1 − TDS) dedifferentiation lift

[stratify by subtype × driver × TDS × dark_matter → candidate ICI-readiness subset]
```

→ 모든 출력은 readiness/vulnerability framing. **response prediction 클레임 절대 불가 (K1)**.

---

## 5. 1차 sanity 분석 (Track B Wk2–4)

DIAL audit 전에도 가능한 sanity check (Module A 만):
- **GSE126698 + GSE76039 + GSE60542 + GSE151179** 풀링 → 7-module score → NMF K=3..6 → ecotype 후보 라벨.
- **microarray 5 코호트** cross-platform 검증 (probe → symbol mapping 후).
- **TCGA-THCA aggregate** (Paper 1 ETL 재사용) → ecotype 라벨 transfer test.

이 단계 결과는 **discovery only** — DIAL audit 통과 후에만 readiness score 에 반영.

---

## 6. Risk / mitigation

| Risk | Mitigation |
|---|---|
| GSE126698 n=28 — 통계적 power 부족 | TCGA-THCA aggressive subset 과 풀링; 단독 발견 X |
| Microarray probe-symbol 매핑 손실 | <50% gene coverage 모듈은 해당 데이터셋에서 drop, supp 에 보고 |
| Riaz Entrez → symbol 변환 손실 | mygene.info bulk conversion (Track B Wk1) |
| GSE76039 가 실제 microarray (RNA-seq 아님) | Landa 2016 supp 재확인; 별도 RNA-seq accession 가 있으면 추가 |
| scRNA RAW 970MB pull 비용 | per-sample selective pull + cellranger summary metric pre-screen |

---

## 7. 산출물 (Track B 진입 시)

- `scoring/bulk_score_matrix.parquet` (N_samples × 7_modules)
- `scoring/scrna_pseudobulk_score_matrix.parquet`
- `scoring/microarray_probe2symbol_mapped_score_matrix.parquet`
- `scoring/per_cohort_zscore.parquet`
- `scoring/ecotype_labels_NMF_K_*.tsv`
- `scoring/dial_verdict.tsv`
- `scoring/integrated_readiness_score.tsv`

---

Track A 동결 유지. Track B 진입 명령 (사용자 명시 "Track B 시작") 대기.
