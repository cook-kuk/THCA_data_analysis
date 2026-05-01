---
title: "Paper 2 brief — cohort raw data access guide"
date: 2026-05-01
purpose: "5 cohorts (TCGA / GSE286332 / GSE213647 / K2 PRJEB11591 / Chu 2018) 의 accession + download method + processing script + 데이터 size. Advisor / reviewer raw 자료 reproducibility 검증용."
scope: "Paper 2 (Autoimmune-overlap PTC, DM2 axis) 5-Pillar 의존 5 cohorts"
---

# Paper 2 brief — cohort raw data access guide

본 brief 의 21 figures + 5 tables 가 의존하는 5 cohorts 의 raw 자료 access 방법.

## 1. Cohort 요약

| # | Cohort | n | 역할 | Access | 상태 |
|---|---|---|---|---|---|
| 1 | TCGA-THCA | 500 | Discovery (Pillar III, V, mediation) | GDC + cBioPortal | ✅ 처리 완료 |
| 2 | GSE286332 | 18 | Pillar II in-cohort, Pillar IV BCR | GEO | ✅ 처리 완료 |
| 3 | GSE213647 | 632 | Pillar III Korean replication | GEO | ✅ 처리 완료 |
| 4 | K2 / PRJEB11591 | 235 typeable | Pillar I HLA + arcasHLA | ENA | ✅ 처리 완료 |
| 5 | Chu 2018 Han Chinese GD | 1,468 GD / 1,490 ctrl | Pillar I forest meta reference | Published table | ✅ 통합 완료 |

## 2. Cohort detail

### 2.1 TCGA-THCA (n=500)

| 항목 | 값 |
|---|---|
| Source | The Cancer Genome Atlas (TCGA) Pan-Cancer |
| Project | TCGA-THCA (Thyroid Carcinoma) |
| Access | GDC Data Portal (https://portal.gdc.cancer.gov/projects/TCGA-THCA) + cBioPortal `thca_tcga` study |
| Data types used | bulk RNA-seq (HTSeq counts + TPM), DNA-seq mutations, clinical |
| Size | RNA-seq ~10 GB, mutations ~100 MB |
| Download | `gdc-client download` with manifest, OR cBioPortal API |
| Processing | `notebooks_or_scripts/v17_*tcga*.py` (multiple) |
| Brief usage | Pillars III (Hashimoto-like signature transfer), V (DM1 sub-A/sub-B), § 3.6 mediation (within-PTC) |
| Metadata | `metadata/sample_master_v3.tsv` (sample manifest with TCGA-THCA samples) |

### 2.2 GSE286332 (Korean Dongguk Univ Lim 2025, n=18)

| 항목 | 값 |
|---|---|
| Source | NCBI GEO + SRA |
| GEO accession | GSE286332 |
| BioProject | PRJNA1208932 |
| Submitter | Lim et al. 2025 Dongguk University |
| Sequencing | Macrogen Seoul, Illumina NovaSeq X |
| Sample composition | 9 PTC + 9 PTC+HT |
| Access | `prefetch SRR12345 && fasterq-dump` (SRA Toolkit) per sample |
| Size | Raw FASTQ ~50 GB total |
| Processing | `notebooks_or_scripts/pillar2_GSE286332/ptc_vs_ptcht_DEG_GSEA.py` (symlink → `v17_P3_GSE286332_ptc_vs_ptcht.py`) |
| Brief usage | **Pillar II** (Fig 3-7 in-cohort mechanism), **Pillar IV** (Fig 12-13 BCR clonal+TLS), § 3.6 mediation primary cohort |
| Output | `results/p3_gse286332/` (DEG + GSEA + 8-gene + HLA module + dm12 predictions) |
| Output also | `results/d5p6_bcr_repertoire/` (BCR), `results/d3p5_pdm1_gradient/` (mediation) |

### 2.3 GSE213647 (Lee 2024 Korean replication, n=632)

| 항목 | 값 |
|---|---|
| Source | NCBI GEO |
| GEO accession | GSE213647 |
| Publication | Lee et al. 2024 |
| Sequencing | Korean PTC bulk RNA-seq |
| Sample composition | 632 PTC samples (large Korean validation cohort) |
| Access | GEO matrix file (already-processed expression) |
| Size | Processed matrix ~500 MB |
| Processing | `notebooks_or_scripts/v17_ULTIMATE_U1C_gse213647.py` |
| Brief usage | **Pillar III Korean arm** (Fig 8 prevalence, Fig 15 sub-B-like rate) |
| Output | `results/d8b_korean_replication/`, `results/d8c_dm1_subB_x_K2_NBNR/` |
| Metadata | `metadata/GSE213647_supplementary_data.xlsx` (sample-level metadata) |

### 2.4 K2 / PRJEB11591 (Yoo 2016 SNU-GMI, n=235 typeable)

| 항목 | 값 |
|---|---|
| Source | EBI ENA |
| ENA accession | PRJEB11591 |
| Publication | Yoo et al. 2016 SNU-GMI Korean PTC |
| Sample composition | 260 manifest samples, 235 valid 4-digit arcasHLA calls |
| Access | `enaBrowserTools` or aspera per run accession |
| Size | Raw FASTQ ~100 GB |
| Processing | arcasHLA 0.6.0 RNA-seq HLA imputation |
| Brief usage | **Pillar I** (n=874 Korean PTC pool 의 K2 arm), § 3.5 K2 NBNR boundary |
| Output | `results/d4p1_panasian_meta/korean_PTC_pool_n908.tsv` (per-sample 4-digit alleles, NOT deprecated) |
| Metadata | `metadata/v3_fusion_anchor_prjeb11591.tsv` |
| Note | K2 ≠ Bundang (memory `v17_K2_vs_bundang_distinction`). Bundang SNUH 은 outreach-stage 별도 cohort. |

### 2.5 Chu X et al. 2018 (Han Chinese GD, n=1,468 GD / 1,490 ctrl)

| 항목 | 값 |
|---|---|
| Citation | Chu X, Pan C-M, Zhao S-X, et al. *J Med Genet* 2018;55(10):685–692 |
| DOI | 10.1136/jmedgenet-2017-105146 |
| PMC | PMC6161647 |
| Sample composition | 1,468 Graves' disease + 1,490 controls (Han Chinese) |
| Methodology | 4-digit SNP2HLA imputation against Pan-Asian reference panel |
| Access | Open published summary statistics (paper Table 2 + supplementary) |
| Brief usage | **Pillar I** reference arm (forest meta vs Korean PTC pool n=874) |
| Integrated | `results/p2_pillar1_forest/chu2018_allele_summary.tsv` (8 focus alleles published values) |
| ⚠️ Citation history | Prior memory + reports listed as "Chen 2018" — INCORRECT. Verified via PMC 6161647 web fetch 2026-05-03. All references corrected. |
| Note | Published table only. No raw genotype download required (summary stats sufficient for forest meta). |

## 3. Reproducibility — End-to-end pipeline

### 3.1 Pillar I (HLA forest meta) reproduce path

```bash
# 1. Korean PTC pool n=874 (already processed)
cat results/d4p1_panasian_meta/korean_PTC_pool_n908.tsv
# 2. Chu 2018 published values (already integrated)
cat results/p2_pillar1_forest/chu2018_allele_summary.tsv
# 3. Re-run forest meta
python notebooks_or_scripts/pillar1_HLA/forest_meta.py
# Outputs → results/p2_pillar1_forest/
```

### 3.2 Pillar II (GSE286332 in-cohort) reproduce path

```bash
# 1. Download GSE286332 raw FASTQ (~50 GB)
prefetch -v SRR... && fasterq-dump SRR...
# 2. Salmon quant (already processed; outputs in cache)
# 3. PyDESeq2 + GSEApy
python notebooks_or_scripts/pillar2_GSE286332/ptc_vs_ptcht_DEG_GSEA.py
# Outputs → results/p3_gse286332/
```

### 3.3 Pillar III (cross-cohort generalization) reproduce path

```bash
# 1. GSE286332 PTC+HT-vs-PTC top DEG signature → TCGA transfer
python notebooks_or_scripts/pillar5_autoimmune/tcga_hashimoto_signature.py
# Outputs → results/d4p2_tcga_hashimoto_signature/
# 2. Korean GSE213647 replication
python notebooks_or_scripts/v17_ULTIMATE_U1C_gse213647.py
# (별도 step, signature transfer in d8b_korean_replication)
```

### 3.4 Pillar IV (BCR clonal + TLS) reproduce path

```bash
# 1. GSE286332 raw FASTQ → BCR repertoire (IgBLAST etc.)
python notebooks_or_scripts/pillar5_autoimmune/bcr_repertoire.py
# Outputs → results/d5p6_bcr_repertoire/
```

### 3.5 § 3.6 Mediation reproduce path

```bash
# Bootstrap 5000-iter mediation analysis
python notebooks_or_scripts/pillar5_autoimmune/pdm1_mediation.py
# Outputs → results/d3p5_pdm1_gradient/
```

## 4. Data size / disk requirement

| Cohort | Raw size | Processed | Required for brief reproduce |
|---|---|---|---|
| TCGA-THCA | ~10 GB | ~500 MB | Already processed, ~500 MB enough |
| GSE286332 | ~50 GB | ~50 MB | Already processed, raw not needed for brief |
| GSE213647 | ~500 MB matrix | ~10 MB | Already processed |
| K2 / PRJEB11591 | ~100 GB | ~5 MB (HLA TSV) | Already processed |
| Chu 2018 | < 1 KB summary | <1 KB | Published table |
| **TOTAL** | **~160 GB raw** | **~570 MB processed** | **Brief reproduce: ~600 MB** |

→ **Reviewer 가 brief 의 모든 figure 를 reproduce 하기 위한 필요 자료 ~600 MB** (raw FASTQ 다운로드 불필요).

## 5. Cohort 사용 권한 / IRB

| Cohort | 사용 권한 |
|---|---|
| TCGA | Open access (controlled tier 별도 신청, 본 분석은 open tier 만 사용) |
| GSE286332 | GEO open submission (Lim 2025) |
| GSE213647 | GEO open submission (Lee 2024) |
| K2 / PRJEB11591 | Open access ENA (Yoo 2016 published) |
| Chu 2018 | Published table (open access summary stats, individual genotypes 별도 신청 필요) |

→ **모든 brief data 는 open access 또는 published**. 개별 환자 IRB 추가 신청 불필요. Bundang SNUH cohort (outreach 단계) 만 별도 IRB.

## 6. Future cohort additions (Phase 1 / scenario A)

| Cohort | 상태 | 계획 |
|---|---|---|
| Bundang SNUH Graves' / PTC+HT | outreach 단계 | n>50 GD + n>30 PTC+HT 가능 시 Pillar I + § 3.5 cookHLA SNP-based 보강 |
| Taiwan CMUH GD (n=2,998 HIBAG) | controlled access | Gap 3 plan 이었으나 Chu 2018 통합으로 deferred — Phase 1 옵션 |
| Japanese Okada 2015 | published summary | Gap 3 plan 이었으나 deferred — Phase 1 옵션 |

---

*Generated 2026-05-01. Reviewer / advisor 가 raw 자료 access 검증 시 본 guide 참고. Brief data 사용 권한 모두 open access — IRB 추가 신청 불필요.*
