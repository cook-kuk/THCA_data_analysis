# Paper 3 ICI — Data Access Blockers

**Date:** 2026-05-06
**Author:** Seungho Cook
**Purpose:** 분석 가능한 데이터와 차단된 데이터를 명시적으로 구분. registry v2 의 `access_type` + `blocker` 컬럼의 종합 view.

---

## 1. 즉시 사용 가능 (public, no application needed)

| Dataset | Modality | n | 이번 sprint 다운로드 |
|---|---|---|---|
| GSE76039 | bulk_microarray (GPL570) | 37 | ✅ series_matrix |
| GSE65144 | microarray | 25 | ✅ series_matrix |
| GSE29265 | microarray | 49 | ✅ series_matrix |
| GSE33630 | microarray | 105 | ✅ series_matrix |
| GSE60542 | microarray | 92 | ✅ series_matrix |
| GSE151179 | gene_expression | 52 | ✅ series_matrix |
| GSE126698 | bulk_RNAseq | 28 | ✅ series_matrix + DE table |
| GSE78220 (Hugo) | bulk_RNAseq + clinical | 28 | ✅ series_matrix + FPKM |
| GSE91061 (Riaz) | bulk_RNAseq + clinical + WES | 109 | ✅ series_matrix + FPKM |
| GSE184362 | scRNA | 11 patients / 158K cells | △ series_matrix only (RAW 970 MB defer) |
| GSE193581 | scRNA + bulk cell line | 32 GSM | ✅ series_matrix + celltype_annotation + cell_line bulk |
| GSE232237 | scRNA | 12 GSM | △ series_matrix only |
| GSE191288 | scRNA | 7 GSM | △ series_matrix only |
| GSE148673 | scRNA | 13 GSM (1 ATC) | △ series_matrix only |

---

## 2. Controlled access (application required)

| Dataset | Authority | Estimated lag | Required for |
|---|---|---|---|
| TCGA-THCA paired tumor-normal BAM | dbGaP | 4–8 weeks | LOHHLA Module C |
| Liu melanoma WES + RNA-seq | dbGaP phs000452 | 4–8 weeks | DIAL audit Tier 2 anchor |
| Yoo SK 2019 Korean ATC | EGA (likely) | 4–8 weeks if controlled | Asian dedifferentiation generalization |

---

## 3. Metadata only (raw data 비공개)

이 데이터는 raw expression / WES 가 publicly 비공개이며, **publication supplementary table + IHC 결과만 사용 가능**. **response prediction 클레임의 evidence가 될 수 없음** — citation 만.

| Trial / Paper | NCT | 사용 |
|---|---|---|
| KEYNOTE-158 thyroid | NCT02628067 | clinical evidence cite |
| KEYNOTE-028 thyroid | NCT02054806 | clinical evidence cite |
| NCT03246958 (nivo+ipi aggressive thyroid) | NCT03246958 | clinical evidence combo cite |
| Dierks 2021 lenva+pembro | — | clinical evidence cite |
| Spartalizumab ATC | NCT02404441 | clinical evidence cite |
| Atezolizumab matched TT ATC | NCT03181100 | clinical evidence cite |
| 다양한 lenva+pembro RAI-r DTC | various | clinical evidence cite |

---

## 4. 영원히 비공개일 수 있는 항목

- **thyroid ICI-treated raw RNA-seq** — KEYNOTE-158 / -028 / Spartalizumab / Sehgal / Dierks / Cabanillas 모두 raw 비공개. 이것이 **K1 (response predictor 금지)** 의 단일 근본 원인.
- 비공개 기관 코호트와의 협력 (연세 / 서울대 / 분당 / 삼성) 만이 K1 해제 경로.

---

## 5. Cross-paper boundary 차단

이번 sprint 에서 사용하지 않는 자산 (Paper 1/2 영역):

| Asset | 이유 | 대신 사용한 자산 |
|---|---|---|
| Paper 1 DM1 cluster 분석 결과 | Paper 1 본진 — 재분석 금지 | dark matter (BRAF/RAS-neg) 정의만 read-only 참조 |
| Paper 2 PTC+HT TLS / AICDA / BCR | Paper 2 Pillar 1/2 본진 | comparator group 라벨로만 참조 |
| Paper 4 GD HLA forest | 별도 paper 영역 | population HLA frequency 참조만 |
| GSE250521 spatial | Paper 1/2 supplement 본진 | TLS overlay 시각화만 read-only |
| PRJEB11591 arcasHLA 결과 | v17_arcasHLA_korean_k2 read-only | Korean HLA frequency anchor 만 인용 |

---

## 6. Sprint 권한 한계 (이번 작업에서 *하지 않은* 것)

- Track A 동결 8 design 파일 수정 안 함.
- Paper 1 / Paper 2 ETL 재실행 안 함.
- NMF / HLA typing / NetMHCpan / LOHHLA / DIAL 실행 안 함.
- scRNA atlas integration (scVI / scANVI) 안 함.
- 실데이터 figure rendering 안 함.
- TCGA-THCA paired BAM dbGaP 신청 안 함 (사용자 명시 결정 대기).
- Yoo 2019 EGA 신청 안 함.

---

## 7. 차단 해제를 위한 사용자 결정 옵션

| 옵션 | 결정 결과 | 시간 |
|---|---|---|
| dbGaP TCGA paired-BAM 신청 | LOHHLA 가능 | 4–8 weeks lag |
| dbGaP Liu phs000452 신청 | DIAL Tier 2 추가 | 4–8 weeks lag |
| EGA Yoo 2019 신청 | Asian ATC anchor | 4–8 weeks lag |
| Yu 교수 → 비공개 thyroid ICI 코호트 협상 | **K1 해제** (response prediction 클레임 가능) | months |
| Track B 진입 명시 | full Module A–E 분석 | Paper 1 출하 + Paper 2 closure 후 |

---

Track A 동결 (Paper 1/2 분석) 유지. 마라톤 비위반.
