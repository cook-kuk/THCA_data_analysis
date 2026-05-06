# Paper 3 ICI — Public Data Download Log

**Date:** 2026-05-06
**Author:** Seungho Cook
**Storage:** `project/results/paper3_ici_public_data/<accession>/`

---

## 1. 다운로드 성공 (총 14건, ~133 MB)

| Accession | File | Size | Source |
|---|---|---|---|
| GSE76039 | series_matrix.txt.gz | 7.0 MB | NCBI GEO ftp |
| GSE65144 | series_matrix.txt.gz | 5.7 MB | NCBI GEO ftp |
| GSE29265 | series_matrix.txt.gz | 15 MB | NCBI GEO ftp |
| GSE33630 | series_matrix.txt.gz | 30 MB | NCBI GEO ftp |
| GSE60542 | series_matrix.txt.gz | 27 MB | NCBI GEO ftp |
| GSE151179 | series_matrix.txt.gz | 7.5 MB | NCBI GEO ftp |
| GSE126698 | series_matrix.txt.gz + DE_Thyroid_totalRNA_all.csv.gz | 8.1 MB | NCBI GEO ftp |
| GSE78220 (Hugo) | series_matrix.txt.gz + PatientFPKM.xlsx | 6.9 MB | NCBI GEO ftp |
| GSE91061 (Riaz) | series_matrix.txt.gz + fpkm.csv.gz | 16 MB | NCBI GEO ftp |
| GSE184362 | series_matrix only (RAW.tar 970 MB skipped) | 4.2 KB | NCBI GEO ftp |
| GSE193581 | series_matrix + celltype_annotation + cell_line_bulk_count | ~1.4 MB | NCBI GEO ftp |
| GSE232237 | series_matrix only (scRNA RAW skipped) | 2.1 KB | NCBI GEO ftp |
| GSE191288 | series_matrix only (scRNA RAW skipped) | 3.4 KB | NCBI GEO ftp |
| GSE148673 | series_matrix only (scRNA RAW skipped) | 2.7 KB | NCBI GEO ftp |

**Total:** ~133 MB (잘 "small/moderate" envelope 안).

---

## 2. 다운로드 시간 및 방법

- **사용 도구:** `curl` (병렬 `xargs -P 8`).
- **시도한 가속:** `aria2c` 없음 / `parallel` 없음 — `xargs -P` 병렬화로 충분.
- **첫 라운드 (6개 microarray 직렬 curl):** ~25 sec.
- **두 번째 라운드 (3개 series_matrix + suppl 직렬):** ~10 sec.
- **세 번째 라운드 (5개 scRNA series_matrix + 3개 suppl, 8-way 병렬):** ~5 sec.

**합계 ~40 sec** for 14 files. 추가 가속 옵션:
- `aria2c -x 8` (multi-connection per file) — `apt-get install aria2c` 필요. 이번 sprint에서는 도입 안 함.
- ENA mirrored downloads (PRJEB* accession은 EBI 가 더 빠를 수 있음) — Phase 2.5 candidate.
- TCGAbiolinks / recount3 R 패키지 — TCGA-THCA 의 표준 처리 매트릭스 직접 pull. 이번 sprint에서 도입 안 함 (Paper 1 ETL outputs 재사용 원칙).

---

## 3. 다운로드 실패 / 스킵 (defer to Track B)

| Accession | 이유 | 대안 |
|---|---|---|
| GSE184362 RAW.tar | 970 MB scRNA — sprint envelope 초과 | series_matrix + filelist만 확보; Track B Wk5에서 selective per-sample mtx pull. |
| GSE193581 RAW.tar | scRNA tar — 큰 파일 | celltype_annotation + cell_line bulk count 확보. scRNA 본체는 Track B. |
| GSE232237 / GSE191288 / GSE148673 RAW | scRNA — 큰 파일 | series_matrix만 확보. |
| GSE151180 (miRNA) | priority 3 — protein-coding module 직접 사용 불가 | Track B 시 miRNA-mRNA 결합 분석에서 검토. |
| TCGA-THCA 처리 매트릭스 | Paper 1 ETL 재사용 원칙 | `project/results/00_qc/` read-only. |
| PRJEB23709 (Gide) | ENA — Track A v1 to_verify 유지 | Track B Wk1 에서 ENA UI 검증. |
| phs000452 (Liu) | dbGaP controlled access | Wk1 application 필요 — `paper3_ici_data_access_blockers.md` 참조. |
| PRJEB25780 (Kim GC) | ENA — pan-cancer Tier 3 | Track B 에서 retrieval. |
| Cho NSCLC | accession 미확정 | Track B Wk1 verification. |

---

## 4. Series matrix 검증 결과 (n_samples 실제값)

| Accession | user 추정 (input prompt) | 실제 (series matrix) | 차이 | 영향 |
|---|---|---|---|---|
| TCGA-THCA | ~500 | (read-only Paper 1) | — | OK |
| GSE76039 | 17 PDTC + 20 ATC = 37 | **37** | — | ✅ user 정확 |
| **GSE126698** | **36 ATC + 18 PDTC + 132 PTC + 55 FTC + 124 normal = 365** | **28** (10 ATC + 6 FTC + 6 NT + 6 PTC) | **차이 -337** | ⚠️ user 추정 매우 부정확 — registry 정정 필요 |
| GSE184362 | 158,577 cells / 11 patients | 23 GSM libraries (11 patients confirmed in series_summary, 158,577 cells confirmed) | — | ✅ |
| GSE193581 | 10 ATC + 7 PTC + 6 normal + 9 cell lines = 32 | **32 GSM** | — | ✅ user 정확 |
| **GSE65144** | **13** | **25** | +12 | ⚠️ user 추정 부정확 |
| GSE33630 | "~105" | **105** | — | ✅ |
| GSE29265 | "~49" | **49** | — | ✅ |
| **GSE60542** | "to_verify" | **92** (PTC primary + nodal mets) | + | ✅ verified now |
| GSE151179 | "to_verify" | **52** | + | ✅ verified now |
| GSE78220 (Hugo) | ~28 | **28** | — | ✅ |
| GSE91061 (Riaz) | ~109 | **109** | — | ✅ |

**결론:** GSE126698 은 user 추정 (n=365) 이 다른 데이터셋과 혼동된 듯. 실제 GSE126698 = "IGF2BP1 is the first positive marker for ATC diagnosis" (n=28). 큰 thyroid TME 코호트는 실제로 **TCGA-THCA (~500) + GSE126698 (28) + GSE76039 (37) + GSE60542 (92) + GSE151179 (52) + microarray들 (GSE33630 105 + GSE29265 49 + GSE65144 25)** = bulk **~388** + 추가 microarray. registry v2 정정 필요.

---

## 5. 추가 검증 필요 사항

- **GSE76039 의 실제 데이터 형식** — 이번 다운로드된 series_matrix는 GPL570 microarray 형식 (54675 probes). user (및 Track A v1) 가 RNA-seq 으로 표기했지만, GEO 의 GSE76039 deposit 은 *microarray*. Landa 2016 paper의 RNA-seq 부분이 별도 series 인지 확인 필요. → **Track B Wk1 verification 항목.**
- **GSE126698 plat = GPL15456 (HiSeq 2500)** — RNA-seq deposit 맞음. 추가 confirm 필요.
- **scRNA series matrix는 sample-level metadata만 — 실제 cell-level 매트릭스는 RAW.tar 또는 cellranger output** — Track B Wk5 selective pull 계획.

---

## 6. 다음 다운로드 후보 (Track B 진입 시)

- TCGA-THCA paired BAM (dbGaP)
- Liu phs000452 (dbGaP)
- Yoo 2019 Korean ATC (EGA application 후)
- GSE184362 / GSE193581 RAW.tar (scRNA atlas 본체)
- IMvigor210 R-package (R 환경 필요)
- PRJEB23709 (Gide melanoma)

---

Track A 동결 (Paper 1/2 분석) 유지. 마라톤 비위반.
