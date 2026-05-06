# Paper 3 ICI — 교수님 보고용 데이터 확장 sprint 결과 (한국어 종합)

**작성일:** 2026-05-06
**저자:** Seungho Cook
**Sprint 권한:** 교수님 관심 표명 → 사용자 명시적 권한 → Track A "no public download" 일부 오버라이드. Track B 의 핵심 분석 (NMF / LOHHLA / NetMHCpan / DIAL) 은 여전히 게이팅 유지.
**Claim guard (binding):** ICI-readiness / immunogenomic vulnerability / candidate immune ecotype / hypothesis-generating prioritization 까지만. **"thyroid 에서 ICI 반응 예측" 절대 불가.**

---

## 1. 한 페이지 요약

> **이번 sprint 로 Paper 3 의 데이터 등록과 metadata 가 v1 (Track A 동결, 8개 design 문서) → v2 로 확장되었다.**
> - 23 개 데이터셋이 단일 마스터 TSV 에 등록되었다.
> - 14 개 공개 GEO series matrix (총 ~133 MB) 가 `project/results/paper3_ici_public_data/` 에 다운로드되어 metadata 추출 + 모듈 gene coverage 검증이 완료되었다.
> - 612 개의 sample-level metadata 가 추출되어 `paper3_ici_public_sample_metadata.tsv` 에 저장되었다.
> - **그러나 페이퍼의 방어 가능한 endpoint 는 여전히 "ICI-readiness / immunogenomic vulnerability prioritization" 이며, "thyroid ICI response prediction" 은 raw thyroid ICI-treated RNA-seq 부재로 인해 데이터상 불가능.**
> - K1 (response predictor 금지) 의 유일한 해제 경로는 **연세 / 서울대 / 분당 / 삼성** 등 비공개 thyroid ICI cohort 와의 협력.

---

## 2. v1 → v2 변경점

| 항목 | v1 (Track A 동결, 5/4) | v2 (이번 sprint, 5/6) |
|---|---|---|
| 데이터셋 등록 형태 | 8개 design 문서 산재 | 단일 마스터 TSV (`paper3_ici_data_registry_v2.tsv`) 23 entry |
| 임상 evidence map | 없음 | `paper3_thyroid_ici_clinical_evidence_map.tsv` 7 entry |
| 공개 thyroid TME registry | design 문서 산재 | `paper3_public_thyroid_tme_dataset_registry.tsv` 14 entry |
| 다운로드 manifest | 없음 | `paper3_ici_public_download_manifest.tsv` (priority + URL + size) |
| sample-level metadata | 없음 | 612 sample TSV (subtype / tissue / RAI / metastasis / treatment / response) |
| dataset QC | 없음 | 14 dataset QC 요약 |
| module gene list | signature_registry 안 산재 | `paper3_ici_module_gene_list.tsv` 7 module × 60 gene |
| 실측 gene coverage | 없음 | `paper3_ici_dataset_gene_coverage.tsv` (RNA-seq sets에서 module hit 검증) |
| 검증된 sample count | 산재 | series matrix header 직접 추출로 정정 (GSE126698 365 → 28 등) |

---

## 3. 발견된 데이터의 분석 가능성

### 3.1 즉시 모듈 스코어링 가능 (Track B 진입 시)

| Dataset | n | 비고 |
|---|---|---|
| TCGA-THCA | ~500 | Paper 1 ETL 재사용 (read-only) |
| GSE126698 | 28 | RNA-seq, IFN 10/10, checkpoint 8/8 module hit 검증 완료 |
| GSE76039 | 37 | microarray (Landa 2016) — Track B 에서 probe→symbol 매핑 |
| GSE60542 | 92 | PTC primary + nodal mets, microarray |
| GSE151179 | 52 | RAI-refractory vs RAI-avid PTC, Clariom D |
| GSE33630 | 105 | PTC + ATC + normal, microarray |
| GSE65144 | 25 | ATC + normal, microarray |
| GSE29265 | 49 | PTC + ATC + normal, microarray |
| **GSE184362 scRNA** | 11 patients / 158K cells / 23 GSM | 1차 atlas 후보 |
| **GSE193581 scRNA + cell lines** | 32 GSM (10 ATC + 7 PTC + 6 normal + 9 ATC cell lines bulk) | ATC transformation trajectory |
| GSE232237 scRNA | 12 GSM (3 normal + 7 PTC + 5 ATC) | atlas 보강 |
| GSE191288 scRNA | 7 GSM (bilateral PTC) | atlas 보강 |
| GSE148673 scRNA | 13 GSM (1 ATC) | ATC fallback |

**분석 합산:** bulk + microarray ~388 sample + scRNA 5 코호트 ~158,000+ cells.

### 3.2 Pan-cancer ICI (DIAL audit 용 — response modeling 은 여기서만)

| Dataset | n | 약물 | 다운로드 여부 |
|---|---|---|---|
| Hugo (GSE78220) | 28 | pembrolizumab | ✅ |
| Riaz (GSE91061) | 109 | nivolumab | ✅ |
| IMvigor210 | 298 | atezolizumab | R-package (Track B) |
| Gide PRJEB23709 | 73 | nivo±ipi/pembro±ipi | Track B |
| Liu phs000452 | 121 | nivo/pembro | dbGaP |
| Kim PRJEB25780 | 45 | pembrolizumab | Track B |

### 3.3 임상 evidence (citation only — 분석 X)

KEYNOTE-158, KEYNOTE-028, NCT03246958 (nivo+ipi), Dierks 2021 (lenva+pembro), Spartalizumab ATC, Atezo+matched TT ATC, lenva+pembro RAI-r DTC.

→ Discussion + Introduction hook + Figure 1A 에서 "thyroid 에서 ICI 임상 시도의 baseline benchmark" 로 인용. **분석 입력 아님.**

---

## 4. 어떤 분석이 어떤 데이터로 가능한가

| 분석 (Module) | 1차 데이터 | 2차 (replication) | 게이팅 |
|---|---|---|---|
| **A. Bulk immune ecotype** | TCGA-THCA + GSE126698 | GSE76039 + GSE60542 + GSE151179 + 4 microarray | G6 (Track B 시작) |
| **B. scRNA atlas + sub-state** | GSE184362 + GSE193581 | GSE232237 + GSE191288 + GSE148673 | G4 (≥3 코호트) + G6 |
| **C. HLA + neoantigen** | TCGA-THCA paired BAM + MC3 MAF | GSE76039 (paired normal 일부) + Yoo 2019 if EGA | G3 (dbGaP) + G6 |
| **D. DIAL audit** | IMvigor210 + Hugo + Riaz | Gide + Liu + Kim GC | G5 (≥3 ICI 코호트) + G6 |
| **E. Integrated readiness score** | Module A–D 통합 | scRNA pseudobulk projection | G3+G4+G5+G6 |

---

## 5. 페이퍼에 가져오는 강화 포인트 (정량화)

1. **GSE126698 (n=28, 진짜 RNA-seq)** — 검증된 모듈 hit 매트릭스 확보. ecotype discovery sanity check 가능.
2. **GSE60542 (n=92 PTC primary + nodal mets)** — Track A v1 에 없던 RAI-refractory adjacent 정보 추가.
3. **GSE151179 (n=52 RAI-r vs RAI-avid)** — aggressive 표현형 proxy 축 신규 확보.
4. **GSE184362 + GSE193581 + GSE232237** scRNA — 11 + 10+7+5 ATC/PTC 환자 + 158K cells. Track A v1 의 to_verify 상태 → 검증 완료.
5. **임상 evidence 7 trial** — Discussion 1.1 hook 에서 "thyroid ICI 시도의 ORR 9–20% 범위" benchmark 로 인용 가능. 페이퍼의 임상 relevance 강화.
6. **Hugo + Riaz processed FPKM 확보** — DIAL audit 의 Tier 1 anchor 2 개 sprint 안에서 즉시 검증 가능 (Track B Wk1).

---

## 6. 변경되지 않은 한계 (정직한 진술)

> Paper 3 가 방어할 수 있는 endpoint 는 sprint v2 이후에도 여전히 **"ICI-readiness / immunogenomic vulnerability prioritization"** 이다. **"thyroid 환자의 ICI 반응 예측"** 은 raw thyroid ICI-treated RNA-seq 의 부재로 인해 데이터상 불가능하다.

이 한계는 단 하나의 이유로 결정된다:
- KEYNOTE-158, KEYNOTE-028, Spartalizumab, Sehgal nivo+ipi, Dierks lenva+pembro, Cabanillas atezo+TT — 모두 raw expression public 비공개.

해제 경로는 단 하나:
- 연세 / 서울대 / 분당 / 삼성 등 한국 기관에서 thyroid ICI-treated 환자의 RNA-seq (또는 처리 매트릭스라도) 을 비공개 협력으로 확보.

---

## 7. 다음 분석 우선순위 (사용자/교수님 결정용)

### A. 마라톤 모드 보존 + Track A 추가 강화 (2026-06-13 까지 유지)
1. 사전등록 (preregistration) 초안 — DIAL verdict lock 시점, weighted scoring 옵션, 1차/민감도 결과 정의.
2. K3 (DIAL all-fail) negative-result 시나리오 사전 abstract.
3. claim guard forbidden words 리스트 강화.
4. cross-paper boundary 매트릭스 1-page mockup.

### B. controlled access 신청 sprint (마라톤 비위반)
1. dbGaP TCGA-THCA paired-BAM 신청 (LOHHLA 위해).
2. dbGaP Liu phs000452 신청.
3. EGA Yoo 2019 상태 확인 + 신청.
4. 모두 4–8 주 lag — Paper 1 출하 후 Track B 시작 시점에 맞춤.

### C. 비공개 thyroid ICI 코호트 협상 (교수님 의사결정 핵심)
- 연세 / 서울대 / 분당 / 삼성 어느 곳이든 thyroid ICI-treated 환자 RNA-seq 1 코호트.
- K1 해제 → 페이퍼 framing 격상 ("vulnerability" → "validated readiness").
- Yu 교수님과의 사전 합의 필요.

### D. Track B 진입 (Paper 1 출하 후 자연 트리거)
- 6/13 Paper 1 bioRxiv + Paper 2 closure 후.
- 사용자 명시 "Track B 시작" 명령 필요.
- 이번 sprint 의 산출물이 Track B 12-주 plan 의 Wk1 데이터 acquisition 의 80% 이상을 미리 처리한 셈.

---

## 8. 산출물 인벤토리

### 8.1 Registry (Phase 1)
- `project/reports/paper3_ici/paper3_ici_data_registry_v2.tsv` (23 entry 마스터)
- `project/reports/paper3_ici/paper3_thyroid_ici_clinical_evidence_map.tsv` (7 trial)
- `project/reports/paper3_ici/paper3_public_thyroid_tme_dataset_registry.tsv` (14 entry)
- `project/reports/paper3_ici/paper3_ici_data_collection_report.md`
- `project/reports/paper3_ici/paper3_ici_professor_data_summary_kr.md` (전 sprint v1)

### 8.2 Download (Phase 2)
- `project/reports/paper3_ici/paper3_ici_public_download_manifest.tsv`
- `project/reports/paper3_ici/paper3_ici_public_download_log.md`
- `project/reports/paper3_ici/paper3_ici_data_access_blockers.md`
- `project/results/paper3_ici_public_data/` (실제 다운로드 자산 ~133 MB)

### 8.3 Metadata (Phase 3)
- `project/results/paper3_ici_data_registry/paper3_ici_public_sample_metadata.tsv` (612 sample)
- `project/results/paper3_ici_data_registry/paper3_ici_public_dataset_qc_summary.tsv` (14 dataset)
- `project/results/paper3_ici_data_registry/extract_metadata.py`

### 8.4 Scoring plan (Phase 4)
- `project/reports/paper3_ici/paper3_ici_module_gene_list.tsv` (7 module × 60 gene)
- `project/results/paper3_ici_data_registry/paper3_ici_dataset_gene_coverage.tsv`
- `project/results/paper3_ici_data_registry/check_gene_coverage.py`
- `project/reports/paper3_ici/paper3_ici_scoring_plan.md`

### 8.5 Professor packet (Phase 5)
- `project/reports/paper3_ici/paper3_ici_data_expansion_professor_packet_kr.md` (이 문서)
- `project/reports/paper3_ici/paper3_ici_one_page_data_summary_kr.md`

---

## 9. 가장 강하게 허용되는 결론 (binding)

> **공개 thyroid TME 데이터셋과 thyroid ICI 임상 evidence 통합으로 Paper 3 를 강화할 수 있다. 그러나 현재 방어 가능한 endpoint 는 여전히 ICI-readiness / immunogenomic vulnerability prioritization 이며, thyroid ICI response prediction 은 아니다.**

---

## 10. 다음 분석 명령 (Track B 진입 시 user 가 발화해야 할 정확한 문구)

> **"Track B 시작 — Paper 3 ICI Module A bulk ecotype on TCGA-THCA + GSE126698 + GSE76039 + GSE60542 + GSE151179."**

이 명령은 G1 (Paper 1 bioRxiv 제출) + G2 (Paper 2 Task A/B/C closure) + G3–G5 (controlled access / scRNA / pan-cancer ICI ≥3 코호트) 가 모두 닫힌 후에만 유효.

---

Track A 동결 (Paper 1/2 분석) 유지. 마라톤 비위반.
