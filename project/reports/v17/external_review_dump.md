# v17 External Review Dump


> ⚠️ **v5.2 retraction notice (2026-04-25):** some DIA-AUC and DIAL numbers below were computed under v5.1 leaky-ComBat protocol (full-pooled-data ComBat before LODO split, information leak). Under proper per-fold ComBat, the numbers shift. Current submission-ready figures are at `reports/html/pages/v17_npj_robustness.html` (manuscript v3 scenario-B reframe). See `reports/v5p2/v5p2_critical_assessment.md`.


## 1. 한 줄 요약

- v17은 TCGA-THCA의 `BRAF/RAS-negative` 178명이 단일 residual bucket이 아니라 **DM1/DM2 두 개의 안정적 expression subtype**으로 갈라짐을 보였고, 이 구조가 external cohort transfer와 trajectory에서도 일정 부분 유지됨을 확인했다.
- v14가 `BRAF_like vs RAS_like` 축의 재현성 검증이었다면, v17은 **그 밖에 남는 표본 자체를 새로운 subtype problem으로 재정의**한 점이 다르다.

## 2. 데이터 인벤토리 (numerical)

### Cohort별 sample n

| Cohort | total n | tumor n | normal n |
|---|---:|---:|---:|
| TCGA-THCA | 572 | 513 | 59 |
| GSE27155 | 99 | 95 | 4 |
| GSE126698 | 28 | 22 | 6 |
| GSE76039 | 37 | 37 | 0 |
| GSE213647 | 632 | 381 | 251 |
| GSE97466 | 141 | 74 | 67 |

### `driver_anchor_v17` 분포

| group | overall n |
|---|---:|
| unknown | 1072 |
| BRAF | 343 |
| RAS | 74 |
| PAX8PPARG | 7 |
| DICER1_EIF1AX_PPM1D | 7 |
| RET | 2 |
| RET_fusion | 1 |
| NTRK_fusion | 1 |
| ALK_fusion | 1 |
| TP53 | 1 |

### cohort별 `driver_anchor_v17` (tumor only)

| Cohort | top distribution |
|---|---|
| TCGA-THCA | `BRAF 284 (55.4%)`, `unknown 164 (32.0%)`, `RAS 54 (10.5%)`, `DICER1/EIF1AX/PPM1D 7 (1.4%)`, rare `NTRK/ALK/RET_fusion/TP53` each `1` |
| GSE27155 | `unknown 45 (47.4%)`, `BRAF 28 (29.5%)`, `RAS 13 (13.7%)`, `PAX8PPARG 7 (7.4%)`, `RET 2 (2.1%)` |
| GSE126698 | `unknown 22 (100%)` |
| GSE76039 | `unknown 37 (100%)` |
| GSE213647 | `unknown 381 (100%)` |
| GSE97466 | `unknown 74 (100%)` |

### Dark Matter / 4-group / fusion / histology / stage

- **Dark Matter n**: `178` (`TCGA-THCA tumor`, `driver_anchor ∉ {BRAF,RAS}`)
- **4-group (degraded)**: `A_braf_only=281`, `B_ras_only=54`, `D_triple_negative=178`, `C_tert_plus=0`
- **TERT caveat**: local WXS MAF에서 promoter hotspot 비관측 → `+TERT arm` 자체가 사라짐
- **Fusion/driver proxy final n**: `RET_fusion=1`, `NTRK_fusion=1`, `ALK_fusion=1`, `PAX8PPARG=7`
- **Proxy source**: `project/metadata/v3_fusion_anchor_tcga.tsv` (TCGA only; `BRAF_V600E / RAS_mutant / other` 성격의 legacy anchor proxy)
- **True Dark Matter (all drivers unresolved)**: `driver_anchor_v17 == unknown = 1072` overall, `164` in `TCGA-THCA tumor`

| Histology | n |
|---|---:|
| cPTC | 1066 |
| FVPTC | 109 |
| PDTC | 25 |
| ATC | 57 |
| FTC | 30 |
| MTC | 2 |
| unknown | 162 |
| normal | 58 |

| TCGA tumor AJCC stage | n |
|---|---:|
| Stage I | 268 |
| Stage II | 48 |
| Stage III | 99 |
| Stage IVA | 47 |
| Stage IVC | 3 |
| Stage IV | 1 |
| Stage 0a | 1 |
| NA | 46 |

## 3. Idea 1 — Dark Matter 결과 (Numerical)

- **Best K**: `2` (consensus stability `0.9788`)
- **Cluster size**: `DM1=109`, `DM2=69`

### Top 5 marker

| Cluster | gene | log2FC | FDR |
|---|---|---:|---:|
| DM1 | DUSP5 | 2.789 | 1.72e-33 |
| DM1 | SLC5A8 | -3.017 | 7.46e-29 |
| DM1 | DUSP6 | 2.382 | 1.20e-26 |
| DM1 | DIO1 | -3.886 | 1.94e-25 |
| DM1 | TPO | -3.910 | 4.07e-25 |
| DM2 | DUSP5 | -2.789 | 1.72e-33 |
| DM2 | SLC5A8 | 3.017 | 7.46e-29 |
| DM2 | DUSP6 | -2.382 | 1.20e-26 |
| DM2 | DIO1 | 3.886 | 1.94e-25 |
| DM2 | TPO | 3.910 | 4.07e-25 |

### Cluster-level biology

| Cluster | n | TDS16 mean | BRS71 mean | RAI mean | age mean |
|---|---:|---:|---:|---:|---:|
| DM1 | 109 | 6.812 | 2.197 | 7.580 | 41.818 |
| DM2 | 69 | 8.126 | 0.753 | 9.201 | 54.708 |

| Cluster | SLC5A5 | TPO | TSHR | TG | DIO1 | DIO2 |
|---|---:|---:|---:|---:|---:|---:|
| DM1 | 1.571 | 7.731 | 8.807 | 14.218 | 4.052 | 6.902 |
| DM2 | 2.319 | 11.641 | 9.256 | 15.502 | 7.938 | 9.025 |

### `v14 molecular_subtype` vs `v17_dark_cluster`

| v17 cluster | BRAF_like | RAS_like | unknown |
|---|---:|---:|---:|
| DM1 | 86 | 18 | 5 |
| DM2 | 25 | 39 | 5 |

## 4. Idea 2 — TERT 4-group (Caveat Heavy)

- **Reality**: local GDC WXS MAF 기반 실행에서 `TERT promoter mutated = 0 / 513`
- **Degradation**: 원래 의도한 `BRAF only / RAS only / +TERT / triple-negative` 4-group이 실제로는 `3-group` (`+TERT arm` 소실)으로 붕괴
- **Fix path**:
  - `PanCanAtlas / cBioPortal THCA mutation supplement` 재수집
  - `Liu et al. 2017 JCO` 계열 TERT promoter supplement 매칭
  - 필요 시 `MUT: promoter hotspot table` 직접 병합
- **Paper impact (quantified)**:
  - 현재 `+TERT` 샘플 수 `0`
  - 공격적 clinical arm 하나를 통째로 잃음
  - survival/KM에서 “aggressive +TERT group” 비교 자체가 불가능
  - 따라서 v17의 임상 narrative는 현재 **driver-negative subtype paper**로는 유지 가능하지만, **driver escalation paper**로는 불완전

## 5. Idea 3 — Fusion Landscape

- **Current frequency (final `driver_anchor_v17`)**: `RET_fusion=1`, `NTRK_fusion=1`, `ALK_fusion=1`, `PAX8PPARG=7`
- **cohort별 fusion-like calls**:
  - `TCGA-THCA`: `RET_fusion 1`, `NTRK_fusion 1`, `ALK_fusion 1`, `PAX8PPARG 0`
  - `GSE27155`: `PAX8PPARG 7`, `RET 2`
  - others: `0`
- **`v3_fusion_anchor_tcga.tsv` 정체**: TCGA sample별 `v3_anchor_6class` 를 담은 internal proxy file (`BRAF_V600E`, `RAS_mutant`, `other`), raw fusion partner annotation 없음
- **TumorFusions / Stransky 2018 raw call과의 차이**:
  - raw call은 fusion partner (`CCDC6-RET`, `NCOA4-RET`, `ETV6-NTRK3` 등)와 breakpoint level 정보 제공
  - 현재 proxy는 partner 정보 없음
  - 따라서 one-to-one 비교 불가, 희귀 fusion under-call 가능성 큼
- **Priority cascade 후 최종 top 분포**: `unknown 1072 > BRAF 343 > RAS 74 > PAX8PPARG 7 = DICER1_EIF1AX_PPM1D 7 > RET 2`
- **“true Dark Matter” final n**: `1072 overall`, `164 TCGA tumor`

## 6. Idea 4 — DIAL Audit (CRITICAL)

| Cohort | Cluster | AUC | DIA-AUC | Identifiability no ComBat | Identifiability ComBat |
|---|---|---:|---:|---:|---:|
| GSE27155 | DM1 | 0.969 | 0.969 | 1.000 | 0.258 |
| GSE27155 | DM2 | 0.969 | 0.969 | 1.000 | 0.258 |
| GSE76039 | DM1 | 0.990 | 0.990 | 1.000 | 0.157 |
| GSE76039 | DM2 | 0.990 | 0.990 | 1.000 | 0.157 |

- **Median DIA-AUC**: `0.979` (8 rows 기준)
- **v4 finding 재현 여부**:
  - exact `1.000 → 0.487`, `AUC 0.989` 는 **수치상 동일하게 재현되지 않음**
  - 그러나 **정성적 방향은 강하게 재현**: `1.000/0.9998 → 0.2577/0.1573`, `DIA-AUC 0.969-1.000`
- **Robust cohort**: `GSE27155`, `GSE76039`
- **깨지는 cohort / 미실행 cohort**:
  - `GSE126698`, `GSE213647`: usable alignment empty
  - `GSE97466`: modality mismatch로 audit 미실행
- **Nature Cancer implication**: 현재는 “2-cohort robust + 3-cohort unresolved” 상태라서 high-impact clinical generalization claim에는 부족

## 7. Idea 5 — Trajectory

- **Usable trajectory subset**: `645` samples
- **Histology projection (usable subset)**:

| Histology | n | pseudotime mean | pseudotime median |
|---|---:|---:|---:|
| cPTC | 431 | 0.739 | 0.786 |
| FVPTC | 102 | 0.490 | 0.499 |
| PDTC | 17 | 0.561 | 0.546 |
| ATC | 20 | 0.823 | 0.852 |

- **Interpretation**:
  - `ATC`가 가장 뒤쪽 tail
  - `PDTC`는 `FVPTC`보다 뒤, `ATC`보다는 앞
  - Dark Matter cluster 중에서는 `DM1 mean pseudotime = 0.736`, `DM2 = 0.295` → **DM1이 고위험 tail에 더 가까움**

### Top 20 dynamic genes

```text
TPO, SLC26A4, DIO1, TG, DUSP5, DIO2, PAX8, MET, LOX, DUOX2,
DUOX1, HLA-DRA, KLK10, FOXE1, DUSP6, CDKN2A, DUSP4, TSHR,
CDKN2B, THRA
```

## 8. Honest Caveats (전체)

- `TERT promoter` 비관측 (`0 / 513`)
- Fusion raw call 부재, `v3_fusion_anchor_tcga.tsv` proxy 재사용
- Wet validation `0`
- `TCGA-CDR` remote fetch는 실패했고, local `tcga_thca_clinical_extended.tsv` fallback으로 `OS`만 사용
- Korean cohort 부재
- Sample size per external cluster validation이 `GSE27155 n=54`, `GSE76039 n=37` 수준으로 작음

## 9. Self-assessment vs Venue

| Item | Score | Note |
|---|---|---|
| Sample size per cluster | 89/100 | DM1 109, DM2 69, 평균 89 |
| Cross-cohort validation | 2/5 | robust는 GSE27155, GSE76039만 |
| Method novelty | medium | DIAL audit를 subtype transfer에 적용했지만 external proof는 부분적 |
| Biological discovery | medium | Dark Matter 2-cluster 구조와 marker set은 새롭지만 driver 확장은 약함 |
| Wet validation | 0 | 부재 |
| Clinical actionability | low | 직접 약물 연결은 약하고 subtype/trajectory 중심 |
| Korean cohort | 0 | 부재 |

## 10. 예상 reviewer 공격 5개

1. `TERT` promoter가 안 보이는데 왜 paper title/narrative에 TERT integration을 넣었는가?
2. Fusion landscape가 raw fusion caller가 아니라 proxy reuse인데, rare driver claim을 얼마나 신뢰할 수 있는가?
3. Cross-cohort validation이 실제로는 5-cohort가 아니라 2-cohort robust 아닌가?
4. DM1/DM2는 unsupervised artifact가 아니라는 강한 orthogonal validation이 있는가?
5. Wet validation 없이 이 subtype이 임상적으로 actionable하다고 말할 수 있는가?

## 11. v15 paper와의 충돌 점검

- **직접 충돌**: 없음. v15는 `BRAF × TROP2` 중심의 therapy/theory narrative, v17은 subtype taxonomy narrative
- **간접 충돌**: 가능. v15가 BRAF/TROP2 축을 더 강하게 전면화하는 반면, v17은 driver-negative space를 강조
- **TROP2 위치**: 현재 v17 산출물에서 TROP2를 별도 중심 marker로 재정의하지 않았고, v15 narrative 쪽이 더 직접적
- **동시 submit 시 self-citation**:
  - v15를 “BRAF-associated therapeutic branch paper”
  - v17을 “taxonomy / subtype expansion paper”
  - 서로 데이터 재사용을 명시하고, claim scope를 분리해야 안전

## 12. 다음 결정 5개

1. `v15 NeurIPS` 일정 유지 vs `v17` 우선 제출
2. `TERT fetch` 1주 추가 재시도 여부
3. Fusion raw reanalysis 여부 (`FASTQ`/raw caller면 4-8주)
4. 분당서울대 / Korean cohort 답변 대기 vs paper 1 먼저 submit
5. Venue 결정: `Bioinformatics short` / `npj Precision Oncology` / `Genome Medicine`
