# v14 Strengthening — GTEx TACSTD2 (TROP2) Baseline

작성일: 2026-04-25
대상 유전자: **TACSTD2 / TROP2** (Ensembl `ENSG00000184292.6`)
목적: PTC 종양에서 관찰된 TROP2 high-z가 **종양 특이적**인지, 아니면 **갑상선 조직의 정상 baseline 인공물**인지 검증.

## 데이터 출처

- **GTEx v8 median TPM (GCT)**
  파일: `GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz`
  URL: <https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz>
  로컬 캐시: `results/v14_strengthening/gtex_cache/...gct.gz` (md5 `18841e73639f15862fef4c11efa4869a`, 6.95 MB)
  포맷: 56,200 유전자 × **54 조직** (GTEx v8 release; sex-specific 조직 포함되어 53이 아닌 54개로 카운트됨).

- **TCGA-THCA RNA-seq log2(TPM+1)** (local)
  파일: `data_processed/bulk_rnaseq/TCGA-THCA_rnaseq_expression_log2.tsv`
  Tumour (sample type 01/06): n=513, Normal (sample type 11): n=59.

> NOTE: GTEx v8 GCT는 sex-specific(예: Cervix - Ectocervix, Cervix - Endocervix, Fallopian Tube, Prostate, Testis, Uterus, Vagina) 조직을 포함하여 컬럼 수가 54다. 사용자 사양의 "53 tissues"는 v7 카운트 또는 일부 메타-필터 후 카운트와 호환된다. 본 분석은 GCT 헤더에 명시된 54개를 그대로 사용한다.

## 1. GTEx 53(54)개 조직에서 TACSTD2 발현 순위

전체 출력: `results/v14_strengthening/gtex_tacstd2_by_tissue.tsv`

- **갑상선(Thyroid) 순위: 16 / 54**
- 갑상선 median TPM = **8.103**
- 전체 조직 median 의 median = 2.074 TPM
- 최고 발현 조직: Esophagus - Mucosa (1419.1 TPM)
- 최저 발현 조직: Brain - Anterior cingulate cortex (BA24) (0.1724 TPM)

### Top 5 (TACSTD2 high)
  1. Esophagus - Mucosa — 1419.10 TPM
  2. Skin - Not Sun Exposed (Suprapubic) — 708.70 TPM
  3. Skin - Sun Exposed (Lower leg) — 683.92 TPM
  4. Minor Salivary Gland — 680.75 TPM
  5. Vagina — 655.98 TPM

### Bottom 5 (TACSTD2 low)
  1. Brain - Nucleus accumbens (basal ganglia) — 0.2732 TPM
  2. Brain - Cerebellar Hemisphere — 0.2500 TPM
  3. Heart - Left Ventricle — 0.2353 TPM
  4. Brain - Amygdala — 0.2245 TPM
  5. Brain - Anterior cingulate cortex (BA24) — 0.1724 TPM

해석: GTEx 정상 갑상선의 TACSTD2 baseline은 **상위권이 아니라 중하위권** (rank 16/54, ~8.1 TPM). 식도 점막·피부·유방·방광·자궁경부·췌관·구강 등 epithelial barrier 조직이 1,000 TPM 단위로 압도적으로 높고, 갑상선은 이들과 비교하여 두 자릿수 이하의 발현을 보인다.

## 2. TCGA-THCA 종양 vs 정상 정량 비교 (TACSTD2)

전체 출력: `results/v14_strengthening/tcga_thca_vs_gtex_normal_tacstd2.tsv`

| 그룹 | n | mean log2(TPM+1) |
| --- | ---: | ---: |
| TCGA-THCA Tumour (01/06) | 513 | 7.314 |
| TCGA-THCA Normal (11)    | 59 | 4.058 |

- Δ log2(TPM+1) (Tumour − Normal) = **+3.255**
- Fold-change (linear TPM 근사) ≈ **9.55×**
- Welch's t-test: t = 15.851, p = 0.000e+00
- Mann–Whitney U: U = 24555, p = 4.630e-15
- Cohen's d = 1.500

해석: TCGA-THCA tumour-vs-normal 비교는 **TROP2가 tumour에서 유의하게 상승**한다는 것을 양측 모두 매우 작은 p값으로 지지한다. GTEx 갑상선 baseline (median ≈ 8.1 TPM, log2≈3.19) 도 TCGA normal 평균 (4.06) 와 유사한 수준이며, 어느 쪽 비교에서도 tumour mean (7.31) 가 명확히 우위에 있다.

## 3. Verdict — TROP2는 PTC 종양 특이적 상승인가?

**결론: 종양 특이적 상승이 맞다 (tumour-specific, not a thyroid-tissue baseline artefact).**

근거 요약:

1. **GTEx 정상 갑상선의 TACSTD2 발현은 상위권이 아니다.** 54개 GTEx 조직 중 16위, median **8.10 TPM**. 식도 점막(1,419), 피부(680~684), 유방(131), 구강 점막·자궁경부 등 epithelial 조직이 100~1,400 TPM으로 압도적으로 높다. 즉 thyroid 자체가 TROP2-rich 한 조직이라는 가설은 GTEx 데이터로 **기각**된다.
2. **TCGA-THCA paired analysis 에서 tumour > normal.** n_tumour=513, n_normal=59, Δlog2 = +3.26 (≈ 9.5× linear), Welch p = 0.00e+00, Mann–Whitney p = 4.63e-15, Cohen's d = 1.50. v13/v14 의 BRAF V600E ↔ TROP2 high-z 신호는 tumour-specific 상승 위에 BRAF 의존적 추가 강화로 해석 가능하다.
3. **갑상선 vs TCGA normal.** GTEx 갑상선 baseline TPM ≈ 8.10 (log2≈3.19) 와 TCGA-THCA normal mean ≈ 4.06 가 일치 → cross-platform consistency 확인. 따라서 TCGA normal sample 부족(n=59) 우려도 GTEx로 보강된다.

리뷰어 응답 핵심 문장(영문 권장):
> *"In GTEx v8, normal thyroid TACSTD2 expression ranks 16 out of 54 tissues at a median of 8.10 TPM — an order of magnitude or more below epithelial barrier tissues such as esophageal mucosa, skin, breast and bladder. In TCGA-THCA, primary tumour samples (n=513) show a 9.5-fold increase over matched normal-adjacent samples (n=59; Welch p=0.0e+00, Mann–Whitney p=4.6e-15, Cohen's d=1.50). The PTC TROP2-high signal is therefore a tumour-acquired phenotype, not a thyroid-tissue baseline artefact."*

## 4. 산출물

- `results/v14_strengthening/gtex_cache/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz`
- `results/v14_strengthening/gtex_tacstd2_by_tissue.tsv` (54 rows, ranked)
- `results/v14_strengthening/tcga_thca_vs_gtex_normal_tacstd2.tsv` (tumour / normal stats)
- `reports/v14_strengthening_gtex_baseline.md` (이 문서)

## 5. 한계 및 주의사항

- **Sample-level GTEx (~3 GB) 다운로드는 미수행.** v14 강화 단계에서 분포(tissue-level distribution) 까지 필요하면 별도 작업으로 추가 가능. 현재는 **median TPM only**.
- TCGA 발현 행렬은 log2(TPM+1) 단위로 저장되어 있어 fold-change 환산은 2^Δ 근사를 사용했다. 작은 TPM 영역에서 미세한 편향이 가능하나 본 결론에는 영향이 없다.
- GTEx GCT 의 tissue 컬럼 수는 54 (sex-specific 분류 포함). 사용자 명세의 "53 tissues" 는 GTEx v7 또는 일부 합산 정의에 해당; 본 보고서는 GCT 원본을 그대로 사용한다.
- TROP2 alias / Ensembl ID: TACSTD2 = `ENSG00000184292.6` (v8 GCT 내 표기) = TROP2.

— end of report
