# 공격적 새 주제 feasibility sprint 요약

Date: 2026-05-06

## 결론

기존 dark-matter taxonomy 재포장은 약합니다. 오늘 빠르게 찔러본 결과, 가장 새롭고 공격적인 주제는 다음입니다.

> **BRAF/RAS/DICER1/EIF1AX/TERT-negative thyroid cancer 안에 존재하는 copy-number-defined residual class**

즉, "driver-negative"라고 부르던 Class6 안에서도 DM2/FVPTC-like 쪽 일부가 **7p/7q gain, 2p/2q loss, 16p/16q gain, 12q gain** 같은 arm-level CNV signature를 갖습니다. 이건 Paper 1의 differentiation score 반복이 아니라, **genomic architecture** 쪽이라 독립성이 더 좋습니다.

## 오늘 테스트한 후보별 판정

| 순위 | 주제 | 오늘 결과 | 판정 |
|---:|---|---|---|
| 1 | Class6 true-driver-negative CNV class | Class6에서 7q gain q~9.6e-5, 7p gain q~5.0e-5, 2p/2q loss q~5-7e-5. Signature는 Class6-DM2 14/45 vs Class6-DM1 1/65 | **최우선 공격 주제** |
| 2 | DM1 hidden fusion/regulatory driver | 기존 f4 fusion 결과는 DM1 fusion+ 76.8%로 강하지만, v17 curated overlay는 ALK/NTRK/RET 각 1개라 충돌 | **고위험 고수익, 재검증 필요** |
| 3 | Epigenetic silencing mechanism | 8-gene methylation과 RAI/TDS score가 강한 음의 상관: rho -0.55~-0.58. TPO methylation DM1 vs DM2 delta beta 0.415, p=1.88e-18 | 강하지만 Paper 1과 겹침 |
| 4 | DICER1/EIF1AX class | TCGA n=7, K2에서도 DM 안 DICER1+EIF1AX n=7 | section으로는 좋지만 단독 논문 약함 |
| 5 | RAI-refractory classifier | GSE151179 직접 재계산. 최고 response AUC는 TPO oriented AUC 0.643, p NS. Primary-only도 약함 | **headline no-go** |
| 6 | TERT-only triple-negative | TCGA n=5, PFI 2/5 | 흥미롭지만 외부 cohort 전에는 hold |

## 제일 새 주제: CNV-defined residual class

기존 driver class:

- Class1 BRAF V600E: n=291
- Class2 RAS hotspot: n=54
- Class4 DICER1/EIF1AX: n=7
- Class5 TERT-only: n=5
- **Class6 true-driver-negative: n=125**

Class6 안에서 arm-CNV가 강하게 enrich된 arm:

| arm | event | Class6 % | others % | Fisher p | BH q |
|---|---|---:|---:|---:|---:|
| 7q | gain | 10.0 | 0.30 | 1.23e-6 | 9.56e-5 |
| 7p | gain | 9.91 | 0.29 | 1.29e-6 | 5.04e-5 |
| 2q | loss | 8.11 | 0.00 | 2.57e-6 | 6.69e-5 |
| 2p | loss | 8.11 | 0.00 | 2.57e-6 | 5.02e-5 |
| 16p | gain | 9.01 | 0.29 | 5.17e-6 | 8.07e-5 |
| 16q | gain | 8.11 | 0.29 | 2.04e-5 | 2.27e-4 |

이 CNV signature는 DM2 쪽에 집중됩니다.

- Class6-DM2: **14/45**
- Class6-DM1: **1/65**

이건 새로운 논문 제목으로 밀 수 있습니다.

> **Copy-number architecture reveals a residual oncogenic class within BRAF/RAS-negative thyroid cancer**

또는 더 공격적으로:

> **True driver-negative thyroid cancers are not driverless: recurrent arm-level copy-number states define a follicular-pattern residual class**

단, "driver"라고 하려면 추가 검증이 필요합니다. 지금은 **copy-number-defined residual class**가 안전합니다.

## RAI-refractory 직접 검증 결과

GSE151179를 GEO에서 직접 파싱했습니다.

- 총 52 samples
- tumor 39
- primary tumor 17
- 8/8 RAI-lineage genes probe 회수 성공

RAI response refractory vs avid:

- best response signal: TPO
- tumor_all n=39, refractory 35, avid 4
- refractory mean 8.34 vs avid 9.08
- oriented AUC 0.643
- p=0.38

결론: 방향은 일부 맞지만 imbalance가 심하고 유의하지 않습니다. **RAI-refractory standalone 논문은 오늘 기준 no-go**입니다.

## 다음 5시간 안에 이어서 할 일

1. Class6 CNV signature를 cBioPortal/GDC GISTIC source에서 독립 확인합니다.
2. CNV+ Class6-DM2가 FVPTC/hurthle/follicular-pattern artifact인지, 진짜 residual genomic class인지 확인합니다.
3. CNV+ vs CNV- Class6에서 expression, methylation, histology, age, stage, PFI를 비교합니다.
4. fusion evidence는 source conflict가 있으므로 cBioPortal structural variant 원표와 v17 overlay를 대조합니다.

## 오늘 산출물

- `aggressive_topic_sprint_report.md`
- `aggressive_sprint_summary.json`
- `gse151179_rai_validation_summary.tsv`
- `gse151179_rai_scores.tsv`
- `driver_class_score_summary.tsv`
- `class6_arm_cnv_enrichment.tsv`
- `class6_dm1_dm2_arm_cnv.tsv`
- `class6_cnv_signature_samples.tsv`
- `methylation_expression_correlations.tsv`

