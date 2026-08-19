# CROSS-Neo Product Demo One-Pager

## 한 줄 포지션
CROSS-Neo product mode는 allele별 작은 모델이 아니라, 여러 HLA에 공통으로 작동하는 pan-allele neoantigen priority ranker입니다.

## 왜 allele-specific으로 가지 않나
- 현재 n에서는 allele별 모델이 biology를 배우기보다 HLA/source shortcut을 외울 가능성이 큽니다.
- 사업 데모에서는 새 병원/새 환자/새 HLA 조합에 버텨야 하므로 하나의 범용 ranker가 더 맞습니다.
- HLA와 supertype은 버리는 정보가 아니라 context, OOD, calibration 신호로 사용합니다.

## 고객에게 보여줄 출력
- 입력: mutant peptide, WT peptide 가능 시, HLA, source/context metadata 가능 시.
- 출력: priority rank, P0/P1/P2 tier, pan-allele core score, predictor-assisted evidence, model agreement, recommended action.
- Customer-safe queue: `/home/seungho/personal/THCA_data_analysis/project/results/cross_neo_product_demo_2026_05_10/product_demo_priority_queue_customer_safe.tsv`

## 내부 리허설 성능
- Fixed pan-allele product score: AUPRC 0.582, top10 precision 0.700, enrichment@10 2.97x over prevalence.
- Testset-aware upper bound: AUPRC 0.634, top10 precision 0.900. This is internal only and not a generalization claim.
- Public predictor assist alone: AUPRC 0.395, top10 precision 0.400, so the demo is not just a wrapper around public predictors.

## 사용할 train/reference source
- NEPdb: core source calibration.
- TESLA_mmc4: low-prevalence false-positive pressure.
- CEDAR: positive-rich prior, heavily downweighted.
- ITSNdb: domain anchor/reference-sensitive support.
- TESLA_mmc7: abstention stress only.

## Claim boundary
- Clean manuscript claim과 product demo claim을 분리합니다.
- Product mode may be predictor-assisted and testset-aware for retrospective demo strategy.
- Do not call this external validation.
- Do not claim quantum advantage.
