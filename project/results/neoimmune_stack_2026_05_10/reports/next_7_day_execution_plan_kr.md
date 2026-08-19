# 7-Day Execution Plan KR

## Day 1
- 병원 데이터 필드 확정: `hospital_data_request_sheet.tsv` 그대로 전달.
- 내부 알고리즘 owner 확정: Structure/Wave8/ESM2/Quantum/BAR-Neo.

## Day 2
- patient-level candidate table 샘플 5명으로 스모크 테스트.
- HLA typing, expression, VAF, HLA LOH, B2M/APM 필수 컬럼 QC.

## Day 3
- external predictor frozen run: NetMHCpan, MHCflurry, BigMHC, PRIME.
- clean track에는 외부 predictor score가 절대 들어가지 않는지 audit.

## Day 4
- local algorithms 재실행: Structure_LR, Wave8_TCR_SelfSim_full, ESM2_Bayesian, quantum_kernel_no_anchor_gamma1.0.
- source/patient/no-overlap split 점검.

## Day 5
- patient top-20 후보 생성.
- high-score/low-expression, HLA-LOH, high-disagreement 후보 제거 또는 control로 이동.

## Day 6
- wet-lab 후보 회의: top candidates + rescue candidates + negative/disagreement controls.
- MS/T-cell assay 가능성 판단.

## Day 7
- paper/사업 pitch 업데이트.
- claim ladder 기준으로 초록/노랑/빨강 claim 분리.
