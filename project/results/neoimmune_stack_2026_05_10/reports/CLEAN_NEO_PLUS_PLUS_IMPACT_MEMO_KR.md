# CLEAN-Neo++ / NeoImmune-Stack Impact Memo KR

## 한 줄 결론
새 foundation model이 아니라, **기존 내 알고리즘 + 공개 predictor + leakage audit + patient top-N handoff**를 묶은 운영 레이어로 가는 것이 가장 임팩트가 큽니다.

## 대박 포인트
- 내 알고리즘이 중심입니다: Structure_LR, Wave8/TCR-self-similarity, ESM2_Bayesian, W7A/W7B, quantum kernel 계열을 clean science track에 배치했습니다.
- 공개 predictor는 보조입니다: NetMHCpan/MHCflurry/BigMHC/PRIME은 production track의 frozen feature/comparator입니다.
- reviewer 방어력이 생겼습니다: clean vs production을 분리했고, strict no-existing-overlap 리더보드를 별도로 만들었습니다.
- 사업화 언어가 생겼습니다: NeoImmune-Stack은 병원/환자별 후보 prioritization operating layer입니다.

## 냉정한 약점
- 현재 public 통합 테이블은 real patient ID가 86명 / 16,298 rows라서 retrospective public patient-level ranking scaffold는 가능해졌습니다.
- 하지만 hospital-grade prospective endpoint, clinical utility, assay-confirmed presentation/immunogenicity는 아직 주장하면 안 됩니다.
- MS immunopeptidomics 없이 presentation claim을 하면 바로 공격받습니다.
- T-cell assay label 없이 immunogenicity claim을 하면 바로 공격받습니다.
- BAR-Neo 쪽 production 숫자가 강해 보여도 clean novel model claim으로 쓰면 안 됩니다.

## 지금 가장 강한 논문 문장
“CLEAN-Neo++ provides a leakage-aware strategy for integrating peptide-HLA presentation, TCR-visible immunogenicity, structure proxies, and quantum-kernel features into patient-level neoantigen candidate ranking.”

## 바로 필요한 병원 데이터
patient_id, sample_id, tumor type, HLA typing, mutation peptide, WT peptide, expression TPM, VAF, clonality, HLA LOH, B2M/APM, candidate source, assay labels, and top-N validation outcome.
