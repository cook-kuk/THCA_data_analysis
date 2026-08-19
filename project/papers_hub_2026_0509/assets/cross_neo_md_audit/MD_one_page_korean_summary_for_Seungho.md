# Seungho용 1페이지 요약: CROSS-Neo MD 결과

현재 결론은 명확합니다. MD는 면역원성 증명이 아니라 구조 audit 레이어입니다. 그래서 논문 메인 claim으로 쓰기보다, wetlab 후보 우선순위와 false positive/false negative 해석에 쓰는 것이 안전합니다.

`GADGVGKSAL / HLA-C*08:02`는 10 ns explicit-solvent CUDA run이 완료됐고, MD evidence는 `MD_MODERATE`입니다. peptide RMSD final은 0.139 nm이고, peptide-MHC anchor와 TCR-peptide contact가 유지되는 쪽으로 보입니다. 이건 wetlab 우선순위 근거로 좋습니다.

`HMTEVVRHC / HLA-A*02:01`는 아직 running입니다. 최신 동기화 기준 10000.0 ps, 온도 299.89 K라서 run 자체는 정상입니다. 하지만 DCD 전체 분석 전에는 contact/RMSD claim을 하지 않는 게 맞습니다.

다음 액션은 3개입니다. 1) HMTEVVRHC 10 ns 끝나면 DCD sync 후 전체 분석 재실행. 2) GADGVGKSAL과 HMTEVVRHC에 WT/decoy/same-HLA positive control을 붙여 3x10 ns batch 실행. 3) 그 다음에만 50-100 ns 장기 MD로 승격합니다.

절대 하면 안 되는 말: MD가 면역원성을 증명했다, 임상 효능을 보였다, SOTA를 검증했다, 짧은 trajectory 하나로 binding을 증명했다. 지금 가능한 말은 `구조적으로 그럴듯해서 wetlab 후보로 올릴 가치가 있다`입니다.
