# CROSS-Neo high-impact one-page summary

## 핵심 결론

- 현재 label 기준 strict preset은 Top-13 후보를 남기고, 이 후보들은 모두 positive label입니다.
- 이 결과는 wetlab 후보 선별용 decision-support evidence이며, 외부검증/임상효용/면역원성 증명은 아닙니다.
- 가장 바로 실험할 후보는 GADGVGKSAL / HLA-C*08:02, HMTEVVRHC / HLA-A*02:01, KLILWRGLK / HLA-A*03:01, ILDKVLVHL / HLA-A*02:01입니다.

## 왜 하이임팩트인가

- cheap DL + Bayesian/dropout uncertainty로 먼저 줄입니다.
- TCR/구조/MD는 비싼 최종 audit layer로만 씁니다.
- 후보마다 required control과 go/no-go 기준이 붙어 wetlab plate로 바로 전환됩니다.

## Source별 현재 성능

- ITSNdb_Val: called=1, TP=1, FP=0, recall=0.200
- ITSNdb_main: called=5, TP=5, FP=0, recall=0.312
- NEPdb: called=7, TP=7, FP=0, recall=0.050

## 실험 우선순위

- 1. GADGVGKSAL / HLA-C*08:02: TIER_1_IMMEDIATE - optimized no-FP preset; paired TCR evidence=13; mainDL=0.770; TCRbranch=0.960; MD=MD_MODERATE
- 2. HMTEVVRHC / HLA-A*02:01: TIER_1_IMMEDIATE - optimized no-FP preset; paired TCR evidence=18; mainDL=0.467; TCRbranch=0.983; MD=MD_VERY_STRONG
- 3. KLILWRGLK / HLA-A*03:01: TIER_1_IMMEDIATE - optimized no-FP preset; paired TCR evidence=0; mainDL=0.803; TCRbranch=0.985; MD=nan
- 4. ILDKVLVHL / HLA-A*02:01: TIER_1_IMMEDIATE - optimized no-FP preset; paired TCR evidence=0; mainDL=0.773; TCRbranch=0.962; MD=nan
- 5. YVDFREYEYY / HLA-A*01:01: TIER_2_BACKUP - optimized no-FP preset; paired TCR evidence=0; mainDL=0.661; TCRbranch=0.839; MD=nan
- 6. GADGVGKSA / HLA-C*08:02: TIER_2_BACKUP - optimized no-FP preset; paired TCR evidence=0; mainDL=0.781; TCRbranch=0.567; MD=nan
- 7. SYLDSGIHF / HLA-A*24:02: TIER_2_BACKUP - optimized no-FP preset; paired TCR evidence=0; mainDL=0.550; TCRbranch=0.836; MD=nan
- 8. ILDTAGREEY / HLA-A*01:01: TIER_2_BACKUP - optimized no-FP preset; paired TCR evidence=0; mainDL=0.604; TCRbranch=0.719; MD=nan

## 금지 claim

- MD가 immunogenicity를 증명한다고 말하지 않습니다.
- 현재 label 최적화 결과를 외부검증이라고 말하지 않습니다.
- TCR/structure가 없는 후보에서 TCR recognition을 증명했다고 말하지 않습니다.
