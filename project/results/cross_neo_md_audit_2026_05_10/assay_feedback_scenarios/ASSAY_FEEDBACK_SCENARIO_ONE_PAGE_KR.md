# assay feedback scenario simulator 요약

이 패키지는 실제 실험값이 들어왔을 때 posterior와 claim이 어떻게 바뀌는지 미리 보여주는 what-if 시뮬레이터입니다.

- mutant activation positive + WT/decoy negative: assay-specific immunogenicity claim까지 상승 가능
- WT positive: mutant-specific claim 차단
- decoy positive: nonspecific/artifact 위험으로 claim 차단
- presentation positive만 있음: presentation plausibility까지만 가능
- recognition positive but activation negative: TCR binding은 있어도 immunogenicity claim은 불가

중요: 여기 rows는 전부 SIMULATED_NOT_REAL_ASSAY_LABEL입니다. 실제 label TSV와 분리되어 있습니다.
