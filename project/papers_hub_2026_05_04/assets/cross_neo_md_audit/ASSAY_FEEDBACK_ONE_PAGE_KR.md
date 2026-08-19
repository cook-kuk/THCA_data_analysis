# CROSS-Neo assay feedback learner 요약

이제 구조가 실제로 data-driven closed loop가 됩니다.

1. CROSS-Neo가 후보 prior를 냅니다.
2. mutant/WT/decoy assay TSV를 넣습니다.
3. posterior readiness, claim state, 다음 실험 우선순위가 자동으로 바뀝니다.
4. WT positive 또는 decoy positive이면 mutant-specific claim이 자동으로 막힙니다.
5. activation positive + WT/decoy negative이면 assay-specific immunogenicity claim으로 올라갈 수 있습니다.

현재 모드: 아직 실제 assay label 없음; 템플릿과 다음 실험 큐 생성

중요: 이 점수는 실험 우선순위용 posterior readiness입니다. 임상 확률이나 면역원성 확정이 아닙니다.
