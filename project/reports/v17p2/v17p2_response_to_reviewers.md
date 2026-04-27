# v17p2 Response To Reviewers

## Reviewer #1 expected concerns
1. Cluster artifact 주장: Layer 1 `bootstrap_consensus.html`, `robustness_method_consistency.html`로 대응.
2. Age confounding 주장: `age_adjusted_cluster_effect.tsv`에서 residual TDS p=2.252117540968525e-18`.
3. Clinical irrelevance 주장: OS event는 적지만 age/stage/histology stratification을 제시.
4. Mechanism 부족 주장: `gsea_top_pathways_bar.html`, `mapk_activity_violin.html`, `thyroid_diff_score_violin.html`.
5. TCGA-only 주장: `external_4cohort_forest.html`, `scrna_dm_signature_umap.html`.

## Reviewer #2 expected concerns
1. ComBat artifact 주장: `dial_combat_threshold.html`에서 identifiability collapse를 정량화.
2. Alternative metric이면 충분하다는 주장: `dial_vs_alternatives_radar.html`.
3. Sample size 우려: `sample_subsample_stability.tsv`, bootstrap persistence.
4. scRNA projection이 bulk signature recycling이라는 주장: per-patient heterogeneity plot로 mixture를 제시.
5. Pan-cancer generalization 부족: 현재 supplementary limitation으로 명시.

## 예상 공격 10개와 답변
1. Random cluster? → bootstrap persistence 0.903
2. Permutation null? → p=0.0033
3. Age only? → residual TDS after age p=2.252117540968525e-18
4. Survival signal weak? → honest limitation, OS events=8
5. No mechanism? → hallmark proxy FDR<0.05 10
6. No external? → robust bulk cohorts=2
7. No single-cell? → scRNA cells scored=66015
8. Batch effect? → lambda threshold=1.0
9. No pan-cancer method proof? → current limitation section으로 후퇴
10. No wet validation? → limitation에서 정직하게 인정
