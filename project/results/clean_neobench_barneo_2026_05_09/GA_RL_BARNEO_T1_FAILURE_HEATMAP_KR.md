# GA/RL BAR-Neo T1 Failure Heatmap KR

## 한 줄 결론

이 패키지는 `GA/RL`에서 살아남은 3개 T1 후보를 **fail / pass로 단순 재분류하지 않고**, exact / near / source / HLA / public overlap 축과 nearest-neighbor 분포를 함께 보여준다. 현재 3개는 clean claim에는 충분히 안전하지만, source/HLA contract와 학습 분포 근접성은 여전히 추적해야 한다.

## Claim boundary

Allowed:

- leakage-aware benchmark analysis
- benchmark-adaptive reliability ranking
- reviewer-safe distribution audit

Forbidden:

- clinical vaccine selection
- new SOTA predictor
- public clean comparator claim without overlap audit
- quantum advantage

## Candidate summary

| candidate_id   | peptide   | hla_allele_4digit   | source_name   |   label |   ga_rl_barneo_priority_score |   barneo_ng_score | leakage_risk_level   | split_source_heldout   | split_hla_heldout   | split_low_prevalence   |   source_positive_prevalence |   source_support_count |   hla_positive_prevalence |   hla_support_count | top_nn_candidate_id   | top_nn_peptide   |   top_nn_similarity |   top_nn_label | top_nn_source_name   | top_nn_hla_allele_4digit   |
|:---------------|:----------|:--------------------|:--------------|--------:|------------------------------:|------------------:|:---------------------|:-----------------------|:--------------------|:-----------------------|-----------------------------:|-----------------------:|--------------------------:|--------------------:|:----------------------|:-----------------|--------------------:|---------------:|:---------------------|:---------------------------|
| CNV0_02407     | KLMNIQQKL | HLA-A*02:01         | ITSNdb_main   |       1 |                      0.729214 |              0.72 | low                  | ITSNdb_main            | HLA-A*02:01         | not_low_prevalence     |                     0.648241 |                    199 |                  0.348348 |                 666 | CNV0_01694            | KLKNKQQQL        |            0.666667 |              0 | TESLA_mmc4           | HLA-B*08:01                |
| CNV0_02410     | MLGEQLFPL | HLA-A*02:01         | ITSNdb_main   |       1 |                      0.72     |              0.72 | low                  | ITSNdb_main            | HLA-A*02:01         | not_low_prevalence     |                     0.648241 |                    199 |                  0.348348 |                 666 | CNV0_00197            | ALGEQVIAL        |            0.555556 |              1 | CEDAR                | HLA-A*02:04                |
| CNV0_02504     | LLVDLAEEL | HLA-A*02:01         | ITSNdb_main   |       1 |                      0.72     |              0.72 | low                  | ITSNdb_main            | HLA-A*02:01         | not_low_prevalence     |                     0.648241 |                    199 |                  0.348348 |                 666 | CNV0_00394            | MLAQLLAEL        |            0.444444 |              1 | CEDAR                | HLA-A*02:04                |

## Contract heatmap

| candidate_id   | axis               | status   |   value | detail                       |
|:---------------|:-------------------|:---------|--------:|:-----------------------------|
| CNV0_02407     | exact_phla_risk    | clear    |       0 | no exact peptide-HLA overlap |
| CNV0_02407     | hla_contract       | heldout  |       1 | held-out HLA=HLA-A*02:01     |
| CNV0_02407     | korean_hla_focus   | focus    |       1 | Korean HLA focus slice       |
| CNV0_02407     | low_prevalence     | ok       |       0 | not low prevalence           |
| CNV0_02407     | near_peptide_risk  | clear    |       0 | no near-peptide cluster hit  |
| CNV0_02407     | public_overlap     | clear    |       0 | no_inhouse_master_overlap    |
| CNV0_02407     | source_contract    | heldout  |       1 | held-out source=ITSNdb_main  |
| CNV0_02407     | supertype_contract | heldout  |       1 | held-out supertype=A02       |
| CNV0_02410     | exact_phla_risk    | clear    |       0 | no exact peptide-HLA overlap |
| CNV0_02410     | hla_contract       | heldout  |       1 | held-out HLA=HLA-A*02:01     |
| CNV0_02410     | korean_hla_focus   | focus    |       1 | Korean HLA focus slice       |
| CNV0_02410     | low_prevalence     | ok       |       0 | not low prevalence           |
| CNV0_02410     | near_peptide_risk  | clear    |       0 | no near-peptide cluster hit  |
| CNV0_02410     | public_overlap     | clear    |       0 | no_inhouse_master_overlap    |
| CNV0_02410     | source_contract    | heldout  |       1 | held-out source=ITSNdb_main  |
| CNV0_02410     | supertype_contract | heldout  |       1 | held-out supertype=A02       |
| CNV0_02504     | exact_phla_risk    | clear    |       0 | no exact peptide-HLA overlap |
| CNV0_02504     | hla_contract       | heldout  |       1 | held-out HLA=HLA-A*02:01     |
| CNV0_02504     | korean_hla_focus   | focus    |       1 | Korean HLA focus slice       |
| CNV0_02504     | low_prevalence     | ok       |       0 | not low prevalence           |

## Nearest training neighbors

| query_candidate_id   | query_peptide   | query_source_name   | query_hla_allele_4digit   | candidate_id   | peptide   | source_name   | hla_allele_4digit   |   label |   similarity | same_source   | same_hla   | same_supertype   | exact_peptide_hla_train_overlap   | near_peptide_train_overlap   | leakage_risk_level   |
|:---------------------|:----------------|:--------------------|:--------------------------|:---------------|:----------|:--------------|:--------------------|--------:|-------------:|:--------------|:-----------|:-----------------|:----------------------------------|:-----------------------------|:---------------------|
| CNV0_02407           | KLMNIQQKL       | ITSNdb_main         | HLA-A*02:01               | CNV0_01694     | KLKNKQQQL | TESLA_mmc4    | HLA-B*08:01         |       0 |     0.666667 | False         | False      | False            | True                              | True                         | high                 |
| CNV0_02407           | KLMNIQQKL       | ITSNdb_main         | HLA-A*02:01               | CNV0_00891     | NLAEAQQKL | CEDAR         | HLA-A*02:02         |       1 |     0.555556 | False         | False      | True             | True                              | True                         | high                 |
| CNV0_02407           | KLMNIQQKL       | ITSNdb_main         | HLA-A*02:01               | CNV0_01901     | KLINSQINL | TESLA_mmc4    | HLA-A*02:01         |       0 |     0.555556 | False         | True       | True             | True                              | True                         | high                 |
| CNV0_02407           | KLMNIQQKL       | ITSNdb_main         | HLA-A*02:01               | CNV0_00051     | GLMATPAKL | CEDAR         | HLA-A*02:01         |       1 |     0.444444 | False         | True       | True             | True                              | True                         | high                 |
| CNV0_02407           | KLMNIQQKL       | ITSNdb_main         | HLA-A*02:01               | CNV0_00109     | RYMNHMQSL | CEDAR         | HLA-A*23:01         |       1 |     0.444444 | False         | False      | False            | True                              | True                         | high                 |
| CNV0_02410           | MLGEQLFPL       | ITSNdb_main         | HLA-A*02:01               | CNV0_00197     | ALGEQVIAL | CEDAR         | HLA-A*02:04         |       1 |     0.555556 | False         | False      | True             | True                              | True                         | high                 |
| CNV0_02410           | MLGEQLFPL       | ITSNdb_main         | HLA-A*02:01               | CNV0_00336     | QLAEMIFNL | CEDAR         | HLA-A*02:01         |       1 |     0.444444 | False         | True       | True             | True                              | True                         | high                 |
| CNV0_02410           | MLGEQLFPL       | ITSNdb_main         | HLA-A*02:01               | CNV0_00338     | MLTERMFNA | CEDAR         | HLA-A*02:01         |       1 |     0.444444 | False         | True       | True             | True                              | True                         | high                 |
| CNV0_02410           | MLGEQLFPL       | ITSNdb_main         | HLA-A*02:01               | CNV0_00350     | QIVERLFSL | CEDAR         | HLA-B*08:01         |       1 |     0.444444 | False         | False      | False            | True                              | True                         | high                 |
| CNV0_02410           | MLGEQLFPL       | ITSNdb_main         | HLA-A*02:01               | CNV0_00394     | MLAQLLAEL | CEDAR         | HLA-A*02:04         |       1 |     0.444444 | False         | False      | True             | True                              | True                         | high                 |
| CNV0_02504           | LLVDLAEEL       | ITSNdb_main         | HLA-A*02:01               | CNV0_00394     | MLAQLLAEL | CEDAR         | HLA-A*02:04         |       1 |     0.444444 | False         | False      | True             | True                              | True                         | high                 |
| CNV0_02504           | LLVDLAEEL       | ITSNdb_main         | HLA-A*02:01               | CNV0_00595     | EFVSRAEER | CEDAR         | HLA-A*33:03         |       1 |     0.444444 | False         | False      | False            | True                              | True                         | high                 |
| CNV0_02504           | LLVDLAEEL       | ITSNdb_main         | HLA-A*02:01               | CNV0_00669     | LLVSGSNVL | CEDAR         | HLA-A*02:01         |       1 |     0.444444 | False         | True       | True             | True                              | True                         | high                 |
| CNV0_02504           | LLVDLAEEL       | ITSNdb_main         | HLA-A*02:01               | CNV0_00767     | EYLDRAEKL | CEDAR         | HLA-A*23:01         |       1 |     0.444444 | False         | False      | False            | True                              | True                         | high                 |
| CNV0_02504           | LLVDLAEEL       | ITSNdb_main         | HLA-A*02:01               | CNV0_01095     | FLLDEAIGL | NEPdb         | HLA-A*02:01         |       1 |     0.444444 | False         | True       | True             | True                              | True                         | high                 |

## Distribution interpretation

- `CNV0_02407`: clean on exact/near/public overlap, but nearest neighbor is a similar `TESLA_mmc4` negative at similarity `0.667`.
- `CNV0_02504`: clean on exact/near/public overlap, but the closest local neighborhood is peptide-similar and mostly `CEDAR` negatives.
- `CNV0_02410`: clean on exact/near/public overlap, and the nearest neighborhood includes a known positive at similarity `0.556`, so this is the strongest distributional survivor of the three.
- All three are `ITSNdb_main` / `HLA-A*02:01` rows, so the analysis is really about source-heldout and HLA-heldout survival, not about removing overlap via trivial memorization.

## Output files

- `/home/seungho/personal/THCA_data_analysis/project/results/clean_neobench_barneo_2026_05_09/ga_rl_t1_failure_heatmap.tsv`
- `/home/seungho/personal/THCA_data_analysis/project/results/clean_neobench_barneo_2026_05_09/ga_rl_t1_nearest_neighbors.tsv`
- `/home/seungho/personal/THCA_data_analysis/project/results/clean_neobench_barneo_2026_05_09/ga_rl_t1_distribution_summary.tsv`
- `/home/seungho/personal/THCA_data_analysis/project/results/clean_neobench_barneo_2026_05_09/figures/fig_ga_rl_t1_failure_heatmap.png`

## Limitations

- The heatmap is a retrospective audit, not external validation.
- It does not support any clinical vaccine recommendation.
- It is intentionally narrow: only the 3 unique GA/RL T1 candidates are shown.
