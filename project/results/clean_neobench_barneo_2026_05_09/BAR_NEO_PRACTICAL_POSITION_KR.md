# BAR-Neo Practical Position Memo

Date: 2026-05-09

## 한 줄 결론

현재 우리 알고리즘의 안전한 위치는 **새 SOTA neoantigen predictor**가 아니라, 여러 predictor를 agent가 호출한 뒤 CLEAN-NeoBench 성능/누수/shift/캘리브레이션 이력으로 가중치를 조절하는 **benchmark-adaptive reliability ranking controller**이다.

## 실용 구조

1. **CLEAN-NeoBench**
   - 후보, label, HLA, source, overlap flag, split contract를 통일한다.
   - exact peptide-HLA, near peptide, source-heldout, HLA-heldout, patient/study-heldout 조건에서 predictor를 평가한다.
   - public pretrained tool은 training-corpus overlap audit 전까지 caveated comparator로만 둔다.

2. **BAR-Neo**
   - 후보 intrinsic feature, method rank percentile, disagreement, overlap/OOD flag를 사용해 reliability score를 만든다.
   - confidence와 abstention reason을 같이 낸다.
   - high leakage, source shift, low metadata, public-only support는 낮은 confidence 또는 abstention으로 처리한다.

3. **BAR-Neo-BMA**
   - 실용 agent layer이다.
   - agent가 Structure_LR, RF/PU/source-balanced models, stacked models, QK fallback, uncertainty branch 등을 expert처럼 호출한다.
   - CLEAN-NeoBench leaderboard와 split robustness에서 posterior method weight를 계산한다.
   - QUBO-style sparse selector로 중복 predictor를 줄이고 context별 expert subset을 고른다.
   - score가 높아도 leakage risk나 expert support 부족이면 abstain한다.

## 현재 score table 위치

- 후보별 실용 ranking: `barneo_bma_candidate_scores.tsv`
- method posterior weight: `barneo_bma_method_weights.tsv`
- selector audit: `barneo_bma_selector_audit.tsv`
- 기본 BAR-Neo 후보 점수: `barneo_candidate_scores.tsv`
- benchmark leaderboard: `clean_neobench_leaderboard.tsv`
- split별 metric: `clean_neobench_split_metrics.tsv`

## 우리 알고리즘 위치

| Layer | 현재 위치 | reviewer-safe claim |
|---|---|---|
| Structure_LR | local anchor | honest internal anchor |
| RF / PU / source-balanced models | strong internal candidates | benchmarked internal comparators |
| W7A / W7B / stacked models | fusion candidates | internal candidate ensemble |
| ESM2 Bayesian | uncertainty/OOD branch | uncertainty-only support |
| QK branch | bounded fallback / selector probe | bounded fallback/fusion component |
| BAR-Neo | reliability ranker | benchmark-adaptive reliability ranking |
| BAR-Neo-BMA | practical agentic ensemble controller | calibrated prioritization with abstention |

## 현재 상위 method weight 해석

`barneo_bma_method_weights.tsv` 기준 최상위권은:

1. `Structure_LR`: anchor라서 AUPRC만 최고가 아니어도 posterior weight 1위.
2. `source_balanced_rf_train_prior_calibrated`: leaderboard상 mean AUPRC 최상위.
3. `source_balanced_plus_pu_rf_train_prior_calibrated`
4. `pu_weighted_rf_train_prior_calibrated`
5. `anchor_lr`
6. `W7B_stacked`
7. `W7A_full`

즉 실용 버전은 “가장 높은 단일 모델 하나”가 아니라, **anchor + calibrated internal models + fusion models를 reliability weight로 섞는 구조**가 맞다.

## 중요한 caveat

- top candidate 일부는 ensemble raw score가 높아도 `High leakage risk invalidates clean benchmark claim` 때문에 low confidence / abstain으로 남는다.
- 현재 BMA 결과는 2,715 후보 중 2,707개가 abstain이다. 이건 실패가 아니라 reviewer-safe mode가 강하게 켜진 상태다.
- fallback-only 후보는 score 상위권을 먹지 않도록 cap을 걸었다.
- clinical vaccine selection이 아니라 research triage이다.

## 다음 실용 개발 순서

1. Public tool training-corpus overlap audit.
2. PAAD/THCA patient-gated demo용 metadata 보강.
3. source-heldout / HLA-heldout collapse를 줄이는 calibration 개선.
4. MHC-II는 별도 benchmark로 분리.
5. QK는 quantum advantage 주장이 아니라 sparse feature/expert selection fallback으로만 유지.
