# Top Quantum Candidate Stress Tests

These are still internal strict-set tests. They are meant to identify failure modes before external benchmarking.

| feature_set | split | n / pos | AUROC | AUPRC |
|---|---|---:|---:|---:|
| quantum_structure_tcr | hla_stratified_group_5fold | 89 / 21 | 0.749 | 0.462 |
| quantum_tcr_no_structure | hla_stratified_group_5fold | 89 / 21 | 0.707 | 0.489 |
| structure_only | hla_stratified_group_5fold | 89 / 21 | 0.674 | 0.386 |
| quantum_structure_no_tcr | hla_stratified_group_5fold | 89 / 21 | 0.659 | 0.361 |
| quantum_only | hla_stratified_group_5fold | 89 / 21 | 0.628 | 0.425 |
| best_single_qk | hla_stratified_group_5fold | 89 / 21 | 0.480 | 0.243 |
| quantum_structure_tcr | repeated_stratified_5x5 | 89 / 21 | 0.785 | 0.603 |
| quantum_tcr_no_structure | repeated_stratified_5x5 | 89 / 21 | 0.716 | 0.523 |
| structure_only | repeated_stratified_5x5 | 89 / 21 | 0.707 | 0.483 |
| quantum_structure_no_tcr | repeated_stratified_5x5 | 89 / 21 | 0.687 | 0.497 |
| quantum_only | repeated_stratified_5x5 | 89 / 21 | 0.650 | 0.426 |
| best_single_qk | repeated_stratified_5x5 | 89 / 21 | 0.518 | 0.294 |

Paper boundary: promote only if the signal survives fixed-feature, fixed-gamma, HLA/study/time-held-out external splits.
