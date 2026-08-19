# NeoBench-Vax Benchmark Challenge Plan

## Purpose
NeoBench-Vax turns the two-paper package into reusable infrastructure.

## Assets
- canonical benchmark registry;
- public predictor score matrix;
- locked split manifests;
- exact/near/source/HLA/study/time holdouts;
- public overlap audit;
- model cards;
- data provenance cards;
- top-k and AUPRC reporting templates;
- failure-mode table.

## Challenge Tracks

### Track A: Clean Comparator Track
No public predictor scores as training features. Public predictors are comparators only.

### Track B: Product-Assisted Track
Public predictor scores allowed, but clearly labeled as assisted triage, not clean scientific validation.

### Track C: Source-Stress Track
Models are ranked by source-heldout top-k enrichment and abstention behavior.

### Track D: OOD/Abstention Track
Models are rewarded for abstaining when source/HLA/length/retrieval evidence is unsafe.

## Leaderboard Metrics
1. AUPRC
2. top5 precision
3. top10 precision
4. enrichment@10
5. recall@10
6. Brier/ECE
7. abstention coverage versus precision
8. source-heldout collapse penalty

## Why This Helps the Papers
Paper 1 proposes the reporting standard. NeoBench-Vax makes the standard concrete. Paper 2 shows CROSS-Neo as the first system evaluated through the benchmark.
