# NeoBench-Vax Leaderboard and Governance

## Purpose
NeoBench-Vax is a benchmark infrastructure for vaccine-ready neoantigen prioritization. It is designed to make models auditable under NEO-PRIOR rather than to reward a single headline AUROC.

## Tracks

| Track | Name | Public predictor scores as features | Intended claim |
|---|---|---:|---|
| A | Clean Comparator | No | Scientific comparison under strict leakage control |
| B | Product-Assisted Triage | Yes, disclosed | Retrospective decision-support utility only |
| C | Source-Stress | No or disclosed by subtrack | Robustness under source/study/prevalence shift |
| D | OOD/Abstention | No or disclosed by subtrack | Reliability and safe non-decision behavior |

## Ranking Metrics
Primary leaderboard rank should be multi-objective:

1. AUPRC over prevalence baseline.
2. top10 precision and enrichment@10.
3. top5 precision for manufacturing-scarce use.
4. retrieval-clean-only AUPRC and top-k.
5. source-heldout collapse penalty.
6. calibration and abstention benefit.

## Model Card Requirements
Every submission should state:

- allowed inputs;
- forbidden inputs;
- public predictor use;
- split contract;
- training data provenance;
- known overlap risks;
- source/HLA/time generalization boundary;
- whether negative labels are hard negatives or unlabeled.

## Governance Principle
A model can lead a product-assisted track and still be ineligible for clean scientific claims. This separation is central to the benchmark.
