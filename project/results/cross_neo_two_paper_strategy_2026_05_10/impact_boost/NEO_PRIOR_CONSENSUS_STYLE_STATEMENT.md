# NEO-PRIOR Consensus-Style Statement

## Full Name
NEO-PRIOR: Neoantigen Prioritization Reporting and Evaluation Standard

## Purpose
NEO-PRIOR defines minimum reporting requirements for computational systems that rank neoantigen candidates for personalized vaccine manufacture, immune testing or clinical prioritization.

## Scope
NEO-PRIOR applies to:
- binding/presentation predictors used in vaccine pipelines;
- immunogenicity predictors;
- pan-allele neoantigen rankers;
- public benchmark papers;
- clinical vaccine candidate-selection pipelines;
- product/demo triage systems.

## Minimum Items

1. State the task: binding, presentation, immunogenicity, vaccine triage or benchmark comparison.
2. Report exact peptide-HLA and exact peptide overlap.
3. Report near-peptide and source-protein-window overlap.
4. Disclose public predictor scores used as features, comparators or product assists.
5. Use train-only retrieval, preprocessing, scaling, calibration and threshold selection.
6. Include exact, near, source/study, HLA allele/supertype and time/assay holdouts when feasible.
7. Treat negative labels as ambiguous unless true non-immunogenicity is experimentally established.
8. Report AUPRC, top-k precision, enrichment over prevalence and recall@k as primary metrics.
9. Report AUROC as secondary.
10. Report calibration, Brier/ECE and OOD/abstention coverage.
11. Publish data provenance and claim-boundary tables.
12. Separate clean scientific validation from product/demo retrospective triage.

## Why It Matters
Neoantigen vaccines depend on a short candidate list. The wrong benchmark can reward memorization, HLA shortcuts or source-specific bias, while the right benchmark can reveal when a system should abstain.

## Intended Use
NEO-PRIOR is a proposed reporting checklist, not a regulatory standard. It is designed to make neoantigen AI papers auditable, comparable and clinically interpretable.
