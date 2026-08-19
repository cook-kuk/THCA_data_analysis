# CROSS-Neo Methods/Benchmarking Manuscript Skeleton

## Working Title
Benchmarking vaccine-ready AI systems for neoantigen prioritization

## Claim Boundary
CROSS-Neo is an internal locked benchmark and prioritization framework. It does not claim external validation, clinical efficacy or quantum advantage.

## Abstract
Personalized neoantigen vaccines require computational systems that select a small number of candidates for manufacture or immune testing, but many predictors are evaluated as retrospective classifiers under overlap-sensitive public benchmarks. We introduce CROSS-Neo, a contamination-controlled benchmark and pan-allele prioritization stack for vaccine-ready neoantigen ranking. The system combines mutant-WT counterfactual features, hard-decoy pressure, fixed sequence/QK fallback signals, explicit retrieval-risk flags, source-aware calibration and OOD-aware abstention. Evaluation uses locked split families, including exact peptide-HLA, near-peptide, HLA/supertype, source/study and public-overlap audits. Primary outcomes are AUPRC, top-k precision, enrichment over prevalence and abstention behavior. Internal retrospective benchmarks show improved practical top-k prioritization over public predictor-only baselines, while source-heldout stress tests identify contexts where abstention remains necessary.

## Results

### A leakage-safe benchmark registry
- master table, split manifests, public overlap audit, retrieval risk flags.

### Public predictor baselines are strong but overlap-sensitive
- BigMHC, MHCflurry, PRIME, DeepImmuno, TransPHLA, NetMHCpan, T-SCAPE, MHCnuggets, NetMHCstabpan.
- Public predictor scores are comparators, not clean training features.

### CROSS-Neo branches have complementary error profiles
- counterfactual RF, QK fixed fallback, structure/geometry, retrieval evidence.
- naive concatenation fails; gated/fixed fusion is safer.

### Hard-decoy and pan-allele fusion improve top-k triage
- product fixed pan-allele score improves over public-only baselines in internal strict demo set.

### Source-heldout stress reveals abstention boundaries
- NEPdb improves, TESLA remains difficult.
- This is a limitation and deployment guardrail, not a failure to hide.

### Reporting standard and model-card outputs
- NEO-PRIOR checklist, data provenance, claim boundary, candidate queue.

## Methods
- data sources and inclusion criteria;
- feature construction;
- split design;
- retrieval index safety;
- calibration and abstention;
- public predictor score handling;
- metrics and statistics;
- software and reproducibility.

## Display Items
See `methods_figure_catalog.tsv`.

## Tables
See `methods_table_catalog.tsv`.
