# CROSS-Neo Two-Paper High-Impact Strategy

## Decision

Yes: this should be **two papers**, not one.

The impact goes up because the papers do different jobs:

1. **Paper 1 defines the field standard.**
2. **Paper 2 implements the standard.**

Trying to put both into one paper weakens both: the Review becomes self-promotional, and the Methods paper becomes underpowered as a universal model claim. Separating them makes the story cleaner and higher-impact.

## Paper 1: Standards Perspective

### Role
Set the field agenda.

### Best Title
**Vaccine-ready AI for personalized cancer immunotherapy**

### Stronger Provocative Title
**Neoantigen prediction needs vaccine-ready benchmarks, not another AUROC leaderboard**

### Core Contribution
NEO-PRIOR: a minimum reporting and evaluation standard for neoantigen prioritization.

### Target Order
1. Nature Medicine Perspective / Comment
2. Nature Cancer Perspective
3. Nature Reviews Immunology
4. Nature Reviews Clinical Oncology
5. Nature Biomedical Engineering Perspective
6. Cancer Cell / Cancer Discovery

### What It Must Not Do
- Do not center CROSS-Neo.
- Do not present unpublished model results as evidence.
- Do not claim public predictors are bad.
- Do not claim external validation.

### Winning Sentence
The next advance in personalized cancer vaccines may come less from predicting more HLA binders than from trusting the shortlist that decides what gets manufactured.

## Paper 2: Benchmark / Validation Implementation

### Role
Show that the standard can be operationalized.

### Best Title
**Benchmarking vaccine-ready AI systems for neoantigen prioritization**

### Alternative Title
**CROSS-Neo: contamination-controlled pan-allele prioritization for personalized neoantigen vaccine triage**

### Core Contribution
A leakage-controlled benchmark and prioritization stack with explicit public-comparator separation, top-k primary metrics, source stress tests and OOD/abstention boundaries.

### Target Order
1. Nature Biomedical Engineering Article
2. Nature Machine Intelligence Analysis
3. Cell Reports Medicine
4. Patterns
5. JITC
6. npj Precision Oncology

### What It Must Not Do
- Do not claim external validation.
- Do not claim quantum advantage.
- Do not hide source-heldout failures.
- Do not present product-demo testset-aware scores as clean generalization.

## Why This Raises Impact

| Problem | One-paper approach | Two-paper approach |
|---|---|---|
| Review seems self-serving | high risk | avoided |
| Method seems underpowered | high risk | framed as implementation/stress-test |
| Public predictor criticism | defensive | becomes field-standard discussion |
| Small-n issue | central weakness | handled through benchmark/abstention framing |
| Business value | product claim only | category leadership + demo |

## Launch Order

1. Send Paper 1 presubmission to Nature Medicine and Nature Cancer.
2. Use the response to sharpen the standards framing.
3. In parallel finalize Paper 2 benchmark package.
4. Submit Paper 2 only after Paper 1 pitch is in motion, so CROSS-Neo can be framed as implementing NEO-PRIOR.

## Shared Claim Boundary

- Review/Perspective: field standard, no new original model claims.
- Methods/Benchmark: internal locked benchmark and implementation, no external validation.
- Product/demo: retrospective triage artifact, not manuscript validation.

## Impact Map

`/home/seungho/personal/THCA_data_analysis/project/results/cross_neo_two_paper_strategy_2026_05_10/figures/two_paper_impact_map.png`
