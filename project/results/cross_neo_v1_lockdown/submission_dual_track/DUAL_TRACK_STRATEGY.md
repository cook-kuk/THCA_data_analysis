# CROSS-Neo Dual-Track Submission Strategy

Status: planning scaffold, not manuscript prose.

## Core Position

Training time is not the bottleneck in neoantigen immunogenicity prediction. The bottleneck is whether a ranking signal survives:

- exact peptide-HLA holdout
- near peptide cluster holdout
- HLA allele/supertype heldout
- source/study heldout
- public-corpus overlap audit
- top-k precision/enrichment rather than AUROC-only reporting

Therefore, CROSS-Neo should be framed around reviewer-defensible evaluation and locked ranking behavior, not fast training.

## Track A: NeurIPS Evaluations & Datasets / Evaluation Science

Best fit if the main contribution is:

> Neoantigen immunogenicity predictors can look strong under conventional internal evaluation but become unstable under contamination-controlled retrieval, public-corpus overlap, HLA, near-neighbor, and source-shift stress tests.

Primary artifact:

- leakage-safe evaluation protocol
- locked split definitions
- public overlap manifest
- QK rescue/harm analysis
- source-heldout collapse diagnosis
- benchmark cards / evaluation cards

Claim boundary:

- Yes: internal locked evaluation framework
- Yes: source-shift stress test
- Yes: top-k/enrichment-first reporting
- No: external validation
- No: quantum advantage
- No: clinical vaccine selection predictor

Most defensible title:

> When Neoantigen Predictors Fail: Leakage-Safe Evaluation and Source-Shift Stress Testing for Immunogenicity Prioritization

## Track B: Bio / AI-for-Science Workshop or Use-Inspired Main Track

Best fit if the main contribution is:

> A contamination-controlled prioritization framework that combines mutant-WT counterfactual encoding with fold-safe QK fallback and explicitly reports where it fails.

Primary artifact:

- CROSS-Neo v1 method
- locked anchor/fusion comparison
- top-k prioritization behavior
- HLA/near-cluster robustness
- source-heldout failure analysis
- biological interpretation of rescued/harmed cases

Claim boundary:

- Yes: stress-tested candidate prioritizer
- Yes: fold-safe fusion improves internal/HLA/near-neighbor ranking
- Yes: source-shift limited
- No: SOTA external predictor
- No: clinical deployment

Most defensible title:

> CROSS-Neo: Contamination-Controlled Counterfactual Ranking for Neoantigen Immunogenicity Prioritization

## Current Decision

Prepare both tracks from the same locked result set.

The ED/evaluation paper is the stronger top-tier route right now because source-heldout collapse becomes evidence, not a weakness. The bio/workshop route is still useful if we add high-compute experiments and a better source-shift mitigation story.

