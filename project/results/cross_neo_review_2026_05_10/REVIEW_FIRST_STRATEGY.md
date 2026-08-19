# Review-First Strategy

## Working Title
From neoantigen discovery to vaccine-ready prioritization: why benchmark contamination, source shift, and top-k utility now matter more than AUROC

## Core Thesis
Personalized neoantigen vaccines are clinically credible again, but the computational bottleneck has shifted. The field no longer only needs better MHC binding prediction. It needs contamination-controlled, source-aware, patient-context-aware ranking systems that optimize the handful of candidates actually manufactured or assayed.

## Why Review First
1. It lets us define the problem before presenting CROSS-Neo.
2. It makes public-overlap and benchmark leakage a field-wide issue, not a defensive caveat about our model.
3. It positions top-k precision, AUPRC, OOD behavior, and abstention as practical clinical/product metrics.
4. It creates the conceptual bridge for the original CROSS-Neo paper: a contamination-controlled pan-allele prioritizer rather than another binding predictor.

## Target Journal Shape
- Fast review / perspective: Nature Reviews Clinical Oncology, Nature Reviews Immunology, Cancer Discovery review, Trends in Cancer, JITC review, Frontiers if speed matters.
- If business timing dominates: preprint + white paper first, then invited-style review.

## Do Not Claim
- Do not claim external validation for CROSS-Neo.
- Do not claim quantum advantage.
- Do not claim public predictors are invalid; say they are powerful but benchmark-overlap-sensitive.
- Do not claim neoantigen vaccines already have universal clinical efficacy.
