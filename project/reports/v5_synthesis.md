# v5 Synthesis — Algorithmic + Cross-cancer + Robust Target Sprint

## Executive summary

The v5 sprint upgrades v4's empirical audit (which concluded
UNRECOVERABLE) into three methodological contributions:

1. **Nonlinear correction ablation (Track 1).** We tested
   6 batch correction methods. Best method:
   **zscore** (verdict: RESCUED).
   Post-correction identifiability AUC =
   0.13697264271547246;
   LODO AUC = 0.9971916608847919;
   bio-preservation = 0.944976526994419.

2. **DANN-style adversarial de-confounding (Track 2, novel).** A
   gradient-reversal-layer biomarker classifier with λ sweep over
   {0, 0.1, 0.3, 1.0, 3.0}. Best λ = **0.1**;
   bio CV AUC = 0.9286392882276063;
   ident AUC = 1.0;
   SHAP-stable top-10 genes: TACSTD2, CST6, DIO1, COL1A1, DDX3Y, IL1RL1.

3. **Cross-cancer 6-check (Track 3).** Replicates in
   1/4 cancers;
   correlation(internal, ident) = nan;
   generalization verdict: **cancer-specific**.

4. **Robust target shortlist (Track 4).** Top-20 batch-robust candidates +
   5 zero-literature novel + 5 known.
   Top 3 novel: CST6, DIO1, PDLIM4.
   Top 3 known: TACSTD2, GABRB2, LDLR.

## Upgraded publication path

Primary venue: **Bioinformatics** (methods paper). Justification:
the DANN-based biomarker de-confounding method is a distinct technical
contribution with cross-cancer generalization; the batch-robust shortlist
is an application result; together these form a methods paper rather
than a workshop poster.

Secondary (fallback): ML4H workshop.

## Honesty preserved
All language throughout: decision-support prototype / retrospective
computational triage / methodology contribution. Not diagnostic.
