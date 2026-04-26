# v5 Upgraded Paper Outline

## Primary venue: Bioinformatics (methods journal)

### Title (working)
"Cohort-invariant biomarker discovery in multi-study bulk RNA-seq via
adversarial domain-adaptation, with cross-cancer framework validation."

### Core novelty
1. **DANN for bulk-RNA-seq biomarker de-confounding.** Existing tools
   (ComBat, Harmony, MNN) correct upstream; we jointly train a biomarker
   classifier adversarially regularized against cohort identity, and
   report a principled Pareto frontier (bio AUC vs 1-ident AUC).
2. **6-check honesty audit framework.** A reusable set of checks
   (internal CV, identifiability, driver-ablation, leakage curve,
   permutation null, LODO) that separate "true signal" from
   "batch-entangled fit". We apply it to 4 cancers to demonstrate the
   fail pattern is general.
3. **Batch-robust target shortlist.** A composite score combining
   adversarial SHAP stability, correction-method survival, external
   replication, and structural druggability. Top-20 + 5 zero-lit novel.

### Figures (proposed)
- Fig 1: pipeline diagram + GRL architecture.
- Fig 2: nonlinear correction radar + PCA grid (Track 1).
- Fig 3: DANN Pareto frontier + SHAP stability (Track 2).
- Fig 4: cross-cancer 6-check grid + internal-vs-ident correlation (Track 3).
- Fig 5: robust-shortlist radar + literature-vs-score scatter (Track 4).
- Fig S1: v4 UNRECOVERABLE baseline (context).

### Ablations (planned)
- Lambda sweep {0, 0.1, 0.3, 1.0, 3.0} on real data.
- Drop-one-cohort ablation (does de-confounding survive when only 2 cohorts?).
- Compare vs. scVI / Harmony on the same biomarker task.

### Justification for venue upgrade (v4 → v5)
- v4 delivered an honest audit. That is useful but not a method.
- v5 delivers: (1) a novel training objective for bulk RNA-seq, (2) a
  reusable audit framework replicated across 4 cancer types, and
  (3) a deployable shortlist with scoring weights and a method card.
- Together these change the contribution class from "empirical audit
  of a single cancer" (workshop fit) to "general method + framework +
  deliverable shortlist" (methods journal fit).

## Secondary venue: ML4H workshop (fallback)
Shorter variant, focused on the DANN novelty only. 4 pages + supplement.

## Gap to JCO-PO / Nature Methods
- JCO-PO: need a prospective RNA-seq cohort + IRB data + decision-curve
  analysis on real clinical operating points. Out of v5 scope.
- Nature Methods: need ≥3 competitor methods on ≥3 real cancer datasets,
  plus a benchmarking harness. Estimated 2-3 months of compute + 2 months
  of writing.
