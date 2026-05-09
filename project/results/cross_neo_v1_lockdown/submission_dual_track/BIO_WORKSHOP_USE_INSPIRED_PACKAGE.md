# Bio / AI-for-Science Workshop Package

Status: scaffold only.

## One-Sentence Thesis

CROSS-Neo is a contamination-controlled neoantigen prioritization framework that combines mutant-WT counterfactual encoding with fold-safe QK fallback, improves internal/HLA/near-neighbor top-k ranking, and explicitly identifies source-shift failure modes.

## Target Venue Shape

Best near-term route:

- NeurIPS AI-for-science / ML-for-biology workshop
- ML for healthcare / computational biology workshop
- Bioinformatics-method workshop track

Longer route:

- NeurIPS main Use-Inspired only if method novelty and cross-dataset robustness become stronger.
- Bioinformatics or Genome Medicine if biological validation/benchmark completion becomes stronger.

## Method Story

Feature groups:

- mutant peptide
- WT peptide if available
- HLA allele/supertype
- frozen counterfactual peptide features
- pMHC geometry features where available
- fixed QK fallback
- OOD/source-shift flags

Locked candidate:

- conservative: `nested_rf_qk_quantum_only_train_selected`
- exact pHLA strongest prespecified: `prespecified_rf_qk_no_anchor_w0.5`
- near-cluster strongest: `rule_gate_rf_qk_fallback_train_selected`

## Result Story

The method improves where a realistic model should improve:

- exact pHLA holdout: top10 0.7
- near-cluster holdout: top10 0.7
- HLA group: top10 0.6
- HLA supertype: top10 0.6

The method fails where the field is genuinely hard:

- NEPdb/TESLA source shift
- low-prevalence top-k brittleness
- label-definition mismatch
- missing WT/source-window features

This should be presented as honest stress testing, not as failure to hide.

## Figures

1. Locked split schema and leakage barriers.
2. Anchor vs fold-safe fusion AUPRC/top-k.
3. QK rescue/harm cases.
4. Source-heldout positive rank distributions.
5. Source-shift heatmap.
6. Public overlap status.
7. Top-k rescue/abstention curve.

## Workshop Abstract Skeleton

Do not use as final prose without author review.

Background:
Neoantigen immunogenicity prediction remains limited by small datasets, public-corpus overlap, HLA imbalance, and source-specific label definitions.

Method:
CROSS-Neo combines mutant-WT counterfactual peptide/HLA encoding with locked, fold-safe QK fallback and a contamination-controlled evaluation contract.

Findings:
Fold-safe fusion improves exact peptide-HLA, near-cluster, and HLA-heldout ranking, with top-10 precision up to 0.7 on locked internal splits. However, source-heldout NEPdb/TESLA stress tests expose collapse, showing that conventional internal metrics are insufficient.

Conclusion:
CROSS-Neo should be treated as a stress-tested prioritization candidate and evaluation framework, not an externally validated clinical predictor.

## Missing Before Workshop Submission

- Add one clean public-overlap audit table.
- Add per-case biological interpretation for top rescued positives and harmed negatives.
- Decide whether QK stays in title; safer to keep it in method, not headline.
- Add reproducibility package and exact command log.

