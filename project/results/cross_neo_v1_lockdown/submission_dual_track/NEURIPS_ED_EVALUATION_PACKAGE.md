# NeurIPS ED / Evaluation Paper Package

Status: scaffold only.

## One-Sentence Thesis

Existing neoantigen immunogenicity predictor comparisons are vulnerable to hidden retrieval overlap, public-corpus contamination, HLA imbalance, near-neighbor leakage, and source-label shift; CROSS-Neo provides a locked evaluation framework that separates real ranking signal from evaluation artifacts.

## Contribution Checklist

- Locked class-I strict set with explicit prevalence.
- Split-safe retrieval audit.
- Exact peptide-HLA and near-peptide cluster holdouts.
- HLA allele/supertype heldouts.
- Source-heldout stress across CEDAR, NEPdb, TESLA_mmc4, TESLA_mmc7.
- Public overlap manifest for MHCflurry, NetMHCpan, BigMHC, PRIME, MixMHCpred, NetMHCstabpan, IEDB, CEDAR, NEPdb, TESLA.
- QK rescue/harm case analysis.
- Top-k precision, enrichment, and AUPRC as primary metrics.
- Explicit forbidden-claim table.

## Result Anchors

Primary locked improvements over `anchor_rf`:

| Split | Best fold-safe method | AUPRC | AUROC | top10 |
|---|---:|---:|---:|---:|
| exact peptide-HLA holdout | prespecified_rf_qk_no_anchor_w0.5 | 0.587 | 0.757 | 0.700 |
| near peptide cluster holdout | rule_gate_rf_qk_fallback_train_selected | 0.519 | 0.723 | 0.700 |
| HLA stratified group | nested_rf_qk_quantum_only_train_selected | 0.555 | 0.761 | 0.600 |
| HLA supertype heldout | prespecified_lr_qk_quantum_only_w0.5 | 0.545 | 0.767 | 0.600 |

Failure evidence:

| Source-heldout | Key finding |
|---|---|
| CEDAR | high prevalence makes top-k less informative |
| NEPdb | QK compact partially rescues top10 to 0.4, but anchor/fusion collapse |
| TESLA_mmc4 | top10 remains 0 across rescue variants |
| TESLA_mmc7 | top10 remains 0 across rescue variants |

## Paper Shape

1. Problem: immunogenicity benchmarks can reward leakage and source shortcuts.
2. Audit design: locked splits, public overlap manifest, retrieval flags.
3. Methods: anchors, QK fallback, fold-safe fusion, source diagnostics.
4. Results: internal/HLA improvements.
5. Failure analysis: source-heldout collapse.
6. Evaluation recommendations: required split/reporting contract.
7. Limitations: small strict set, unresolved public overlap, no external validation.

## Reviewer-Safe Claim

CROSS-Neo v1 is an internal locked-split evaluation framework for neoantigen immunogenicity prioritization. It shows that fold-safe counterfactual/QK fusion improves internal/HLA/near-neighbor ranking, while source-heldout stress exposes major source-shift limitations that should prevent external-validation or clinical-use claims.

## Must Not Claim

- external validation
- quantum advantage
- clean public comparator status before row-level public overlap audit
- clinical vaccine selection readiness
- SOTA predictor status based only on internal splits

## Missing Before Submission

- Finish public training-corpus downloads/parsing for MHCflurry, NetMHCpan, BigMHC, PRIME, MixMHCpred, NetMHCstabpan.
- Convert overlap manifest into exact row-level audit.
- Add evaluation cards and dataset cards.
- Add bootstrap CIs to all headline tables.
- Package code/reproduction script with deterministic runner.

