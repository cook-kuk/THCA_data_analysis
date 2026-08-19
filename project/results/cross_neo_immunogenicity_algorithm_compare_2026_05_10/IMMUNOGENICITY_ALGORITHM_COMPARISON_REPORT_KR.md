# CROSS-Neo immunogenicity algorithm comparison

Generated: 2026-05-10T22:30:42

## Bottom line

This run grafts an external immunogenicity algorithm lane into CROSS-Neo and compares it against the existing CROSS-Neo/Kaggle-style scores on the same retrospective known-answer label set.

- Candidates evaluated: 2,715
- Best top-96 score: `immunogenicity_discovery_score` with top96 precision 1.000 (96/96)
- BigMHC IM comparator: top96 precision 0.875, AUPRC 0.672, AUROC 0.722
- Existing `stress_guarded_discovery_score`: top96 precision 0.948, AUPRC 0.872, AUROC 0.873

## Interpretation

BigMHC IM is the strongest locally runnable public immunogenicity comparator because it models class-I presentation and then transfer-learns immune-response labels. It should be treated as an external comparator and discovery feature, not as a clean manuscript feature until public training-overlap is audited.

The new `immunogenicity_discovery_score` is a fixed-weight integration:

`0.30 CROSS-Neo + 0.25 BigMHC IM + 0.15 BigMHC EL + 0.15 TCR expert + 0.10 mutant/WT foreignness proxy + 0.05 MD/control readiness`

The companion `immunogenicity_claim_safe_score` applies existing validity/patient/overlap caps. It is expected to be more conservative than discovery ranking.

## Action counts

| immunogenicity_action               |    n |
|:------------------------------------|-----:|
| deprioritize_or_control_only        | 2346 |
| watchlist_or_active_learning        |  350 |
| wetlab_priority_with_claim_boundary |   19 |

## Claim boundary

These are retrospective known-label metrics and assay-prioritization scores. They do not establish prospective wetlab immunogenicity, clinical vaccine selection, or clean SOTA status. The next claim upgrade requires the preregistered 96-well assay result interpreter.
