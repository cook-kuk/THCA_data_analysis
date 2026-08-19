# DL Funnel Threshold Optimization

## Purpose

This searches slider thresholds against the available labels to create operating presets instead of stopping at manual threshold tuning.

## Recommended Presets

| preset_name              | call_rule              |   called_positive |   TP |   TN |   FP |   FN |   precision |    recall |        F1 |   specificity |      FPR |
|:-------------------------|:-----------------------|------------------:|-----:|-----:|-----:|-----:|------------:|----------:|----------:|--------------:|---------:|
| NO_FALSE_POSITIVE_MAX_TP | md_escalation          |                13 |   13 |  488 |    0 |  148 |    1        | 0.0807453 | 0.149425  |      1        | 0        |
| HIGH_PRECISION_MIN_FP    | md_escalation          |                13 |   13 |  488 |    0 |  148 |    1        | 0.0807453 | 0.149425  |      1        | 0        |
| BALANCED_F1              | non_culled             |               347 |  133 |  274 |  214 |   28 |    0.383285 | 0.826087  | 0.523622  |      0.561475 | 0.438525 |
| RECALL_PRESERVING        | non_culled             |               112 |   63 |  439 |   49 |   98 |    0.5625   | 0.391304  | 0.461538  |      0.89959  | 0.10041  |
| WETLAB_ULTRA_STRICT      | wetlab_shortlist       |                 2 |    2 |  488 |    0 |  159 |    1        | 0.0124224 | 0.0245399 |      1        | 0        |
| STRUCTURE_MD_STRICT      | structure_md_supported |                 2 |    2 |  488 |    0 |  159 |    1        | 0.0124224 | 0.0245399 |      1        | 0        |

## Claim Boundary

These presets are optimized on the current labeled candidate table and must be treated as decision-support settings, not externally validated clinical thresholds.
