# Phase 2 — CLAM attention-MIL (UNI features) on TCGA-THCA

- 54 slides (25 DM1 / 29 DM2)
- 5-fold CV mean AUC: 0.951 ± 0.081
- Per-fold: [1.0, 1.0, 0.9642857142857143, 1.0, 0.7916666666666667]
- Pooled cross-fold AUC: 0.874

- Closure ResNet50 baseline AUC ~ 0.55 (negative)
- Kill-switch: PASS if held-out AUC > 0.70; MARGINAL 0.60-0.70; FAIL ≤ 0.60

## VERDICT

**PASS** — foundation+CLAM exceeds closure ResNet50; Paper 2 launch unlocked