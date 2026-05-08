# Phase 2 — CLAM attention-MIL (UNI features) on TCGA-THCA

- 18 slides (8 DM1 / 10 DM2)
- 5-fold CV mean AUC: 0.833 ± 0.211
- Per-fold: [1.0, 1.0, 0.6666666666666667, 1.0, 0.5]
- Pooled cross-fold AUC: 0.763

- Closure ResNet50 baseline AUC ~ 0.55 (negative)
- Kill-switch: PASS if held-out AUC > 0.70; MARGINAL 0.60-0.70; FAIL ≤ 0.60

## VERDICT

**PASS** — foundation+CLAM exceeds closure ResNet50; Paper 2 launch unlocked