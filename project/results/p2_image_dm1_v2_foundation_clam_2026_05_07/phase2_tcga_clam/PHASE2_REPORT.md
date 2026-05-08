# Phase 2 — CLAM attention-MIL (UNI features) on TCGA-THCA

- 59 slides (29 DM1 / 30 DM2)
- 5-fold CV mean AUC: 0.830 ± 0.139
- Per-fold: [0.7142857142857143, 0.7222222222222222, 1.0, 0.7142857142857143, 1.0]
- Pooled cross-fold AUC: 0.746

- Closure ResNet50 baseline AUC ~ 0.55 (negative)
- Kill-switch: PASS if held-out AUC > 0.70; MARGINAL 0.60-0.70; FAIL ≤ 0.60

## VERDICT

**PASS** — foundation+CLAM exceeds closure ResNet50; Paper 2 launch unlocked