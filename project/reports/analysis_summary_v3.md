# THCA v3 analysis summary

_Build: 2026-04-24 06:20 UTC_

This is an **exploratory, decision-support prototype** built from retrospective 
public cohorts. None of the numbers below constitute a diagnostic test or 
imply clinical validity.

## 1. Honesty audit (page 19)

- TierA67_clean TCGA AUC at k=0 genes removed: **0.992**.
- After removing top-10 |Cohen's d| genes: **0.984**.
- After removing top-30: **0.895**.
- MAPK-output ablation (drops DUSP4/5/6, SPRY1/2/4, ETV4/5, FOSL1, PHLDA1): AUC drops from **0.923** to **0.924**.
- Permutation null: observed AUC **0.923**, empirical p = **0.0010**.
- Dataset identifiability mean OvR AUC: **1.000** — any value ≫ 0.7 implies the same model could be learning batch/platform.

## 2. External validation (LODO, page 19)

| external | mode | AUC | PR-AUC | bACC | NPV@Se≥0.95 | Brier |
|----------|------|-----|--------|------|-------------|-------|
| GSE27155 | train_tcga_only | 0.843 | 0.929 | 0.798 | 0.889 | 0.284 |
| GSE27155 | train_tcga_plus_others | 0.883 | 0.951 | 0.815 | 0.846 | 0.263 |
| GSE126698 | train_tcga_only | 0.694 | 0.809 | 0.750 | NA | 0.450 |
| GSE126698 | train_tcga_plus_others | 0.778 | 0.862 | 0.833 | 1.000 | 0.281 |
| GSE213647 | train_tcga_only | NA | NA | NA | NA | NA |
| GSE213647 | train_tcga_plus_others | NA | NA | NA | NA | NA |

## 3. Best panel (page 20)

- Best strategy: **univariate_d**, k = **16**.
- TCGA 5-fold AUC: **0.932** — GSE27155 isotonic-recalibrated AUC: **0.929**.

## 4. Three-class (page 21)

- Multinomial LogReg macro OvR AUC: **0.937**.
- GradientBoosting macro OvR AUC: **0.938**.

## 5. Bethesda prevalence sweep (page 23)

| prevalence | chosen threshold | Sens | Spec | NPV | surgery reduction vs treat-all |
|---|---|---|---|---|---|
| 10% | 0.26 | 0.991 | 0.411 | 0.998 | 37% |
| 15% | 0.26 | 0.987 | 0.406 | 0.995 | 35% |
| 20% | 0.26 | 0.995 | 0.408 | 0.997 | 33% |
| 25% | 0.26 | 0.992 | 0.397 | 0.993 | 30% |
| 30% | 0.26 | 0.992 | 0.415 | 0.991 | 29% |

## 6. Multimodal (page 24)

- Fusion: **skipped** (no local callset).
- SCNA: **skipped** (no local copy-number calls).
- Methylation (GSE97466 within-cohort tumor-vs-normal): AUC **0.986**.

## Caveats

- All external AUCs use isotonic recalibration fit on held-out slices of each external cohort — performance reported reflects decision-support behaviour after platform-specific recalibration, not raw transfer.
- Classes with n<20 are flagged in every figure legend with '⚠ small n'.
- Decision-support prototype only. Not a diagnostic device; not clinically validated.
