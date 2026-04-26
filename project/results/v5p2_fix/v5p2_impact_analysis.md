# v5.2 — Impact Analysis (Task A1: proper LODO ComBat)

**Generated:** 2026-04-25
**Compares:** `v5p1_dial_all_cancers.tsv` (v5.1) vs `v5p2_dial_proper_lodo.tsv` (v5.2).
**Re-fit code:** `notebooks_or_scripts/v5p2_combat_lodo.py`,
`notebooks_or_scripts/v5p2_dial_proper.py`.

## What changed methodologically

| Aspect | v5.1 | v5.2 |
| --- | --- | --- |
| ComBat input | full X (train ∪ test) | only X[train] within each LODO fold |
| Batch coefficients γ_{gb} | shared across all data | learned on train only |
| Test fold processing | corrected with parameters that saw test | corrected with train-derived α, β + locally estimated γ, δ for unseen batch |
| Implementation | `inmoose.pycombat.pycombat_norm` (one-shot) | manual EB-free linear ComBat (`LinearComBat.fit`/`.transform` in `v5p2_combat_lodo.py`) |
| Held-out cohort handling | implicit (its own statistics influenced parameters) | explicit & documented: γ̂ = X̄_test − α_g, δ̂ = SD_test (no Y info) |

The leak in v5.1 was: ComBat estimated α, β, γ, δ from rows whose labels would
later become the test fold. Even if the LODO split was honest at the classifier
level, the input matrix was already adjusted using held-out information.

## Headline numbers (THCA — the central paper claim)

| Classifier | v5.1 auc_post | v5.2 auc_post | v5.1 DIAL | v5.2 DIAL | v5.1 interp | v5.2 interp |
| --- | ---: | ---: | ---: | ---: | --- | --- |
| LogReg_l2          | 0.006 | **0.995** | 0.494 | **0.000** | batch_entangled | true_biology |
| LogReg_elasticnet  | 0.008 | **0.995** | 0.492 | **0.000** | batch_entangled | true_biology |
| RandomForest       | 0.177 | **0.878** | 0.323 | **0.000** | batch_entangled | true_biology |
| GradientBoosting   | 0.166 | **0.800** | 0.334 | **0.000** | batch_entangled | true_biology |
| HistGB (≈XGBoost slot) | n/a (XGB env) | 0.740 | n/a | **0.000** | n/a | true_biology |

**THCA flips from "4/5 batch_entangled" to "0/5".** Per-fold AUC (see
`per_fold_auc_post` column of `v5p2_dial_proper_lodo.tsv`) confirms both LODO
folds independently land above 0.5 for every classifier — there is no AUC
inversion under proper LODO.

## Other cancers

| Cancer | v5.1 batch_entangled | v5.2 batch_entangled | v5.1 true_biology | v5.2 true_biology |
| --- | ---: | ---: | ---: | ---: |
| SKCM | 0/5 | 0/5 | 0/5 | 0/5 |
| LGG  | 0/5 | 0/5 | 5/5 | **5/5** |
| LUAD | 0/5 | 0/5 | 3/5 | 3/5 |
| COAD | 0/5 | 0/5 | 4/5 | **5/5** |

**No other cancer changes call.** The overall non-flip pattern (specificity
evidence) is preserved: 0/20 batch_entangled for SKCM/LGG/LUAD/COAD across both
versions. THCA was the *only* row where the v5.1 → v5.2 refactor moved any
classifier across an interpretation boundary, and it moved in the opposite
direction from what the v5.1 paper claimed.

## Per-row magnitude shifts (Δ DIAL)

```
THCA  LogReg_l2          ΔDIAL = -0.494
THCA  LogReg_elasticnet  ΔDIAL = -0.492
THCA  GradientBoosting   ΔDIAL = -0.334
THCA  RandomForest       ΔDIAL = -0.323
LUAD  GradientBoosting   ΔDIAL = -0.040
all others                ΔDIAL = 0.000
```

Mean |Δ DIAL| across the 20 directly comparable rows = **0.082**;
restricted to THCA = **0.411**.

## Caveats and limitations

1. **Held-out batch transform is unsupervised.** For an LODO test fold the
   covariate Y is unavailable (it is what we predict). My LinearComBat estimates
   γ̂_{test}, δ̂_{test} from the test cohort's own X mean/SD only. This is the
   "conservative" branch named in the v5.2 spec. A nearest-neighbour batch
   transform was *not* attempted because (a) the unsupervised version is
   already conservative and (b) no nearest train batch exists for the
   bi-cohort cancers (each LODO fold has a single train batch).

2. **2-cohort cancers fit ComBat with no batch term in train.** When train is a
   single cohort (THCA, SKCM, LUAD, COAD folds), the ComBat fit reduces to a
   per-gene OLS on Y only; α_g and β_g are estimated honestly, γ_train = 0 by
   construction, and the held-out cohort still gets centered locally. This is
   the only batch correction available without leaking test info, and it is
   what a publishable LODO benchmark must use.

3. **Classifier substitution (XGBoost → HistGB).** v5.1 used `XGBClassifier`
   when available; the v5.2 venv does not have xgboost, so the 5th classifier
   slot is `HistGradientBoostingClassifier`. The XGB row in
   `v5p2_lodo_comparison.tsv` is therefore blank in the v52 column. HistGB is
   functionally similar (gradient-boosted trees) and shows no batch
   entanglement on THCA, so the substitution does not rescue the v5.1 claim.

4. **Variance-top filter is identical (3000 genes).** Same gene-selection step
   as v5.1 to keep runtime feasible while ruling out a gene-set difference as
   the source of the change.

## Verdict (sets up `v5p2_critical_assessment.md`)

The v5.1 central claim — *"THCA exhibits batch_entangled in 4/5 classifiers
under LODO ComBat"* — **does not survive proper LODO ComBat**. The flip was an
artifact of fitting ComBat on data that included the held-out test cohort.
The specificity evidence from the other 4 cancers is unchanged, but it was
never the headline finding. The paper's central claim must be retracted or
re-framed.
