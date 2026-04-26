# v5.1 Change Log

_Date: 2026-04-24 · Author: THYRAI Research · Contact: kukshomr@gmail.com_

## What changed from v5 to v5.1

The v5 DIAL cross-cancer audit has been retracted and superseded by v5.1. v5 reported a cross-cancer-general batch-entanglement regime based on StratifiedKFold cross-validation over a pooled matrix that mixed multiple cohorts, with semi-synthetic augmentations added where cohort coverage was thin. Both choices were wrong for the claim being made. v5.1 re-ran the audit under leave-one-dataset-out (LODO) on real multi-cohort expression matrices only, and the conclusion changed qualitatively: the label-flip regime is not cross-cancer-general — it is specific to THCA BRAF-vs-RAS. Four of five THCA classifiers (LogReg_l2, LogReg_elasticnet, RandomForest, GradientBoosting) flip direction post-ComBat with DIAL in [0.323, 0.494]; all 20 non-THCA classifier-cancer combinations are null. This is a specificity finding, not a generalization failure of the DIAL method.

## StratifiedKFold -> LODO (why)

StratifiedKFold over a pooled cross-cohort matrix leaks cohort identity across folds. Pre-ComBat AUC inflates because the classifier learns cohort-specific expression offsets as a proxy for subtype; post-ComBat AUC then deflates disproportionately because ComBat removes the very axis the classifier was exploiting, and the "flip" is partly an artifact of that leak. LODO — train on n-1 cohorts, fit ComBat on training cohorts only, evaluate on the held-out cohort — is the correct regime for a cross-cohort claim. Under LODO the pre-ComBat AUC no longer benefits from cohort-identity leak, and any post-ComBat flip that remains reflects a genuine property of the data generating process. In THCA the flip survives LODO; in the other four cancers it does not appear at all.

## Real data vs semi-synthetic

v5 augmented thin cohorts with semi-synthetic samples to reach sample-size thresholds. This conflated real biology with the generating process of the augmentation and invalidated inferences about batch-entanglement in those cohorts. v5.1 uses only real cohort expression matrices. Sample counts are modest (THCA 392, SKCM 212, LGG 142, LUAD 236, COAD 372), and the shared-gene intersection ranges from 9,179 (SKCM) to 20,964 (LGG) before variance filtering. We train classifiers on the variance-top-3,000 genes per cancer.

## Interpretation rule update (dial>0.3 primary, no ident gate under LODO)

v5 required both DIAL > 0.3 AND a low post-correction batch identifiability AUC to call a pair batch_entangled. Under LODO, post-correction batch identifiability is ill-defined (the held-out cohort is by construction the identifiability target), so the ident gate is dropped. The v5.1 rule is: DIAL > 0.3 -> batch_entangled; DIAL <= 0.1 and auc_post > 0.7 -> true_biology; DIAL <= 0.1 and auc_post < 0.6 -> no_signal; otherwise ambiguous.

## Retracted claims vs Valid v5.1 claims

| Claim | v5 (retracted) | v5.1 (valid) |
| --- | --- | --- |
| Cross-cancer generalization of batch entanglement | "cross-cancer general label-flip regime under ComBat" | **Retracted.** THCA-specific only. 4/25 pairs flip, all in THCA. |
| Best correction method (zscore + DANN rescue) | post-ident 0.137, LODO 0.997 reported as RESCUED | **Retracted.** Rescue was inflated by StratifiedKFold leak; not reported in v5.1. |
| DIAL as diagnostic probe | Generic diagnostic usable on any cancer | Still valid. DIAL correctly fires on THCA and correctly stays silent on SKCM/LGG/LUAD/COAD. |
| Semi-synthetic augmentation | Used to boost cohort N | **Retracted.** v5.1 uses real data only. |
| THCA BRAF/RAS flip magnitude | DIAL ~= 0.31 (LogReg_l2) under StratifiedKFold | **Revised.** DIAL = 0.494 (LogReg_l2) under LODO on real data; cleaner complete-flip geometry. |

All v5.1 findings are computational. Prospective cohort validation is pending through ongoing collaboration with Prof. 유형원 at Seoul National University Bundang Hospital (유형원 교수 분당서울대 협력 진행 중).
