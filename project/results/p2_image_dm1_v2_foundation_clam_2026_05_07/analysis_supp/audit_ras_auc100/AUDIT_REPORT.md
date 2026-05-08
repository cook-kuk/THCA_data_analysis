# RAS-like AUC=1.000 audit — data-split / shortcut investigation

**Question:** Paper 2 image-DM1 CLAM (UNI features) reports RAS-like / FVPTC subgroup AUC=1.000 (n=16, n_pos=3 / n_neg=13). Is this real biological signal, a data-split artifact, or a histology shortcut?

## Verdict TL;DR

(see Section 6 — filled after retrain finishes)

## 1. Are RAS-like and FVPTC the same 16 slides?

| RAS_like | FVPTC | RAS ∩ FVPTC | union |
|---|---|---|---|
| 16 | 16 | 16 | 16 |

**Result.** RAS_like ∩ FVPTC = 16/16 — the two subgroup AUC=1.000 rows are **not independent corroboration**, they are the same 16 slides labelled twice.

## 2. Raw OOF prediction distribution within the 16 RAS-like slides

| rank | submitter_id | label | prob_DM1 | fold |
|---:|---|:---:|---:|---:|
| 1 | TCGA-EM-A3FL | DM2 | 0.000 | 1 |
| 2 | TCGA-EM-A1YB | DM2 | 0.001 | 4 |
| 3 | TCGA-EM-A1YD | DM2 | 0.019 | 1 |
| 4 | TCGA-DE-A2OL | DM2 | 0.028 | 3 |
| 5 | TCGA-EM-A1YE | DM2 | 0.035 | 4 |
| 6 | TCGA-EM-A4FH | DM2 | 0.043 | 2 |
| 7 | TCGA-EM-A2OY | DM2 | 0.097 | 3 |
| 8 | TCGA-EM-A1YC | DM2 | 0.103 | 2 |
| 9 | TCGA-FY-A3WA | DM2 | 0.108 | 3 |
| 10 | TCGA-EM-A1CW | DM2 | 0.151 | 2 |
| 11 | TCGA-EM-A3O9 | DM2 | 0.191 | 4 |
| 12 | TCGA-BJ-A0ZG | DM2 | 0.374 | 3 |
| 13 | TCGA-EM-A3FP | DM2 | 0.532 | 5 |
| 14 | TCGA-DJ-A13W | **DM1** | 0.549 | 3 |
| 15 | TCGA-DJ-A2PX | **DM1** | 0.807 | 5 |
| 16 | TCGA-FK-A3S3 | **DM1** | 0.834 | 5 |

**Separation gap.** lowest DM1 prob = **0.549**, highest DM2 prob = **0.532**. The 'perfect ranking' is achieved by **0.017** margin. If the borderline DM2 (TCGA-EM-A3FP, 0.532) had been ranked above the borderline DM1 (TCGA-DJ-A13W, 0.549), AUC would drop to 38/39 = 0.974.

## 3. Permutation null + bootstrap CI

| group | n (pos/neg) | observed AUC | perm p (2-sided) | null p99 | boot 95% CI |
|---|---:|---:|---:|---:|---|
| molecular_subtype=RAS_like | 16 (3/13) | 1.000 | 0.0030 | 0.923 | [1.000, 1.000] |
| histology_subtype=FVPTC | 16 (3/13) | 1.000 | 0.0039 | 0.923 | [1.000, 1.000] |
| molecular_subtype=BRAF_like | 41 (26/15) | 0.592 | 0.3405 | 0.715 | [0.400, 0.758] |
| histology_subtype=cPTC | 41 (26/15) | 0.592 | 0.3437 | 0.715 | [0.416, 0.767] |

**Combinatorial baseline.** With n_pos=3 / n_neg=13, P(perfect ranking | random) = 1/C(16,3) = **0.0018** (0.18%). Observed permutation p=0.0030 — model-produced ranking is **better than random label assignment** but only barely (perfect by chance every 1 in 333 random label shuffles).

**Bootstrap CI [1.000, 1.000] is misleading.** It reflects deterministic resampling of an already-perfect ranking — not external generalization uncertainty. With n_pos=3, the *effective* sample size is tiny.

## 4. Histology / molecular shortcut probe

| Probe | AUC / corr |
|---|---:|
| Histology-only LogReg → label (DM1 vs DM2) | 0.680 |
| Molecular-only LogReg → label (DM1 vs DM2) | 0.680 |
| corr(FVPTC indicator, CLAM prob_DM1) | -0.359 |

**Reading.** Knowing only the histology label (FVPTC/cPTC) lets a logistic regression hit AUC=0.68 on the DM1/DM2 task — purely from base-rate asymmetry (81% of FVPTC → DM2, 63% of cPTC → DM1). CLAM's overall AUC=0.746 is only +0.07 over this trivial baseline. corr(FVPTC, prob_DM1)=−0.36 indicates CLAM has internalised histology pattern as a feature.

## 5. Per-fold composition (unstratified KFold seed=42)

| fold | n_total | n_pos (DM1) | n_neg (DM2) |
|---:|---:|---:|---:|
| 1 | 2 | 0 | 2 |
| 2 | 3 | 0 | 3 |
| 3 | 5 | 1 | 4 |
| 4 | 3 | 0 | 3 |
| 5 | 3 | 2 | 1 |

**Issue.** Folds 1, 2, 4 contained ZERO RAS-like positives. All 3 RAS-like DM1 cases concentrate in folds 3 (1) and 5 (2). Models in folds 1/2/4 saw all 3 RAS-like DM1 in training and thus produced confident low-prob_DM1 for their test-set RAS-like DM2 cases — locking in 11 of the 13 DM2 ranks. Strict OOF accounting still holds (each test slide was held out from its own fold), but the apparent *separation* is easier to achieve than under a stratified split.

## 6. Split-stress retrain (Strategies A / B / C)

(retraining still running — re-run `audit_make_report.py` after `RETRAIN_SUMMARY.json` exists)

## 7. What this audit shows

- **Not a hard train/test leak.** OOF predictions come from 5 separately-trained models that each held out their fold's slides. Patient/case-level deduplication confirms 59 unique cases = 59 unique slides (no patient leakage).
- **The two AUC=1.000 rows are one finding, not two.** RAS_like ∩ FVPTC = 16/16; same slides labelled twice.
- **The perfect ranking is fragile.** 0.017 prob margin between lowest DM1 (0.549) and highest DM2 (0.532). n_pos=3 → effective sample size is tiny.
- **Histology shortcut is plausible.** Histology-only LogReg AUC=0.68 already solves much of the task; corr(FVPTC, CLAM prob)=−0.36.
- **Unstratified KFold(seed=42) concentrated all 3 positives into 2 of 5 folds**, leaving the other 3 folds to confidently downrank their RAS-like DM2 test slides.

**For Paper 2 reporting:** treat overall AUC=0.746 + Korean K2 prospective validation as the load-bearing result. Move 'RAS-like AUC=1.000' from main figure to an honest caveat box — report alongside the 0.017 separation gap and the n_pos=3 limit so reviewers can't weaponize the apparent perfection.
