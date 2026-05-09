# RAS-like AUC=1.000 audit — data-split / shortcut investigation

**Question:** Paper 2 image-DM1 CLAM (UNI features) reports RAS-like / FVPTC subgroup AUC=1.000 (n=16, n_pos=3 / n_neg=13). Is this real biological signal, a data-split artifact, or a histology shortcut?

## Verdict TL;DR

**The signal survives, but the 1.000 was an unstable upper outlier — not a reproducible result.**

| split strategy | RAS-like AUC | overall AUC |
|---|---:|---:|
| Original (GPU, seed=42, KFold 5) | **1.000** | 0.746 |
| CPU rerun seed=42, KFold 5 | 0.795 | 0.722 |
| Multi-seed median (n=6 seeds) | **0.744** [0.436–0.923] | 0.739 |
| StratifiedKFold(3) on label×subtype | **0.923** | 0.793 |
| LOO within 16 RAS-like slides | **0.923** | n/a |

**Interpretation.**
1. **The exact 1.000 is not reproducible.** Even at seed=42, CPU rerun gives 0.795 (GPU↔CPU numerical drift). Original 1.000 was one realization of a high-variance estimator.
2. **Multi-seed range [0.44, 0.92]** — at seed=2026, RAS-like AUC drops to 0.44 (worse than random). With n_pos=3 in n=16, plain KFold(5) is unstable.
3. **Under more rigorous splits (stratified, LOO), RAS-like AUC stabilises at 0.923.** This is the honest headline number: RAS-like is genuinely easier to classify than BRAF-like (≈0.65), but not perfectly.
4. **Histology-only LogReg AUC = 0.68** — half the RAS-like advantage is base-rate (FVPTC 81% DM2). CLAM adds +0.24 over the trivial baseline.

**Recommendation for Paper 2:** report **RAS-like AUC = 0.92 (stratified-CV)** in main text; relegate the original 1.000 to a methods footnote with the seed-sensitivity disclosure. The clinical claim ("FVPTC/RAS-like sub-cohort is more separable than cPTC/BRAF-like") is preserved without overclaiming.

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

### Strategy A — multi-seed KFold(5)

| seed | overall AUC | RAS-like AUC | BRAF-like AUC | seconds |
|---:|---:|---:|---:|---:|
| 42 | 0.722 | 0.795 | 0.646 | 89.5 |
| 1 | 0.664 | 0.590 | 0.600 | 81.5 |
| 7 | 0.795 | 0.692 | 0.736 | 73.7 |
| 13 | 0.759 | 0.821 | 0.736 | 73.6 |
| 99 | 0.747 | 0.923 | 0.572 | 74.4 |
| 2026 | 0.730 | 0.436 | 0.741 | 73.1 |

### Strategy B — StratifiedKFold(3) on label × molecular_subtype

- Overall AUC: **0.793**
- RAS-like AUC: **0.923**
- BRAF-like AUC: **0.738**

Per-fold breakdown:

| fold | n_val | val AUC | RAS pos | RAS neg |
|---:|---:|---:|---:|---:|
| 1 | 20 | 0.889 | 1 | 5 |
| 2 | 20 | 0.860 | 1 | 4 |
| 3 | 19 | 0.900 | 1 | 4 |

### Strategy C — Leave-one-out on the 16 RAS-like slides

- n = 16 (3 DM1 / 13 DM2)
- LOO AUC = **0.923**

Each held-out slide was predicted by a model trained on the other 58 slides (other 15 RAS-like slides remained in training). Most stringent test of the AUC=1.000 claim.

## 7. Summary of evidence

**A. What is real (signal survives split changes):**
- Stratified-CV(3) and LOO-on-RAS-like both give **RAS-like AUC = 0.923** — robust to split choice.
- BRAF-like AUC stays in 0.57–0.74 across all strategies — consistent moderate signal.
- Overall AUC stays in 0.66–0.80 across 6 random seeds.
- Permutation test p=0.003 — the model's ranking within the 16 RAS-like slides is not from random noise.

**B. What is not real (the 1.000 itself is fragile):**
- CPU rerun at the same seed=42 gives RAS-like AUC=0.795, not 1.000 — the original number does not reproduce on CPU.
- Multi-seed range [0.436, 0.923] — at one seed (2026), RAS-like AUC is 0.436, *worse than random*. With only 3 positives in n=16, KFold(5) is too unstable to interpret any one realization.
- Bootstrap CI [1.000, 1.000] reflects deterministic resampling of an already-perfect ranking — not generalization uncertainty.
- The two rows (RAS_like AUC=1.000 + FVPTC AUC=1.000) are the same 16 slides labelled twice.
- Half the RAS-like advantage is histology base-rate (FVPTC-only LogReg AUC = 0.68; corr(FVPTC, CLAM prob_DM1) = −0.36).

**C. Recommended Paper 2 reporting change:**
- Replace headline `RAS-like AUC = 1.000 (boot CI [1.0, 1.0])` with `RAS-like AUC = 0.92 (StratifiedCV-3, LOO-confirmed; n=16, n_pos=3)`
- Footnote: "Original report 1.000 from a single KFold(seed=42) realization is not reproducible across GPU/CPU or alternative seeds; the stratified-CV / LOO estimate of 0.92 is more reliable."
- Reviewer-defense: "RAS-like advantage is not pure histology shortcut: histology-only baseline gives 0.68, CLAM adds +0.24."
