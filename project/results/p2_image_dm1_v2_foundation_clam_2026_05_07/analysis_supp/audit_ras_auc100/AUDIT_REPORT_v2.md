# Paper 2 image-DM1 — RAS-like AUC=1.000 audit (v2 final)

## Verdict

**RAS-like AUC=1.000 is fragile AND the broader image-DM1 claim has multiple paper-blocking confounds.**

| Issue | Severity | Evidence |
|---|:---:|---|
| RAS_like ∩ FVPTC = 16/16 (same slides, not independent) | **MAJOR** | C1 |
| RAS-like 1.000 separation gap = 0.017 prob | **MAJOR** | C2 |
| Multi-seed RAS-like AUC range [0.44, 0.92] | **MAJOR** | A multi-seed |
| CPU rerun seed=42 gives 0.79 not 1.00 (irreproducible) | **MAJOR** | A seed=42 |
| Male AUC = 0.35 (worse than random) | **PAPER-KILLING** | E1 |
| All 3 RAS-like DM1 cases from 2 TSS centers (DJ, FK); all 10 TSS=EM cases are DM2 | **PAPER-KILLING** | E2 |
| 0 off-diagonal (cPTC, RAS_like) cases — cannot disentangle histology from molecular subtype | **PAPER-KILLING** | E4 |
| Clinical-only LR (histology + sex + TSS) AUC=0.768 **> CLAM 0.746** | **PAPER-KILLING** | E5 |

## A. Data structure (no retraining)

- 59 slides = 59 unique cases (no patient-level leakage).
- Histology × molecular subtype crosstab (in our 59):

| | RAS_like | BRAF_like | unknown |
|---|---:|---:|---:|
| FVPTC | 16 | 0 | 0 |
| cPTC | 0 | 41 | 0 |
| unknown | 0 | 0 | 2 |

→ **Histology and molecular subtype are 100% co-linear in our subset.** RAS_like ∩ FVPTC = 16/16. The two subgroup AUC=1.000 rows are one finding labelled twice.

## B. Raw OOF predictions for the 16 RAS-like slides

| rank | submitter_id | TSS | label | prob_DM1 | fold |
|---:|---|:---:|:---:|---:|---:|
| 1 | TCGA-EM-A3FL | EM | DM2 | 0.000 | 1 |
| 2 | TCGA-EM-A1YB | EM | DM2 | 0.001 | 4 |
| 3 | TCGA-EM-A1YD | EM | DM2 | 0.019 | 1 |
| 4 | TCGA-DE-A2OL | DE | DM2 | 0.028 | 3 |
| 5 | TCGA-EM-A1YE | EM | DM2 | 0.035 | 4 |
| 6 | TCGA-EM-A4FH | EM | DM2 | 0.043 | 2 |
| 7 | TCGA-EM-A2OY | EM | DM2 | 0.097 | 3 |
| 8 | TCGA-EM-A1YC | EM | DM2 | 0.103 | 2 |
| 9 | TCGA-FY-A3WA | FY | DM2 | 0.108 | 3 |
| 10 | TCGA-EM-A1CW | EM | DM2 | 0.151 | 2 |
| 11 | TCGA-EM-A3O9 | EM | DM2 | 0.191 | 4 |
| 12 | TCGA-BJ-A0ZG | BJ | DM2 | 0.374 | 3 |
| 13 | TCGA-EM-A3FP | EM | DM2 | 0.532 | 5 |
| 14 | TCGA-DJ-A13W | DJ | **DM1** | 0.549 | 3 |
| 15 | TCGA-DJ-A2PX | DJ | **DM1** | 0.807 | 5 |
| 16 | TCGA-FK-A3S3 | FK | **DM1** | 0.834 | 5 |

**Separation gap = 0.017** (lowest DM1 0.549 − highest DM2 0.532). 1 misrank → AUC = 0.974. Bootstrap CI [1.000, 1.000] is deterministic resampling, not generalization.

## C. Permutation null + bootstrap CI

| group | n (pos/neg) | observed AUC | perm p (2-sided) | boot 95% CI |
|---|---:|---:|---:|---|
| molecular_subtype=RAS_like | 16 (3/13) | 1.000 | 0.0030 | [1.000, 1.000] |
| histology_subtype=FVPTC | 16 (3/13) | 1.000 | 0.0039 | [1.000, 1.000] |
| molecular_subtype=BRAF_like | 41 (26/15) | 0.592 | 0.3405 | [0.400, 0.758] |
| histology_subtype=cPTC | 41 (26/15) | 0.592 | 0.3437 | [0.416, 0.767] |

Combinatorial baseline P(perfect | random)=1/C(16,3)=0.0018. Perm p=0.003 means model ranking is barely above the random-shuffle null (1 in 333).

## D. Split-stress retrain

| split strategy | overall AUC | RAS-like AUC | BRAF-like AUC |
|---|---:|---:|---:|
| original (GPU, KFold seed=42) | 0.746 | **1.000** | 0.592 |
| KFold seed=42 (CPU) | 0.722 | 0.795 | 0.646 |
| KFold seed=1 (CPU) | 0.664 | 0.590 | 0.600 |
| KFold seed=7 (CPU) | 0.795 | 0.692 | 0.736 |
| KFold seed=13 (CPU) | 0.759 | 0.821 | 0.736 |
| KFold seed=99 (CPU) | 0.747 | 0.923 | 0.572 |
| KFold seed=2026 (CPU) | 0.730 | 0.436 | 0.741 |
| **StratifiedKFold(3)** | 0.793 | **0.923** | 0.738 |
| **LOO 16 RAS-like** | n/a | **0.923** | n/a |

- Original 1.000 does not reproduce on CPU (seed=42 gives 0.795). Multi-seed range [0.44, 0.92].
- Stratified-CV / LOO converge to **0.92** — the honest estimate.

## E. Phase-2 confounders — paper-killing findings

### E1. Sex-stratified AUC

| sex | n (pos/neg) | AUC | 95% boot CI |
|---|---:|---:|---|
| sex=Female | 46 (24/22) | 0.833 | [0.705, 0.942] |
| sex=Male | 13 (5/8) | 0.350 | [0.000, 0.733] |

→ Model is **worse than random in n=13 Males**. CI lower bound 0.000. Sex-stratified review will flag this.

### E2. TSS distribution among 16 RAS-like slides

| TSS | n | n DM1 | mean prob_DM1 |
|---|---:|---:|---:|
| BJ | 1 | 0 | 0.374 |
| DE | 1 | 0 | 0.028 |
| DJ | 2 | 2 | 0.678 |
| EM | 10 | 0 | 0.117 |
| FK | 1 | 1 | 0.834 |
| FY | 1 | 0 | 0.108 |

→ **All 3 RAS-like DM1 cases concentrate in TSS=DJ (×2) + FK (×1).** All 10 RAS-like cases from TSS=EM are DM2. The model could be classifying **scanning center / staining batch**, not biology.

### E5. Clinical-covariate baseline beats CLAM

| Predictor | AUC |
|---|---:|
| Histology only (FVPTC indicator) | 0.680 |
| Histology + sex + TSS dummies (5-fold OOF LR) | **0.768** |
| CLAM-UNI overall (5-fold OOF) | 0.746 |

→ **CLAM image gain over clinical = -0.022** — image model **does not add over a 3-feature logistic regression** on histology + sex + TSS.

## F. Phase-3 decisive retraining

### E6. Init-only variance (StratifiedKFold split fixed)

| init seed | overall AUC | RAS-like AUC | BRAF-like AUC |
|---:|---:|---:|---:|
| 42 | 0.793 | 0.923 | 0.738 |
| 1 | 0.679 | 0.795 | 0.623 |
| 7 | 0.805 | 0.949 | 0.692 |
| 13 | 0.795 | 1.000 | 0.751 |
| 99 | 0.814 | 0.769 | 0.803 |
| 2026 | 0.762 | 0.615 | 0.703 |

RAS-like across 6 init seeds at fixed split: median 0.923, range [0.615, 1.000]

### E7. Label-shuffle null (3 reps)

| rep | overall AUC | RAS-like AUC |
|---:|---:|---:|
| 1 | 0.534 | 0.317 |
| 2 | 0.554 | 0.698 |
| 3 | 0.648 | 0.635 |

→ Shuffled-label AUC ≈ 0.5 = no memorization shortcut. Healthy null (model isn't fitting noise structure).

### E8. Leave-one-TSS-out — TSS batch-effect probe

| TSS held out | n test | n DM1 | LOTO AUC |
|---|---:|---:|---:|
| BJ | 13 | 3 | 0.933 |
| CE | 5 | 5 | n/a |
| DE | 3 | 1 | 1.000 |
| DJ | 7 | 6 | 0.667 |
| EL | 4 | 4 | n/a |
| EM | 13 | 2 | 0.955 |
| ET | 4 | 2 | 0.500 |
| FK | 3 | 2 | 1.000 |
| FY | 3 | 0 | n/a |
| H2 | 1 | 1 | n/a |
| J8 | 1 | 1 | n/a |
| KS | 2 | 2 | n/a |

- Pooled OOF overall AUC under LOTO = **0.602**
- Pooled OOF RAS-like AUC under LOTO = **0.308**

→ If pooled LOTO RAS-like AUC drops sharply vs the 0.92 from random / stratified split, the AUC was TSS-batch-driven. If it stays ≥ 0.85, the signal survives center-holdout.

## G. Phase-4 — solutions attempted

Tried 4 fixes for the TSS / sex / image-gain failures:

### S1. TSS-ComBat on UNI features (per-tile mean centering by TSS)

| init seed | overall | RAS-like | BRAF-like | Male | Female |
|---:|---:|---:|---:|---:|---:|
| 42 | 0.707 | 0.615 | 0.672 | 0.950 | 0.688 |
| 1 | 0.626 | 0.923 | 0.487 | 0.750 | 0.589 |
| 7 | 0.661 | 0.692 | 0.682 | 0.625 | 0.672 |
| 13 | 0.629 | 0.821 | 0.536 | 0.525 | 0.661 |

- RAS-like median 0.821 (was 0.92 stratified, 1.000 original) — TSS effect partially removed → AUC drops to honest range.
- Male median 0.750 (was 0.35 raw) — sex asymmetry partially fixed.
- Overall median ~0.66 (was 0.74) — image gain smaller after batch correction.

### S2. TSS-balanced StratifiedKFold on (label × TSS-bucket)

| init seed | overall | RAS-like | BRAF-like |
|---:|---:|---:|---:|
| 42 | 0.731 | 0.795 | 0.651 |
| 1 | 0.800 | 0.923 | 0.721 |
| 7 | 0.844 | 0.897 | 0.813 |
| 13 | 0.801 | 0.923 | 0.754 |

- RAS-like median 0.923 — explicit TSS stratification preserves the apparent RAS-like advantage. Combined with the LOTO 0.308 result, this confirms: RAS-like AUC stays high **only when training and test contain the same TSS centers**.

### S3. Multimodal LR (CLAM_prob + clinical features)

| Predictor | OOF AUC |
|---|---:|
| clinical_only | 0.768 |
| clam_only | 0.691 |
| multimodal | 0.756 |

→ **Image channel net gain over clinical = -0.012**. Adding CLAM to clinical features does not improve AUC; the image model is **redundant with clinical metadata**.

### S5. EM-out subsample (drop 10 RAS-like-EM-DM2 anchor)

| init seed | overall | RAS-like | BRAF-like |
|---:|---:|---:|---:|
| 42 | 0.667 | 0.667 | 0.633 |
| 1 | 0.688 | 0.667 | 0.682 |
| 7 | 0.712 | 0.778 | 0.690 |
| 13 | 0.671 | 0.778 | 0.633 |

- After dropping the 10 (TSS=EM, RAS_like, all DM2) slides: RAS-like median = 0.778 (was 0.92). Signal degrades meaningfully without the EM-DM2 anchor — confirms that the apparent RAS-like advantage was driven by the EM-DM2 majority.

## H. Bottom line for Paper 2

**Cannot ship 'image-DM1 AUC=0.83 (UNI)' as the headline as currently framed.** Three independent confounds:

1. **Sex asymmetry** — Female 0.83 vs Male 0.35. Reviewers will ask for sex-balanced training.
2. **TSS batch confound** — molecular subtype is co-linear with scanning center. Need leave-one-TSS-out as primary reporting metric, OR ComBat-equivalent feature normalization.
3. **No image gain over clinical** — histology + sex + TSS dummies already give AUC=0.768. CLAM at 0.746 is *below* this baseline; the image channel is not the load-bearing signal in TCGA.

**Path forward options:**

- **Defer Paper 2 launch** until Korean K2 H&E (or FFPE multi-site) cohort is in. K2 will have different TSS structure → if RAS-like AUC=0.92 holds on K2 with TSS-balanced split, the signal is real.
- **Reframe to clinical-augmented model**: report (histology + sex + TSS + CLAM-prob) as the predictor; compare against (clinical only) baseline; main claim becomes "CLAM adds X% on top of clinical". Currently X = 0%.
- **Restrict to TSS-balanced subsample** for the main result; report TSS-stratified analysis as primary.

**Do NOT publish RAS-like AUC=1.000 as currently framed.**

---

## I. Final synthesis — what each fix bought us

| Fix | What it solves | What it doesn't |
|---|---|---|
| S1 ComBat | TSS center batch effect; sex asymmetry partially | n_pos=3 RAS-like instability remains |
| S2 TSS-balanced split | confirms TSS-balanced AUC ≈ 0.91 | RAS-like AUC still seed-dependent [0.80, 0.92]; **LOTO 0.308 still the truth** |
| S3 multimodal LR | quantifies image gain | image gain = NEGATIVE; clinical alone (0.768) > clinical + image (0.756) |
| S5 EM-out | confirms EM-DM2 anchor is load-bearing | RAS-like AUC drops to 0.72; signal partially survives |

**No fix recovers the 1.000.** The honest TCGA-only ceiling is overall AUC ≈ 0.65–0.80, RAS-like AUC ≈ 0.70–0.90 (highly seed-dependent), and the image channel adds ZERO over (histology + sex + TSS) logistic regression.

**The fundamental fix is more data, not more analysis.**

- TCGA-THCA full WSI = ~500 slides (10× current 59)
- K2 cohort (planned) = ~320 slides with different TSS structure → independent test of any signal that survives TCGA training
- FFPE multi-site (planned) = 600 slides → adequate power for stratified analysis

Until then: **defer Paper 2 main-text image-DM1 claim**, or reframe as "clinical-only + image as auxiliary" with image gain reported truthfully (0% on TCGA).
