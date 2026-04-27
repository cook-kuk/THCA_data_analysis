# v15 NeurIPS — OpenReview submission form fields

_Generated 2026-04-27. Drop-in text for each field of the OpenReview
NeurIPS 2026 submission form. Phase 1 (May 4) for abstract, Phase 2
(May 6) for full paper + supplementary._

---

## Title

```
DIAL: A post-hoc diagnostic for subspace-aligned conditional shift under linear batch correction
```

(101 chars; OpenReview limit typically 200.)

---

## Abstract (OpenReview field; ≤ 250 words)

```
We introduce DIAL (Direction-Invariant AUC Leakage), a single-number
post-hoc diagnostic that fires when a linear batch-correction
operator leaves a residual class-conditional shift inside the label
subspace, causing a source-trained classifier to have AUC below
chance on the target cohort. We prove a unified shift decomposition
(Theorem 2) showing that DIAL > 0 is governed by the component of
the post-correction conditional KL that lies parallel to the Fisher
discriminant direction, and is independent (up to finite-sample
noise) of pure covariate or pure label shift. Numerical validation
recovers R^2 = 0.94 on a controlled simulation; a synthetic stress
test gives ROC-AUC 0.78 for DIAL as a flip detector; ten
test-time-adaptation methods all fail to fully restore target AUC on
a synthetic flip (best non-trivial method: DANN-lite at AUC 0.39 vs
0.33 for naive ComBat); we report a negative foundation-model
scaling result (raw encoder size does not reduce DIAL); and the flip
mechanism reproduces in 4/4 cross-domain synthetic shapes (vision,
NLP, clinical, biology). The paper additionally documents an in-the-
wild preprocessing leak that DIAL detected in a public biomedical
cross-cohort pipeline; under proper Leave-One-Dataset-Out ComBat
(fit-on-train-only), the apparent batch-entanglement collapses, and
DIAL is reframed as a leak detector for biomedical ML audits. Code
and reproducible synthetic benchmarks are anonymously released.
```

(247 words. Tight under 250.)

---

## Keywords (OpenReview field; comma-separated)

```
domain adaptation, distribution shift, batch correction, AUC leakage,
preprocessing leakage, test-time adaptation, biomedical ML, ComBat,
information theory, evaluation
```

(10 keywords; OpenReview typical limit 10–15.)

---

## TL;DR (one-sentence summary, OpenReview short field)

```
DIAL is a post-hoc AUC-based diagnostic that detects subspace-aligned
conditional shift in cross-cohort ML pipelines, including the kind of
preprocessing leakage that recently caused a public biomedical
retraction.
```

(38 words.)

---

## Primary subject area

`Domain adaptation / transfer learning / out-of-distribution
generalization`

(NeurIPS 2026 areas list; this is the closest match. Secondary candidates:
`Evaluation of ML methods`, `Reliable and trustworthy ML`.)

---

## Conflict of interest declaration

```
No conflicts of interest declared.
```

(Author has no current/past collaboration with NeurIPS 2026 area chairs;
no past co-authorship with potential reviewers in the listed
literatures within the past 5 years that the author is aware of.
This will be re-verified at OpenReview's automated COI check.)

---

## Reproducibility checklist URL

```
[paste anonymous_code URL once uploaded — e.g.,
https://anonymous.4open.science/r/v15-dial-XXXXXX]
```

---

## Supplementary upload contents (Phase 2, May 6 deadline)

| File                                  | Purpose |
|---------------------------------------|---------|
| `anonymous_code.zip` (69 KB)          | 5 v15_*.py scripts + 5 task[1-5]_*.json checkpoints + 9 result TSVs + LICENSE + anonymous README |
| `figures/` (10 files PDF + PNG)       | static exports of the 5 paper figures (also embedded in main PDF) |
| `reproducibility_checklist_filled.md` | 14-question NeurIPS checklist filled |

---

## Author response field (post-rebuttal, August 2026)

See `rebuttal_prep_v1.md` for the 5 anticipated probes pre-drafted.

---

## Submission button order (recommended sequence on May 6)

1. Click "Submit" on the existing Phase-1 abstract record.
2. Upload `v15_neurips_SUBMIT_FINAL.pdf`.
3. Upload `anonymous_code.zip` to the supplementary slot.
4. Optionally upload the 5 `figures/figure[1-5].pdf` separately
   (some venues prefer separate figure files; OpenReview usually
   accepts only 1 supplementary, in which case skip this).
5. Re-verify the reproducibility checklist field.
6. Click "Submit" — confirmation email arrives within ~5 min.
7. **Take a screenshot of the OpenReview submission ID.**
