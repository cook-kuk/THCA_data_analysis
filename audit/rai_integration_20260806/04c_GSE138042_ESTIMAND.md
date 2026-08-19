# AUDIT 04c — GSE138042 estimand mismatch and an annotation-batch signal

Generated 2026-08-06.

## The problem

The first pass reported two numbers side by side as if they described one analysis:

- primary contrast, refractory 13 versus radiosensitive 10: Cohen's d = −0.42, P = 0.34
- "composition-adjusted" logistic: OR = 0.288 (0.101–0.819), **P = 0.020**

The adjusted model was fitted on **all 95 libraries** with the outcome defined as "is this
library named `RAIR-*`". That is a different question from the primary contrast, and the
significant adjusted odds ratio belongs to the confounded comparison, not the clean one.

| Model | Population | Estimand |
|---|---|---|
| A (reported as "adjusted") | 13 RAIR vs 82 others — including 17 benign follicular adenomas, 3 medullary carcinomas, and 62 carcinomas of which only 10 have a recorded radioiodine outcome | "is a library named RAIR distinguishable from the rest of this series" |
| B (the paper's actual contrast) | 13 refractory vs 10 radiosensitive, all carcinomas with a documented radioiodine outcome | "among treated patients, does the score differ by radioiodine outcome" |

## The same covariates applied to the correct contrast

| Model | OR (panel z) | 95% CI | P |
|---|---|---|---|
| y ~ panel z | 0.442 | 0.091–2.140 | 0.310 |
| y ~ panel z + immune | 0.204 | 0.024–1.699 | 0.142 |
| y ~ panel z + immune + stroma | 0.205 | 0.023–1.838 | 0.157 |

**None is significant.** The point estimate moves in the hypothesised direction and is
similar in magnitude to Model A, but with 23 patients the interval is uninformative.

**Correction required:** the sentence "adjusted for immune and stromal content, the panel
retains OR 0.288, P = 0.020" must not be attached to the refractory-versus-radiosensitive
comparison. Either report it explicitly as the whole-series contrast with its comparator
contamination stated, or drop it.

## A second problem: the refractory group is annotated differently

Cross-tabulating histology label against radioiodine status:

| `cancer_type` string | no RAI outcome recorded | radiosensitive | refractory |
|---|---|---|---|
| Follicular adenoma | 17 | 0 | 0 |
| Follicular cancer | 17 | 1 | 0 |
| Medullar cancer | 3 | 0 | 0 |
| Papillary cancer | 35 | 9 | 0 |
| **Papillary thyroid cancer** | 0 | 0 | **7** |
| **Follicular thyroid cancer** | 0 | 0 | **5** |
| **Poorly differentiated thyroid cancer** | 0 | 0 | **1** |

Every refractory sample carries a histology string that **no non-refractory sample uses**
("Papillary **thyroid** cancer" rather than "Papillary cancer"). The two groups were
evidently annotated in separate batches, and the `RAIR-*` versus `TC-*` library naming
points the same way.

This does not by itself invalidate the comparison, but it means group membership is
perfectly confounded with annotation batch, and therefore plausibly with collection period,
processing run and sequencing batch. Combined with the poorly differentiated case appearing
only in the refractory arm, the earlier reading stands and strengthens: **a substantial part
of the apparent effect reflects advanced-versus-early disease and batch, not radioiodine
behaviour.**

## Consequence

GSE138042 should be presented as **supportive but uninformative**: the direction is
consistent, no comparison in it reaches significance once the comparator is correct, the
cohort detects only d ≥ 1.24 at 80% power, and group membership is confounded with
annotation batch. It cannot validate or falsify the hypothesis.
