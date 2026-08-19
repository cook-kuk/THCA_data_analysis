# AUDIT 04d — BLOCKING: manuscript Figure 6D is a tissue-type contrast, not a radioiodine contrast

**Priority: P0 fatal. This must be resolved before either integration strategy proceeds.**
Generated 2026-08-06.

## What the manuscript currently claims

`NATURE_COMMUNICATIONS_FULL_DRAFT_v2_2026_07_08.md:124`

> The axis aligns with external radioiodine-refractory biology: in GSE151179, post-RAI
> refractory tumours shifted toward a DM1-like thyroid-differentiation state
> (Cohen's d ≈ −1.0; Mann-Whitney p ≈ 10⁻⁴), indicating that the DM1 programme resembles the
> transcriptional state observed after clinical RAI failure (Fig. 6d).

`05_figure_captions_NC.md:125`

> **(D) DM1-like state in post-radioiodine refractory disease (GSE151179).** Boxplot of the
> thyroid-differentiation score in pre- versus post-RAI samples … post-RAI thyroid_diff
> Cohen's d = −1.01, Mann-Whitney p = 0.0001. The post-RAI dedifferentiation profile
> recapitulates the DM1 transcriptional state, providing an external clinical anchor that the
> DM1 axis names the same biology that defines RAI-refractory progression.

## The contrast is perfectly confounded with tissue type

Cross-tabulating the `collection before/after rai` field against `tissue type` in GSE151179:

| Tissue type | before RAI | after RAI |
|---|---|---|
| non-neoplastic thyroid | **13** | 0 |
| primary tumour | **17** | 0 |
| synchronous lymph node metastasis | 5 | 0 |
| lymph node metastasis post RAI | 0 | **9** |
| lymph node metastasis_1 post RAI | 0 | **4** |
| lymph node metastasis_2 post RAI | 0 | **4** |

The "after" arm is **100% lymph node metastasis**. The "before" arm is 37% normal thyroid and
49% primary tumour. Timing and tissue type are not merely correlated — they are perfectly
nested. No sample is both.

## Decomposition

Using the eight-gene panel score on the same series:

| Contrast | n (after vs before) | Cohen's d | P |
|---|---|---|---|
| All samples — as published | 17 vs 35 | **−0.859** | 3.8 × 10⁻⁴ |
| Tumours only (drop 13 non-neoplastic thyroid) | 17 vs 22 | **−0.377** | 0.043 |
| **Lymph node metastases only (tissue held constant)** | **17 vs 5** | **−0.132** | **0.49** |

And the contrast that actually dominates this dataset:

| Contrast | n | Cohen's d | P |
|---|---|---|---|
| **non-neoplastic thyroid vs lymph node metastasis** | 13 vs 22 | **+1.923** | 1.2 × 10⁻⁵ |
| lymph node metastasis vs primary tumour | 22 vs 17 | −0.411 | 0.098 |

Removing the normal thyroids halves the effect. Holding tissue type constant removes it
entirely. The published d ≈ −1.0 is largely **normal thyroid versus lymph node metastasis** —
a tissue-type difference that would appear in any thyroid differentiation score regardless of
whether radioiodine had ever been given.

(The published figure used a `thyroid_diff` module rather than the eight-gene panel, so the
magnitudes are not expected to match to the decimal. The direction, the significance, and the
decomposition pattern are what matter, and all three reproduce.)

## Why this is fatal rather than merely imprecise

The sentence claims the DM1 axis "names the same biology that defines RAI-refractory
progression". What the data show is that the DM1 axis distinguishes normal thyroid from
lymph node metastasis. Those are different statements, and the second one supports no
radioiodine claim at all.

A reviewer can reproduce this in ten minutes: GSE151179 is public, the tissue-type field is in
the sample characteristics, and the confound is visible in a single cross-tabulation.

## Required action before any manuscript edit

1. **Do not use the whole-series contrast.** If the panel is to be shown in GSE151179 at all,
   the comparison must hold tissue type constant, and that comparison is null
   (d = −0.13, P = 0.49, 17 vs 5).
2. **Rewrite `NC_v2:124` and `05_figure_captions_NC.md:125`.** Neither "aligns with external
   radioiodine-refractory biology" nor "names the same biology that defines RAI-refractory
   progression" survives.
3. **Reconsider whether Figure 6D should exist.** With tissue held constant the panel is only
   powered to detect d ≥ 1.4 in this series (17 versus 5). A null on five samples is not worth
   a main-figure panel.
4. **Check the generating script** `fig6d_post_rai_box.py` (listed in
   `05_figure_captions_NC.md:180`) to confirm which samples it includes, and whether the
   13 non-neoplastic thyroid samples are in the "pre-RAI" box as the metadata implies.

## Relationship to the other GSE151179 result

This is a **different** contrast from the one analysed earlier today. Today's primary analysis
used `rai uptake at the metastatic site` (Yes 25 / No 27), the field that actually records
whether radioiodine reached disease, and found d = +0.37, P = 0.53, collapsing to β = +0.035
under purity adjustment (`04b_GSE151179_INDEPENDENCE.md`).

So GSE151179 now has two independent null results against the panel — one on uptake, one on
timing with tissue held constant — and one strongly positive result that turns out to be
normal thyroid versus metastasis. The dataset does not support the manuscript's use of it.
