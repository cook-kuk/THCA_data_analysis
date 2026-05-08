# Cover Letter — Paper 2 (H&E → DM1 image classifier)

**Target journal:** Cell Reports Medicine
**Date:** 2026-05-08
**Corresponding author:** Seungho Cook

---

## Para 1 — VOICE-PROTECTED (do NOT generate prose)

```
[Para 1 — author keyboard only. State the paper, the discovery, why it matters
to the journal, in author's voice. NOT to be auto-generated.]
```

---

## Para 2 — Significance vs concurrent work

The H&E-to-molecular space has advanced rapidly in the past year, most notably
with GigaTIME [@valanarasu_2026_gigatime], a foundation-model framework that
predicts cell-level multiplex immunofluorescence (mIF) marker expression from
H&E across pan-cancer cohorts. Our work occupies an orthogonal layer. GigaTIME
operates at the single-cell level and outputs virtual mIF channels for tumor-
immune microenvironment phenotyping; in contrast, the present study performs
slide-level classification of a molecularly defined transcriptomic
sub-population — DM1, an 8-gene RNA-defined sub-stratifier identified within
BRAF/RAS-negative papillary thyroid carcinoma (PTC) and enriched for a
Hashimoto-overlap immune phenotype. The disease scope is also distinct: rather
than pan-cancer TIME mapping, we focus on a single endocrine indication where
the molecular sub-class carries direct downstream implications for
radioactive-iodine refractoriness and immune-checkpoint candidacy. The two
approaches are therefore complementary rather than competitive.

---

## Para 3 — Method, result, and clinical translation

We pair a frozen tile-level foundation-model encoder (UNI when access is
gated, ImageNet-21k ViT-L as fallback) with a CLAM gated-attention multiple-
instance learning (MIL) aggregator to produce slide-level DM1 probabilities.
In Phase 1 spatial validation on GSE250521, attention maps correlated with
ground-truth molecular maps in 8 of 16 slides at |ρ|>0.3 (max ρ = -0.61),
clearing the pre-registered kill-switch threshold. In Phase 2 classification
on TCGA-THCA, the model reached pooled cross-fold AUC 0.76 with mean fold AUC
0.83 ± 0.21 (N=89 final). Robustness checks include bootstrap 95% CI
[0.611-0.862], stage-stratified subgroup analysis, calibration assessment with
Hosmer-Lemeshow goodness-of-fit, and tile-count sensitivity, all pre-specified.
A closure baseline using a ResNet50 encoder on the same MIL head yielded AUC
0.55, confirming that encoder architecture, not biology, was the prior
limitation. Clinically, this enables a low-cost chain:
a $0 H&E slide flags DM1-likely cases, triggers a reflex 8-gene RNA panel,
and surfaces patients at elevated RAI-refractory risk and ICI-candidate
immune phenotype.

---

## Para 4 — Editor ask

We believe Cell Reports Medicine is the appropriate venue: the manuscript
sits squarely at the intersection of clinical translation, AI/digital
pathology, and thyroid endocrinology that the journal regularly publishes. We
respectfully suggest reviewers spanning the three required perspectives.

- `[Reviewer 1 Name, Institution, thyroid imaging / digital pathology]`
- `[Reviewer 2 Name, Institution, ML pathology / weakly-supervised MIL]`
- `[Reviewer 3 Name, Institution, thyroid molecular endocrinology / RAI-refractory PTC]`

We have no conflicts of interest with the suggested reviewers. Thank you for
considering this work.
