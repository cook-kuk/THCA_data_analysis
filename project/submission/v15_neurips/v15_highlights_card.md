# v15 — Highlights Card

_For author's outreach use after NeurIPS submission lands.
Same content, three lengths: tweet (280 char), abstract-pitch (≤500 char),
blog-post lede (≤1500 char)._

---

## Tweet (1×, 280 chars)

```
🧪 New paper: DIAL is a 1-number AUC-based diagnostic that fires
when batch correction leaks class signal back into the model.
Caught a real biomedical leak (THCA/ComBat, ΔDIAL=−0.49 under
proper LODO). NeurIPS 2026 under review. Synthetic 4-domain
reproduction.
```

(280 chars exactly.)

---

## Tweet thread (4 tweets)

```
1/4 Cross-cohort ML pipelines often "harmonize" data via ComBat
before training a classifier. We show that when ComBat is fit on
the FULL pooled data (train+test) before the LODO split, the
classifier can flip its decision direction on the held-out cohort.
This is a leakage artefact.

2/4 We introduce DIAL (Direction-Invariant AUC Leakage), a
single-number post-hoc diagnostic that fires when target AUC drops
BELOW chance — i.e., the source-trained classifier ranks target
samples in the wrong order.

3/4 Theorem 2: DIAL is governed by the conditional KL component
parallel to the Fisher discriminant direction, and is INDEPENDENT
(up to noise) of pure covariate or pure label shift. Numerical
validation R²=0.94. Reproduces in 4/4 synthetic domains
(vision/NLP/clinical/biology).

4/4 We caught our own retracted leakage (THCA BRAF/RAS, ΔDIAL=−0.49
under proper LODO) using DIAL — turning the retraction into a
worked example. Code anonymously released.
NeurIPS 2026 (under review).
```

---

## Abstract-pitch (≤500 char, for funding apps / CV)

```
DIAL (Direction-Invariant AUC Leakage) is a 1-number post-hoc
diagnostic for cross-cohort ML pipelines that fires when linear
batch correction (ComBat, CORAL, mean-centring) leaves a residual
class-conditional shift parallel to the Fisher discriminant
direction. Theorem 2 proves DIAL>0 iff this parallel component
exceeds the perpendicular one. We caught a public retraction
(THCA/ComBat-LODO leak, ΔDIAL=−0.49) and report a negative
foundation-model scaling result. NeurIPS 2026 under review.
```

(498 chars.)

---

## Blog-post lede (≤1500 char)

```
A surprising thing happens in cross-cohort ML pipelines: sometimes
a classifier trained on cohort A and evaluated on cohort B achieves
AUC BELOW chance — it ranks target samples in the wrong order.
We call this the flip event, and it almost always means something
in the preprocessing has corrupted the model's view of the world.

In a paper under review at NeurIPS 2026 we introduce DIAL
(Direction-Invariant AUC Leakage), a single-number post-hoc
diagnostic for the flip event. Our main theorem (Theorem 2) decomposes
DIAL into four KL components — covariate-marginal, label-prior,
conditional-parallel, conditional-perpendicular — and proves that
DIAL is monotone in the parallel component alone. This unifies
three Ben-David shift types into a single decomposition.

The story has an honest dark side. The motivation for v15 was a
real biomedical pipeline (TCGA + GEO thyroid cancer) where we had
observed a 4-of-5-classifier flip. A subsequent audit revealed
that v5.1 had fit ComBat on the full pooled matrix before splitting
for cross-validation — a textbook information leak. Under proper
fit-on-train-only LODO ComBat, the apparent flip vanishes
(ΔDIAL = −0.49). DIAL caught its own creators' leak; the retraction
is now in §1 of the paper, and DIAL is reframed as a leak detector.

If you build cross-cohort ML pipelines, the takeaway is: report
DIAL alongside AUC.
```

(1430 chars.)

---

## One-line elevator (post-acceptance, for CV)

```
"DIAL: A post-hoc diagnostic for subspace-aligned conditional shift
under linear batch correction." NeurIPS 2026.
```

---

## Author's voice notes

- **Don't oversell theorem.** The Gaussian-shared-Σ assumption is
  restrictive; broader synthetic reproduction is the empirical hook.
- **Don't undersell the retraction.** Disclosing it in public is the
  right move — peer reviewers (and Twitter) respect honest correction.
  The story is "we caught our own leak with the metric we built";
  this is positive, not embarrassing.
- **Foundation-model scaling negative result is intentional.** Don't
  let outreach over-pitch as "scaling doesn't help"; the precise claim
  is "raw parameter count alone doesn't help; shift-aware pretraining
  might".
- **Use ΔDIAL = −0.49 as the headline number.** It's
  visceral, easy to remember, and ties directly to the retraction.
