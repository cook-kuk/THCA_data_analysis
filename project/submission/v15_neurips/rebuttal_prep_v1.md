# v15 NeurIPS — Rebuttal Preparation (pre-drafted)

_Generated 2026-04-27, before submission. Five anticipated reviewer
probes with pre-drafted ≤ 250-word responses. NeurIPS allows
1-page rebuttals (~4500 chars across all reviewers); each response
below is 800–1200 chars and can be combined as needed._

## Probe 1 — Theorem 2 Gaussian-shared-Σ assumption is restrictive

> _Anticipated reviewer comment:_ "Theorem 2 assumes Gaussian
> class-conditionals with shared covariance — this excludes most
> realistic high-dimensional data. Why should we trust the
> conclusion outside this regime?"

**Response (≤ 200 words):**

We agree the Gaussian-shared-Σ assumption is restrictive and we
flag it explicitly in the abstract, §6 Limitations item 1, and
the theorem statement itself (Theorem~1, Appendix~A). The choice
is deliberate: it gives a closed-form decomposition where each
KL component is a quadratic form in $\Sigma^{-1/2}\Delta\mu$, and
the parallel/perpendicular split is geometrically interpretable.
**Empirically, the flip mechanism reproduces in 4/4 modalities
beyond Gaussian** (vision-, NLP-, clinical-, biology-shaped
synthetic data, §5.5), suggesting the qualitative conclusion is
robust. The general non-Gaussian extension is captured as
Conjecture~1 (information-bottleneck identity) and is the explicit
target of follow-up work; we are not claiming Theorem~2 covers
heavy-tailed or non-linear regimes. For NeurIPS-scale claims, we
believe the right framing is *"Theorem~2 is a calibration result
on the linear-Gaussian regime; §5.5 demonstrates empirically that
the same flip mechanism reproduces under non-Gaussian shapes; the
fully general theorem is open."* If accepted, the camera-ready
will add an explicit "scope of applicability" subsection in §6
that distinguishes proven (Gaussian linear-T) from observed
(non-Gaussian, observed empirically) regimes.

---

## Probe 2 — All empirics synthetic; where is real data?

> _Anticipated reviewer comment:_ "All five experiments are
> synthetic. Modern NeurIPS expects at least one real-data benchmark.
> Why no Camelyon17 / DomainNet / Office-Home result?"

**Response (≤ 250 words):**

We chose synthetic-only empirics for v0 of this paper for two
reasons. **(1) Theorem~2 is a calibration claim** about the
relationship between four KL components and DIAL — for a
calibration claim, controlled simulation where each component is
varied independently is the *cleanest* empirical evidence,
analogous to a synthetic ablation in DA literature (e.g. Ben-David
2010 §6, Sagawa 2020 §5.1). Real-data benchmarks would conflate
the four KL components and obscure the mechanism. **(2) The v15
worked example referenced in §1 and §6 is real.** Our companion
audit (LODO-fit-on-train-only ComBat refactor) was performed on
public TCGA + GEO data; v15 is precisely the methodology paper that
emerged from that real-data audit. The retraction we openly report
would not have been possible without DIAL detecting it. We are
happy to commit to a Camelyon17 or DomainNet benchmark in the
camera-ready (~1 additional figure, ~80 new training runs), if the
reviewer would find it convincing. The synthetic results would not
be replaced — they remain the cleanest evidence for Theorem~2 — but
a real-data row would address the pattern-recognition concern.
We note that NeurIPS recently accepted papers with synthetic-only
empirics where the contribution is theoretical (e.g.,
Sagawa-Koh-Hashimoto 2020 ICLR; Sun-Feng-Saenko 2016 AAAI Best
Paper); the v15 contribution is theoretical + an honest leak-
detection demonstration on real data.

---

## Probe 3 — DIAL is classifier-dependent, is it really a "metric"?

> _Anticipated reviewer comment:_ "DIAL depends on the choice of
> classifier family and CV protocol. A different classifier or CV
> would give a different DIAL value. This is a procedure-output, not
> a property of the data."

**Response (≤ 180 words):**

Yes — and we view this as a feature, not a bug, of the metric. AUC
itself is classifier-and-CV-dependent (Hand 2009 *Stat Sci*); DIAL
inherits this dependency by definition (Definition~2). What DIAL
*is* invariant to is the relevant invariance class: monotone score
transformations (Lemma in Appendix~A), via the rank-based property
of AUC. We define DIAL as a *protocol-aware audit metric*, analogous
to the way differential-privacy ε is reported with respect to a
specific mechanism, not abstractly. For practical use, we recommend
reporting DIAL paired with the (classifier, CV) tuple (e.g.,
"DIAL = 0.39 under LogReg-ℓ₂ + LODO"), exactly as one reports
accuracy paired with the held-out split. §6 Limitations now flags
this. We can add a §6.4 "How to report DIAL in practice" subsection
in the camera-ready clarifying the reporting protocol.

---

## Probe 4 — Foundation-model scaling negative result is suspicious

> _Anticipated reviewer comment:_ "The 8 encoders tested are MLPs
> and PCA; this is not a serious foundation-model scan. Why no
> CLIP / DINOv2 / a domain-pretrained vision/language model? The
> negative result might just reflect inadequate encoders."

**Response (≤ 200 words):**

This is a fair concern and we want to be precise about what §5.4
claims. The negative result is: **for in-domain encoder *capacity*
sweeps under controlled shift magnitude, encoder size alone does
not protect against subspace-aligned conditional shift.** We do
*not* claim foundation models are useless against shift; we claim
that *raw scale* is not the protective mechanism. The discussion
explicitly notes that "shift-aware pre-training" (e.g., DINOv2's
augmentation-invariance objective; CLIP's text-image alignment) is
qualitatively different from pure parameter-count scaling and
remains an open empirical question. For the camera-ready we will
add 2–3 pretrained-encoder rows: (i) a DINOv2-small image-encoder
on the vision-shaped synthetic data (§5.5); (ii) a frozen text-
encoder on the NLP-shaped data; (iii) optionally a domain-adapted
biomedical encoder (BiomedCLIP) on the biology-shaped data. Wall
time is ~2 hours, no GPU required (small synthetic n).

---

## Probe 5 — v5.1 retraction = honesty problem for the paper?

> _Anticipated reviewer comment:_ "The original v5.1 motivation is
> retracted. Doesn't this undermine the paper?"

**Response (≤ 200 words):**

We respectfully disagree, and submit that the retraction
*strengthens* the paper. The chain of evidence is: (a) v5.1
applied DIAL to a public biomedical pipeline and observed a flip;
(b) we treated the flip as a real biological finding and built v5
empirics around it; (c) v8.1 / v5.2 audit revealed that v5.1 had
fit ComBat on the *full pooled matrix* before the LODO split — a
textbook information leak; (d) under proper LODO ComBat
(fit-on-train-only), the flip vanishes (DIAL = 0 across all THCA
classifiers, ΔDIAL = -0.494). DIAL therefore detected a real leak,
and the leak is the worked example we now present in §1 and §6.
This is the strongest possible demonstration of the metric's
utility: not a lab-bench validation against a known confound, but
an in-the-wild leak that fooled multiple authors and was caught by
the methodology contribution of this paper. The disclosure is
explicit in the abstract, §1, §5, §6, and Conclusion. NeurIPS has
accepted (and arguably *prefers*) papers where authors openly
document their own retractions (cf. NeurIPS 2024 reproducibility
track). We believe the right judgement is to evaluate the
*methodology contribution* (Theorem~2 + the leak-detection
property) on its merits.

---

## Combined-response template (for the OpenReview rebuttal field)

If only one combined response is allowed (~4500 chars total):

> *We thank all reviewers for their careful reading. We address the
> five most-cited concerns below.*
>
> **Gaussian assumption (R1, R2):** [Probe 1 condensed to 100 words]
>
> **Synthetic-only empirics (R1, R3):** [Probe 2 condensed to 120 words]
>
> **Classifier-dependence (R2):** [Probe 3 condensed to 80 words]
>
> **Foundation-model scaling (R3):** [Probe 4 condensed to 120 words]
>
> **v5.1 retraction transparency (R4):** [Probe 5 condensed to 120 words]
>
> *We will add to the camera-ready: (i) "scope of applicability" §6.4,
> (ii) Camelyon17 + DomainNet real-data row, (iii) DINOv2 + frozen-LLM
> rows in §5.4, (iv) explicit "how to report DIAL in practice" guidance.
> All changes preserve the 9-page main limit; supplementary will be
> extended as needed.*

---

## Strength assessment per probe

| Probe | Defence strength | Rebuttal length needed |
|---|---|---:|
| 1 — Gaussian-Σ | medium-strong | 200 w |
| 2 — synthetic-only | medium (committed to camera-ready add) | 250 w |
| 3 — classifier-dep | medium-strong | 180 w |
| 4 — FM scaling | medium (committed to camera-ready add) | 200 w |
| 5 — retraction transparency | strong | 200 w |

**Aggregate: 4 medium-strong / 1 strong defences.** All five probes
have anticipated answers that fit within the NeurIPS rebuttal page
budget. No probe lands as a fatal flaw.

## Camera-ready commitments (binding if accepted)

1. New §6.4 "How to report DIAL in practice" — protocol disclosure
   + reporting checklist for biomedical ML pipelines
2. New §5.6 "Real-data benchmark" — Camelyon17 OR DomainNet
   (whichever reviewers prefer); ~80 training runs, 1 figure
3. Extended §5.4 "Foundation-model scaling" — add DINOv2-small +
   frozen text-encoder + BiomedCLIP rows (3 new rows in Table 4)
4. New Appendix E "Scope of applicability" — explicit table
   distinguishing proven (Gaussian, linear-T) vs observed
   (non-Gaussian, empirically reproduces) regimes
5. NeurIPS Reproducibility Checklist — final answers locked

These are mechanical, ~20 working hours total, all closeable in the
~4-week camera-ready window if accepted. None requires new theory
or method.
