# v15 — 3-Reviewer Self-Review Simulation

_Generated 2026-04-27. Simulates three NeurIPS-style reviewers
(R1: theorist, R2: empiricist, R3: applied-ML practitioner) using
the NeurIPS 2026 review form. Goal is to anticipate scores and
weak points before submission, NOT to pass judgement._

NeurIPS 2026 expected rating axes:
- **Soundness** (1–4): is the methodology / theory rigorous?
- **Presentation** (1–4): is it clear and well written?
- **Contribution** (1–4): is it novel / impactful?
- **Overall** (1–10): NeurIPS recommended scale

---

## Reviewer 1 — Theorist

> *Reviewer profile: theory-track NeurIPS regular; areas covariate
> shift, learning theory, information geometry. Expects rigorous
> proof and clear assumptions.*

**Soundness: 3/4.**
Theorem 1's proof in Appendix A is well-structured (6 steps + 3
lemmas + 1 supporting Proposition). Each TODO has been replaced
with a citation or supporting lemma. The Gaussian-shared-Σ
assumption is explicitly stated in the theorem and abstract. The
deferred Conjecture 1 (information-bottleneck identity) is honestly
flagged. Two concerns: (a) the H-divergence sandwich (Proposition 1)
proof is a sketch ("by Pinsker's + Le Cam's"); a NeurIPS reviewer
may want one extra page of explicit calculation; (b) the McDiarmid
bound constant (Lemma 3) is cited from Agarwal 2005 but not
re-derived — fine for this venue, just possible nit.

**Presentation: 3/4.**
Theorem 1 is stated precisely; corollaries are clean; the proof
sketch is structured. The main paper compresses the proof to ~1
page (§4.2) which is appropriate. Notation is consistent: $w_Y$,
$\Sigma$, $P_Y$ used throughout. Minor: $\Delta_{\text{cond}}^{\parallel}$
notation in math is fine but spell-out "parallel component" in
prose to aid readers from non-information-geometry backgrounds.

**Contribution: 3/4.**
The unified-shift-decomposition theorem is genuinely novel: it
gives a closed-form connection between four KL components and
the AUC-flip event, in the linear-Gaussian regime. The connection
to Ben-David's $\mathcal{H}\Delta\mathcal{H}$ and Sugiyama's KLIEP
is correct and useful. Honest negative result on foundation-model
scaling is intellectually honest. Theorem 1's assumption regime
(Gaussian-shared-Σ + linear T) is restrictive but not stranger
than typical NeurIPS theory papers.

**Overall: 6/10.** Solid theoretical contribution. Above the
borderline. Would accept with revisions.

**Anticipated requests:**
- Probe the "monotone in $\Delta^{\parallel}$ alone" claim under
  shared-Σ relaxation (Σ_s ≠ Σ_t).
- Add an explicit non-Gaussian counter-example showing the theorem
  *fails* outside the assumption regime — this would strengthen
  scope of claim.

---

## Reviewer 2 — Empiricist

> *Reviewer profile: builds DA / TTA benchmarks, runs hundreds of
> experiments, expects real-data results. Less interested in proofs.*

**Soundness: 2/4.**
The synthetic-only empirical strategy is the main weakness from
this perspective. The 5 experiments are well-controlled and
internally consistent (claims verified against checkpoints), but
"v15 tested DIAL on 4 modality-shaped synthetic datasets" is a
weaker empirical claim than "v15 tested DIAL on 4 real benchmarks".
The retraction case study in §1 / §6 is a real-data demonstration
but is described qualitatively, not quantitatively in the
experiments section. The "negative result on foundation-model
scaling" rests on 8 encoders that are MLPs / random projections /
PCA — these are not the foundation models reviewers have in
mind.

**Presentation: 3/4.**
Tables 1, 3, 4 are well-organised. Figures 2–5 use consistent
error-bar notation. The TTA benchmark §5.3 is the cleanest
section: 10 methods, 5 seeds, target AUC + DIAL columns, easy to
compare. The cross-domain heatmap §5.5 is intuitive.

**Contribution: 3/4.**
DIAL itself is genuinely useful as an audit metric and the
metric-theory unification is a contribution. But the empirical
contribution is shallower than typical NeurIPS standards — no
new dataset, no new architecture, no SOTA claim. The honest
v5.1 retraction is a contribution to community norms but does
not substitute for empirical depth.

**Overall: 5/10.** Borderline. Would request major revisions:
add Camelyon17 + DomainNet + one foundation-model encoder for
the rebuttal phase.

**Anticipated requests:**
- Run DIAL on Camelyon17 (medical imaging, n ≈ 400, real cross-
  hospital shift). Wall time: ~2 hr GPU.
- Run DIAL on DomainNet (n ≈ 600k, 6 domains). Wall time: ~12 hr GPU.
- Add at least one bona-fide pretrained encoder row to §5.4.
  CLIP-ViT-B/32 or DINOv2-small, frozen.

---

## Reviewer 3 — Applied-ML practitioner

> *Reviewer profile: works in industry; cares about practical
> utility, code release, deployment. Less concerned with proofs.*

**Soundness: 3/4.**
The methodology is sound: DIAL is well-defined, has explicit
properties, and the limit cases (covariate-only / label-only / 
orthogonal-conditional all give DIAL = 0) match intuition. The
v5.2 retraction story is exactly the kind of audit-fix narrative
that practitioners care about: "we built a metric, applied it to
public data, caught our own leak." This is genuinely useful.

**Presentation: 4/4.**
Reads cleanly. Code release is anonymous-ready. Repro checklist
filled. Hyperparameter table in appendix is comprehensive.
Figures are static PDFs, not Plotly HTML — practical for printing.

**Contribution: 4/4.**
For my use case: I need a 1-number diagnostic I can drop into
existing pipelines to flag preprocessing leakage. DIAL is exactly
that. The Theorem 1 result tells me when to expect DIAL to fire.
The leak-detection demonstration on real biomedical data is
gold — I can show it to my CTO.

**Overall: 7/10.** Above the bar. Would accept.

**Anticipated requests:**
- Provide a *deployed example*: a regulatory-style audit log
  format (CSV or JSON) that the metric can output for a given
  pipeline. Makes it usable in FDA / EMA pipelines.
- Add a "How to interpret DIAL = 0.X" guide in §6.4.
  E.g., "DIAL ∈ [0, 0.05] = pipeline likely safe, DIAL ∈ [0.05, 0.1]
  = ambiguous, DIAL ≥ 0.1 = audit required."

---

## Aggregate prediction

| Reviewer | Soundness | Presentation | Contribution | Overall |
|---|---:|---:|---:|---:|
| R1 (theorist)     | 3 | 3 | 3 | 6 |
| R2 (empiricist)   | 2 | 3 | 3 | 5 |
| R3 (practitioner) | 3 | 4 | 4 | 7 |
| **Mean**          | **2.7** | **3.3** | **3.3** | **6.0** |

**Predicted NeurIPS verdict: borderline accept (P ≈ 30 %).**

NeurIPS 2025 main-track acceptance rate was ~25 %. A predicted
mean Overall of 6.0 with one reviewer at 5 puts the paper at the
borderline. The R2 (empiricist) probe — synthetic-only empirics —
is the single largest risk and is exactly what the rebuttal_prep
Probe 2 prepares for.

## Strongest defences

- R1's Soundness can move 3 → 4 if we add the explicit non-Gaussian
  counter-example (~15 hours of additional work, post-submission)
- R2's Soundness can move 2 → 3 if we commit Camelyon17 +
  DomainNet rows in the rebuttal (~24 hours of additional work,
  post-submission)
- R3's Overall is robust at 7+; little risk on this axis

## Sensitivity analysis

If the empiricist reviewer (R2) is replaced by a second theorist
or practitioner, mean Overall climbs to ~6.5 and acceptance
probability climbs to ~40%. If R2 is replaced by another
empiricist, mean Overall drops to ~5.5 and acceptance drops to
~20%. Bottom line: **the paper's fate depends on getting at
most 1 empiricist reviewer in the random-3 draw.**

## Decision

**Submit anyway.** The paper is honestly within the NeurIPS-
acceptance band, and the upside (signal value + camera-ready
opportunity to add real-data) outweighs the downside (re-target
to TMLR or workshop if rejected). NeurIPS rebuttal phase allows
for camera-ready commitments that can move R2's Soundness up.
