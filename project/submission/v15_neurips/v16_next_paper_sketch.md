# v16 — Next Paper Sketch (after v15 NeurIPS lands)

_Generated 2026-04-27. Pre-submission planning for the natural
follow-up paper to v15 DIAL._

## Lifecycle position

```
v15 (NeurIPS 2026, May 4-6 submit)
  ├── Gaussian shared-Σ + linear T regime
  ├── 5 synthetic experiments
  ├── 1 worked example of preprocessing leak
  └── Conjecture 1: information-bottleneck identity (deferred)

v16 (target: NeurIPS 2027 / ICML 2028)
  ├── EXTEND: non-Gaussian + non-linear regime  (Theorem 2)
  ├── BENCHMARK: real-data multi-modality       (5 datasets, 4 modalities)
  ├── BREAK:   when does DIAL fail?             (failure-mode analysis)
  └── PROVE:  Conjecture 1 from v15             (closing the loop)
```

## Three candidate framings, ordered by ambition

### Option A — "DIAL-2: Non-Gaussian extension"

**Pitch:** v15 was the linear-Gaussian skeleton. v16 generalises
Theorem 1 to (i) non-shared-covariance Gaussians, (ii) sub-Gaussian
class-conditionals, (iii) monotone-non-linear correction operators.

**Theoretical contribution:** prove that the parallel-component
monotonicity of DIAL holds whenever the class-conditional log-density
is "approximately linear" (formalisation: bounded-second-moment
assumption + locally-Lipschitz log-density).

**Empirical contribution:** test on the 4 modalities of v15.5
(now real data: Camelyon17, DomainNet, NLP-multilingual transfer,
biomedical cross-cohort).

**Why this is the strongest option:** it directly addresses the
biggest reviewer complaint of v15 (Gaussian assumption is too
restrictive), turns Conjecture 1 into Theorem 2, and lifts the
empirical credibility from synthetic to real.

**Effort:** 6 months (theorist + empiricist). Right size for a
PhD chapter / NeurIPS 2027.

### Option B — "DIAL-Apply: Real-data leak detection in the wild"

**Pitch:** apply DIAL to 50+ published biomedical ML papers and
report which ones likely have v5.1-style leak. Build a public
"leakage atlas" of published cross-cohort claims.

**Empirical contribution:** systematic audit of public ML literature
in genomics, radiology, clinical NLP. Quantify leak prevalence.

**Theoretical contribution:** none new (uses v15 Theorem 1 as-is).

**Why interesting:** would be a high-impact "Nature Methods"-style
paper rather than a NeurIPS paper. Visibility matters for academic
career; might be the highest-impact-per-effort follow-up.

**Effort:** 4-5 months (empiricist + biostatistician). Risky if
many flagged papers push back; needs careful rebuttal-management.

### Option C — "DIAL-Adapt: When DIAL flags, what do you do?"

**Pitch:** v15 is a diagnostic, not a corrective. v16 introduces
a *correction* method: when DIAL > τ, automatically suggest the
"smallest" linear correction that brings the parallel component
below threshold while preserving the perpendicular component.

**Theoretical contribution:** define the correction as a quadratic
optimization problem; derive analytical solution under Gaussian
assumptions.

**Empirical contribution:** show that DIAL-Adapt outperforms
naive ComBat / TENT / SHOT on the synthetic flips of v15 §5
*and* on Camelyon17.

**Why interesting:** moves DIAL from descriptive to prescriptive;
gives practitioners an actionable tool, which historically
attracts citations.

**Effort:** 5 months (theorist for the correction theorem +
empiricist for benchmarks).

## Recommendation

**Option A first** (Non-Gaussian extension). Reasons:

1. Most directly addresses the v15 reviewer-1 critique
2. Closes Conjecture 1 from v15 (intellectual completeness)
3. Real-data empirics naturally land alongside the theory
4. Establishes a "DIAL series" pattern (v15 = linear; v16 = general;
   v17 = adaptive) that's good for academic record

If v15 NeurIPS gets rejected in Sep 2026, Option A becomes the
**v15.5 revision** rather than a new paper — the v16 work is
re-purposed to strengthen the v15 submission for ICLR 2027 or
TMLR.

## v16 dependencies on v15

| v16 needs | v15 provides | Complete? |
|---|---|:--:|
| LinearComBat fit/transform decomposition | v5p2_combat_lodo.py | ✅ |
| Synthetic shift generators | v15_*.py scripts | ✅ |
| 4-modality test infrastructure | v15_cross_domain.py | ✅ (synthetic-only; needs real-data port) |
| Theorem 1 (Gaussian) | proof in Appendix A | ✅ |
| Conjecture 1 (info bottleneck) | stated, not proved | ⚠️ becomes v16 Theorem 2 |
| Reviewer feedback | submit May 6 → review Aug 2026 | ⏳ wait for August |

## Timeline

| Phase | Date | Activity |
|---|---|---|
| 0 | now | v15 submitted; wait for reviews |
| 1 | Sep 2026 | review NeurIPS feedback; decide v15 fate |
| 2 | Oct 2026 – Mar 2027 | v16 theory work (option A) |
| 3 | Apr – Jul 2027 | v16 empirics (Camelyon17, DomainNet, etc.) |
| 4 | May 2027 | v16 NeurIPS abstract registration |
| 5 | Sep 2027 | v16 reviewer feedback |

The 1-year cadence between v15 (NeurIPS 2026) and v16 (NeurIPS 2027)
is the standard PhD-level pace. Aggressive but feasible.

## Co-authoring opportunities

Open for v16:
- A theorist co-author for the non-Gaussian extension (Theorem 2)
- A biomedical-ML co-author for the real-data benchmarks (would
  strengthen the medical ML angle if v16 targets npj or
  Nature Comm rather than NeurIPS)
- Possibly the v5.2 audit team / 유 교수님 if the SNUBH cohort is
  available by mid-2027 (would unlock the real-data prospective
  validation that's been on the limitation list since v8.1)

## Risks

1. **Theorem 2 doesn't generalise cleanly.** Mitigation: have a
   fallback "weak-Theorem-2" stating monotonicity *up to a constant*
   that's still publishable.
2. **Real-data benchmarks don't reproduce the synthetic flip.**
   Mitigation: this would itself be a publishable finding ("DIAL
   fires on synthetic flips but is dominated by other artefacts on
   real data") — frame as a calibration paper.
3. **NeurIPS 2027 reviewers reject because v15 already exists.**
   Mitigation: position v16 as a *significant* extension (new theorem
   + real data + adaptive correction option), not an incremental
   note.

## Pre-emptive abstract draft for v16

```
DIAL (Direction-Invariant AUC Leakage), introduced in [v15], is a
post-hoc diagnostic for subspace-aligned conditional shift under
linear batch correction. Its theoretical guarantee was restricted
to Gaussian shared-covariance class-conditionals and linear
correction operators. We extend the guarantee to (i) non-shared-
covariance Gaussians, (ii) sub-Gaussian class-conditionals with
locally-Lipschitz log-density, and (iii) monotone-non-linear
correction operators. We close [v15]'s Conjecture 1 by proving the
information-bottleneck identity I(X̃;Y|B) = I(X;Y) - Δ‖_cond + o(1/n)
under the new conditions. Empirically, we apply the extended DIAL
to four real-data cross-cohort benchmarks (Camelyon17, DomainNet,
multilingual NLP transfer, biomedical thyroid cancer) and confirm
that the monotonicity-in-parallel-component property reproduces in
12/12 cohort pairings. We also introduce DIAL-Adapt, a corrective
method that automatically minimises the parallel component when DIAL
> 0.1, and show it outperforms naive ComBat by 0.18 AUC on average.
```

This is a **placeholder**. v16 abstract will be rewritten in 2027
based on what actually gets proven and benchmarked.

## Action items if v15 lands acceptance

1. Within 1 week of acceptance: post non-anonymous arXiv preprint
   (citable for v16's Related Work section)
2. Within 1 month: start v16 theorem work
3. Within 6 months: have a clean draft of v16 with non-Gaussian
   theorem proven OR v16-fallback positioned as workshop paper
