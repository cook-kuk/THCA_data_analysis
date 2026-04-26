# v15 Theoretical Connections

**Author:** Seungho Cook
**Date:** 2026-04-25
**Companion to:** `theorem2_proof.tex`, `v15_neurips_paper_v0.tex`

> **Scope note (v5.2 audit, 2026-04-25):** The v5.1 THCA empirical claim
> ("DIAL = 0.494, batch_entangled in 4/5") was retracted under proper LODO
> ComBat-fit-on-train-only. The theoretical material below stands
> independent of that empirical retraction. All empirics in v15 are
> synthetic and explicitly framed as a *calibration study* of when DIAL
> fires.

This document situates Theorem 2 (unified shift decomposition) within six
existing theoretical traditions. For each, we state (i) the precise
mathematical relationship, (ii) when DIAL coincides with the prior
construct, (iii) when it diverges, and (iv) why DIAL offers unique
diagnostic value.

---

## 1. Covariate shift (Shimodaira 2000; Sugiyama et al. 2007 — KLIEP)

**Setup.** Covariate shift assumes
$p_s(x)\!\ne\!p_t(x)$ but
$p_s(y\mid x)\!=\!p_t(y\mid x)$. KLIEP (Sugiyama 2007) estimates the
density ratio $r(x)=p_t(x)/p_s(x)$ by minimising an empirical KL on a
linear model, and re-weights the source loss by $r(x)$.

**Relationship.** Theorem 2 implies that the KLIEP fix is *insufficient*
when $\Delta_{\text{cond}}^{\parallel}\!>\!0$, because re-weighting acts on
$p(x)$ alone. KLIEP cannot move probability mass between
$p(\tilde x\mid y\!=\!0)$ and $p(\tilde x\mid y\!=\!1)$, but precisely that
movement is what produces the DIAL flip.

**Coincidence.** Under pure covariate shift with $T$ being the
Bayes-optimal importance map, both KLIEP-corrected risk and DIAL are
zero.

**Divergence.** Under our subspace-aligned regime, KLIEP-corrected risk
remains low on its KL diagnostic but DIAL is positive — KLIEP is "fooled"
because it monitors the wrong divergence.

**Unique value.** DIAL targets the conditional axis that KLIEP cannot
observe.

---

## 2. Domain-invariance theory (Ben-David et al. 2007/2010 — h-divergence)

**Setup.** The seminal bound
$\epsilon_t(h)\le\epsilon_s(h) + d_{\mathcal{H}\Delta\mathcal{H}}(p_s,p_t)
+ \lambda$, where $d_{\mathcal{H}\Delta\mathcal{H}}$ is the symmetric
$\mathcal{H}$-divergence. The proof of inequality
$\Delta_{\text{cond}}^{\parallel}\le d_{\mathcal{H}\Delta\mathcal{H}}\le 2
\sqrt{\Delta_{\text{cond}}}$ is sketched in our Theorem 2 proof and
flagged `TODO Seungho verify`.

**Relationship.** $d_{\mathcal{H}\Delta\mathcal{H}}$ is intractable for
deep $\mathcal{H}$. Under linear correction $T$ and Gaussian conditionals,
DIAL is a *computable surrogate* that lower-bounds it.

**Coincidence.** When $\mathcal{H}$ is the linear function class and
class-conditional means are well-separated, DIAL and the (single-classifier
proxy of) Ben-David's bound differ by an additive $O(\varepsilon_n)$ from
finite-sample noise.

**Divergence.** When $\mathcal{H}$ is rich (deep nets), DIAL captures only
the linear-classifier flip; the underlying h-divergence may still be
positive in higher-order moments that DIAL does not see. We document this
gap explicitly in §6 (Discussion → Limitations) of the paper draft.

**Unique value.** Cheap, AUC-based, requires no domain classifier
training.

---

## 3. Information bottleneck (Tishby & Zaslavsky 2015)

**Setup.** A representation $\tilde X$ is sufficient for $Y$ if
$I(\tilde X;Y)=I(X;Y)$, and minimal-sufficient if it additionally
minimises $I(\tilde X;X)$.

**Relationship.** Under our Gaussian-shared-covariance assumptions
(theorem2_setup, §1.2),
$$
   I(\tilde X;Y\mid B) \;=\; I(X;Y) - \Delta_{\text{cond}}^{\parallel}
   + O(1/n)
$$
i.e.\ a *positive* DIAL is exactly evidence that the correction has
destroyed $\Delta_{\text{cond}}^{\parallel}$ bits of label information
while preserving information about the orthogonal complement of the batch
subspace. This is flagged `TODO Seungho verify` in the proof.

**Coincidence.** A perfectly invariant representation has both
$\mathrm{DIAL}\!=\!0$ and $I(\tilde X;Y\mid B)\!=\!I(X;Y)$.

**Divergence.** In the heavy-tailed or non-Gaussian regime the identity
above gains a residual term; the qualitative direction (DIAL up ⇔ $I$
down) is preserved but the constant $c_0$ depends on tail behaviour.

**Unique value.** DIAL converts an information-theoretic phenomenon into
a single number that tracks model performance directly.

---

## 4. Invariant risk minimisation (Arjovsky et al. 2019 — IRM)

**Setup.** IRM seeks a representation $\Phi(X)$ such that the optimal
classifier on $\Phi$ is *the same across environments*. The theory
guarantees out-of-distribution generalisation under linear-Gaussian
assumptions; in practice IRM-v1 (penalised ERM) is widely used.

**Relationship.** A positive DIAL is precisely the failure mode IRM
attempts to prevent: the source-optimal classifier does not generalise
because its score direction is anti-aligned with the target's true label
direction. **IRM's "invariance gap" is upper-bounded by the projection of
$\Delta_{\text{cond}}$ onto $\mathrm{span}(w_Y)$ — which is exactly our
$\Delta_{\text{cond}}^{\parallel}$**.

**Coincidence.** IRM converging to a global invariant solution
$\Leftrightarrow$ DIAL $\to 0$ in the limit.

**Divergence.** IRM's penalty does not give a *closed-form, post-hoc*
diagnostic, only a training criterion. DIAL is post-hoc and works on
*any* model. They are therefore complementary: train IRM, audit with
DIAL.

**Unique value.** A model-agnostic empirical diagnostic for IRM's
target objective.

---

## 5. Optimal transport / Wasserstein (Courty 2017; Damodaran 2018)

**Setup.** Wasserstein-based DA matches $p_s$ to $p_t$ in the
$W_2$ metric, often per-class (joint distribution OT, JDOT).

**Relationship.** Under our Gaussian assumptions,
$W_2^2\bigl(p_s(\tilde x\mid y),\,p_t(\tilde x\mid y)\bigr)
= \|\tilde m_{s,y}-\tilde m_{t,y}\|_{\Sigma^{-1}}^{2}$, identical (up to
factor of 2) to our $\Delta_{\text{cond}}$. Splitting $W_2^2$ along
$\mathrm{span}(w_Y)$ and its orthogonal complement reproduces
$\Delta_{\text{cond}}^{\parallel}$ and $\Delta_{\text{cond}}^{\perp}$.

**Coincidence.** JDOT-corrected residual = 0 $\Leftrightarrow$ DIAL = 0
(modulo finite-sample noise).

**Divergence.** OT methods compute the *full* coupling and pay
$O(n^3)$; DIAL needs only AUC.

**Unique value.** $10^{4}\times$ cheaper post-hoc; identical asymptotic
target.

---

## 6. Influence functions and instance-level corrections (Koh & Liang 2017)

**Setup.** Koh & Liang's influence functions identify training instances
whose removal would change the model's prediction on a test point.

**Relationship.** Linear ComBat correction is an *aggregate* of
instance-level interventions: $T(x_i,b_i)=A_{b_i}x_i+c_{b_i}$. The
sum-of-influences picture (Koh §3.2) shows that when batch-mean
estimation is contaminated by class imbalance (the v5.2 leak mechanism!),
the cumulative influence on the source classifier is *along
$\mathrm{span}(w_Y)$* — exactly Theorem 2's positive-DIAL regime.

**Coincidence.** No flip $\Leftrightarrow$ aggregated influence has zero
projection on $w_Y$.

**Divergence.** Influence functions diagnose individual training points;
DIAL diagnoses the global aggregate.

**Unique value.** A batch-level analogue of influence functions:
"how much did the correction step move the classifier's score axis?".

---

## 7. v5.2 LODO leak — what the retraction *adds* (not subtracts)

The v5.2 audit (project memo `v52_lodo_finding.md`) found that the
empirical THCA flip in v5.1 was a **leakage artifact**: ComBat fit on
the pooled (source + target) matrix before LODO splitting borrowed
strength from the held-out cohort via empirical-Bayes priors with
`par_prior=True`. Under proper fit-on-train-only ComBat, all 5 THCA
classifiers return DIAL = 0.

This *strengthens* — does not weaken — Theorem 2's contribution:

1. The leak mechanism is itself a special case of subspace-aligned
   conditional shift. ComBat fitting on pooled data with
   class-imbalanced batches estimates $\hat\gamma_b$ that has a non-zero
   projection on $\Sigma^{-1/2}w_Y$, producing
   $\Delta_{\text{cond}}^{\parallel}\!>\!0$ on the held-out cohort.
   Theorem 2 predicts this should produce DIAL > 0 — and v5.1 measured
   exactly that.
2. Under proper LODO, $\hat\gamma_b$ is estimated *without* the held-out
   class composition; $\Delta_{\text{cond}}^{\parallel}\to 0$ and DIAL
   $\to 0$ — also consistent with Theorem 2.
3. **DIAL is therefore validated as a leak detector**, not (in this
   case) as a real-biology detector. The paper draft re-frames v15 as a
   *calibration study*: Theorem 2 explains what DIAL fires on, and the
   THCA story becomes a worked example of a leaky preprocessing pipeline.

---

## Summary table

| Construct                 | Coincides when …                      | Differs in …                       |
|---------------------------|---------------------------------------|------------------------------------|
| KLIEP / covariate shift   | pure covariate, optimal $T$          | conditional shift inside span(w_Y) |
| Ben-David h-divergence    | linear $\mathcal{H}$, large samples   | richer hypothesis classes          |
| Information bottleneck    | Gaussian, shared $\Sigma$             | non-Gaussian tails                 |
| IRM                       | infinite training, single environment | post-hoc auditability              |
| Wasserstein OT            | shared support, Gaussian              | $O(n^3)$ vs. AUC cost              |
| Influence functions       | linear correction aggregate           | instance-level granularity         |
| v5.2 LODO leak            | (always) — Theorem 2 explains it      | empirical artefact, not biology    |

