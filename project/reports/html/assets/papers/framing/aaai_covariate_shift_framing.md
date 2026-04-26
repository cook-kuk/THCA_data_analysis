# DIAL in the Covariate-Shift Taxonomy
*(v10 AAAI 2027 -- 2026-04-24)*

## Standard taxonomy

Let $P_{\text{tr}}$ and $P_{\text{te}}$ be the joint train and test
distributions over $(X,Y)$.  The canonical classification (Quionero-
Candela, Sugiyama, Schwaighofer \& Lawrence 2009) is:

- **Covariate shift:** $P_{\text{tr}}(X) \ne P_{\text{te}}(X)$,
  $P(Y\mid X)$ invariant. Correctable with importance weighting.
- **Label shift (prior shift):** $P(X\mid Y)$ invariant,
  $P(Y)$ differs.  Corrected by reweighting per-class priors.
- **Concept shift:** $P(Y\mid X)$ differs. The hardest case.
- **Conditional shift:** general $P(X\mid Y)$ and $P(Y)$ both differ.

## Our setting: batch-induced pseudo-conditional shift

In multi-cohort biomarker development, each cohort $B\in\{1,\dots,K\}$
has its own conditional
\[
   P_B(X\mid Y) \;=\; P(X\mid Y, B).
\]
The marginal $P(X\mid Y) = \sum_B P_B(X\mid Y)\Pr(B)$ is intermediate.
Standard batch correction (ComBat, Harmony, scVI) aims to remove
$P(X\mid B)$, not $P(X\mid Y, B)$.

### The failure mode

Let $w_Y, w_B$ be the dominant label and batch directions (per
Theorem~1). When the subspaces $\mathrm{span}(w_Y)$ and
$\mathrm{span}(w_B)$ are \emph{sufficiently aligned} (cosine-squared
$\rho^2 > 1/2$), linear correction operators cannot factor
$P(X\mid Y,B)$ cleanly. They instead produce a corrected
$\tilde P(X\mid Y)$ whose support-density ordering along $w_Y$
\emph{inverts} relative to the pre-correction conditional density.

We call this \textbf{subspace-aligned pseudo-conditional shift}. It is
not covariate shift (marginals are equalised by construction). It is
not label shift. It is a flavour of \emph{concept shift} that is
induced \emph{by} the attempted harmonisation, not a pre-existing
property of the data.

## Relation to existing methods

\begin{center}
\begin{tabular}{lll}
\toprule
Method & Target shift removed & Audit metric \\
\midrule
IWA  (Shimodaira 2000)      & $P(X)$ shift               & density ratio \\
ComBat (Johnson 2007)       & per-feature batch mean/var & none \\
SVA  (Leek 2007)            & unknown surrogate factors  & none \\
scVI (Lopez 2018)           & latent-variable batch      & ARI / KBET \\
Harmony (Korsunsky 2019)    & per-PC batch clustering    & LISI \\
scIB (Büttner 2019)         & benchmark suite            & integration metrics \\
\textbf{DIAL (ours)}        & \emph{subspace alignment}  & direction-invariant AUC \\
\bottomrule
\end{tabular}
\end{center}

Every method to the left of DIAL measures \emph{how completely} batch
structure is removed. None measures \emph{what collateral damage} the
removal has done to the \emph{label} signal. DIAL closes that gap.

## Theoretical connection to importance weighting

ComBat with Y covariate can be derived as a constrained importance
reweighting problem. Write the empirical risk as
\[
  R(\theta) \;=\; \E_{X\sim P_{\text{tr}}(\cdot)}\!\left[
     r(X,Y;\theta)\cdot\frac{P_{\text{te}}(X,Y)}{P_{\text{tr}}(X,Y)}\right].
\]
When $P_{\text{tr}}$ and $P_{\text{te}}$ differ only in per-batch
location, IWA with a linear log-density ratio reproduces the ComBat
update for the additive-Gaussian model~\eqref{eq:additive-model}.
\TODO{Make this derivation rigorous: match the ComBat eBayes-shrunken
update to a KLIEP-style LASSO solution for the density ratio; see
Sugiyama, Suzuki \& Kanamori 2012 eq.~3.4.}

Under this lens:

- The \emph{batch-pseudo-conditional} structure we describe is exactly
  the case where the assumed parametric form of the density ratio is
  \emph{misspecified} in a specific way -- it cannot disentangle the
  $w_Y$ and $w_B$ contributions when they are collinear.
- DIAL flags this misspecification empirically.
- Proposition~2 (quantum bound) then shows how lifting to an
  exponentially-large feature space mitigates the misspecification in
  a probabilistic sense.

## Takeaway

The pseudo-conditional shift we study is a specific, detectable failure
mode of linear importance-weighting / batch-correction hybrids in
cohort-rich biomedical ML. DIAL is the audit metric for this failure
mode; Theorem~1 characterises when it fires; Proposition~2 quantifies
the gain from quantum feature-space lifting.

Word count: ~620.
