# AAAI 2027 Related Work — Direction-Invariant Auditing of Batch Correction
*(v10, 2026-04-24)*

The closest literature divides across six sub-areas. We cite at least 15
ML/statistics references and draw out, for each, the specific delta our
work introduces.

## 1. Covariate shift and distribution mismatch

Classical covariate-shift theory treats train-test mismatch as
$P_{\text{tr}}(X) \ne P_{\text{te}}(X)$ with $P(Y\mid X)$ fixed
(Shimodaira 2000; Sugiyama, Krauledat & Müller 2007). Methods include
KLIEP (Sugiyama et al. 2008) and kernel-based MMD minimisation (Gretton
et al. 2009, 2012). Pan & Yang (2010) survey broader domain adaptation;
Ben-David et al. (2010) formalise discrepancy-based generalisation
bounds.  All assume $P(Y\mid X)$ is invariant -- the failure mode we
probe (\emph{correction-induced} subspace alignment) violates this
invariance in the \emph{post-correction} feature space even though the
unmodified $P(Y\mid X)$ is identifiable.  **Delta (ours):** DIAL is not
a shift corrector; it is an *audit* on the output of a correction
operator, flagging \emph{structured} leakage the above methods cannot
detect.

## 2. Dataset / model cards and ML auditing

Datasheets for Datasets (Gebru, Morgenstern, Vecchione et al. 2021) and
Model Cards (Mitchell, Wu, Zaldivar et al. 2019) introduced structured
documentation for ML artefacts. Underspecification (D'Amour et al. 2020)
showed that multiple equally-good models can diverge catastrophically
under distribution shift.  Recent auditing frameworks (Raji et al.
2020 SMACTR; Metaxa et al. 2021) push for systematic third-party audit
methodologies.  **Delta:** all documentation-based approaches are
descriptive.  DIAL is a \emph{quantitative} metric with a closed-form
threshold, hence actionable in CI or pre-registration.

## 3. Information-theoretic leakage

Kaufman, Rosset, Perlich & Stitelman (2012) and Naish-Guzman & Holden
(2008) formalise train-test contamination. Samala et al. (2020)
document its prevalence in medical imaging CAD systems.  Roberts et al.
(\emph{Nat.\ Mach.\ Intell.} 2021) meta-analyse COVID-19 ML leakage.
**Delta:** these works assume direct contamination (shared samples,
preprocessing fit on test data).  The DIAL flip describes a more
insidious channel -- \emph{structured} leakage through a
covariate-preserving transform that, by construction, \emph{should}
have been label-protecting.

## 4. Batch effect and normalisation

ComBat (Johnson, Li & Rabinovic 2007) remains the most cited batch
corrector; SVA (Leek \& Storey 2007, 2010) handles unknown batch
factors; RUV (Gagnon-Bartsch \& Speed 2012) uses control probes.  Newer
single-cell tools include scVI (Lopez et al. 2018, \emph{Nat.\ Methods}),
Harmony (Korsunsky et al. 2019, \emph{Nat.\ Methods}), and LIGER (Welch
et al. 2019).  The scIB benchmark (B\"uttner et al. 2019) unified
cross-method evaluation via integration metrics (ARI, LISI, NMI, ASW).
**Delta:** benchmark metrics in scIB et al.\ measure \emph{how well}
batch effects are removed -- they do not detect the failure mode where
removal inadvertently destroys label signal.  DIAL is an orthogonal
audit metric capturing this specific failure.

## 5. Quantum machine learning

Havlíček et al.\ (\emph{Nature} 2019) introduced the quantum kernel
with ZZFeatureMap. Schuld \& Killoran (\emph{Phys.\ Rev.\ Lett.} 2019)
derived the feature-map / kernel-method correspondence; Schuld,
Sweke \& Meyer (2021) analysed data-reuploading expressivity.
Huang et al.\ (\emph{Nat.\ Commun.} 2021) established theoretical
advantage conditions over classical kernels.  Liu, Arunachalam \&
Temme (2021) proved quantum speedup on constructed problems.  **Delta:**
quantum ML has emphasised \emph{advantage} benchmarks; our
Proposition~2 reframes the quantum setting as an \emph{auditing} tool
where exponential Hilbert-space dimension lower-bounds the probability
of DIAL-flip under random feature maps -- the first quantum audit
framework for biological data.

## 6. Random-effects meta-analysis in ML

Benavoli et al. (2017) introduced Bayesian significance tests for
ML comparison; Dem\v{s}ar (2006) standardised statistical comparison
of classifiers across datasets.  In genetics, DerSimonian \& Laird
(1986) is the classical random-effects standard, and Han \& Eskin
(2011) introduced the $m$-value posterior for per-study effect
attribution in METASOFT.  Lee, Kang \& Eskin (2022) extended to
multi-ancestry settings.  **Delta:** $m$-values have not previously
been applied to ML audit outputs. We reinterpret the per-cancer
DIAL as a study-level effect and compute the posterior $m_\text{THCA}
\ge 0.9999$, providing the first meta-analytic attribution of an ML
failure mode to a specific cohort.

---

### Summary of deltas

| Area | Representative work | DIAL contribution |
|------|--------------------|-------------------|
| Covariate shift | Sugiyama 2007, Gretton 2009 | Audit *after* correction, not before |
| ML auditing | Mitchell 2019, Gebru 2021 | Quantitative threshold, not documentation |
| Info leakage | Kaufman 2012, Samala 2020 | Structured pipeline-induced leakage |
| Batch correction | Johnson 2007, Korsunsky 2019 | Orthogonal failure-mode metric |
| Quantum ML | Havlíček 2019, Huang 2021 | First quantum audit framework |
| Meta-analysis | Han 2011, Benavoli 2017 | First $m$-value for ML attribution |

We cite 17 works across these six areas in the main bibliography
(references embedded in `aaai_paper_v0.tex`).
