# Presubmission Inquiry Pack

## Universal Subject Line
Presubmission inquiry: reporting standard for vaccine-ready neoantigen AI

## Nature Medicine / Nature Cancer Version
Dear Editors,

We are preparing a Perspective proposing NEO-PRIOR, a reporting and evaluation standard for AI systems that rank neoantigen candidates for personalized cancer vaccine manufacture, immune testing and clinical prioritization.

The motivation is that candidate selection for personalized vaccines now depends on computational shortlists, yet published evaluations often mix binding, presentation, immunogenicity and triage tasks, while public-data overlap, near-peptide memorization, source/study shift, negative-label ambiguity, calibration and top-k utility are not reported consistently.

We propose a minimum reporting checklist and a clean separation between scientific validation and product-assisted triage. A companion benchmark framework, NeoBench-Vax, illustrates how models can be evaluated using locked splits, overlap audits, model cards and abstention-aware metrics. We do not claim external validation or universal model superiority; instead, we argue that the field needs an auditable evaluation contract before candidate-ranking AI can be interpreted clinically.

We would be grateful for your view on whether this topic would be suitable as a Perspective or Comment.

## Nature Biomedical Engineering Version
Dear Editors,

We are preparing a manuscript on NeoBench-Vax, a benchmark and validation infrastructure for vaccine-ready neoantigen prioritization systems, together with CROSS-Neo as a first contamination-controlled implementation.

The work treats neoantigen prioritization as a biomedical decision-support system rather than a single-score predictor. The benchmark includes locked exact, near-peptide, source/study, HLA and retrieval-clean splits; public-overlap audits; model-card requirements; AUPRC/top-k/enrichment as primary metrics; calibration and OOD/abstention reporting; and a clean separation between public predictors as comparators versus product-assisted features.

Our claim boundary is conservative: internal locked retrospective evaluation only, no external validation claim, and no quantum-advantage claim. The central contribution is a reproducible evaluation framework for computational systems that may influence assay and manufacturing shortlists in personalized cancer vaccine workflows.

We would welcome guidance on whether this is better suited as an Article or Perspective.

## Nature Machine Intelligence Version
Dear Editors,

We are preparing an analysis of AI evaluation failure modes in neoantigen prioritization for personalized cancer vaccines. The manuscript focuses on leakage-sensitive benchmarks, public predictor reuse, ambiguous negative labels, source/study shift, uncertainty, abstention and top-k decision metrics.

We introduce NEO-PRIOR as a reporting checklist and NeoBench-Vax as a benchmark framework, with CROSS-Neo as a first implementation. The broader AI contribution is an evaluation contract for small-n biomedical ranking problems where the practical objective is a short, high-precision candidate list rather than a calibrated population-level classifier.

We would appreciate your advice on fit for an Analysis, Comment or related format.
