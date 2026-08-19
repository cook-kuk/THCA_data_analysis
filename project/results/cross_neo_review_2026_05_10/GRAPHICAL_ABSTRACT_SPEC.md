# Graphical Abstract Spec

## Concept

The figure should show that the field is moving from "predict binders" to "prioritize manufacturable vaccine candidates."

## Layout

### Left: Clinical Funnel

Patient tumor sequencing -> mutation calling -> candidate peptides -> HLA presentation -> immunogenicity evidence -> manufacturable vaccine shortlist.

Label the bottleneck:

> Thousands of candidates become 5-20 practical vaccine slots.

### Middle: Why Current Benchmarks Break

Four warning icons:
- exact peptide-HLA overlap;
- near-peptide/source-protein overlap;
- HLA/source/study shortcut;
- incomplete negative labels.

Caption:

> A high AUROC leaderboard can still fail the vaccine shortlist.

### Right: Vaccine-Ready Evaluation Contract

A checklist panel:
- AUPRC
- top5/top10 precision
- enrichment over prevalence
- exact/near/source/HLA/study/time holdouts
- train-only retrieval audit
- calibration and OOD abstention

Final badge:

> Rank the candidates worth manufacturing.

## Style

- Clinical-translational, not tech-dashboard.
- Avoid quantum imagery.
- Avoid implying external validation.
- Use restrained colors: blue/teal for clinical workflow, amber for benchmark risks, green for evaluation contract.

## Caption Draft

Personalized neoantigen vaccine design is a constrained ranking problem. Candidate peptides pass through biological filters, but only a small shortlist can be manufactured or assayed. Evaluation therefore needs contamination-controlled splits, source-aware stress tests and top-k metrics that reflect vaccine decision-making rather than AUROC alone.
