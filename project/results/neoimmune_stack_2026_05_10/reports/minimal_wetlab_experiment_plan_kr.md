# Minimal Wet-Lab Experiment That Makes This Paper Real KR

## Cohort
5-10 patients first, then 30-50 for paper-grade patient-level metrics.

## Per patient
- HLA class I typing
- somatic mutation table
- mutant and WT peptides
- expression TPM
- VAF/clonality
- HLA LOH
- B2M/TAP/APM context

## Candidate selection
- Top 20 by production stack
- Top 10 clean-local rescue candidates
- 5 binding-only high-score negative-risk controls
- 5 model-disagreement controls

## Assays
- MS immunopeptidomics if feasible for presentation evidence
- T-cell assay for immunogenicity evidence
- Track synthesis success and assay failure separately

## Endpoints
- Primary: patient-level hit rate@20 and Recall@20
- Secondary: strict no-overlap AUPRC and rescued-positive fraction

## Claim after success
Patient-context-aware ranking improves wet-lab candidate prioritization. Still not a clinical efficacy claim.
