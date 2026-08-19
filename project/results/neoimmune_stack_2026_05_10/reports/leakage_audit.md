# Leakage audit

## Controls implemented
- exact peptide leakage
- peptide-HLA pair leakage
- mutant-WT pair leakage
- source protein/window leakage where source windows are available
- study leakage
- patient leakage
- HLA allele leakage
- external tool training contamination risk flags when provided by existing artifacts

## Summary
- existing_overlap_flags: 2,327/22,090 (0.105)
- source_protein_window_available: 0/22,090 (0.000)

## Reviewer boundary
Random split performance is smoke-test evidence only. Clean-science claims should be read from strict, source-heldout, patient-heldout, or no-overlap views.
