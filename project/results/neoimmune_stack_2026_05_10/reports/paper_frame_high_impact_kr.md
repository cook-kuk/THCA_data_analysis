# High-Impact Paper Frame KR

## Title
CLEAN-Neo++: Leakage-aware integration of TCR-visible immunogenicity, structure proxies, and public presentation predictors for patient-level cancer vaccine candidate ranking

## Figure Flow
1. Problem: binding/presentation alone is not immunogenicity.
2. System: clean local algorithm track vs production frozen-predictor stack.
3. Leakage: exact peptide, peptide-HLA, study, patient, HLA allele, public-tool training risk.
4. Clean result: local algorithms under strict no-existing-overlap.
5. Production result: practical stack compared with BigMHC/MHCflurry/PRIME/NetMHCpan artifacts.
6. Failure audit: false positives and rescued positives.
7. Patient handoff: top-20 report format and metadata requirements.

## Reviewer-Proof Main Claim
We do not claim clinical efficacy. We show that candidate ranking needs a leakage-aware, patient-context-aware operating layer and that local immunogenicity branches can be evaluated separately from public presentation predictors.

## Must Not Say
- “validated vaccine candidates”
- “predicts clinical response”
- “presentation confirmed” without MS
- “immunogenicity confirmed” without T-cell assay
- “LLM predicts vaccine efficacy”
