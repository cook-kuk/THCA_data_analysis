# Failure case audit

## False-positive annotations
- low expression
- low VAF/subclonal risk
- weak or missing presentation evidence
- high WT/self similarity when available
- poor TCR/self-similarity branch evidence when available
- HLA LOH/APM/B2M risk
- model disagreement
- source-study/leakage risk
- external predictor conflict

## Rescued positive annotations
- missed by binding-only/public comparators when available
- rescued by TCR/self-similarity branch
- rescued by structure branch
- rescued by quantum kernel branch
- rescued by production stack
- rescued by patient context

- Failure/rescue rows written: 139

## Boundary
This audit explains ranking behavior; it does not assert biological validation.
