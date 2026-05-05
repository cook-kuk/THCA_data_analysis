# P3/P9 Execution Boundary Summary

## Paper 3 result summary
- Immune score rows: 815
- Joined lineage/immune rows: 815
- PASS signatures: ANTIGEN_PROCESSING, CXCL13_AXIS, CYTOLYTIC, DEDIFFERENTIATED_EPI, IFNG_TIS, IFN_GAMMA_RESPONSE, M2_TAM, MDSC_LIKE, MHC1_CORE, MHC2_CORE, STROMAL_EXCLUSION, TLS_CABRITA, TREG, T_CELL_EXHAUSTION
- Output is an ICI-readiness atlas, not an ICI response predictor.

## Paper 9 result summary
- Thyroid DepMap models: 25
- Candidate dependency rows: 26
- Drug association rows: 140
- Output is a cell-intrinsic vulnerability map, not clinical utility evidence.

## Overlap genes/pathways
- Direct signature-target overlap: none
- Shared pathway themes: interferon/JAK-STAT context, dedifferentiation/stromal programs, chromatin and survival dependencies.

## Difference
- Paper 3 = immune readiness and tumor microenvironment state.
- Paper 9 = cell-intrinsic dependency and drug-response nomination from cell-line data.

## Claim guard
- No treatment recommendation.
- No patient treatment selection.
- No ICI response predictor claim.
- No TROP2 synthetic-lethality claim.

## Next gates
- Protected-access planning only if local credentials/files are explicitly supplied.
- Wet-lab validation planning for Paper 9 candidates with non-high artifact risk.
- Manuscript prose remains untouched.
