# Paper 3 ICI Full-Public Report

## 1. Data used
- Local/public expression signature score rows: 815
- Lineage-immune joined rows: 815
- TCGA immune scores used only from already-local processed immune signature table when available.
- HLA LOH and neoantigen data were not imputed.

## 2. Signature registry
- Registry output: `/home/seungho/personal/THCA_data_analysis/project/results/p3_p9_full_execution/paper3/signature_registry.tsv`
- Signatures registered: 14

## 3. DIAL audit result
- PASS: 14
- FLIP: 0
- COLLAPSE: 0
- AMBIGUOUS: 0

## 4. ICI readiness score
- Terms used after PASS filter: IFNG_TIS, CYTOLYTIC, MHC1_CORE, MHC2_CORE, TLS_CABRITA, CXCL13_AXIS, TREG, M2_TAM, MDSC_LIKE, STROMAL_EXCLUSION
- Scored samples: 778

## 5. Immune ecotypes
- Ecotype rows: 815
- ecotype_4: 354
- ecotype_2: 241
- unassigned: 96
- ecotype_1: 71
- ecotype_3: 53

## 6. HLA/neoantigen blockers
- See `03_paper3_protected_access_blockers.md`.

## 7. Claim boundary
- Allowed: ICI readiness, vulnerability, and research prioritization language.
- Forbidden: ICI response predictor, treatment recommendation, patient treatment selection.

## 8. Figure list
- `/home/seungho/personal/THCA_data_analysis/project/results/p3_p9_full_execution/paper3/fig_p3_immune_heatmap.png`
- `/home/seungho/personal/THCA_data_analysis/project/results/p3_p9_full_execution/paper3/fig_p3_lineage_vs_ici_readiness.png`
- `/home/seungho/personal/THCA_data_analysis/project/results/p3_p9_full_execution/paper3/fig_p3_dial_audit.png`

## 9. Next gates
- Add protected HLA LOH or neoantigen analyses only if local/credentialed inputs are explicitly provided.
- Treat weak or negative immune associations as honest negative results.
