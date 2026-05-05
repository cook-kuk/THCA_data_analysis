# Paper 9 Thyroid Model Coverage Audit

## Coverage Summary
- thyroid_flagged_models: 25 
- strict_thyroid_lineage_models: 22 OncoTree primary/subtype/lineage thyroid; excludes thyroid-site non-thyroid cancers
- crispr_overlap: 14 
- strict_thyroid_crispr_overlap: 11 
- ccle_expression_overlap: 25 
- strict_thyroid_expression_overlap: 22 
- prism_overlap: 16 all PRISM 24Q2 models, not GLS-specific
- gdsc_ctrp_overlap: 15 any GDSC1/GDSC2/CTRP AUC model
- gls_expression_dependency_drug_complete: 9 GLS expression + GLS dependency + matched glutaminase inhibitor response
- strict_gls_expression_dependency_drug_complete: 9 strict thyroid subset with GLS expression + GLS dependency + matched glutaminase inhibitor response

## Histology / Subtype Annotation
### OncotreePrimaryDisease
- Anaplastic Thyroid Cancer: 11
- Well-Differentiated Thyroid Cancer: 8
- Medullary Thyroid Cancer: 2
- Head and Neck Squamous Cell Carcinoma: 2
- Sarcoma, NOS: 1
- Poorly Differentiated Thyroid Cancer: 1
### OncotreeSubtype
- Anaplastic Thyroid Cancer: 11
- Follicular Thyroid Cancer: 5
- Papillary Thyroid Cancer: 3
- Medullary Thyroid Cancer: 2
- Head and Neck Squamous Cell Carcinoma: 2
- Sarcoma, NOS: 1
- Poorly Differentiated Thyroid Cancer: 1
### OncotreeCode
- THAP: 11
- THFO: 5
- THPA: 3
- THME: 2
- HNSC: 2
- SARCNOS: 1
- THPD: 1
### SampleCollectionSite
- thyroid: 13
- pleural_cavity: 6
- lymph_node: 2
- pleural_effusion: 2
- soft_tissue: 1
- lung: 1

## Boundary
- Thyroid-only estimates remain underpowered and are used as coverage/audit evidence only.
- No protected data or new large downloads were used.
