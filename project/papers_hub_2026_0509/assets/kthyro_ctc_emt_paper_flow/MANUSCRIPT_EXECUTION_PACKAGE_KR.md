# Execution Package


## Main/Supp Figure Layout

| paper_slot           | content                              | asset                          | manuscript_position         | takeaway                                                         |
|:---------------------|:-------------------------------------|:-------------------------------|:----------------------------|:-----------------------------------------------------------------|
| Main Fig 1           | Study logic and data provenance      | PF01 + PF02 combined           | Introduction/overview       | public data as prior, Yu as feasibility, hospital as future test |
| Main Fig 2           | Mutation/state heterogeneity         | F03 + F04                      | Results 1                   | BRAF/mutation alone is insufficient                              |
| Main Fig 3           | Spatial tissue-state organization    | F05 + source spatial coherence | Results 2                   | tissue state is organized, not random                            |
| Main Fig 4           | Marker context and module support    | F06 + F07 + PF10               | Results 3-4                 | marker panel is biologically motivated but not validated         |
| Main Fig 5           | CTC-EMT-NGS interpretation framework | PF03 + F08                     | Results 5                   | phenotype must be clone-anchored                                 |
| Main Fig 6           | Dataset gap and future unlock        | PF04 + PF12                    | Discussion                  | current paper stops before discovery claims                      |
| Main Fig 7           | Claim boundary and reviewer defense  | PF05 + PF07                    | Discussion/Reviewer defense | explicit allowed/forbidden claim boundary                        |
| Supplementary Fig S1 | CV2 visual atlas                     | PF08                           | Supplement/meeting          | all visual assets in one board                                   |
| Supplementary Fig S2 | Publishability decision              | F09                            | Internal supplement         | paper route decision                                             |
| Supplementary Tables | All data dictionaries and manifests  | TSV tables                     | Methods/Supplement          | reproducible provenance                                          |


## Methods Reproducibility Checklist

| analysis_unit        | input_files                                                    | primary_output   | method_or_metric                     | status        | manuscript_qc_note                             |
|:---------------------|:---------------------------------------------------------------|:-----------------|:-------------------------------------|:--------------|:-----------------------------------------------|
| TCGA BRAF split      | tcga_braf_only_label_distribution.tsv                          | F03              | BRAF tumors across state labels      | computed      | confirm exact n/percent before manuscript      |
| TCGA driver variance | tcga_axis_variance_models.tsv                                  | F04              | driver-only model boundary           | computed      | report as descriptive only                     |
| Spatial coherence    | spatial_coherence_statistics.tsv                               | F05              | same-niche z-score across slides     | computed      | cite GSE250521 and spot count                  |
| scRNA marker context | scrna_module_celltype_attribution.tsv                          | F06              | module-cell context                  | computed      | state source dataset in Methods                |
| Proteomics proxy     | bulk_proteomics_kthyro_module_contrasts.tsv; dediff_trends.tsv | F07              | protein-direction support            | computed      | call it bulk proxy, not spatial/CTC proteomics |
| Marker support score | marker_module_cross_layer_support_matrix.tsv                   | PF10             | 0/1/2 planning support scores        | new synthesis | qualitative planning score only                |
| Claim boundary       | claim_boundary_matrix.tsv                                      | PF05             | allowed vs forbidden statements      | computed      | use in Discussion and reviewer memo            |
| Future unlock        | future_data_unlock_map.tsv                                     | PF12             | hospital data unlocks future figures | new synthesis | do not mix with current Results                |


## Hospital Data Request CRF

| data_block       | fields                                                               | scope                      | why_needed                    | priority                        |
|:-----------------|:---------------------------------------------------------------------|:---------------------------|:------------------------------|:--------------------------------|
| patient_core     | patient_id, age, sex, diagnosis date, surgery date                   | all patients               | de-identified linkage key     | required                        |
| tumor_pathology  | tumor size, multifocality, LVI, LNM, ETE, margin, ATA risk           | resected tumor             | pathology report fields       | required                        |
| blood_timepoints | T0 pre-op, T1 2 weeks, T2 3 months, T3 6-12 months/recurrence        | serial blood               | exact draw date and tube type | required                        |
| CTC_counts       | total CTC, epithelial, hybrid, mesenchymal counts/fractions          | each blood timepoint       | assay method and QC flag      | required                        |
| CTC_marker_qc    | marker intensity, threshold, batch, image/cell QC if available       | cell-level or sample-level | needed for reproducibility    | high priority                   |
| tumor_normal_NGS | variant list, VAF, coverage, CNV/fusion, normal filter               | tumor-normal pair          | clone anchor                  | required                        |
| cfDNA_NGS        | variant, VAF, depth, fragment/QC, timepoint                          | serial plasma              | low signal expected           | high priority                   |
| followup         | Tg, anti-Tg, ultrasound, RAI, recurrence/persistence, last follow-up | post-op clinical follow-up | date-stamped                  | required for future association |
| sample_logistics | blood volume, processing time, storage, shipping, freeze-thaw        | all biospecimens           | failure-mode audit            | required                        |


## Reviewer Risk Heatmap

| reviewer_risk                              |   severity_1_5 | mitigation                                                   | support                  | boundary                          |
|:-------------------------------------------|---------------:|:-------------------------------------------------------------|:-------------------------|:----------------------------------|
| No raw CTC data                            |              5 | Frame current paper as public prior/framework                | PF02, PF11, PF12         | Do not claim transition discovery |
| Public pilot not CTC biology               |              5 | State explicitly that it supports marker/state priors only   | PF05, PF07               | Strong boundary                   |
| Qualitative support scores look subjective |              4 | Label as planning scores and show source layers              | PF10, marker support TSV | Not validation metrics            |
| Spatial relevance challenged               |              3 | Use only as tissue organization/ROI context                  | F05, spatial table       | No CTC origin claim               |
| Mutation-alone reviewer                    |              3 | Show BRAF state split and driver variance boundary           | F03, F04                 | Descriptive, not causal           |
| Too broad / platform-like                  |              4 | Keep paper scoped to PTC, thyroidectomy, CTC-EMT, NGS anchor | PF11                     | Remove drug/AI platform language  |
| Clinical utility overreach                 |              5 | Use research-grade and future-validation language            | PF05                     | No diagnostic claim               |
| Low-input CTC NGS feasibility              |              4 | Stage-gate as CTC-high subset feasibility                    | PF12, CRF                | No universal success claim        |


## Safe Sentence Bank

| claim_area     | safe_sentence                                                                                       | forbidden_sentence                                   |
|:---------------|:----------------------------------------------------------------------------------------------------|:-----------------------------------------------------|
| Mutation/state | BRAF-mutant PTC remains heterogeneous across public tissue-state labels.                            | BRAF status predicts postoperative CTC state.        |
| Spatial        | Spatial coherence supports organized tissue-state context for marker selection.                     | Spatial data prove the origin of CTC shedding.       |
| scRNA          | Single-cell marker context supports a biologically organized panel prior.                           | scRNA validates the circulating CTC phenotype.       |
| Proteomics     | Bulk proteomics proxy supports directional consistency of selected modules.                         | Proteomics proves CTC protein state.                 |
| Marker support | Cross-layer support differs by module and should guide prospective assay priorities.                | Support score is assay validation.                   |
| NGS anchor     | Phenotypic blood signals require matched genetic anchoring before residual-risk interpretation.     | Marker-positive cells are residual disease.          |
| Paper identity | The current manuscript is a public-data prior/framework paper.                                      | The current manuscript discovers CTC-EMT transition. |
| Future data    | Serial hospital CTC/cfDNA/tissue NGS data will be required to test transition and risk association. | Public data already support recurrence prediction.   |


## Next 72h Execution

| window   | task                                                    | input_output                                             | done_definition           |
|:---------|:--------------------------------------------------------|:---------------------------------------------------------|:--------------------------|
| 0-12h    | Freeze title, abstract scaffold, and main figure layout | MANUSCRIPT_BLUEPRINT_KR.md + main_supp_figure_layout.tsv | ready for Prof. Yu        |
| 12-24h   | Send hospital/Yu CRF request                            | hospital_data_request_crf.tsv                            | data request sent         |
| 24-36h   | Convert figure cards into manuscript outline            | figure_caption_and_results_scaffold.tsv                  | Results skeleton ready    |
| 36-48h   | Audit exact dataset citations and methods               | methods_reproducibility_checklist.tsv                    | Methods provenance locked |
| 48-60h   | Prepare Samsung appendix deck pages from PF11/PF12/PF10 | web figures                                              | proposal appendix ready   |
| 60-72h   | Decide preprint/framework route vs wait for Yu raw data | submission_strategy_matrix.tsv                           | go/no-go decision         |
