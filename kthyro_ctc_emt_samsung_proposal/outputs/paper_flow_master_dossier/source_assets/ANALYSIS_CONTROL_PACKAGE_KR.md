# Analysis Control Package


## Key Number Ledger

| metric                   | value   | unit   | source_table                                | why_it_matters                                   | boundary                             |
|:-------------------------|:--------|:-------|:--------------------------------------------|:-------------------------------------------------|:-------------------------------------|
| TCGA BRAF-mutant subset  | 274     | tumors | tcga_braf_only_label_distribution.tsv       | BRAF-mutant tumors used for mutation/state split | tissue only; no CTC                  |
| BRAF tissue-state labels | 7       | labels | tcga_braf_only_label_distribution.tsv       | BRAF tumors distribute across multiple states    | label scheme is public-pilot derived |
| Largest BRAF state       | 33.2    | %      | tcga_braf_only_label_distribution.tsv       | largest state is only about one-third            | do not infer recurrence              |
| BRAF ambiguity           | 66.8    | %      | tcga_braf_only_label_distribution.tsv       | remaining BRAF cases are outside largest state   | descriptive only                     |
| Spatial slides           | 16      | slides | spatial_coherence_statistics.tsv            | spatial coherence evaluated across slides        | not CTC origin proof                 |
| Spatial spots            | 57144   | spots  | spatial_spot_vulnerability_scores.tsv       | spot-scale tissue-state context                  | not cell-level pathology             |
| Spatial z>2 slides       | 16/16   | slides | spatial_coherence_statistics.tsv            | same-niche coherence supports organization       | not shedding evidence                |
| scRNA module rows        | 12      | rows   | scrna_module_celltype_attribution.tsv       | marker-context prior                             | not blood validation                 |
| Proteomics contrast rows | 40      | rows   | bulk_proteomics_kthyro_module_contrasts.tsv | protein-direction support                        | bulk proxy only                      |
| Driver-only median R2    | 17.2    | %      | tcga_axis_variance_models.tsv               | driver-only model boundary                       | descriptive not causal               |


## Statistical Analysis Plan

| analysis_stage         | question                                                   | input_data                                     | planned_method                                             | output                      | allowed_interpretation                   | boundary                                  |
|:-----------------------|:-----------------------------------------------------------|:-----------------------------------------------|:-----------------------------------------------------------|:----------------------------|:-----------------------------------------|:------------------------------------------|
| Current in silico RQ1  | Does driver mutation collapse PTC into one state?          | TCGA BRAF label distribution                   | fraction by state label; descriptive distribution          | barplot + state split table | Mutation alone is insufficient           | No CTC or outcome claim                   |
| Current in silico RQ2  | Are tissue-state labels spatially organized?               | GSE250521 spatial scores                       | same-neighbor coherence z-score                            | coherence barplot           | Tissue context supports marker/ROI logic | No CTC origin proof                       |
| Current in silico RQ3  | Which marker modules have public support?                  | TCGA/spatial/scRNA/proteomics/Yu anchor coding | 0/1/2 cross-layer support score                            | marker support heatmap      | Prioritizes marker modules               | Not assay validation                      |
| Current in silico RQ4  | What claims are currently writeable?                       | claim/evidence scorecard                       | evidence grade + support figure mapping                    | claim scorecard             | Framework paper is feasible              | No discovery claim                        |
| Future prospective RQ1 | Does thyroidectomy perturb CTC E/EM/M state?               | serial CTC phenotype table                     | paired change; mixed-effects or paired nonparametric model | transition plot             | Primary biological endpoint              | Requires raw serial data                  |
| Future prospective RQ2 | Are persistent postoperative signals genetically anchored? | tumor-normal NGS + cfDNA/CTC NGS               | variant concordance/VAF kinetics/QC pass                   | clone map                   | Genetic map of persistent signals        | No universal success assumption           |
| Future prospective RQ3 | Do EM/M persistent signals associate with risk features?   | CTC/cfDNA + pathology + Tg/US/follow-up        | association model; exploratory FDR control                 | risk association panel      | Hypothesis-generating risk biology       | No clinical prediction without validation |


## Stage-Gate Kill Criteria

| gate   | stage                     | pass_rule                                     | kill_or_downgrade_rule                        | action_if_failed                   | status         |
|:-------|:--------------------------|:----------------------------------------------|:----------------------------------------------|:-----------------------------------|:---------------|
| G1     | Public-data manuscript    | all source tables present and linked          | missing or undocumented input                 | do not write unsupported Methods   | currently pass |
| G2     | CTC assay reproducibility | repeatable CTC detection/subtyping            | CTC detection/subtyping not reproducible      | stop transition claim              | future         |
| G3     | Tissue-normal NGS         | success rate >=85%                            | success <85%                                  | clone-anchor claim weakened        | future         |
| G4     | cfDNA feasibility         | usable signal in high-risk or CTC-high subset | unusable cfDNA even in enriched subset        | cfDNA becomes exploratory/optional | future         |
| G5     | CTC-enriched NGS          | QC pass in CTC-high subset                    | low-input CTC NGS fails QC                    | do not claim CTC genotype map      | future         |
| G6     | Follow-up capture         | >=80% key follow-up fields                    | follow-up capture <80%                        | no outcome/risk association        | future         |
| G7     | Claim discipline          | all claims map to evidence table              | unsupported clinical/prediction claim appears | remove sentence/figure             | currently pass |


## Reproducibility Manifest

| artifact                           | purpose                                     | location                  | use_in_manuscript     | qc_action                          |
|:-----------------------------------|:--------------------------------------------|:--------------------------|:----------------------|:-----------------------------------|
| build_no_yu_in_silico_paper.py     | generates no-Yu figures/tables/reports      | in_silico_paper_no_yu     | source analysis layer | rerun before manuscript freeze     |
| build_no_yu_cv2_storyboard.py      | generates CV2 storyboards                   | visual_qc                 | visual overview       | optional rerun                     |
| build_no_yu_in_silico_web.py       | deploys no-Yu dossier                       | papers hub                | web evidence page     | deployed                           |
| build_paper_flow_master_dossier.py | generates master manuscript control dossier | paper_flow_master_dossier | current page          | rerun after edits                  |
| key_number_ledger.tsv              | locks numbers used in text                  | tables                    | fact checking         | must match final manuscript        |
| statistical_analysis_plan.tsv      | separates current and future analyses       | tables                    | reviewer defense      | must not overclaim future analyses |
| hospital_data_request_crf.tsv      | defines hospital/Yu fields                  | tables                    | data request          | send to collaborator               |


## Scenario Decision Tree

| scenario                            | paper_route                                         | what_to_do                                         | boundary                                        |
|:------------------------------------|:----------------------------------------------------|:---------------------------------------------------|:------------------------------------------------|
| If Yu raw data arrives              | upgrade to serial CTC transition manuscript         | add patient-level transition and clone map figures | still avoid clinical prediction until follow-up |
| If only aggregate Yu data arrives   | keep current prior/framework paper                  | cite as feasibility anchor only                    | no reanalysis claim                             |
| If no hospital data before deadline | submit/prepare framework paper and Samsung appendix | use PF11/PF13 execution route                      | make limitations explicit                       |
| If CTC NGS fails                    | retain phenotype-transition biology only            | downgrade genetic-map Aim to cfDNA/tissue anchor   | do not hide failure                             |
| If follow-up is immature            | publish biological transition/gap paper             | defer recurrence/risk association                  | no prediction language                          |
