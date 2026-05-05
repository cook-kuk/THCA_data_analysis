# Paper 9 GLS Sprint 1 Report

## 1. Executive verdict
- GLS status: KEEP_AS_CANDIDATE (MEDIUM_MINUS).
- Interpretation: GLS/glutamine axis remains a candidate vulnerability hypothesis, but Sprint 1 does not prove synthetic lethality.
- Main limitation: thyroid-only CRISPR/drug overlap is small, so pan-cancer artifact correction carries most of the evidence weight.

## 2. Thyroid model coverage
- thyroid_flagged_models: 25
- strict_thyroid_lineage_models: 22
- crispr_overlap: 14
- strict_thyroid_crispr_overlap: 11
- ccle_expression_overlap: 25
- strict_thyroid_expression_overlap: 22
- prism_overlap: 16
- gdsc_ctrp_overlap: 15
- gls_expression_dependency_drug_complete: 9
- strict_gls_expression_dependency_drug_complete: 9

## 3. GLS dependency result
- GLS dependency vs lineage-silenced score: Spearman r=-0.0881, p=0.0029, n=1140.
- Negative direction means higher lineage-silenced score is associated with more negative GLS gene effect.

## 4. Artifact-control result
- GLS dependency vs common-essential proxy: Spearman r=0.0972, p=0.00102, n=1140.
- Control-gene distributions and nonessential controls are reported in `gls_control_gene_correlations.tsv` and `fig_gls_artifact_controls.png`.

## 5. Glutamine pathway result
- Glutamine dependency module vs lineage-silenced score: r=0.125, p=2.29e-05.
- OXPHOS/glycolysis comparator modules are included in `glutamine_pathway_module_associations.tsv`.

## 6. Drug concordance result
- GDSC1_AUC BPTES GDSC1:1425 DPC-007404: r=-0.0153, p=0.694, concordance=discordant_or_no_support.

## 7. Updated target ranking
- GLS: KEEP_AS_CANDIDATE (MEDIUM_MINUS); artifact risk=low_to_moderate.
- glutamine_pathway_module: KEEP_AS_AXIS (HYPOTHESIS_SUPPORT); artifact risk=moderate; pathway overlaps MYC/LDHA metabolic proliferation biology.
- CHEK2: UNCHANGED_LOW (drug_only_nominal); artifact risk=moderate_or_low.
- MCL1: UNCHANGED_LOW (drug_only_nominal); artifact risk=moderate_or_low.
- DNMT1: UNCHANGED_LOW (artifact_flagged); artifact risk=HIGH_common_essential_or_panlethal.
- STAT3: UNCHANGED_LOW (weak_or_negative); artifact risk=moderate_or_low.
- JAK1: UNCHANGED_LOW (drug_only_nominal); artifact risk=moderate_or_low.
- HDAC1: UNCHANGED_LOW (drug_only_nominal); artifact risk=moderate_or_low.
- HDAC2: UNCHANGED_LOW (drug_only_nominal); artifact risk=moderate_or_low.

## 8. Claim boundary
Allowed:
- GLS/glutamine axis is a candidate vulnerability.
- Research-use hypothesis.

Forbidden:
- Validated synthetic lethality.
- Clinical treatment recommendation.
- Patient selection.

## 9. Wet-lab validation plan
- GLS inhibitor sensitivity in DM1-like thyroid models.
- Rescue with glutamine pathway modulation.
- Organoid/cell-line validation with lineage-state stratification.
