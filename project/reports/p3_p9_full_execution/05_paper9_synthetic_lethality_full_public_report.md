# Paper 9 Synthetic Lethality Full-Public Report

## 1. Data downloaded/used
- model_metadata: `Model.csv` (DepMap Public 26Q1, 697.5 kB, existing)
- crispr_gene_effect: `CRISPRGeneEffect.csv` (DepMap Public 26Q1, 440.6 MB, existing)
- ccle_expression: `OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv` (DepMap Public 26Q1, 305.0 MB, existing)
- common_essential: `AchillesCommonEssentialControls.csv` (DepMap Public 26Q1, 17.0 kB, existing)
- nonessential_controls: `AchillesNonessentialControls.csv` (DepMap Public 26Q1, 11.5 kB, existing)
- prism_lfc_collapsed: `Repurposing_Public_24Q2_LFC_COLLAPSED.csv` (PRISM Primary Repurposing DepMap Public 24Q2, 150.1 MB, existing)
- prism_treatment_metadata: `Repurposing_Public_24Q2_Treatment_Meta_Data.csv` (PRISM Primary Repurposing DepMap Public 24Q2, 1.5 MB, existing)
- prism_cell_metadata: `Repurposing_Public_24Q2_Cell_Line_Meta_Data.csv` (PRISM Primary Repurposing DepMap Public 24Q2, 157.3 kB, existing)
- gdsc1_auc: `GDSC1AUCMatrix.csv` (Harmonized GDSC 25Q2, 4.3 MB, existing)
- gdsc1_ic50: `GDSC1Log2IC50Matrix.csv` (Harmonized GDSC 25Q2, 1.9 MB, existing)
- gdsc1_conditions: `GDSC1Log2ViabilityCollapsedConditions.csv` (Harmonized GDSC 25Q2, 118.1 kB, existing)
- gdsc2_auc: `GDSC2AUCMatrix.csv` (Harmonized GDSC 25Q2, 3.9 MB, existing)
- gdsc2_ic50: `GDSC2Log2IC50Matrix.csv` (Harmonized GDSC 25Q2, 1.3 MB, existing)
- gdsc2_conditions: `GDSC2Log2ViabilityCollapsedConditions.csv` (Harmonized GDSC 25Q2, 155.7 kB, existing)
- ctrp_auc: `CTRPAUCMatrix.csv` (Harmonized CTD^2 25Q2, 6.1 MB, existing)
- ctrp_ic50: `CTRPLog2IC50Matrix.csv` (Harmonized CTD^2 25Q2, 3.6 MB, existing)
- ctrp_conditions: `CTRPLog2ViabilityCollapsedConditions.csv` (Harmonized CTD^2 25Q2, 432.3 kB, existing)

## 2. Model/cell-line coverage
- DepMap expression-scored models: 1719
- Thyroid-flagged models: 25

## 3. Lineage-state scoring result
- RAI/thyroid differentiation and TF-collapse scores were computed from public processed CCLE/DepMap expression.
- DM1-like/lineage-silenced score is the inverse of RAI_8 within the DepMap expression context.

## 4. Dependency association result
- Candidate targets tested: 26
- Nominal pan-cancer dependency associations with p < 0.05: 5
- Negative dependency correlations indicate stronger dependency as lineage silencing increases.

## 5. Drug-response result
- Drug association rows: 140
- Drug metrics are direction-normalized where lower viability/AUC/LFC indicates greater sensitivity.

## 6. Top candidate targets
- GLS: MEDIUM (dependency_plus_expression_or_thyroid), dependency association=-0.0881
- CHEK2: LOW (drug_only_nominal), dependency association=0.153
- MCL1: LOW (drug_only_nominal), dependency association=0.126
- BCL2L1: LOW (artifact_flagged), dependency association=0.11
- DNMT1: LOW (artifact_flagged), dependency association=-0.0645
- HDAC2: LOW (drug_only_nominal), dependency association=0.0489
- STAT3: LOW (weak_or_negative), dependency association=0.0487
- IL6R: LOW (weak_or_negative), dependency association=-0.0443
- JAK1: LOW (drug_only_nominal), dependency association=0.0391
- HDAC1: LOW (drug_only_nominal), dependency association=0.0376

## 7. What failed
- Thyroid-only dependency associations are underpowered if few thyroid models overlap the CRISPR matrix.
- Drug metadata matching is keyword-based and therefore candidate-generating, not confirmatory.
- Protected patient-level data and raw sequencing were not used.

## 8. Claim boundary
- Allowed: candidate dependency nomination and research-use vulnerability prioritization.
- Forbidden: clinical treatment selection, proven drug response, or patient treatment recommendation.

## 9. Wet-lab validation priority
- Prioritize targets with negative lineage-dependency association, non-high artifact risk, expression support, and drug/dependency concordance.
- Treat common-essential or proliferation-linked candidates as artifact-risk until validated with orthogonal assays.
