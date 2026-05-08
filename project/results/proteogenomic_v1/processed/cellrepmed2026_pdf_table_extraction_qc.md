# CRM PDF table extraction QC
- input: `processed/cellrepmed2026_mmc1_supp_info_full.txt`
- out-dir: `processed`
- detected 5 'Table S{N}.' blocks

## Table S1 — Key clinicopathologic data summary of the Discovery Cohort.
- expected: Discovery Cohort 113 patients × clinical + CC subtype
- header (14 cols): ['Patient ID', 'Age', 'Sex', 'ECOG performance status', 'Immune status', 'Status of advanced disease', 'CC subtype', 'Proteomics', 'NGS', 'Digital pathology', 'Radioiodine therapy', 'RAI sensitivity', 'Biochemical response to RAI', 'Histological Type']
- data rows: 121
- warnings (18): ['single-cell line treated as continuation: Eastern Cooperative', 'single-cell line treated as continuation: ECOG', 'single-cell line treated as continuation: Oncology Group']
- wrote: `cellrepmed2026_table_S1_clinical.tsv`

## Table S2 — Next-generation sequencing data summary of tumors in the Discovery Cohort.
- expected: NGS mutation calls per patient
- header (5 cols): ['Patient ID', 'Gene Symbol', 'Mutation Description', 'Mutation Site (reference protein)', 'Variant Allele Frequency']
- data rows: 175
- wrote: `cellrepmed2026_table_S2_ngs.tsv`

## Table S3 — Single-cell RNA-seq and spatial RNA-seq cohort (External Cohort 1).
- expected: External Cohort 1 — scRNA + spatial RNA-seq
- header (15 cols): ['No.', 'Patient Id', 'Age', 'Sex', 'ECOG performance status', 'Immune status', 'scRNA-seq\u200b\u200b', 'spRNA-seq\u200b\u200b', 'CC subtype', 'Tissue', 'TDS', '\u200bSpecific_Mutation_Types', 'Stromal_Score Immune_Score', 'ESTIMATE_Score', 'Tumor_Purity']
- data rows: 15
- wrote: `cellrepmed2026_table_S3_scrna_spatial_cohort.tsv`

## Table S4 — Real-world systemic treatment cohort of advanced DTC (External Cohort 2).
- expected: External Cohort 2 — real-world systemic therapy
- header (3 cols): ['Type of sample acquisition (1.', 'Treatment purpose (1.', 'Condition (1. locally advanced; 2. distant']
- data rows: 120
- wrote: `cellrepmed2026_table_S4_realworld_treatment_cohort.tsv`

## Table S5 — Signature gene sets used to calculate TDS, stromal and immune scores by AddModuleScore in spRNA-seq data.
- expected: TDS / stromal / immune signature gene lists
- header (3 cols): ['TDS signature genes', 'Stromal score signature genes', 'Immune score signature genes']
- data rows: 20
- wrote: `cellrepmed2026_table_S5_signature_gene_sets.tsv`

