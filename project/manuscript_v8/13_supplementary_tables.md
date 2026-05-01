---
title: "Paper 1 manuscript v8 — Supplementary Tables S1-S10 (index + source TSVs)"
date: 2026-05-01
status: v1 — existing TSV files linked. Final formatting (Excel multi-sheet workbook) W4-W5.
---

# Supplementary Tables Index

Paper 1 manuscript v8 supplementary tables, with source TSV file paths in `project/results/`. Final submission packaging will combine these into a single multi-sheet Excel workbook (`THCA_paper1_supplementary_tables.xlsx`).

---

## S1 — TIERA67 67-gene panel definition (7 categories)

**Description.** Curated candidate gene pool spanning thyroid biological categories used for 8-gene panel selection.

**Source.** `project/metadata/tierA67_genes.txt`

**Format (10 columns):**

| Column | Type | Description |
|---|---|---|
| `gene_symbol` | str | HUGO symbol |
| `category` | str | One of 7 categories: TDS_core / MAPK_output_ERK / Driver_anchor / Aggressive_marker / Dediff_invasion / Immune_stromal_light / Thyroid_lineage_extra |
| `n_genes_in_category` | int | Count |
| `Yoo2016_inclusion` | bool | In Yoo 2016 BRS 273-gene |
| `TCGA2014_BRS` | bool | In TCGA 2014 BRS |
| `8gene_panel` | bool | In final 8-gene panel |
| `entrez_id` | int | NCBI Entrez |
| `ensembl_id` | str | GENCODE v44 |
| `chromosome` | str | hg38 |
| `notes` | str | biological role |

n=67 rows.

---

## S2 — Cohort overview

**Description.** Multi-cohort sample-level descriptor table.

**Source.** Composed from:
- `project/results/tables/tcga_thca_clinical_extended.tsv` (TCGA n=504)
- `project/results/audit_2026_04_30/round3/cbio_sv_thca.tsv` (TCGA SV)
- `project/results/audit_2026_04_30/round5/r5_2_sample_methylation_8gene.tsv` (TCGA HM450)
- `project/results/v17_korean/` (K2 + Lee + GSE286332)

**Format (~12 columns):**

| Column | Description |
|---|---|
| `cohort` | TCGA / MSK-IMPACT / K2 / Lee / GSE286332 / GSE184362 / Lu2023 |
| `n_total` | sample count |
| `n_with_OS` | overall survival annotation |
| `n_with_RNAseq` | bulk RNA-seq |
| `n_with_WGS` | WGS / WES |
| `n_with_HM450` | methylation |
| `n_with_SV` | structural variant |
| `assay_platform` | RNA-seq / array / scRNA-seq |
| `tissue_type` | FF / FFPE / mixed |
| `ancestry_majority` | EUR / EA / mixed |
| `event_rate_pct` | OS event rate |
| `cite` | primary publication |

---

## S3 — TCGA-THCA per-sample 8-gene scores + DM call + clinical

**Description.** Per-sample 8-gene panel score, P(DM1) classification probability, DM cluster call, and clinical metadata.

**Source.** `project/results/tables/tcga_thca_8gene_per_sample.tsv` (compose from existing v17 outputs)

**Format (~25 columns):**

| Column | Description |
|---|---|
| `sample_id` | TCGA-XX-XXXX-01 (primary tumor) |
| `SLC5A5_log2tpm` ... `DIO1_log2tpm` | 8 panel gene log2(TPM+1) |
| `panel_mean_z` | within-cohort z-mean |
| `P_DM1` | logistic regression class probability |
| `DM_call` | DM1 / DM2 |
| `BRS_273gene` | Yoo 2016 BRAF-RAS Score |
| `BRAF_V600E` | mutation status |
| `RAS_hotspot` | NRAS/HRAS/KRAS hotspot |
| `TERT_promoter` | promoter mutation |
| `fusion_status` | RET/NTRK/ALK/BRAF/none |
| `age_at_diagnosis` | years |
| `sex` | M/F |
| `stage` | I-IV |
| `OS_status` | dead/alive |
| `OS_days` | follow-up |

n=504 rows.

---

## S4 — TIERA67 univariate Cohen's d ranking (full 67 genes)

**Description.** Per-gene Cohen's d (DM1 vs DM2) ranking across all 67 TIERA67 genes, including drivers.

**Source.** `project/results/p1_driver_mrna_audit/top20_by_d_full_tiera67.tsv` (extend to full 67)

**Format:** gene_symbol, category, mean_DM1, mean_DM2, Cohen_d, MW_p, BH_FDR_q, rank.

Highlights (per memory): 8-gene at ranks 4 (TPO), 5 (DIO1), 18 (TG), 19 (PAX8), 25 (FOXE1), 38 (NKX2-1), 40 (TSHR), 50 (SLC5A5); drivers at ranks 15 (CDKN2B), 22 (RET), 24 (CDKN2A), 52 (BRAF), 56 (TERT), 63 (KRAS), 66 (NRAS), 67 (HRAS).

n=67 rows.

---

## S5 — Pan-genome cluster ARI ladder

**Description.** Cluster reproducibility across panel sizes 8, 16, 67 (TIERA67), 200, 1000, 5000 (pan-genome MAD-ranked).

**Source.** `project/results/p4_pangenome_vs_tiera67/ari_comparison.tsv`

**Format:** panel_definition, n_genes, ARI_vs_8gene, NMI_vs_8gene, hypergeometric_p_in_top100.

Highlights: 8-gene 0.49 / TIERA67 0.90 / Pan top-200 0.86 / top-1000 0.90 / top-5000 0.92 / Driver_anchor 12 only −0.007.

---

## S6 — DM1 fusion details (per-sample fusion class + partner)

**Description.** Per-sample fusion partner annotations within TCGA + MSK SV-tested cohort.

**Source.** Compose from:
- `project/results/audit_2026_04_30/round3/cbio_sv_thca.tsv` (TCGA n=542)
- `project/results/audit_2026_04_30/round5/msk_sv_thca.tsv` (MSK n=12)

**Format:** sample_id, cohort, fusion_class (RET/NTRK/ALK/BRAF/PAX8-PPARG/none), fusion_partner_5p, fusion_partner_3p, exon_breakpoint, DM_call, DM1_subcluster (sub-A/sub-B).

Highlights: TCGA RET fusion CCDC6-RET 17, NCOA4-RET 3, other 13. MSK RET fusion CCDC6-RET 3, NCOA4-RET 2.

---

## S7 — DM1 vs DM2 per-gene HM450 methylation β-values

**Description.** Per-gene HM450 promoter β-value comparison across 8-gene panel + DIO2 + SLC26A4 + DUOX2.

**Source.** `project/results/audit_2026_04_30/round5/r5_2_per_gene_methylation_DM.tsv`

**Format:** gene_symbol, n_DM1, n_DM2, n_notDM, mean_DM1, mean_DM2, mean_notDM, Cohen_d, MW_p, BH_FDR_q, fusion_independence_p_within_DM1.

Highlights: TPO d=2.30 p=1.9e-18; DIO1 d=1.24; TSHR d=1.20; SLC5A5 d=0.22 NS.

---

## S8 — Meta-analysis raw inputs (TCGA + MSK Cox + DerSimonian-Laird)

**Description.** Per-cohort hazard ratio inputs and pooled meta-analysis output.

**Source.** Compose from:
- `project/results/audit_2026_04_30/round2/n1_meta.json`
- `project/results/p6_multi_cohort_meta.json`

**Format:** cohort, n_total, n_events, log_HR_DM1_vs_DM2, SE_logHR, HR, HR_lower_95, HR_upper_95, p_value, weight_in_meta.

Highlights: TCGA HR 2.30 [0.77, 6.88], n=504, 16 events. MSK HR 2.67 [1.17, 6.10], n=117, 38 events. Pooled HR 2.53 [1.31, 4.89], I²=0%.

---

## S9 — Korean Pan-Asian arcasHLA per-sample genotypes

**Description.** Per-sample 4-digit HLA allele calls for Korean PTC pool (n=874) and GSE286332 (n=18).

**Source.** Compose from:
- `project/results/d4p1_panasian_meta/korean_PTC_pool_n908.tsv`
- `project/results/v17_korean/arcasHLA_GSE286332/{18 × genotype.json}`

**Format:** sample_id, cohort, HLA_A_1, HLA_A_2, HLA_B_1, HLA_B_2, HLA_C_1, HLA_C_2, HLA_DRB1_1, HLA_DRB1_2, HLA_DPB1_1, HLA_DPB1_2, HLA_DQB1_1, HLA_DQB1_2.

★ **Paper 2 backbone** — Paper 1 supp only.

---

## S10 — Reviewer Q&A pre-empt 12 items

**Description.** 12 anticipated reviewer questions with prepared responses (revision-ready text).

**Source.** `project/manuscript_v8/09_reviewer_qa.md`

---

# Final packaging plan (W5)

1. Execute compose scripts to generate each TSV from existing v17 outputs (where missing).
2. Convert each TSV to Excel sheet using `pandas.to_excel`.
3. Combine into single multi-sheet workbook `THCA_paper1_supplementary_tables.xlsx`.
4. Write Supplementary Table caption `.docx` per Cell Press submission template.

Submission package:
- `manuscript.docx` (main text + figures inline + references)
- `figures.zip` (Fig 1-8 PDF/PNG + S1-S9 PDF/PNG)
- `THCA_paper1_supplementary_tables.xlsx`
- `cover_letter.docx`
- `reviewer_qa.docx` (S10 standalone)
