---
dir: results/p3_gse286332
pillar: Paper 2 Pillar II — In-cohort mechanism (GSE286332 PTC vs PTC+HT)
status: STRONG
generator: notebooks_or_scripts/v17_P3_GSE286332_ptc_vs_ptcht.py
---

# p3_gse286332 — Paper 2 Pillar II source data

GSE286332 (Korean Dongguk Univ Lim 2025, n=18 = 9 PTC + 9 PTC+HT) bulk RNA-seq DEG + GSEA + 8-gene RAI score + HLA module score.

## Files

| File | Description | Used by |
|---|---|---|
| `P3_summary.json` | Top-level structured results: DEG counts, panel_8gene_compare, hla_per_gene_compare, deg_top10_up | Brief Fig 3, 5, 6; full structured anchor |
| `8gene_panel_per_sample.tsv` | Per-sample 8-gene RAI score (n=18) | Pillar II / Pillar III mediation input |
| `8gene_per_gene_compare.tsv` | Per-gene PTC+HT vs PTC Cohen's d + MW p | Brief **Fig 3** (8-gene Cohen d bar) |
| `deg_ptcht_vs_ptc.tsv` | PyDESeq2 DEGs (10,380 padj<0.05) | Brief **Fig 5** (top 10 up DEG) + GSEA input |
| `dm12_predictions.tsv` | TCGA-trained classifier P_DM1 / P_DM2 per sample | Pillar III mediation, § 3.6 |
| `gsea_MSigDB_Hallmark_2020.tsv` | Hallmark gene set NES + FDR | Brief **Fig 4** |
| `gsea_KEGG_2021_Human.tsv` | KEGG pathway NES + FDR | Brief **Fig 4** (Type I diabetes NES=+1.92) |
| `gsea_Reactome_2022.tsv` | Reactome pathway NES + FDR | Brief **Fig 4** (TCR/CD3/ZAP-70 cascades) |
| `hla_module_scores.tsv` | HLA-I + HLA-II module mean expr per sample | Brief **Fig 6 ★★★** (HLA-II d=+3.65 paper-defining) |

## Brief usage cross-reference

| Brief Fig | Source file (relative to this dir) |
|---|---|
| Fig 3 — 8-gene RAI Cohen's d | `8gene_per_gene_compare.tsv` + `P3_summary.json:panel_8gene_compare` |
| Fig 4 — GSEA top pathways NES | `gsea_MSigDB_Hallmark_2020.tsv` + `gsea_KEGG_2021_Human.tsv` + `gsea_Reactome_2022.tsv` |
| Fig 5 — Top 10 up DEGs | `deg_ptcht_vs_ptc.tsv` + `P3_summary.json:deg_top10_up` |
| **Fig 6 ★★★** — HLA-I/II module Cohen d | `hla_module_scores.tsv` |

## Key claims

- 10,380 DEGs (PyDESeq2 BH-FDR < 0.05; up 6,004 / down 4,376)
- 8-gene RAI panel Cohen's d = **−1.602** (MW p=4.1e−4); per-gene TF backbone (PAX8/NKX2-1/FOXE1) d=−1.7 ~ −2.3, SLC5A5/NIS d=+0.25 (NS preserved)
- HLA-II module Cohen's d = **+3.65** (★★★ paper-defining magnitude; d=+2.34 for HLA-I)
- GSEA: Hallmark Allograft Rejection NES=+2.12 (FDR=0); KEGG Type I Diabetes NES=+1.92; Reactome TCR/CD3 cascades all FDR≤2e−4

## Cohort

- BioProject PRJNA1208932 (Dongguk Univ, Macrogen Seoul, Illumina NovaSeq X)
- n=18 (9 PTC + 9 PTC+HT), Korean
- arcasHLA outputs in `../d4p1_panasian_meta/GSE286332_arcasHLA_genotypes.tsv`
