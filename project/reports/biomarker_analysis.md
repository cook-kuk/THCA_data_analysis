# THCA Biomarker Analysis: Known vs Novel

## TL;DR

Tested 51,711 genes on TCGA-THCA BRAF_like (n=392) vs RAS_like (n=111). 148 of the 221 hardcoded + panel reference genes were present; 70 replicated (same direction, external p<0.05) in ≥1 external cohort. 2773 genes passed the novel-validated filter (FDR<0.05, |d|>0.5, replication≥1, median log2>1, present in ≥2 cohorts). Top novel candidates: **TACSTD2** (d=+2.52, replicated 2/2), **TMPRSS4** (d=+2.09, replicated 2/2), **PLEKHA6** (d=+2.22, replicated 2/2). Genes with large TCGA effects that fail to replicate in BOTH externals are **not** called biomarkers.

## Methods summary

- **Discovery**: Welch t-test on log2 expression between BRAF_like and RAS_like TCGA samples, BH-FDR across all tested genes. log2FC is defined as mean(BRAF) − mean(RAS).
- **Replication**: same test re-run on GSE27155 microarray and GSE126698 RNA-seq. A gene replicates in a cohort if its external log2FC has the same sign as TCGA and the external p-value is <0.05 (per-cohort, uncorrected).
- **External labels**: GSE27155 uses `molecular_subtype` from `sample_master.tsv` (cPTC/tall-cell→BRAF_like, FVPTC/FTC→RAS_like as already encoded by rerun_v2). GSE126698 sample codes are mapped by prefix: P*→BRAF_like (PTC), F*→RAS_like (FTC); A* (ATC) and N* (normal) are excluded.
- **Known set**: hardcoded reference table (below) plus genes in `tds16_genes.txt`, `tierA67_genes.txt`, `brs71_genes.txt`.
- **Novel filter**: NOT in any known list AND TCGA FDR<0.05 AND |d|>0.5 AND replication_rate≥1 AND median TCGA log2>1 AND present in ≥2 cohorts.
- **Novelty score** = |d|_TCGA × replication_rate × −log10(FDR) (used for ranking, not as a hard threshold).

## Known biomarkers: what replicated, what didn't

| gene | category | log2FC_TCGA | d | FDR_TCGA | repl_GSE27155 | repl_GSE126698 | rep_rate |
|------|----------|-------------|---|----------|-----|-----|---|
| SYT12 | - | +4.88 | +2.47 | 1.07e-66 | 0 | 1 | 1 |
| FN1 | - | +4.78 | +2.45 | 4.85e-58 | 1 | 0 | 1 |
| SERPINA1 | BRS_literature | +4.19 | +2.38 | 1.68e-37 | 1 | 1 | 2 |
| KCNN4 | - | +3.78 | +2.37 | 4.47e-59 | 1 | 1 | 2 |
| DCSTAMP | - | +5.35 | +2.24 | 1.75e-82 | 1 | 1 | 2 |
| ERBB3 | BRS_literature | +2.29 | +2.20 | 2.42e-38 | 1 | 0 | 1 |
| TMPRSS6 | - | +4.31 | +2.14 | 1.56e-74 | 1 | 0 | 1 |
| KRT19 | Clinical_IHC | +3.29 | +2.07 | 1.47e-31 | 1 | 1 | 2 |
| STAC | - | +2.13 | +2.07 | 2.68e-64 | 1 | 0 | 1 |
| LGALS3 | Clinical_IHC | +2.45 | +2.04 | 2.37e-28 | 1 | 1 | 2 |
| CREB5 | - | +2.09 | +2.03 | 3.59e-49 | 1 | 1 | 2 |
| MGAT3 | - | +2.84 | +1.95 | 1.09e-47 | 0 | 1 | 1 |
| FAM111A-DT | - | +1.27 | +1.93 | 1.77e-48 | 0 | 0 | 0 |
| MET | Aggressiveness_EMT | +2.04 | +1.92 | 1.47e-36 | 1 | 1 | 2 |
| AP002358.1 | - | +2.34 | +1.92 | 1.01e-93 | 0 | 0 | 0 |
| KLK7 | - | +3.94 | +1.91 | 1.74e-63 | 1 | 0 | 1 |
| WNT10A | - | +2.47 | +1.87 | 5.16e-65 | 0 | 0 | 0 |
| FCHO1 | - | +1.68 | +1.86 | 2.34e-47 | 0 | 0 | 0 |
| CLDN10 | - | +3.16 | +1.84 | 2.30e-67 | 1 | 0 | 1 |
| TIMP1 | BRS_literature | +2.53 | +1.83 | 4.95e-28 | 1 | 1 | 2 |
| LRRC52-AS1 | - | +2.82 | +1.82 | 3.26e-69 | 0 | 0 | 0 |
| KLK10 | - | +3.50 | +1.81 | 1.40e-37 | 1 | 0 | 1 |
| SPOCK2 | - | +2.75 | +1.79 | 9.42e-49 | 1 | 1 | 2 |
| AC009549.1 | - | +1.27 | +1.78 | 6.84e-66 | 0 | 0 | 0 |
| TMEM92 | - | +2.12 | +1.77 | 4.61e-78 | 0 | 0 | 0 |
| DRAXIN | - | +1.37 | +1.73 | 1.64e-41 | 0 | 0 | 0 |
| SLC25A47P1 | - | +1.58 | +1.71 | 1.06e-86 | 0 | 0 | 0 |
| ST6GALNAC5 | - | +2.91 | +1.71 | 1.49e-61 | 1 | 1 | 2 |
| INAVA | - | +1.89 | +1.71 | 3.77e-82 | 0 | 0 | 0 |
| DUSP5 | MAPK_output | +2.25 | +1.68 | 2.81e-37 | 1 | 0 | 1 |
| SLC5A8 | Thyroid_lineage_TDS | -3.06 | -1.68 | 8.53e-35 | 0 | 0 | 0 |
| WARS1P1 | - | +1.24 | +1.66 | 1.75e-84 | 0 | 0 | 0 |
| BNC1 | - | +1.65 | +1.63 | 4.99e-66 | 1 | 1 | 2 |
| CYP1B1-AS1 | - | +0.90 | +1.62 | 8.99e-46 | 0 | 0 | 0 |
| DIO1 | Thyroid_lineage_TDS | -3.89 | -1.62 | 8.74e-36 | 1 | 1 | 2 |
| TPO | Thyroid_lineage_TDS | -4.14 | -1.61 | 5.51e-40 | 1 | 1 | 2 |
| RPS29P11 | - | +1.25 | +1.59 | 7.29e-80 | 0 | 0 | 0 |
| LINC00607 | - | +0.99 | +1.55 | 5.88e-38 | 0 | 0 | 0 |
| PNPLA5 | - | +1.53 | +1.54 | 2.27e-82 | 0 | 0 | 0 |
| DUSP6 | MAPK_output | +1.59 | +1.52 | 1.08e-21 | 1 | 1 | 2 |
| IGFL2 | - | +2.22 | +1.50 | 1.49e-61 | 0 | 1 | 1 |
| CEACAM6 | - | +2.60 | +1.50 | 5.58e-73 | 1 | 0 | 1 |
| ADAMTS14 | - | +1.76 | +1.49 | 5.92e-60 | 0 | 0 | 0 |
| AL137026.1 | - | +1.18 | +1.46 | 4.31e-49 | 0 | 0 | 0 |
| CRLF2 | - | +1.39 | +1.44 | 1.93e-74 | 1 | 0 | 1 |
| LY6G6C | - | +1.24 | +1.44 | 2.10e-61 | 0 | 0 | 0 |
| LOX | Aggressiveness_EMT | +1.91 | +1.43 | 6.33e-41 | 1 | 0 | 1 |
| ACTBL2 | - | +0.95 | +1.43 | 4.99e-71 | 0 | 0 | 0 |
| AC002401.4 | - | +1.07 | +1.42 | 4.92e-63 | 0 | 0 | 0 |
| SYT1 | - | +2.06 | +1.41 | 4.95e-51 | 1 | 1 | 2 |
| AL096865.1 | - | +0.93 | +1.39 | 3.79e-46 | 0 | 0 | 0 |
| TM4SF4 | - | +1.81 | +1.39 | 5.02e-42 | 1 | 0 | 1 |
| DSC3 | - | +2.10 | +1.38 | 3.15e-51 | 1 | 0 | 1 |
| RNF183 | - | +1.27 | +1.36 | 6.47e-52 | 0 | 0 | 0 |
| BRINP2 | - | +1.25 | +1.36 | 2.66e-47 | 0 | 0 | 0 |
| CDK5RAP2 | - | +1.17 | +1.35 | 4.38e-42 | 1 | 0 | 1 |
| LINC02408 | - | +1.48 | +1.33 | 3.68e-48 | 0 | 0 | 0 |
| DIO2 | Thyroid_lineage_TDS | -1.65 | -1.32 | 2.58e-23 | 1 | 1 | 2 |
| FAM155B | - | -2.14 | -1.32 | 2.74e-35 | 1 | 0 | 1 |
| CFB | - | +0.96 | +1.31 | 2.04e-37 | 1 | 0 | 1 |
| IVL | - | +1.92 | +1.30 | 1.63e-53 | 1 | 0 | 1 |
| HLA-G | - | +1.91 | +1.29 | 4.14e-43 | 1 | 0 | 1 |
| BEND6 | - | +1.12 | +1.28 | 1.48e-48 | 0 | 0 | 0 |
| AC004847.1 | - | +1.77 | +1.27 | 1.38e-39 | 0 | 0 | 0 |
| CRYBG2 | - | +0.93 | +1.24 | 6.73e-47 | 0 | 0 | 0 |
| VGLL1 | - | +1.14 | +1.24 | 2.52e-54 | 0 | 0 | 0 |
| SIGLEC6 | - | +2.26 | +1.22 | 4.82e-61 | 1 | 0 | 1 |
| HLA-DRA | - | +1.68 | +1.21 | 4.79e-18 | 1 | 0 | 1 |
| SDR16C5 | - | +1.07 | +1.19 | 8.56e-60 | 0 | 0 | 0 |
| TMPRSS11E | - | +1.47 | +1.19 | 6.86e-53 | 1 | 0 | 1 |
| SPOCD1 | - | +1.05 | +1.18 | 2.21e-51 | 0 | 0 | 0 |
| CDKN2B | Aggressiveness_EMT | +1.07 | +1.18 | 4.70e-31 | 0 | 0 | 0 |
| SLC26A4 | Thyroid_lineage_TDS | -2.34 | -1.18 | 3.13e-26 | 0 | 0 | 0 |
| DMBX1 | - | +1.60 | +1.18 | 3.97e-52 | 0 | 0 | 0 |
| RAB27B | - | +1.76 | +1.17 | 8.34e-47 | 0 | 0 | 0 |
| FGFBP1 | - | +1.80 | +1.17 | 3.30e-43 | 0 | 1 | 1 |
| HMGA2 | BRS_literature | +1.96 | +1.16 | 1.26e-13 | 1 | 1 | 2 |
| ELFN2 | - | +1.61 | +1.16 | 3.44e-41 | 0 | 0 | 0 |
| ANGPTL4 | BRS_literature | +1.33 | +1.16 | 6.18e-17 | 0 | 1 | 1 |
| PODNL1 | - | +1.20 | +1.16 | 2.08e-48 | 0 | 0 | 0 |
| ARSI | - | +0.94 | +1.14 | 8.40e-49 | 0 | 0 | 0 |
| FOXP3 | - | +0.93 | +1.13 | 1.75e-36 | 0 | 0 | 0 |
| VTCN1 | - | +1.77 | +1.11 | 2.01e-51 | 1 | 0 | 1 |
| DUSP4 | MAPK_output | +1.09 | +1.11 | 7.33e-13 | 1 | 1 | 2 |
| MYBPH | - | +1.10 | +1.09 | 2.27e-35 | 0 | 1 | 1 |
| CST5 | - | +0.85 | +1.07 | 6.88e-48 | 1 | 0 | 1 |
| KLK6 | - | +1.28 | +1.07 | 3.11e-48 | 1 | 0 | 1 |
| TG | Thyroid_lineage_TDS | -1.48 | -1.07 | 1.79e-21 | 1 | 0 | 1 |
| DUOX2 | Thyroid_lineage_TDS | -1.50 | -1.06 | 7.42e-17 | 1 | 0 | 1 |
| SLC6A20 | - | +1.82 | +1.05 | 4.29e-27 | 1 | 0 | 1 |
| CD274 | - | +0.94 | +1.04 | 1.76e-13 | 0 | 0 | 0 |
| EREG | - | +1.14 | +1.01 | 1.17e-42 | 1 | 0 | 1 |
| SLC6A14 | - | +1.32 | +1.00 | 8.73e-36 | 1 | 0 | 1 |
| IYD | Thyroid_lineage_TDS | -1.61 | -0.94 | 3.37e-16 | 0 | 1 | 1 |
| MMP9 | Aggressiveness_EMT | +1.51 | +0.92 | 7.23e-14 | 0 | 0 | 0 |
| PAX8 | Thyroid_lineage_TDS | -0.58 | -0.91 | 2.10e-10 | 1 | 0 | 1 |
| CST2 | - | +1.26 | +0.91 | 1.18e-20 | 1 | 0 | 1 |
| FOSL1 | MAPK_output | +1.05 | +0.88 | 2.51e-16 | 0 | 0 | 0 |
| PIK3CA | Aggressiveness_EMT | +0.48 | +0.87 | 2.81e-12 | 1 | 0 | 1 |
| DUOX1 | Thyroid_lineage_TDS | -0.76 | -0.85 | 3.17e-10 | 1 | 0 | 1 |
| CDKN2A | Aggressiveness_EMT | +0.91 | +0.84 | 2.23e-18 | 0 | 0 | 0 |
| CTNNB1 | - | +0.43 | +0.77 | 4.43e-08 | 1 | 1 | 2 |
| LRP4 | BRS_literature | +1.19 | +0.76 | 4.49e-07 | 1 | 1 | 2 |
| ETV4 | MAPK_output | +0.84 | +0.75 | 4.61e-07 | 0 | 0 | 0 |
| SPRY4 | MAPK_output | -0.52 | -0.72 | 5.49e-10 | 0 | 0 | 0 |
| RET | Driver_mutation_fusion_proxy | +1.07 | +0.71 | 1.67e-12 | 1 | 1 | 2 |
| FOXE1 | Thyroid_lineage_TDS | -0.56 | -0.70 | 7.03e-07 | 1 | 0 | 1 |
| CITED1 | Clinical_IHC | +0.78 | +0.67 | 2.52e-08 | 1 | 1 | 2 |
| MSH2 | - | +0.29 | +0.67 | 6.63e-07 | 1 | 0 | 1 |
| MSLN | Clinical_IHC | +0.76 | +0.67 | 6.09e-16 | 0 | 0 | 0 |
| CDH1 | Aggressiveness_EMT | +0.32 | +0.58 | 3.69e-05 | 0 | 0 | 0 |
| TWIST1 | Aggressiveness_EMT | +0.38 | +0.49 | 8.92e-07 | 1 | 0 | 1 |
| ALK | Driver_mutation_fusion_proxy | +0.53 | +0.49 | 1.87e-05 | 0 | 1 | 1 |
| SNAI2 | Aggressiveness_EMT | +0.40 | +0.48 | 3.07e-05 | 0 | 1 | 1 |
| CEACAM5 | Clinical_IHC | +0.23 | +0.45 | 2.18e-11 | 0 | 0 | 0 |
| IDO1 | - | +0.60 | +0.45 | 1.48e-04 | 0 | 0 | 0 |
| TSHR | Thyroid_lineage_TDS | -0.30 | -0.41 | 3.37e-04 | 1 | 0 | 1 |
| ATM | - | +0.22 | +0.41 | 2.20e-03 | 0 | 0 | 0 |
| TP53 | Aggressiveness_EMT | +0.14 | +0.39 | 5.21e-03 | 0 | 0 | 0 |
| ZEB1 | Aggressiveness_EMT | -0.29 | -0.37 | 4.79e-04 | 0 | 0 | 0 |
| CDH2 | Aggressiveness_EMT | -0.68 | -0.35 | 6.76e-03 | 0 | 0 | 0 |
| THRA | Thyroid_lineage_TDS | -0.14 | -0.34 | 8.98e-04 | 0 | 0 | 0 |
| GLIS3 | Thyroid_lineage_TDS | -0.23 | -0.31 | 4.37e-02 | 0 | 0 | 0 |
| NTRK1 | Driver_mutation_fusion_proxy | +0.21 | +0.31 | 2.40e-06 | 0 | 0 | 0 |
| CD8A | - | +0.34 | +0.28 | 2.50e-02 | 1 | 0 | 1 |
| NKX2-1 | Thyroid_lineage_TDS | -0.18 | -0.27 | 8.35e-02 | 1 | 0 | 1 |
| APC | - | -0.15 | -0.26 | 4.72e-02 | 1 | 0 | 1 |
| SLC5A5 | Thyroid_lineage_TDS | -0.44 | -0.26 | 6.33e-02 | 0 | 0 | 0 |
| PHLDA1 | MAPK_output | -0.22 | -0.23 | 7.23e-02 | 0 | 0 | 0 |
| THRB | Thyroid_lineage_TDS | +0.16 | +0.23 | 7.45e-02 | 0 | 0 | 0 |
| ZEB2 | Aggressiveness_EMT | +0.17 | +0.21 | 9.08e-02 | 1 | 1 | 2 |
| ETV5 | MAPK_output | +0.16 | +0.21 | 2.44e-01 | 1 | 0 | 1 |
| KRAS | Driver_mutation_fusion_proxy | +0.07 | +0.20 | 1.72e-01 | 0 | 0 | 0 |
| TERT | Driver_mutation_fusion_proxy | +0.05 | +0.18 | 3.68e-02 | 0 | 0 | 0 |
| SPRY2 | MAPK_output | +0.13 | +0.18 | 2.80e-01 | 1 | 0 | 1 |
| NRAS | Driver_mutation_fusion_proxy | +0.07 | +0.17 | 2.15e-01 | 0 | 0 | 0 |
| PPARG | Driver_mutation_fusion_proxy | -0.16 | -0.17 | 3.17e-01 | 1 | 0 | 1 |
| HRAS | Driver_mutation_fusion_proxy | -0.09 | -0.15 | 3.02e-01 | 0 | 0 | 0 |
| SPRY1 | MAPK_output | +0.10 | +0.11 | 5.07e-01 | 0 | 0 | 0 |
| VIM | Aggressiveness_EMT | +0.06 | +0.11 | 5.62e-01 | 0 | 0 | 0 |
| THADA | - | +0.04 | +0.09 | 7.17e-01 | 0 | 0 | 0 |
| AKT1 | Aggressiveness_EMT | -0.03 | -0.09 | 6.43e-01 | 0 | 0 | 0 |
| BRAF | Driver_mutation_fusion_proxy | -0.02 | -0.07 | 6.84e-01 | 1 | 0 | 1 |
| CALCA | Clinical_IHC | +0.08 | +0.06 | 7.27e-01 | 0 | 0 | 0 |
| EIF1AX | Driver_mutation_fusion_proxy | -0.01 | -0.02 | 9.21e-01 | 0 | 0 | 0 |
| SNAI1 | Aggressiveness_EMT | -0.02 | -0.02 | 9.28e-01 | 0 | 0 | 0 |
| PTEN | Aggressiveness_EMT | -0.00 | -0.01 | 9.78e-01 | 0 | 0 | 0 |
| NTRK3 | Driver_mutation_fusion_proxy | -0.01 | -0.01 | 9.78e-01 | 0 | 0 | 0 |

## Novel candidate biomarkers

| gene | log2FC_TCGA | d | FDR_TCGA | repl_GSE27155 | repl_GSE126698 | novelty_score | median_log2_TCGA |
|------|-------------|---|----------|-----|-----|---|---|
| TACSTD2 | +4.87 | +2.52 | 6.23e-52 | 1 | 1 | 258.50 | 8.16 |
| TMPRSS4 | +4.13 | +2.09 | 3.25e-59 | 1 | 1 | 244.58 | 4.90 |
| PLEKHA6 | +2.35 | +2.22 | 3.04e-54 | 1 | 1 | 237.60 | 3.87 |
| CYP1B1 | +3.47 | +2.22 | 8.99e-53 | 1 | 1 | 231.15 | 5.74 |
| LDLR | +2.24 | +2.02 | 2.76e-53 | 1 | 1 | 211.87 | 4.39 |
| GABRB2 | +3.70 | +2.24 | 2.90e-47 | 1 | 1 | 208.65 | 6.34 |
| B3GNT3 | +3.78 | +2.02 | 8.01e-52 | 1 | 1 | 205.95 | 4.54 |
| PTPRE | +2.38 | +2.26 | 3.83e-45 | 1 | 1 | 200.82 | 6.03 |
| KCNQ3 | +2.48 | +1.96 | 4.95e-51 | 1 | 1 | 197.30 | 4.55 |
| LY6E | +2.13 | +2.16 | 5.84e-43 | 1 | 1 | 182.27 | 7.42 |
| TAGLN2 | +1.22 | +2.18 | 2.35e-42 | 1 | 1 | 181.58 | 8.55 |
| PDLIM4 | +3.44 | +2.25 | 3.07e-40 | 1 | 1 | 178.00 | 6.86 |
| ITGA3 | +1.52 | +2.14 | 5.01e-42 | 1 | 1 | 177.13 | 9.56 |
| BID | +1.35 | +2.13 | 1.74e-41 | 1 | 1 | 173.25 | 5.21 |
| CST6 | +3.87 | +2.04 | 1.17e-41 | 1 | 1 | 167.38 | 5.82 |
| MAMLD1 | +1.95 | +1.94 | 6.65e-42 | 1 | 1 | 160.11 | 3.99 |
| COL8A2 | +2.74 | +2.12 | 2.85e-38 | 1 | 1 | 159.21 | 6.75 |
| MICAL2 | +1.61 | +1.84 | 1.04e-43 | 1 | 1 | 157.98 | 5.45 |
| SFN | +2.92 | +1.78 | 6.43e-45 | 1 | 1 | 157.67 | 3.73 |
| TGFBR1 | +1.44 | +1.81 | 7.17e-44 | 1 | 1 | 156.59 | 6.61 |

### Top 5 novel candidates — short context

- **TACSTD2** — TCGA d=+2.52, log2FC=+4.87, FDR=6.23e-52; replicates in both externals; direction: up in BRAF_like. Gene symbol only, no literature claim in this pipeline.
- **TMPRSS4** — TCGA d=+2.09, log2FC=+4.13, FDR=3.25e-59; replicates in both externals; direction: up in BRAF_like. Gene symbol only, no literature claim in this pipeline.
- **PLEKHA6** — TCGA d=+2.22, log2FC=+2.35, FDR=3.04e-54; replicates in both externals; direction: up in BRAF_like. Gene symbol only, no literature claim in this pipeline.
- **CYP1B1** — TCGA d=+2.22, log2FC=+3.47, FDR=8.99e-53; replicates in both externals; direction: up in BRAF_like. Gene symbol only, no literature claim in this pipeline.
- **LDLR** — TCGA d=+2.02, log2FC=+2.24, FDR=2.76e-53; replicates in both externals; direction: up in BRAF_like. Gene symbol only, no literature claim in this pipeline.

### Top 5 known markers with clearest signal (replicated)

- **SYT12** () — d=+2.47, FDR=1.07e-66, replicates in GSE126698 only.
- **FN1** () — d=+2.45, FDR=4.85e-58, replicates in GSE27155 only.
- **SERPINA1** (BRS_literature) — d=+2.38, FDR=1.68e-37, replicates in both externals.
- **KCNN4** () — d=+2.37, FDR=4.47e-59, replicates in both externals.
- **DCSTAMP** () — d=+2.24, FDR=1.75e-82, replicates in both externals.

## Caveats

- **Label proxies in external cohorts**: GSE27155 subtype labels come from histology-derived mappings already baked into `sample_master.tsv`; GSE126698 labels are inferred from sample-ID prefix (P→PTC→BRAF_like, F→FTC→RAS_like). These are not molecularly-defined subtypes.
- **Small N in GSE126698**: only ~6 BRAF-like and ~6 RAS-like tumors. Per-cohort p-values are low-power; replicate-in-126698 flags should be read as directional support, not independent statistical validation.
- **Microarray coverage**: GSE27155 carries ~13k probes mapped to genes; ~60% of TCGA-tested genes are simply not present there, which mechanically caps `replication_rate` at 1 for many otherwise-valid candidates.
- **No cross-cohort joint test**: FDR is computed per cohort. There is no joint meta-analysis FDR, so the novel-validated set is a conjunction of per-cohort filters, not a joint-model discovery.
- **Novelty is operational, not literature-exhaustive**: 'novel' here strictly means 'not in any hardcoded reference panel or gene list used by this pipeline'. Many genes flagged novel likely have prior thyroid-cancer literature that simply is not encoded here.
- A gene with |d|≈2 in TCGA that fails to replicate in BOTH externals is NOT reported as a biomarker; the novel-validated filter explicitly requires replication_rate≥1.
