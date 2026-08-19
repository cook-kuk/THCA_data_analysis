# Dataset inventory — RAI response genomics atlas

This is the human-readable companion to `results/tables/dataset_inventory.csv`. Each row maps a dataset to its label tier and to the manuscript's evidence ladder. **Always cite the tier** in derived figures and tables.

| Dataset | Accession | Cancer | n | Omics | RAI label | Tier | Access | Priority | Next action |
|---------|-----------|--------|---|-------|-----------|------|--------|----------|-------------|
| GSE151179 | GEO | PTC | ~50 (verify) | bulk mRNA | RAI-avid vs refractory; sample type | 2 | public | 1 | download + validate panel |
| GSE151180 | GEO | PTC | matched | miRNA | same as GSE151179 | 2 | public | 3 | defer until mRNA done |
| GSE151181 | GEO | — | — | SuperSeries | — | — | public | — | navigation only |
| GSE299988 | GEO | PTC | small (~10–30) | bulk expression | RAI-avid vs non-avid; +normal | 2 | public | 2 | download supportive |
| Boucai 2023 CCR | request | metastatic DTC | exceptional responders cohort | RNA-seq + WES | RECIST exceptional vs non-responder | 1 | supplement + request | 1 | data_request_boucai.md |
| Mu 2024 JCEM | HRA004166 | distant metastatic DTC | 220 | targeted NGS | 4-class RAI uptake patterns | 2 (gray-zone) | controlled | 1 | data_request_mu_hra004166.md |
| TCGA-THCA | GDC | PTC majority | ~504 | RNA-seq + WES + HM450 + CNV + clinical | none (proxies only) | 4 | public | 1 | discovery context, no overclaim |
| GSE112202 | GEO | thyroid (cells/tumors) | ? | bulk expression | digoxin redifferentiation | 5 | public | 2 | check panel directionality |
| GSE184362 | GEO | PTC | scRNA, ~149k cells | scRNA-seq | primary / paratumor / LN met / RAI-refractory distant met | 2 + cell-state | public | 3 | plan scRNA workflow first |
| GSE76039 (Landa) | GEO | PDTC + ATC | ~85 | bulk expression + HM450 | aggressive vs not | 4 | public | reserve | already used in Paper 1 R17 layer reserve |
| GSE286332 | GEO | Korean PTC vs PTC+HT | 18 | bulk RNA | HT-overlap (Paper 1 standard) | 4 | public | reserve | already used in Paper 1 R17 L3 |
| Mun 2025 proteome | published | DTC + PDTC + ATC | 336 | proteome | dedifferentiation modules | 4 | public | reserve | already used in Paper 1 R17 L4 |
| Lu 2023 GSE193581 | GEO | PTC + ATC | scRNA, ~67k | scRNA | per-sample histology, no direct RAI | 4 | public | reserve | already used in Paper 1 R17 L2 |
| Redifferentiation trials | various | DTC RAIR | ~10–20 each | bulk RNA + RAI imaging | restored uptake post-MAPKi | 5 | supplement only | 2 | scan trial supplements |

## Notes per dataset

### GSE151179 (priority 1, label-anchored RAI cohort)

This is the **first validation target**. Verify on download:

- exact sample count (announced as PTC primary + lymph-node metastasis with RAI-refractory annotation in the source paper, but counts vary by GEO record version).
- normalization status (already log2 / quantile-normalized vs raw).
- platform — affects gene-to-probe mapping and 8-gene panel coverage.
- label fields: search `characteristics_ch1` for `iodine`, `i-131`, `rai`, `avid`, `refractory`, `remission`, `persistence`, `uptake`.

### GSE299988 (priority 2, small supportive)

Treat as **supportive validation only**, not primary proof. Sample size likely < 30. If panel directionality holds, include as supplementary; if it doesn't, document honestly and probe for confounders (platform, batch, normal-contamination).

### Mu 2024 / HRA004166 (priority 1, gray zone)

The four-class uptake pattern data is **the single best argument** for the manuscript's "molecular gray zone" framing. Even if controlled access blocks per-patient genomics, **the supplementary tables** (likely PDF/XLSX) usually publish mutation-frequency tables by uptake class — these alone enable a driver-overlap argument. See `docs/data_request_mu_hra004166.md`.

### Boucai 2023 CCR (priority 1, Tier 1 label)

The only public source with **RECIST-based RAI response** in metastatic DTC at meaningful sample size. Pursue both supplementary tables (for definitions + gene lists) and a data-sharing request for raw expression (`docs/data_request_boucai.md`).

### TCGA-THCA (priority 1, Tier 4 proxy only)

Discovery cohort. Use for:

- panel score distribution across drivers.
- methylation × panel × driver three-way analysis (already done in Paper 1).
- weak survival proxies (PFI / DFI / recurrence) with honest under-powered disclosure.

**Never** claim TCGA-THCA validates the panel against RAI response.

### GSE184362 (priority 3, scRNA mechanism)

Single-cell. The Pu 2021 cohort includes per-sample tissue annotation (primary tumor / paratumor / LN met / **RAI-refractory distant metastasis**). The RAI-refractory distant met sample is a unique cell-state anchor. Heavy compute (~150k cells); plan before downloading raw.

### Redifferentiation trials (priority 2, Tier 5)

Often only patient-level summary tables in supplements. Useful for **direction-of-effect** demonstration: panel score should *rise* under successful redifferentiation. Selumetinib (Ho 2013 NEJM), dabrafenib/trametinib (Rothenberg 2015 + MERAIODE), vemurafenib (Brose 2016 Lancet Oncol).

## Cross-reference to Paper 1 R17 work (existing)

Datasets that have already been characterized in Paper 1's R17 two-axis reconciliation work and live in `project/results/r17_tcga_panel_d4p2_reconciliation/`:

- TCGA bulk (R17 L1).
- Lu 2023 sc (R17 L2).
- GSE286332 Korean (R17 L3).
- Mun 2025 proteome (R17 L4).
- TCGA TERT × R17 zone (R17 L5).
- TCGA Cox per zone (R17 L6).
- TCGA per-zone 8-gene profile (R17 L7).
- TCGA HM450 β × R17 zone (R17 L8).

These are reused as the **discovery / mechanism** half of the RAI atlas manuscript. The new contribution of this workspace is the **Tier 1/2 label-anchored validation** against GSE151179 + Boucai + Mu + GSE299988.
