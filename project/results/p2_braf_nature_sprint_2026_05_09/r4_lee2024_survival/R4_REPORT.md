# R4 — Lee 2024 GSE213647 survival deep-mine

## Verdict: NOT_RECOVERED (per-sample survival absent)

H8's "blocked" upgraded to "blocked-with-evidence". Lee 2024 (Nat Commun 2024;15:1163, DOI 10.1038/s41467-024-45366-0, PMID 38331894) does not release per-sample OS/PFI/DSS/time-to-event for the Korean SNUBH/CNUH/KRIBB cohort. Reporting Summary (MOESM3) explicitly checks **Clinical data = n/a**. Supplementary Tables 1–3 (MOESM1, Korean N=369 primary tumors) report only **Recurrence** and **Distant metastasis** as group-stratified contingency totals — no time-to-event, no OS, no DSS, no DFI, no per-sample TERT or BRAF/RAS.

## What was checked

1. **GEO** (`series_matrix.txt.gz`, `family.soft.gz`): only 5 sample_characteristics fields = tissue / cell type / cell subtype / **genotype: NA** / **treatment: NA** for all 632 samples.
2. **Local Lee2024_MOESM5.xlsx**: 15 sheets — metabolomics, DEG, KEGG, scRNA, STR. ST13 RNAseq sheet has only Sample ID / GEO ID / Tissue / Fixation / Library kit / Dx_year / Age / Sex.
3. **Springer static-content probe** all MOESM 1–12 × {xlsx,pdf,zip,csv,tsv,txt}: public = MOESM1.pdf (21 pp, Suppl Figs + ST1–ST5), MOESM2.pdf (Peer Review File), MOESM3.pdf (Reporting Summary), MOESM5.xlsx, MOESM6.xlsx (35 sheets, figure source data only). MOESM4 + 7–12 = 403. No source-data zip.
4. **MOESM1 Tables 1–3** (Korean N=369): Recurrence 53/369, Distant met 22/369 — group-totals only across SHMT2-low/high, MTHFD2-low/high, TDS<0/>0. Cannot re-map to our DM1/DM2 calls.
5. **MOESM1 Tables 4–5**: OS/DSS/DFI rows are TCGA-THCA N=500 re-analysis, not Korean. Redundant with our `h6_merged_clinical.tsv`.
6. **cBioPortal API**: no Lee 2024 study ingested (no thyroid / korean / snubh / lee / kribb / cnuh study).

## Implication

- Survival N stays at TCGA-THCA only (h6, n≈513). Lee 2024 cannot double the survival N.
- Lee 2024 still contributes **expression-side Korean replication** (h3 n=632) and 8-gene panel × histology — unaffected.
- **Manuscript recommendation**: cite Lee 2024 for expression-side replication; do NOT claim survival replication; Limitations note Korean OS replication requires author data sharing (Yea Eun Kang, Bon Seok Koo, Seon Kyu Kim) or another deposited Korean RNA-seq + OS cohort (none currently public).

## Alternative paths (out of sprint budget)

1. **Author contact** (IRB-controlled SNUBH/CNUH/KRIBB), 4–12 weeks; not paper-v8-blocking.
2. **Mutect2 BRAF/TERT re-calling from BAMs** (FASTQ already on `/data/thca/v17_korean/GSE213647/per_sample/`) feasible, but survival data still missing → half-useful. Skip unless author contact succeeds.
3. **Korean substitutes**: GSE286332 (n=18) too small; PRJEB11591/K2 (n=180) has BRAF/RAS but no TERT and no public OS. No public Korean cohort currently combines RNA-seq + survival + driver calls.

## Files

- `r4_clinical_table.tsv` (10 rows: Korean ST1–ST3 + TCGA ST4–ST5 totals; no per-sample)
- `r4_survival_replication.tsv` (5 rows, all status=BLOCKED_* with reason codes)
- `_supp/41467_2024_45366_MOESM{1,2,3,5,6}_ESM.*` — full supplementary cache

## Wall-clock

~20 min agent time. Exited within 30-min budget.
