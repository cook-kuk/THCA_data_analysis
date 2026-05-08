# Neoantigen dataset inventory

Generated: 2026-05-07 17:42:17
Total master records: **105,058**

## Source coverage

```
                      source_dataset  n_records  n_with_peptide  n_with_hla  n_positive  n_negative
                     BigMHC_example1          9               9           9           4           5
                               CEDAR       5000            4985        4941        4837         163
                       IEDB_tcell_v3     100000          100000       88577       38314       61686
    NeoRanking_Gartner_nmers_ranking         46              46           0          46           0
        clinical_vaccine_BNT122_PDAC          1               0           0           1           0
       clinical_vaccine_ELI002_mKRAS          1               0           0           1           0
clinical_vaccine_KEYNOTE942_melanoma          1               0           0           1           0
```

## Validation level distribution

| Level | Records |
|---|---|
| LEVEL_4_TCELL_ASSAY | 100,037 |
| LEVEL_2_MS_ELUTED_LIGAND | 4,794 |
| LEVEL_5_PATIENT_MATCHED_TCELL | 182 |
| LEVEL_1_BINDING_ASSAY | 27 |
| LEVEL_7_TCR_PMHC_STRUCTURE | 15 |
| LEVEL_6_CLINICAL_VACCINE_RESPONSE | 3 |


## Per-source × validation level

```
                      source_dataset                  validation_level      n
                     BigMHC_example1             LEVEL_1_BINDING_ASSAY      9
                               CEDAR             LEVEL_1_BINDING_ASSAY     18
                               CEDAR          LEVEL_2_MS_ELUTED_LIGAND   4794
                               CEDAR               LEVEL_4_TCELL_ASSAY     37
                               CEDAR     LEVEL_5_PATIENT_MATCHED_TCELL    136
                               CEDAR        LEVEL_7_TCR_PMHC_STRUCTURE     15
                       IEDB_tcell_v3               LEVEL_4_TCELL_ASSAY 100000
    NeoRanking_Gartner_nmers_ranking     LEVEL_5_PATIENT_MATCHED_TCELL     46
        clinical_vaccine_BNT122_PDAC LEVEL_6_CLINICAL_VACCINE_RESPONSE      1
       clinical_vaccine_ELI002_mKRAS LEVEL_6_CLINICAL_VACCINE_RESPONSE      1
clinical_vaccine_KEYNOTE942_melanoma LEVEL_6_CLINICAL_VACCINE_RESPONSE      1
```

## Status by source

| Source | Status | Count |
|---|---|---|
| CEDAR | SUCCESS | {per_source.get('CEDAR',0):,} |
| IEDB | PARTIAL (deduplicated against CEDAR) | (overlapping, not loaded separately) |
| BigMHC GitHub | PARTIAL (example only) | {per_source.get('BigMHC_example1',0)} |
| BigMHC Mendeley bulk | MANUAL_REQUIRED | manual download from `data.mendeley.com/datasets/dvmz6pkzvb/4` |
| TESLA | MANUAL_REQUIRED (Cell suppl) | 0 (R scripts only on GitHub) |
| NeoRanking Gartner | SUCCESS | {per_source.get('NeoRanking_Gartner_nmers_ranking',0)} |
| NEPdb | MANUAL_REQUIRED | site reachable, no bulk URL |
| dbPepNeo | BROKEN_URL | DNS down |
| TSNAdb | MANUAL_REQUIRED | needs download-button click |
| Neodb | BROKEN_URL | DNS down |
| PRIME / MixMHCpred | BLOCKED_BY_LICENSE | predictor binaries gated |
| VDJdb | SUCCESS (cloned, not yet ingested) | TCR layer, separate table |
| McPAS-TCR | MANUAL_REQUIRED | email-gated form |
| NeoTCR | SUCCESS (cloned, not yet ingested) | TCR layer |
| Clinical vaccine evidence | SUCCESS | 3 (curated paper list) |

## Manual download list

See `data_sources/manual_downloads.md` for exact procedures for each.

## Next actions (in priority order)

1. **TESLA Cell supplementary** (Wells 2020) — required for canonical patient-matched ground truth.
2. **BigMHC Mendeley bulk** — needed for el_train/im_train predictor benchmarking.
3. **NEPdb form-gated download** — needed for negative labels (positive + negative T-cell-validated).
4. **TSNAdb download button** — needed for SNV/INDEL/fusion-derived neoantigens.
5. **PRIME license application** — needed for production-grade prediction.
6. **Yonsei (Korean PDAC) cohort** — drop-in slot already prepared at scripts/21_yonsei_dropin.py.

## Thyroid relevance

Records with thyroid driver gene match (BRAF/RAS/TERT/RET/NTRK/ALK/TP53/EIF1AX/DICER1/PAX8/PPARG/PTEN/SMAD4/CDKN2A): **{master_status['n_thyroid_relevant']}**

Output: `data_processed/thyroid_relevance_candidates.csv`.

⚠ NONE of these has clinical-vaccine-validated efficacy in thyroid cancer as of 2026-05-07. Treat as exploratory candidates.
