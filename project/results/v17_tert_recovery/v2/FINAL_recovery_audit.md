# v17 TERT promoter recovery — FINAL audit (v2, 9-source)

_Run date: 2026-04-27_
_Working dir: `project/results/v17_tert_recovery/v2/`_

## Headline: Scenario A (best case) — 36 TCGA-THCA TERT promoter mutations recovered

The v1 attempt (`recovery_status.md`) checked only `thca_tcga` and
`thca_tcga_pan_can_atlas_2018` cBioPortal studies and found 4 TERT mutations,
all non-promoter. The v1 conclusion was that TCGA-THCA WES does not capture
TERT promoter (true for those two studies).

**v2 finding:** the cBioPortal `thca_tcga_pub` study — built from the TCGA
Thyroid Cancer 2014 *Cell* paper publication MAF — contains **36 TERT
promoter mutations** at canonical hotspots (chr5:1295228 = C228T, chr5:1295250
= C250T) for TCGA-THCA samples. The 2014 publication MAF added Sanger-validated
TERT promoter calls that the WES-derived legacy MAFs lacked. The v1 sweep
did not check this study.

## Headline metrics

| Metric | Value |
|---|---|
| Sources attempted | 9 |
| TCGA-THCA TERT promoter mutations recovered | **36** (7.0% of 513) |
| External (non-TCGA) TERT promoter records | **1,844** (MSK-IMPACT 1361, MSK-CHORD 342, MSK PDTC/ATC 81, MSK 2016 60, ODG 21) |
| Canonical positions | C228T (chr5:1295228) and C250T (chr5:1295250), GRCh37 |
| Wall time | 10.4 s (asyncio orchestrator) |

## Per-source result

| Source | Records | Promoter | TCGA matched | Status |
|---|---:|---:|---:|---|
| S1 Liu 2017 JCO supplement | 0 | 0 | 0 | All 8 candidate URLs returned 0-byte / 404 / non-data HTML |
| S2 Landa 2016 PDTC/ATC | 0 | 0 | 0 | JCI URL guesses failed; cBioPortal `thca_mskcc_2016` covered via S6 instead |
| S3 Yoo 2019 Korean cohort | 0 | 0 | 0 | Springer ESM URLs returned non-XLSX bodies; manual download required |
| S4 Pozdeyev 2024 / GENIE | 61 | 60 | 0 | Pozdeyev 2024 paper not located via DOI guess; recovered TERT promoter from MSK PDTC/ATC and `odg_msk_2017` |
| S5 COSMIC | 0 | 0 | 0 | Public pages returned no machine-readable promoter data; bulk requires academic license (guide written) |
| **S6 cBioPortal all-thyroid** | **2,534** | **1,820** | **36** | Iterated 9 thyroid-relevant studies, including the previously-missed `thca_tcga_pub` |
| S7 GDC controlled access | 10,000 | 0 | 0 | Listed 470 unique TCGA-THCA cases with WGS BAMs (controlled); dbGaP application guide written |
| S8 bioRxiv / EuropePMC preprints | 50 | 0 | 0 | Found 50 TERT/thyroid papers; sample-level extraction needs manual triage |
| S9 PubMed publication mining | 0 | 0 | 0 | 68 PubMed IDs searched, no auto-extracted TCGA-barcode/TERT pairing |

## TERT+ patient profile (n=36)

```
Quad group (out of 281 BRAF, 54 RAS, 178 triple-neg in cohort):
  BRAF only            25  (8.9% of BRAF)
  RAS only              6  (11.1% of RAS)
  Triple negative       5  (2.8% of triple-neg)

AJCC stage:
  Stage I               4
  Stage II              4
  Stage III            10
  Stage IVA             8
  Stage IVC             3
  Stage IV              1
  unknown/blank         6

TDS16 differentiation score (lower = more dedifferentiated):
  TERT+ mean    6.238
  WT mean       6.879
  delta        -0.641   (TERT+ less differentiated, consistent with biology)
```

22 / 30 staged TERT+ patients are stage III/IV (73%) — strong stage skew vs
the cohort baseline.

## Cross-tabs

### v17 dark-matter (DM1 vs DM2) cluster

DM clustering only covers 178 of 513 patients (the dark-matter-positive
subcohort). Among those:

```
                  mutated  wildtype
DM1                     4       106
DM2                     1        68
unassigned (non-DM)    31       303
```

Fisher DM1 vs DM2: odds = 2.57, **p = 0.65** — under-powered with only 5
TERT+ falling inside the DM subcohort.

### Quad-group (full cohort)

```
                  mutated  wildtype
A_braf_only            25       256
B_ras_only              6        48
D_triple_negative       5       173
```

### driver_anchor

```
        mutated  wildtype
BRAF         25       256
RAS           6        48
NTRK          1         0
TP53          0         1
unknown       4       172
```

### molecular_subtype

```
          mutated  wildtype
BRAF_like      29       363
RAS_like        7       104
unknown         0        10
```

## Survival (overall TCGA-THCA cohort, n=504 with OS)

```
TERT+ vs WT logrank        p = 4.92e-06
  TERT+ events:    6 / 36   (16.7%)
  WT events:      10 / 468  (2.1%)
  Test stat:      20.87

4-group (BRAF / RAS / TERT+ / Triple-neg) logrank
  p = 3.78e-05
  Test stat:     23.14
  Group sizes:   BRAF_only 250, RAS_only 48, TERT+ 36, Triple_neg 170
```

These are strong p-values driven by the elevated event rate in TERT+
(~8x baseline). The absolute event count (6) is small, so hazard ratios will
have wide confidence intervals — present with caution and include CI bands.

## Caveats and audit trail

- All 9 sources logged URL-by-URL in `logs/*_attempts.json` (success and
  failure, with status codes, byte counts, timings).
- Data provenance: 36 TCGA-THCA promoter calls all originate from
  `thca_tcga_pub` (cBioPortal mirror of the 2014 *Cell* paper). These were
  produced by Sanger validation in the original publication; we did not
  re-call from raw BAMs.
- The Liu 2017 *JCO* follow-up paper (potentially 100+ additional Sanger
  calls) was **not** recovered — all 8 candidate supplement URLs failed.
  Manual outreach to Liu R or an institutional library remains the highest-
  leverage outstanding source.
- 1,844 external TERT-promoter records exist (MSK-IMPACT 1361, MSK-CHORD 342,
  MSK 2016 60, PDTC/ATC 81, ODG 21). These cannot be merged into TCGA-THCA at
  the patient level (different cohorts) but are valid as **external validation
  cohorts** for any TERT-related signature derived from TCGA.
- Pozdeyev 2024 *Nature Cancer* paper not located via automated DOI guesses.
- Yoo 2019 Korean cohort URLs returned non-data bodies; sample-level Korean
  TERT validation remains unavailable from automated retrieval.

## Files

| File | Contents |
|---|---|
| `FINAL_tert_status_integrated.tsv` | 36-row integrated TERT-mutated patient list |
| `FINAL_promoter_records_all_sources.tsv` | 1,880 promoter mutation records (incl. external) |
| `sample_master_v17_tert_v2.tsv` | v17 sample_master with `tert_promoter_v2` and `tert_promoter_integrated` columns + merged `v17_dark_cluster` |
| `FINAL_crosstab_quad_group.tsv` | Quad-group cross-tab |
| `FINAL_crosstab_4group.tsv` | BRAF / RAS / TERT+ / Triple-neg |
| `FINAL_crosstab_aggressive.tsv` | aggressive_flag cross-tab |
| `FINAL_crosstab_molecular_subtype.tsv` | molecular subtype cross-tab |
| `FINAL_crosstab_driver_anchor.tsv` | driver anchor cross-tab |
| `FINAL_extended_summary.json` | Survival + Fisher + crosstab counts in JSON |
| `parsed/S6_cbioportal_all_studies_summary.tsv` | Per-study TERT counts |
| `parsed/S6_cbioportal_all_promoter_mutations.tsv` | All 1,820 promoter calls |
| `S5_cosmic_application_guide.md` | COSMIC academic-access workflow |
| `S7_gdc_controlled_application_guide.md` | dbGaP TCGA-THCA WGS workflow |
| `orchestrator_summary.json` | Per-source orchestrator output |

## Verdict

**Scenario A** — TERT promoter status recovered for 36 TCGA-THCA patients
(7.0% of cohort) with strong, biologically plausible association with stage
and survival (logrank p = 4.9e-6). The 4-group analysis (BRAF / RAS / TERT+ /
Triple-neg) is now feasible in the v17 paper main text. External validation
cohort (MSK-IMPACT etc.) available for sensitivity analysis at the cohort
(not patient-matched) level.
