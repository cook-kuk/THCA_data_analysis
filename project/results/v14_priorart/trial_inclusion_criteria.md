# Task 3 — Sacituzumab thyroid trial inclusion criteria (NCT06235216, NCT07521670)

_Date: 2026-04-25 · Source: ClinicalTrials.gov v2 API, retrieved live._

## NCT06235216 — SETHY (Sacituzumab govitEcan in THYroid Cancers)

| Field | Value |
|---|---|
| Sponsor | Grupo Espanol de Tumores Neuroendocrinos |
| Status | **RECRUITING** |
| Phase | 2 |
| Started | 2024-09-13 |
| Primary completion | 2027-12 |
| Drug | Sacituzumab govitecan, 10 mg/kg IV |
| Primary outcome | Objective Response Rate (ORR) |
| Conditions | Differentiated Thyroid Cancer (cohort A); Anaplastic Thyroid Cancer (cohort B) |

### BRAF stratification?  **NO.**
Inclusion criteria specify radioactive-iodine refractory differentiated TC (cohort A, must have progressed on Sorafenib/Lenvatinib/Cabozantinib, ≤3 prior lines) or ATC (cohort B, ≤1 prior line). **No BRAF V600E genotype requirement; no TROP2 IHC requirement.** Archival tumour tissue is requested for translational studies but not gated on biomarker status.

### TROP2 expression requirement?  **NO** (archival sample requested for translational analysis only).

## NCT07521670 — STRAP (Sacituzumab Tirumotecan in R/M Adenoid Cystic Carcinoma and PTC)

| Field | Value |
|---|---|
| Sponsor | National Cancer Centre, Singapore |
| Status | **NOT_YET_RECRUITING** |
| Phase | 2 |
| Started (planned) | 2026-05 |
| Primary completion | 2028-10-31 |
| Drug | Sacituzumab Tirumotecan, 4 mg/kg IV (D1, D15 of 28-day cycle) |
| Primary outcome | Objective Response Rate (ORR), investigator-assessed |
| Conditions | Recurrent/Metastatic ACC; PTC |

### BRAF stratification?  **NO.**
No BRAF genotype requirement. Cohort B is "papillary thyroid carcinoma" without further molecular stratification.

### TROP2 expression requirement?  **EXPLICITLY WAIVED.**
Verbatim from inclusion criteria #2: *"Trophoblast cell-surface antigen 2 (TROP2) expression testing by immunohistochemistry or other methods is not required for enrollment. Provision of archival tumour tissue (where available) will be requested to support retrospective analysis."*

This is the strongest design fact in our favour — STRAP **explicitly enrols TROP2-unselected**, with retrospective TROP2 analysis only.

## Implication for v14 contribution

Both trials enrol **without BRAF or TROP2 stratification** and use ORR as the primary endpoint. This means:

- If sacituzumab govitecan/tirumotecan responds preferentially in BRAF-like × TROP2-high subset (our hypothesis), the unselected ORR will be **diluted** by RAS-like or TROP2-low non-responders.
- A **prospective enrichment** or a **stratified secondary analysis** (with biomarker-positive subgroup ORR as a key supporting endpoint) would be the v14-aligned trial-design proposal.
- Both trials request archival tissue for translational analysis — this is the natural insertion point for the BRAF-like + TACSTD2-high two-gate retrospective overlay.

## Recommended language for v14 paper Discussion §6.X

> "Two ongoing or planned phase-2 trials of TROP2-directed antibody–drug conjugates in thyroid cancer (NCT06235216, NCT07521670) currently enrol without BRAF genotype stratification, and STRAP (NCT07521670) explicitly waives TROP2 IHC at enrollment. Our v14 cell-line and pharmacogenomic data provide a pre-specified hypothesis for retrospective biomarker-positive subgroup analysis (ORR in BRAF-mutant ∩ TACSTD2-high vs. complement) using the archival tissue collected per protocol; if borne out, a prospective enrichment design would be the natural successor study."
