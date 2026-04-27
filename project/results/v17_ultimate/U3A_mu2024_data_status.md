# U3A — Mu Z et al. JCEM 2024 data accessibility status

_Date: 2026-04-27_  •  _Author: Seungho Cook_  •  _v17 ULTIMATE sprint_

## TL;DR

**Mu 2024 supplementary / raw data accessible? NO (controlled).**

The Oxford Academic / PMC version of the paper is open-access for *reading*,
but the underlying NGS data is deposited under **NGDC accession
`HRA004166`** with **controlled access** — a Data Access Committee (DAC)
application is required. The published supplement does not contain a
patient-level mutation × RAI-avidity table that would allow us to bypass
the DAC.

## Citation

- Mu Z, Zhang X, Sun D, Sun Y, Shi C, et al.. *Characterizing Genetic Alterations Related to Radioiodine Avidity in Metastatic Thyroid Cancer*.
  **J Clin Endocrinol Metab** 2024;109(5):1231-1240.
- PMID: 38060243  •  PMCID: PMC11031230  •  DOI: 10.1210/clinem/dgad697

## Cohort

- N patients: **214** (DM-DTC)
- N samples : **281**
- RAI labels: I-RAIA (n=134) vs I-RAIR (n=80); subgroups C-RAIA / P-RAIR / G-RAIR

## Access verdict

- Data accession: **HRA004166 (NGDC GSA-Human)**
- Access tier  : **CONTROLLED — DAC application required**
- DAC contact  : Lin Yansong, linys@pumch.cn (Peking Union Medical College Hospital)
- Reason       : Data Availability: 'Restrictions apply ... to preserve patient confidentiality'. NGS deposited under HRA004166 at NGDC requires Data Access Committee approval; no openly downloadable variant table or supplementary mutation matrix on Oxford Academic / PMC.

## What we did NOT do (and why)

- We did **not** attempt Sci-Hub or other paywall-bypass mirrors.
  The paper itself is already open-access on PMC; the controlled item is
  patient-level NGS, where bypass would be both unethical and illegal.
- We did **not** scrape NGDC. A DAC request is the correct path.

## Recommended next step

1. Lean on the open RAI-avid/refractory cohorts in `U3A_rai_landscape.tsv`
   (GSE151181 SuperSeries, n=99) as the primary external validation set —
   already in flight in U1B.
2. **Optionally** submit a DAC request to Lin Yansong (linys@pumch.cn) for
   HRA004166 if reviewers ask for clinical-grade RAI labels beyond GEO.
3. Treat Mu 2024 as **citation + concordance reference**, not as a
   downloadable test set, in the manuscript.
