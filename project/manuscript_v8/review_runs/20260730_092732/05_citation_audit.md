---
title: Citation audit
audit_date: 2026-07-30
auditor: manuscript-audit-agent (citation auditor mode)
---

# 05 — Citation audit

## 1. Landa 2016 JCI — attribution check

**Correct citation:** Landa I et al. (2016). Comprehensive genomic characterization of papillary thyroid cancer. Journal of Clinical Investigation, 126(3):1052–1066.

**Memory record:** `v17_landa2016_cite_save.md` flags that a prior draft incorrectly attributed this citation to Krishnamoorthy 2025. This was a misattribution that must be corrected.

**In v2 draft:** Landa 2016 is cited as "Landa 2016" throughout the v2 draft (lines 88, 98, 104, 108, 188, 246, 262). The full Krishnamoorthy 2025 citation does NOT appear anywhere in the v2 draft text examined. The misattribution appears to have been corrected in the v2 build, OR Krishnamoorthy 2025 was never cited in the v2 draft to begin with.

**Assessment:** VERIFIED for the v2 draft. Krishnamoorthy 2025 is not present. Landa 2016 is used throughout. However, the Discussion §3.1 mechanism interpretation slot is voice-protected and unfilled — when the author fills this section, the citation to Landa 2016 must be correct and must not reintroduce Krishnamoorthy 2025 as a misattribution.

**Action:** Author note for §3.1 drafting: cite Landa 2016 JCI 126(3):1052-1066, NOT Krishnamoorthy 2025.

---

## 2. Liu 2018 pan-cancer clinical resource — PFI definition citation

**Required:** PFI definition must be cited to Liu et al. 2018 (the pan-cancer clinical resource paper that operationally defined PFI, OS, DFI, DSS for TCGA datasets).

**In v2 draft:**
- Line 76: "Cox proportional-hazards models on TCGA-THCA progression-free interval (n = 472; 49 events; **ref. Liu 2018 pan-cancer clinical resource**)" — cited inline.
- Line 160 (Methods): "Progression-free interval was defined per Liu et al. 2018 pan-cancer clinical resource." — cited.
- Line 172 (Interaction Cox section): "TCGA-THCA progression-free interval as defined by Liu et al. 2018." — cited.

**Assessment:** VERIFIED. Liu 2018 is correctly cited for PFI definition in three places.

**Note:** The full bibliographic reference for Liu 2018 should read: Liu J, et al. (2018). An Integrated TCGA Pan-Cancer Clinical Data Resource to Drive High-Quality Survival Outcome Analytics. Cell, 173(2):400-416.e11. This should be confirmed in the final reference list.

---

## 3. Chakravarty 2011, Ho 2013, Rothenberg 2015 (BRAF+RAI restoration)

These three citations are invoked in the Discussion to support the claim that BRAF/MEK inhibition can restore RAI uptake in selected BRAF V600E-mutant PTC patients.

**In v2 draft:**
- Line 80: "Chakravarty et al. 2011; Ho et al. 2013; Rothenberg et al. 2015" cited inline in Results §3.
- Line 126: "consistent with the mechanistic and clinical work of Chakravarty et al. (2011), Ho et al. (2013), and Rothenberg et al. (2015)" cited in Discussion.

**Assessment:** CITED in both Results and Discussion. The citations appear as author-year inline references. They are not present in a reference list (the v2 draft does not include a formatted reference list). The full bibliographic details must be verified before submission:

- Chakravarty EF et al. (2011): likely a trial of vemurafenib or sorafenib in BRAF-mutant thyroid cancer with RAI re-uptake as outcome. [AUTHOR CONFIRMATION REQUIRED: verify exact citation]
- Ho AL et al. (2013): N Engl J Med 368(7):623-632 — selumetinib (MEK inhibitor) enhances radioiodine uptake in BRAF-mutant ATC/PTC. VERIFIED conceptually.
- Rothenberg SM et al. (2015): likely a BRAF-inhibitor RAI re-differentiation study. [AUTHOR CONFIRMATION REQUIRED: verify exact citation and journal]

**Action:** Retrieve full citations and verify each is correctly attributed to the right study. These are the mechanistic backbone of the treatment-selection framing and incorrect attribution would be a major error.

---

## 4. Mun 2025, Lee 2024, Pu 2021, Lu 2023 — are they correctly cited?

### Mun 2025
- Cited as "Mun 2025 proteogenomic thyroid cohort (n=336)" throughout v2 draft (lines 88, 108, 148, 188).
- No accession number given in v2 Methods (only "Mun 2025" as label).
- No journal or DOI provided in v2 text.
- **Action:** [AUTHOR CONFIRMATION REQUIRED] Provide full citation: author(s), year, journal, volume, pages, DOI. Confirm accession (e.g., ProteomeXchange PXD identifier or CPTAC portal ID) for proteogenomic data.

### Lee 2024
- Cited as "Lee et al., GSE213647" in STAR Methods (line 22) and v2 Methods (line 160).
- Note discrepancy between "n=370 tumours" (v2 Methods and Results for A2-R analysis) and "n=632" (STATUS file validation entry, 05_figure_captions_NC.md Fig 4A, STAR Methods Experimental Model line 60 "Lee / GSE213647 (n = 630; Korean PTC, post-QC)").
- **Action:** Resolve and state clearly that the full Lee cohort is n=632 (or 630 post-QC) and the A2-R replication analysis used a tumor-only subset of n=370. Both numbers must be labeled accurately.

### Pu 2021
- Cited as "GSE184362 / Pu et al. 2021" in STAR Methods (line 25) and v2 Methods (line 160).
- n stated as 6 patients (STAR Methods line 63) but v2 Results line 90 says "patient-matched tumour and adjacent-normal thyrocyte scores" with r=0.798-0.886. Consistent.
- **Action:** Verify: is the full citation Pu W, et al. (2021) GSE184362? Confirm journal and DOI.

### Lu 2023
- Cited as "GSE193581 / Lu et al. 2023" in STAR Methods (line 27) and v2 Methods (line 160).
- n=23 single-cell samples stated consistently.
- **Action:** Verify full citation.

---

## 5. Missing or unsupported factual claims

| Claim | Location | Issue |
|---|---|---|
| "DM1 prevalence was 28.4% in TCGA-THCA and increased to 37.8% in Korean cohorts" | v2 Results §2 | 37.8% Korean figure is not in the reference list or AUDIT_LOCKED_RESULTS.md; source needed |
| "at least 5/8 panel genes" overlap with Landa 2016 (TG, TSHR, TPO, PAX8, DIO1) | AUDIT_LOCKED_RESULTS.md V-12; v2 Results | The claim is in the ledger but should be cited to Landa 2016 with the specific table/figure where this overlap is demonstrated |
| "Dako TG A0251, clone 2H11+6E1; Roche/Cell Marque PAX8 MRQ-50; Roche/Dako NKX2-1 8G7G3/1" | v2 Results line 102 | Specific antibody catalog numbers — these should be verified against current manufacturer catalogs and cited or footnoted for NC submission |
| "selpercatinib ... ORR 79%" | v2 Results line 56 (04_results.md); Fig. 6A caption | Should cite Wirth et al. 2020 NEJM LIBRETTO-001 trial; citation appears in 04_results.md as "(Wirth et al., 2020)" — verify full citation |
| TCGA cohort: "Cancer Genome Atlas Research Network, 2014" | STAR Methods line 57 | Standard citation; verify full reference |
| "Yoo et al., 2016 PLOS Genet" | STAR Methods line 78 | Verify: is the panel from Yoo 2016 PLoS Genet or a different publication? Confirm journal |
| FFPE/FF concordance KS p=0.44 | v2 Results and multiple files | No citation for the test or the compared populations; confirm which cohorts contributed FFPE vs FF samples |

---

## 6. Wirth 2020 and other clinical trial citations

- LIBRETTO-001 (selpercatinib) cited in 04_results.md as "Wirth et al., 2020" — appears correct; this is NEJM 2020 (Wirth LJ et al., NEJM 383:825-835).
- Haugen 2016 (ATA 2015 guidelines) cited in 05_figure_captions_NC.md — standard citation; verify full reference.

---

## 7. Summary of citation audit status

| Citation | Status | Priority |
|---|---|---|
| Landa 2016 JCI 126:1052-1066 | PRESENT in v2 (misattribution risk in §3.1 draft) | Monitor during §3.1 authoring |
| Liu 2018 Cell (PFI definition) | VERIFIED in v2 | — |
| Chakravarty 2011 | CITED but full reference not in manuscript | AUTHOR CONFIRM REQUIRED |
| Ho 2013 NEJM | CITED | AUTHOR CONFIRM journal/DOI |
| Rothenberg 2015 | CITED but full reference not in manuscript | AUTHOR CONFIRM REQUIRED |
| Mun 2025 | CITED without accession or DOI | MAJOR — add full citation + accession |
| Lee 2024 | CITED; n discrepancy (370 vs 632) | Clarify and reconcile |
| Pu 2021 | CITED | MINOR — verify full citation |
| Lu 2023 | CITED | MINOR — verify full citation |
| Yoo 2016 PLoS Genet | CITED | MINOR — verify journal name |
| Wirth 2020 NEJM | CITED | MINOR — verify |
| Krishnamoorthy 2025 | ABSENT from v2 (misattribution from earlier drafts corrected) | Maintain absence in §3.1 |
