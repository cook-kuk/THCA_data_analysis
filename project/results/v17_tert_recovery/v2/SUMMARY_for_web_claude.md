# v17 TERT promoter recovery — Summary for strategic judgment

_Run: 2026-04-27 · 9-source asyncio sweep · Wall time 10.4 s_

---

## TL;DR

A second-pass exhaustive 9-source sweep recovered **36 TCGA-THCA TERT
promoter mutations** that the first attempt missed. The first attempt
correctly concluded TCGA-THCA WES does not capture TERT promoter, but
checked only `thca_tcga` and `thca_tcga_pan_can_atlas_2018` on cBioPortal.
The second attempt found that the cBioPortal `thca_tcga_pub` study (mirror
of the 2014 *Cell* paper publication MAF) contains **Sanger-validated TERT
promoter calls** at canonical positions chr5:1295228 (C228T) and chr5:1295250
(C250T) for 36 TCGA-THCA patients.

**Survival signal is strong:** TERT+ vs WT logrank **p = 4.9e-6** (events
6/36 = 16.7% vs 10/468 = 2.1%, ~8x event rate).

This makes the 4-group narrative (BRAF / RAS / +TERT / Triple-neg) defensible
in the v17 paper main text, which was the goal of the sprint.

---

## Sprint scope

The sprint was framed as a "fire-and-forget last attempt" before declaring
TERT a fundamental limitation. Three planned scenarios:

| Scenario | Expected outcome | What actually happened |
|---|---|---|
| A — best (≥30 patients) | Paper figure, npj P 70%→78% | **MET. 36 patients.** |
| B — realistic (10–29) | Supplementary figure, P 70%→73% | n/a |
| C — worst (<10) | Limitations strengthened, no P shift | n/a |

---

## What recovered (per source)

| # | Source | TCGA promoter | External promoter | Notes |
|---|---|---:|---:|---|
| S1 | Liu 2017 JCO supplement (8 URL guesses) | 0 | 0 | All URLs dead — needs manual library access |
| S2 | Landa 2016 PDTC/ATC (JCI + cBioPortal) | 0 | covered by S6 | JCI URL guesses failed |
| S3 | Yoo 2019 Korean cohort (Nat Comm) | 0 | 0 | Springer ESM URLs returned non-data — needs manual download |
| S4 | Pozdeyev 2024 Nat Cancer (DOI guess + GENIE) | 0 | 60 | Paper not located via DOI guess; recovered TERT from MSK substitute panels |
| S5 | COSMIC v100 (public scrape) | 0 | 0 | Bulk requires academic license; application guide written |
| **S6** | **cBioPortal all-thyroid (9 studies)** | **36** | **1,820** | **Found the missing `thca_tcga_pub` study + MSK-IMPACT/CHORD** |
| S7 | GDC controlled-access scout (TCGA WGS) | 0 | 0 | 470 TCGA-THCA WGS cases listed; dbGaP DAR needed (4–8 weeks) |
| S8 | bioRxiv / EuropePMC preprints | 0 | 0 | 50 papers found; sample-level extraction needs manual triage |
| S9 | PubMed publication mining (PMC OA) | 0 | 0 | 68 PMIDs searched, no auto-extracted TCGA/TERT pairs |

External promoter records (1,844 total): **MSK-IMPACT 1,361 · MSK-CHORD 342 ·
MSK-PDTC-ATC 81 · MSK-thyroid-2016 60 · ODG-MSK-2017 21**. Cannot patient-
match to TCGA (different cohorts) but valid as **cohort-level external
validation** for prevalence sanity checks.

---

## TERT+ patient profile (n=36)

```
Quad group (% of group that is TERT+):
  BRAF only            25 / 281   (8.9%)
  RAS only              6 /  54   (11.1%)
  Triple negative       5 / 178   (2.8%)

AJCC stage:
  Stage I               4
  Stage II              4
  Stage III            10
  Stage IVA             8
  Stage IVC             3
  Stage IV              1
  unknown / blank       6

  → 22 / 30 staged TERT+ are Stage III/IV (73%)

Differentiation (TDS16, lower = more dedifferentiated):
  TERT+ mean    6.238
  WT mean       6.879
  delta        -0.641
```

These distributions match published thyroid TERT biology (TERT+ enriched in
advanced stage, less differentiated, higher event rate).

---

## Survival analysis

```
TERT+ vs WT logrank:        p = 4.92e-06    test stat = 20.87
  TERT+ events:    6 / 36   (16.7%)
  WT events:      10 / 468  (2.1%)

4-group (BRAF / RAS / TERT+ / Triple-neg):
  logrank         p = 3.78e-05    test stat = 23.14
  Group sizes:    BRAF_only 250, RAS_only 48, TERT+ 36, Triple_neg 170
```

**Honest caveat:** absolute event count in TERT+ is 6 — small. Hazard
ratios will have wide 95% CIs. Recommended to use Cox with Firth correction
or bootstrap CIs, and present n + event count next to every HR.

---

## Cross-tabs

### v17 dark-matter (DM1 vs DM2) cluster

DM clustering only covers 178 of 513 patients (the dark-matter-positive
subcohort). Among them:

|  | mutated | wildtype |
|---|---:|---:|
| DM1 | 4 | 106 |
| DM2 | 1 | 68 |
| (non-DM) | 31 | 303 |

Fisher DM1 vs DM2: odds = 2.57, **p = 0.65** — under-powered.

The DM-cluster narrative cannot be salvaged from this — only 5 TERT+ fall in
the DM subcohort. The strong signal lives in the **whole-cohort 4-group
analysis**, not the DM split.

### Driver anchor

|  | mutated | wildtype |
|---|---:|---:|
| BRAF | 25 | 256 |
| RAS | 6 | 48 |
| NTRK | 1 | 0 |
| TP53 | 0 | 1 |
| unknown | 4 | 172 |

---

## Honest limitations

1. **Provenance:** All 36 calls come from `thca_tcga_pub` (cBioPortal mirror
   of TCGA Cell 2014 publication MAF). These were Sanger-validated by the
   original publication; we did not re-call from BAMs. This is the same
   level of trust as any cBioPortal-derived analysis — but worth noting
   plainly.

2. **Liu 2017 *JCO* supplement (PMID 27979994)** — potentially **100+
   additional Sanger-validated TERT calls** for TCGA-THCA, but all 8
   candidate URLs returned 0-byte / 404 / non-data. This is the single
   highest-leverage outstanding source. Manual library / interlibrary loan
   or author outreach (Liu R, Xing M) would unlock it.

3. **Pozdeyev 2024 *Nature Cancer*** — not located via DOI guesses. The
   paper claims 1,408 PTC + 423 PDTC + 248 ATC with TERT assayed — recovery
   would substantially extend external validation.

4. **Yoo 2019 *Nat Commun* Korean cohort** — Springer ESM URLs returned
   non-data bodies. Manual download is straightforward; would unlock a
   non-TCGA cohort cross-validation.

5. **DM cluster Fisher is underpowered** — only 5 TERT+ in DM subcohort.
   Cannot make the "TERT enriches in DM1" claim from this data.

6. **External cohorts (MSK 1,844 records) cannot be patient-matched** to
   TCGA — they validate prevalence at the cohort level, not at the
   per-patient signature level.

7. **The recovery is reproducible audit-wise:** every URL attempted is
   logged in `logs/*_attempts.json` with status, bytes, timing — closes
   the "why didn't you check TERT?" reviewer attack vector regardless of
   outcome.

---

## Paper-impact judgment (informal)

- **npj Precision Oncology (current target):** ~70% → ~78%
  - 4-group narrative now defensible
  - Survival signal strong
  - Single-cohort caveat remains (no patient-matched external validation)
- **Genome Medicine:** ~40% → ~50%
  - Meaningful but not uniquely positioning
- **Reviewer attack vector closed**, regardless of venue, by the 9-source
  audit log

---

## Question for strategic judgment

1. **Worth re-running the survival figure as a main-text panel** of the
   v17 paper, or keep it supplementary? Pro: p = 4.9e-6 is striking and
   the 4-group split is clean. Con: 6 events is small and HR CIs will be
   wide; reviewers may push back on event count.

2. **How aggressively to chase the 3 outstanding sources** (Liu 2017,
   Pozdeyev 2024, Yoo 2019) before submission?
   - Manual outreach to Liu R / Xing M for the JCO supplement is the
     highest-EV next move (~hours of effort, +100 patients possible).
   - Pozdeyev / Fagin email to MSK for 2024 supplement — moderate effort,
     uncertain reply.
   - Yoo SK email or KOBIC request for Korean cohort — moderate effort,
     opens a Korean-cohort cross-validation that aligns with the
     Bundang/SNU narrative.

3. **DM1/DM2 narrative**: the 5-patient DM subcohort cross-tab is too small.
   Should the paper drop the "TERT enriches in DM1" claim entirely and
   pivot to "TERT defines a separate aggressive axis orthogonal to DM
   clustering"? The latter is supported by the data; the former isn't.

4. **External validation framing:** is it acceptable to report MSK-IMPACT
   prevalence (1,361 TERT+ in advanced/aggressive thyroid) as a "cohort-
   level prevalence consistency check" without patient matching? Or does
   that read as window dressing?

---

## Files (in `project/results/v17_tert_recovery/v2/`)

| File | Contents |
|---|---|
| `FINAL_recovery_audit.md` | Full per-source audit |
| `FINAL_paper_implication.md` | Concrete next actions |
| `FINAL_tert_status_integrated.tsv` | 36-row mutated patient list |
| `FINAL_promoter_records_all_sources.tsv` | 1,880 raw records (incl. external) |
| `sample_master_v17_tert_v2.tsv` | 513-row sample master with `tert_promoter_integrated` |
| `FINAL_crosstab_*.tsv` | Per-axis cross-tabs |
| `FINAL_extended_summary.json` | Survival + Fisher in JSON |
| `parsed/S6_cbioportal_all_studies_summary.tsv` | Per-study TERT counts |
| `parsed/S6_cbioportal_all_promoter_mutations.tsv` | All 1,820 promoter calls |
| `S5_cosmic_application_guide.md` | COSMIC academic-access workflow |
| `S7_gdc_controlled_application_guide.md` | dbGaP TCGA-THCA WGS workflow |
| `logs/*_attempts.json` | Per-URL audit (bytes, status, timing) |
| `orchestrator_summary.json` | Orchestrator run output |
