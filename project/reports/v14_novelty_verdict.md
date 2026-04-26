# v14 Novelty Verdict — BRAF × TROP2 × Thyroid

_Date: 2026-04-25 · Author: Seungho Cook · 30-min reality check before paper writing._
_Sources: PubMed E-utilities (12 prior-art queries, 94 records; +
new thyroid-TROP2 temporal sweep, 38 records 1981–2026),
ClinicalTrials.gov v2 (NCT06235216, NCT07521670), v14 internal results
(CCLE 13-line + PRISM 19Q4 — both independent of v13/DIAL evidence,
which v5.2 LODO finding ruled out as leak artifact)._

---

## A. Timing analysis

### Thyroid-TROP2 publication rate (PubMed, query `(TROP2 OR TACSTD2) AND thyroid`)

| Bucket | n | Rate (papers/yr) |
|---|---:|---:|
| Pre-2010          | 0  | 0   |
| 2010–2014         | 0  | 0   |
| 2015–2019         | 12 | 2.4 |
| 2020–2023         | 12 | 3.0 |
| **2024–2026**     | **14** | **4.7** |
| Total             | **38** | — |

Per-year drilldown of the recent window: 2024 → 5, 2025 → 6, 2026 → 3 (YTD).

> **"Underexplored → emerging" 정확한가?**
> **Half-true.** The field has had a sustained ~3 papers/yr since 2017,
> ramping ~1.5× to 4.7/yr in 2024–2026. By the rate-acceleration rule
> (>1.8× ⇒ BURST, 0.8–1.8× ⇒ SUSTAINED), this is **SUSTAINED with mild
> acceleration**, not "newly hot." Re-frame the v14 narrative from
> *"underexplored → emerging"* to **"steadily warm, accelerating."**

### Where the niche is genuinely thin

- **BRAF × TROP2 × thyroid intersection: only 2 of 38 papers** mention BRAF in title/abstract.
  - PMID 31949805 (Liu 2018, *Int J Clin Exp Pathol*) — IHC: BRAF V600E ↔ TROP2-high in PTC.
  - PMID 38696857 (2024, miRNA paper) — incidental BRAF mention, not a TROP2-stratification study.
  - PMID 29228520 (Bychkov 2018) is at the same intersection but did not surface in the title/abstract sweep — captured in `braf_trop2_thyroid_priorart.tsv` query A.
- **Cell-line transcriptomic restatement (CCLE/PRISM): 0 prior reports.**
- **Two-gate stratification (BRAF × TACSTD2): 0 prior reports.**

### Trial timing (the actionable window)

| Trial | NCT | Phase | Status | n | BRAF in eligibility? | TROP2 IHC required? |
|---|---|---|---|---:|---|---|
| SETHY (sacituzumab govitecan, DTC + ATC) | NCT06235216 | 2 | RECRUITING (since 2024-09) | 42 | **No** (0 mentions) | **No** (0 mentions) |
| STRAP (sacituzumab tirumotecan, PTC + ACC) | NCT07521670 | 2 | NOT_YET_RECRUITING | 68 | **No** (0 mentions) | **No** — explicitly waived in inclusion criterion §2 |

Primary completion 2027-12 (SETHY) / 2028-10 (STRAP). **Both trials run
BRAF-blind and TROP2-blind**, so a 2026-mid-year publication is
pre-emptively informative for retrospective stratification.

> **Timing assessment: *adequate window, not golden.*** With the rate
> climbing, the BRAF × TROP2 thyroid niche could be claimed by another
> group in 18–24 months. Move quickly; do not over-claim novelty.

---

## B. Genuinely novel in v14 (post-v5.2)

1. **First cell-line transcriptomic corroboration of the BRAF-axis TROP2
   up-regulation in PTC.** The 2018 IHC papers (Liu PMID 31949805,
   Bychkov PMID 29228520) established the link at the protein level in
   primary tumours. Our v14 result restates it at the **CCLE 13-line
   transcriptomic level** (top-3 TACSTD2 lines all V600E). N is small
   but orthogonal to IHC, and necessary to enable downstream cell-line
   pharmacology. — *Lit support: strong (re-confirmation).*
   - **Strengthening pass 2026-04-25 (Tier 1):**
     - **GTEx baseline** — normal thyroid TACSTD2 = mid-pack rank 16/54
       tissues (8.10 TPM); TCGA-THCA tumour vs normal Δlog2 = +3.26
       (~9.55× FC), Cohen's d = 1.50, Welch p ≈ 0. **TROP2 high in PTC
       is a tumour-acquired phenotype, not a thyroid-tissue baseline.**
     - **4th-cohort GEO replication (GSE58545, n=45, Bauer 2014)** — never
       used in v8/v8p1/v13/v14. BRAF V600E PTC mean log2 TACSTD2 = 10.98
       vs WT 9.35 vs normal 5.05; PTC-internal z = +0.326 vs −0.652.
       Mann-Whitney p = 0.0020 (significant), Welch p = 0.086 (borderline,
       small WT arm n=9). Direction matches; **partial replication**.

2. **First thyroid-specific PRISM read-out of BRAF-like SN-38-class
   selectivity.** No PubMed-indexed record (52 records reviewed, queries
   G/H/T2c/T2d) reports topotecan/irinotecan ΔLFC stratified by BRAF in
   thyroid lines. The closest neighbour PMID 34560251 is melanoma TOP1
   in BRAF-i resistance — different tumour, different direction.
   — *Lit support: moderate (concept adjacent in melanoma; first
   specific thyroid report).*
   - **Strengthening pass 2026-04-25 (Tier 1):**
     - **PRISM 24Q2 independent REP.1M screen — confirmed.**
       *belotecan* (new topo-I poison only in 24Q2, not in 19Q4) shows
       ΔLFC = **−1.63** under BRS52 grouping. **9 of 9 topo-I poisons**
       in 24Q2 panel show negative ΔLFC (range −0.71 to −2.19). Class
       signal is coherent across two independent screens.
     - **24Q2 caveat — must disclose in paper.** The 24Q2
       `Extended_Primary_Data_Matrix` re-uses 19Q4 LFC verbatim (4 d.p.
       match) for compounds not re-screened (topotecan, irinotecan,
       SN-38, exatecan, simvastatin, kaempferol). Genuine independent
       replication is limited to compounds re-screened in REP.1M
       (belotecan adds; atorvastatin null below).
     - **Grouping caveat re-confirmed.** Effect requires BRS52
       transcriptomic grouping (`v14_prediction`); under raw mutation
       grouping topotecan/irinotecan ΔLFC flips positive — already
       acknowledged in v14, still applies.

3. **Pre-specified, stratification-ready biomarker hypothesis (BRAF-like ×
   TACSTD2-high) for two BRAF-unselected sacituzumab thyroid trials.**
   The proposal — not the underlying biology — is novel; both trials
   verified BRAF-blind and TROP2-blind above.
   — *Lit support: strong (trial-design facts unambiguous).*

4. **Mechanistic coherence of the BRAF-selective PRISM hit list (with
   one weakening).** Of the 5 BRAF-selective PRISM compounds, 3 land on
   independently established axes — topotecan (TROP2-ADC payload class;
   Liu 2018, Bychkov 2018, TROPiCS-02), simvastatin / atorvastatin
   (LDLR / cholesterol axis, well-precedented in cancer pharmacology) —
   and a 4th (kaempferol) is a flavonoid with documented CYP1B1 activity
   in the literature. The coherence is read out from the PRISM screen
   alone (subtype-stratified ΔLFC), independent of any v13/DIAL-derived
   prior. — *Lit support: moderate (each axis has external citations).*
   - **Strengthening pass 2026-04-25 weakens the LDLR axis specifically.**
     Atorvastatin re-screened in PRISM REP.1M (24Q2 independent screen)
     gives ΔLFC = **−0.038** (vs **−0.645** in 19Q4) — does NOT replicate.
     Rosuvastatin drifts +0.004 → +0.40 (small magnitude, RAS-selective
     direction). **Treat the LDLR / statin axis as hypothesis-level
     only**; do not headline it in the paper. Topotecan / SN-38 /
     belotecan class signal is unaffected.

## C. Already established (must cite)

| Established claim | Citation | Year |
|---|---|---|
| BRAF V600E ↔ TROP2 over-expression in PTC (IHC, primary tumour) | Liu PMID 31949805; Bychkov PMID 29228520 | 2018 |
| TROP-2 / 5hmC in follicular-patterned thyroid neoplasms | PMID 35121239 | 2022 |
| TROP2 as ADC target in anaplastic thyroid carcinoma | PMID 35158847; PMID 33289434 | 2021–2022 |
| Pan-tumour TROP2 IHC landscape (n=18,563, incl. thyroid) | PMID 35477165 | 2022 |
| Sacituzumab govitecan clinical activity (TNBC, urothelial; TROPiCS-02) | PMID 39067901 | 2024 |
| BRAF-mut mCRC patients exposed to irinotecan as standard chemo backbone | PMID 41852303, 41505697 | 2024–2026 |
| TOP1 modulation in melanoma BRAF-i resistance | PMID 34560251 | 2021 |

The two 2018 IHC papers **must appear in the abstract or Introduction**
— failing to cite them while claiming any BRAF→TROP2 link in PTC would
fail review on novelty grounds.

## D. Recommended framing

### Title options (ranked by realism)

| | Title | Posture |
|---|---|---|
| **A. 보수적 (recommended fallback)** | *Cell-line and pharmacogenomic corroboration of BRAF-driven TROP2 expression in papillary thyroid cancer: orthogonal evidence for stratified sacituzumab-govitecan use.* | Confirmatory + applied |
| **B. 중도적 (recommended primary)** | *A BRAF-like × TROP2-high two-gate stratification hypothesis for TROP2-directed antibody–drug conjugate trials in papillary thyroid cancer.* | Translational |
| **C. 야심적 (avoid)** | *BRAF-axis primary sensitivity to topoisomerase-I poisons in papillary thyroid cancer.* | Discovery — would need xenograft we lack |

### Venue probabilities (acceptance odds)

| Venue | Title A | Title B | Title C |
|---|---:|---:|---:|
| **Bioinformatics (OUP)** | high | high | low |
| **JCO Precision Oncology** | high | high | medium |
| **Endocrine-Related Cancer** | medium | high | low |
| **Cancer Medicine** | high | medium | low |
| **Cancers / Frontiers in Oncology** | high (fallback) | high (fallback) | medium |
| **Clinical Cancer Research** | low | medium | low |
| **Nature Communications** | n/a | n/a | n/a (no wet-lab) |
| **Cancer Discovery / Cancer Cell** | n/a | n/a | n/a |

> **Recommendation:** **Title B → JCO Precision Oncology primary,**
> Bioinformatics secondary. Title A as graceful retreat if reviewers
> challenge the trial-stratification framing. Skip Title C and any
> top-tier venue without wet-lab follow-up.

## E. Strengthening priorities

### Tier 1 — computational, weeks (before submission)

- Add explicit citation of **Liu 2018 (PMID 31949805) and Bychkov 2018
  (PMID 29228520)** in the v14 §5.Y opening sentence. Frame the v14
  result as **"CCLE/PRISM observations consistent with the 2018 IHC
  link"** — anchored on the IHC priors directly, NOT on internal
  v13/DIAL inheritance (v5.2 LODO finding rules out citing v13 DIAL
  evidence downstream).
- Pull **TROPiCS-02 sacituzumab efficacy numbers** (PMID 39067901) into
  the discussion for clinical context.
- ~~**Independent PRISM Repurposing 23Q4 / 24Q2 re-analysis**~~
  **DONE 2026-04-25** — 24Q2 downloaded (Figshare 25917643). Class
  signal confirmed via independent REP.1M (belotecan ΔLFC = −1.63;
  9/9 topo-I poisons negative). Two caveats surfaced and must appear
  in paper: (a) 24Q2 reuses 19Q4 LFC for non-rescreened compounds;
  (b) atorvastatin null re-screen (LDLR axis demoted to hypothesis-
  level). See `reports/v14_strengthening_prism_newer.md`.
- ~~**GTEx normal-thyroid TROP2 baseline**~~ **DONE 2026-04-25** —
  GTEx v8 median TPM downloaded; thyroid TACSTD2 = mid-pack rank
  16/54 (8.10 TPM); TCGA tumour vs normal Δlog2 = +3.26, Cohen's d
  = 1.50, p ≈ 0. **Tumour-acquired, not baseline artefact.** See
  `reports/v14_strengthening_gtex_baseline.md`.
- ~~**External validation in ICGC-THCA / CPTAC-THCA**~~
  **DONE 2026-04-25, with negative result on the named cohorts.**
  ICGC THCA-CN/SA expose mutations only (no expression); ICGC THCA-US
  is the same TCGA cohort already used. CPTAC has no thyroid study at
  all. **Pivot:** GSE58545 (n=45, Bauer 2014, never used in v8/v8p1/
  v13/v14) selected as 4th cohort. BRAF V600E PTC TACSTD2: Mann-Whitney
  p = 0.0020 vs WT, Welch p = 0.086 (small WT arm n=9). Direction
  matches; **partial replication**. See
  `reports/v14_strengthening_external_cohort.md`.
- ~~**bioRxiv / medRxiv / Europe PMC preprint scan**~~ **DONE 2026-04-25**
  — 10 queries × Europe PMC `SRC:PPR`, total 5 preprint hits, all
  triaged as **no conflict** (lung NET, lung fibroblast senolytics,
  pan-cancer MPZL1, CAR T, cfDNA). See
  `results/v14_priorart/preprint_scoop_scan.md`. **Re-run within 48 h
  of submission** (preprint indexing lag 1–14 days).
- ~~**Patent search** for "BRAF-stratified TROP2 ADC thyroid"~~
  **DONE 2026-04-25** via WebSearch + WebFetch on top candidates. **No
  patent claims BRAF-stratified TROP2 ADC use in thyroid.** Closest
  neighbours: **US20250152731A1** (UT System, 2021, TROP2-low + DNMT-i,
  thyroid covered but BRAF-blind — *complementary*, not blocking) and
  **JP7525633B2** (Immunomedics, 2020, sacituzumab govitecan biomarkers
  = DDR genes; **no thyroid, no BRAF**). Composition-of-matter held by
  Daiichi Sankyo / Immunomedics is normal commercial context. **Add one
  Discussion sentence citing US20250152731A1 as orthogonal stratification.**
  See `results/v14_priorart/patent_landscape_findings.md`.
- ~~**ASCO / ESMO / AACR / ATA 2025–2026 abstracts**~~
  **DONE 2026-04-25** via WebSearch on each venue. **No conference
  scoop risk.** SETHY trial protocol (ASCO 2024 `TPS6130`) is the only
  thyroid-TROP2-ADC public abstract; BRAF-blind design re-confirmed.
  TROP2-ADC AACR/ESMO activity is dense but breast/lung-focused.
  Citations to add (Discussion clinical-context paragraph): 2025 ATA
  DTC Guidelines, ASCO Thyroid Cancer Guideline (`JCO-26-00235`),
  ESMO 2024 dabrafenib + trametinib BRAF V600E DTC real-world,
  AACR ADC-biomarker review (`CCR 32:4:661`). Final 5-min manual
  check on ATA 2025 abstract title list morning of submission. See
  `results/v14_priorart/conference_abstract_scan.md`.
- Add a caveat box noting only 2 of 38 thyroid-TROP2 papers mention
  BRAF — i.e. the intersection is genuinely thin.

### Tier 2 — collaboration required, months

- **Retrospective TROP2 IHC + BRAF genotype overlay** on an institutional
  PTC cohort with sacituzumab-eligible (RAI-refractory) patients. The
  clinical co-author with archival pathology access is the bottleneck.
  Highest-leverage Tier-2 ask; converts the verdict from "computational
  corroboration" to "computational + IHC translational evidence."
- **NCT06235216 / NCT07521670 PI contact** for a pre-specified
  retrospective biomarker subgroup analysis collaboration (SETHY's
  archival-tissue requirement makes this concretely tractable).
- **Single-institution TMA** (~50–80 PTC) with TROP2 IHC + BRAF V600E
  genotyping and cross-cohort prediction.
- **scRNA-seq subtyping** to check whether the BRAF→TROP2 link is
  driven by a specific subpopulation within tumours.

### Tier 3 — wet-lab, year+

- **Sacituzumab govitecan vs free topotecan in BRAF-like (BCPAP, 8505C)
  and RAS-like (FTC-133, KTC-1) PTC xenografts** — directly answers the
  payload-vs-delivery question. Required to justify Title C.
- **MEK-i + statin combination test** in paired BRAF-like / RAS-like
  lines to disentangle the LDLR-axis observation from a generic
  BRAF-pathway metabolic vulnerability.
- **Investigator-initiated prospective trial** of sacituzumab govitecan
  in BRAF-V600E PTC (RAI-refractory, post-lenvatinib) with TROP2 IHC
  enrichment — converts the stratification proposal into a prospective
  biomarker-driven study.
- **Pan-cancer extension** to other BRAF-mutant tumours where TROP2 is
  already clinically measured (NSCLC, mCRC).

---

## Bottom line, four lines (post-strengthening 2026-04-25)

> **(1) Timing.** Thyroid TROP2 is **sustained, not bursting** —
> 3.0 → 4.7 papers/yr. The "underexplored → emerging" framing is
> half-correct; rephrase to "steadily warm, accelerating."
>
> **(2) Novelty.** Our **BRAF × TROP2 × thyroid intersection** is
> genuinely thin (2 of 38 papers). The **cell-line / PRISM
> restatement** and the **two-gate stratification proposal** are novel;
> the underlying biology (BRAF V600E ↔ TROP2-high in PTC, IHC) is
> already established (Liu 2018, Bychkov 2018) and must be cited.
>
> **(3) Strengthening evidence (Tier 1 closed).** GTEx confirms PTC
> TROP2 is tumour-acquired (Cohen's d=1.50, p≈0). 4th cohort GSE58545
> partially replicates BRAF↔TROP2 (Mann-Whitney p=0.0020). PRISM 24Q2
> belotecan REP.1M independently confirms the topo-I class signal
> (ΔLFC=−1.63), but atorvastatin null re-screen demotes the LDLR axis
> to hypothesis-level. ICGC/CPTAC THCA expression cohorts do not exist
> as available datasets — note this in the Methods. The PRISM 24Q2
> data-reuse caveat must be disclosed transparently. **Patent + ASCO/
> ESMO/AACR/ATA scoop scans closed: no IP block, no scoop. The
> BRAF-stratified TROP2-ADC thyroid niche is open.**
>
> **(4) Submission.** Submit **Title B to JCO Precision Oncology**
> within 4–6 weeks, anchor the Discussion in 2018 IHC priors + GTEx
> baseline + GSE58545 4th-cohort partial replication, hook the two
> BRAF-blind sacituzumab trials, headline belotecan as independent
> topo-I class confirmation, demote the statin axis to a hypothesis,
> and avoid Title C's causal language. **No "first to show" claims.
> Honest corroboration framing.**
