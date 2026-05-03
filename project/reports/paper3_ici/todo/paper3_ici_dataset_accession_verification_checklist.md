# Paper 3 ICI — Dataset Accession Verification Checklist

**Status:** Verification = web-UI lookup only (GEO / SRA / EGA / ENA / ArrayExpress portals). Marathon-compatible because no download / no compute. **Defer all `verify` actions to Track B kickoff** unless an item is needed to unblock controlled-access application drafting (Section B / C of `paper3_ici_controlled_access_checklist.md`).
**Author:** Seungho Cook
**Created:** 2026-05-04

**Output target:** `paper3_ici_dataset_registry_VERIFIED.tsv` produced at Track B Wk 1.

---

## A. Bulk RNA-seq / WES — discovery cohorts

For each, verify (a) accession resolves, (b) raw FASTQ or processed counts available, (c) clinical/sample metadata accessible, (d) license/access mode, (e) sample count matches `paper3_ici_dataset_registry.md` claim.

- [ ] **TCGA-THCA** — internal mirror confirmed; reuse Paper 1 ETL
- [ ] **GSE76039** — Landa 2016 PDTC/ATC; verify n=37 RNA-seq + paired WES per-sample availability
- [ ] **GSE33630** — verify Affy U133 Plus 2 platform, n≈105
- [ ] **GSE65144** — verify n≈13 ATC + normal
- [ ] **GSE54958** — verify PTC stages cohort
- [ ] **GSE53157** — verify n≈27 PTC
- [ ] **GSE29265** — verify n≈49 PTC + ATC + normal
- [ ] **Yoo SK 2019 Nat Commun Korean ATC/aDTC** — pin accession (likely EGA); verify access mode
- [ ] **PRJKA210106 Korean PTC n=282** — verify ENA presence + access mode + paper boundary care (already touched by Paper 1/2 lineage)
- [ ] **PRJEB11591 Yoo 2016 SNU-GMI** — internal arcasHLA results confirmed (`v17_arcasHLA_korean_k2`); read-only reference

---

## B. Single-cell RNA-seq — atlas

- [ ] **GSE184362** — verify thyroid PTC + ATC scRNA presence; confirm 10x platform; n samples
- [ ] **GSE193581** — verify PTC scRNA; n samples
- [ ] **GSE232237** — verify; n samples
- [ ] **GSE191288** — verify; n samples
- [ ] **GSE148673** — verify; n samples
- [ ] **Han 2024 JCI Insight** — pull accession from supp; verify

For each cohort: record cell count (post-QC estimate), tumor subtype distribution (PTC vs PDTC vs ATC vs normal), 10x version (v2 / v3 / v3.1).

---

## C. Spatial transcriptomics

- [ ] **GSE250521** — internal (Paper 1/2 spatial); confirm read-only access
- [ ] **GSE230424** — verify Visium presence; n samples
- [ ] **GSE248205** — verify; n samples
- [ ] **Ning / Liao / Zheng spatial** — locate accessions via PubMed; verify

---

## D. Pan-cancer ICI reference (DIAL)

For each, verify accession + RECIST / response endpoint metadata + pre-treatment sample availability.

- [ ] **IMvigor210** — `IMvigor210CoreBiologies` R-package access confirmed; n≈298
- [ ] **Hugo GSE78220** — verify n≈28 melanoma anti-PD-1 pre-treatment
- [ ] **Riaz GSE91061** — verify n≈109; isolate pre-treatment subset
- [ ] **Gide PRJEB23709** — verify n≈73 melanoma anti-PD-1 ± anti-CTLA-4
- [ ] **Liu melanoma phs000452** — verify dbGaP accession current
- [ ] **Kim gastric PRJEB25780** — verify n≈45 pembrolizumab
- [ ] **Cho NSCLC** — pin accession from publication; verify

---

## E. Per-accession verification record (fill at Track B Wk 1)

| Accession | Resolves? | Raw available? | Metadata? | Access mode | n confirmed | Verified date | Verifier |
|---|---|---|---|---|---|---|---|
| TCGA-THCA | yes | yes | yes | open + controlled | ~500 | 2026-05-04 | internal | 
| GSE76039 | — | — | — | — | — | — | — |
| GSE33630 | — | — | — | — | — | — | — |
| GSE65144 | — | — | — | — | — | — | — |
| GSE54958 | — | — | — | — | — | — | — |
| GSE53157 | — | — | — | — | — | — | — |
| GSE29265 | — | — | — | — | — | — | — |
| Yoo 2019 ATC | — | — | — | — | — | — | — |
| PRJKA210106 | — | — | — | — | — | — | — |
| PRJEB11591 | yes | n/a (internal) | yes | open | 260 | 2026-04-29 | internal |
| GSE184362 | — | — | — | — | — | — | — |
| GSE193581 | — | — | — | — | — | — | — |
| GSE232237 | — | — | — | — | — | — | — |
| GSE191288 | — | — | — | — | — | — | — |
| GSE148673 | — | — | — | — | — | — | — |
| Han 2024 | — | — | — | — | — | — | — |
| GSE250521 | yes | n/a (internal) | yes | open | — | 2026-05-03 | internal |
| GSE230424 | — | — | — | — | — | — | — |
| GSE248205 | — | — | — | — | — | — | — |
| Ning/Liao/Zheng spatial | — | — | — | — | — | — | — |
| IMvigor210 | — | — | — | open R-pkg | ~298 | — | — |
| GSE78220 | — | — | — | — | ~28 | — | — |
| GSE91061 | — | — | — | — | ~109 | — | — |
| PRJEB23709 (Gide) | — | — | — | — | ~73 | — | — |
| Liu phs000452 | — | — | — | controlled | ~121 | — | — |
| PRJEB25780 (Kim) | — | — | — | — | ~45 | — | — |
| Cho NSCLC | — | — | — | — | — | — | — |

---

## F. Verification rules

- Verification is **portal lookup only**. No `wget` / `prefetch` / `gdc-client` / dbGaP-token download during verification.
- If an accession does not resolve → flag in registry as `MISSING`, not silently dropped.
- If an accession resolves but n disagrees with claim by > 20% → flag as `N_DISCREPANCY`, document, do not silently proceed.
- If a paper-supp lookup is needed: read paper supp via existing PDF tooling — do not download supplementary tables to disk during marathon.

---

## G. Verification cadence

- **During marathon (now → Paper 1 ship):** verification work blocked unless trivially fast (<5 min per accession) AND not displacing Paper 1/2 writing. Default: defer.
- **At Track B Wk 1:** full sweep across all unchecked items; produce `paper3_ici_dataset_registry_VERIFIED.tsv`.
- **At any Paper 3 prose drafting step:** verify accessions cited, even before Track B compute.

---

## H. Pre-Track-B output

When all items checked → produce `paper3_ici_dataset_registry_VERIFIED.tsv` with columns:
`accession, source_portal, n_confirmed, modality, access_mode, license, verified_date, verifier, notes`

Drop into `project/reports/paper3_ici/` (separate from frozen design bundle).

---

Paper 3 Track A frozen. Return to Paper 1/2 marathon.
