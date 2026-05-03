# Paper 3 ICI — Dataset Registry

**Working title:** HLA loss, neoantigen architecture, and immune ecotypes define ICI vulnerability in molecularly dark thyroid cancer
**Track:** A — design only. No download, no analysis. All accessions marked `to_verify` until Track B kickoff.
**Author:** Seungho Cook
**Status date:** 2026-05-04

**Claim guard:** Datasets are catalogued for ICI vulnerability / immunogenomic prioritization, never for "ICI response prediction" in thyroid unless thyroid ICI-treated raw RNA-seq is acquired.

---

## 0. Verification policy

For each accession:
- **`to_verify`** — accession string is from prior literature / user input; existence, version, file types, and access conditions must be checked at Track B kickoff via GEO/SRA/EGA/dbGaP web UI (not by Claude during Track A).
- **`internal`** — already on local mirror or processed for Paper 1/2; reuse read-only with paper-boundary discipline.
- **`controlled`** — known to be controlled-access (dbGaP / EGA / institutional). Application timing must be added to the 12-week plan.
- **`tier_x`** — see priority tiers below.

Priority tiers:
- **Tier 1** — required for headline figure. Track B blocks if unavailable.
- **Tier 2** — replication / cross-platform. Track B proceeds if at least 50% of Tier 2 are usable.
- **Tier 3** — supplementary / sensitivity. Loss tolerable.

---

## 1. Bulk RNA-seq ± WES — discovery & dedifferentiation cohorts

| Accession | Modality | n (claimed) | Subtype | Role | Tier | Status | Notes |
|---|---|---|---|---|---|---|---|
| TCGA-THCA | RNA-seq + WES + clinical | ~500 | PTC, FVPTC, TCV; ATC limited | Primary discovery: ecotype NMF, HLA LOH, neoantigen, dark-matter stratification | 1 | internal | Driver/fusion calls already curated in Paper 1 ETL. WES BAM available via GDC; LOHHLA needs raw BAM not just MAF. |
| GSE76039 | RNA-seq | ~37 (PDTC + ATC; Landa 2016) | PDTC, ATC | Dedifferentiation axis anchor | 1 | to_verify | Landa 2016 JCI (cite already saved in `v17_landa2016_cite_save`). Small n — pair with TCGA aggressive subset for power. |
| GSE33630 | microarray (Affy U133 Plus 2) | ~105 | PTC, ATC, normal | Cross-platform replication | 2 | to_verify | Microarray, so signature scoring only — no neoantigen, no HLA LOH. |
| GSE65144 | microarray | ~13 | ATC + normal | ATC enrichment | 2 | to_verify | Tiny n; suitable as ATC-vs-normal contrast in figure inset. |
| GSE54958 | microarray | small | PTC stages | Stage-axis sensitivity | 3 | to_verify | |
| GSE53157 | microarray | ~27 | PTC | Replication | 3 | to_verify | |
| GSE29265 | microarray | ~49 | PTC, ATC, normal | Cross-platform | 2 | to_verify | |
| Yoo SK et al. 2019 *Nat Commun* | RNA-seq | Korean ATC + advanced DTC | ATC, advanced DTC | Asian dedifferentiation generalization | 1 | to_verify | Accession not yet pinned; check supp; may be EGA-controlled. If controlled, add EGA application to Wk 1. |
| PRJKA210106 | RNA-seq | n=282 Korean PTC | PTC | Asian baseline | 2 | to_verify | Cross-paper care — already touched by Paper 1/2 lineage. Use for HLA-typing reference, not for ecotype discovery, to keep Paper 3 distinct. |
| PRJEB11591 (Yoo 2016 SNU-GMI) | RNA-seq | n=260 Korean PTC | PTC | HLA-typing reference (arcasHLA already done) | 3 | internal (read-only) | Per `v17_arcasHLA_korean_k2`. Do not re-run; do not duplicate Paper 1/2 analyses. Use only as Korean HLA-frequency anchor. |

---

## 2. Single-cell RNA-seq — atlas construction

| Accession | Platform | n samples (claimed) | Subtype focus | Role | Tier | Status | Notes |
|---|---|---|---|---|---|---|---|
| GSE184362 | 10x 3' | thyroid tumor + normal | PTC + ATC | Primary atlas — epithelial + immune + stromal | 1 | to_verify | Per literature, ATC/PDTC content. Verify n PDTC/ATC vs PTC. |
| GSE193581 | 10x | PTC | PTC | T cell / myeloid sub-states | 1 | to_verify | |
| GSE232237 | 10x | thyroid | mixed | Atlas integration | 2 | to_verify | |
| GSE191288 | 10x | thyroid | PTC | Replication | 2 | to_verify | |
| GSE148673 | 10x | thyroid | mixed | Replication | 3 | to_verify | |
| Han et al. 2024 *JCI Insight* | 10x | thyroid | dedifferentiation | Atlas | 2 | to_verify | Accession to be pulled from supp at Track B kickoff. |

**Atlas-side risks:** ATC/PDTC scRNA samples are scarce. Combined ATC + PDTC scRNA cells across all six accessions are likely <5K — viable for state annotation by reference projection but not for de novo trajectory inference. The atlas plan in `paper3_ici_scRNA_reference_atlas_summary.md` accounts for this.

---

## 3. Spatial transcriptomics — TLS / niche localization

| Accession | Platform | n (claimed) | Subtype | Role | Tier | Status | Notes |
|---|---|---|---|---|---|---|---|
| GSE250521 | Visium | DM1-spatial validation set | PTC + DM1 | TLS / B-cell niche localization (already used in Paper 1/2 supplement) | 2 | internal (read-only) | Paper-boundary care. Paper 3 may reference niche maps but must NOT duplicate Paper 1's DM1 spatial analysis. |
| GSE230424 | Visium | thyroid | mixed | Spatial replication | 3 | to_verify | |
| GSE248205 | Visium | thyroid | mixed | Spatial replication | 3 | to_verify | |
| Ning / Liao / Zheng spatial | Visium / VisiumHD | thyroid | mixed | Asian spatial | 3 | to_verify | Accessions to be located via PubMed → supp at Track B kickoff. |

---

## 4. Pan-cancer ICI reference — DIAL audit + transfer

ICI-treated raw RNA-seq cohorts used for direction-invariance testing of the signature suite. Thyroid ICI-treated raw RNA-seq is *not* expected to be available, hence the transfer-only framing.

| Accession | Disease | n (claimed) | ICI agent | Role | Tier | Status | Notes |
|---|---|---|---|---|---|---|---|
| IMvigor210 (EGAS00001002556 / SRP172670) | UC | ~298 | atezolizumab | DIAL anchor 1 | 1 | controlled / partial-public | Pre-treatment RNA-seq + RECIST. Public package available via `IMvigor210CoreBiologies` R data drop. |
| Hugo GSE78220 | melanoma | ~28 | anti-PD-1 | DIAL anchor 2 | 1 | to_verify | Small but heavily reused in pan-cancer ICI papers; pre-treatment. |
| Riaz GSE91061 | melanoma | ~109 | nivolumab | DIAL anchor 3 | 1 | to_verify | Pre + on-treatment; use pre only for DIAL. |
| Gide melanoma (PRJEB23709) | melanoma | ~73 | anti-PD-1 ± anti-CTLA-4 | DIAL anchor 4 | 2 | to_verify | |
| Liu melanoma (dbGaP phs000452) | melanoma | ~121 | anti-PD-1 | DIAL anchor 5 | 2 | controlled | dbGaP application required. If timeline blocks, drop to Tier 3. |
| Kim gastric (PRJEB25780) | GC | ~45 | pembrolizumab | DIAL extension | 3 | to_verify | Pan-cancer breadth. |
| Cho NSCLC | NSCLC | varies | anti-PD-1 / anti-PD-L1 | DIAL extension | 3 | to_verify | Specific accession to be pinned at Track B kickoff. |

**Thyroid ICI gap (key risk).** No ICI-treated thyroid raw RNA-seq cohort is reliably public as of design date. Implications:
- Headline claim must be "ICI vulnerability / readiness", not "response predictor".
- Validation strategy = pan-cancer transfer (DIAL audit identifies which signatures retain direction in non-thyroid → apply only those to thyroid).
- A single-arm thyroid ICI report (small, e.g. lenvatinib + pembrolizumab in ATC) may exist as supplementary tables only; treat as case-level anecdote, not validation.

---

## 5. Internal / Paper-1 / Paper-2 reuse — boundary discipline

Reusable read-only artifacts (no re-analysis allowed during Track B unless explicitly opened):

| Path | Origin | Allowed reuse |
|---|---|---|
| `project/results/00_qc/` | Paper 1 ETL | Sample-level QC tables — read-only reference. |
| `project/results/01_spatial_score/` | Paper 1/2 spatial pipeline | Spatial niche maps for TLS visualization only. Do not recompute. |
| `project/results/03_pathology_poc/` | Paper 2 H&E projection | Pathology features as comparator — *not* as Paper 3 input. |
| PRJEB11591 arcasHLA results | per `v17_arcasHLA_korean_k2` | Korean HLA frequency anchor. Do not re-run. |
| TCGA-THCA driver/fusion calls | Paper 1 ETL | Dark-matter (BRAF/RAS-negative, fusion-negative) stratification. |

**Boundary rules (binding for Paper 3 prose):**
- Paper 1 (DM1 molecular) — Paper 3 may *cite* the dark-matter definition. May not duplicate DM1 cluster analyses.
- Paper 2 (H&E-DM1 / HT-overlap PTC) — Paper 3 may *cite* the immune-rich PTC+HT phenotype as one comparator subgroup. May not duplicate TLS/AICDA/BCR analyses.
- Paper 4 (GD HLA backlog) — Paper 3 is cancer-side; reference at most at the population HLA-frequency level.

---

## 6. Verification checklist (to execute at Track B kickoff, NOT now)

1. For each `to_verify` accession: open GEO/SRA/EGA/ENA UI; confirm (a) records exist, (b) raw FASTQ or processed counts available, (c) clinical/sample metadata accessible, (d) license/access mode.
2. For each `controlled` accession: identify required application (dbGaP, EGA, institutional), determine 4–8 week lag, place in 12-week plan accordingly.
3. Resolve Yoo SK 2019 Korean ATC accession; if EGA-controlled, decide whether to apply or drop to Tier 2 with TCGA aggressive subset substitution.
4. Resolve Han 2024 JCI Insight accession from paper supp.
5. Resolve Ning / Liao / Zheng spatial accessions.
6. Resolve Cho NSCLC ICI accession.

Verification must be logged into `paper3_ici_dataset_registry_VERIFIED.tsv` (TSV form) at kickoff.

---

## 7. Out-of-scope datasets (explicit exclusions)

- ICI-treated cohorts in non-relevant diseases (e.g., RCC, CRC) unless DIAL audit specifically benefits — added only via amendment.
- Non-thyroid scRNA atlases used as reference (e.g., pan-cancer T-cell atlases like Zheng 2021, Andreatta 2022) — these are *reference annotations*, not data inputs, and will be cited as method dependencies in `paper3_ici_scRNA_reference_atlas_summary.md`.
- Animal models — out of scope for this paper.
- Patient-derived organoids without paired tumor RNA-seq — out of scope.

---

Track A completed. No marathon violation.
