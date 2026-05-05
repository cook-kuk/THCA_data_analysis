# GPL570 thyroid dedifferentiation validation — access check
**Date:** 2026-05-04
**Scope:** Confirm public availability and processed-matrix usability of three GPL570 thyroid cancer datasets selected to extend the Landa 2016 / GSE76039 first-pass replication.
**Discipline:** processed Series Matrix only; no raw CEL; no other GEO datasets touched.

---

## 0. Allowed dataset list (per 2026-05-04 decision update)

```
GSE33630   GSE29265   GSE65144
```
No additional GEO acquisitions in this session.

---

## 1. Per-dataset access summary

### GSE33630 — Tomás et al. — Normal/PTC/ATC reference cohort

| Field | Value |
|---|---|
| Title | "Normal thyrocytes vs papillary vs anaplastic thyroid carcinomas" |
| Platform | GPL570 (Affymetrix HG-U133 Plus 2.0) — all samples |
| n samples | **105** |
| Histology composition | 49 PTC, 11 ATC, 45 patient-matched non-tumor (Normal) |
| Series Matrix | YES — `GSE33630_series_matrix.txt.gz` |
| Series Matrix size | ~30 MB gzipped |
| Series Matrix URL | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE33nnn/GSE33630/matrix/GSE33630_series_matrix.txt.gz |
| Histology metadata source | `!Sample_characteristics_ch1` — `pathological diagnostic: {anaplastic, papillary, non-tumor}` |
| Mutation / survival | none in GEO record |
| Verdict | **GO** — strongest cohort: full Normal / PTC / ATC ladder, all panel genes (22/22) mappable on GPL570. |

### GSE29265 — Tomás et al. — Sporadic vs post-Chernobyl PTC + ATC

| Field | Value |
|---|---|
| Title | "Sporadic vs. Post-Chernobyl Papillary vs. Anaplastic Thyroid Cancers" |
| Platform | GPL570 — all samples |
| n samples | **49** |
| Histology composition | 20 PTC (10 sporadic + 10 Chernobyl), 9 ATC, 20 patient-matched nontumor |
| Series Matrix | YES — `GSE29265_series_matrix.txt.gz` |
| Series Matrix size | ~15 MB gzipped |
| Series Matrix URL | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE29nnn/GSE29265/matrix/GSE29265_series_matrix.txt.gz |
| Histology metadata source | `!Sample_title` — strings include "Anaplastic thyroid carcinoma", "Papillary thyroid carcinoma", "Patient-matched non-tumor control" |
| Additional metadata | `origin` (sporadic / Chernobyl / NA), `dominant variant` (papillary / follicular / solitary), capsule invasion, extra-thyroid invasion, regional metastasis, age |
| Mutation / survival | none in GEO record |
| Verdict | **GO** — independent ATC vs PTC vs Normal cohort with origin annotation; same 22/22 panel coverage. |

### GSE65144 — von Roemeling et al. — ATC vs Normal

| Field | Value |
|---|---|
| Title | "Gene array analysis of anaplastic thyroid carcinoma tissue versus matched and unmatched normal thyroid tissue" |
| Platform | GPL570 — all samples |
| n samples | **25** |
| Histology composition | 12 ATC, 13 Normal thyroid (1 matched, 11 unmatched ATC; 1 matched, 12 unmatched Normal) |
| Series Matrix | YES — `GSE65144_series_matrix.txt.gz` |
| Series Matrix size | ~5.7 MB gzipped |
| Series Matrix URL | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE65nnn/GSE65144/matrix/GSE65144_series_matrix.txt.gz |
| Histology metadata source | `!Sample_characteristics_ch1` — `tissue type: {anaplastic thyroid carcinoma (ATC), normal thyroid}`. ATC subtype (spindle / squamoid / spindle+giant / etc.) recorded but not used for primary contrast. |
| Mutation / survival | none in GEO record |
| Verdict | **GO** — third independent ATC vs Normal cohort; small but clean. 22/22 panel coverage. |

---

## 2. Reused asset

GPL570 platform annotation (`GPL570_full.soft`, ~82 MB) is reused from
`project/results/p_landa_2016/raw/GPL570_full.soft`. **Not re-downloaded.** Already
covered by `.gitignore` line for `project/results/p_landa_2016/raw/`.

---

## 3. New gitignore line

```
project/results/p_gpl570_validation/raw/
```
appended to repo root `.gitignore`. The three Series Matrix .txt.gz files (~50 MB combined) are reproducible from the URLs in §1 and must not be committed.

---

## 4. Out-of-scope / explicitly NOT acquired

- GSE33538, GSE6004, GSE53157, GSE60542, GSE82208, any other GEO dataset
- TCGA methylation
- DepMap / CCLE / PRISM
- TCGA WSI / pathology
- Any RNA-seq cohort (this pack is GPL570-only by design)

---

## 5. Timebox check

- Access check: 30 min (within 2 h budget).
- Data download: 3 datasets in parallel, total ~3 min.
- Series Matrix usable for all three datasets — no FAIL/skip outcomes; no raw CEL needed.
