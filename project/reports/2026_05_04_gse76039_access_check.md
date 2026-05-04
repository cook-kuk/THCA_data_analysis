# GSE76039 — access check (Step 1)

**Date:** 2026-05-04 · **Owner:** Seungho Cook
**Mode:** Read-only access verification → small download. **No raw FASTQ alignment, no CEL re-normalisation, no GPU/RunPod, no TCGA WSI, no H&E-DM1 retry.**
**Scope:** Confirm GSE76039 has a usable processed expression matrix before further pipeline runs.
**Verdict:** **GO.** gcRMA-normalised Series Matrix available; all 22 needed panel genes mappable through GPL570 annotation.

---

## 1. Landing-page summary

| Field | Value |
|---|---|
| Title | "Genomic and Transcriptomic Hallmarks of Poorly-Differentiated and Anaplastic Thyroid Cancers" |
| Series | GSE76039 |
| Status | Public (Feb 15 2016) |
| Submission / last update | Dec 15 2015 / Mar 17 2021 |
| Citation | PubMed 26878173 (Landa 2016 JCI 126(3):1052–1066) + PubMed 33654411 (follow-up) |
| Contributors | Iñigo Landa, Tihana Ibrahimpasic, Laura Boucai, Rileen Sinha, Jeffrey A. Knauf, Snjezana Dogan, Ronak H. Shah, Julio C. Ricarte-Filho, Gnana P. Krishnamoorthy, Bin Xu, Nikolaus Schultz, Michael F. Berger, Chris Sander, Barry S. Taylor, Ronald Ghossein, Ian Ganly, James A. Fagin |
| Total samples | 37 |
| Histology breakdown | 17 poorly-differentiated thyroid carcinoma (PDTC) + 20 anaplastic thyroid carcinoma (ATC) |
| **Important caveat** | **No PTC, no normal samples in this series.** Series Matrix is the advanced-disease subset of the larger Landa 2016 work (n=117 PDTC+ATC for targeted DNA-seq; n=37 used for transcriptomics). |
| Platform | GPL570 — Affymetrix Human Genome U133 Plus 2.0 Array (54,675 probes) |
| Modality | Microarray (NOT RNA-seq, contrary to strategy-memo §3.1 description; corrected here) |
| Pre-processing | **gcRMA-normalised** (per `!Sample_data_processing` field) — provided as Series Matrix log2 expression |
| Sex | 28 female / 9 male |
| Tumor type detail | 31 primary; 4 recurrent in neck; 1 metastasis; 1 recurrent/persistent metastatic; 1 primary residual |
| BRAF/RAS/TERT mutation status | **NOT in GEO characteristics.** Available only via the Landa 2016 paper supplementary (out of scope for this acquisition). |

---

## 2. Supplementary file inventory (per GEO record)

| Filename | Size | Content | Action |
|---|---|---|---|
| `GSE76039_RAW.tar` | 164.7 MB | Raw CEL files (pre-gcRMA) | **NOT downloaded** — gcRMA-normalised Series Matrix supersedes |
| `GSE76039_series_matrix.txt.gz` | 6.9 MB | gcRMA log2 + per-sample metadata in one file | **DOWNLOADED** |
| GPL570 platform annotation (separate, via GEO query CGI) | 85.6 MB | Probe → Gene Symbol → Entrez mapping | **DOWNLOADED** (one-time, cacheable for any GPL570 study) |

---

## 3. Files acquired (CPU-only, no GPU/RunPod)

| Path | Size | md5 |
|---|---:|---|
| `project/results/p_landa_2016/raw/GSE76039_series_matrix.txt.gz` | 7,240,015 B (6.9 MB) | `b81615a5f5cffe27b10ba136111f8efa` |
| `project/results/p_landa_2016/raw/GPL570_full.soft` | 85,601,440 B (85.6 MB) | `ab5b044c4795d307e9201cda9a0e6a4f` |

Gzip integrity verified (`gunzip -t` OK on Series Matrix).

---

## 4. Access-check verdict

| Criterion | Status |
|---|---|
| Processed expression matrix exists (not raw FASTQ / not raw CEL only) | ✓ gcRMA log2 in Series Matrix |
| Per-sample metadata present | ✓ histology / sex / tumor type via `!Sample_characteristics_ch1` |
| Probe → gene mapping resolvable | ✓ GPL570 SOFT annotation parsed (54,675 probes; "Gene Symbol" column populated) |
| All 22 panel-union genes mappable | ✓ (8/8 RAI panel + 8/8 NONOVERLAP + all mechanism genes; see Step 4 in `2026_05_04_gse76039_first_pass_report.md`) |
| User-defined STOP condition (raw-FASTQ-only) | ✓ NOT triggered |
| Time-to-matrix budget (6 hours) | ✓ ≪ budget (~30 min total: download + parse + first scores) |

**→ Proceed to Step 4–7 (parse + score + tests + figures).** First-pass results in `2026_05_04_gse76039_first_pass_report.md`.

---

## 5. Marathon-discipline attestation

| Constraint | Status |
|---|---|
| New dataset (not in user's authorised list) | ✓ none — GSE76039 only |
| Raw FASTQ alignment / STAR / Salmon / kallisto | ✓ NOT performed |
| CEL re-normalisation (RMA / gcRMA from raw) | ✓ NOT performed (used Landa lab's gcRMA output) |
| GPU / RunPod / Azure burst | ✓ none |
| TCGA WSI download | ✓ none |
| H&E-DM1 retry | ✓ none |
| GSE33630 / TCGA methylation / Paper 3/4 / Paper 2 touch | ✓ none |
| Manuscript prose modification | ✓ none |
| Voice-protected sections (Hook / Aim / Discussion / Limitations / Cover ¶1 / Q9) | ✓ untouched |

---

## 6. Cross-references

- `2026_05_04_minimal_external_data_plan.md` §3.1 — acquisition slot for GSE76039 (the only authorised acquisition this turn)
- `2026_05_04_8gene_curated_vs_denovo_final_strategy.md` §12.6 — Nature Cancer reach gate (independent external molecular validation)
- `2026_05_04_paper1_dm1_full_molecular_only_lock.md` §2 — manuscript-safe envelope (16 claims) into which any GSE76039 result must fit
- `v17_landa2016_cite_save` — bib correction memo (`Landa2016JCI` is already in `2026_05_03_references.bib`)
- `v17_korean_K2_calibration` — within-sample-centered profile rationale (used here to handle microarray vs RNA-seq cross-platform robustness)
