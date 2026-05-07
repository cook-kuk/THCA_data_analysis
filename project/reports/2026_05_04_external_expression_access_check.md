# External expression validation — Step 1 access check

**Date:** 2026-05-04
**Scope:** Paper 1 driver-orthogonal RAI-lineage / DM1 differentiation axis — external thyroid expression replication (processed-only, no raw FASTQ/CEL alignment).
**Working directory:** `project/results/p_external_expression_validation/`

GEO landing pages were fetched for the user-supplied candidate list only. **No additional GEO search/download performed.**

---

## Access summary table

| # | Tier | GEO | Title (short) | Organism | Platform | Modality | n total | Histology groups (n) | Series Matrix | Processed in matrix | Raw form | Verdict |
|---|------|-----|---------------|----------|----------|----------|---------|----------------------|---------------|---------------------|----------|---------|
| 1 | T1 | GSE126698 | IGF2BP1 ATC marker | Hs | GPL15456 (Illumina HiScanSQ) | RNA-seq | 28 | ATC 10 / PTC 6 / FTC 6 / Normal 6 | Yes (TXT) | **Sample-level expression NOT confirmed in matrix** — supp files are DE tables (`DE_*.csv.gz`); raw=SRA (PRJNA523137 / SRP186236) | SRA FASTQ | **CONDITIONAL** — open Series Matrix to confirm whether per-sample counts exist; if matrix only carries DE summary, mark FAIL (no raw FASTQ alignment allowed) |
| 2 | T2 | GSE33630 | Normal vs PTC vs ATC | Hs | GPL570 | Affy U133+2 array | 105 | ATC 11 / PTC 49 / Normal 45 | Yes (TXT) | Yes ("Processed data included within Sample table") | CEL TAR (849 MB) | **GO** — processed-only path viable |
| 3 | T2 | GSE29265 | Sporadic vs Chernobyl PTC + ATC | Hs | GPL570 | Affy U133+2 array | 49 | ATC 9 / PTC 20 / Paired-normal 20 | Yes (TXT) | Yes (sample table) | CEL TAR (213 MB) | **GO** |
| 4 | T2 | GSE65144 | ATC vs matched/unmatched normal | Hs | GPL570 | Affy U133+2 array | 25 | ATC 12 / Normal 13 | Yes (TXT) | Yes (sample table) | CEL TAR (107 MB) | **GO** |
| 5 | T2 | GSE53157 | PDTC progression series | Hs | GPL570 | Affy U133+2 array | 27 | PDTC 5 / cPTC 7 / fvPTC 8 / FTC 4 / Normal 2 / Commercial pool 1 | Yes (TXT) | Yes (sample table) | CEL TAR (126 MB) | **GO** (drop commercial pool sample) |
| 6 | T3 hold | GSE53072 | Tiny ATC vs normal | Hs | GPL6244 (HuGene-1.0-ST) | Array | 9 | ATC 5 / Normal 3 / Pool 1 | Yes (TXT) | Yes | CEL TAR (40 MB) | **HOLD** per user (tiny + non-GPL570) |
| 7 | T3 hold | GSE60542 | PTC primary vs nodal mets + LN | Hs | GPL570 | Affy U133+2 array | 92 | PTC N+ 17 / PTC N0 11 / LN-met 17 / Paired-normal 24 / Normal LN 4 | Yes (TXT) | Yes | CEL TAR (646 MB) | **HOLD** per user (lymphoid confounding) |
| 8 | T3 hold | GSE120177 | CDK7/THZ1 in ATC cell lines | Hs | GPL23227 (BGISEQ-500) | RNA-seq + ChIP | 12 | BCPAP / CAL-62 cell lines × DMSO/THZ1/Input/H3K27ac | Yes (TXT) | Yes (small TAR) | SRA | **HOLD** per user (Paper 9 scope only; cell line perturbation, not patient histology) |

---

## Per-dataset notes

### GSE126698 — Tier 1, conditional
- RNA-seq, only ATC/PDTC/PTC/FTC/normal patient series in this set; no PDTC arm despite the user mention of "ATC/PDTC/PTC/FTC/normal" — **PDTC is absent**.
- Supplementary CSVs listed are DE tables (ATC vs others, total). Per-sample counts must come from the Series Matrix itself or from `GSE126698_*` if a counts file exists; will verify at Step 2 download time.
- If the Series Matrix carries no per-sample expression (DE-only summary), this dataset will be marked FAIL at Step 2 — **raw FASTQ alignment is forbidden**.

### GSE33630, GSE29265, GSE65144, GSE53157 — Tier 2 GPL570 pack
- All four are GPL570 Affy U133+2 with processed values inside the Series Matrix. Standard log2-RMA-like values are expected; absolute scaling will not be cross-pooled (per spec, per-dataset z-scoring only).
- GSE29265 has paired-normal design (matched controls per patient).
- GSE53157 has a "Commercial RNA pool" sample to drop.
- GSE60542 (held) overlaps with GSE29265 / GSE33630 lab (Maenhaut/Detours group); cross-study sample re-use risk noted but irrelevant here since GSE60542 is held.

### GSE53072, GSE60542, GSE120177 — Tier 3, hold
- Per user instructions, no download for this sweep. Recorded only.

---

## What this step did NOT do
- No additional GEO search beyond the 8 listed accessions.
- No FTP downloads yet.
- No raw CEL/FASTQ acquisition.
- No alignment, no probe-level reprocessing.

## Step 2 plan (processed-only)
- Download Series Matrix `*_series_matrix.txt.gz` for: GSE126698, GSE33630, GSE29265, GSE65144, GSE53157.
- Download SOFT platform annotation for: GPL570, GPL15456 (probe → gene mapping).
- If GSE126698 Series Matrix lacks per-sample expression, mark FAIL and skip — will not pull SRA.

---

## Forbidden-action audit (this step)
- ✅ No new GEO accessions discovered or queued.
- ✅ No raw download.
- ✅ No alignment.
- ✅ No GPU/RunPod.
- ✅ No manuscript prose touched.
- ✅ No voice-protected section drafted.
