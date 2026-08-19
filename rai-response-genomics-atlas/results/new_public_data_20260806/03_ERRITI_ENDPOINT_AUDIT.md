# 03 — ERRITI endpoint audit and the radioiodine endpoint hierarchy

**Grade C · ENDPOINT_REFERENCE_ONLY.** No analysis of ERRITI data is possible or attempted.

## What was obtained

Both ERRITI supplementary documents (Weber et al., *Clin Cancer Res* 2022, ccr-22-0437),
via the AACR figshare mirror:

| File | Size | Content |
|---|---|---|
| `erriti_39940925.docx` | 35 KB | Supplemental Table S1 — solid-tumour sequencing panel gene list; full eligibility criteria |
| `erriti_39940928.docx` | 13 KB | figure legends |

**There is no per-patient data table.** The individual uptake, dosimetry and thyroglobulin
values exist only as graphics inside main-text Figure 2. Per the user's instruction —
*"ERRITI에는 expression/proteomics가 없으면 8-gene score 분석을 시도하지 마라"* —
no eight-gene analysis was attempted, and none is possible: there is no molecular matrix
of any kind, only a panel gene list.

## Why ERRITI still matters

Two things are worth keeping.

**1. The mutation panel sequences the iodide machinery.** ERRITI's solid-tumour panel
includes **SLC5A5 (exons 2, 3, 6, 11)** and all of **TSHR**. Very few clinical panels do.
This is a useful precedent to cite when arguing that iodide-handling genes deserve
routine capture — most panels would not detect variation in the gene with the most direct
mechanistic claim on iodide transport.

**2. It measured the full endpoint ladder.** ERRITI is the only source in our inventory
that captured every tier at once: I-123 scintigraphy uptake, **I-124 PET dosimetry with
absorbed dose in Gy**, I-131 therapy, thyroglobulin trajectory, RECIST, and FDG-PET.

That makes it the **reference definition** against which the vaguer endpoints elsewhere
can be positioned.

## The hierarchy

Written to `rai-response-genomics-atlas/results/rai_endpoint_hierarchy.tsv`.
Seven tiers, ordered by how directly each measures the radioiodine exposure itself.

| Tier | Endpoint | Who has it |
|---|---|---|
| 1 | **Absorbed dose to lesion** (I-124 PET dosimetry, Gy) | ERRITI; some redifferentiation arms. **No public omics dataset.** |
| 2 | **Qualitative uptake at a lesion** (I-123/I-131 scan) | ERRITI; GSE151179; Liu 2024; E-MTAB-12837/12900; all 9 redifferentiation arms |
| 3 | **Biochemical response** (thyroglobulin trajectory) | ERRITI; Zhang 2026; Borowczyk 2024 (n=646, largest) |
| 4 | **Structural response by RECIST after RAI** | ERRITI; Zhang 2026; redifferentiation arms |
| 5 | **Clinician-assigned refractory status** (ATA 2015) | GSE138042; Zhang 2026; Liu 2024; DECISION; E-MTAB |
| 6 | **Registry treatment-response code** | TCGA-THCA only |
| 7 | **Downstream clinical outcome after RAI** | TCGA-THCA; Zhang 2026 (PFS) |

## What the hierarchy implies for our claims

**The tier gap is the whole problem.** Our only significant results (TCGA
OR = 0.221 and OR = 0.224) sit at **tier 7**, the tier furthest from the exposure and the
one that cannot distinguish a radioiodine effect from prognosis, because the
treatment-by-score interaction is unestimable there (2 events among non-treated patients).

Everything we tested at **tier 2**, the tier that actually records whether radioiodine
reached disease, was null: GSE151179 uptake d = +0.37 (P = 0.53, collapsing to β = +0.035
under purity adjustment), E-MTAB redifferentiation response d = +0.34 (P = 0.80),
Liu 2024 uptake d = −0.51 (P = 0.088, trend only).

**No dataset anywhere gives us tier 1 with a transcriptome.** That is the single
structural hole in the field, and it is worth stating plainly in the manuscript rather
than leaving a reviewer to notice it.

## Cross-tier heterogeneity is a stated reason for our meta-analysis I²

The driver-axis meta-analysis pooled Siraj 2022, Boucai 2023 and Zhang 2026 — all tier 5,
but with three different local definitions of refractoriness. I² = 63% is consistent with
definition heterogeneity, not only with sampling variation. The hierarchy table is the
evidence for saying so.

## Files

- `rai-response-genomics-atlas/results/rai_endpoint_hierarchy.tsv`
- `external_data/erriti_redifferentiation/source/erriti_39940925.docx`
- `external_data/erriti_redifferentiation/source/erriti_39940928.docx`
