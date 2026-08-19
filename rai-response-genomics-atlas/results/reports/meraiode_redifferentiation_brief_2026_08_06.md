# E-MTAB-12837 / 12900 — redifferentiation + I-131 response, and why the panel cannot be tested here

**Date** 2026-08-06 · **Script** `scripts/meraiode_redifferentiation_panel_2026_08_06.py`
**Source** ArrayExpress E-MTAB-12837 (BRAF arm, n = 13; ENA PRJEB61187) and E-MTAB-12900
(RAS arm, n = 8; ENA PRJEB61819). Trial paper: Leboulleux S et al., *Clin Cancer Res*
2023;29(13):2401-9 — phase II dabrafenib + trametinib + ¹³¹I redifferentiation.
Raw FASTQ (21 runs, 6.5 GB) downloaded and re-quantified locally with kallisto 0.50.1
against a 93-transcript index covering the eight panel genes;
`/data/rai_atlas/raw/E-MTAB/`.

## Why we went after this

These are the only public per-patient datasets pairing a tumour molecular profile with the
response to an intervention whose explicit purpose is to restore iodine handling. Response
is recorded per sample as RECIST and PERCIST responder / non-responder, and — unusually —
`tumor cell content` is recorded too, so the composition covariate that decided the
GSE151179 and TCGA analyses is available. Neither submission deposits a processed expression
matrix, which is presumably why nobody has reused them.

## Hard finding: the assay cannot measure three of the eight genes

The trial used the **HTG EdgeSeq Oncology Biomarker Panel**, a 2,549-gene pan-cancer probe
set. Re-quantification gives zero counts across all 21 libraries for three genes:

| Measurable (5) | Median CPM | Absent from probe set (3) |
|---|---|---|
| PAX8 | 588,375 | **SLC5A5 (NIS)** — 0/21 nonzero |
| FOXE1 | 126,846 | **TG** — 0/21 nonzero |
| TSHR | 113,534 | **DIO1** — 0/21 nonzero |
| NKX2-1 | 95,912 | |
| TPO | 54,803 | |

This was confirmed independently against a published copy of the HTG panel gene list, where
NKX2-1 appears under its legacy symbol TTF1. The three missing genes are the iodide-handling
machinery — including the sodium-iodide symporter, the single most mechanistically relevant
transcript for radioiodine. **Any score computed here is a 5-gene transcription-factor-heavy
surrogate, not the panel.**

## Result: directionally correct, firmly null

| Endpoint | n responder / non-responder | Cohen's d | AUC | P |
|---|---|---|---|---|
| RECIST | 9 / 12 | +0.34 | 0.54 | 0.80 |
| PERCIST | 11 / 7 | +0.60 | 0.64 | 0.37 |
| Clinical history | 9 / 12 | +0.34 | 0.54 | 0.80 |

Responders sit higher on the 5-gene score in every endpoint, which is the hypothesised
direction, and nothing approaches significance. With 9 versus 12 the study detects
d ≈ 1.3 at 80% power, so a real effect of the size seen in GSE138042 (d ≈ 1.0 on the full
panel) would have been missed here more often than not. This is a power failure, not
evidence against the hypothesis, and it should be reported as such or not at all.

## Verdict

Do not put this in the manuscript as a validation attempt. It is worth one sentence in the
data-availability or limitations discussion: the only public redifferentiation-trial
transcriptomes use a targeted panel lacking NIS, TG and DIO1, so the field currently has no
public dataset in which a thyroid differentiation score can be tested against restored
iodine uptake. That is a genuine gap statement and it is defensible because we did the work.

## Files

- `results/tables/meraiode_panel_per_sample_2026_08_06.tsv`
- `results/tables/meraiode_response_test_2026_08_06.tsv`
- `results/figures/figure_meraiode_redifferentiation_2026_08_06.{png,pdf}`
- Raw FASTQ + kallisto output retained at `/data/rai_atlas/raw/E-MTAB/`

Note on provenance: E-MTAB-12900 is labelled the "KRAS cohort" but the samples carry NRAS
and HRAS mutations, and no publication is linked to it. Treat it as a sibling RAS arm of the
same programme, unverified.
