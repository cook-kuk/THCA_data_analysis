# GSE138042 — 8-gene panel and radioiodine refractoriness, with the clean contrast

**Date** 2026-08-06 (revised same day after obtaining the clinical annotation)
**Script** `scripts/gse138042_rair_panel_2026_08_06.py`
**Source** GEO GSE138042, Heliyon 2021, PMID 33748479 / PMC7970325, doi 10.1016/j.heliyon.2021.e06408.
Count matrix `GSE138042_mRNA_seq_thyroid.csv.gz` (33,513 expressed genes × 95 libraries) and
per-patient clinical annotation `mmc1.xlsx`, both public and both now local at
`/data/rai_atlas/raw/GSE138042/`.

> **How the clinical file was obtained.** Supplementary File 1 sits behind a JavaScript
> download gate on PMC that defeats scripted retrieval. The Europe PMC supplementary-files
> endpoint returns the whole set as a zip with no gate:
> `https://www.ebi.ac.uk/europepmc/webservices/rest/PMC7970325/supplementaryFiles`
> This is worth remembering — it works for any PMC article.

## The correction that matters

The first pass compared the 13 `RAIR-*` libraries against all 82 remaining libraries and
gave d = −1.03, P = 1.7 × 10⁻⁴. The clinical file shows that comparator was wrong for the
question: it contained 17 benign follicular adenomas, 3 medullary carcinomas, and 62
carcinomas of which only 10 have a recorded radioiodine outcome. The annotation column
`Radiation/Iodine therapy status` gives the contrast the paper actually defines:
**radioresistance 13 versus radiosensitivity 10**, all of them carcinomas that received
radioiodine.

| Contrast | n | Cohen's d | 95% CI | P |
|---|---|---|---|---|
| **PRIMARY — radioresistant vs radiosensitive** (all radioiodine-treated) | 13 / 10 | **−0.42** | −1.44 to +0.42 | **0.34** |
| Sensitivity — refractory vs all malignant (no adenoma, no medullary) | 13 / 62 | −0.94 | — | 4.1 × 10⁻⁴ |
| Context — refractory vs whole 82-library cohort | 13 / 82 | −1.03 | −1.80 to −0.52 | 1.7 × 10⁻⁴ |

**The head-to-head test is null.** Direction is as hypothesised — refractory tumours score
lower — but the interval spans zero and P = 0.34. With 13 versus 10 the study detects
d ≥ 1.24 at 80% power, and power at the observed effect is 0.16. This cannot adjudicate the
hypothesis either way.

## What the three rows together actually say

Removing the benign adenomas barely moves the effect (−1.03 → −0.94), so the large result is
**not** an artefact of benign tissue in the comparator. But the effect halves when the
comparator is narrowed to the ten patients with a documented radioiodine outcome
(−0.94 → −0.42). The most economical reading is that a substantial part of the large effect
reflects **advanced versus early disease** rather than radioiodine behaviour specifically —
the refractory samples are metastatic/advanced carcinomas including a poorly differentiated
case, while the 62-sample malignant comparator is dominated by primary resections.

This is the same lesson as the TCGA and GSE151179 analyses, arriving a third time: the
differentiation axis tracks disease stage and dedifferentiation robustly, and the extra step
to "radioiodine-specific" is exactly where the evidence thins.

## Per-gene, on the primary contrast

| Gene | Cohen's d | P | q (BH) |
|---|---|---|---|
| TG | −1.18 | 0.014 | 0.057 |
| TSHR | −1.10 | 0.012 | 0.096 |
| SLC5A5 (NIS) | −0.67 | 0.18 | 0.36 |
| TPO | −0.63 | 0.28 | 0.44 |
| NKX2-1 | −0.01 | 0.93 | 0.93 |
| FOXE1 | +0.30 | 0.60 | 0.68 |
| DIO1 | +0.45 | 0.34 | 0.45 |
| PAX8 | +0.72 | 0.13 | 0.34 |

TG and TSHR are nominally significant and do not survive correction. On the broad contrast
NIS was significant (q = 0.006); on the head-to-head it is not (q = 0.36). Reporting only the
broad contrast would therefore have overstated the mechanistic specificity, and that is worth
saying plainly rather than quietly dropping.

## Composition control (broad contrast)

Radioiodine-refractory tumours have lower immune content (d = −0.65, P = 0.049). In the
logistic model on the full 95 libraries, panel z retains OR = 0.288 (0.101–0.819), P = 0.020
with immune and stromal transcripts held constant, and immune content is itself significant
(OR 0.285, P = 0.024). The composition adjustment does not explain the broad effect — unlike
GSE151179, where the equivalent adjustment removed it entirely.

## Verdict

Report GSE138042 as a **supportive but underpowered** cohort:

> In an independent cohort with documented radioiodine outcome (GSE138042, n = 23 treated
> patients), radioiodine-refractory carcinomas scored lower on the differentiation panel than
> radiosensitive carcinomas, but the difference was not significant (Cohen's d = −0.42,
> 95% CI −1.44 to +0.42, P = 0.34), and the cohort is powered only for d ≥ 1.24. The larger
> effect seen against the full surgical cohort (d = −0.94 restricted to malignant samples,
> P = 4.1 × 10⁻⁴) is confounded by disease stage and should not be presented as a
> radioiodine-specific result.

The design target from the GSE151179 power calculation stands unchanged: **≈226 patients**
with a lesion-level radioiodine outcome. No public dataset reaches it. The nearest candidates
remain Zhang 2026 (n = 113 advanced DTC, GSA HRA011340), Mu 2024 (n = 214, HRA004166) and
Siraj 2022 / EGAS00001001788 (n = 158).

## Files

- `results/tables/gse138042_rair_panel_main_2026_08_06.tsv`
- `results/tables/gse138042_rair_per_gene_2026_08_06.tsv`
- `results/tables/gse138042_rair_adjusted_2026_08_06.tsv`
- `results/tables/gse138042_rair_per_sample_2026_08_06.tsv` (panel z, composition scores, RAI status, histology, BRAF, age, sex)
- `results/figures/figure_gse138042_rair_panel_2026_08_06.{png,pdf}`

Counts → CPM → log2(x+1) → per-gene z across all 95 libraries → panel score = mean z.
Seed 20260806, bootstrap 5,000 draws.
