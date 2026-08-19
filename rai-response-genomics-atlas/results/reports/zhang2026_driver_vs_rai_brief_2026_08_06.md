# Zhang 2026 (n = 113 advanced DTC) — molecular subtype predicts radioiodine refractoriness, driver mutation does not

**Date** 2026-08-06 · **Script** `scripts/zhang2026_extract_and_test_2026_08_06.py`
**Source** Zhang T, Cui W, Tang H, … Shi X. "Proteogenomic characterization delineates
clinically relevant subtypes of advanced differentiated thyroid cancer." *Cell Reports
Medicine* 2026;7(3):102661, doi 10.1016/j.xcrm.2026.102661, PMID 41794039, PMCID PMC13006413.
Fudan University Shanghai Cancer Center.

## How the data was obtained without contacting anyone

Both deposited accessions are gated — proteomics at iProX IPX0011848000 (directory returns
403, no public manifest) and sequencing at GSA-Human HRA011340 (controlled, DAC HDAC005903).
The per-patient tables, however, are printed in the open supplemental PDF:

```
https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13006413/supplementaryFiles   → mmc1.pdf
pdftotext -layout mmc1.pdf -                                                       → one line per patient
```

Table S1 renders as 113 fixed-width rows with fourteen fields, and Table S2 as 175
per-patient mutation rows. Both parsed cleanly. Local copies:
`results/tables/zhang2026_patient_table_2026_08_06.tsv` and `…_mutations_….tsv`.

## Cohort

113 patients with advanced differentiated thyroid carcinoma, and the balance is what every
previous cohort lacked:

| Field | Distribution |
|---|---|
| **RAI sensitivity** | **Refractory 58 / Avid 55** |
| Biochemical response to RAI | G1 24 · G2 16 · G3 63 · NA 10 |
| Proteomic consensus subtype | CC1 43 · CC2 40 · CC3 30 |
| Histology | cPTC 79 · fvPTC 18 · FTC 16 |
| Driver class (from Table S2) | BRAF V600E 46 · BRAF/RAS-negative 56 · RAS hotspot 11 |

A 51% refractory fraction is the regime in which a discriminator has room to work — the
opposite of TCGA-THCA, where 86% achieve complete response.

## Result

| Classifier | Refractory fraction by class | χ² P |
|---|---|---|
| **Proteomic consensus subtype** | **CC1 21% → CC2 57% → CC3 87%** | **1.44 × 10⁻⁷** |
| Driver mutation class | BRAF V600E 61% · BRAF/RAS-negative 48% · RAS hotspot 27% | 0.109 |
| Histology | FTC 56% · cPTC 54% · fvPTC 33% | 0.248 |

CC1 versus CC3 alone: **odds ratio 24.6, Fisher P = 2.0 × 10⁻⁸**.

The biochemical response grade moves with subtype in the same direction: CC1 is 18/40 G1
(≥50% thyroglobulin decline), while CC3 is 29/30 G3 (≥10% thyroglobulin rise).

## Why this matters for our manuscript

This is an independent, adequately sized, advanced-disease cohort that makes both halves of
our argument at once:

1. **Driver mutation class does not identify radioiodine-refractory disease** (P = 0.109).
   This replicates the Siraj 2022 result exactly (n = 158, P = 0.67,
   `siraj2022_driver_vs_rai_brief_2026_08_06.md`). Two independent cohorts, 271 patients
   total, both null. The BRAF/RAS axis is not the axis that governs iodine handling.
2. **An expression-derived molecular class does identify it, strongly** (P = 1.4 × 10⁻⁷,
   CC1 vs CC3 OR 24.6). Their CC subtypes come from consensus clustering of the global
   proteome — conceptually the same manoeuvre as deriving DM1/DM2 from an eight-gene
   transcript panel, but at protein level and in advanced disease.

Together these say the differentiation/expression axis is the right axis and the driver axis
is not, which is the premise the manuscript rests on, now demonstrated in someone else's
data with our own re-analysis behind it rather than borrowed from their conclusions.

## What it does not do

It does not test our panel. Their subtypes are proteome-wide consensus clusters, not the
eight-gene score, and the per-patient protein quantities live behind iProX. Testing the
panel here requires either the iProX matrix or the GSA sequencing, both of which need a
request to Xiao Shi (`xshi11@fudan.edu.cn`), who is simultaneously the paper's lead contact
and the GSA data-access committee contact — one email covers both.

Also note the direction of the non-significant driver trend: BRAF V600E tumours are
*numerically more* refractory here (61%), while in the manuscript's framing the
BRAF/RAS-negative compartment is the one of interest. That compartment sits at 48%, close to
the cohort mean. Do not present the driver-negative compartment as intrinsically
radioiodine-refractory — neither this cohort nor Siraj supports that.

## Files

- `results/tables/zhang2026_patient_table_2026_08_06.tsv` (113 rows, 15 fields)
- `results/tables/zhang2026_mutations_2026_08_06.tsv` (175 rows)
- `results/tables/zhang2026_driver_vs_rai_2026_08_06.tsv`
- `results/figures/figure_zhang2026_driver_rai_2026_08_06.{png,pdf}`
- Source PDF retained at `/data/rai_atlas/external/supp_bulk/mmc1.pdf`

Driver calls from Table S2: BRAF V600E; RAS = NRAS/HRAS/KRAS at codons 12, 13 or 61;
everything else BRAF/RAS-negative. The paper's own genomic analyses use the same panel.
