---
title: "MSK-IMPACT thyroid cohort — enrichment bias caveat"
date: 2026-04-29
purpose: PROMPT E deliverable. Document MSK referral bias for paper Methods + Discussion.
---

# Bias quantification (MSK n=117 vs TCGA-THCA n=513)

| Variable | MSK-IMPACT | TCGA-THCA | Test | p |
|----------|-----------|-----------|------|---|
| % PTC | 0.0% | 94.0% | chi² | 6.6e-131 |
| % PDTC | 71.8% | 0.0% | — | — |
| % ATC | 28.2% | 0.0% | — | — |
| Median age | 61.0 | 46.3 | Mann-Whitney | 8.5e-12 |
| M1 (distant met) frac | 37.6% | — | — | — |
| Stage IV frac (TCGA) | — | 9.9% | — | — |

The MSK 2016 thyroid panel is by design enriched for **advanced / refractory thyroid cancer** (PDTC and ATC referrals to a tertiary cancer center). It is **not a representative thyroid cancer cohort** in the epidemiological sense.

# Implication for the 8-gene panel paper

The MSK cohort is appropriate for testing the panel's behavior in **dedifferentiated / aggressive thyroid cancer** (the right tail of the differentiation continuum) but is **not a generalizability validation set** for the broader PTC population. Findings in MSK should be reported as:

> "MSK-IMPACT thyroid samples (n=117, all PDTC/ATC referrals; median age 61 vs TCGA 46; chi-squared p=6.6×10⁻¹³¹ vs TCGA histology distribution) constitute an enriched advanced-disease cohort and were used here to test whether the 8-gene differentiation signature behaves consistently in the dedifferentiated tail. They are not intended as an epidemiologically representative thyroid cancer validation set."

# Suggested Methods + Discussion text

**Methods § Cohorts.** "We deliberately included MSK-IMPACT (Landa 2016) as an enriched advanced-disease comparator because PDTC and ATC are systematically underrepresented in the TCGA-THCA discovery cohort. We do not compute generalizability metrics on MSK without weighting; instead we use MSK to characterize signature behavior at the dedifferentiated extreme of the continuum."

**Discussion § Limitations.** "The MSK-IMPACT thyroid panel reflects tertiary referral bias and overrepresents PDTC/ATC by ~70-fold and ~unbounded ratio respectively versus TCGA-THCA. This bias is intentional in our study design (we specifically wanted dedifferentiated cases) but cautions against generalizing MSK-derived effect sizes to community PTC."

# Files

- `histology_comparison.tsv` — MSK vs TCGA harmonized histology bucket counts
- `msk_bias_summary.json` — programmatic summary
- `msk_bias_panel.pdf/png` — 2-panel figure (histology + age)

Reproduce: `python project/notebooks_or_scripts/v17_msk_bias_doc.py`
