# v4 Synthesis

**Generated:** 2026-04-24T07:52:18.328477+00:00

## Executive summary
- Track A verdict: **UNRECOVERABLE**
- Track B: 6 prevalences x 20,000 simulated patients; 18 operating points.
- Track C: PRJEB11591 alt-URL retry status: **PARTIAL_DOWNLOAD**
- Track D: 6/6 audit checks pass.

## Track verdicts
| Track | Status |
|------|--------|
| A - ComBat rescue | UNRECOVERABLE (pre_ident=1.000, post_ident=0.487, pre_lodo=0.997, post_lodo=0.011) |
| B - Bethesda decision | 18 operating points; KRW cost-utility NOT generalizable |
| C - PRJEB11591 retry | PARTIAL_DOWNLOAD |
| D - Honest audit | 6/6 checks pass |

## 4x4 Decision matrix (venue x feasibility/time/impact/risk)
| venue | feasibility | time | impact | risk(inv) | total |
|------|---|---|---|---|---|
| ML4H workshop (fallback) | 5 | 5 | 3 | 5 | 21 |
| MLCB workshop (fallback) | 5 | 5 | 2 | 5 | 19 |
| Bioinformatics (methods) | 4 | 3 | 4 | 3 | 18 |
| JCO Precision Oncology (clinical) | 3 | 2 | 5 | 2 | 17 |

Winner: **ML4H workshop (fallback)** (total=21).

## 14-day daily task list
- D01: freeze v4 verdicts; circulate honest_audit.md to co-authors
- D02: finalize mini-paper skeleton outline; draft abstract & methods
- D03: re-run Track A with non-linear correction (BBKNN / Harmony on PCA) as ablation
- D04: wet-lab contact for PRJEB11591 author email; send outreach
- D05: extend Bethesda sim with miscalibration sensitivity analysis
- D06: replicate operating-point analysis with real TCGA probability outputs
- D07: literature scan for comparable commercial panels (Afirma, ThyGenX)
- D08: consolidate all v4 figures into supplementary PDF
- D09: prepare 1-slide and 5-slide pitch decks for Gemma hackathon
- D10: draft submission cover letter (JCO PO primary; Bioinformatics backup)
- D11: IRB / data-use documentation review for each cohort
- D12: freeze LaTeX main text draft (v0.1)
- D13: internal review pass; address 3 reviewers' honest-concerns list
- D14: submit to JCO PO; parallel preprint to bioRxiv
