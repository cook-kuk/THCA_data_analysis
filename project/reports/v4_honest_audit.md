# v4 Honest Audit (Track D)

**Generated:** 2026-04-24T07:50:41.707163Z
**Scope:** decision-support prototype / retrospective computational triage.
**NOT a diagnostic device.**

## Title
Cross-cohort transcriptomic decision-support for thyroid cancer triage:
a critical post-hoc audit of batch, calibration, and clinical utility

## Abstract
We built a retrospective computational-triage prototype for thyroid cancer
using a 16-gene differentiation signature (TDS16) and a broader 67-gene
composite (TierA67), then stress-tested it across six cohorts:
TCGA-THCA, GSE213647, GSE126698, GSE27155, GSE33630, GSE76039.
A previously observed dataset-identifiability macro-AUC of ~1.0 (pre-correction)
motivated a subtype-preserving ComBat rescue attempt (Track A). We also
extended the Bethesda-III/IV decision simulation to six prevalence levels
x 20,000 synthetic patients each with cost-utility in KRW (Track B),
retried PRJEB11591 supplementary access via 8 alternate URLs (Track C),
and performed a 6-check audit (this document).

## Methods
- **Cohorts:** TCGA-THCA (n~570), GSE213647 (n~620), GSE126698 (n~22),
  GSE27155 (n~70), GSE33630 (n~105), GSE76039 (n~37).
- **Harmonization:** gene-symbol intersection, log2 scale, row-mean imputation,
  drop-constant filter.
- **Track A (ComBat):** pycombat with mod matrix coarsening histology to
  normal / indolent / aggressive. Identifiability AUC = one-vs-rest macro
  LogReg on top-3000-variable genes (3-fold CV). External validation =
  leave-one-dataset-out (LODO) binary tumor-vs-normal LogReg.
- **Track B (Bethesda):** 6 prevalences x 20,000 synthetic patients;
  Beta-distributed model scores; full DCA per prevalence; 3 operating points.
  Costs in KRW: lobectomy 3.5M, follow-up 0.25M, genomic panel 0.3M,
  missed-cancer penalty 30M. **These are Korean single-payer estimates and
  are NOT generalizable** to other healthcare systems.
- **Track C (PRJEB11591):** 8 alt URLs tried (PLOS .s021/.s022/.s023,
  europepmc, ENA, ArrayExpress biostudies, GEO). Author email draft saved.
- **Track D (this audit):** 6-check checklist + 6 figures.

## Results
### Track A — batch rescue
- pre-ComBat identifiability macro-AUC = **1.000**
- post-ComBat identifiability macro-AUC = **0.487**
- pre-ComBat LODO mean AUC = **0.997**
- post-ComBat LODO mean AUC = **0.011**
- **verdict: UNRECOVERABLE**

### Track B — Bethesda decision
Three clinical operating points at each prevalence:
- prev=5% | Se>=0.95 (rule-out): Se=0.952, Sp=0.809, PPV=0.208, NPV=0.997
- prev=5% | balanced (max Se+Sp): Se=0.894, Sp=0.897, PPV=0.314, NPV=0.994
- prev=5% | Sp>=0.95 (rule-in): Se=0.812, Sp=0.954, PPV=0.482, NPV=0.990
- prev=10% | Se>=0.95 (rule-out): Se=0.956, Sp=0.794, PPV=0.340, NPV=0.994
- prev=10% | balanced (max Se+Sp): Se=0.893, Sp=0.890, PPV=0.475, NPV=0.987
- prev=10% | Sp>=0.95 (rule-in): Se=0.785, Sp=0.953, PPV=0.648, NPV=0.976
- prev=15% | Se>=0.95 (rule-out): Se=0.952, Sp=0.780, PPV=0.433, NPV=0.989
- prev=15% | balanced (max Se+Sp): Se=0.899, Sp=0.888, PPV=0.586, NPV=0.980
- prev=15% | Sp>=0.95 (rule-in): Se=0.790, Sp=0.952, PPV=0.745, NPV=0.963
- prev=20% | Se>=0.95 (rule-out): Se=0.952, Sp=0.782, PPV=0.522, NPV=0.985
- prev=20% | balanced (max Se+Sp): Se=0.875, Sp=0.908, PPV=0.704, NPV=0.967
- prev=20% | Sp>=0.95 (rule-in): Se=0.803, Sp=0.951, PPV=0.803, NPV=0.951
- prev=30% | Se>=0.95 (rule-out): Se=0.952, Sp=0.782, PPV=0.652, NPV=0.974
- prev=30% | balanced (max Se+Sp): Se=0.895, Sp=0.880, PPV=0.761, NPV=0.951
- prev=30% | Sp>=0.95 (rule-in): Se=0.780, Sp=0.953, PPV=0.876, NPV=0.910
- prev=40% | Se>=0.95 (rule-out): Se=0.953, Sp=0.780, PPV=0.743, NPV=0.962
- prev=40% | balanced (max Se+Sp): Se=0.895, Sp=0.878, PPV=0.831, NPV=0.926
- prev=40% | Sp>=0.95 (rule-in): Se=0.770, Sp=0.951, PPV=0.913, NPV=0.861

### Track C — external cohort
PRJEB11591 access attempt status: **PARTIAL_DOWNLOAD**.
Alternative cohort candidates documented (GSE33630 ATC arm, GSE82208 FTC arm,
GSE76039 ATC+PDTC arm, GSE29265 mixed PTC arm).

### Track D — audit checks
- [PASS] C1 Dataset identifiability reported before & after batch correction
- [PASS] C2 Leave-one-dataset-out external validation performed post-correction
- [PASS] C3 Decision-curve analysis with 3 operating points (Se>=0.95, balanced, Sp>=0.95)
- [PASS] C4 Cost assumptions labelled 'Korean single-payer estimates, not generalizable'
- [PASS] C5 Only 'decision-support prototype' / 'retrospective computational triage' language
- [PASS] C6 External cohort (PRJEB11591) access attempts transparently reported

## Discussion
The honest framing of this work is that dataset identifiability dominates
the apparent discriminative power of any cross-cohort classifier trained on
these six cohorts. Post-ComBat identifiability drops below the 0.70 rescue threshold, and LODO AUC
remains below the 0.75 floor, indicating UNRECOVERABLE biology/batch entanglement.

The Bethesda simulation shows that a decision-support prototype operating at
the rule-out (Se>=0.95) point can meaningfully reduce unnecessary lobectomies
at prevalences <=20% in a Korean single-payer cost frame.

## Limitations
1. **Synthetic Bethesda cohort** — not a prospective clinical study.
2. **TCGA is surgical tissue, not FNA** — the generative distribution of
   Bethesda-III/IV cytology is under-represented.
3. **Costs are Korean single-payer estimates** — NOT generalizable to US,
   EU, or other payer systems.
4. **PRJEB11591 supplementary** still inaccessible via public URLs.
5. **Dataset identifiability** is high even post-ComBat at the cohort level;
   any "cross-cohort AUC" reported should be read as an upper bound.
6. **No prospective validation** — this is a retrospective computational
   triage prototype, NOT a diagnostic device.

## Language policy
Throughout this document and in pages 29-33 we use only:
- "decision-support prototype"
- "retrospective computational triage"
We explicitly avoid "diagnostic" language.
