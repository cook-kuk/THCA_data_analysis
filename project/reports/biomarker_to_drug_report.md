# Biomarker → Drug Discovery Report

- Generated: 2026-04-24 11:31:09
- Biomarker source: `biomarker_validated.tsv`
- Targets analysed: **8**
- Compounds collected: **15**
- Literature records: **27**

## Top 8 targets

| Gene | Novel | Classification | #Cmpds | Best pChEMBL | Thyroid papers | Top PMID | Top Compound |
|---|---|---|---|---|---|---|---|
| TACSTD2 | True | novel_target | 0 | — | 4 | 41413112 | — |
| TMPRSS4 | True | novel_target | 0 | — | 5 | 41656803 | — |
| PLEKHA6 | True | novel_target | 0 | — | 0 | — | — |
| CYP1B1 | True | validated_target | 10 | 8.75 | 5 | 40900788 | CHEMBL3132932 |
| LDLR | True | emerging_target | 4 | 6.22 | 5 | 41340871 | CHEMBL115992 |
| GABRB2 | True | emerging_target | 1 | 6.89 | 5 | 40900788 | CANNABIDIOL |
| B3GNT3 | True | novel_target | 0 | — | 0 | — | — |
| PTPRE | True | novel_target | 0 | — | 3 | 36654463 | — |

## Pharma-quality highlights (3 novel targets)
### TACSTD2 — novel_target
**TACSTD2** surfaces from biomarker discovery with high novelty_score and no known high-confidence chemical starting points in ChEMBL. For a pharma portfolio this is a **first-in-class opportunity** — the target hypothesis is de-risked by replicated differential expression in thyroid cancer cohorts, while the chemical matter is greenfield (fragment screening / DNA-encoded library campaigns would be the logical next step). Thyroid-specific literature count = 4.

### TMPRSS4 — novel_target
**TMPRSS4** surfaces from biomarker discovery with high novelty_score and no known high-confidence chemical starting points in ChEMBL. For a pharma portfolio this is a **first-in-class opportunity** — the target hypothesis is de-risked by replicated differential expression in thyroid cancer cohorts, while the chemical matter is greenfield (fragment screening / DNA-encoded library campaigns would be the logical next step). Thyroid-specific literature count = 5.

### PLEKHA6 — novel_target
**PLEKHA6** surfaces from biomarker discovery with high novelty_score and no known high-confidence chemical starting points in ChEMBL. For a pharma portfolio this is a **first-in-class opportunity** — the target hypothesis is de-risked by replicated differential expression in thyroid cancer cohorts, while the chemical matter is greenfield (fragment screening / DNA-encoded library campaigns would be the logical next step). Thyroid-specific literature count = 0.


## Caveat
This is computational triage only — not a drug recommendation. All compound activities come from heterogeneous public assays; thyroid-specific efficacy is not guaranteed. Confirmatory in-vitro / in-vivo work required before any claim.

## Run log
    [11:30:43] PHASE 0: polling for biomarker_validated.tsv (max 1200s)
    [11:30:43] PHASE 0: found validated biomarkers after 0s
    [11:30:44] PHASE 0: loaded 2843 biomarker rows from biomarker_validated.tsv
    [11:30:44] PHASE 0: selected targets: ['TACSTD2', 'TMPRSS4', 'PLEKHA6', 'CYP1B1', 'LDLR', 'GABRB2', 'B3GNT3', 'PTPRE']
    [11:30:44] PHASE 1-3: TACSTD2
    [11:30:45]     pubmed thyroid hits: 4 (cancer total≈14)
    [11:30:47]     chembl target=CHEMBL3856163 pref=Tumor-associated calcium signal transducer 2 compounds=0
    [11:30:48] PHASE 1-3: TMPRSS4
    [11:30:49]     pubmed thyroid hits: 5 (cancer total≈60)
    [11:30:50]     chembl target=CHEMBL2331048 pref=Transmembrane protease serine 4 compounds=0
    [11:30:51] PHASE 1-3: PLEKHA6
    [11:30:51]     pubmed thyroid hits: 0 (cancer total≈0)
    [11:30:51]     chembl target=None pref=None compounds=0
    [11:30:52] PHASE 1-3: CYP1B1
    [11:30:53]     pubmed thyroid hits: 5 (cancer total≈242)
    [11:30:57]     chembl target=CHEMBL4878 pref=Cytochrome P450 1B1 compounds=10
    [11:30:57] PHASE 1-3: LDLR
    [11:30:59]     pubmed thyroid hits: 5 (cancer total≈45)
    [11:31:01]     chembl target=CHEMBL3311 pref=Low-density lipoprotein receptor compounds=4
    [11:31:01] PHASE 1-3: GABRB2
    [11:31:02]     pubmed thyroid hits: 5 (cancer total≈5)
    [11:31:04]     chembl target=CHEMBL1920 pref=Gamma-aminobutyric acid receptor subunit beta-2 compounds=1
    [11:31:04] PHASE 1-3: B3GNT3
    [11:31:05]     pubmed thyroid hits: 0 (cancer total≈13)
    [11:31:06]     chembl target=CHEMBL3325305 pref=N-acetyllactosaminide beta-1,3-N-acetylglucosaminyltransferase 3 compounds=0
    [11:31:06] PHASE 1-3: PTPRE
    [11:31:07]     pubmed thyroid hits: 3 (cancer total≈3)
    [11:31:08]     chembl target=CHEMBL4850 pref=Receptor-type tyrosine-protein phosphatase epsilon compounds=0
    [11:31:09] PHASE 4: payload size = 16.8 KB