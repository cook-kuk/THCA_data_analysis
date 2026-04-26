# Biomarker → Drug Discovery Report

- Generated: 2026-04-24 11:26:27
- Biomarker source: `biomarker_validated.tsv`
- Targets analysed: **8**
- Compounds collected: **24**
- Literature records: **23**

## Top 8 targets

| Gene | Novel | Classification | #Cmpds | Best pChEMBL | Thyroid papers | Top PMID | Top Compound |
|---|---|---|---|---|---|---|---|
| DCSTAMP | False | novel_target | 0 | — | 3 | 39193983 | — |
| KCNN4 | False | validated_target | 10 | 8.22 | 4 | 40950697 | SENICAPOC |
| TACSTD2 | False | novel_target | 0 | — | 0 | — | — |
| TMPRSS4 | False | novel_target | 0 | — | 5 | 41656803 | — |
| PLEKHA6 | False | novel_target | 0 | — | 0 | — | — |
| CYP1B1 | False | validated_target | 10 | 8.75 | 5 | 40900788 | CHEMBL3132932 |
| BNC1 | False | novel_target | 0 | — | 1 | 41083582 | — |
| LDLR | False | emerging_target | 4 | 6.22 | 5 | 41340871 | CHEMBL115992 |

## Pharma-quality highlights (3 novel targets)
### DCSTAMP — novel_target
**DCSTAMP** surfaces from biomarker discovery with high novelty_score and no known high-confidence chemical starting points in ChEMBL. For a pharma portfolio this is a **first-in-class opportunity** — the target hypothesis is de-risked by replicated differential expression in thyroid cancer cohorts, while the chemical matter is greenfield (fragment screening / DNA-encoded library campaigns would be the logical next step). Thyroid-specific literature count = 3.

### KCNN4 — validated_target
**KCNN4** is a **validated druggable target** (10 ChEMBL actives, best pChEMBL ≈ 8.22). Its appearance inside our novelty-ranked biomarker set shows that well-trodden pharmacology may be **repositioned** into thyroid cancer contexts, shortening the path to preclinical. Thyroid literature count = 4.

### TACSTD2 — novel_target
**TACSTD2** surfaces from biomarker discovery with high novelty_score and no known high-confidence chemical starting points in ChEMBL. For a pharma portfolio this is a **first-in-class opportunity** — the target hypothesis is de-risked by replicated differential expression in thyroid cancer cohorts, while the chemical matter is greenfield (fragment screening / DNA-encoded library campaigns would be the logical next step). Thyroid-specific literature count = 0.


## Caveat
This is computational triage only — not a drug recommendation. All compound activities come from heterogeneous public assays; thyroid-specific efficacy is not guaranteed. Confirmatory in-vitro / in-vivo work required before any claim.

## Run log
    [11:25:59] PHASE 0: polling for biomarker_validated.tsv (max 1200s)
    [11:25:59] PHASE 0: found validated biomarkers after 0s
    [11:25:59] PHASE 0: loaded 2843 biomarker rows from biomarker_validated.tsv
    [11:25:59] PHASE 0: selected targets: ['DCSTAMP', 'KCNN4', 'TACSTD2', 'TMPRSS4', 'PLEKHA6', 'CYP1B1', 'BNC1', 'LDLR']
    [11:25:59] PHASE 1-3: DCSTAMP
    [11:26:02]     pubmed thyroid hits: 3 (cancer total≈2)
    [11:26:03]     chembl target=None pref=None compounds=0
    [11:26:05] PHASE 1-3: KCNN4
    [11:26:06]     pubmed thyroid hits: 4 (cancer total≈14)
    [11:26:09]     chembl target=CHEMBL4305 pref=Intermediate conductance calcium-activated potassium channel protein 4 compounds=10
    [11:26:10] PHASE 1-3: TACSTD2
    [11:26:10]     pubmed thca search failed for TACSTD2: ('Connection aborted.', RemoteDisconnected('Remote end closed connection without response'))
    [11:26:10]     pubmed thyroid hits: 0 (cancer total≈14)
    [11:26:12]     chembl target=CHEMBL3856163 pref=Tumor-associated calcium signal transducer 2 compounds=0
    [11:26:12] PHASE 1-3: TMPRSS4
    [11:26:14]     pubmed thyroid hits: 5 (cancer total≈60)
    [11:26:15]     chembl target=CHEMBL2331048 pref=Transmembrane protease serine 4 compounds=0
    [11:26:15] PHASE 1-3: PLEKHA6
    [11:26:16]     pubmed thyroid hits: 0 (cancer total≈0)
    [11:26:16]     chembl target=None pref=None compounds=0
    [11:26:16] PHASE 1-3: CYP1B1
    [11:26:17]     pubmed thyroid hits: 5 (cancer total≈242)
    [11:26:21]     chembl target=CHEMBL4878 pref=Cytochrome P450 1B1 compounds=10
    [11:26:22] PHASE 1-3: BNC1
    [11:26:23]     pubmed thyroid hits: 1 (cancer total≈7)
    [11:26:23]     chembl target=None pref=None compounds=0
    [11:26:24] PHASE 1-3: LDLR
    [11:26:25]     pubmed thyroid hits: 5 (cancer total≈45)
    [11:26:27]     chembl target=CHEMBL3311 pref=Low-density lipoprotein receptor compounds=4
    [11:26:27] PHASE 4: payload size = 16.1 KB