# Moderna/Merck Blinded Validation Protocol Scaffold

## Objective

Test whether CROSS-Neo KG-GA improves fixed-slot individualized neoantigen selection versus sponsor incumbent scoring on blinded candidate pools.

## Data Flow

1. Sponsor provides candidate table without outcome labels or incumbent rank if possible.
2. CROSS-Neo team freezes parser, model version, config hash, and score columns before processing.
3. CROSS-Neo returns per-candidate score, per-patient top-34 list, reason codes, uncertainty flags, and assay-control annotations.
4. Sponsor unblinds labels and calculates pre-specified endpoints.
5. Joint review decides whether to expand to wetlab or license option.

## Minimum Input Schema

| column | required | note |
|---|---|---|
| patient_id | yes | de-identified |
| mutant_peptide | yes | 8-15 aa for MHC-I path; longer class-II can be separated |
| hla_allele | yes | high-resolution preferred |
| wt_peptide | preferred | enables mutant-vs-WT specificity filters |
| expression / RNA support | preferred | optional but valuable |
| variant context | preferred | SNV/indel/fusion/source metadata if available |
| sponsor_score | optional | held out until after lock if possible |
| wetlab_label | yes, withheld | binary/ordinal response, presentation, or selection outcome |

## Success Criteria

| module                  | required_from_partner                                                                                                   | our_commitment                                              | pass_signal                                                                             |
|:------------------------|:------------------------------------------------------------------------------------------------------------------------|:------------------------------------------------------------|:----------------------------------------------------------------------------------------|
| Blinded candidate input | patient_id, mutant peptide, optional WT peptide, HLA, expression, variant context, current internal scores if available | Freeze parser and score generation before label unblinding  | >=95% candidate parse success and reproducible score file hash                          |
| Primary endpoint        | binary or ordinal immunogenicity / selection / wetlab labels withheld until lock                                        | Top-34 ranking and calibrated score per patient             | AUPRC or top-34 hit-rate lift versus incumbent sponsor baseline                         |
| Slot-value endpoint     | per-patient vaccine slot budget, ideally 20-34                                                                          | Rank candidates under a V940-like slot constraint           | More validated candidates in fixed 34-slot payload without loss of HLA/source diversity |
| False-positive control  | decoys, WT peptides, benign/self-like candidates, failed assay records                                                  | Report claim-safe and production-priority scores separately | Lower false-positive burden at matched sensitivity                                      |
| Explainability          | no extra confidential biology required                                                                                  | Per-candidate reason codes, uncertainty, and audit flags    | Sponsor scientist can adjudicate top candidates without opening model internals         |
| Wetlab translation      | assay format, HLA cell line/APC constraints, readout thresholds                                                         | Mutant/WT/decoy control map and top-N plate-ready export    | Assay-ready list accepted by wetlab operations                                          |

## Statistical Readout

- Primary: per-candidate AUPRC and per-patient top-34 hit rate versus incumbent baseline.
- Secondary: AUROC, top-10/top-20/top-34 precision, calibration, false-positive burden at matched sensitivity.
- Stratification: tumor type, HLA family, mutation class, expression availability, payload slot budget.
- Hard rule: no model fitting, threshold tuning, or weight update on sponsor labels before primary analysis.

## Decision Gates

| stage                      | status       | deliverable                                                    | exit_criterion                                                                                      | share_level              |
|:---------------------------|:-------------|:---------------------------------------------------------------|:----------------------------------------------------------------------------------------------------|:-------------------------|
| 0. Internal lock           | done         | Frozen local KG-GA comparison package                          | AUPRC 0.933 on 2,715-row board with explicit caveats                                                | internal only            |
| 1. IP hygiene              | next         | Provisional patent / invention disclosure                      | Claims cover KG-guided evolutionary neoantigen selection, audit gates, and assay-control generation | attorney / tech transfer |
| 2. Non-confidential teaser | ready        | One-page BD teaser with masked method details                  | No code, no full feature list, no exact candidate tables                                            | pre-NDA                  |
| 3. CDA/NDA                 | next         | Mutual NDA and blinded benchmark data-use terms                | Sponsor cannot use submitted scores to train without written license                                | legal                    |
| 4. Blinded validation      | design-ready | Frozen scoring container/API run on sponsor blinded candidates | Pre-specified top-34 enrichment/AUPRC lift vs incumbent selection stack                             | NDA                      |
| 5. Option/license          | conditional  | Field-limited oncology neoantigen selection license            | Payment tied to validation pass, clinical program use, and milestones                               | BD/legal                 |

