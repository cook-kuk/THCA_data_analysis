# Neoantigen Method Blueprint - 2026-05-09

## Bottom Line

The current result is not a clean SOTA claim. It is a leakage/generalization story that can become a paper-grade method if the method is redesigned around contamination control, epitope-shift robustness, and class I/class II coverage.

Proposed method: **CLEAN-Neo** — contamination-controlled, leakage-aware, evidence-calibrated neoantigen prioritization under epitope shift.

## Current Evidence To Preserve

Class I ITSNdb no-overlap:

| method | n / pos | AUROC | AUPRC | interpretation |
|---|---:|---:|---:|---|
| Wave8_TCR_SelfSim_full | 106 / 33 | 0.735 | 0.598 | apparent leader, exact TCR-reference sensitive |
| MHCflurry | 106 / 33 | 0.668 | 0.553 | public-pretrained upper-bound comparator, training overlap unresolved |
| Structure_LR | 103 / 32 | 0.653 | 0.489 | cleanest local anchor |
| BigMHC_IM | 106 / 33 | 0.626 | 0.380 | public-pretrained, in-master inflated |
| ESM2_Bayesian | 103 / 32 | 0.411 | 0.274 | leakage-failure baseline |

Strict class I subset after removing TCR/self exact hits and requiring HLA pseudo-sequence:

| method | n / pos | AUROC | AUPRC | interpretation |
|---|---:|---:|---:|---|
| Structure_LR | 89 / 21 | 0.679 | 0.435 | strict-set leader |
| MHCflurry | 89 / 21 | 0.657 | 0.488 | still caveated |
| ESMFold min peptide pLDDT | 89 / 21 | 0.644 | 0.313 | weak univariate structure proxy |
| Wave8_TCR_SelfSim_full | 89 / 21 | 0.632 | 0.403 | drops after exact-reference removal |
| ESMFold 3D LR CV | 89 / 21 | 0.569 | 0.279 | not enough as replacement |

Decision: **do not claim ESMFold/AlphaFold proxy replaces missing TCR evidence yet.** Use structure as a fallback feature and uncertainty signal; keep Structure_LR as the clean local baseline.

Class II status: no local comparable class-II peptide immunogenicity benchmark exists yet. Class II must be built as a separate task, not mixed into the ITSNdb class-I table.

## What Recent Work Says

1. **The field is leakage-sensitive.** IMMREP23 found reasonable performance for seen pMHC targets but near-random behavior for unseen peptides, and explicitly identified negative-set leakage/bias as a benchmark problem.
2. **Presentation models are moving beyond binding.** BigMHC and HLApollo use eluted-ligand presentation pretraining and immunogenicity transfer; HLApollo adds source-protein language features and negative-set switching.
3. **Structure/multimodal is now credible, but needs calibration.** ImmunoStruct integrates sequence, structure, and biochemical features for class-I pMHC immunogenicity. NeoaPred-style mutant-vs-WT surface/structure deltas are aligned with the biological problem.
4. **Class II is no longer optional.** ImmuScope uses weakly supervised MIL for MHC-II/CD4 antigen presentation and immunogenicity over broad DR/DQ/DP coverage.
5. **TCR prediction remains the weak link.** Recent TCR-pMHC challenges and models show that unseen-epitope generalization is still unsolved; exact TCR retrieval must be treated as evidence with flags, not a hidden feature.

## CLEAN-Neo Architecture

### Inputs

- Mutant peptide and wild-type counterpart.
- HLA class I allele and class II DR/DQ/DP alleles.
- Source protein window, mutation position, flanking residues, expression, clonality/VAF when available.
- Optional patient-level bag of candidates for vaccine or ICI response cohorts.
- Optional TCR repertoire or reference TCR-pMHC retrieval evidence, always leakage-flagged.

### Branches

| branch | role | implementation target |
|---|---|---|
| Delta sequence branch | learn `mut - WT` foreignness and anchor/non-anchor changes | frozen ESM2/SaProt/ProSST-style encoder + small head |
| Presentation branch | predict MHC-I/MHC-II presentation without using benchmarked tools as features | HLA pseudo-sequence cross-attention, flanks/source-protein context |
| Structure branch | model TCR-facing geometry and pMHC confidence | ESMFold/Boltz/AlphaFold-derived ensemble features, peptide exposure, bulge, anchor geometry, pLDDT, WT-mut deltas |
| Retrieval branch | use evidence without silent memorization | HLA/motif/embedding/source-gene neighbors from train-only index; exact/near matches become separate evidence flags |
| TCR latent branch | learn recognition potential when TCR data are absent | auxiliary pretraining on public TCR-pMHC pairs, evaluated with unseen-pMHC split |
| Patient MIL branch | translate peptide scores to patient/vaccine ranking | attention MIL over candidate bags with clinical outcome weak labels |
| Calibration branch | prevent overconfident deployment | group-wise temperature scaling + conformal abstention by HLA, assay, study, peptide length, OOD distance |

### Training

- Multi-task objectives: MHC-I presentation, MHC-II presentation, immunogenicity, TCR-pMHC specificity, patient-level response.
- Positive-unlabeled risk or reliable-negative selection for immunogenicity, because most negatives are assay-limited unobserved candidates.
- GroupDRO/domain-balanced sampling by study, HLA supertype, assay type, peptide length, tumor type, and source database.
- No MHCflurry/NetMHCpan/BigMHC scores as main-model features when those tools are also benchmark comparators.

## Benchmark Contract

Report every result under these splits:

| split | purpose |
|---|---|
| exact peptide-HLA holdout | removes direct memorization |
| near peptide similarity holdout | removes motif-level memorization |
| source-protein window holdout | removes same-gene/window leakage |
| study/patient holdout | removes cohort-specific shortcuts |
| HLA allele/supertype holdout | tests rare-allele transfer |
| unseen pMHC/TCR holdout | tests the IMMREP failure mode |
| time-split holdout | protects against public-database training overlap |
| strict public-tool overlap audit | excludes rows found in available MHCflurry/NetMHCpan/BigMHC/PRIME/CEDAR/IEDB/TESLA/NEPdb training-like corpora where possible |

Primary metrics:

- AUPRC, top-k precision, enrichment at vaccine-relevant k.
- AUROC as secondary.
- Calibration error and abstention coverage.
- Per-HLA and per-study tables, not only pooled metrics.

## Paper-Grade Experiments

| experiment | output | go/no-go criterion |
|---|---|---|
| E0 locked comparison | current class I/class II table | proves why old wave claims are unsafe |
| E1 contamination audit | overlap matrix against public corpora | every comparator labeled clean/caveated |
| E2 CLEAN-Neo class I | strict ITSNdb/TESLA/CEDAR/NEPdb benchmark | beats Structure_LR and caveated public tools on strict AUPRC/top-k |
| E3 structure ablation | ESMFold/Boltz/AlphaFold branch ablation | improves strict external split, not only CV |
| E4 retrieval ablation | exact retrieval on/off, near retrieval on/off | gains survive exact-match removal |
| E5 class II module | DR/DQ/DP benchmark table | establishes class-II companion, not pooled with class I |
| E6 patient MIL | vaccine/ICI patient-level ranking | improves top-k candidate selection and calibration |

## Immediate Implementation Order

1. Freeze current wave comparison and strict ESMFold results as the negative/diagnostic baseline.
2. Build a public-corpus overlap index for CEDAR, IEDB, TESLA, NEPdb, MHCflurry-available training files, NetMHCpan/DTU docs where downloadable, BigMHC release data, PRIME/MixMHCpred release data.
3. Build class-I `strict_v2` and class-II `strict_v1` benchmark tables with explicit overlap flags.
4. Implement CLEAN-Neo v0 with interpretable tabular branches first: delta sequence features, structure proxy, leakage-safe retrieval, calibration/abstention.
5. Only after v0 beats Structure_LR on strict external splits, add frozen PLM/geometric encoders.

## Sources Checked

- IMMREP23 challenge: https://discovery.ucl.ac.uk/id/eprint/10205429/
- IMMREP23 data/Kaggle pointer: https://github.com/justin-barton/IMMREP23/
- MHCflurry docs: https://openvax.github.io/mhcflurry/intro.html
- BigMHC: https://www.nature.com/articles/s42256-023-00694-6
- HLApollo: https://www.nature.com/articles/s41467-024-54887-7
- ImmunoStruct: https://www.nature.com/articles/s42256-025-01163-y
- ImmuScope class-II MIL: https://www.nature.com/articles/s42256-025-01073-z
- ProSST NeurIPS 2024: https://papers.nips.cc/paper_files/paper/2024/hash/3ed57b293db0aab7cc30c44f45262348-Abstract-Conference.html
- AlphaFold 3: https://www.nature.com/articles/s41586-024-07487-w
- RFdiffusion: https://www.nature.com/articles/s41586-023-06415-8
- MHCXAI/MHC-Bench overlap discussion: https://www.nature.com/articles/s42003-024-05968-2
