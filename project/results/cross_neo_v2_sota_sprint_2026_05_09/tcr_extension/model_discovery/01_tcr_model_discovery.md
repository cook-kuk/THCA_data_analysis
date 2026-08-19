# CROSS-Neo-TCR Model Discovery Sprint

Date: 2026-05-09

## Bottom Line

We found enough additional model surface to keep pushing. The immediate model queue should not be structure-first; it should be **pMTnet + TEPCAM first**, then **ERGO-II adapter**, then **NetTCR-2.2/TSpred after CDR1/CDR2 reconstruction**.

## Local Runtime Discovery

- NetTCR-2.2 local repo: `/data/thca/_tmp/relocated_2026_05_09/NetTCR-2.2`
- NetTCR-2.2 pretrained TFLite weights found: `True`
- ColabFold batch executable: `/home/seungho/personal/THCA_data_analysis/.venv/bin/colabfold_batch`
- TensorFlow available: `True`
- PyTorch available: `True`
- Transformers available: `True`
- ESM available: `True`
- Stitchr/ANARCI/tcrdist3 available: `stitchr=False`, `anarci=True`, `tcrdist3=False`

## NetTCR-2.2 Sanity Run

- Small example predictions completed: n=100, positives=18
- Example AUPRC=0.932, AUROC=0.980
- Output: `/home/seungho/personal/THCA_data_analysis/project/results/cross_neo_v2_sota_sprint_2026_05_09/tcr_extension/model_discovery/nettcr_sanity/nettcr_small_example_predictions.csv`
- This is only a runtime sanity check on repository example data, not a CROSS-Neo validation result.

## Same-Pilot External Model Comparison

| Model | Inputs | n | AUPRC | AUROC | Read |
|---|---|---:|---:|---:|---|
| pMTnet | TCRbeta + peptide + HLA | 355 | 0.7486 | 0.7167 | P0 runnable; pilot claim-limited |
| PanPep | TCRbeta + peptide | 374 | 0.5248 | 0.4855 | P1 runnable; no HLA in score |
| TEPCAM | TCRbeta + peptide | 374 | 0.7532 | 0.7156 | P1 runnable; no HLA in score |
| ERGO-II-vdjdb | TCRbeta + peptide + coarse MHC | 374 | 0.5587 | 0.4870 | P1 runnable; adapter is beta-only/unknown-VJ |

These are fast shuffled-decoy pilot numbers, not clean external benchmark claims.

## Fast Multi-Panel Challenge Benchmark

| Model | n | AUPRC | AUROC |
|---|---:|---:|---:|
| mean_pMTnet_TEPCAM | 1191 | 0.7841 | 0.7652 |
| logreg_groupcv_pMTnet_TEPCAM | 1177 | 0.7803 | 0.7683 |
| pMTnet | 1177 | 0.7479 | 0.7255 |
| TEPCAM | 1191 | 0.6878 | 0.6970 |

The no-training mean pMTnet/TEPCAM ensemble is currently the best practical diagnostic score.

## P0 Models

| Model | Why P0 | Immediate blocker | Local status |
|---|---|---|---|
| NetTCR-2.2 | TCR-available subset baseline/expert only | CROSS-Neo registry lacks CDR1/CDR2; need V-gene CDR reconstruction or full-chain TCR import | local_ready:/data/thca/_tmp/relocated_2026_05_09/NetTCR-2.2 |
| ERGO-II | strong next practical model because CDR3-only path matches current registry | CROSS-Neo input adapter and strict split evaluation | local_ready:/data/thca/tcr_model_repos/ERGO-II |
| pMTnet | diagnostic peptide-HLA-TCR score; useful for beta-only rows | strict leakage-filtered benchmark and score calibration | local_ready:/data/thca/tcr_model_repos/pMTnet |

## P1 Models

| Model | Role | Immediate blocker |
|---|---|---|
| TSpred | paired TCR subset benchmark; seen/unseen epitope split baseline | requires CDR1/CDR2 reconstruction for registry-scale use |
| PanPep | unseen-peptide stress test and negative-control/diagnostic baseline | weak zero-shot pilot; try few-shot only if compatible support labels exist |
| TEPCAM | sequence-only TCRbeta-peptide expert; compare to pMTnet | no HLA input; needs strict source/epitope-heldout validation |
| TCR-H | explainable low-compute hard-split baseline | clone and map feature script |
| TCRen | structure-guided case diagnostics, especially unseen epitopes | requires modeled/experimental TCR-pMHC structures |
| tFold-TCR | structure feature generator, not recognition classifier alone | not installed; full-chain sequences missing for many jobs |
| TCRdock | geometry/contact diagnostics | not installed locally |
| AlphaFold-Multimer/ColabFold | structure feature generator and parser validation | full-chain sequences missing for most registry rows |
| UniPMT | unified pMHC/TCR expert benchmark | install and verify checkpoint/data format |

## Expanded Model Inventory

| Model | Family | Input | HLA modeled | Paired alpha/beta | Priority | Local status |
|---|---|---|---|---|---|---|
| NetTCR-2.2 | sequence paired-chain CNN | peptide + CDR1/2/3 alpha + CDR1/2/3 beta | no explicit HLA | required/preferred | P0 | local_ready:/data/thca/_tmp/relocated_2026_05_09/NetTCR-2.2 |
| ERGO-II | sequence LSTM/attention TCR-peptide | peptide + CDR3 beta; optional alpha, V/J, MHC, T-cell type | optional MHC feature | not required | P0 | local_ready:/data/thca/tcr_model_repos/ERGO-II |
| pMTnet | sequence beta-chain peptide-HLA model | CDR3 beta + peptide + HLA | yes | no, beta-only | P0/P1 | local_ready:/data/thca/tcr_model_repos/pMTnet |
| TSpred | paired-chain CNN + reciprocal attention | paired alpha/beta TCR sequence + epitope | no explicit HLA in headline task | yes | P1 | local_ready:/data/thca/tcr_model_repos/TSpred |
| PanPep | meta-learning peptide-TCR | peptide + TCR sequence, commonly CDR3 beta in public usage; extensions reported for alpha/alpha-beta | no explicit HLA in original core task | not required in original core task | P1 | local_ready:/data/thca/tcr_model_repos/PanPep |
| TEPCAM | cross-attention + multi-channel convolution | TCR sequence + peptide | no explicit HLA | verify code; likely sequence-level TCR/epitope | P1/P2 | local_ready:/data/thca/tcr_model_repos/TEPCAM |
| TEINet | deep TCR-epitope sequence model | CDR3 beta + epitope | no explicit HLA | no | P2 | not_found |
| TCR-H | SVM physicochemical CDR3beta-epitope | CDR3 beta + epitope | no explicit HLA | no | P1/P2 | local_ready:/data/thca/tcr_model_repos/TCR-H |
| MixTCRpred | dual-alpha-aware TCR-epitope model | paired/dual alpha TCR context + beta + epitope depending dataset | not central | yes/optional depending mode | P2 | not_found |
| TCRen | structure/statistical potential | TCR-pMHC model/structure or homology model + peptide context | yes via pMHC structure | yes for TCR complex | P1/P2 | not_found |
| TCRLens | structure-aware EGNN | TCR-pMHC-I structure/interface graph | yes, class I focus | yes | P2 | not_found |
| tFold-TCR | structure prediction | full TCR alpha/beta + MHC + beta2M + peptide sequences | yes via structure | yes | P1 | not_found |
| TCRdock | TCR-pMHC structure/docking geometry | peptide-HLA + TCR alpha/beta V/J/CDR3 or structures | yes | yes | P1 | not_found |
| AlphaFold-Multimer/ColabFold | general complex structure prediction | full chain FASTA: TCR alpha/beta, MHC, beta2M, peptide | yes | yes for TCR-pMHC | P1 | local_ready:/home/seungho/personal/THCA_data_analysis/.venv/bin/colabfold_batch |
| UniPMT | unified peptide-MHC-TCR graph/multitask | peptide + MHC + TCR features depending task | yes | verify, often CDR3 beta heavy | P1/P2 | not_found |
| UnifyImmun | unified cross-attention pHLA/pTCR | peptide + HLA and/or peptide + TCR | yes for pHLA branch | verify | P2 | not_found |
| ImRex | interaction-map CNN | CDR3 beta + epitope | no explicit HLA | no | P2 | not_found |
| TITAN | bimodal attention | TCR sequence + epitope | no explicit HLA | not central | P2 | not_found |

## Hard Decision

1. **Run pMTnet and TEPCAM before structure-first work**, because current rows already have peptide and CDR3 beta, and pMTnet additionally uses HLA.
2. **Move ERGO-II next**, because runtime is now fixed but CROSS-Neo input mapping still needs a clean adapter.
3. **Do not throw away NetTCR-2.2/TSpred**. They are concrete local assets, but need a CDR1/CDR2 reconstruction layer from V genes or full-chain TCR imports.
4. **Do not make structure the next bottleneck**. tFold/TCRdock/TCRen/TCRLens become valuable after full-chain sequence recovery; today they are mostly blocked.
5. **Use all external models as optional TCR expert scores**, not as replacements for the pMHC ranker.

## Sources

- NetTCR-2.2: https://github.com/mnielLab/NetTCR-2.2
- TSpred: https://academic.oup.com/bioinformatics/article/40/8/btae472/7721043
- PanPep: https://www.nature.com/articles/s42256-023-00619-3
- TEPCAM: https://pubmed.ncbi.nlm.nih.gov/37983648/
- TCR-H: https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2024.1426173/full
- TCRen: https://www.nature.com/articles/s43588-024-00653-0
- TCRLens: https://www.lifescience.net/publications/1944865/tcrlens-structure-aware-equivariant-graph-learning/
- tFold-TCR: https://github.com/TencentAI4S/tfold
- Comparative model dependency warning: https://pmc.ncbi.nlm.nih.gov/articles/PMC10152969/
