# Web Search Prompt: Better Neoantigen Algorithms

You are a research agent. Find better algorithmic directions for a paper-grade neoantigen immunogenicity predictor, using recent web literature and primary sources. Do not generate manuscript prose. Produce a method-design report with citations, tables, and concrete next experiments.

## Project Context

We are building a neoantigen prioritization method for MHC class I and eventually class II. The current local benchmark is ITSNdb MHC-I peptide immunogenicity.

Current key finding:

- Full no-overlap class I: `Wave8_TCR_SelfSim_full` looks best, AUROC ~0.735, but is exact-reference sensitive.
- Strict class I after removing TCR/self exact hits: `Structure_LR` is the clean locked leader, AUROC ~0.679.
- `MHCflurry` is useful but public-pretrained / training-overlap unresolved, not a clean baseline.
- ESMFold pseudo-complex structure features are weak alone; structure is useful as fallback/uncertainty, not a full reference substitute.
- Quantum should not be abandoned: internal strict-set pilot found `quantum_kernel_no_anchor_gamma1.0` AUROC ~0.785, AUPRC ~0.603, but it is small-n internal CV and needs external/time-split validation.

## Algorithms Already Tried Locally

Compare new methods against these categories and explain what they improve:

| category | examples tried | current read |
|---|---|---|
| Public MHC-I predictors | MHCflurry, NetMHCpan, MHCnuggets, NetMHCstabpan | useful comparators; public-training contamination risk |
| Public immunogenicity/presentation models | BigMHC_IM, PRIME, DeepImmuno, TransPHLA, TSCAPE/TITANiAN | mixed; often inflated on in-master, weaker on no-overlap |
| Local sequence/PLM model | ESM2_Bayesian | fails no-overlap; memorization/leakage failure baseline |
| Local biophysics baseline | Structure_LR, RF_biophys | Structure_LR is cleanest strict anchor |
| Stacking/ensembles | Stack_mean, Stack_median, Stack_LR, OOF stacks | often no robust no-overlap win |
| Test-time adaptation | TENT | weak on no-overlap |
| Multitask/deep variants | MultiTask A/B/C/D/full, LoRA/GroupDRO/MIRO/MoLE attempts | not yet convincing under strict split |
| Quantum | GP_quantum, VQC, W7A_QK_only, W7A_full, quantum-kernel fallback | keep; best new pilot is quantum-kernel + structure/no-anchor |
| Retrieval/reference | Wave8 TCR motif, TCR+SelfSim, self-similarity | strong when references exist, exact-hit sensitive |
| Structure proxy | ESMFold peptide+linker+HLA-pseudo pseudo-complex | folded successfully, weak standalone, useful as uncertainty/geometry branch |

## Search Goals

Find algorithm ideas that can beat or improve on:

1. `Structure_LR` as the clean no-reference baseline.
2. `Wave8` when reference evidence exists but must be leakage-flagged.
3. `quantum_kernel + structure` as a no-reference research candidate.
4. Public-pretrained tools while controlling for training contamination.

Prioritize methods from 2023-2026. Use primary sources whenever possible: conference proceedings, official benchmark/challenge pages, PubMed/Nature/Cell/Science pages, arXiv only if no peer-reviewed source exists.

## Venues And Sources To Search

Search at least these areas:

- ML conferences: NeurIPS, ICML, ICLR, CVPR, KDD, AAAI, IJCAI.
- Bio/medical ML: RECOMB, ISMB, MICCAI, MLSB, ML4H, CASP.
- Immunology/neoantigen benchmarks: IMMREP23, IMMREP25, TESLA/Wells benchmark, IEDB automated benchmarks, CEDAR, NEPdb, CAGI, DREAM, Kaggle TCR/pMHC challenges.
- Protein modeling: AlphaFold 3, Boltz, ESM3, SaProt, ProSST, RFdiffusion, protein structure diffusion/flow matching, protein retrieval-augmented models.
- Neoantigen/pMHC methods: BigMHC, HLApollo, PRIME, MixMHCpred/MixMHC2pred, NetMHCpan/NetMHCIIpan, ImmunoStruct, NeoaPred, ImmuScope, pVACtools/pVACview.
- TCR-pMHC recognition: EPACT, TITAN/TITANiAN, TCR-ESM, TouCAN, PanPep, unseen-pMHC benchmarks.

## Search Queries

Use and expand these queries:

```text
2026 neoantigen immunogenicity prediction benchmark leakage contamination
2025 neoantigen prediction structure sequence multimodal immunogenicity
2024 peptide MHC immunogenicity model structure deep learning benchmark
IMMREP23 unseen pMHC TCR specificity negative bias leakage
IMMREP25 unseen peptides TCR pMHC challenge methods
MHC-Bench training overlap peptide MHC predictor benchmark
HLApollo neoantigen presentation transformer source protein features
ImmunoStruct sequence structure biochemical peptide MHC immunogenicity
ImmuScope MHC class II weak supervision multiple instance learning
SaProt ProSST protein language model structure tokens immunology application
protein retrieval augmented language model neoantigen peptide MHC
positive unlabeled learning immunogenicity neoantigen negative labels
conformal prediction OOD calibration biomedical protein model
quantum kernel protein peptide immunogenicity classification
quantum machine learning peptide MHC neoantigen prediction
geometric deep learning peptide MHC TCR-facing surface immunogenicity
AlphaFold3 Boltz peptide MHC structure immunogenicity prediction
```

## Required Output

Produce a Markdown report with these sections:

### 1. Executive Decision

Give a ranked list of 5-8 algorithm directions:

- `adopt_now`
- `pilot_next`
- `watch_only`
- `reject`

Each item must say why it is better than our current attempts, what data it needs, and what leakage risk it introduces.

### 2. Literature Table

Create a table:

| year | method/paper | source URL | task | class I/II | model idea | benchmark | leakage/overlap handling | relevance to our repo |

Do not include uncited claims.

### 3. Competition/Benchmark Lessons

Summarize lessons from IMMREP23/25, TESLA, IEDB benchmarks, CEDAR/NEPdb, MHC-Bench or related datasets:

- seen vs unseen pMHC behavior
- negative-set bias
- public training contamination
- time-split feasibility
- allele/study/patient splits
- AUPRC/top-k vs AUROC

### 4. Better Algorithm Designs

Propose concrete designs that could become a paper method:

1. **Leakage-safe retrieval branch**
   - train-only retrieval index
   - exact/near match flags
   - no hidden reference leakage

2. **Quantum-structure fallback branch**
   - angle-encoded fidelity kernel or other quantum kernel
   - structure proxy features
   - no public-pretrained score features
   - evaluate with HLA/time/study holdout

3. **Multimodal PLM branch**
   - mut peptide, WT peptide, HLA pseudo-sequence, source-protein window
   - frozen ESM/SaProt/ProSST/ESM3-style embeddings
   - small calibrated head

4. **pMHC geometry branch**
   - AlphaFold3/Boltz/ESMFold-derived ensemble
   - TCR-facing exposure, anchor/non-anchor separation, peptide bulge, confidence
   - explicit uncertainty, not single-structure truth

5. **PU/MIL learning branch**
   - positive-unlabeled immunogenicity labels
   - patient/vaccine candidate bag learning
   - class-II CD4 module

6. **Calibration/OOD branch**
   - conformal abstention
   - HLA/study/assay/peptide-length OOD flags
   - top-k precision at vaccine-relevant k

### 5. Evaluation Contract

Design a benchmark table that includes:

- exact peptide-HLA holdout
- near peptide similarity holdout
- source-protein window holdout
- patient/study holdout
- HLA allele/supertype holdout
- unseen pMHC split
- time-split holdout
- public-corpus overlap audit against MHCflurry, NetMHCpan, BigMHC, PRIME/MixMHCpred, CEDAR, IEDB, TESLA, NEPdb

Metrics:

- AUPRC
- top-k precision/enrichment
- AUROC as secondary
- calibration error
- abstention coverage
- per-HLA and per-study breakdown

### 6. Implementation Plan

Give a fast plan:

- 24-hour pilot
- 3-day pilot
- 1-week paper-grade validation

For each, list exact artifacts to produce:

- input TSVs
- overlap audit tables
- method comparison tables
- figures
- failure-mode tables

## Strict Rules

- Do not recommend using MHCflurry/NetMHCpan/BigMHC scores as training features if those tools are benchmark comparators.
- Do not call MHCflurry a clean external baseline unless training overlap is proven absent.
- Do not pool class I and class II AUROCs unless the endpoint and benchmark are comparable.
- Do not accept exact TCR/self retrieval gains without exact-match and near-match removal stress tests.
- Do not optimize only AUROC; prioritize AUPRC and top-k precision.
- Mark every claim as one of:
  - `confirmed by source`
  - `inference from source`
  - `hypothesis for pilot`

## Desired Final Recommendation

End with a single recommended method name and architecture. The recommendation should be conservative enough for reviewers but novel enough for top-tier ML/bioinformatics:

Example target framing:

> A contamination-controlled, class-aware, retrieval-flagged, quantum-structure fallback model for neoantigen prioritization under unseen-pMHC shift.

Also provide a one-paragraph reviewer defense explaining why this is not just another MHC binding predictor.
