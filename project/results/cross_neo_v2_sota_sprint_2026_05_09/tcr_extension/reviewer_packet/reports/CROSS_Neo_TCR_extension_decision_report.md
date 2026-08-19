# CROSS-Neo-TCR Extension Decision Report

Decision: **PROMOTE_TO_DIAGNOSTIC_CASE_STUDY** and **PROMOTE_TO_WETLAB_PRIORITIZATION_TOOL** for the current sprint. Hold main-method promotion until paired TCR benchmarks, stronger TCR sequence models, and structure validation are complete.

## Required Answers

1. Paired TCR alpha/beta rows: **64,037** in the local TCR registry.
2. Peptide-HLA-TCR labeled rows: **98,718** with any TCR sequence; **56,329** with paired alpha/beta.
3. Cancer neoantigen overlap: **2,016/2,715** CROSS-Neo rows have some TCR-resource overlap; **101** exact paired TCR-pMHC rows; **1,876** rows have cancer-context evidence.
4. Top-k/ranking improvement: internal/exact/near pilots improve with TCR evidence. Conservative exact holdout: AUPRC 0.588 -> 0.666 (delta +0.078), AUROC delta +0.094; near holdout: AUPRC 0.393 -> 0.558 (delta +0.165), AUROC delta +0.133.
5. Structure over sequence-only: **not yet supported**. Structure files currently contain missingness/QC features unless PDB/template evidence exists. However, external sequence model escalation is now supported by pMTnet and TEPCAM pilots.
6. Mutant-WT interface delta: **not yet supported**. Needs parsed mutant/WT TCR-pMHC structures.
7. Source-heldout rescue: partial. Conservative NEPdb source-heldout readout: AUPRC 0.233 -> 0.259 (delta +0.026), AUROC delta +0.050. Raw TCR evidence is stronger but treated as leakage-prone.
8. False positives: case audit files identify main-high/TCR-low and TCR-high-only negatives, but explanation remains diagnostic until structures or external TCR assays support it.
9. Main claim or supplement: **supplement/diagnostic branch now**, not the main CROSS-Neo ranking claim.
10. Wetlab candidates: prioritize rows in `tcr_wetlab_candidate_prioritization.tsv` and the de-duplicated `tcr_wetlab_candidate_prioritization_unique_pmhc.tsv` with high conservative TCR-augmented score, positive delta, exact paired or cancer-context TCR evidence, and low pathogen-only dependence.

## Decoy Recognition Pilot

- Handcrafted paired positive versus shuffled-TCR decoy pilot best result: pmhc_only AUPRC 0.500, AUROC 0.500. This indicates the simple handcrafted TCR sequence features are not sufficient for de novo cognate-recognition prediction.
- External pMTnet pilot on positive public TCR-pMHC rows versus same peptide-HLA shuffled-TCR decoys: AUPRC 0.749, AUROC 0.717. This supports escalating to external TCR-specific predictors, still under synthetic-decoy/public-source claim limits.
- Same-pilot external comparison:
  - pMTnet, TCRbeta + peptide + HLA: AUPRC 0.7486, AUROC 0.7167.
  - TEPCAM, TCRbeta + peptide: AUPRC 0.7532, AUROC 0.7156.
  - ERGO-II-vdjdb beta-only/unknown-VJ adapter, TCRbeta + peptide + coarse MHC: AUPRC 0.5587, AUROC 0.4870.
  - PanPep zero-shot, TCRbeta + peptide: AUPRC 0.5248, AUROC 0.4855.
- Runtime conclusion: use pMTnet as P0 HLA-aware expert and TEPCAM as P1 sequence-only expert. Keep ERGO-II for richer alpha/beta/VJ rows rather than beta-only unknown-VJ mode, and keep PanPep as diagnostic/negative-control unless few-shot improves.

## Fast External Expert Challenge Benchmark

Six 100-positive challenge panels were built from public positives plus same-panel shuffled-TCR decoys: all public, cancer-context, paired-TCR, VDJdb, McPAS-TCR, and 10x-derived rows.

- pMTnet alone: AUPRC 0.7479, AUROC 0.7255.
- TEPCAM alone: AUPRC 0.6878, AUROC 0.6970.
- Mean pMTnet/TEPCAM diagnostic ensemble: AUPRC 0.7841, AUROC 0.7652.
- Peptide-HLA GroupKFold logistic pMTnet/TEPCAM ensemble: AUPRC 0.7803, AUROC 0.7683.

This supports a practical two-expert TCR branch. It still does not support a main-method or clinical claim because the benchmark uses synthetic decoys and public-source data.

## Wetlab Candidate External Expert Rescoring

The top 25 wetlab PMHC candidates were checked for exact peptide-HLA TCR rows and rescored with pMTnet/TEPCAM where exact TCR evidence was available.

- HMTEVVRHC / HLA-A*02:01: 4 exact TCR rows scored; external mean 0.7121, max 0.8029; remains the top external-adjusted candidate.
- GADGVGKSAL / HLA-C*08:02: 5 exact TCR rows scored; external mean 0.5703, max 0.6790; remains second after external adjustment.
- Most other top-priority PMHC candidates lacked exact modelable TCR rows, so their ranking remains based on pMHC/TCR-evidence metadata rather than external TCR expert scores.

Wetlab prioritization implication: prioritize HMTEVVRHC/HLA-A*02:01 first for TCR-aware follow-up; treat GADGVGKSAL/HLA-C*08:02 as a second external-supported candidate; keep unscored candidates as pMHC-led candidates, not TCR-recognition-supported candidates.

## MD Escalation Queue

Prediction-fragile, high-value candidates can be escalated to MD, but only as a diagnostic layer after credible structures are available.

- P0 TCR-pMHC MD after chain recovery: HMTEVVRHC / HLA-A*02:01 and GADGVGKSAL / HLA-C*08:02.
- P0 existing structure breakthrough: experimentally solved PDB complexes were found for both P0 peptides. HMTEVVRHC / HLA-A*02:01 has ready pilots from 6VRN, 6VRM, 6VQO, and 7RM4; GADGVGKSAL / HLA-C*08:02 has ready pilots from 6UON.
- OpenMM dry repair/minimization sanity check succeeded for representative P0 complexes 6VRN and 6UON.
- OpenMM local smoke dynamics succeeded for both representatives using short vacuum smoke tests. Final peptide RMSD was 0.0091 nm for 6VRN/HMTEVVRHC and 0.0144 nm for 6UON/GADGVGKSAL. This verifies the runner/trajectory/analysis plumbing, not production stability.
- P1 pMHC bulge/stability MD: long class-I peptides such as MAWSLGVLVALPFPL, FSLVFLVYSVFKNNV, and MSSLAATTFHWKKCR across their candidate HLA contexts.
- P1 pMHC uncertainty MD plus TCR search: high-priority candidates without exact modelable TCR rows, including MSFSHLFYL / HLA-A*02:01 and VLMMPFSIV / HLA-A*02:01.
- Prepared MD job stubs: top-20 queue, 20 jobs, 3,300 ns planned aggregate production time if every stub is eventually run.
- Current blocker: local OpenMM/MDTraj/PDBFixer are now available, but only CPU/Reference OpenMM platforms are present. Production MD should run on a CUDA-enabled GPU pod. GROMACS and MDAnalysis remain unavailable locally.

Claim boundary: MD can support interface diagnostics, conformational-stability triage, and mutant-WT cross-reactivity hypotheses. It cannot prove immunogenicity by itself.

## Figure Package

Figure package generated under `project/results/cross_neo_v2_sota_sprint_2026_05_09/tcr_extension/figures/`:

- `fig_tcr1_recognition_gap`: presentation vs recognition gap.
- `fig_tcr2_architecture`: pMHC, counterfactual, TCR sequence, structure, and MoE gate.
- `fig_tcr3_structure_pipeline`: TCR-pMHC structure workflow and bottleneck.
- `fig_tcr4_mutant_wt_interface_delta`: mutant vs WT interface-delta schematic.
- `fig_tcr5_tcr_available_subset_performance`: pMTnet/TEPCAM/ensemble performance.
- `fig_tcr6_case_studies`: wetlab candidate rescoring.
- `fig_tcr7_claim_boundary`: allowed versus forbidden claims.

## Claim Boundary

- Do not replace the main pMHC neoantigen model with the TCR branch.
- Do not transfer pathogen epitope labels to cancer neoantigens.
- Do not claim clinical utility or universal TCR-aware prediction.
- Safe current claim: optional TCR evidence/diagnostic layer plus wetlab prioritization scaffold.
