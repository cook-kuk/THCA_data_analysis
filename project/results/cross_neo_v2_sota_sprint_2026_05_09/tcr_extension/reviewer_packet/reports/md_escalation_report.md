# CROSS-Neo-TCR MD Escalation Queue

Purpose: identify prediction-fragile, high-value neoantigen candidates that justify expensive MD or enhanced MD.

Important boundary: MD is an escalation diagnostic, not proof of immunogenicity. It should be run only after structure inputs are credible.

## Runtime Reality

- Local GROMACS/OpenMM/MDTraj/MDAnalysis were not detected in the current Python environment.
- Most TCR-pMHC structure jobs are currently blocked by missing full TCR and MHC sequences.
- Therefore the immediate action is a prioritized MD queue plus chain-recovery requirements, not blind production MD.

## Top MD Escalation Candidates

| Rank | Peptide | HLA | Tier | MD score | Fragility | External TCR mean | Protocol |
|---:|---|---|---|---:|---:|---:|---|
| 1 | HMTEVVRHC | HLA-A*02:01 | P0_MD_TCR_pMHC | 0.791 | 0.360 | 0.712 | chain_recovery_then_TCR-pMHC_MD_3x100ns |
| 2 | GADGVGKSAL | HLA-C*08:02 | P0_MD_TCR_pMHC | 0.571 | 0.316 | 0.570 | chain_recovery_then_TCR-pMHC_MD_3x100ns |
| 3 | MAWSLGVLVALPFPL | HLA-B*40:01 | P1_MD_pMHC_bulge | 0.477 | 0.860 | NA | pMHC_bulge_stability_MD_3x50ns |
| 4 | MAWSLGVLVALPFPL | HLA-A*02:01 | P1_MD_pMHC_bulge | 0.477 | 0.860 | NA | pMHC_bulge_stability_MD_3x50ns |
| 5 | FSLVFLVYSVFKNNV | HLA-A*11:01 | P1_MD_pMHC_bulge | 0.450 | 0.840 | NA | pMHC_bulge_stability_MD_3x50ns |
| 6 | MAWSLGVLVALPFPL | HLA-B*18:01 | P1_MD_pMHC_bulge | 0.449 | 0.850 | NA | pMHC_bulge_stability_MD_3x50ns |
| 7 | MAWSLGVLVALPFPL | HLA-A*25:01 | P1_MD_pMHC_bulge | 0.429 | 0.842 | NA | pMHC_bulge_stability_MD_3x50ns |
| 8 | MSSLAATTFHWKKCR | HLA-B*07:02 | P1_MD_pMHC_bulge | 0.428 | 0.835 | NA | pMHC_bulge_stability_MD_3x50ns |
| 9 | MSSLAATTFHWKKCR | HLA-A*68:01 | P1_MD_pMHC_bulge | 0.428 | 0.835 | NA | pMHC_bulge_stability_MD_3x50ns |
| 10 | MSSLAATTFHWKKCR | HLA-B*35:03 | P1_MD_pMHC_bulge | 0.425 | 0.834 | NA | pMHC_bulge_stability_MD_3x50ns |

## Recommended Use

1. Start with P0 exact paired TCR candidates: HMTEVVRHC/HLA-A*02:01 and GADGVGKSAL/HLA-C*08:02.
2. Recover full TCR alpha/beta and MHC sequences before TCR-pMHC MD.
3. For long peptides without exact TCR evidence, run pMHC stability MD first and keep the recognition claim off.
4. Only run mutant-WT counterfactual MD when WT peptide is available.

## Outputs

- `md_escalation_queue.tsv`
- `md_escalation_queue_top20.tsv`
- `md_simulation_tiers.tsv`
