# One-Page Paper Summary: CROSS-Neo MD Audit Layer

We added an explicit-solvent OpenMM MD audit layer for selected CROSS-Neo/TCR-aware neoantigen candidates. The layer is designed to test structural plausibility and diagnostic consistency, not to replace the main pMHC ranking model and not to prove immunogenicity.

`GADGVGKSAL / HLA-C*08:02` completed a 10 ns CUDA explicit-solvent run and one 1 ns alternate screen. The primary run reached an MD evidence label of `MD_MODERATE`, with final peptide RMSD 0.139 nm and persistent TCR-peptide contact signal. This supports wetlab prioritization and case interpretation, while remaining below the threshold for any immunogenicity or clinical claim.

`HMTEVVRHC / HLA-A*02:01` remains a high-priority active run. The latest synced state reached 10000.0 ps at 299.89 K, but full DCD-based RMSD/contact analysis is pending. Claims for this candidate should wait until the 10 ns trajectory and follow-up replicate screens are parsed.

Integration with CROSS-Neo predictions found 182 MD-supported model-high rows and 48 model-low/MD-high disagreement rows. These cases are useful for model audit and wetlab triage. No WT or decoy trajectories are available yet, so mutant-specific recognition and WT cross-reactivity claims remain blocked.
