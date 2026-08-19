# CROSS-Neo-TCR Reviewer Q&A Fact Blocks

Status: factual defense scaffold for the TCR-aware extension. This is not final manuscript voice.

## Q1. Does CROSS-Neo-TCR replace the main pMHC neoantigen model?

No. The TCR branch is optional and only used when TCR evidence is available. The main CROSS-Neo ranker remains the pMHC/presentation-centered model because most neoantigen immunogenicity datasets lack paired TCR alpha/beta sequences.

## Q2. How much paired TCR data exists locally?

The local TCR registry contains 537,618 rows. Paired alpha/beta TCR rows total 64,037. Any TCR sequence is available for 123,035 rows.

## Q3. How many rows contain peptide-HLA-TCR labels?

There are 98,718 peptide-HLA-TCR labeled rows with any TCR sequence and 56,329 with paired alpha/beta.

## Q4. How much of CROSS-Neo overlaps TCR resources?

Among 2,715 CROSS-Neo rows, 2,016 have some TCR-resource overlap, 101 are exact paired TCR-pMHC matches, and 1,876 have cancer-context evidence. There are 699 rows with no TCR match.

## Q5. Are peptide-only TCR matches treated as ground truth?

No. Peptide-only and near-peptide TCR matches are diagnostic evidence only. Labels are not transferred from pathogen epitopes to cancer neoantigens without compatible annotation.

## Q6. Which external TCR models currently work locally?

Runnable local models include pMTnet, TEPCAM, ERGO-II, PanPep, NetTCR-2.2, TSpred, and TCR-H assets. pMTnet and TEPCAM currently have the strongest practical fit to available registry fields.

## Q7. What is the best current TCR expert score?

The no-training mean pMTnet/TEPCAM ensemble is the best practical diagnostic score in the fast challenge benchmark: AUPRC 0.7841 and AUROC 0.7652 on 1,191 rows across six panels.

## Q8. Why not promote ERGO-II or PanPep now?

ERGO-II runtime works, but the beta-only/unknown-VJ adapter is weak in the fast pilot. ERGO-II should be restricted to richer rows with alpha/beta, V/J, peptide, and MHC fields. PanPep zero-shot was near random in the current pilot and should remain diagnostic or negative-control unless few-shot mode improves.

## Q9. Does the source-heldout benchmark support a broad claim?

No. Source-heldout performance is source-dependent: VDJdb heldout AUPRC 0.6523 / AUROC 0.6040, McPAS-TCR heldout AUPRC 0.8025 / AUROC 0.8232, and VDJdb 10x heldout AUPRC 0.8408 / AUROC 0.8340. This supports diagnostic utility, not universal TCR-aware prediction.

## Q10. Does structure improve over sequence-only?

Not yet. The structure pipeline is ready as QC/parser infrastructure, but most rows lack full TCR and MHC chains. Structure-derived improvement cannot be claimed until modelable mutant and wildtype complexes are generated and evaluated.

## Q11. Which wetlab candidates are most supported by external TCR experts?

The top external-supported candidate is HMTEVVRHC / HLA-A*02:01, with four exact modelable TCR rows, external mean score 0.7121, and max 0.8029. The second is GADGVGKSAL / HLA-C*08:02, with five exact modelable TCR rows, external mean 0.5703, and max 0.6790.

## Q12. What is the safest current manuscript positioning?

Promote CROSS-Neo-TCR as a diagnostic case-study and wetlab prioritization branch. Hold the main TCR-aware method claim until strict external paired-TCR benchmarks, stronger structure validation, and ideally wetlab validation are complete.

## Q13. Can prediction-fragile neoantigens be escalated to MD?

Yes, but as a targeted diagnostic escalation, not a broad default. The current queue promotes HMTEVVRHC / HLA-A*02:01 and GADGVGKSAL / HLA-C*08:02 to P0 TCR-pMHC MD after full TCR/MHC chain recovery. Long or conformationally uncertain class-I peptides are P1 pMHC-stability MD candidates. The prepared top-20 MD stubs total 3,300 ns planned production time if all are eventually run. MD evidence should support interface plausibility, peptide stability, and mutant-WT cross-reactivity hypotheses, not prove immunogenicity.

## Q14. Did any P0 structure pilot actually pass local MD preparation?

Yes. Existing solved TCR-pMHC PDBs were found for both P0 candidates. HMTEVVRHC / HLA-A*02:01 has ready pilot complexes from 6VRN, 6VRM, 6VQO, and 7RM4; GADGVGKSAL / HLA-C*08:02 has ready pilot complexes from 6UON. Representative complexes 6VRN and 6UON passed PDBFixer repair plus OpenMM dry minimization. This is still a preparation/QC result, not production MD and not wetlab validation.

## Q15. Did the OpenMM runner create trajectories?

Yes. Local CPU vacuum smoke tests succeeded for both P0 representative complexes and produced DCD trajectories, state logs, final PDBs, and MDTraj contact summaries. Final peptide RMSD was 0.0091 nm for 6VRN/HMTEVVRHC and 0.0144 nm for 6UON/GADGVGKSAL. These smoke tests validate the runner and analysis plumbing only; the real next step is CUDA explicit-solvent 10 ns pilot MD.
