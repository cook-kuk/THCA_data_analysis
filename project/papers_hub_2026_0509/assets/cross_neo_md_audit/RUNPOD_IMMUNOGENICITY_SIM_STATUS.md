# RunPod Immunogenicity Simulation Status

Checked: 2026-05-10 12:22 KST

## Active queue

- Pod: `thca-neo-bayesian-aux`
- SSH endpoint used: `135.84.176.142:20878`
- Queue watcher PID: `3290`
- Queue script: `queue_after_pid.sh 2213 run_immediate_ready_batch.sh`
- Status: the watcher detected PID `2213` had finished and started `run_immediate_ready_batch.sh`.

## Currently running job

- PID: `3757`
- Job: `tier1_full_10ns_replicate_or_control_GADGVGKSAL_HLA-C_08_02_mutant_TCR-pMHC_6UON_GADGVGKSAL_HLA-C0802_chainF_rep01_10p0ns`
- Module meaning: additional mutant TCR-pMHC explicit-solvent replicate for the lead `GADGVGKSAL / HLA-C*08:02` wetlab candidate.
- Latest synced remote state at check:
  - time: `3899.999999921601 ps`
  - temperature: `300.7206488292319 K`
  - speed: `44.2 ns/day`

## Newly completed and synced

- `HMTEVVRHC / HLA-A*02:01` primary 10 ns run synced from RunPod.
- Latest parsed state: `10000.000002 ps`, `299.891412 K`.
- QC: `ok`, peptide RMSD final `0.065110 nm`.
- MD evidence: `MD_VERY_STRONG`, score `0.821789`.
- Top-13 decision package updated from `PENDING_MD_COMPLETION` to `TIER_A_EXPERIMENT_NOW_WITH_WT_CURATION`.

## Boundary

This is simulation support for immunogenicity testing. It is not immunogenicity proof, clinical validation, or a replacement for IFN-gamma / ICS / killing assays.
