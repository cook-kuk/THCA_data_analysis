# CROSS-Neo-TCR MD Job Stubs

These are execution stubs for candidates already selected by the MD escalation queue.

- Planned jobs: 20
- P0 TCR-pMHC jobs: 2
- Planned aggregate production time if all top-20 jobs run: 3300 ns
- GROMACS available now: False
- OpenMM available now: True

Do not launch these blindly. First supply QC-passed PDB structures, repair missing atoms/protonation, choose force field and water model, and confirm peptide/HLA/TCR chain identities.

Outputs:

- `md_job_manifest.tsv`
- `md_engine_availability.tsv`
- `submit_md_jobs.template.sh`
