# P3/P9 Full Public Execution Manifest

- Full execution started: 2026-05-06T02:25:36+09:00
- Branch: `exec/p3-p9-full-public-2026-05-06`
- Workspace: `/home/seungho/personal/THCA_data_analysis`
- Scope: Paper 3 ICI readiness atlas and Paper 9 synthetic-lethality vulnerability map.
- Paper 1 manuscript edits: prohibited; no manuscript prose, Hook, Aim, Discussion, Limitations, Cover, or Reviewer Q9 edits are in scope.
- Paper 2/4 HLA work: prohibited; no Paper 2/4 HLA execution is in scope.
- Protected data policy: protected BAM/WES/dbGaP/raw FASTQ or controlled-access files will not be downloaded. HLA LOH and neoantigen execution require already-local WES/BAM/VCF/HLA inputs; otherwise they are blocked.
- Data policy: public processed releases and already-local processed files may be used. Large raw BAM/WES/FASTQ/WSI/GPU workflows are excluded.
- Expected outputs: local inventory, shared lineage-state master table, Paper 3 signature scores/readiness/ecotypes/figures/report, Paper 9 dependency/drug-response vulnerability rankings/figures/report, and boundary summary.
- Disk budget: stop and ask before any single download over 10 GB; stop and ask if total new data would exceed 50 GB.
- Wall-time checkpoints: checkpoint and report if runtime exceeds 24 hours.
- Compute policy: CPU only unless GPU is explicitly approved later.
- Claim guards: no ICI response predictor claim, no clinical treatment recommendation, no patient treatment selection claim, no TROP2 synthetic-lethality claim.
