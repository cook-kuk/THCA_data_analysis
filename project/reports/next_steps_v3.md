# THCA v3 — prioritised next steps

## P0 — immediate (data already on disk)

- [ ] Re-run panel-k curve with QUBO-neal on a dense k grid (every 2) to get a smoother plateau estimate.
- [ ] Bootstrap 95% CIs over external AUCs (currently only TCGA CV has CIs).
- [ ] Re-check sample_master labels for GSE27155 `unknown` samples — ~29/99 are unlabeled; partial-label semi-supervised may help.

## P1 — pulls that require extra raw data

- [ ] Fetch TCGA-THCA Illumina 450k beta values from GDC → enables inter-cohort methylation LODO with GSE97466.
- [ ] Fetch STAR-Fusion / Arriba callset for TCGA + GSE213647 → Step 7 Task A becomes real.
- [ ] Fetch GISTIC2 focal / arm-level SCNA for TCGA → Step 7 Task B.

## P2 — methodological upgrades

- [ ] Replace isotonic recalibration with an explicit domain-adaptation step (e.g. CORAL or subsampled DANN) and compare.
- [ ] Explore a compact 'MAPK-free' panel derived by QUBO on the leakage-clean feature pool, with cross-platform constraints.
- [ ] Stratify survival analysis by ATA risk tier instead of TDS tertile, once sample_master carries complete ATA labels.

## P3 — product / governance

- [ ] Publish the honesty audit as a standalone supplementary note.
- [ ] Add a per-prediction uncertainty field to the decision-support prototype (Platt + isotonic disagreement as a simple proxy).
- [ ] Draft a data-sheet-for-datasets style card for each cohort (coverage, platform, labelling provenance, known biases).
