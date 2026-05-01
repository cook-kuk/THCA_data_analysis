# Dark Matter of Thyroid Cancer

> An 8-gene transcriptional axis sub-stratifies BRAF/RAS-negative thyroid cancer
> at single-cell resolution and identifies DICER1/EIF1AX as the FVPTC-like genomic anchor.

**Author**: Seungho Cook (corresponding); collaborators TBD.
**Status**: Phase 1 + Phase 2 sprint completed 2026-04-29 (single-day intensive).
**Target venue**: Cell Reports Medicine (1차) / Nature Communications (reach).

## TL;DR

- Driver-negative thyroid cancer (BRAF V600E−, RAS hotspot−) = **28% of TCGA-THCA**, currently un-stratified by molecular tests.
- Our 8-gene panel sub-stratifies this "Dark Matter" into:
  - **DM1** (n=89, mean age 42, cPTC-architectured, mechanism unknown — true mystery)
  - **DM2** (n=55, mean age 54, FVPTC-like, **DICER1/EIF1AX/PPM1D 9× enriched**, p=0.0175)
- Reproduced at single-cell resolution across 4 sc cohorts (r ≥ 0.86).
- Korean Yoo 2016 K2 cohort (n=180): **NBNR ↔ DM concordance 93.5%**, DICER1+EIF1AX 10.3% ≈ TCGA 10.9%.
- Pan-cancer specificity confirmed (LUAD r=0.235 vs THCA r=0.99).
- 8-gene → cPTC/FVPTC predictive AUC = 0.71 (5-fold CV).
- Bootstrap cluster ARI = 0.994 (random null = 0).

## Repository structure

```
project/
├── data_processed/        # bulk RNA-seq + microarray expression matrices
├── data_raw/              # GDC + GEO + ENA raw downloads
├── results/
│   ├── dark_matter_phase1/  # Phase 1 (Discovery)
│   └── dark_matter_phase2/  # Phase 2 + 3 + 4 (Validation + Deepening + ULTRA)
│       └── web/             # Live HTML dashboard
├── notebooks_or_scripts/  # analysis scripts
└── reports/html/         # public dashboard at http://40.82.129.113:8012/
```

## Live dashboard
http://40.82.129.113:8012/reports/html/dark_matter/v2.html

## Key scripts (reproducible analyses)
```
project/results/dark_matter_phase1/
├── run_steps_1_to_6.py      # Phase 1 viability test
├── step2_redo_pfi.py        # Liu 2018 TCGA-CDR Cox
├── sc_analysis.py           # GSE241184 sc pipeline

project/results/dark_matter_phase2/
├── p2a_gse193581.py         # P2-A1 external sc
├── p2a2_gse184362.py        # P2-A2 adult multi-patient
├── p2a2_multisite.py        # T/P/LN trajectory
├── p2b_yoo2016_parse.py     # Korean K2 mutation parse
├── p2c_dm1_deep_dive.py     # DM1 mechanism mining
├── p2d_driver_map.py        # 7-class driver landscape
├── p3_deepening.py          # DEG + pathway + AUC + sc + K2 fix
├── p4_final_batch.py        # Pan-cancer + ATC + Stage I + DEG table + BibTeX + Methods
├── p5_speed.py              # Age-stratified + DICER1+ PFI
└── p6_more.py               # Mutual exclusivity + cell composition + master table
```

## Software requirements
- Python 3.12, scanpy 1.12.1, lifelines 0.27, pandas 2.3, scipy 1.13, scikit-learn 1.5
- See `requirements.txt`

## Citation
If you use this work, cite:
> Cook S et al. (2026) *Dark matter of thyroid cancer: an 8-gene transcriptional axis identifies DICER1/EIF1AX-anchored FVPTC-like sub-cluster within driver-negative tumors* (manuscript in preparation).

BibTeX entry: see `data/references.bib`.

## License
MIT (see `LICENSE`). Data sources retain their original licenses.

## Acknowledgments
Yu Hyeong-Won (Seoul National University Bundang Hospital) — clinical guidance.
Yoo SK et al. (PLoS Genet 2016) — K2 Korean cohort dataset.
Original TCGA, GSE241184/GSE193581/GSE184362 publication teams.
