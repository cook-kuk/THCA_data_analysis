# RAI Response Genomics Atlas — thyroid cancer

Reproducible workspace for identifying and re-analyzing public datasets that contain **radioactive iodine (RAI) response, RAI avidity, RAI refractoriness, remission/persistence, or redifferentiation** labels in thyroid cancer, and for validating an **eight-gene thyroid differentiation / iodide-handling panel** against those labels.

> **Manuscript framing**: "An eight-gene thyroid differentiation / iodide-handling silence axis identifies patients likely to fail radioactive iodine therapy, including a molecular gray zone not fully explained by BRAF/RAS/TERT or standard histologic categories."

This workspace is sibling-level inside the existing `THCA_data_analysis` repository. It does **not** modify Paper 1 or any other ongoing track.

## Quick layout

```
rai-response-genomics-atlas/
├── README.md
├── config/
│   ├── datasets.yaml              # priority datasets + accession + access conditions
│   └── eight_gene_panel.yaml      # SLC5A5/TPO/TG/TSHR/PAX8/NKX2-1/FOXE1/DIO1 (Paper 1 confirmed)
├── docs/
│   ├── literature_review.md
│   ├── dataset_inventory.md
│   ├── rai_response_label_taxonomy.md
│   ├── manuscript_strategy.md
│   ├── reviewer_risk_register.md
│   ├── data_request_boucai.md
│   └── data_request_mu_hra004166.md
├── scripts/
│   ├── 00_setup_environment.py
│   ├── 01_search_geo_metadata.py
│   ├── 02_download_geo_dataset.py
│   ├── 03_parse_geo_metadata.py
│   ├── 04_score_gene_panel.py
│   ├── 05_validate_rai_labels.py
│   ├── 06_make_figures.py
│   └── 07_build_dataset_inventory.py
├── notebooks/
│   ├── 01_gse151179_validation.ipynb
│   ├── 02_gse299988_validation.ipynb
│   ├── 03_tcga_thca_discovery_plan.ipynb
│   └── 04_redifferentiation_dataset_scan.ipynb
├── data → /data2/rai_atlas        # symlink, large outputs land on Premium SSD /data2
│   ├── raw/        # GEO matrices, supplementary tables as fetched
│   ├── interim/    # parsed metadata, gene-mapped expression
│   ├── processed/  # panel scores, label-joined tables
│   └── external/   # supplementary tables from papers (PDF/XLSX → TSV)
├── results/
│   ├── tables/
│   ├── figures/
│   └── reports/
└── logs/
```

## Eight-gene panel (canonical, Paper 1 confirmed)

`SLC5A5, TPO, TG, TSHR, PAX8, NKX2-1, FOXE1, DIO1` — see `config/eight_gene_panel.yaml`. **DIO1**, not DUOX2 (the placeholder in some specs). Panel score = mean z-score of the 8 genes within the cohort.

## How to run

```bash
cd /home/seungho/personal/THCA_data_analysis/rai-response-genomics-atlas
python3 scripts/00_setup_environment.py            # check libs, write logs/env.txt
python3 scripts/01_search_geo_metadata.py          # list priority accessions
python3 scripts/02_download_geo_dataset.py GSE151179
python3 scripts/03_parse_geo_metadata.py GSE151179
python3 scripts/04_score_gene_panel.py GSE151179
python3 scripts/05_validate_rai_labels.py GSE151179
python3 scripts/06_make_figures.py GSE151179
python3 scripts/07_build_dataset_inventory.py
```

## Priority datasets

| Tier | Dataset | Accession | Why |
|------|---------|-----------|-----|
| 1 (direct response) | Boucai 2023 CCR exceptional responders | request | RECIST RAI response; expression + genomics |
| 1 (direct response) | Mu 2024 JCEM RAI uptake patterns | NGDC HRA004166 | 220 patients, 4 RAI-uptake patterns, targeted NGS |
| 2 (avidity) | GSE151179 / GSE151180 / GSE151181 | GEO | RAI-avid vs RAI-refractory PTC, mRNA + miRNA |
| 2 (avidity) | GSE299988 | GEO | small PTC, RAI-avid vs non-avid |
| 4 (proxy) | TCGA-THCA | GDC | discovery cohort, no direct RAI labels |
| 5 (redifferentiation) | GSE112202 + MERAIODE trial supplements | GEO + trial | redifferentiation after MAPK inhibition |
| 4 (mechanism) | GSE184362 | GEO | scRNA, malignant cell states |

See `config/datasets.yaml` and `docs/dataset_inventory.md` for details.

## Scientific style rules (enforced in scripts and docs)

- Never say "predicts RAI response" unless a dataset has true response labels.
- Use: *associated with, stratifies, consistent with RAI refractoriness, captures a differentiation-linked RAI failure axis, supports clinical relevance*.
- Always document: exact sample count · exact label source · number of usable genes · missing genes · platform · normalization status.
- Each dataset is mapped to a **Tier** in the label taxonomy (`docs/rai_response_label_taxonomy.md`).

## Storage convention

`data/` is a symlink to `/data2/rai_atlas/` (Premium SSD, 438 G free). All large GEO matrices land there; the code path stays portable. `results/` and `logs/` are kept on root disk (small).
