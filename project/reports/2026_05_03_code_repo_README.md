# Code Repository Prep — github.com/seunghocook/thyca-paper-2026

**Date:** 2026-05-03 (marathon scaffolding/infra)
**Purpose:** GitHub repo structure + README + requirements for paper submission.

---

## Proposed repo structure

```
thyca-paper-2026/
├── README.md                          (this file as template)
├── requirements.txt                   (pinned dependencies)
├── LICENSE                            (MIT or CC-BY-4.0)
├── CITATION.cff                       (citation metadata)
├── .gitignore                         (data/ excluded, only metadata committed)
├── notebooks_or_scripts/              (renamed from current)
│   ├── pillar1_HLA/
│   │   ├── arcasHLA_GSE286332.sh       (burst pipeline)
│   │   ├── parse_genotype_json.py
│   │   └── forest_meta.py              (= v17_paper2_pillar1_forest.py)
│   ├── pillar2_GSE286332/
│   │   └── ptc_vs_ptcht_DEG_GSEA.py    (= v17_P3_GSE286332_ptc_vs_ptcht.py)
│   ├── pillar3_driver/
│   │   └── driver_mrna_audit.py        (= v17_P1_driver_mrna_audit.py)
│   ├── pillar4_robustness/
│   │   └── pangenome_vs_tiera67.py     (= v17_P4_pangenome_vs_tiera67.py)
│   ├── pillar5_autoimmune/
│   │   ├── pdm1_mediation.py           (= v17_D3P5_pdm1_gradient.py)
│   │   ├── tcga_hashimoto_signature.py (= v17_D4P2_tcga_hashimoto_signature.py)
│   │   ├── bcr_repertoire.py           (= v17_D5P6_bcr_repertoire.py)
│   │   └── dm1_subcluster.py           (= v17_D6P7_dm1_subcluster.py)
│   ├── replication/
│   │   ├── korean_GSE213647_hashimoto.py (= v17_D8B_korean_replication.py)
│   │   └── subB_x_K2_NBNR.py             (= v17_D8C_dm1_subB_x_K2_NBNR.py)
│   ├── figures/
│   │   ├── F1_panasian_HLA_forest.py
│   │   ├── F2_gse286332_multipanel.py
│   │   ├── F3_driver_neutrality.py
│   │   ├── F4_pangenome_robustness.py
│   │   └── F5_autoimmune_PTC_mechanism.py
│   ├── suppl_tables/
│   │   └── build_S1_S8.py              (= v17_suppl_tables_build.py)
│   └── helpers/
│       ├── ensembl_to_symbol.tsv         (62,700 mappings)
│       └── tierA67_genes.txt
├── submission_data/                    (TSVs ready for submission)
│   ├── S1_cohort_assembly.tsv
│   ├── S2_TIERA67_candidate_pool.tsv
│   ├── S3_gse286332_top300_DEGs.tsv
│   ├── S4_panasian_HLA_forest_full.tsv
│   ├── S5_mediation_analysis.tsv
│   ├── S5b_OLS_decomposition.tsv
│   ├── S6_tcga_hashimoto_DM_crosstab.tsv
│   ├── S7_cross_cohort_generalization.tsv
│   └── S8a/b/c_dm1_subcluster_*.tsv
├── figures/                            (final main + suppl PDFs)
│   └── (F1-F5 PDF + PNG)
├── docs/
│   ├── methods_M1_M11.md               (= 2026_05_03_methods_M1_M11_scaffold.md)
│   ├── reviewer_QA.md
│   ├── ATA_2015_cheatsheet.md
│   └── STAR_Methods_KeyResources.md
└── tests/                              (smoke tests for reproducibility)
    └── test_signature_score.py
```

---

## README.md template

```markdown
# THCA Paper 2026 — Phase 0 Cancer Paper Code & Data

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.TBD.svg)](https://doi.org/10.5281/zenodo.TBD)
[![Code License](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

This repository accompanies the manuscript **"[paper title — voice-protected]"** by Cook et al. (2026), establishing a transcriptional differentiation axis (DM1/DM2) and autoimmune-PTC mechanism layer for papillary thyroid carcinoma stratification.

## Overview

Five-pillar evidence structure:
1. **Korean Pan-Asian HLA cohort (n=874)** — arcasHLA RNA-seq imputation + Chu 2018 forest meta
2. **GSE286332 PTC vs PTC+HT molecular dissection** — 10,380 DEGs, IFN-γ + HLA-II
3. **Driver mRNA neutrality** — BRAF transcript Cohen d=−0.04 vs WT
4. **Pan-genome cluster robustness** — TIERA67 ARI=0.90 ≈ unbiased top-5000
5. **Autoimmune-PTC mechanism layer** — TCGA Hashimoto-DM2 OR=0.20 (p=6e-10), HLA-II 140% mediation, BCR clonal + TLS d=+1.96, sub-B NBNR cluster

## Quick start

```bash
git clone https://github.com/seunghocook/thyca-paper-2026.git
cd thyca-paper-2026
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
```

## Reproduce main figures

```bash
# Pillar 1
python notebooks_or_scripts/pillar1_HLA/forest_meta.py

# Pillar 2
python notebooks_or_scripts/pillar2_GSE286332/ptc_vs_ptcht_DEG_GSEA.py

# Pillar 3
python notebooks_or_scripts/pillar3_driver/driver_mrna_audit.py

# Pillar 4
python notebooks_or_scripts/pillar4_robustness/pangenome_vs_tiera67.py

# Pillar 5
python notebooks_or_scripts/pillar5_autoimmune/pdm1_mediation.py
python notebooks_or_scripts/pillar5_autoimmune/tcga_hashimoto_signature.py
python notebooks_or_scripts/pillar5_autoimmune/bcr_repertoire.py
python notebooks_or_scripts/pillar5_autoimmune/dm1_subcluster.py

# Replication
python notebooks_or_scripts/replication/korean_GSE213647_hashimoto.py
python notebooks_or_scripts/replication/subB_x_K2_NBNR.py

# Figures
for f in notebooks_or_scripts/figures/F*.py; do python $f; done

# Suppl tables
python notebooks_or_scripts/suppl_tables/build_S1_S8.py
```

## Data

- **Raw RNA-seq:** TCGA-THCA (GDC), GSE213647 + GSE286332 (NCBI GEO), PRJEB11591 (ENA)
- **HLA:** arcasHLA imputation outputs at `submission_data/HLA/`
- **Chu 2018 reference:** Published table values from PMC6161647

## Citation

```bibtex
@article{Cook2026thca,
  author = {Cook, Seungho and ...},
  title = {[paper title]},
  journal = {[journal]},
  year = {2026}
}
```

## License

Code: MIT License. Documentation: CC-BY-4.0.

## Contact

Lead: Seungho Cook (kukshomr@gmail.com)

```

---

## requirements.txt

```
# Core scientific stack
numpy>=1.26
pandas>=2.0
scipy>=1.11
statsmodels>=0.14
scikit-learn>=1.3
matplotlib>=3.7

# Differential expression + GSEA
pydeseq2==0.5.4
gseapy==1.1.13

# Data fetching
biopython>=1.81
requests>=2.31
pyyaml>=6.0

# (External tools — install separately)
# arcasHLA v0.6.0
# kallisto v0.46.1
# bowtie2 v2.4.4
# samtools
# sra-tools 3.2.1
# IPD-IMGT/HLA v3.44.0
```

---

## .gitignore

```
# Raw data (never committed)
data/
/data/thca/
*.fastq.gz
*.bam
*.sra

# Compiled python
__pycache__/
*.py[cod]
*.so

# Distribution
build/
dist/
*.egg-info/

# Virtual env
.venv/
venv/

# IDE
.vscode/
.idea/

# Logs
*.log
/tmp/

# OS
.DS_Store
Thumbs.db

# Submission cache (already in submission_data/, don't double-commit)
results/figures/*/*.png
```

---

## LICENSE (MIT)

```
MIT License

Copyright (c) 2026 Seungho Cook

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

---

## CITATION.cff

```yaml
cff-version: 1.2.0
title: "THCA Paper 2026 — Phase 0 Cancer Paper Code & Data"
authors:
  - family-names: "Cook"
    given-names: "Seungho"
    orcid: "https://orcid.org/0000-0000-0000-0000"
type: software
version: "1.0"
doi: "10.5281/zenodo.TBD"
date-released: 2026-06-13
license: MIT
repository-code: "https://github.com/seunghocook/thyca-paper-2026"
keywords:
  - thyroid cancer
  - PTC
  - HLA
  - autoimmune
  - Hashimoto's
  - transcriptional differentiation
```

---

## Pre-submission checklist

- [ ] GitHub repo public
- [ ] Zenodo doi reserved (release v1.0)
- [ ] All scripts re-tested with `pip install -r requirements.txt`
- [ ] `data/` paths abstracted (config file or environment variable)
- [ ] README.md final (with paper title from voice-protected section)
- [ ] LICENSE + CITATION.cff committed
- [ ] Smoke tests pass: `pytest tests/`
- [ ] All figures regeneratable from script + submission_data/

---

## Voice-protected items (NOT in this scaffold)

- Paper title (Decision 1 — user voice)
- Abstract first paragraph (user voice)
- README.md "Overview" deeper framing tone (user voice for hook + claim language)
- Author list ordering & ORCID (user + Yu professor decision)
