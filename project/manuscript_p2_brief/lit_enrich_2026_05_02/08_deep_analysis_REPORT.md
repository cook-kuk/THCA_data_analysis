# 08 lit_enrich DEEP analysis report

_Generated: 2026-05-02 15:46_

Builds on raw `data/*.json` (5.9 MB) — does dedupe + cross-source merge + ranking,
outputs paper-ready tables. Does NOT replace `01-07_*.md` (those stay as the
narrative summary); this is the **structured drill-down**.

---

## 1. Reference verification — triangulation matrix

_All 17 entries × 4 sources (CrossRef / PubMed / OpenAlex / EuropePMC). Sorted by_
_(incomplete-first, then n_confirm-desc). **8 incomplete entries**_
_remain; **9/17 have ≥1 confirming source**._

| key                | incomplete   |   n_confirm | crossref   | pubmed   | openalex   | europepmc   | crossref_title                                                                   |
|:-------------------|:-------------|------------:|:-----------|:---------|:-----------|:------------|:---------------------------------------------------------------------------------|
| Ringel2025         | True         |           0 | False      | False    | False      | False       |                                                                                  |
| Liu2017            | True         |           0 | False      | False    | False      | False       |                                                                                  |
| Wang2024           | True         |           0 | False      | False    | False      | False       |                                                                                  |
| SEER_PTC           | True         |           0 | False      | False    | False      | False       |                                                                                  |
| Lu2023             | True         |           0 | False      | False    | False      | False       |                                                                                  |
| Bradley2010        | True         |           0 | False      | False    | False      | False       |                                                                                  |
| Krishnamoorthy2025 | True         |           0 | False      | False    | False      | False       |                                                                                  |
| Chu2018            | True         |           0 | False      | False    | False      | False       |                                                                                  |
| Haugen2016         | False        |           1 | True       | False    | False      | False       | 2015 American Thyroid Association Management Guidelines for Adult Patients with  |
| TCGA2014           | False        |           1 | True       | False    | False      | False       | Integrated Genomic Characterization of Papillary Thyroid Carcinoma               |
| Yoo2016            | False        |           1 | True       | False    | False      | False       | Comprehensive Analysis of the Transcriptional and Mutational Landscape of Follic |
| Landa2016          | False        |           1 | True       | False    | False      | False       | Genomic and transcriptomic hallmarks of poorly differentiated and anaplastic thy |
| Pu2021             | False        |           1 | True       | False    | False      | False       | Single-cell transcriptomic analysis of the tumor ecosystems underlying initiatio |
| Wirth2020          | False        |           1 | True       | False    | False      | False       | Efficacy of Selpercatinib in             <i>RET</i>             -Altered Thyroid |
| Xing2014           | False        |           1 | True       | False    | False      | False       | Molecular pathogenesis and mechanisms of thyroid cancer                          |
| Cibas2017          | False        |           1 | True       | False    | False      | False       | The 2017 Bethesda System for Reporting Thyroid Cytopathology                     |
| Pozdeyev2018       | False        |           1 | True       | False    | False      | False       | Genetic Analysis of 779 Advanced Differentiated and Anaplastic Thyroid Cancers   |

---

## 2. Competitive landscape — Top-5 per claim (cross-source dedupe)

_Mined 8 claims × 5 sources (OpenAlex/SemScholar/E-PMC/bioRxiv/arXiv);_
_dedupe collapsed 279 raw rows → 40 top-5 papers across 8 claims._
_Composite rank = `n_sources*100 + min(cited_by,500) + influential*5 + (year-2000)`._

### 8-gene / driver-excluded PTC stratification

|   n_sources |   year | first_author   | title                                                                            |   cited_by |   influential | doi                          |
|------------:|-------:|:---------------|:---------------------------------------------------------------------------------|-----------:|--------------:|:-----------------------------|
|           1 |   2022 | Anand          | Cancer chemotherapy and beyond: Current status, drug candidates, associated risk |       1578 |             0 | 10.1016/j.gendis.2022.02.007 |
|           1 |   2021 | McGrail        | High tumor mutation burden fails to predict immune checkpoint blockade response  |       1039 |             0 | 10.1016/j.annonc.2021.02.006 |
|           1 |   2020 | Jiang          | Role of PI3K/AKT pathway in cancer: the framework of malignant behavior          |        570 |             0 | 10.1007/s11033-020-05435-1   |
|           1 |   2020 | Tran           | Advances in bladder cancer biology and therapy                                   |        779 |             0 | 10.1038/s41568-020-00313-1   |
|           1 |   2019 | Yuan           | Mechanisms underlying the activation of TERT transcription and telomerase activi |        460 |             0 | 10.1038/s41388-019-0872-9    |

### BCR repertoire / TLS in thyroid cancer

|   n_sources |   year | first_author   | title                                                                   |   cited_by |   influential | doi                        |
|------------:|-------:|:---------------|:------------------------------------------------------------------------|-----------:|--------------:|:---------------------------|
|           1 |   2023 | Sun            | T cells in health and disease                                           |        874 |             0 | 10.1038/s41392-023-01471-y |
|           1 |   2023 | Pisetsky       | Pathogenesis of autoimmune disease                                      |        573 |             0 | 10.1038/s41581-023-00720-1 |
|           1 |   2022 | Wang           | Therapeutic peptides: current applications and future directions        |       1890 |             0 | 10.1038/s41392-022-00904-4 |
|           1 |   2021 | Zhu            | Combination strategies to maximize the benefits of cancer immunotherapy |        645 |             0 | 10.1186/s13045-021-01164-5 |
|           1 |   2020 | Jin            | The updated landscape of tumor microenvironment and drug repurposing    |       1218 |             0 | 10.1038/s41392-020-00280-x |

### BRAF/RAS-negative ('dark matter') PTC subtypes

|   n_sources |   year | first_author   | title                                                                            |   cited_by |   influential | doi                        |
|------------:|-------:|:---------------|:---------------------------------------------------------------------------------|-----------:|--------------:|:---------------------------|
|           1 |   2023 | Bahar          | Targeting the RAS/RAF/MAPK pathway for cancer therapy: from mechanism to clinica |        723 |             0 | 10.1038/s41392-023-01705-z |
|           1 |   2023 | Glaviano       | PI3K/AKT/mTOR signaling transduction pathway and targeted therapies in cancer    |       1787 |             0 | 10.1186/s12943-023-01827-6 |
|           1 |   2020 | Jiang          | Role of PI3K/AKT pathway in cancer: the framework of malignant behavior          |        570 |             0 | 10.1007/s11033-020-05435-1 |
|           1 |   2020 | Yoo            | Glutamine reliance in cell metabolism                                            |        961 |             0 | 10.1038/s12276-020-00504-8 |
|           1 |   2020 | Degirmenci     | Targeting Aberrant RAS/RAF/MEK/ERK Signaling for Cancer Therapy                  |        531 |             0 | 10.3390/cells9010198       |

### HLA-II PTC autoimmunity (DPB1*05:01 etc.)

|   n_sources |   year | first_author   | title                                                                            |   cited_by |   influential | doi                          |
|------------:|-------:|:---------------|:---------------------------------------------------------------------------------|-----------:|--------------:|:-----------------------------|
|           1 |   2025 | Bryliński      | Effects of Trace Elements on Endocrine Function and Pathogenesis of Thyroid Dise |         24 |             0 | 10.3390/nu17030398           |
|           1 |   2021 | Mikosch        | Hashimoto’s thyroiditis and coexisting disorders in correlation with HLA status— |         26 |             0 | 10.1007/s10354-021-00879-x   |
|           1 |   2023 | Xiong          | How does SARS‐CoV‐2 infection impact on immunity, procession and treatment of pa |         16 |             0 | 10.1002/jmv.28487            |
|           1 |   2019 | Lu             | The Major Histocompatibility Complex Class II–CD4 Immunologic Synapse in Alcohol |         17 |             0 | 10.1016/j.ajpath.2019.09.019 |
|           1 |   2023 | Ding           | Unveiling the mystery of Riehl's melanosis: An update from pathogenesis, diagnos |         13 |             0 | 10.1111/pcmr.13108           |

### Hashimoto-like signature & PTC outcomes

|   n_sources |   year | first_author       | title                                                                            |   cited_by |   influential | doi                        |
|------------:|-------:|:-------------------|:---------------------------------------------------------------------------------|-----------:|--------------:|:---------------------------|
|           1 |   2022 | Kłubo-Gwieździńska | Hashimoto thyroiditis: an evidence-based guide: etiology, diagnosis and treatmen |        181 |             0 | 10.20452/pamw.16222        |
|           1 |   2023 | Zeng               | Understanding tumour endothelial cell heterogeneity and function from single-cel |        172 |             0 | 10.1038/s41568-023-00591-5 |
|           2 |   2021 | Pan                | Papillary Thyroid Carcinoma Landscape and Its Immunological Link With Hashimoto  |         42 |             0 | 10.3389/fcell.2021.758339  |
|           1 |   2021 | Romei              | A Narrative Review of Genetic Alterations in Primary Thyroid Epithelial Cancer   |         89 |             0 | 10.3390/ijms22041726       |
|           1 |   2021 | Menicali           | Immune Landscape of Thyroid Cancers: New Insights                                |         87 |             0 | 10.3389/fendo.2020.637826  |

### Pan-Asian HLA fine-mapping (Graves' / autoimmune thyroid)

|   n_sources |   year | first_author      | title                                                                            |   cited_by |   influential | doi                        |
|------------:|-------:|:------------------|:---------------------------------------------------------------------------------|-----------:|--------------:|:---------------------------|
|           1 |   2018 | Rojas             | Molecular mimicry and autoimmunity                                               |        597 |             0 | 10.1016/j.jaut.2018.10.012 |
|           1 |   2017 | Matzaraki         | The MHC locus and genetic susceptibility to autoimmune and infectious diseases   |        597 |             0 | 10.1186/s13059-017-1207-1  |
|           1 |   2016 | Yazdani           | Selective IgA Deficiency: Epidemiology, Pathogenesis, Clinical Phenotype, Diagno |        197 |             0 | 10.1111/sji.12499          |
|           1 |   2018 | Tye–Din           | Celiac Disease: A Review of Current Concepts in Pathogenesis, Prevention, and No |        193 |             0 | 10.3389/fped.2018.00350    |
|           1 |   2023 | Vargas‐Uricoechea | Molecular Mechanisms in Autoimmune Thyroid Disease                               |        178 |             0 | 10.3390/cells12060918      |

### Single-cell PTC progression

|   n_sources |   year | first_author     | title                                                                            |   cited_by |   influential | doi                          |
|------------:|-------:|:-----------------|:---------------------------------------------------------------------------------|-----------:|--------------:|:-----------------------------|
|           1 |   2023 | Glaviano         | PI3K/AKT/mTOR signaling transduction pathway and targeted therapies in cancer    |       1787 |             0 | 10.1186/s12943-023-01827-6   |
|           1 |   2022 | Anand            | Cancer chemotherapy and beyond: Current status, drug candidates, associated risk |       1578 |             0 | 10.1016/j.gendis.2022.02.007 |
|           1 |   2022 | Li               | Lactate metabolism in human health and disease                                   |       1243 |             0 | 10.1038/s41392-022-01151-3   |
|           1 |   2022 | Tong             | Targeting cell death pathways for cancer therapy: recent developments in necropt |        803 |             0 | 10.1186/s13045-022-01392-3   |
|           1 |   2021 | Dhatchinamoorthy | Cancer Immune Evasion Through Loss of MHC Class I Antigen Presentation           |       1003 |             0 | 10.3389/fimmu.2021.636568    |

### TIERA-like Korean PTC molecular cohort

|   n_sources |   year | first_author   | title                                                                            |   cited_by |   influential | doi                        |
|------------:|-------:|:---------------|:---------------------------------------------------------------------------------|-----------:|--------------:|:---------------------------|
|           1 |   2020 | Aaltonen       | Pan-cancer analysis of whole genomes                                             |       3240 |             0 | 10.1038/s41586-020-1969-6  |
|           1 |   2020 | Yang           | Brief introduction of medical database and data mining technology in big data er |        592 |             0 | 10.1111/jebm.12373         |
|           1 |   2018 | Si             | The roles of metallothioneins in carcinogenesis                                  |        442 |             0 | 10.1186/s13045-018-0645-x  |
|           1 |   2018 | Patel          | Performance of a Genomic Sequencing Classifier for the Preoperative Diagnosis of |        375 |             0 | 10.1001/jamasurg.2018.1153 |
|           1 |   2019 | Abdullah       | Papillary Thyroid Cancer: Genetic Alterations and Molecular Biomarker Investigat |        354 |             0 | 10.7150/ijms.29935         |

---

## 3. Bibliography GAP candidates — top-3 per claim NOT in our bib

_These are high-rank papers from the dedupe that don't appear in our 17-entry bib._
_Triage candidates for Discussion §3 expansion._

### 8-gene / driver-excluded PTC stratification

|   n_sources |   year | first_author   | title                                                                            |   cited_by | doi                          |
|------------:|-------:|:---------------|:---------------------------------------------------------------------------------|-----------:|:-----------------------------|
|           1 |   2022 | Anand          | Cancer chemotherapy and beyond: Current status, drug candidates, associated risk |       1578 | 10.1016/j.gendis.2022.02.007 |
|           1 |   2021 | McGrail        | High tumor mutation burden fails to predict immune checkpoint blockade response  |       1039 | 10.1016/j.annonc.2021.02.006 |
|           1 |   2020 | Jiang          | Role of PI3K/AKT pathway in cancer: the framework of malignant behavior          |        570 | 10.1007/s11033-020-05435-1   |

### BCR repertoire / TLS in thyroid cancer

|   n_sources |   year | first_author   | title                                                            |   cited_by | doi                        |
|------------:|-------:|:---------------|:-----------------------------------------------------------------|-----------:|:---------------------------|
|           1 |   2023 | Sun            | T cells in health and disease                                    |        874 | 10.1038/s41392-023-01471-y |
|           1 |   2023 | Pisetsky       | Pathogenesis of autoimmune disease                               |        573 | 10.1038/s41581-023-00720-1 |
|           1 |   2022 | Wang           | Therapeutic peptides: current applications and future directions |       1890 | 10.1038/s41392-022-00904-4 |

### BRAF/RAS-negative ('dark matter') PTC subtypes

|   n_sources |   year | first_author   | title                                                                            |   cited_by | doi                        |
|------------:|-------:|:---------------|:---------------------------------------------------------------------------------|-----------:|:---------------------------|
|           1 |   2023 | Bahar          | Targeting the RAS/RAF/MAPK pathway for cancer therapy: from mechanism to clinica |        723 | 10.1038/s41392-023-01705-z |
|           1 |   2023 | Glaviano       | PI3K/AKT/mTOR signaling transduction pathway and targeted therapies in cancer    |       1787 | 10.1186/s12943-023-01827-6 |
|           1 |   2020 | Jiang          | Role of PI3K/AKT pathway in cancer: the framework of malignant behavior          |        570 | 10.1007/s11033-020-05435-1 |

### HLA-II PTC autoimmunity (DPB1*05:01 etc.)

|   n_sources |   year | first_author   | title                                                                            |   cited_by | doi                        |
|------------:|-------:|:---------------|:---------------------------------------------------------------------------------|-----------:|:---------------------------|
|           1 |   2025 | Bryliński      | Effects of Trace Elements on Endocrine Function and Pathogenesis of Thyroid Dise |         24 | 10.3390/nu17030398         |
|           1 |   2021 | Mikosch        | Hashimoto’s thyroiditis and coexisting disorders in correlation with HLA status— |         26 | 10.1007/s10354-021-00879-x |
|           1 |   2023 | Xiong          | How does SARS‐CoV‐2 infection impact on immunity, procession and treatment of pa |         16 | 10.1002/jmv.28487          |

### Hashimoto-like signature & PTC outcomes

|   n_sources |   year | first_author       | title                                                                            |   cited_by | doi                        |
|------------:|-------:|:-------------------|:---------------------------------------------------------------------------------|-----------:|:---------------------------|
|           1 |   2022 | Kłubo-Gwieździńska | Hashimoto thyroiditis: an evidence-based guide: etiology, diagnosis and treatmen |        181 | 10.20452/pamw.16222        |
|           1 |   2023 | Zeng               | Understanding tumour endothelial cell heterogeneity and function from single-cel |        172 | 10.1038/s41568-023-00591-5 |
|           2 |   2021 | Pan                | Papillary Thyroid Carcinoma Landscape and Its Immunological Link With Hashimoto  |         42 | 10.3389/fcell.2021.758339  |

### Pan-Asian HLA fine-mapping (Graves' / autoimmune thyroid)

|   n_sources |   year | first_author   | title                                                                            |   cited_by | doi                        |
|------------:|-------:|:---------------|:---------------------------------------------------------------------------------|-----------:|:---------------------------|
|           1 |   2018 | Rojas          | Molecular mimicry and autoimmunity                                               |        597 | 10.1016/j.jaut.2018.10.012 |
|           1 |   2017 | Matzaraki      | The MHC locus and genetic susceptibility to autoimmune and infectious diseases   |        597 | 10.1186/s13059-017-1207-1  |
|           1 |   2016 | Yazdani        | Selective IgA Deficiency: Epidemiology, Pathogenesis, Clinical Phenotype, Diagno |        197 | 10.1111/sji.12499          |

### Single-cell PTC progression

|   n_sources |   year | first_author   | title                                                                            |   cited_by | doi                          |
|------------:|-------:|:---------------|:---------------------------------------------------------------------------------|-----------:|:-----------------------------|
|           1 |   2023 | Glaviano       | PI3K/AKT/mTOR signaling transduction pathway and targeted therapies in cancer    |       1787 | 10.1186/s12943-023-01827-6   |
|           1 |   2022 | Anand          | Cancer chemotherapy and beyond: Current status, drug candidates, associated risk |       1578 | 10.1016/j.gendis.2022.02.007 |
|           1 |   2022 | Li             | Lactate metabolism in human health and disease                                   |       1243 | 10.1038/s41392-022-01151-3   |

### TIERA-like Korean PTC molecular cohort

|   n_sources |   year | first_author   | title                                                                            |   cited_by | doi                       |
|------------:|-------:|:---------------|:---------------------------------------------------------------------------------|-----------:|:--------------------------|
|           1 |   2020 | Aaltonen       | Pan-cancer analysis of whole genomes                                             |       3240 | 10.1038/s41586-020-1969-6 |
|           1 |   2020 | Yang           | Brief introduction of medical database and data mining technology in big data er |        592 | 10.1111/jebm.12373        |
|           1 |   2018 | Si             | The roles of metallothioneins in carcinogenesis                                  |        442 | 10.1186/s13045-018-0645-x |

---

## 4. HLA AFND — sample-weighted East Asian matrix (Pillar 1 input)

_From 123 (allele × population × study) raw records, computed sample-_
_weighted mean frequency per (allele, country) using AFND `sample_size` as weights._

### Weighted matrix

| allele     | population   |   weighted_freq |   total_n |   n_studies |
|:-----------|:-------------|----------------:|----------:|------------:|
| B*46:01    | China        |       0.117921  |      2584 |          18 |
| B*46:01    | Japan        |       0.0390204 |       538 |           3 |
| B*46:01    | Korea        |       0.044     |       485 |           1 |
| B*46:01    | Taiwan       |       0.118773  |      1583 |           7 |
| DPB1*05:01 | China        |       0.35543   |      1969 |          20 |
| DPB1*05:01 | Japan        |       0.385741  |      1235 |           6 |
| DPB1*05:01 | Korea        |       0.366737  |       680 |           3 |
| DPB1*05:01 | Taiwan       |       0.521446  |       704 |           5 |
| DQB1*06:02 | China        |       0.0667075 |      2759 |          15 |
| DQB1*06:02 | Japan        |       0.0770119 |       421 |           2 |
| DRB1*04:05 | China        |       0.0537887 |      2645 |          12 |
| DRB1*04:05 | Japan        |       0.117969  |       421 |           2 |
| DRB1*04:05 | Korea        |       0.06      |       201 |           1 |
| DRB1*15:01 | China        |       0.0924874 |      2795 |          13 |
| DRB1*15:01 | Korea        |       0.09      |       201 |           1 |

### Pillar 1 forest input — DPB1*05:01 East Asian baseline

| Country | Weighted freq | Total n | Studies | vs Korean PTC pool 53.2% |
|---|---|---|---|---|
| China | 0.355 (35.5%) | 1969 | 20 | Δ = +0.177 |
| Japan | 0.386 (38.6%) | 1235 | 6 | Δ = +0.146 |
| Korea | 0.367 (36.7%) | 680 | 3 | Δ = +0.165 |
| Taiwan | 0.521 (52.1%) | 704 | 5 | Δ = +0.011 |

_**Interpretation:** Korean PTC pool 53.2% > Korean baseline 36.7% (Δ=+16.5pp,_
_pooled n=680). Forest meta input ready._

---

## 5. Clinical trials — structured drill-down

_Total 65 active/recruiting trials across 5 queries; **38 use targeted_
_drugs** (RET/BRAF/IO/TKI). Phase + intervention + sponsor extracted from CT.gov v2._

### Phase distribution (all 65 trials)

| phase         |   n_trials |
|:--------------|-----------:|
| PHASE2        |         32 |
| PHASE1        |          5 |
| PHASE1,PHASE2 |          5 |
| PHASE4        |          3 |
| PHASE3        |          3 |
| EARLY_PHASE1  |          2 |
| PHASE2,PHASE3 |          1 |

### Sponsor class (all 65 trials)

| class    |   n_trials |
|:---------|-----------:|
| OTHER    |         46 |
| INDUSTRY |         13 |
| NIH      |          5 |
| NETWORK  |          1 |

### Top targeted-drug trials (sorted: phase desc, recruiting first)

| nct_id      | phase         | status                |   enrollment | intervention_names                                                               | lead_sponsor                           |
|:------------|:--------------|:----------------------|-------------:|:---------------------------------------------------------------------------------|:---------------------------------------|
| NCT04940052 | PHASE3        | ACTIVE_NOT_RECRUITING |          153 | Dabrafenib; Trametinib; Trametinib Placebo; Dabrafenib placebo                   | Novartis Pharmaceuticals               |
| NCT02465060 | PHASE2        | ACTIVE_NOT_RECRUITING |         6452 | Adavosertib; Afatinib; Afatinib Dimaleate; Binimetinib; Biopsy Procedure; Biospe | National Cancer Institute (NCI)        |
| NCT02834013 | PHASE2        | ACTIVE_NOT_RECRUITING |          818 | Biospecimen Collection; Computed Tomography; Echocardiography Test; Ipilimumab;  | National Cancer Institute (NCI)        |
| NCT02834013 | PHASE2        | ACTIVE_NOT_RECRUITING |          818 | Biospecimen Collection; Computed Tomography; Echocardiography Test; Ipilimumab;  | National Cancer Institute (NCI)        |
| NCT04280081 | PHASE2        | ACTIVE_NOT_RECRUITING |           77 | Selpercatinib                                                                    | Eli Lilly and Company                  |
| NCT05696548 | PHASE2        | ACTIVE_NOT_RECRUITING |           51 | Lenvatinib; Nivolumab                                                            | National Cancer Center Hospital East   |
| NCT03181100 | PHASE2        | ACTIVE_NOT_RECRUITING |           50 | Atezolizumab; Bevacizumab; Cobimetinib; Nab-paclitaxel; Paclitaxel; Vemurafenib  | M.D. Anderson Cancer Center            |
| NCT03181100 | PHASE2        | ACTIVE_NOT_RECRUITING |           50 | Atezolizumab; Bevacizumab; Cobimetinib; Nab-paclitaxel; Paclitaxel; Vemurafenib  | M.D. Anderson Cancer Center            |
| NCT04439292 | PHASE2        | ACTIVE_NOT_RECRUITING |           35 | Dabrafenib Mesylate; Trametinib Dimethyl Sulfoxide                               | National Cancer Institute (NCI)        |
| NCT04675710 | PHASE2        | ACTIVE_NOT_RECRUITING |           30 | Conventional Surgery; Dabrafenib; Intensity-Modulated Radiation Therapy; Pembrol | M.D. Anderson Cancer Center            |
| NCT03449108 | PHASE2        | ACTIVE_NOT_RECRUITING |           30 | Aldesleukin; Autologous Tumor Infiltrating Lymphocytes LN-145; Autologous Tumor  | M.D. Anderson Cancer Center            |
| NCT04759911 | PHASE2        | ACTIVE_NOT_RECRUITING |           30 | Quality-of-Life Assessment; Questionnaire Administration; Selpercatinib; Therape | M.D. Anderson Cancer Center            |
| NCT04675710 | PHASE2        | ACTIVE_NOT_RECRUITING |           30 | Conventional Surgery; Dabrafenib; Intensity-Modulated Radiation Therapy; Pembrol | M.D. Anderson Cancer Center            |
| NCT04171622 | PHASE2        | ACTIVE_NOT_RECRUITING |           25 | Lenvatinib; Pembrolizumab                                                        | M.D. Anderson Cancer Center            |
| NCT04061980 | PHASE2        | ACTIVE_NOT_RECRUITING |           24 | Binimetinib; Encorafenib; Nivolumab                                              | Providence Health & Services           |
| NCT04238624 | PHASE2        | ACTIVE_NOT_RECRUITING |           16 | Dabrafenib; Trametinib                                                           | Memorial Sloan Kettering Cancer Center |
| NCT04238624 | PHASE2        | ACTIVE_NOT_RECRUITING |           16 | Dabrafenib; Trametinib                                                           | Memorial Sloan Kettering Cancer Center |
| NCT05119296 | PHASE2        | ACTIVE_NOT_RECRUITING |           12 | Pembrolizumab (Keytruda)                                                         | Stanford University                    |
| NCT03914300 | PHASE2        | ACTIVE_NOT_RECRUITING |           11 | Biospecimen Collection; Cabozantinib S-malate; Computed Tomography; Ipilimumab;  | National Cancer Institute (NCI)        |
| NCT03899792 | PHASE1,PHASE2 | ACTIVE_NOT_RECRUITING |           36 | Selpercatinib                                                                    | Eli Lilly and Company                  |

---

## 6. Open-access full-text discovery (Unpaywall)

_10/17 entries have OA copy. PMC IDs + PDF URLs ready for Methods/Discussion citations._

| key                | is_oa   | oa_status   | pmcid       | pdf_url                                                                          |
|:-------------------|:--------|:------------|:------------|:---------------------------------------------------------------------------------|
| Haugen2016         | True    | bronze      | PMC4739132  | https://www.liebertpub.com/doi/pdf/10.1089/thy.2015.0020                         |
| Ringel2025         | True    | green       |             |                                                                                  |
| TCGA2014           | True    | bronze      | PMC4243044  | http://www.cell.com/article/S0092867414012380/pdf                                |
| Yoo2016            | True    | gold        | PMC4986964  | https://journals.plos.org/plosgenetics/article/file?id=10.1371/journal.pgen.1006 |
| Landa2016          | True    | bronze      | PMC4767360  | http://www.jci.org/articles/view/85271/files/pdf                                 |
| Pu2021             | True    | gold        | PMC8523608  | https://www.nature.com/articles/s41467-021-26343-3.pdf                           |
| Wirth2020          | True    | bronze      | PMC10777663 | https://www.nejm.org/doi/pdf/10.1056/NEJMoa2005651?articleTools=true             |
| Xing2014           | False   | closed      | PMC3791171  |                                                                                  |
| Liu2017            |         | no_doi      |             |                                                                                  |
| Wang2024           |         | no_doi      |             |                                                                                  |
| Cibas2017          | True    | bronze      |             | https://www.liebertpub.com/doi/pdf/10.1089/thy.2017.0500                         |
| SEER_PTC           |         | no_doi      |             |                                                                                  |
| Lu2023             |         | no_doi      |             |                                                                                  |
| Bradley2010        |         | no_doi      |             |                                                                                  |
| Krishnamoorthy2025 | True    | green       | PMC11849553 |                                                                                  |
| Chu2018            |         | no_doi      |             |                                                                                  |
| Pozdeyev2018       | True    | green       | PMC6030480  |                                                                                  |

---

## 7. Cross-reference: dedup ∩ verified_bib

_2 papers appear in BOTH our bib AND the mined top-ranked list — these are_
_**self-validated citations** (independently surfaced as top by 5-source mining)._

| claim                                       | doi                        | title                                                                            | first_author   |   year |
|:--------------------------------------------|:---------------------------|:---------------------------------------------------------------------------------|:---------------|-------:|
| 8-gene / driver-excluded PTC stratification | 10.1038/s41467-021-26343-3 | Single-cell transcriptomic analysis of the tumor ecosystems underlying initiatio | Pu             |   2021 |
| Single-cell PTC progression                 | 10.1038/s41467-021-26343-3 | Single-cell transcriptomic analysis of the tumor ecosystems underlying initiatio | Pu             |   2021 |

---

## 8. Action items for manuscript_v8

**Immediate (paste-ready):**
1. Pillar 1 forest: replace placeholder Korean baseline with **36.7% (n=680, 3 studies)** — §4.
2. Pan-Asian gradient panel: DPB1*05:01 Korea 36.7% / China 35.5% / Japan 38.6% / Taiwan 52.1%
3. Discussion §3.x translational outlook: cite NCT06458036 (RAISE selpercatinib pre-RAI),
   NCT06475989 (Phase 3 targeted vs chemo), NCT04675710 (pembro+dabra+trame neoadjuvant).

**Triage (review gap_candidates_top3 per claim):**
- 24 high-rank papers not yet in bib — accept/reject per claim.

**Pending (next pass):**
- Fetch full abstracts for 10 OA Unpaywall PDFs (currently link-only).
- Sem Scholar 0-hit on 6/8 claims (rate-limited) — retry with backoff.
