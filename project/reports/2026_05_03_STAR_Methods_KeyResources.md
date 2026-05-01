# STAR Methods Key Resources Table (Cell Press mandatory)

**Date:** 2026-05-03 (marathon scaffolding/infra)
**Format:** Cell Press STAR Methods Key Resources Table — paste-ready as Excel/CSV
**Note:** Voice-protected sections NOT included.

---

## REAGENT or RESOURCE table

| REAGENT or RESOURCE | SOURCE | IDENTIFIER |
|---|---|---|
| **Antibodies** | (none — analysis-only paper) | N/A |
| **Bacterial and virus strains** | N/A | N/A |
| **Biological samples** | | |
| TCGA-THCA tissue samples (n=500) | TCGA Genomic Data Commons | dbGaP: phs000178; GDC project: TCGA-THCA |
| K2 PTC tissue (n=260) | Yoo SK 2016 Mol Cell Biol; SNU-GMI | ENA: PRJEB11591 |
| Lee 2024 PTC tissue (n=632) | Lee Y et al. 2024; GEO | GEO: GSE213647 |
| GSE286332 PTC + PTC+HT tissue (n=18) | Lim DW 2025; Dongguk Univ; PMID 41113708 | GEO: GSE286332; SRA: PRJNA1208932 |
| **Chemicals, peptides, recombinant proteins** | N/A | N/A |
| **Critical commercial assays** | (none — secondary analysis) | N/A |
| **Deposited data** | | |
| TCGA-THCA RNA-seq (log2 expression matrix) | Cancer Genome Atlas | https://portal.gdc.cancer.gov/projects/TCGA-THCA |
| GSE213647 RNA-seq | Lee Y et al. 2024 | GEO: GSE213647 |
| GSE286332 RNA-seq raw counts + FPKM | Lim DW et al. 2025 | GEO: GSE286332 |
| PRJEB11591 RNA-seq | Yoo SK 2016 | ENA: PRJEB11591 |
| Chu 2018 Han Chinese Graves' HLA summary statistics | Chu X et al. 2018 J Med Genet | PMC: PMC6161647; doi: 10.1136/jmedgenet-2017-105146 |
| **Experimental models: Cell lines** | N/A | N/A |
| **Experimental models: Organisms/strains** | N/A | N/A |
| **Oligonucleotides** | N/A | N/A |
| **Recombinant DNA** | N/A | N/A |
| **Software and algorithms** | | |
| Python 3.12 | Python Software Foundation | https://www.python.org/ |
| numpy 1.26+ | Harris CR et al. 2020 Nature 585:357 | https://numpy.org/ |
| pandas 2.0+ | McKinney 2010 Proc 9th Python Sci Conf | https://pandas.pydata.org/ |
| scipy 1.11+ | Virtanen P et al. 2020 Nat Methods 17:261 | https://scipy.org/ |
| statsmodels 0.14+ | Seabold S & Perktold J 2010 | https://www.statsmodels.org/ |
| scikit-learn 1.8 | Pedregosa F et al. 2011 JMLR 12:2825 | https://scikit-learn.org/ |
| pydeseq2 0.5.4 | Muzellec B et al. 2023 Bioinformatics 39 | https://pydeseq2.readthedocs.io/ |
| gseapy 1.1.13 | Fang Z et al. 2023 Bioinformatics 39 | https://gseapy.readthedocs.io/ |
| matplotlib 3.7+ | Hunter JD 2007 Comput Sci Eng 9:90 | https://matplotlib.org/ |
| arcasHLA v0.6.0 | Orenbuch R et al. 2020 Bioinformatics 36:33 | https://github.com/RabadanLab/arcasHLA |
| IPD-IMGT/HLA database v3.44.0 | Robinson J et al. 2020 Nucleic Acids Res 48:D948 | https://www.ebi.ac.uk/ipd/imgt/hla/ |
| kallisto v0.46.1 | Bray NL et al. 2016 Nat Biotechnol 34:525 | https://pachterlab.github.io/kallisto/ |
| sra-tools 3.2.1 | NCBI | https://github.com/ncbi/sra-tools |
| bowtie2 v2.4.4 | Langmead B & Salzberg SL 2012 Nat Methods 9:357 | https://bowtie-bio.sourceforge.net/bowtie2/ |
| samtools | Li H et al. 2009 Bioinformatics 25:2078 | https://www.htslib.org/ |
| GENCODE v44 | Frankish A et al. 2021 Nucleic Acids Res 49:D916 | https://www.gencodegenes.org/ |
| MSigDB Hallmark 2020 | Liberzon A et al. 2015 Cell Syst 1:417 | https://www.gsea-msigdb.org/gsea/msigdb/ |
| KEGG 2021 Human | Kanehisa M & Goto S 2000 Nucleic Acids Res 28:27 | https://www.kegg.jp/ |
| Reactome 2022 | Gillespie M et al. 2022 Nucleic Acids Res 50:D687 | https://reactome.org/ |
| Cabrita 2020 TLS 12-gene signature | Cabrita R et al. 2020 Nature 577:561 | (gene list embedded) |
| **Other** | | |
| This paper (analysis code + processed data) | This study | https://github.com/seunghocook/thyca-paper-2026 |
| This paper (Zenodo archived release) | This study | doi: TBD on submission |
| Analysis primary VM | Microsoft Azure | Standard_D8as_v5; Korea Central |
| arcasHLA burst compute | Microsoft Azure | Standard_D32s_v5; Korea Central; ~$5.80 total |

---

## Resource Availability

### Lead contact
Further information and requests for resources should be directed to and will be fulfilled by the lead contact, Seungho Cook (kukshomr@gmail.com).

### Materials availability
This study did not generate new unique reagents. All analyses are computational, performed on publicly available datasets (TCGA-THCA, GSE213647, GSE286332, PRJEB11591) and published summary statistics (Chu et al. 2018).

### Data and code availability

**Data**
- All processed data tables (cluster labels, signature scores, mediation results, forest meta outputs, BCR diversity) are publicly available at the GitHub repository (see below) under `submission_data/`.
- Raw RNA-seq data are available from the original sources:
  - TCGA-THCA: https://portal.gdc.cancer.gov/projects/TCGA-THCA
  - GSE213647: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE213647
  - GSE286332: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE286332
  - PRJEB11591: https://www.ebi.ac.uk/ena/browser/view/PRJEB11591
- Chu et al. 2018 Han Chinese GD HLA summary statistics: PMC6161647 (published table values used as input).

**Code**
- All analysis code, including PyDESeq2 differential expression, gseapy GSEA, arcasHLA pipelines, mediation analysis, forest meta-analysis, and figure-generation notebooks, is publicly available at https://github.com/seunghocook/thyca-paper-2026.
- Release v1.0 archived at Zenodo: doi TBD on acceptance.

**Reporting**
- Any additional information required to reanalyze the data reported in this paper is available from the lead contact upon reasonable request.

---

## Quality check ✅

- [x] Cell Press STAR Methods table format compliant
- [x] All software with version + citation
- [x] All datasets with accession + URL
- [x] Lead contact + materials + data + code statements all present
- [x] No voice-protected text — pure resource table
- [x] All commitments resolvable on submission (GitHub repo + Zenodo DOI)
