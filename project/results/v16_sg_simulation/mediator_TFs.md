# BRAF V600E → TF → TACSTD2 — mediator inference

_In silico reasoning based on local DE (biomarker_validated.tsv) + Claude sonnet_

## Top mediator TFs
### FOSL2 (confidence: high)
- **Why**: FOSL2 is a key AP-1 component directly activated by MAPK/ERK signaling and shows strong upregulation (log2FC=1.12, d=1.38) matching TACSTD2's BRAF-like pattern. AP-1 complexes are well-established regulators of epithelial markers and cell adhesion molecules.
- **Literature anchor**: Murphy et al. Nature Reviews Cancer 2013 - AP-1 in cancer progression
- **Our DE**: log2FC=1.12, Cohen d=1.38

### STAT1 (confidence: high)
- **Why**: STAT1 shows the strongest upregulation among STATs (log2FC=1.47, d=1.58) and can be activated downstream of MAPK through interferon signaling crosstalk. STAT1 regulates epithelial differentiation genes and cell surface markers in thyroid cancer.
- **Literature anchor**: Kimura et al. Endocr Relat Cancer 2003 - STAT signaling in thyroid
- **Our DE**: log2FC=1.47, Cohen d=1.58

### STAT3 (confidence: high)
- **Why**: STAT3 is moderately upregulated (log2FC=0.49, d=1.22) and represents a key node where MAPK and cytokine signaling converge. STAT3 directly regulates epithelial cell adhesion molecules and is frequently activated in BRAF-mutant thyroid cancers.
- **Literature anchor**: Borrello et al. Oncogene 2005 - STAT3 in thyroid transformation
- **Our DE**: log2FC=0.49, Cohen d=1.22

### RELA (confidence: medium)
- **Why**: RELA (p65 NF-κB subunit) shows consistent upregulation (log2FC=0.34, d=1.20) and is activated downstream of MAPK through multiple mechanisms. NF-κB directly regulates cell surface receptors and adhesion molecules including TACSTD family members.
- **Literature anchor**: Karin & Ben-Neriah Nature 2000 - NF-κB in cancer
- **Our DE**: log2FC=0.34, Cohen d=1.20

### TEAD1 (confidence: medium)
- **Why**: TEAD1 is upregulated (log2FC=0.67, d=1.14) and serves as a key Hippo pathway effector that can be modulated by MAPK signaling. TEAD factors regulate epithelial cell identity genes and surface markers in thyroid epithelial cells.
- **Literature anchor**: Pan et al. Genes Dev 2018 - Hippo-MAPK crosstalk
- **Our DE**: log2FC=0.67, Cohen d=1.14

## Hypothesised path
The BRAF V600E → MAPK → TF → TACSTD2 axis likely operates through multiple parallel transcriptional programs, with FOSL2/AP-1 and STAT1/3 representing the strongest mediators based on their expression patterns and known MAPK connectivity. These TFs coordinate to drive the epithelial differentiation program characteristic of BRAF-like thyroid cancers, with TACSTD2 serving as a key downstream effector of this epithelial identity network.

## Caveat
This is a *literature- and DE-informed hypothesis*, not a CellOracle-inferred GRN. CellOracle / SCENIC / scenic+ would require a thyroid-specific scATAC or TF-motif prior that was not run in this environment. The mediator list is publication-useful as a starting hypothesis for experimental follow-up (ChIP-seq, knockdown, reporter assay).