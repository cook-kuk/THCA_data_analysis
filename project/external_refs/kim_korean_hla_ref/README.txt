###########################################################
# Korean Reference Panel v1.0 for imputing HLA variants
# Contact: 
#   Kwangwoo Kim (kkim@khu.ac.kr)
#   Sang-Cheol Bae (scbae@hanyang.ac.kr)
# revised on Aug 23, 2025
###########################################################

Thank you for downloading the Korean imputation reference panel.

Our HLA reference panel is designed for use with SNP2HLA. 
Please visit http://www.broadinstitute.org/mpg/snp2hla/ to learn more about SNP2HLA.

A total of 413 unrelated Korean subjects were analyzed for MHC SNPs 
within the extended MHC locus and for classical alleles of six HLA genes: 
HLA-A, -B, -C, -DRB1, -DPB1, and -DQB1. 

The HLA reference panel was constructed by phasing:
  - 5,858 MHC SNPs
  - 233 classical HLA alleles
  - 1,387 amino acid residue markers from 1,025 amino acid positions 
    (as binary variables).

The Korean HLA reference panel is highly applicable and suitable 
for various genome-wide array datasets from East Asian populations, 
including Han Chinese, Japanese, and Koreans.

-----------------------------------------------------------
Citation
-----------------------------------------------------------
Kwangwoo Kim, So-Young Bang, Hye-Soon Lee, Sang-Cheol Bae. (2014)
Construction and Application of a Korean Reference Panel for Imputing 
Classical Alleles and Amino Acids of Human Leukocyte Antigen Genes. 
PLoS ONE. 9(11), e112546. doi: 10.1371/journal.pone.0112546. [November 2014] 

-----------------------------------------------------------
Usage
-----------------------------------------------------------
./SNP2HLA.csh DATA (.bed/.bim/.fam) REFERENCE (.bgl.phased/.markers) OUTPUT plink {optional: max_memory[mb] window_size}

-----------------------------------------------------------
Example
-----------------------------------------------------------
./SNP2HLA.csh /path_to/KOR_REF_1.0/HapMap3_CHB_JPT/hapmap3_r2_b36_chr6.MHC._plink.id61 /path_to/KOR_REF_1.0/KOR_REF/Kim_KOR_HLA OUTPUT/imputed plink

