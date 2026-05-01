============================================================
[1] TCGA-THCA cbio mutation table — K601E + V600E count
   total rows: 11857
   BRAF rows: 304, unique patients: 287

   per-patient BRAF class:
braf_class
V600E                282
other                  2
K601E                  2
V600_syn_or_other      1

   K601E patient IDs: ['TCGA-EM-A3O8', 'TCGA-ET-A4KQ']

   K601E in sample_master:
     TCGA-ET-A4KQ  TERT=wildtype  age=45.7  os_event=0.0  os_days=955.0  histology=FVPTC
     TCGA-EM-A3O8  TERT=wildtype  age=34.0  os_event=0.0  os_days=381.0  histology=FVPTC

============================================================
[2] MSK external cohorts — TERT promoter + BRAF positions
    (fetching mutation tables from cBioPortal datahub)

   Known external promoter records (from v2 sprint):
     total: 1880

   external TERT prevalence by cohort + position:
pos_label           C228T  C250T  other
cohort_short                           
nan                  1280    352    188
thyroid_mskcc_2016     49      9      2

   fetching thyroid_mskcc_2016 full mutation table...
     fetching https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/thyro...
     loaded: (564, 46), cols sample: ['Hugo_Symbol', 'Entrez_Gene_Id', 'Center', 'NCBI_Build', 'Chromosome', 'Start_Position', 'End_Position', 'Strand', 'Consequence', 'Variant_Classification']
     MSK thyroid BRAF rows: 44
     MSK BRAF HGVSp_Short:
HGVSp_Short
p.V600E    43
p.I300V     1
     saved to /opt/thyroid-dash/project/results/v17_tert_recovery/v3/v3_thyroid_mskcc_2016_mutations.tsv.gz

   fetching thyroid_mskcc_2016 clinical...
     fetching https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/thyro...
     fetching https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/thyro...
     sample clinical: (117, 17), cols: ['PATIENT_ID', 'SAMPLE_ID', 'OTHER_SAMPLE_ID', 'SAMPLE_CLASS', 'CANCER_TYPE', 'CANCER_TYPE_DETAILED', 'ONCOTREE_CODE', 'SAMPLE_TYPE', 'SITE_OF_RECURRENCE', 'SAMPLE_TYPE_DETAIL']
     patient clinical: (117, 9), cols: ['PATIENT_ID', 'SEX', 'OS_STATUS', 'OS_MONTHS', 'AGE', 'PDTC_DEFINITION', 'M_STAGE', 'PATH_N_STAGE', 'PATH_T_STAGE']
     OS_STATUS distribution:
       OS_STATUS: {'0:LIVING': 68, '1:DECEASED': 47}
       OS_MONTHS: {10.3: 2, 2.83: 2, 2.27: 1, 3.22: 1, 1.81: 1}

============================================================
[3] Combined BRAF+TERT analysis across cohorts
   MSK TERT rows: 61
   MSK TERT samples (mutation table + ext promoter records): 59

   COMBINED table: (545, 8)
   cohort × braf_class × tert breakdown:
tert                                mutated  wildtype
cohort           braf_class                          
MSK-thyroid-2016 V600E                   26        17
TCGA-THCA        K601E                    0         2
                 RAS_only                 6        48
                 TripleNeg                4       157
                 V600E                   26       256
                 V600_syn_or_other        0         1
                 other                    0         2

   wrote /opt/thyroid-dash/project/results/v17_tert_recovery/v3/v3_external_combined_BRAF_TERT.tsv

============================================================
[4] Combined BRAF×TERT survival
   n with OS: 540, events: 38

   V600E: n=322, ev=30
     V600E+TERT+: n=52, ev=23
     V600E+TERT-: n=270, ev=7
     logrank V600E+TERT vs V600E-only: chi2=89.70, p=2.77e-21

   K601E (combined cohorts): n=2, ev=0
     K601E+TERT+: n=0
     K601E+TERT-: n=2

   RAS_only: n=54, ev=1
     RAS+TERT+: n=6, ev=0
     RAS+TERT-: n=48, ev=1