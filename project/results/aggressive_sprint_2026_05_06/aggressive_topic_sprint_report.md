# Aggressive-topic feasibility sprint

Date: 2026-05-06

Scope: rapid local + public-download tests; no Paper 1 edits.

## Bottom line

The strongest aggressive new angle is **not** generic dark-matter taxonomy. The best candidate from today's test is a **copy-number-defined residual class inside true driver-negative / DM2-like tumors**, with a secondary high-risk option of hidden fusion/regulatory mechanism after fusion-call reconciliation. The direct RAI-refractory paper is currently weak after re-scoring GSE151179.

## 1. Direct RAI-refractory validation: weak / no-go as headline

GSE151179 parsed from GEO: 52 samples, 39 tumors, 17 primary tumors. Probe recovery: {'SLC5A5': 1, 'TPO': 1, 'TSHR': 1, 'TG': 1, 'PAX8': 1, 'NKX2-1': 1, 'FOXE1': 1, 'DIO1': 1}.

Best refractory-response tests by oriented AUC:

| subset | score | n | n_pos | n_neg | mean_pos | mean_neg | cohens_d_pos_vs_neg | mw_p | auc_pos_high | auc_best_direction | lower_score_predicts_pos |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| tumor_all | TPO | 39 | 35 | 4 | 8.344 | 9.084 | -0.4563 | 0.3812 | 0.3571 | 0.6429 | True |
| pre_rai_tumor | TSHR | 22 | 18 | 4 | 9.228 | 9.128 | 0.06532 | 0.4342 | 0.6389 | 0.6389 | False |
| tumor_all | DIO1 | 39 | 35 | 4 | 5.234 | 5.345 | -0.105 | 0.4066 | 0.3643 | 0.6357 | True |
| primary_only | TSHR | 17 | 13 | 4 | 9.37 | 9.128 | 0.2066 | 0.4773 | 0.6346 | 0.6346 | False |
| tumor_all | rai6_raw_mean | 39 | 35 | 4 | 8.447 | 8.684 | -0.2579 | 0.4328 | 0.3714 | 0.6286 | True |
| tumor_all | rai8_raw_mean | 39 | 35 | 4 | 7.928 | 8.146 | -0.2515 | 0.4328 | 0.3714 | 0.6286 | True |
| tumor_all | TG | 39 | 35 | 4 | 10.73 | 11.06 | -0.1883 | 0.4328 | 0.3714 | 0.6286 | True |
| primary_only | PAX8 | 17 | 13 | 4 | 9.141 | 9.338 | -0.1841 | 0.5487 | 0.3846 | 0.6154 | True |


Best metastatic RAI-uptake tests by oriented AUC:

| subset | score | n | n_pos | n_neg | mean_pos | mean_neg | cohens_d_pos_vs_neg | mw_p | auc_pos_high | auc_best_direction | lower_score_predicts_pos |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| primary_only | NKX2-1 | 17 | 6 | 11 | 7.683 | 7.446 | 0.7092 | 0.149 | 0.7273 | 0.7273 | False |
| pre_rai_tumor | TPO | 22 | 8 | 14 | 8.345 | 9.351 | -0.6048 | 0.1653 | 0.3125 | 0.6875 | True |
| primary_only | TPO | 17 | 6 | 11 | 8.449 | 9.287 | -0.5328 | 0.3011 | 0.3333 | 0.6667 | True |
| tumor_all | TPO | 39 | 20 | 19 | 7.971 | 8.893 | -0.5887 | 0.1189 | 0.3526 | 0.6474 | True |
| primary_only | PAX8 | 17 | 6 | 11 | 9.16 | 9.202 | -0.03938 | 0.4043 | 0.6364 | 0.6364 | False |
| pre_rai_tumor | TG | 22 | 8 | 14 | 9.926 | 11.17 | -0.644 | 0.3301 | 0.3661 | 0.6339 | True |
| metastatic_tumor | DIO1 | 22 | 14 | 8 | 4.777 | 5.206 | -0.5687 | 0.402 | 0.3839 | 0.6161 | True |
| primary_only | rai6_z_mean | 17 | 6 | 11 | -0.1344 | -0.09258 | -0.07182 | 0.5249 | 0.6061 | 0.6061 | False |


Interpretation: if lower thyroid-lineage/RAI score predicts refractory/no-uptake, the signal is at best modest and sample-imbalanced. This can support a limitation or pilot validation, not an aggressive standalone RAI-refractory classifier.

## 2. True driver-negative / hidden mechanism: most interesting, but validation-gated

Driver-class score summary:

| driver_class | n | DM1_n | DM2_n | PFI_events | PFI_rate | rai_mean | tds16_mean | methyl8_mean | age_mean |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Class1_BRAF_V600E | 291 | 1 | 0 | 34 | 0.1168 | 7.035 | 6.327 | 0.3719 | 48.28 |
| Class2_RAS_hotspot | 54 | 0 | 0 | 7 | 0.1296 | 8.914 | 7.884 | 0.2715 | 46.05 |
| Class4_DICER1_EIF1AX | 7 | 1 | 6 | 1 | 0.1429 | 9.031 | 7.961 | 0.2365 | 52.59 |
| Class5_TERT_only | 5 | 4 | 1 | 2 | 0.4 | 6.791 | 6.163 | 0.3556 | 54.37 |
| Class6_True_driver_neg | 125 | 76 | 48 | 6 | 0.048 | 8.327 | 7.41 | 0.3359 | 46.35 |


Methylation-expression correlations in DM1/DM2:

| x | y | n | spearman_r | p |
| --- | --- | --- | --- | --- |
| mean_8g_beta | rai_score_v17 | 135 | -0.5516 | 4.102e-12 |
| mean_8g_beta | tds16_score_v17 | 135 | -0.583 | 1.184e-13 |


Existing methylation result: DM1 is hypermethylated across the 8 thyroid-lineage/RAI genes versus DM2; TPO shows delta beta 0.415, Cohen d 2.299, p=1.88e-18. This is a real mechanism signal, but it overlaps heavily with the Paper 1 differentiation axis.

Class6 arm-level CNV enrichment top rows:

| arm | call | class6_call_pct | other_call_pct | odds_ratio | fisher_p | bh_q |
| --- | --- | --- | --- | --- | --- | --- |
| 7q | Gain | 10 | 0.295 | 37.56 | 1.225e-06 | 9.556e-05 |
| 7p | Gain | 9.91 | 0.2941 | 37.29 | 1.292e-06 | 5.04e-05 |
| 2q | Loss | 8.108 | 0 |  | 2.574e-06 | 6.692e-05 |
| 2p | Loss | 8.108 | 0 |  | 2.574e-06 | 5.019e-05 |
| 16p | Gain | 9.009 | 0.2941 | 33.56 | 5.174e-06 | 8.071e-05 |
| 12q | Gain | 7.207 | 0 |  | 1.107e-05 | 0.0001439 |
| 16q | Gain | 8.108 | 0.2941 | 29.91 | 2.039e-05 | 0.0002272 |
| 21q | Loss | 7.207 | 0.2941 | 26.33 | 7.905e-05 | 0.0007707 |
| 12p | Gain | 5.455 | 0 |  | 0.000192 | 0.001664 |
| 17q | Gain | 7.207 | 0.8824 | 8.725 | 0.000913 | 0.007122 |


Interpretation: this is a real signal, not just noise. The top Class6 arm-level CNV enrichments remain significant after the simple BH screen, and the pattern is concentrated in Class6/DM2 rather than Class6/DM1. A seven-arm quick signature using 7q gain, 7p gain, 2q loss, 2p loss, 16p gain, 16q gain, and 12q gain marks **14/45 evaluable Class6-DM2** versus **1/65 evaluable Class6-DM1** samples. This is the strongest "new" aggressive topic because it is genomic, not just the Paper 1 lineage score.

## 3. Fusion angle: high upside, current evidence inconsistent

PanCancer/SV-style local fusion result reports DM1 fusion-positive 76.8%, DM2 30.9%, not-DM 12.8%. RET/NTRK/ALK/BRAF fusion classes are concentrated in DM1 in that table.

But v17 curated fusion overlay is sparse: DM1 has ALK 1, NTRK 1, RET 1; DM2 has 0. Therefore this is **the most aggressive testable new topic**, but it is not publishable until raw RNA/SV calls are reconciled and externally validated.

## 4. DICER1/EIF1AX class: real but small

TCGA class counts: {'Class1_BRAF_V600E': 291, 'Class6_True_driver_neg': 125, 'Class2_RAS_hotspot': 54, 'Class4_DICER1_EIF1AX': 7, 'Class5_TERT_only': 5}. K2/Yoo: BRAF/RAS-negative 68/180 = 37.78%; DICER1+EIF1AX in K2 dark matter = 7 cases.

Interpretation: biologically clean, replicated frequency, but n=7 in TCGA and n=7 in K2 dark-matter. Good as a class inside a broader paper, weak as a standalone aggressive paper.

## 5. TERT-only class: too small right now

Local class map: TERT-only n=5, PFI events 2/5. This is provocative but too small for a standalone paper without external TERT-rich cohort access.

## Ranked recommendation

| rank | topic | why | blocker | decision |
| --- | --- | --- | --- | --- |
| 1 | Copy-number-defined true-driver-negative DM2 class | Strongest new signal; genomic; less Paper 1 overlap | Needs validation against cBioPortal/GDC GISTIC source and external cohort | Attack first |
| 2 | DM1 hidden fusion/regulatory driver | Highest conceptual novelty | Fusion evidence conflicts by source; needs raw fusion/SV reconciliation | Attack second |
| 3 | True driver-negative epigenetic silencing | Strong methylation signal | Overlaps Paper 1 lineage axis | Only if framed as mechanism, not biomarker |
| 4 | DICER1/EIF1AX FVPTC-like class | Replicated rare-driver enrichment | Small n | Use as section, not standalone |
| 5 | RAI-refractory classifier | Direct public labels available | GSE151179 signal weak/imbalanced | No-go as headline today |
| 6 | TERT-only triple-negative | Provocative event rate | n=5 | Hold until external cohort |

## Immediate next work if continuing today

1. Validate the Class6/DM2 CNV signature from independent GDC/cBioPortal GISTIC calls and check whether it corresponds to whole-chromosome/arm aneuploidy rather than focal artifacts.  
2. Reconcile fusion evidence: compare cBioPortal PanCancer structural variants, v17 fusion overlay, and raw GSE184362/GSE232237 metadata.  
3. If raw FASTQ or STAR-Fusion-compatible files exist for K2/GSE184362, run fusion calling on BRAF/RAS-negative or DM1-like samples.  
4. For epigenetic mechanism, test whether promoter methylation explains expression after adjusting for histology and BRAF/RAS class.  
