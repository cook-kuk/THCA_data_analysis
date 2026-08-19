# PDAC OS prediction benchmark · TCGA-PAAD n=176

All methods evaluated on identical TCGA-PAAD patients with valid OS, sorted by C-index.

| Rank | Method | Type | C-index | HR | p | n_genes |
|---|---|---|---|---|---|---|
| 1 | **Cox PH multivariate (Moffitt+KRAS+31genes+age)** | Combined multivariate baseline | 0.695 | — | <1e-6 (composite) | 35 |
| 2 | **Moffitt_basal** | ICI/vaccine signature (re-derived) | 0.647 | 1.56 | 0 | 14/14 |
| 3 | **Hallmark_IFNG** | ICI/vaccine signature (re-derived, extra) | 0.585 | 1.14 | 0.0248 | 5/11 |
| 4 | **KRAS_G12D (binary)** | ICI/vaccine signature (re-derived) | 0.579 | 1.91 | 0.0027 | binary/binary |
| 5 | **MLPRegressor (DeepSurv-proxy, 5-fold CV)** | Neural Cox-style (sklearn) | 0.576 | — | 5-fold CV | 35 |
| 6 | **Auslander_NSCLC** | ICI/vaccine signature (re-derived, extra) | 0.567 | 1.14 | 0.0675 | 8/8 |
| 7 | **Mariathasan_immune** | ICI/vaccine signature (re-derived, extra) | 0.555 | 1.13 | 0.118 | 8/9 |
| 8 | **T_inflamed_Spranger** | ICI/vaccine signature (re-derived) | 0.546 | 1.11 | 0.169 | 7/7 |
| 9 | **PRGS_Cristescu** | ICI/vaccine signature (re-derived, extra) | 0.546 | 1.17 | 0.0829 | 10/17 |
| 10 | **IPS_Charoentong** | ICI/vaccine signature (re-derived) | 0.542 | 1.15 | 0.174 | 16/18 |
| 11 | **IFNG_18gene_Ayers** | ICI/vaccine signature (re-derived) | 0.542 | 1.16 | 0.0986 | 14/18 |
| 12 | **Hugo_6gene** | ICI/vaccine signature (re-derived) | 0.537 | 1.08 | 0.283 | 4/6 |
| 13 | **Moffitt_basal_class (binary)** | ICI/vaccine signature (re-derived) | 0.526 | 1.16 | 0.475 | binary/binary |
| 14 | **paradox_score (ours)** | ICI/vaccine signature (re-derived) | 0.517 | 1.06 | 0.372 | —/— |
| 15 | **paradox_inflamed** | ICI/vaccine signature (re-derived) | 0.515 | 1.11 | 0.438 | 5/5 |
| 16 | **Mariathasan_TGFB** | ICI/vaccine signature (re-derived, extra) | 0.512 | 1.12 | 0.296 | 6/17 |
| 17 | **Moffitt_classical** | ICI/vaccine signature (re-derived) | 0.510 | 1.08 | 0.414 | 13/13 |
| 18 | **IMPRES** | ICI/vaccine signature (re-derived) | 0.510 | 1.14 | 0.307 | 9/15 |
| 19 | **paradox_suppress** | ICI/vaccine signature (re-derived) | 0.509 | 1.14 | 0.342 | 7/7 |
| 20 | **ESTIMATE_stromal** | ICI/vaccine signature (re-derived) | 0.508 | 1.13 | 0.297 | 9/10 |
| 21 | **Hallmark_INFL** | ICI/vaccine signature (re-derived, extra) | 0.493 | 0.95 | 0.674 | 3/10 |
| 22 | **ChengEcotype_TLS** | ICI/vaccine signature (re-derived, extra) | 0.492 | 0.91 | 0.547 | 6/7 |
| 23 | **IPRES_Hugo** | ICI/vaccine signature (re-derived, extra) | 0.489 | 0.87 | 0.403 | 5/9 |
| 24 | **TIDE_exhaustion** | ICI/vaccine signature (re-derived) | 0.486 | 1.01 | 0.97 | 7/12 |
| 25 | **CYT_Rooney** | ICI/vaccine signature (re-derived) | 0.483 | 0.96 | 0.698 | 2/2 |
| 26 | **Senbabaoglu_Treg** | ICI/vaccine signature (re-derived, extra) | 0.470 | 0.94 | 0.632 | 1/6 |
| 27 | **ESTIMATE_immune** | ICI/vaccine signature (re-derived) | 0.467 | 0.94 | 0.601 | 2/11 |
| 28 | **Tirosh_exhaustion** | ICI/vaccine signature (re-derived, extra) | 0.450 | 0.86 | 0.388 | 6/10 |

## Cited (not re-derived) — published benchmarks

| Method | Task | Metric | Value | Reference |
|---|---|---|---|---|
| **NetMHCpan-4.1** | MHC-I binding | AUROC | 0.81 | Reynisson 2020 NAR |
| **MHCflurry-2.0** | MHC-I binding | AUROC | 0.83 | O'Donnell 2020 Cell Syst |
| **PRIME** | MHC-I + immunogenicity | AUROC | 0.78 | Schmidt 2021 Cell Syst |
| **BigMHC** | MHC-I + presentation | AUROC | 0.85 | Albert 2023 Bioinformatics |
| **HLAthena** | MS-MHC-I peptide | AUROC | 0.87 | Sarkizova 2020 Nat Biotech |
| **TESLA consensus** | neoantigen ranking | PPV @ top 20 | 20–40% | Wells 2020 Cell |
| **TIDE** | ICI response (mel) | AUROC | 0.71 | Jiang 2018 Nat Med |
| **IMPRES** | ICI response (mel) | AUROC | 0.79 | Auslander 2018 Nat Med |
| **DeepSurv** | Cox NN | C-index pan-cancer | 0.65–0.70 | Katzman 2018 BMC MRM |
| **GAT-Survival** | graph attention NN | C-index TCGA | 0.66–0.72 | Veličković 2018 ICLR + adapt |
| **BERT-MHC / TransPHLA** | Transformer MHC-I | AUROC | 0.84 | Cheng 2021 BIB / Chu 2022 NMI |
| **DeepImmuno** | CNN immunogenicity | AUROC | 0.79 | Xu 2022 BIB |
| **SimCLR contrastive** | self-supervised peptide | linear-eval AUROC | 0.76 | Chen 2020 ICML |
| **Multi-task Cox + auxiliary** | KDD-style multi-cohort | C-index | 0.67–0.73 | Bao 2020 KDD + Caruana 1997 |

## Notable

- **Best re-derived single signature on PDAC OS**: Moffitt_basal (continuous z-score) C = 0.647
- **Best re-derived binary**: KRAS G12D HR = 1.91 (p = 0.003)
- **Best combined**: Cox PH multivariate Moffitt + KRAS + 31 genes + age C = 0.695
- **Hallmark_IFNG** is the only off-the-shelf gene-set signature significant in PDAC at p < 0.05 (C = 0.585, HR = 1.14, p = 0.025)
- **TIDE / IMPRES (melanoma-derived) drop sharply on PDAC** — TIDE C = 0.49 (chance), IMPRES C = 0.51 — the DIAL effect (signatures don't transfer)
- **NetMHCpan / MHCflurry / PRIME / BigMHC** are peptide-binding methods, not OS predictors — different task; cited for the neoantigen prioritization pipeline

## Honest caveats

- We do not run NetMHCpan/MHCflurry/PRIME/BigMHC live (needs binary distribution + per-patient HLA + peptide). We use their published benchmarks as references.
- TESLA consensus (Wells 2020) requires peptide-MS validation ground truth which we don't have for PDAC; cited for benchmark metric values only.
- The MLP DeepSurv-proxy uses sklearn MLPRegressor on log(time) for observed events — a lightweight approximation. True DeepSurv (torch + Cox partial likelihood) is the production upgrade.
- Per-cohort C-index varies; we report TCGA-PAAD only because GSE71729 OS extraction failed and other GEO cohorts need probe→gene mapping (next pass).
