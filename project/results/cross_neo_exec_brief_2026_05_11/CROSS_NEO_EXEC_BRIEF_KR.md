# CROSS-Neo Executive Brief

## 결론

- `KG_GA_evolved` is the local champion on the 2,715-candidate CROSS-Neo board: AUPRC 0.933, AUROC 0.943, top96 96/96.
- Frozen validation-like subset: AUPRC 0.978, AUROC 0.999, top24 11/24.
- Low/medium leakage subset: AUPRC 0.757, top96 36/96.
- Paperability reward: 0.922, above `CROSS_claimsafe` 0.750 and `CROSS_integrated` 0.703.

## Same-Board Comparator Ranking

| algorithm | all_AUPRC | all_AUROC | frozen_validation_AUPRC | low_medium_leakage_AUPRC | v6_smoke_AUPRC |
|---|---:|---:|---:|---:|---:|
| KG_GA_evolved | 0.933 | 0.943 | 0.978 | 0.757 | 0.812 |
| CROSS_integrated | 0.878 | 0.896 | 0.714 | 0.573 | 0.790 |
| CROSS_stress | 0.872 | 0.873 | 0.651 | 0.601 | 0.777 |
| CROSS_claimsafe | 0.861 | 0.864 | 0.846 | 0.617 | 0.788 |
| CROSS_BMA | 0.861 | 0.860 | 0.660 | 0.636 | 0.767 |
| CROSS_finetuned | 0.841 | 0.863 | 0.606 | 0.669 | 0.681 |
| BigMHC_IM | 0.672 | 0.722 | 0.334 | 0.374 | 0.672 |
| BigMHC_EL | 0.648 | 0.719 | 0.308 | 0.306 | 0.866 |

## Strict No-Overlap Reliability

| algorithm | neo_strict_AUPRC | neo_strict_AUROC |
|---|---:|---:|
| BAR-Neo_confidence | 0.857 | 0.980 |
| BAR-Neo | 0.825 | 0.979 |
| BAR-Neo_patient_gated | 0.820 | 0.978 |
| BAR-Neo-X | 0.800 | 0.966 |
| BAR-Neo-X_claim_safe | 0.800 | 0.966 |
| RF_biophys | 0.752 | 0.685 |
| MHCflurry_2.0_presentation | 0.726 | 0.621 |
| KG_GA_evolved_controller | 0.624 | 0.897 |

## Broader Public Context

- `RF_biophys` is the strongest local public-context baseline: best AUROC 0.973, mean AUROC 0.891.
- `PRIME`: best AUROC 0.805.
- `MHCflurry`: best AUROC 0.756.
- `BigMHC_IM`: best AUROC 0.748.
- `TransPHLA`: best AUROC 0.728.
- `NetMHCpan_4.1`: best AUROC 0.567.

## Latest Literature SOTA Watchlist

- `NetMHCpan-4.2`: official 2026 server update for MHC-I presentation / neoepitope mode, with transfer learning + structural features.
- `ImmugenX`: pMHC immunogenicity, AUROC 0.619, AP 0.514.
- `PepFore / NeoaPred`: structure/surface-aware foreignness model, AUROC 0.81, AUPRC 0.54; BigMHC 0.70/0.30 in the same paper.
- `NeoGuider`: advanced feature engineering, benchmarked on 7 cohorts / 113 patients / 635 immunogenic candidates.
- `TrambaHLApan`: Transformer + Mamba, presentation + immunogenicity joint predictor.
- `NeoTImmuML`: best model AUROC 0.8707 on independent test set.

## Claim Boundary

- Allowed: retrospective prioritization champion, blinded validation proposal, assay-control generation, selection-layer collaboration.
- Not allowed: prospective immunogenicity proof, clinical utility, or universal SOTA over every public method without reruns.
- Clean boundary: `BAR-Neo_confidence` remains the strict no-overlap reliability winner; `KG_GA_evolved` is the discovery/prioritization champion.

## BioDarwin / NeoICI bridge

- `BioDarwin_mode_bank_v4` is the industrial locked winner in the neoantigen lane: AUPRC 0.628, AUROC 0.639 on TG4050 locked v0.
- `NeoICI_GA_RL_combo_v1` is the combo-therapy readiness lane: AUPRC 0.863 on the local gauntlet, with open combo replay support.
- Open combo replay: 9 paired samples across GSE222011 / GSE255830, median response-like index 0.033.
- Proxy comparison: NeoICI-style proxy Spearman 0.817 vs NeoPrecis-style proxy 0.417 on the same replay; CTMS1_followup4 is the clearest case lane.
- Practical split: use `KG_GA_evolved` for same-board discovery / prioritization, `BAR-Neo_confidence` for strict no-overlap reliability, `BioDarwin` for neoantigen selection, and `NeoICI` for combo-readiness framing.

## BD Packaging

- Non-confidential teaser: `project/results/cross_neo_moderna_bd_package_2026_05_11/MODERNA_NONCONFIDENTIAL_TEASER.md`
- Outreach draft: `project/results/cross_neo_moderna_bd_package_2026_05_11/MODERNA_OUTREACH_EMAIL_AFTER_IP.md`
- Pre-NDA packet: `project/results/cross_neo_moderna_bd_package_2026_05_11/pre_nda_send_packet.zip`
- Contact identity: `kukshomr@snu.ac.kr`, `@Ho15421542`

## Official Contact Routing

- Merck BD&L: `https://www.merck.com/research/business-development-and-licensing/`
- Merck fallback contact: `https://www.merck.com/contact-us/`
- Moderna routing page: `https://www.modernatx.com/en-US/contact-moderna`
