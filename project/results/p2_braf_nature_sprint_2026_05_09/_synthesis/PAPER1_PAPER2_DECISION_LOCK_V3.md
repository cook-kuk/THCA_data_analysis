# Paper 1/2 Decision Lock v3

Build: 2026-05-09 KST. This is factual scaffolding only. It is not
manuscript prose and does not edit Hook, Aim, Discussion 3.1, Limitations,
Cover Letter paragraph 1, or Q9.

## Executive Disposition

| Item | Verdict | Why | Safe placement |
|---|---|---|---|
| Paper 1 BRAF clinical axis | Lead / absorb | DM2 BRAF-like PFI HR 4.67; BRAF-cPTC DM2 PFI HR 5.68; DM1 x TERT pooled OS HR 19.2 | Paper 1 figure scaffold + supplement + reviewer Q factual blocks |
| Paper 2 image-DM1 | Pilot survives | UNI LOTO overall AUC 0.852 and BRAF-like AUC 0.855, but RAS-like is 0.718 and TCGA TSS risk remains | Paper 2 pilot / supplement / future validation |
| ViT-L RAS-like AUC=1.000 | Reject as headline | TSS-only RAS-like LR AUC 1.000; ViT-L LOTO RAS-like AUC 0.308 | Methods/audit cautionary example |
| K2 full-transcriptome re-quant | Negative-control reserve | 9/9 local FASTQs complete, full panel coverage, but all DM2 and no significant tumor-normal panel signal | Reviewer reserve only |
| H15v2 foundation rescue | Defer | No usable output; current image conclusion already bounded | Future work only |
| RunPod state | Stopped | All pods EXITED after result recovery | Operational note |

## Claim Tiers

### Tier 1: Safe Main-Analysis Claims

| Claim | Exact factual support | Source |
|---|---|---|
| DM2 is the BRAF-context adverse arm | TCGA BRAF-like DM2-vs-not-DM PFI HR 4.67 [1.41, 15.5], p=0.012; BRAF-cPTC DM2 PFI HR 5.68 | `r1_meta_survival/R1_REPORT.md`, `h12_rai_response/H12_REPORT.md` |
| DM1 x TERT is the replicated aggressive escape | TCGA+Landa OS REML HR 19.2 [7.16, 51.3], p=4.2e-9, I2=26%; Landa-only HR 12.78 and 14/14 deaths | `r1_meta_survival/R1_REPORT.md`, `h8_tert_replication/H8_REPORT.md` |
| Lee 2024 is expression-side Korean replication only | n=632 RNA with HT/FA/MAPK coverage; no public per-sample survival, driver, or TERT | `h30_lee2024_full/H30_REPORT.md`, `r4_lee2024_survival/R4_REPORT.md` |
| Immune axis is not just purity | HT-13 d remains +0.88 after leukocyte residualization | `r6_purity_audit/R6_REPORT.md` |
| Compartment model is supported | HT-13 peaks in myeloid/B/T; FA-12 in malignant cells; MAPK-9 in BRAF malignant cells | `r7_singlecell/R7v2_REPORT.md` |

### Tier 2: Safe Supplement / Reviewer-Reserve Claims

| Claim | Exact factual support | Source |
|---|---|---|
| UNI image signal survives center holdout overall | UNI LOTO overall AUC 0.852, BRAF-like AUC 0.855, RAS-like AUC 0.718 | `analysis_supp/audit_uni_loto/UNI_LOTO_SUMMARY.json` |
| K2 local re-quant was technically successful | 9/9 local FASTQs quantified, 61,228 gene symbols, HT-13/FA-12/MAPK-9 coverage 100% | `h_k2_full_quant/k2_full_summary.json` |
| K2 local re-quant does not validate DM1 vs DM2 | 9/9 samples inherit DM2 calls; HT-13 tumor-normal AUC 0.40, MW p=0.730; FA-12/MAPK-9 AUC 0.75, p=0.286 | `h_k2_full_quant/K2_PILOT_VERDICT.md` |
| Strict ESMFold outputs were recovered before pod stop | 89 rows in `strict_esmfold_features.tsv`; pod stopped after recovery | `p_neo_bayesian_2026_05_09/wave8b_strict_structure_2026_05_09/` |

### Tier 3: Blocked / Do-Not-Claim

| Blocked claim | Reason |
|---|---|
| K2 full-cohort Korean transcriptome replication | Only 9 local FASTQs are available out of 262 metadata runs. |
| K2 DM1-vs-DM2 full-transcriptome validation | All 9 quantified samples inherit DM2 calls. |
| Lee 2024 survival / driver / TERT replication | Public data do not expose those per-sample fields. |
| TCGA-only H&E image-DM1 main claim | TSS confounding and subgroup sparsity remain attackable. |
| H15v2 result improves the story | No usable local H15v2 result was present. |
| Causal methylation mechanism | Q9 boundary remains correlative: motivates rather than confirms. |

## Figure Allocation Queue

| Slot | Proposed factual content | Inputs | Status |
|---|---|---|---|
| Paper 1 main / clinical panel | BRAF-context DM2 adverse arm plus DM1 x TERT escape | R1, H8, H12, H24 | Ready as figure scaffold |
| Paper 1 biology panel | HT immune x FA thyroid-lipid compartment split | H2, H5, H10, R7 | Ready as figure scaffold |
| Paper 1 supplement | Purity and panel-disentanglement defense | R3, R6 | Ready |
| Paper 1 supplement / reviewer reserve | K2 n=9 technical pilot, non-deployable | `h_k2_full_quant/` | Ready but reserve only |
| Paper 2 supplement | UNI LOTO and TSS audit | `audit_uni_loto/`, `audit_ras_auc100/` | Ready |
| Paper 2 future-work gate | K2 H&E / external WSI validation | K2 H&E outreach, external slides | Blocked |

## Reviewer Q Mapping

| Reviewer risk | Factual answer | Artifact |
|---|---|---|
| "Is DM2 risk just stage/TSS?" | H24 sensitivity: 9/9 deployable subgroups HR>1; TSS-cluster robust HR 5.91 [2.95, 11.8] | `h24_survival_sensitivity/H24_REPORT.md` |
| "Is HT-13 just leukocyte contamination?" | Residual d remains +0.88 after leukocyte residualization | `r6_purity_audit/R6_REPORT.md` |
| "Is image-DM1 robust?" | ViT-L no; UNI pilot yes overall but not subgroup-main. Report as pilot. | `paper2_image_dm1_audit_dossier.html` |
| "Does K2 validate this?" | No. Completed local n=9 pilot is all-DM2 and non-deployable. | `h_k2_full_quant/K2_PILOT_VERDICT.md` |
| "Is this causal methylation?" | No. Keep Q9 boundary: correlative mechanism, not causal proof. | `r9_hm450_immune/R9_REPORT.md`, `h13_mediation/H13_REPORT.md` |

## Weakness Mitigation Lock

| Weakness | Reviewer-safe answer | Promotion status |
|---|---|---|
| DM1 is not generically adverse in BRAF/cPTC | Split the clinical story: DM2 is the BRAF-context adverse/RAI-refractory arm; DM1 x TERT is the replicated aggressive escape niche. | Promoted only as a two-gate claim |
| Lee 2024 lacks public survival / driver / TERT | Use Lee for expression-side Korean replication and HT/FA/MAPK coverage only. | Keep out of survival/driver claims |
| K2 local re-quant is n=9 and all DM2 | Treat as completed technical negative-control, not Korean DM1-vs-DM2 validation. | Reviewer reserve only |
| UNI image signal survives overall but not enough for subgroup deployment | Use UNI LOTO as cautious Paper 2 pilot support: overall AUC 0.852, BRAF-like 0.855, RAS-like 0.718. | Pilot / supplement |
| ViT-L RAS-like AUC=1.000 was TSS-driven | Retain as audit evidence showing why the image paper is bounded. | Cautionary methods example |
| HM450/RNA mechanism is correlative | Preserve the protected Q9 "motivates rather than confirms" boundary. | Supportive mechanism layer, not causal proof |
| MAPK forest has high heterogeneity | Explain using pre-specified two-route model: TCGA/Lee MAPK route, GSE286332 HT route, GSE76039 saturation. | Q14 / supplement |
| H15v2 has no usable rescue output | Do not restart GPU rescue without new WSI/K2 H&E or thyroid-pretrained model rationale. | Defer |

## Operational Closeout

| Resource | Final state |
|---|---|
| `thca-neo-bayesian-aux` | EXITED |
| `thca-img-dm1-a6000` | EXITED |
| `thca-spark-dm-a6000-v4` | EXITED after strict ESMFold output recovery |
| K2 kallisto runner | Completed 9/9; aggregate exit=0 |
| Live dossier | `project/papers_hub_2026_05_04/paper1_paper2_braf_axis_takeover_dossier.html` |

## Next Safe Actions

1. Add factual cross-references to non-voice-protected figure/caption scaffolds only.
2. Keep K2 and H15v2 out of main claims.
3. Use UNI LOTO as Paper 2 pilot support, with TSS caveat adjacent to the result.
4. Do not generate protected manuscript prose from this file.
