# Weakness Reinforcement Lock v1

Build: 2026-05-09 KST. Factual scaffold only. This file does not generate or
edit Hook, Aim, Discussion 3.1, Limitations 3.4, Cover Letter paragraph 1, or
Q9 prose.

## Executive Lock

The strongest manuscript posture is not to hide weak results. The strongest
posture is to assign every weak result to a narrow, traceable role and prevent
it from contaminating the main claim.

| Reviewer attack | Current weak point | Reinforcement move | Safe manuscript role | Source |
|---|---|---|---|---|
| "This is just a broad DM1-is-bad story." | DM1 is not generically adverse in BRAF/cPTC; DM2 is the adverse arm except for DM1 x TERT. | Split the clinical claim into two gates: DM2 = BRAF-context risk; DM1 x TERT = aggressive escape niche. | Main clinical scaffold / figure caption. | `r1_meta_survival/R1_REPORT.md`, `h8_tert_replication/H8_REPORT.md`, `h12_rai_response/H12_REPORT.md` |
| "Korean validation is overclaimed." | Lee lacks public survival, driver, and TERT fields; K2 full-transcriptome pilot is only n=9 and all DM2. | Use Lee for expression-side replication only; use K2 as completed negative-control reserve only. | Reviewer reserve / supplement caveat. | `h30_lee2024_full/H30_REPORT.md`, `r4_lee2024_survival/R4_REPORT.md`, `h_k2_full_quant/K2_PILOT_VERDICT.md` |
| "Image-DM1 is TSS/site artifact." | ViT-L RAS-like AUC=1.000 collapses under LOTO; RAS-like UNI LOTO is 0.718. | Keep UNI as cautious pilot because overall AUC=0.852 and BRAF-like AUC=0.855 survive center holdout; do not make subgroup-main image claims. | Paper 2 pilot / supplement. | `analysis_supp/audit_uni_loto/UNI_LOTO_SUMMARY.json`, `paper2_image_dm1_audit_dossier.html` |
| "Mechanism is not causal." | HM450/RNA/cohort evidence is correlative; no CRISPRi, ChIP-seq, DNMT perturbation, or PDX rescue. | Preserve Q9 boundary: methylation and MAPK/HT routes motivate, not confirm, an upstream silencing mechanism. | Q9 protected author voice; figure caveat only. | `r9_hm450_immune/R9_REPORT.md`, `h13_mediation/H13_REPORT.md`, `p_deconv_2026_05_08/v13_*.tsv`, `v14_cross_cohort_forest.tsv` |
| "Panel is cherry-picked." | Eight genes selected from a known thyroid differentiation program. | Point to TDS-16 lock: Panel-8 behaves like canonical TDS-16 across MAPK correlation, driver-class d, sub-A/sub-B convergence, and MAPK-high AUC. | Q13 / Supplementary Fig SX_v13. | `project/results/p_deconv_2026_05_08/v13_*.tsv` |
| "Cross-cohort heterogeneity weakens the MAPK story." | I2 is high in the five-cohort forest. | Treat heterogeneity as predicted by the two-route model: TCGA/Lee MAPK route, GSE286332 HT route, GSE76039 saturation. | Q14 / Supplementary Fig SX_v14. | `project/results/p_deconv_2026_05_08/v14_cross_cohort_forest.tsv` |
| "Purity explains immune axes." | HT-13 includes immune-lineage biology. | Do not call it contamination; call it a tissue-state component. Residual d remains +0.88 after leukocyte adjustment. | Main biology scaffold + supplement. | `r6_purity_audit/R6_REPORT.md`, `r7_singlecell/R7v2_REPORT.md` |
| "GPU rescue failed." | H15v2 has no usable output and should not be restarted casually. | Make the failed/absent rescue an operational boundary: all pods stopped, no more GPU spend unless K2 H&E or external WSI becomes available. | Operations note only. | `PAPER1_PAPER2_DECISION_LOCK_V3.md` |

## Claim Wording Guardrails

| Unsafe wording | Replace with factual wording |
|---|---|
| "K2 validates the Korean transcriptome result." | "The local K2 re-quant is a completed n=9 negative-control pilot and is not used as Korean DM1/DM2 validation." |
| "Lee validates survival." | "Lee/GSE213647 supports expression-side replication; public fields do not permit per-sample survival, driver, or TERT validation." |
| "Image model reads DM1 robustly." | "UNI image-DM1 signal survives center holdout overall, but subgroup and TSS limitations keep it in pilot/supplementary status." |
| "DM1 is the aggressive subtype." | "In BRAF-context disease, DM2 is the adverse arm; DM1 becomes high-risk chiefly through the replicated DM1 x TERT escape niche." |
| "Methylation explains the mechanism." | "Promoter methylation is a strong correlative layer that motivates, rather than confirms, upstream silencing." |
| "The 8-gene panel is the best subset." | "The 8-gene panel is a compact deployable readout of the canonical thyroid differentiation program; TDS-16 controls show no material gain." |

## Placement Map

| Artifact | Add / preserve | Do not add |
|---|---|---|
| `05_figure_captions.md` | Factual cross-references to K2 negative-control reserve, UNI LOTO pilot, TDS-16 lock, and cross-cohort forest. | No new causal mechanism interpretation beyond Q9 boundary. |
| `09_reviewer_qa.md` | Preserve existing Q1-Q14 structure unless the author explicitly opens new reviewer-Q numbers. Existing Q13/Q14 already cover panel and cross-cohort mechanism defenses. | No Q9 rewrite and no new Q-number expansion by default. |
| `paper1_paper2_braf_axis_takeover_dossier.html` | Dedicated weakness-lock section with reviewer attack -> factual defense -> disposition. | No claim promotion for H15v2 or K2. |
| Paper 1 main Results | DM2 BRAF-context risk and DM1 x TERT escape as separate claims. | Generic "DM1 adverse" phrasing. |
| Paper 2 | Pilot / supplement frame. | Standalone clinical deployment claim. |

## Final Disposition

Paper 1 is stronger if the weak edges are visible and bounded. The manuscript
should lead with the strong BRAF clinical axis, use mechanism data as correlative
support, keep Paper 2 as a pilot, and keep K2/H15v2 as reserve/negative-control
artifacts rather than rescue evidence.
