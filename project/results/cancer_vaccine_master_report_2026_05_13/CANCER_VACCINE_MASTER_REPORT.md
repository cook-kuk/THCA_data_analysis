---
title: "Cancer vaccine master report"
date: 2026-05-13
status: working scaffold
purpose: long-form customer-facing synthesis of the cancer vaccine / neoantigen / BioDarwin stack
---

# Cancer vaccine master report

## 1. Executive thesis

This project is not a single model. It is a **decision system** for neoantigen vaccine selection that combines:

1. a literature map,
2. a cleaned public corpus,
3. leakage-aware ranking,
4. structural support,
5. a router / ensemble layer,
6. a 96-well validation scaffold,
7. and a customer-facing deliverable set.

The useful claim is therefore narrower and stronger than "one model wins everywhere":

- the platform is research-grade,
- the shortlist is executable,
- the comparison is honest,
- the prospective path is explicit.

## 2. What this report is for

This report is the long-form version for a first-time customer MD.

It answers:

1. what data were used,
2. how the ranking logic evolved,
3. what figures matter,
4. which algorithms were compared,
5. why the project is high impact,
6. what is still blocked,
7. and which page to open next.

## 3. Fast scoreboard

| metric | value | meaning |
|---|---:|---|
| sources ingested | 13 | broad public benchmark base |
| cleaned records | 119,237 | master corpus scale |
| leakage-free benchmark | 19,817 | safe evaluation pool |
| RF 5-fold AUROC | 0.854 +/- 0.008 | strongest simple baseline |
| 7-peptide coverage | 66.3% | executable set-cover candidate |
| top KRAS candidate | GADGVGKSAL | practical PAAD lead |
| ESM2 LOSO rescue | +0.073 on NEPdb-out | partial domain-shift fix |
| public-router v1 -> v4 | 0.7104 -> 0.9205 mean AUPRC | regime-aware routing gain |
| industrial locked v0 | 0.6298 AUPRC / 0.6458 AUROC | diagnostic stress test |
| 96-well locked control | 0.6430 AUPRC / 0.6806 AUROC | locked TG4050 control |

## 4. Presentation flow

This is the slide order that works best for a customer MD.

| slide | headline | speaker line | why it lands |
|---|---|---|---|
| 1 | Problem and data spine | Neoantigen selection is noisy and easy to overclaim. | Starts with the right tension |
| 2 | What we used | We audited the data instead of hiding the mess. | Gives the page credibility |
| 3 | What we found | RF anchors the ranking, ESM2 rescues some shift rows. | Makes the output concrete |
| 4 | Why the stack wins | BioDarwin + NeoQ combines ranking, audit, structure, and deliverable. | Explains the comparison |
| 5 | How we validate | The 96-well plate turns ranking into a live readout schema. | Bridges to wetlab |
| 6 | What we need | Give us a locked prospective batch. | Ends with a clear ask |

## 5. How to read the project

Read it in this order:

1. clean brief,
2. this master report,
3. timeline index,
4. full dossier,
5. BioDarwin router dossier,
6. 96-well interpreter,
7. external comparison,
8. partner pack.

That order prevents the common mistake of reading a benchmark table before the data spine.

## 6. Chronology

| date | phase | key artifact | what changed |
|---|---|---|---|
| 2026-05-07 | literature map | `post_tesla_key_papers.md` | shifted the field view from binding-only prediction to presentation / recognition / clinical-response logic |
| 2026-05-08 | master corpus | `cancer_vaccine_full_dossier.html` | consolidated the data spine, leakage control, and structural support |
| 2026-05-09 | execution surface | `cancer_vaccine_agent.html` | moved from a report to a working tool-like surface |
| 2026-05-10 | search build | `biodarwin_pan_vaccine_ga_rl.html` | started GA/RL search and RF-anchor expansion |
| 2026-05-11 | routing layer | `biodarwin_public_router_impact_dossier.html` | moved lookup -> GA/RL -> LOBO blend -> contextual MoE |
| 2026-05-11 | validation layer | `biodarwin_96well_interpreter_v1.html` | added hit/fail/helper/probe logic for 96-well reading |
| 2026-05-11 | external comparison | `biodarwin_external_mega_comparison_2026_05_11.html` | added public winners and reviewer-safe boundary handling |
| 2026-05-12 | pitch front door | `cancer_vaccine_clean_brief_2026_05_12.html` | turned the work into a deck-like customer page |
| 2026-05-13 | chronology layer | `cancer_vaccine_timeline_2026_05_13.html` | taught the order of the stack |
| 2026-05-13 | long synthesis | this page | combines the chronology, figures, algorithms, and claim boundary |

## 7. Figure wall

| figure | what it says | why it matters |
|---|---|---|
| Source / leakage wall | what came in, what was removed, and why the benchmark is honest | proves the shortlist is not inflated |
| NeoQ decision rule | why the shortlist is rule-based rather than highest-score-wins | makes the logic legible |
| Algorithm winner map | floor, rescue, challenger, and front-door stack | shows the role of each comparator |
| TESLA reference panel | benchmark context, not whole claim | keeps the field reference clear |
| 96-well interpretation | retrospective ranking to prospective lab decision surface | connects to validation |
| Ask / boundary panel | next step and claim boundary | ends with an action |

## 8. What data were used

| layer | used data or evidence | why it matters |
|---|---|---|
| public corpus | 13 sources, 119,237 cleaned records | gives the project scale and traceability |
| safe evaluation | 19,817 leakage-free benchmark rows | keeps the score honest |
| benchmark context | TESLA, CEDAR, NEPdb, ITSNdb, PRIME, MHCflurry, DeepImmuno, NetMHCpan, TransPHLA | shows the external landscape and why winners differ |
| structural support | 115 ESMFold structures + real pMHC / TCR-pMHC examples | anchors the shortlist in biology |
| deliverable layer | 7-peptide set cover, KRAS lead `GADGVGKSAL` | makes the output executable |
| customer surface | clean brief, timeline, partner pack, 96-well interpreter | makes the stack readable outside the lab |

## 9. Important figures

| figure family | what it shows | why it matters |
|---|---|---|
| source / leakage figure | 13 sources and overlap audit | proves the shortlist is not inflated by leakage |
| model leaderboard figure | RF vs ESM2 vs XGBoost and public comparators | tests the ranking claim instead of assuming it |
| set-cover figure | 7-peptide world coverage | turns scores into a concrete deliverable |
| KRAS figure | `GADGVGKSAL` lead | makes the output actionable |
| structure figure | 115 ESMFold + real pMHC / TCR-pMHC examples | supports biological plausibility |
| router figure | lookup -> GA/RL -> LOBO -> MoE | explains why the system improved over time |
| 96-well figure | hit / fail / helper / probe lanes | connects retrospective ranking to prospective validation |
| external comparison figure | public comparator landscape | prevents a narrow or inflated claim |

## 10. Algorithm comparison

| algorithm / comparator | role | what it tested |
|---|---|---|
| RF biophysics | honest floor | whether simple biophysical features already solve the task |
| ESM2 + LogReg | rescue branch | whether protein-language features help under domain shift |
| XGBoost | challenger | whether a nonlinear tabular model beats the anchor |
| NetMHCpan / MHCflurry / PRIME / DeepImmuno / TransPHLA | public comparators | where the project sits in the broader literature |
| BAR-Neo / KG_GA / CROSS-Neo | same-board champions | how strong board-specific systems behave |
| BioDarwin + NeoQ | front-door stack | whether ranking, audit, structure, and deliverable can be unified |

The useful conclusion is not "one model wins everything."  
The useful conclusion is:

- RF is the anchor,
- ESM2 is the rescue branch,
- XGBoost is the challenger,
- public methods are the context,
- same-board champions prove the platform idea,
- BioDarwin + NeoQ is the only stack that combines ranking, audit, structure, and output in one front door.

## 11. BioDarwin evolution

| stage | old state | new state | why it matters |
|---|---|---|---|
| v1 | lookup policy | frozen public pool baseline | shows a naive policy is not enough |
| v2 | GA/RL | feature-gated MoE | shows the router can learn |
| v3 | LOBO blend | freeze point | gives a claim-safe operating line |
| v4 | contextual MoE | operating champion | shows same-pool saturation explicitly |
| locked TG4050 | diagnostic test | industrial stress test | shows the route still beats a strong comparator |
| 96-well control | score table | wetlab-facing interpreter | makes the system prospective-ready |

Current read:

- the public-pool search is strong,
- the same-pool ceiling is reached,
- the next real gain needs either a new expert pool or a new batch.

## 12. External landscape

The external public landscape is not dominated by one universal winner.

| benchmark family | current winner / strong method | takeaway |
|---|---|---|
| CEDAR partial | MHCflurry | this dataset rewards a different inductive bias |
| ITSNdb | DeepImmuno | this dataset has its own winner logic |
| NEPdb | PRIME | public benchmark winners are dataset-specific |
| TESLA mmc4 | PRIME | low positive-rate setting changes the ceiling |
| TESLA mmc7 | PRIME / BigMHC_IM context | benchmark sensitivity remains high |

This is the reason the project should not claim universal SOTA.

Instead, the claim is:

- the project learns a router,
- the router chooses the right expert family,
- and the system is honest about when a public comparator is better.

## 13. 96-well decision layer

The 96-well plate is the prospective-ready bridge.

| lane | meaning | use |
|---|---|---|
| hit | high-confidence candidate | priority for wetlab follow-up |
| fail | low-confidence / low-consensus row | down-rank or exclude |
| helper | ambiguous but informative | useful for calibration |
| probe | exploratory row | used to test generalization |

The important point is that the plate is not just a list.
It is a decision surface.

## 14. What the project already does

1. Ingests public neoantigen data.
2. Audits leakage instead of hiding it.
3. Ranks candidates with both simple and LM-based models.
4. Keeps structure and TCR support attached to the ranking stack.
5. Produces executable deliverables instead of only scores.
6. Separates retrospective tables from prospective scaffolds.

## 15. What is still blocked

| blocked item | why |
|---|---|
| universal predictor claim | assay and distribution shift remain real |
| wetlab superiority claim | no new locked prospective batch yet |
| one-number SOTA claim | dataset-specific winners remain strong |
| clinical deployment claim | the work is still a research platform |

## 16. Why this is high impact

| reason | short explanation |
|---|---|
| scale | 13 sources and 119,237 cleaned records is a real corpus |
| honesty | leakage is audited instead of hidden |
| utility | the output is a shortlist, not just a score |
| structure | ranking is tied to structural evidence |
| delivery | the project is packaged for a customer MD |
| next step | a locked prospective batch has a clear path |

## 17. Figure / page map

| open order | page | what it does |
|---|---|---|
| 1 | `cancer_vaccine_clean_brief_2026_05_12.html` | pitch-deck front door |
| 2 | `cancer_vaccine_timeline_2026_05_13.html` | chronology and learning order |
| 3 | `cancer_vaccine_full_dossier.html` | deep evidence and methods |
| 4 | `biodarwin_public_router_impact_dossier.html` | router evolution and validation layer |
| 5 | `biodarwin_96well_interpreter_v1.html` | wetlab-ready scoring |
| 6 | `biodarwin_external_mega_comparison_2026_05_11.html` | comparator landscape |
| 7 | `cross_neo_partner_pack.html` | customer-side delivery logic |

## 18. What to ask next

| ask | needed data | why |
|---|---|---|
| prospective locked batch | fresh peptides and labels | turns the platform into a live validation story |
| independent cohort | same pipeline, new population | checks external reproducibility |
| HLA / patient metadata | provenance and typing | reduces ambiguity in ranking |
| wetlab readout schema | hit / fail / helper / probe | lets the interpreter update cleanly |

## 19. Bottom line

This work is strongest when it is presented as a chronology of increasingly strict claims:

1. literature map,
2. corpus build,
3. leakage-aware ranking,
4. structural support,
5. router evolution,
6. external comparator landscape,
7. wetlab scaffold,
8. customer-facing brief.

That is the story the project can defend now.
