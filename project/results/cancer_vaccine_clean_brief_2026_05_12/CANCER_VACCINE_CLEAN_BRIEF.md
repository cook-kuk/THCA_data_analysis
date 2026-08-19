---
title: "Cancer vaccine clean brief"
date: 2026-05-12
status: working scaffold
purpose: clean executive brief for the post-TESLA neoantigen platform
---

# Cancer vaccine clean brief

## 1. One-line thesis

The platform is a research-grade neoantigen vaccine selection system: it ingests public neoantigen data, audits leakage, ranks candidates with simple and LM-based models, and turns the best hits into executable vaccine outputs.

## 2. Current scoreboard

| metric | value | meaning |
|---|---:|---|
| sources ingested | 13 | broad public benchmark base |
| cleaned records | 119,237 | master project scale |
| leakage-free benchmark | 19,817 | safe evaluation pool |
| 5-fold RF AUROC | 0.854 +/- 0.008 | strongest simple baseline |
| 7-peptide coverage | 66.3% | population set-cover candidate |
| top KRAS PAAD candidate | GADGVGKSAL | proba 0.871 |
| ESM2 LOSO rescue | +0.073 on NEPdb-out | partial domain-shift fix |

## 3. What matters now

| lane | current status | why it matters |
|---|---|---|
| Selection | public-data ranking + leakage audit | makes the shortlist defensible |
| Structure | 115 ESMFold + 6 real structures | prevents score-only overclaiming |
| Translation | 7-peptide set cover + KRAS lead | gives an executable vaccine output |
| Boundary | no prospective locked batch yet | keeps the claim research-grade |

## 4. Pitch deck flow

1. Problem: neoantigen selection is noisy and easy to overclaim.
2. Proof: we audited the data instead of hiding leakage.
3. Model: RF anchors the ranking, ESM2 rescues some LOSO rows.
4. Impact: the output is a real shortlist, not a theory slide.
5. Ask: a locked prospective batch turns this into validation.

## 5. What the page should say first

| lane | best current result | note |
|---|---:|---|
| biophysics baseline | RF AUROC 0.854 +/- 0.008 | best simple model |
| LM augmentation | ESM2 improves some LOSO rows | not universal |
| delivery | 7-peptide set cover at 66.3% | clear executable output |
| structure | 115 ESMFold + 6 real pMHC/TCR-pMHC | visual and structural support |
| top KRAS lead | GADGVGKSAL, proba 0.871 | PAAD / actionable shortlist |

## 6. What the platform already does

1. Ranks neoantigen candidates from public data.
2. Flags training overlap instead of hiding it.
3. Keeps TCR and structure layers attached to the ranking stack.
4. Produces a compact vaccine deliverable set, not just scores.

## 7. What is still blocked

| blocked item | why |
|---|---|
| universal predictor claim | dataset and assay shift remain real |
| wetlab superiority claim | no prospective locked validation yet |
| one-number SOTA claim | overlap and LOSO failure make it unsafe |
| clinical deployment claim | this is still a research platform |

## 8. Evidence stack that should survive review

| layer | evidence | short read |
|---|---|---|
| corpus | 13 sources, 119,237 cleaned records | public benchmark base |
| leakage control | 19,817 leakage-free benchmark | overlap-filtered pool |
| ranking | RF 0.854 +/- 0.008, ESM2 rescue +0.073 | strongest current models |
| structure | 115 ESMFold + 6 real structures | structural anchoring |
| deliverable | 7-peptide set cover, 66.3% coverage | executable shortlist |
| lead | GADGVGKSAL at 0.871 | most concrete KRAS item |

## 9. Important figures

| figure | what it shows | why it matters |
|---|---|---|
| source/leakage figure | 13 sources and the overlap audit | proves the shortlist is not inflated by leakage |
| model leaderboard figure | RF vs ESM2 vs XGBoost and public comparators | shows the ranking claim is tested, not assumed |
| set-cover figure | 7-peptide world coverage | turns score into a concrete deliverable |
| KRAS figure | GADGVGKSAL at 0.871 | makes the lead actionable for PAAD follow-up |
| structure figure | 115 ESMFold + 6 real structures | anchors the shortlist in biology, not score only |
| LOSO figure | NEPdb-out +0.073 rescue | shows the LM adds value under domain shift |

All main figures are click-to-enlarge in the HTML brief, so the advisor can inspect the plots directly instead of reading them as thumbnails.

## 10. TESLA benchmark in one slide

| point | summary |
|---|---|
| What TESLA is | a canonical public neoantigen benchmark used to compare selection tools |
| Why it matters | it exposes generalization gaps and lets us compare against prior work |
| Scale | 5 subjects, 535 validated peptides, 34 positives in the benchmark description used by validation papers |
| How we use it | as a stress-test source and benchmark reference, not as the only evidence |
| Why it is not enough | it is small and benchmark-heavy, so prospective batch data is still needed |

**TESLA visual takeaway**

- Neo-intline was reported above older baselines on TESLA in the validation-paper comparison.
- The key point is not the absolute number; it is that the benchmark can separate methods.
- The page therefore treats TESLA as a filled reference slide, not as the whole claim.

## 11. Comparison figure

| panel | what it says | why it helps us |
|---|---|---|
| front-door stack | BioDarwin + NeoQ sits at the top | combines ranking, audit, structure, and output in one view |
| RF anchor | strongest simple baseline | honest floor, but not the whole package |
| ESM2 rescue | source-shift helper | useful branch, not the main story |
| public comparators | MHCflurry / PRIME / DeepImmuno | context only, overlap-sensitive |

| system | ranking | leakage audit | structure | deliverable | best use |
|---|---|---|---|---|---|
| **BioDarwin + NeoQ** | high | explicit | yes | 7-peptide cover | front door |
| RF biophysics | high floor | yes | no | no | baseline |
| ESM2 + LogReg | shift rescue | yes | no | no | fallback branch |
| XGBoost | challenger | yes | no | no | ceiling check |
| MHCflurry / PRIME / DeepImmuno | context only | caveated | no | no | literature context |
| BAR-Neo / KG_GA / CROSS-Neo | other boards | board-specific | some | some | cross-page winners |

The comparison figure is claim-weighted, not raw-AUROC weighted. It is meant to show why our stack should occupy the front door of the page.

## 12. Algorithm comparison

| algorithm / comparator | what it compared | current read |
|---|---|---|
| RF biophysics | 5-fold CV on leakage-aware benchmark | strongest simple baseline and honest floor |
| ESM2-35M + LogReg | sequence-LM feature stack | rescue branch for some LOSO / source-shift rows |
| XGBoost | nonlinear tabular learner | tests whether a booster beats the RF ceiling |
| NetMHCpan / MHCflurry / PRIME / DeepImmuno / TransPHLA | public benchmark landscape | caveated comparators, not universal SOTA |
| BAR-Neo / KG_GA / CROSS-Neo | other repo pages / same-board champions | different board, different winner logic |
| structural layer | ESMFold + real pMHC/TCR-pMHC | supports biological plausibility |
| BioDarwin + NeoQ | ranked shortlist with audit and output | front-door stack with the widest claim surface |

**Comparison framing**

- RF is the anchor card.
- ESM2 is the rescue card.
- XGBoost is the challenger card.
- Public comparators remain comparators, not co-winners.
- Cross-page champions can win on their own boards; this page should not flatten them into one number.
- BioDarwin + NeoQ is the only row that simultaneously combines ranking, audit, structure, and a concrete translational output.

## 13. Why this is high impact

| reason | short explanation |
|---|---|
| scale | 13 sources and 119,237 cleaned records is not a toy dataset |
| honesty | leakage is audited instead of hidden |
| utility | the output is a real shortlist, not only a score |
| structure | ranking is tied to structural evidence |
| boundary | the page is explicit about what still needs prospective validation |

## 14. Clean paper shape

| section | job |
|---|---|
| Intro | why neoantigen selection is hard |
| Data | what sources were ingested and filtered |
| Ranking | what model stack performs best |
| Deliverables | what concrete vaccine candidates come out |
| Caveats | what fails and why |
| Methods | how to reproduce the pipeline |

## 15. Executive decision matrix

| decision | recommendation |
|---|---|
| main framing | neoantigen selection platform |
| strongest simple baseline | RF biophysics model |
| strongest practical output | 7-peptide set cover + KRAS lead |
| best readability | clean brief + separate full dossier |
| next unlock | new locked prospective batch |

## 16. What to ask next

| ask | needed data | why |
|---|---|---|
| prospective locked batch | fresh peptides + hit/fail labels | turns the project into a validation story |
| FFPE / institutional validation | independent cohort with same pipeline | checks external reproducibility |
| HLA / patient metadata | typing, peptide provenance, assay context | reduces ambiguity in shortlist ranking |
| wetlab readout schema | hit / fail / helper / control | lets the interpreter update without redesign |

## 17. Reproducibility anchors

- Project root: `/data/neoantigen_vaccine_hub`
- Source record: `project/papers_hub_2026_05_04/CANCER_VACCINE_PROJECT_FULL.md`
- Detailed dossier: `project/papers_hub_2026_05_04/cancer_vaccine_full_dossier.html`
- Agent page: `project/papers_hub_2026_05_04/cancer_vaccine_agent.html`
