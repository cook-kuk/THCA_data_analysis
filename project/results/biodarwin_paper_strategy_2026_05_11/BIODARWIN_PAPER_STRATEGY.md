---
title: "BioDarwin paper strategy dossier"
date: 2026-05-11
status: working scaffold
purpose: turn the current BioDarwin results into a stronger, reviewer-safe paper structure without touching voice-protected manuscript prose
---

# BioDarwin paper strategy dossier

## 1. Core thesis

The paper is strongest when framed as a **regime-aware routing paper**, not as a universal SOTA claim.

- What it does: route each sample to the expert family that best fits its regime.
- What it proves well: the router improves over lookup and over the best single expert on the frozen public bundles, and it wins the locked industrial mini-set in diagnostic mode.
- What it does not prove: universal dominance across all benchmark families, or prospective wetlab superiority without a new locked batch.

## 1.1 What makes the paper better

The paper gets stronger when it is split by job, not by dataset size.

- One paper should be the **platform/router paper**.
- A second paper can be the **application/validation paper** if new locked data arrive.
- Do not force a universal SOTA story into a single manuscript.

## 2. Current evidence ledger

| layer | result | interpretation |
|---|---:|---|
| public-router v1 lookup | mean AUPRC 0.7104 | baseline lookup policy |
| public-router v2 GA/RL | mean AUPRC 0.9186 | feature-gated MoE over public metadata |
| public-router v3 LOBO blend | mean AUPRC 0.9205 | frozen blend point |
| public-router v4 contextual MoE | mean AUPRC 0.9205 | same-pool ceiling reached |
| industrial locked v0 | AUPRC 0.629766 / AUROC 0.645833 | router beats BigMHC_IM, but diagnostic only |
| 96-well locked control | AUPRC 0.643014 / AUROC 0.680556 | public-anchor / mode-bank wins locked TG4050 control |
| 96-well GA/RL search | AUPRC 0.669359 / AUROC 0.701389 | best retrospective meta-score champion |

## 3. Claim ladder

1. **Supported.** BioDarwin improves over lookup on the frozen public pool.
2. **Supported.** BioDarwin mode-bank / public-anchor wins the 25-row locked TG4050 control.
3. **Supported.** GA/RL search can improve the 96-well meta-score over the manual router.
4. **Supported with caveat.** The current expert pool is saturated; v4 equals v3 blend.
5. **Blocked.** Universal all-dataset SOTA.
6. **Blocked.** Prospective wetlab superiority before a new locked batch.

## 4. Recommended paper shape

### Main paper theme

The best paper is a **platform / router paper**:

- define the routing problem,
- show that a regime-aware controller outperforms lookup,
- show that the controller remains useful on a locked industrial stress test,
- show the 96-well decision surface,
- state the claim boundary explicitly.

### Not recommended

- A single "we beat everything everywhere" claim.
- A paper that hides overlap-heavy and low-power datasets.
- A paper that treats dataset-specific winners as row-wise generalization.

## 5. Figure order

| figure | title | main message |
|---|---|---|
| Fig 1 | Problem, regimes, and benchmark map | why a router is needed at all |
| Fig 2 | Router progression v1 -> v4 | lookup -> GA/RL -> LOBO blend -> contextual MoE |
| Fig 3 | 96-well plate and interpreter | how hit/fail/helper/probe lanes are defined |
| Fig 4 | External family landscape | dataset-specific winners are real and must be routed |
| Fig 5 | Genome / ablation / feature weights | what the GA/RL actually learned |
| Fig 6 | Claim boundary and failure modes | what is supported, what is not, and why |
| Supplement | overlap audit, reviewer-safe comparator table, frozen trace | reviewer defense and reproducibility |

## 6. Ablation matrix

| arm | role | keep in main text? | why |
|---|---|---:|---|
| v1 lookup | baseline | yes | establishes the ceiling of a naive policy |
| v2 public GA/RL | search upgrade | yes | proves the controller can learn from metadata |
| v3 LOBO blend | freeze point | yes | best claim-safe current router |
| v4 contextual MoE | operating champion | yes | same-pool saturation shown explicitly |
| RF_biophys anchor | disagreement feature | supplement | useful anchor, not standalone claim model |
| public anchor | external family expert | yes | explains why public-family routing improves |
| mode-bank | internal-feature expert | yes | best internal-feature route |
| reviewer-safe baselines | comparator landscape | supplement | essential for fairness, not row-wise claim |

## 7. Reviewer attack map

| attack | best answer |
|---|---|
| "This is just tuned on TG4050." | the public family table shows dataset-specific winners, and the industrial lock is separated as diagnostic until frozen |
| "Overlap made the result easy." | overlap-heavy rows are called out explicitly; reviewer-safe tables are separated from performance-max tables |
| "Why not claim universal SOTA?" | because the task and objective differ across bundles; a router is the correct model class |
| "Where is prospective validation?" | the 96-well plate/interpreter is the prospective-ready scaffold; a new locked batch is the next unlock |
| "RF_biophys should be the main model." | RF_biophys is an anchor/disagreement feature, not the final model |

## 8. Decision matrix

| decision | recommendation | rationale |
|---|---|---|
| main paper framing | regime-aware router | strongest claim, cleanest evidence, reviewer-safe |
| current freeze point | v3 LOBO blend | simplest stable benchmark story |
| operating champion | v4 contextual MoE | best current public-pool performance |
| next performance jump | new expert pool or new data | same-pool search is saturated |
| comparator policy | dataset-specific winner tables | avoids false universal claim |

## 9. Paper split plan

| paper | purpose | keep / delay |
|---|---|---|
| Platform paper | show routing logic, regime-aware gains, and reviewer-safe comparator handling | keep now |
| Application paper | focus on a single locked industrial use case or wetlab decision lane | delay until new locked batch |
| Prospective paper | use the 96-well interpreter plus a new untouched batch | delay until the next freeze |

### Platform paper

- Main story: regime-aware routing beats lookup and stays honest across family-specific winners.
- Main claim: improved selection quality under explicit regime control.
- Main risk: reviewer asks for universal dominance. Answer: not the right claim class.

### Application paper

- Main story: one frozen application lane, one comparator board, one locked decision outcome.
- Main claim: a narrow deployment setting with a clear operating envelope.
- Main risk: too many benchmark families. Answer: keep the bench landscape in supplement only.

## 10. Journal ladder

| level | target | fit | current status |
|---|---|---|---|
| floor | npj / JCI Insight | conservative, claim-safe, strong validation story | feasible now |
| reach | Cell Reports Medicine / Nature Communications | needs locked prospective evidence or a cleaner application lane | not yet locked |
| stretch | Cell / Nature Medicine | needs a new prospective batch plus stronger cross-site proof | blocked |

## 11. Next unlocks

1. Add a new expert pool instead of further squeezing the same pool.
2. Re-run the 96-well search after the next locked batch.
3. Keep the reviewer-safe benchmark table separate from the performance-max table.
4. Use the claim-boundary table as the manuscript editing guardrail.
5. Decide whether the next manuscript is platform-only or platform + application split.

## 12. Manuscript section map

| section | core job | evidence to place there |
|---|---|---|
| Title / Abstract | say it is a router paper, not a universal-SOTA paper | v1->v4 progression, 5/5 bundle lift, industrial locked v0, 96-well scaffold |
| Introduction | define the routing problem and why a regime-aware controller is needed | dataset-specific winners, public-family heterogeneity, reviewer-safe comparator split |
| Results 1 | show the lookup baseline is not enough | v1 lookup vs v2 GA/RL |
| Results 2 | show the learning-based router improves the public pool | v2 GA/RL vs v3 LOBO blend vs v4 contextual MoE |
| Results 3 | show the industrial locked stress test | BigMHC_IM vs regime_router / mode-bank |
| Results 4 | show the 96-well interpretation layer | locked-control table, hit/fail/helper/probe split |
| Results 5 | show the external landscape is dataset-specific | MHCflurry / PRIME / DeepImmuno / RF_biophys / BigMHC_IM split winners |
| Discussion | define the claim boundary and why the paper is still strong | same-pool saturation, new pool needed, prospective validation left open |
| Methods | make the router reproducible | data filters, overlap audit, genome encoding, GA/RL reward, score aggregation |

## 13. Writing order

1. Freeze the figure order.
2. Freeze the claim boundary wording.
3. Draft Results section headers only.
4. Draft Methods module headers.
5. Assemble the comparator tables.
6. Write the discussion only after the claim surface is fixed.

## 14. Completeness check

- One main claim family.
- One freeze point.
- One operating champion.
- One locked stress test.
- One external comparator landscape.
- One prospective-ready 96-well scaffold.
- No universal dominance claim.
