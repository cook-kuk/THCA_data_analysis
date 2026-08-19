---
title: "BioDarwin submission blueprint"
date: 2026-05-11
status: working scaffold
purpose: convert the BioDarwin strategy dossier into a concrete submission-ready manuscript blueprint
---

# BioDarwin submission blueprint

## 1. Submission thesis

BioDarwin is a **regime-aware routing system** for biomedical predictors. The paper should show that a learned router beats lookup on a frozen public pool, remains useful on a locked industrial mini-set, and can be interpreted prospectively with a 96-well plate scaffold.

This is a paper about **where to route a sample**, not a claim that one model dominates every dataset.

## 2. Title candidates

| rank | title | why it works |
|---|---|---|
| 1 | BioDarwin: regime-aware routing for biomedical expert selection | clean, accurate, broad |
| 2 | BioDarwin improves biomedical expert selection with regime-aware routing | emphasizes result and method |
| 3 | BioDarwin learns when to trust different biomedical experts | more readable, slightly less technical |
| 4 | BioDarwin: GA/RL-guided mixture-of-experts routing for biomedical prediction | stronger method signal, more technical |

## 3. Abstract contract

The abstract should contain only these elements:

1. Problem: biomedical benchmarks are heterogeneous and no single expert wins everywhere.
2. Method: BioDarwin uses public metadata, GA/RL search, LOBO blending, and contextual MoE routing.
3. Result 1: frozen public bundles improve from lookup to the current operating champion.
4. Result 2: locked industrial TG4050 control beats BigMHC_IM in diagnostic mode.
5. Result 3: 96-well plate/interpreter defines hit/fail/helper/probe lanes for prospective use.
6. Boundary: not universal SOTA; not prospective superiority until a new locked batch is tested.

## 4. Main-text skeleton

| section | job | key content |
|---|---|---|
| Title / Abstract | say what the paper is | router paper, not universal SOTA |
| Introduction | explain why the problem exists | dataset-specific winners, heterogeneity, need for routing |
| Results 1 | establish baseline failure | lookup is not enough |
| Results 2 | show learning works | GA/RL -> LOBO -> contextual MoE |
| Results 3 | show locked industrial stress test | BigMHC_IM comparator, diagnostic claim only |
| Results 4 | show prospective-ready plate | 96-well interpreter and candidate lanes |
| Results 5 | show external landscape | reviewer-safe vs performance-max comparison |
| Discussion | close the claim boundary | saturation, new pool needed, split paper logic |
| Methods | reproduce the system | filters, overlap audit, genome, reward, scoring |

## 5. Figure sequence

| figure | role | evidence |
|---|---|---|
| Fig 1 | routing problem | benchmark map and regimes |
| Fig 2 | learning curve | v1 lookup -> v4 contextual MoE |
| Fig 3 | plate reader | 96-well interpreter |
| Fig 4 | external landscape | public benchmark winners by dataset |
| Fig 5 | genome/weights | what GA/RL learned |
| Fig 6 | failure modes | what is not claimed |

## 6. Results wording contract

Use the following phrasing discipline:

- say "improves selection quality"
- say "wins the locked industrial mini-set in diagnostic mode"
- say "same-pool saturation"
- say "new expert pool or new data required"
- do not say "universal dominance"
- do not say "prospective superiority"

## 7. Supplement contract

Keep these out of the main text:

- overlap audit details
- reviewer-safe comparator landscape tables
- low-power / all-positive bundles
- stress tests that only support caveats
- RF_biophys anchor diagnostics

## 8. Reviewer defense map

| reviewer question | answer |
|---|---|
| Is this just a tuned benchmark? | no, the main claim is routing under explicit regime control |
| Why not claim SOTA everywhere? | because dataset-specific winners differ and the correct model class is a router |
| What about overlap? | handled explicitly in reviewer-safe tables |
| Where is prospective validation? | the 96-well interpreter is prospective-ready; a new locked batch is still needed |
| Why keep RF_biophys? | as an anchor/disagreement feature, not as the final model |

## 9. Manuscript order of operations

1. Freeze the title.
2. Freeze the claim boundary.
3. Lock the figure order.
4. Draft section headers only.
5. Build the comparator tables.
6. Write the discussion last.

## 10. Hard framing rules

- Keep the paper in the selection class, not the universal predictor class.
- Keep the 96-well plate as a prospective scaffold, not as completed wetlab proof.
- Keep dataset-specific winners in a separate comparator landscape.
- Keep RF_biophys as an anchor feature, not the final model.
- Freeze v3 as the claim-safe point and v4 as the operating champion.

## 11. Abstract sentence map

| sentence | job |
|---|---|
| 1 | define the heterogeneous routing problem |
| 2 | name the BioDarwin method stack |
| 3 | state the frozen public-pool improvement |
| 4 | state the industrial locked diagnostic win |
| 5 | state the 96-well prospective scaffold |
| 6 | close the claim boundary |

## 12. Results headers

1. `Results 1 | A lookup baseline is not enough`
2. `Results 2 | GA/RL and contextual routing improve the frozen public pool`
3. `Results 3 | The locked industrial mini-set is solved in diagnostic mode`
4. `Results 4 | The 96-well interpreter converts ranking into execution lanes`
5. `Results 5 | Public benchmark winners are dataset-specific`
6. `Results 6 | The current expert pool is saturated`

## 13. Figure-to-section map

| figure | manuscript section | what it should prove |
|---|---|---|
| Fig 1 | Introduction + Results 1 | why a router is necessary |
| Fig 2 | Results 2 | why learning beats lookup |
| Fig 3 | Results 4 | how the 96-well lanes are interpreted |
| Fig 4 | Results 5 | why external winners must be routed |
| Fig 5 | Methods + Results 2 | what the genome learned |
| Fig 6 | Discussion | what the paper is not allowed to claim |

## 14. Forbidden sentence patterns

- one model wins every dataset
- industrial diagnostic results are already prospective proof
- the external landscape can be collapsed into one number
- the 96-well plate is a completed wetlab validation
- same-pool search will keep improving indefinitely

## 15. Supplement placement map

| supplement block | content |
|---|---|
| overlap audit | public-training and comparator leakage checks |
| comparator landscape | dataset-specific winner tables |
| stress tests | low-power, all-positive, or diagnostic-only rows |
| anchor diagnostics | RF_biophys / disagreement traces |
| prospective scaffold | 96-well helper/probe lanes and frozen trace |

## 16. Reviewer kill-switch

If any reviewer-facing sentence sounds like universal SOTA, replace it with a routing sentence:

- best single expert on this frozen pool
- dataset-specific winner
- operating champion in the current pool
- claim-safe freeze point
- new pool needed for the next jump

## 17. Execution checklist

- [ ] Title selected
- [ ] Abstract contract approved
- [ ] Section headers frozen
- [ ] Fig-to-section map frozen
- [ ] Comparator landscape separated
- [ ] Supplement boundary frozen
- [ ] Reviewer kill-switch inserted
