# CROSS-Neo supertask v12

Generated: 2026-05-11T07:51:45

## Headline

v12 expands the backlog to 60 tasks and keeps them split into analysis, execution, reviewer, and next-step lanes.

## Core counts

- Known-answer top96: 91/96
- Top10: 10/10
- Risk sum: 0.149
- Execution wells: 96
- Control count: 12

## Task board

| track     | task_code   | task                                                                             | group               |   priority_score | task_id      | status   | claim_boundary         | dependency          | why_it_matters                    |
|:----------|:------------|:---------------------------------------------------------------------------------|:--------------------|-----------------:|:-------------|:---------|:-----------------------|:--------------------|:----------------------------------|
| analysis  | t12-02      | split the reviewer rebuttal into evidence, boundary, and response sublists       | reviewer defense    |            0.960 | T12-E75CA1E7 | ready    | analysis / claim-safe  | v11 board           | reduces rebuttal latency          |
| analysis  | t12-06      | split the claim-boundary table into allowed, locked, and forbidden columns       | claim safety        |            0.950 | T12-22879F1F | ready    | analysis / claim-safe  | v11 board           | increases task granularity        |
| reviewer  | t12-19      | turn reviewer objections into a slide-by-slide response matrix                   | reviewer defense    |            0.950 | T12-2890FFED | ready    | reviewer / claim-safe  | v10 reviewer packet | reduces rebuttal latency          |
| analysis  | t12-01      | split the external-validation work into cohort, sample, and assay sublists       | external validation |            0.940 | T12-70A009AA | ready    | analysis / claim-safe  | v11 board           | connects demo to prospective work |
| reviewer  | t12-23      | turn the claim boundary table into a locked / unlocked matrix                    | claim safety        |            0.940 | T12-56941F7C | ready    | reviewer / claim-safe  | v10 reviewer packet | increases task granularity        |
| next      | t12-35      | write a single-page what-is-not-claimed checklist                                | QA                  |            0.940 | T12-EEF8A7F2 | ready    | next-step / claim-safe | v7-v11 artifacts    | increases task granularity        |
| next      | t12-50      | prepare a reviewer-facing summary of the boundary between demo and evidence      | claim safety        |            0.940 | T12-E2565D3D | ready    | next-step / claim-safe | v7-v11 artifacts    | increases task granularity        |
| analysis  | t12-07      | split the task board into dependencies, sequencing, and deliverable tags         | project mgmt        |            0.930 | T12-BBC88895 | ready    | analysis / claim-safe  | v11 board           | increases task granularity        |
| reviewer  | t12-22      | turn the control manifest into a review-safety board                             | reviewer defense    |            0.930 | T12-34E88860 | ready    | reviewer / claim-safe  | v10 reviewer packet | reduces rebuttal latency          |
| next      | t12-34      | write a single-page what-is-claimed-now checklist                                | QA                  |            0.930 | T12-1FB9B26E | ready    | next-step / claim-safe | v7-v11 artifacts    | increases task granularity        |
| analysis  | t12-03      | split the figure-caption pack into main text, supplement, and deck captions      | figures             |            0.920 | T12-D53A08C9 | ready    | analysis / claim-safe  | v11 board           | increases task granularity        |
| reviewer  | t12-21      | turn the execution packet into a one-screen checklist                            | reviewer defense    |            0.920 | T12-332C6F73 | ready    | reviewer / claim-safe  | v10 reviewer packet | reduces rebuttal latency          |
| next      | t12-33      | write a single-page do-not-overclaim checklist                                   | QA                  |            0.920 | T12-3CE3CDB7 | ready    | next-step / claim-safe | v7-v11 artifacts    | increases task granularity        |
| analysis  | t12-04      | split the false-positive audit into likely cause, action, and owner              | failure registry    |            0.910 | T12-280BDF70 | ready    | analysis / claim-safe  | v11 board           | increases task granularity        |
| execution | t12-14      | prepare a bench printout for the 96-well map with blank entry columns            | lab handoff         |            0.910 | T12-3E603809 | ready    | execution / claim-safe | v6 execution packet | improves operational traceability |
| next      | t12-27      | queue a minimal external validation pilot for the top three candidates           | next step           |            0.910 | T12-91177D28 | ready    | next-step / claim-safe | v7-v11 artifacts    | increases task granularity        |
| next      | t12-46      | prepare a reviewer-facing summary of how controls block overclaim drift          | reviewer defense    |            0.910 | T12-70260FBF | ready    | next-step / claim-safe | v7-v11 artifacts    | reduces rebuttal latency          |
| next      | t12-60      | split the task explosion into a reusable template for future packets             | QA                  |            0.910 | T12-C476614C | ready    | next-step / claim-safe | v7-v11 artifacts    | increases task granularity        |
| analysis  | t12-08      | split the control manifest into bench, reviewer, and pitch controls              | controls            |            0.900 | T12-A026472D | ready    | analysis / claim-safe  | v11 board           | increases task granularity        |
| execution | t12-11      | build a cohort-specific validation queue for the highest-priority candidates     | lab handoff         |            0.900 | T12-124187FD | ready    | execution / claim-safe | v6 execution packet | improves operational traceability |
| execution | t12-15      | prepare a bench printout for the reagent order list                              | lab handoff         |            0.900 | T12-90663B06 | ready    | execution / claim-safe | v6 execution packet | improves operational traceability |
| reviewer  | t12-20      | turn the five red wells into a visible failure narrative                         | reviewer defense    |            0.900 | T12-E5B7F298 | ready    | reviewer / claim-safe  | v10 reviewer packet | reduces rebuttal latency          |
| next      | t12-36      | write a single-page table of file anchors for every figure                       | QA                  |            0.900 | T12-9A41D4FF | ready    | next-step / claim-safe | v7-v11 artifacts    | increases task granularity        |
| next      | t12-41      | prepare a reviewer-facing summary of why 91/96 is impressive but bounded         | reviewer defense    |            0.900 | T12-9CBEF5B0 | ready    | next-step / claim-safe | v7-v11 artifacts    | reduces rebuttal latency          |
| next      | t12-42      | prepare a reviewer-facing summary of why 10/10 matters but is retrospective      | reviewer defense    |            0.900 | T12-4F280198 | ready    | next-step / claim-safe | v7-v11 artifacts    | reduces rebuttal latency          |
| next      | t12-56      | split the task board into now / next / later lanes                               | project mgmt        |            0.900 | T12-6350AD5C | ready    | next-step / claim-safe | v7-v11 artifacts    | increases task granularity        |
| analysis  | t12-05      | split the known-answer board into top10, top24, top48, and top96 views           | known-answer        |            0.890 | T12-39E5ECC4 | ready    | analysis / claim-safe  | v11 board           | increases task granularity        |
| execution | t12-13      | build an assay-level handoff queue with control requirements                     | lab handoff         |            0.890 | T12-84760E52 | ready    | execution / claim-safe | v6 execution packet | improves operational traceability |
| reviewer  | t12-24      | turn the score-board into a ranking comparison figure                            | reviewer defense    |            0.890 | T12-7D58A4F4 | ready    | reviewer / claim-safe  | v10 reviewer packet | reduces rebuttal latency          |
| next      | t12-37      | write a single-page table of file anchors for every task                         | QA                  |            0.890 | T12-4A739311 | ready    | next-step / claim-safe | v7-v11 artifacts    | increases task granularity        |
| next      | t12-43      | prepare a reviewer-facing summary of why 0.149 risk sum matters                  | reviewer defense    |            0.890 | T12-3929B9C7 | ready    | next-step / claim-safe | v7-v11 artifacts    | reduces rebuttal latency          |
| next      | t12-55      | split the reviewer objection pack into standard and deep-dive responses          | reviewer defense    |            0.890 | T12-A6A265A0 | ready    | next-step / claim-safe | v7-v11 artifacts    | reduces rebuttal latency          |
| next      | t12-57      | split the task board into must-do / should-do / reserve lanes                    | project mgmt        |            0.890 | T12-41E253A8 | ready    | next-step / claim-safe | v7-v11 artifacts    | increases task granularity        |
| analysis  | t12-10      | split the endpoint slide into power, null risk, and unlock rules                 | preregistration     |            0.880 | T12-CA041206 | ready    | analysis / claim-safe  | v11 board           | increases task granularity        |
| execution | t12-12      | build a sample-level handoff queue with expected readout types                   | lab handoff         |            0.880 | T12-5CD843A6 | ready    | execution / claim-safe | v6 execution packet | improves operational traceability |
| execution | t12-16      | prepare a bench printout for the candidate response entry template               | lab handoff         |            0.880 | T12-EE6F0B15 | ready    | execution / claim-safe | v6 execution packet | improves operational traceability |
| reviewer  | t12-26      | turn the mosaic into a single screenshot page                                    | visual bundle       |            0.880 | T12-03E544AB | ready    | reviewer / claim-safe  | v10 reviewer packet | increases task granularity        |
| next      | t12-28      | queue a minimal rescue pilot for the five false-positive wells                   | next step           |            0.880 | T12-FAC05859 | ready    | next-step / claim-safe | v7-v11 artifacts    | increases task granularity        |
| next      | t12-38      | write a single-page dependency list for all backlog tasks                        | QA                  |            0.880 | T12-AA306166 | ready    | next-step / claim-safe | v7-v11 artifacts    | increases task granularity        |
| next      | t12-44      | prepare a reviewer-facing summary of why 96 wells matters operationally          | reviewer defense    |            0.880 | T12-A1A87E3E | ready    | next-step / claim-safe | v7-v11 artifacts    | reduces rebuttal latency          |
| next      | t12-45      | prepare a reviewer-facing summary of why 114 order lines matter operationally    | reviewer defense    |            0.880 | T12-046AC4CB | ready    | next-step / claim-safe | v7-v11 artifacts    | reduces rebuttal latency          |
| next      | t12-54      | split the control matrix into active and locked controls                         | controls            |            0.880 | T12-5078767D | ready    | next-step / claim-safe | v7-v11 artifacts    | increases task granularity        |
| next      | t12-58      | split the task board into bench / analysis / figure lanes                        | project mgmt        |            0.880 | T12-3581E157 | ready    | next-step / claim-safe | v7-v11 artifacts    | increases task granularity        |
| analysis  | t12-09      | split the figure mosaic into editor, reviewer, and lab views                     | visual bundle       |            0.870 | T12-D2C19BD2 | ready    | analysis / claim-safe  | v11 board           | increases task granularity        |
| execution | t12-18      | prepare a bench printout for the claim boundary ledger                           | lab handoff         |            0.870 | T12-9A9B7BB9 | ready    | execution / claim-safe | v6 execution packet | improves operational traceability |
| next      | t12-29      | queue a minimal QC pilot for the positive-control arm                            | next step           |            0.870 | T12-D89D496D | ready    | next-step / claim-safe | v7-v11 artifacts    | increases task granularity        |
| next      | t12-40      | write a single-page dependency-to-deliverable crosswalk                          | QA                  |            0.870 | T12-9F01D698 | ready    | next-step / claim-safe | v7-v11 artifacts    | increases task granularity        |
| next      | t12-47      | prepare a reviewer-facing summary of how failure wells improve the story         | reviewer defense    |            0.870 | T12-9CF41658 | ready    | next-step / claim-safe | v7-v11 artifacts    | reduces rebuttal latency          |
| next      | t12-59      | split the review response into concise and detailed views                        | reviewer defense    |            0.870 | T12-28D63A1F | ready    | next-step / claim-safe | v7-v11 artifacts    | reduces rebuttal latency          |
| execution | t12-17      | prepare a bench printout for the failure-action registry                         | lab handoff         |            0.860 | T12-85667E1F | ready    | execution / claim-safe | v6 execution packet | improves operational traceability |
| next      | t12-30      | queue a minimal specificity pilot for the hard-negative arm                      | next step           |            0.860 | T12-8BC285E0 | ready    | next-step / claim-safe | v7-v11 artifacts    | increases task granularity        |
| next      | t12-39      | write a single-page owner map for all backlog tasks                              | QA                  |            0.860 | T12-8EF92825 | ready    | next-step / claim-safe | v7-v11 artifacts    | increases task granularity        |
| next      | t12-48      | prepare a reviewer-facing summary of how task granularity improves actionability | reviewer defense    |            0.860 | T12-9B767A11 | ready    | next-step / claim-safe | v7-v11 artifacts    | reduces rebuttal latency          |
| next      | t12-31      | queue a minimal provenance-audit pilot for the label-noise arm                   | next step           |            0.850 | T12-7DA9C7AF | ready    | next-step / claim-safe | v7-v11 artifacts    | increases task granularity        |
| next      | t12-49      | prepare a reviewer-facing summary of how the packet can be executed next week    | reviewer defense    |            0.850 | T12-FCE351AF | ready    | next-step / claim-safe | v7-v11 artifacts    | reduces rebuttal latency          |
| reviewer  | t12-25      | turn the impact timeline into a sequencing chart                                 | reviewer defense    |            0.840 | T12-A76EFB1F | ready    | reviewer / claim-safe  | v10 reviewer packet | reduces rebuttal latency          |
| next      | t12-32      | queue a minimal model-boundary pilot for the TCR/structure arm                   | next step           |            0.840 | T12-26F2E0A7 | ready    | next-step / claim-safe | v7-v11 artifacts    | increases task granularity        |
| next      | t12-51      | split the deck into title, metrics, controls, and appendix sections              | figures             |            0.830 | T12-D27BE925 | ready    | next-step / claim-safe | v7-v11 artifacts    | increases task granularity        |
| next      | t12-52      | split the deck into editor, reviewer, and lab versions                           | figures             |            0.820 | T12-A866C439 | ready    | next-step / claim-safe | v7-v11 artifacts    | increases task granularity        |
| next      | t12-53      | split the figure bundle into 1-slide, 3-slide, and 6-slide cuts                  | figures             |            0.810 | T12-A7A871C0 | ready    | next-step / claim-safe | v7-v11 artifacts    | increases task granularity        |

## Control matrix

| control_id               | purpose                                       | anchor_task   | boundary_state   |   linked_tasks |
|:-------------------------|:----------------------------------------------|:--------------|:-----------------|---------------:|
| claim_safety_ledger      | allowed vs forbidden claims side by side      | t12-23        | locked           |              1 |
| bench_printouts          | bench-readable files for the next run         | t12-14        | active           |              1 |
| reviewer_response_matrix | slide-by-slide objection responses            | t12-19        | active           |              1 |
| failure_registry         | five false positives turned into action items | t12-20        | active           |              1 |
| task_splitter            | break work into smaller, owner-tagged pieces  | t12-06        | active           |              1 |
| qa_one_pagers            | single-page QA checklists                     | t12-33        | active           |              1 |
| deck_cuts                | editor / reviewer / lab deck variants         | t12-52        | active           |              1 |
| control_lock             | controls clearly labeled active vs locked     | t12-54        | locked           |              1 |
| external_validation      | prospective pilot queue                       | t12-27        | active           |              1 |
| mosaic_bundle            | single screenshot bundle                      | t12-53        | active           |              1 |
| risk_sum_lock            | confirmatory risk sum remains fixed           | t12-43        | locked           |              1 |
| order_line_trace         | 114-line execution trace                      | t12-15        | active           |              1 |

## Reviewer objections

| objection                                              | response                                                                                                | anchor                      | task_code   | risk             |
|:-------------------------------------------------------|:--------------------------------------------------------------------------------------------------------|:----------------------------|:------------|:-----------------|
| This is just known-label cherry-picking.               | Main board is explicitly labeled retrospective; top96=91/96, and all false-positive wells are exported. | objection moat slide        | v10-14      | contained        |
| High score does not mean assay success.                | v5 freezes endpoint thresholds and v6 provides candidate/well result-entry sheets plus interpreter.     | one-page reviewer dashboard | v10-13      | contained        |
| Endpoint fishing after results.                        | Confirmatory thresholds are predeclared; null-risk sum=0.149.                                           | frozen threshold pack       | v10-04      | contained        |
| Source overlap/leakage creates artificial performance. | Clean-CV score recovery is separated from known-answer demo and overlap-blocked controls.               | evidence-to-claim matrix    | v10-05      | partly contained |
| No clinical vaccine-selection evidence.                | The clinical claim is explicitly locked in the evidence matrix.                                         | claim boundary table        | v10-08      | locked           |
| The task list is still too coarse.                     | v12 expands to 60 tasks and keeps them tagged by track and dependency.                                  | task board                  | t12-06      | contained        |
| The reviewer defense is still not granular enough.     | v12 adds slide-by-slide responses and claim-safety ledgers.                                             | reviewer response matrix    | t12-19      | contained        |
| The bench handoff is too abstract.                     | v12 adds bench printouts for map, reagent, response, and failure registries.                            | bench printouts             | t12-14      | contained        |
| You still need more next steps.                        | v12 adds 33 next-step and QA microtasks.                                                                | next steps                  | t12-33      | contained        |

## Boundary

All tasks are scaffolds. They do not upgrade the retrospective known-answer board into prospective validation, and they do not create a clinical vaccine-selection claim.
