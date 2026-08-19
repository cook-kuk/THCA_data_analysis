# CROSS-Neo task expansion v11

Generated: 2026-05-11T07:38:09

## Headline

v11 expands the backlog to 32 tasks and keeps them split across analysis, execution, reviewer, and next-step lanes.

## Core counts

- Known-answer top96: 91/96
- Top10: 10/10
- Risk sum: 0.149
- Execution wells: 96
- Control count: 10

## Task board

| track     | task_code   | task                                                                            | group               |   priority_score | why_it_matters                      | task_id      | status   | claim_boundary                     | dependency          |
|:----------|:------------|:--------------------------------------------------------------------------------|:--------------------|-----------------:|:------------------------------------|:-------------|:---------|:-----------------------------------|:--------------------|
| analysis  | t11-06      | split the 20-task backlog into mini-tasks by owner and dependency               | project mgmt        |            0.950 | more granular execution             | T11-0AD48FA4 | ready    | analysis / reviewer-safe           | v7-v10 artifacts    |
| reviewer  | t11-17      | show which claims are allowed now versus locked for later                       | claim safety        |            0.940 | clean boundary discipline           | T11-03357044 | ready    | reviewer defense / claim-safe      | v10 reviewer packet |
| analysis  | t11-03      | write a rebuttal matrix for the five most likely reviewer objections            | reviewer defense    |            0.930 | pre-answered objections             | T11-178AF38A | ready    | analysis / reviewer-safe           | v7-v10 artifacts    |
| next      | t11-27      | add a reproducibility checklist for every table in the packet                   | QA                  |            0.930 | lowers reviewer friction            | T11-D4E196E5 | queued   | task backlog only                  | v7-v10 artifacts    |
| analysis  | t11-01      | write an external-validation task list for K2/Bundang-style follow-up           | external validation |            0.920 | retrospective to prospective bridge | T11-A913B5AA | ready    | analysis / reviewer-safe           | v7-v10 artifacts    |
| reviewer  | t11-15      | add an objection-response table with exact evidence file anchors                | reviewer defense    |            0.920 | faster rebuttal                     | T11-F7A117BC | ready    | reviewer defense / claim-safe      | v10 reviewer packet |
| next      | t11-30      | create a reviewer-facing summary of what is deliberately not claimed            | QA                  |            0.920 | tightens scope                      | T11-DDA3836D | queued   | task backlog only                  | v7-v10 artifacts    |
| analysis  | t11-02      | write a figure-caption expansion pack for the strongest visuals                 | figures             |            0.910 | editor-ready caption workflow       | T11-82C00CEF | ready    | analysis / reviewer-safe           | v7-v10 artifacts    |
| reviewer  | t11-18      | separate positive controls from discovery lanes in a visible matrix             | reviewer defense    |            0.910 | prevents benchmark confusion        | T11-463705D6 | ready    | reviewer defense / claim-safe      | v10 reviewer packet |
| next      | t11-31      | create a single-slide narrative of why 91/96 matters and what it does not prove | reviewer defense    |            0.910 | high-impact but bounded             | T11-94291DFC | queued   | task backlog only                  | v7-v10 artifacts    |
| analysis  | t11-05      | add a claim-boundary ledger mapping every allowed claim to a forbidden upgrade  | claim safety        |            0.900 | prevents overclaim drift            | T11-FF393155 | ready    | analysis / reviewer-safe           | v7-v10 artifacts    |
| execution | t11-11      | prepare a per-candidate checklist that maps assay conditions to readouts        | lab handoff         |            0.900 | less ambiguity at bench             | T11-B4D861B9 | ready    | execution scaffold / reviewer-safe | v6 execution packet |
| reviewer  | t11-16      | build a red-flag registry for leakage, fishing, and claim inflation             | reviewer defense    |            0.900 | transparent risk handling           | T11-4E375345 | ready    | reviewer defense / claim-safe      | v10 reviewer packet |
| next      | t11-25      | prepare a decision-tree for what changes if the next assay is negative          | next step           |            0.900 | pre-commit response plan            | T11-C6DA9A7B | queued   | task backlog only                  | v7-v10 artifacts    |
| execution | t11-09      | create a public-safe result-entry template for future assay readout             | lab handoff         |            0.890 | operational cleanliness             | T11-DA97013B | ready    | execution scaffold / reviewer-safe | v6 execution packet |
| reviewer  | t11-19      | map the top 5 false-positive wells to action items                              | failure registry    |            0.890 | uses misses constructively          | T11-CE320A6C | ready    | reviewer defense / claim-safe      | v10 reviewer packet |
| next      | t11-26      | prepare a decision-tree for what changes if the next assay is positive          | next step           |            0.890 | unlock path remains visible         | T11-9748CB6E | queued   | task backlog only                  | v7-v10 artifacts    |
| analysis  | t11-04      | build a candidate triage table from the top96 false-positive wells              | failure registry    |            0.880 | turn misses into action             | T11-BD48DAB7 | ready    | analysis / reviewer-safe           | v7-v10 artifacts    |
| execution | t11-10      | create a control-sequence list for positive, negative, and provenance checks    | lab handoff         |            0.880 | ordering without redesign           | T11-4137B67E | ready    | execution scaffold / reviewer-safe | v6 execution packet |
| next      | t11-32      | build a ‘do not overcall’ checklist for every future update                     | QA                  |            0.880 | keeps claims honest                 | T11-EC2E4DE7 | queued   | task backlog only                  | v7-v10 artifacts    |
| analysis  | t11-08      | expand the endpoint plan into a slide with power and risk thresholds            | preregistration     |            0.870 | reviewer-friendly stats             | T11-C0AFCA0C | ready    | analysis / reviewer-safe           | v7-v10 artifacts    |
| reviewer  | t11-21      | make a dense mosaic for a single-screenshot review                              | visual bundle       |            0.870 | fast scan path                      | T11-9C983253 | ready    | reviewer defense / claim-safe      | v7-v10 artifacts    |
| execution | t11-12      | generate a per-well task map that can be printed for the bench                  | lab handoff         |            0.860 | well-level traceability             | T11-732604FF | ready    | execution scaffold / reviewer-safe | v6 execution packet |
| next      | t11-23      | prepare a mini-experiment backlog for top three high-priority candidates        | next step           |            0.860 | focus on likely wins                | T11-86715E7A | queued   | task backlog only                  | v7-v10 artifacts    |
| reviewer  | t11-20      | draft a one-slide summary that compares v3, v5, v7, v8, v10                     | summary             |            0.850 | trajectory clarity                  | T11-3B0A4DDB | ready    | reviewer defense / claim-safe      | v7-v10 artifacts    |
| analysis  | t11-07      | expand the known-answer board into a score-by-source breakdown                  | known-answer        |            0.840 | more informative demo board         | T11-249540E2 | ready    | analysis / reviewer-safe           | v7-v10 artifacts    |
| execution | t11-14      | produce a standalone figure index that points to every key image                | packaging           |            0.840 | fast navigation                     | T11-4F8A5CD3 | ready    | execution scaffold / reviewer-safe | v7-v10 artifacts    |
| next      | t11-24      | prepare a rescue backlog for the false-negative / rescue wells                  | next step           |            0.840 | convert uncertainty to a queue      | T11-FB54C844 | queued   | task backlog only                  | v7-v10 artifacts    |
| reviewer  | t11-22      | add a version timeline from score recovery to reviewer response                 | visual bundle       |            0.830 | shows momentum                      | T11-C9D7EE22 | ready    | reviewer defense / claim-safe      | v7-v10 artifacts    |
| execution | t11-13      | bundle the execution packet and reviewer packet into one delivery zip           | packaging           |            0.820 | single handoff artifact             | T11-4A48FFBC | ready    | execution scaffold / reviewer-safe | v7-v10 artifacts    |
| next      | t11-29      | crosswalk the task board to existing v6/v7/v8/v9 assets                         | QA                  |            0.810 | ensures traceability                | T11-5BC1ED8E | queued   | task backlog only                  | v7-v10 artifacts    |
| next      | t11-28      | add a glossary of claim-boundary language used across packets                   | QA                  |            0.800 | prevents inconsistent wording       | T11-D57FC140 | queued   | task backlog only                  | v7-v10 artifacts    |

## Control matrix

| control_id                | purpose                                                          | anchor_task   | boundary_state   |   linked_tasks |
|:--------------------------|:-----------------------------------------------------------------|:--------------|:-----------------|---------------:|
| external_validation_block | separate prospective external validation from retrospective demo | t11-01        | locked           |              1 |
| figure_caption_pack       | attach explicit boundary language to every high-impact figure    | t11-02        | active           |              1 |
| reviewer_rebuttal_map     | pre-answer objections with file anchors                          | t11-03        | active           |              1 |
| failure_to_action         | use false positives as a queue, not a discarded appendix         | t11-04        | active           |              1 |
| claim_boundary_ledger     | keep allowed claims and forbidden upgrades side by side          | t11-05        | locked           |              1 |
| mini_task_split           | split the backlog into granular owner/dependency items           | t11-06        | active           |              1 |
| score_by_source           | show how the score behaves by source slice                       | t11-07        | active           |              1 |
| power_risk_slide          | freeze the endpoint risk slide                                   | t11-08        | locked           |              1 |
| result_entry_template     | standardize the future assay readout entry                       | t11-09        | active           |              1 |
| control_sequence_list     | track positive / negative / provenance controls                  | t11-10        | active           |              1 |

## Reviewer objection pack

| objection                                              | response                                                                                                | anchor                      | task_code   | risk             |
|:-------------------------------------------------------|:--------------------------------------------------------------------------------------------------------|:----------------------------|:------------|:-----------------|
| This is just known-label cherry-picking.               | Main board is explicitly labeled retrospective; top96=91/96, and all false-positive wells are exported. | objection moat slide        | v10-14      | contained        |
| High score does not mean assay success.                | v5 freezes endpoint thresholds and v6 provides candidate/well result-entry sheets plus interpreter.     | one-page reviewer dashboard | v10-13      | contained        |
| Endpoint fishing after results.                        | Confirmatory thresholds are predeclared; null-risk sum=0.149.                                           | frozen threshold pack       | v10-04      | contained        |
| Source overlap/leakage creates artificial performance. | Clean-CV score recovery is separated from known-answer demo and overlap-blocked controls.               | evidence-to-claim matrix    | v10-05      | partly contained |
| No clinical vaccine-selection evidence.                | The clinical claim is explicitly locked in the evidence matrix.                                         | claim boundary table        | v10-08      | locked           |
| Need more tasks and finer granularity.                 | v11 splits the backlog into 32 tasks with dependencies and boundary tags.                               | reviewer task board         | t11-06      | contained        |
| This is still too optimistic.                          | v11 includes a do-not-overcall checklist and a what-is-not-claimed summary.                             | QA packet                   | t11-30      | contained        |

## Boundary

All tasks are reviewer-safe scaffolds. They are not prospective wetlab validation, and they do not upgrade the clinical vaccine-selection claim.
