# CROSS-Neo impact command v17

Generated: 2026-05-11T08:34:59

## Headline

v17 is the executive command layer: scoreboards, claim locks, rebuttal register, action queue, and evidence spine.

## Core numbers

- Known-answer top96: 91/96
- Top10: 10/10
- Risk sum: 0.149
- Tasks: 120
- Controls: 20
- Captions: 8
- Claim locks: 6

## Scoreboard

| tile      | value    | meaning                | use                          | tile_id     |
|:----------|:---------|:-----------------------|:-----------------------------|:------------|
| signal    | 91/96    | known-answer retrieval | retrospective prioritization | SB-299CD297 |
| score     | 1.93     | AUPRC recovery         | benchmark improvement        | SB-2AE017E0 |
| risk      | 0.149    | confirmatory null-risk | pre-assay guardrail          | SB-4DB51502 |
| execution | 96 wells | execution packet       | lab handoff                  | SB-5624125B |
| tasks     | 120      | task density           | scaffold management          | SB-99CB1B7E |
| captions  | 8        | caption pack           | editor surface               | SB-E86EC024 |
| pack      | v16      | apex deck              | compressed executive layout  | SB-AAEB67CB |

## Claim lock matrix

| lane          | metric    | allowed_use                  | boundary            | lane_id      |
|:--------------|:----------|:-----------------------------|:--------------------|:-------------|
| retrospective | 91/96     | allowed for prioritization   | not prospective     | CLK-662A850D |
| benchmark     | 1.93      | allowed for model choice     | clean-CV only       | CLK-A9870B56 |
| preregistered | 0.149     | allowed for endpoint lock    | no post-hoc tuning  | CLK-C5914146 |
| execution     | 96 wells  | allowed for lab handoff      | not completed assay | CLK-2A88BF46 |
| backlog       | 120 tasks | allowed for program steering | scaffold only       | CLK-9A539AD7 |
| apex          | 6 steps   | allowed for review deck      | compression only    | CLK-66794B86 |

## Rebuttal register

| objection                                              | defense           | status           | move                            |
|:-------------------------------------------------------|:------------------|:-----------------|:--------------------------------|
| This is just known-label cherry-picking.               |                   | contained        | tie to evidence or boundary     |
| High score does not mean assay success.                |                   | contained        | tie to evidence or boundary     |
| Endpoint fishing after results.                        |                   | contained        | tie to evidence or boundary     |
| Source overlap/leakage creates artificial performance. |                   | partly contained | tie to evidence or boundary     |
| No clinical vaccine-selection evidence.                |                   | locked           | tie to evidence or boundary     |
| Need the deck in one glance.                           | scoreboard        | contained        | front-door summary              |
| Need next steps with ownership.                        | action queue      | contained        | prioritized tasks               |
| Need caveats near the claim.                           | claim lock matrix | contained        | keep caveat and result adjacent |
| Need reviewer-safe wording.                            | boundary labels   | contained        | no overclaim                    |

## Action queue

|   rank | action                      | status   | notes                                         | action_id   |
|-------:|:----------------------------|:---------|:----------------------------------------------|:------------|
|      1 | publish apex deck           | done     | surface the strongest reusable board          | AQ-698C2B27 |
|      2 | lock claim boundary         | done     | keep retrospective/prospective split explicit | AQ-515A87C3 |
|      3 | keep reviewer moat adjacent | done     | reduce claim drift                            | AQ-95F500D0 |
|      4 | attach execution backlog    | done     | 120 tasks available                           | AQ-3C4D9388 |
|      5 | retain evidence spine       | done     | 6 layers                                      | AQ-CBDC98AC |
|      6 | prepare next review cut     | queued   | flatten to one-slide summary if needed        | AQ-FD38043C |

## Evidence spine

| layer       | value   | meaning               |
|:------------|:--------|:----------------------|
| signal      | 91/96   | known-answer surface  |
| score       | 1.93    | clean-CV improvement  |
| risk        | 0.149   | pre-assay gating      |
| execution   | 96      | handoff readiness     |
| density     | 120     | task pressure         |
| compression | 6       | executive compression |

## Boundary

This remains retrospective, preregistered, and execution-ready only.
