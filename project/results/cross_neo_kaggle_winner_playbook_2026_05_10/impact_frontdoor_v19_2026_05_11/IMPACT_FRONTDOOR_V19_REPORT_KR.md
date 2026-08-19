# CROSS-Neo impact frontdoor v19

Generated: 2026-05-11T13:17:32

## Headline

v19 is the front-door memo. It compresses the war room into a board snap-brief.

## Core numbers

- Known-answer top96: 91/96
- Top10: 10/10
- Risk sum: 0.149
- Tasks: 120
- Controls: 20
- Claim locks: 6
- Launch checks: 6

## Headline strip

| tile      | value    | meaning          | use               | tile_id      |
|:----------|:---------|:-----------------|:------------------|:-------------|
| signal    | 91/96    | known-answer     | signal lock       | HDR-9532EAFD |
| score     | 1.93     | AUPRC recovery   | benchmark lift    | HDR-2AE017E0 |
| risk      | 0.149    | risk sum         | hard boundary     | HDR-FBFDAA9B |
| execution | 96 wells | execution packet | lab handoff       | HDR-5624125B |
| locks     | 6        | claim locks      | reviewer boundary | HDR-657BFEDA |
| launch    | 6        | launch checks    | board readiness   | HDR-32AC131C |

## Verdict matrix

| verdict   | dimension            | metric   | status   | verdict_id   |
|:----------|:---------------------|:---------|:---------|:-------------|
| GO        | retrospective signal | 91/96    | allowed  | VDR-98F3707A |
| GO        | benchmark lift       | 1.93     | allowed  | VDR-47A4B429 |
| LOCK      | confirmatory risk    | 0.149    | boundary | VDR-3FEDD46B |
| GO        | execution            | 96 wells | allowed  | VDR-98C73968 |
| LOCK      | claim locks          | 6        | boundary | VDR-DF39630B |
| GO        | launch readiness     | 6        | allowed  | VDR-75C47BFB |

## Risk triage

| risk                                                   | mitigation                  | severity   | owner      | risk_id      |
|:-------------------------------------------------------|:----------------------------|:-----------|:-----------|:-------------|
| This is just known-label cherry-picking.               |                             | contained  | front door | TRI-C18E4414 |
| High score does not mean assay success.                |                             | contained  | front door | TRI-965117D4 |
| Endpoint fishing after results.                        |                             | contained  | front door | TRI-32418712 |
| Source overlap/leakage creates artificial performance. |                             | contained  | front door | TRI-6646BC4F |
| No clinical vaccine-selection evidence.                |                             | contained  | front door | TRI-6F23BC45 |
| Overclaiming the boundary                              | claim lock + verdict matrix | contained  | narrative  | TRI-FB3863BB |
| Deck is too wide                                       | front-door compression      | contained  | design     | TRI-5B1B3C92 |
| Action queue loses urgency                             | ranked next steps           | contained  | program    | TRI-AEF9228F |

## Action queue

|   rank | action                     | status   | notes                                 | action_id    |
|-------:|:---------------------------|:---------|:--------------------------------------|:-------------|
|      1 | publish front-door memo    | done     | one page, one read                    | ACT-FE1856EE |
|      2 | keep claim boundary locked | done     | no prospective overclaim              | ACT-08B29C1F |
|      3 | attach risk triage         | done     | reviewer objections adjacent          | ACT-1D123FDE |
|      4 | retain execution packet    | done     | 120 tasks preserved                   | ACT-39EC6E22 |
|      5 | retain launch checklist    | done     | 6 checks preserved                    | ACT-DA69F0FF |
|      6 | prepare board summary      | queued   | compress to single slide if requested | ACT-F5C10537 |

## Evidence spine

| layer     | value   | meaning         |
|:----------|:--------|:----------------|
| signal    | 91/96   | known-answer    |
| score     | 1.93    | benchmark lift  |
| risk      | 0.149   | risk lock       |
| execution | 96      | lab readiness   |
| locks     | 6       | claim controls  |
| launch    | 6       | board readiness |

## Boundary

This remains retrospective, preregistered, and execution-ready only.
