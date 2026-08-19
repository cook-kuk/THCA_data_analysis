# CROSS-Neo impact war room v18

Generated: 2026-05-11T08:59:17

## Headline

v18 is the board-facing war room. It is designed to be scanned once and understood.

## Core numbers

- Known-answer top96: 91/96
- Top10: 10/10
- Risk sum: 0.149
- Tasks: 120
- Controls: 20
- Claim locks: 6
- Launch checks: 6

## Scoreboard

| tile      | value    | meaning                | use               | tile_id      |
|:----------|:---------|:-----------------------|:------------------|:-------------|
| signal    | 91/96    | known-answer retrieval | signal lock       | SBW-299CD297 |
| score     | 1.93     | AUPRC recovery         | benchmark lift    | SBW-2AE017E0 |
| risk      | 0.149    | confirmatory risk      | hard boundary     | SBW-297E34DE |
| execution | 96 wells | execution packet       | lab handoff       | SBW-5624125B |
| tasks     | 120      | task density           | program load      | SBW-99CB1B7E |
| locks     | 6        | claim locks            | reviewer boundary | SBW-657BFEDA |
| deck      | war room | board-facing surface   | front door        | SBW-4CEE70AB |

## Risk register

| risk                                                   | mitigation          | status           | owner         | risk_id     |
|:-------------------------------------------------------|:--------------------|:-----------------|:--------------|:------------|
| This is just known-label cherry-picking.               |                     | contained        | board surface | RR-C18E4414 |
| High score does not mean assay success.                |                     | contained        | board surface | RR-965117D4 |
| Endpoint fishing after results.                        |                     | contained        | board surface | RR-32418712 |
| Source overlap/leakage creates artificial performance. |                     | partly contained | board surface | RR-21E4DC32 |
| No clinical vaccine-selection evidence.                |                     | locked           | board surface | RR-7A6A4F42 |
| Claims drift beyond retrospective boundary             | claim lock matrix   | contained        | narrative     | RR-93BE45FE |
| Reviewer asks for prospective proof                    | execution packet    | contained        | lab           | RR-5C125B4C |
| Deck becomes too dense                                 | war room cut        | contained        | design        | RR-DD60A004 |
| Action backlog loses priority                          | ranked action queue | contained        | program       | RR-CBF871BC |

## Lane decision board

|   lane | headline             | metric    | best_use                 | lane_id     |
|-------:|:---------------------|:----------|:-------------------------|:------------|
|      1 | retrospective signal | 91/96     | reviewer-safe demo       | LB-4750B606 |
|      2 | benchmark lift       | 91        | comparison lane          | LB-F6E96AE1 |
|      3 | preregistered lock   | 0.149     | pre-assay boundary       | LB-26989684 |
|      4 | execution ready      | 96 wells  | handoff lane             | LB-A42647EF |
|      5 | task pressure        | 120 tasks | program lane             | LB-718F163F |
|      6 | board surface        | war room  | front-door communication | LB-C9B86417 |

## Launch checklist

|   step | item                      | status   | notes                            | step_id     |
|-------:|:--------------------------|:---------|:---------------------------------|:------------|
|      1 | claim lock matrix visible | done     | boundary is explicit             | LC-A942F051 |
|      2 | risk register adjacent    | done     | reviewer objections are mapped   | LC-FB05406F |
|      3 | action queue ranked       | done     | 120 tasks remain runnable        | LC-8C63B84A |
|      4 | evidence bridge present   | done     | signal to boundary spine visible | LC-AB23AC2A |
|      5 | mosaic attached           | done     | dense visual summary is ready    | LC-37D35656 |
|      6 | one-screen board verified | done     | 6 claim locks remain active      | LC-D31E34FA |

## Evidence bridge

| layer     | value   | meaning               |
|:----------|:--------|:----------------------|
| signal    | 91/96   | known-answer          |
| score     | 1.93    | benchmark lift        |
| risk      | 0.149   | confirmatory boundary |
| execution | 96      | lab readiness         |
| density   | 120     | workload              |
| locks     | 6       | claim controls        |

## Boundary

This deck remains retrospective, preregistered, and execution-ready only.
