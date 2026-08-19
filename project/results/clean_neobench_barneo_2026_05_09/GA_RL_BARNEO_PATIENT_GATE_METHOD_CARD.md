# GA/RL Patient Gate Method Card

## What It Does

This pack takes the GA/RL unique T1 handoff candidates and passes them through the PAAD/THCA patient-gated research triage matrix.

## Claim Boundary

It does not make clinical recommendations. It ranks research triage priority only.

## Rule

- `clinical_use=false` remains the default.
- Best scenario is only a research handoff context.
- Missing disease timing, expression, clonality, and safety metadata keep confidence capped.
