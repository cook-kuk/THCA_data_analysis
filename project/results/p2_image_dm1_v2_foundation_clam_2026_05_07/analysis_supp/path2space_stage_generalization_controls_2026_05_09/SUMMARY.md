# GSE250521 Stage-Generalization Controls

## Verdict

This is a caveat-bearing control.

| Layer | Result | Interpretation |
|---|---:|---|
| Existing LOSO DM1/RAI stage-centered rho | 0.392 | positive after stage centering |
| Existing LOSO DM1/RAI PT/PTC/LPTC/ATC rho | 0.246 / 0.487 / 0.471 / 0.191 | ATC is the weak stage |
| Leave-one-stage-out raw DM1/RAI pooled rho | 0.335 | raw cross-stage transfer support |
| Leave-one-stage-out residual DM1/RAI pooled rho | 0.140 | residual transfer is weaker caveat |

Safe wording: stage-centered local signal remains and raw leave-one-stage-out transfer is positive, but coord+QC residual transfer and ATC stability are weaker.
