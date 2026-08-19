# Post-NDA Validation Protocol

## Phase 1: Freeze
- receive blinded candidate file
- confirm slot budget
- confirm endpoint definitions

## Phase 2: Run
- execute frozen scorer
- generate top-34 and optional top-20, top-30 views
- emit uncertainty and reason codes

## Phase 3: Unblind
- compare to sponsor labels or wetlab outputs
- calculate AUPRC, AUROC, precision@N, and false-positive burden

## Phase 4: Decision
- go to option/license if lift is meaningful
- otherwise archive the run and keep route warm

