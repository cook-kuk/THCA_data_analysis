# RunPod Queue Status

Date: 2026-05-10

## Remote

- host: `135.84.176.142:20878`
- path: `/runpod-volume/cross_neo_md_diverse_2026_05_10`
- queued script: `queue_after_pid.sh 2213 run_immediate_ready_batch.sh`
- queue watcher PID: `3290`

## Current Dependency

- active HMTEVVRHC 10 ns OpenMM PID: `1956`
- existing HMTEVVRHC replicate-screen watcher PID: `2213`
- the GADGVGKSAL diverse batch waits for PID `2213`, so it starts only after the existing HMTEVVRHC current run and its already queued short screens finish.

## Batch Waiting To Start

- `run_immediate_ready_batch.sh`
- 5 jobs, 50 ns total requested trajectory length
- content: 2 additional GADGVGKSAL mutant TCR-pMHC 10 ns replicates plus 3 pMHC same-HLA positive-control 10 ns replicates

## Boundary

This is structural robustness and control simulation work only. It is not immunogenicity, vaccine efficacy, or clinical validation evidence.
