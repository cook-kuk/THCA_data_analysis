# Fallback Algorithm Sweep

Strict class-I set: n=89, positives=21. Multiprocessing workers: 7.

All trained rows are internal repeated-CV pilots, not external paper claims.

| rank | method | family | n / pos | AUROC | AUPRC | note |
|---:|---|---|---:|---:|---:|---|
| 1 | quantum_kernel_no_anchor_gamma1.0 | quantum_kernel | 89 / 21 | 0.785 | 0.603 | angle-encoded fidelity kernel; internal repeated CV |
| 2 | quantum_kernel_compact_gamma0.5 | quantum_kernel | 89 / 21 | 0.707 | 0.457 | angle-encoded fidelity kernel; internal repeated CV |
| 3 | Structure_LR | fixed | 89 / 21 | 0.679 | 0.435 | locked existing score |
| 4 | lr_compact_reference_replacement | cv | 89 / 21 | 0.660 | 0.412 | small v0 candidate |
| 5 | MHCflurry | fixed | 89 / 21 | 0.657 | 0.488 | locked existing score |
| 6 | rf_compact_reference_replacement | cv | 89 / 21 | 0.654 | 0.446 | small v0 candidate |
| 7 | quantum_kernel_compact_gamma0.25 | quantum_kernel | 89 / 21 | 0.638 | 0.374 | angle-encoded fidelity kernel; internal repeated CV |
| 8 | Wave8_TCR_SelfSim_full | fixed | 89 / 21 | 0.632 | 0.403 | locked existing score |
| 9 | Wave8_TCR_SelfSim_no_exact | fixed | 89 / 21 | 0.629 | 0.472 | locked existing score |
| 10 | Wave8_TCR_motif_only | fixed | 89 / 21 | 0.625 | 0.462 | locked existing score |
| 11 | Stack_mean_E1 | fixed | 89 / 21 | 0.610 | 0.409 | locked existing score |
| 12 | BigMHC_IM | fixed | 89 / 21 | 0.596 | 0.288 | locked existing score |
| 13 | lr_anchor_quantum_structure | cv | 89 / 21 | 0.590 | 0.342 | internal repeated CV |
| 14 | extratrees_anchor_quantum_structure | cv | 89 / 21 | 0.588 | 0.309 | internal repeated CV |
| 15 | quantum_kernel_compact_gamma1.0 | quantum_kernel | 89 / 21 | 0.588 | 0.274 | angle-encoded fidelity kernel; internal repeated CV |
| 16 | W7A_QK_only | fixed | 89 / 21 | 0.586 | 0.352 | locked existing score |
| 17 | W7B_stacked | fixed | 89 / 21 | 0.578 | 0.294 | locked existing score |
| 18 | VQC | fixed | 89 / 21 | 0.573 | 0.309 | locked existing score |

## Decision

- Primary no-reference substitute for now: `Structure_LR` if we require a locked, non-CV score.
- Best v0 research candidate: `quantum_kernel_no_anchor_gamma1.0`, because it uses quantum scores + structure proxy + soft TCR motif only, with no public-pretrained MHCflurry/BigMHC and no `Structure_LR` anchor.
- Keep quantum as an active branch: the top quantum-kernel pilot beats the locked strict baselines in internal repeated CV.
- Do not use public-pretrained scores as clean evidence; keep them as upper-bound/caveated comparators.
- Structure proxy is useful as an uncertainty/fallback feature, not a standalone reference substitute.

## Permutation Check

`quantum_kernel_no_anchor_gamma1.0` was checked with 100 label permutations.

| metric | observed | permutation null mean | null 95% range | empirical p |
|---|---:|---:|---:|---:|
| AUROC | 0.785 | 0.495 | 0.310-0.696 | 0.010 |
| AUPRC | 0.603 | 0.276 | 0.176-0.477 | 0.010 |

This supports continuing the quantum branch. It still does not replace a leakage-controlled external/time-split benchmark.
