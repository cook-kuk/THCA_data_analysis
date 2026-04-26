# Quantum ML full benchmark (publication-style)

Classical LogReg_l2 reference AUC = **1.0000**.

| Algorithm | Kind | CV AUC [95% CI] | PR-AUC | Wall (s) | Verdict vs classical |
|---|---|---|---|---|---|
| LogReg_l2 (classical baseline) | classical | 1.000 [1.000, 1.000] | 1.000 | 0.0 | TIE or WIN (CI above classical baseline) |
| MLP (classical baseline) | classical | 0.981 [0.955, 0.997] | 0.995 | 0.2 | LOSE |
| VQC (ZZ feature map) | quantum | 0.488 [0.410, 0.571] | 0.841 | 88.0 | LOSE |
| VQC (PAULI feature map) | quantum | 0.782 [0.729, 0.833] | 0.955 | 44.6 | LOSE |
| VQC (AMP feature map) | quantum | 0.844 [0.782, 0.902] | 0.964 | 34.2 | LOSE |
| QSVM (fidelity kernel) | quantum | 0.754 [0.685, 0.818] | 0.935 | 76.6 | LOSE |
| QNN (EstimatorQNN + 3 ansatz layers) | quantum | 0.594 [0.516, 0.672] | 0.881 | 25.6 | LOSE |
| Quantum Kitchen Sinks | quantum-inspired | 0.491 [0.420, 0.565] | 0.835 | 23.3 | LOSE |
| QAOA feature selection -> LogReg | quantum | 1.000 [1.000, 1.000] | 1.000 | 0.0 | TIE or WIN (CI above classical baseline) |
| dwave-neal SA feature selection -> LogReg | classical-SA | 1.000 [1.000, 1.000] | 1.000 | 0.0 | TIE or WIN (CI above classical baseline) |
| Grover-inspired AA (demo, sqrt(N) asymptotic) | quantum-demo | 1.000 [1.000, 1.000] | 1.000 | 0.0 | TIE or WIN (CI above classical baseline) |
| Quantum k-means (swap-test, k=2) | quantum-unsupervised | — | nan | 1.6 | N/A (unsupervised) |
| QBoost (16 weak stumps, QUBO w/ neal) | quantum-inspired-ensemble | 0.976 [0.962, 0.989] | 0.993 | 0.2 | LOSE |

## Honesty clause
- `dwave-neal` is a **classical** simulated annealer; it is NOT quantum hardware.
- QAOA + Aer runs on a noiseless classical simulator (`qiskit_aer.primitives.Sampler`).
- PennyLane `lightning.qubit` is also a classical statevector simulator.
- All AUC CIs are 95% percentile bootstraps with 800 resamples on held-out CV predictions.
- Grover-inspired demo has O(sqrt(N)) asymptotic speedup but at n=55 genes that is irrelevant.
- Quantum k-means is unsupervised, so AUC is not reported; ARI vs label is the metric.

## qPCA vs classical PCA (first 8 eigenvalues)

| Rank | Classical PCA | Quantum-style PCA (rho eig) |
|---|---|---|
| 1 | 13.1760 | 0.9075 |
| 2 | 7.8659 | 0.0242 |
| 3 | 5.2459 | 0.0109 |
| 4 | 4.2053 | 0.0090 |
| 5 | 3.0079 | 0.0059 |
| 6 | 1.9111 | 0.0049 |
| 7 | 1.6845 | 0.0033 |
| 8 | 1.4480 | 0.0030 |

## Grover demo
- qubits: 3, iterations: 2
- target index: 4
- target probability after AA: 0.9453124999999959
