# dl_quantum_summary

## Environment

- torch: ok 2.11.0+cpu
- dwave-neal: ok
- qiskit: ok 2.4.0
- pennylane: ok 0.44.1
- plotly: ok

## Setup

- Training cohort: TCGA-THCA, BRAF_like vs RAS_like (MAF-anchored)
- Samples (training): 333, class balance: {1: 279, 0: 54}
- External cohorts evaluated: ['GSE27155', 'GSE126698']
- QUBO feature selection algorithm: QUBO_SA_dwave-neal (classical simulated annealing; QAOA not installed)
- QUBO selected genes (k=14): DUSP6, KLK10, DIO1, IYD, CDH1, FOSL1, SPRY4, PIK3CA, DUSP4, MMP9, TSHR, IDO1, AKT1, NKX2-1
- Quantum classifier backend: pennylane_VQC (lightning.qubit)

## Internal CV AUC (5-fold except VQC=3-fold)

| Model | Feature set | AUC | 95% CI | Wall seconds |
|---|---|---|---|---|
| LogReg_l2 (classical) | TierA67 | 1.000 | [1.000, 1.000] | 0.505 |
| LogReg_l2 (classical) | TierA67_clean | 1.000 | [1.000, 1.000] | 0.675 |
| MLP_torch | TDS16 | 0.990 | [0.982, 0.997] | 6.578 |
| MLP_torch | TierA67 | 1.000 | [0.999, 1.000] | 3.713 |
| MLP_torch | TierA67_clean | 1.000 | [1.000, 1.000] | 3.417 |
| MLP_torch | variance_top50 | 1.000 | [1.000, 1.000] | 3.972 |
| pennylane_VQC (lightning.qubit) | TierA67_clean | 0.732 | [0.656, 0.808] | 21.013 |

## CI overlap

- MLP (TierA67_clean) CI vs LogReg_l2 CI: overlap = yes
- VQC (TierA67_clean) CI vs LogReg_l2 CI: overlap = no

## Honest verdict

- MLP did NOT improve over LogReg on internal CV
- Quantum VQC did NOT improve over LogReg (AUC 0.732 vs 1.000)

## Caveats

- The classical LogReg AUC on TierA67 is already at or near the ceiling (AUC ~1.0)
  on TCGA-THCA. There is essentially no room for a more complex model to 'win'
  on internal CV. Any observed improvement is likely noise.
- Quantum feature selection here is solved by classical simulated annealing
  (dwave-neal), not a real quantum device. The QUBO formulation is valid, but
  the quantum advantage claim would require a QAOA run on real hardware.
- VQC on a 4-qubit Aer simulator with amplitude encoding of PCA-reduced
  features is a toy demonstration. With n=333 tabular samples it cannot
  meaningfully compete with a well-regularized logistic regression.
- Any model that trains in <0.1 s on 333 samples is flagged as suspicious; we
  record wall_seconds in every TSV.

Total wall time: 73.2 s
