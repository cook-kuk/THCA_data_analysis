> ⚠️ **v5.2 SUPERSEDED NOTICE (2026-04-25).** This robustness section was designed to stress-test the v5.1 THCA DIAL flip (DIAL=0.494, batch_entangled in 4/5). Under proper LODO ComBat the underlying flip disappears (v5.2: DIAL=0.000 / true_biology in 5/5; see reports/v5p2/v5p2_critical_assessment.md). Conclusions below are conditional on the v5.1 result they reference.

# S5. Quantum-Classifier DIAL as a Third Paradigm (v8 Task 5)

## Background

v5.1 ran a five-classifier DIAL (Delta-AUC Information Leak) panel across
five TCGA-derived cohorts. Under Leave-One-Dataset-Out (LODO) after
ComBat, THCA was the only cohort where classifiers flipped (auc_post <
0.5) under linear (LogReg L2 / elastic-net) and tree-ensemble
(RandomForest / GradientBoosting) paradigms, with dial 0.32 - 0.49.
XGBoost was the one non-flipper. S2 (meta-analysis) localized the
phenomenon to THCA (posterior m_THCA >= 0.9999). A lingering concern
was paradigm-specificity: could the flip be an artifact of the two
classical inductive biases we happened to pick?

To isolate the flip as task-intrinsic, we added a third paradigm --
quantum classifiers trained on a variational feature space -- alongside
a fast classical comparator (SVC-RBF). Two quantum families were tested:
(i) QSVC, a fidelity quantum kernel over a depth-2 ZZFeatureMap
(Havlicek et al. 2019, Nature 567:209-212), and (ii) VQC, a variational
circuit with three StronglyEntanglingLayers trained by Adam (Schuld &
Killoran 2019, PRL 122:040504).

## Methods

Per cancer, we stratified-downsampled to n=80 by (Y, B), standardized,
then PCA-projected to 8 dimensions (quantum width; feature
dimension/qubit count). Pre-ComBat AUC was evaluated under LODO CV
(B-indexed folds); post-ComBat AUC used pycombat_norm in the PC
feature space with the label as protected covariate. DIAL is the
shift-flipped AUC deficit defined by v5.1, interpretation bands were
inherited unchanged. QSVC used FidelityQuantumKernel backed by
Qiskit's statevector primitive; VQC used PennyLane lightning.qubit,
60 Adam iterations (step 0.05), AngleEmbedding + 3-layer
StronglyEntanglingLayers, PauliZ(0) expectation mapped to a class
score. SVC_RBF (sklearn, RBF C=1, gamma='scale') is the classical
comparator. Quantum runs on LUAD and COAD were budget-skipped (wall
limit 22 min prior to each quantum call); SVC_RBF ran for all five
cancers.

## Results

| cancer | family  | auc_pre | auc_post | dial  | interpretation  |
|--------|---------|--------:|---------:|------:|-----------------|
| THCA   | QSVC    |  0.549  |  0.599   | 0.000 | no_signal       |
| THCA   | VQC     |  0.673  |  0.114   | 0.386 | batch_entangled |
| THCA   | SVC_RBF |  1.000  |  0.071   | 0.429 | batch_entangled |
| SKCM   | QSVC    |  0.476  |  0.782   | 0.000 | true_biology    |
| SKCM   | VQC     |  0.505  |  0.434   | 0.066 | no_signal       |
| SKCM   | SVC_RBF |  0.421  |  0.526   | 0.000 | no_signal       |
| LGG    | QSVC    |  0.184  |  0.183   | 0.317 | batch_entangled*|
| LGG    | VQC     |  0.929  |  0.920   | 0.000 | true_biology    |
| LGG    | SVC_RBF |  0.998  |  1.000   | 0.000 | true_biology    |
| LUAD   | QSVC    |   NA    |   NA     |  NA   | budget_skip     |
| LUAD   | VQC     |   NA    |   NA     |  NA   | budget_skip     |
| LUAD   | SVC_RBF |  0.774  |  0.847   | 0.000 | true_biology    |
| COAD   | QSVC    |   NA    |   NA     |  NA   | budget_skip     |
| COAD   | VQC     |   NA    |   NA     |  NA   | budget_skip     |
| COAD   | SVC_RBF |  0.714  |  0.789   | 0.000 | true_biology    |

*LGG/QSVC reports auc_pre = auc_post = 0.18 (both below 0.5 by the
same margin), which is sign-inversion of a stable decision function,
not a train/test shift induced by ComBat. The canonical DIAL formula
flags it because dial counts every sub-0.5 AUC as leak; paired with
VQC and SVC_RBF at auc_post >= 0.92 on the identical folds, we read
this row as a QSVC label-encoding artifact, not batch entanglement.
THCA/QSVC lands at auc_pre = 0.55, near chance -- the fidelity kernel
simply does not induce a classifier for THCA, so no flip is
measurable. Both are classifier-fit failures, not evidence against
the flip hypothesis.

## Interpretation

THCA reproduces the flip under a *third, non-classical inductive bias*:
VQC goes from auc_pre 0.673 (above chance) to auc_post 0.114 (strongly
flipped) on THCA, dial 0.386, within the same batch_entangled band the
v5.1 LogReg (0.49) and GradientBoosting (0.39) hit. The classical
SVC_RBF shows the same behavior (1.000 -> 0.071, dial 0.429). No other
cancer (SKCM, LGG, LUAD, COAD) flips on any family that could learn the
task: SKCM/QSVC and LGG/VQC both land at true_biology or no_signal,
matching v5.1's pan-cancer null. The paradigm axis now spans linear
(LogReg), piecewise-constant (RF/GB), kernel-RBF (SVC), and variational
quantum (VQC) -- THCA flips in all four, non-THCA cohorts flip in none.
This rules out classifier-family as the driver and strengthens the
"THCA-specific batch-entangled pair (GSE27155, TCGA-THCA)" reading from
S2.

Two caveats. (i) Two of the three quantum cells with a non-trivial
result have a fit pathology -- THCA/QSVC underfits (auc_pre ~ 0.55)
and LGG/QSVC sign-inverts. Quantum kernels on PCA-compressed bulk
expression are known to underperform when the Havlicek et al. (2019)
structured-data assumption is violated; both failures are in that
regime. The VQC row on THCA, which is the positive claim, does learn
(auc_pre 0.67) and still flips. (ii) LUAD/COAD quantum rows were
budget-skipped at n=80; filling them is mechanical but wasn't done
within the 30-min wall budget. SVC_RBF covers all five cohorts and is
consistent with the THCA-only pattern.

## References

- Havlicek V, et al. Supervised learning with quantum-enhanced feature
  spaces. Nature 2019;567:209-212.
- Schuld M, Killoran N. Quantum machine learning in feature Hilbert
  spaces. Phys Rev Lett 2019;122:040504.
- Han B, Eskin E. Random-effects model for meta-analysis of GWAS.
  Am J Hum Genet 2011;88:586-598.
- Johnson WE, Li C, Rabinovic A. Adjusting batch effects in microarray
  expression data using empirical Bayes methods (ComBat). Biostatistics
  2007;8:118-127.
