
# BAR-Neo Method Card

## Algorithm

BAR-Neo is a Benchmark-Adaptive Reliability Neoantigen Ranking layer. It is not a giant deep model. It combines candidate intrinsic features, method rank percentiles, predictor disagreement, overlap/leakage flags, HLA/source support, historical source-heldout behavior, and transparent patient gates.

## Feature Groups

- Candidate intrinsic features: peptide length, hydrophobic/aromatic/charge proxy, cysteine count, glycine/proline count, anchor residue proxies.
- HLA features: HLA gene and supertype one-hot, allele support count, Korean-HLA focus flags.
- Method-derived features: score percentiles, number of methods available, clean internal score, caveated public score, public-vs-internal disagreement.
- Reliability features: exact/near/source/study/patient overlap flags, leakage risk level, source prevalence, source-heldout historical top-k behavior.
- Uncertainty features: predictor disagreement, low method availability, underrepresented HLA, uncertainty-only model support where available.

## Calibration Logic

Method scores are direction-normalized and calibrated with fold-safe isotonic/logistic calibration when feasible. Small methods fall back to rank-percentile calibration and are marked through abstention logic when unstable.

## Abstention Logic

BAR-Neo abstains or lowers confidence for high leakage risk, underrepresented HLA, source-shift risk, high method disagreement, caveated-public-only evidence, MHC-II rows in a Class-I benchmark, failed presentation gates, low-prevalence source collapse, and sparse method coverage.

## Patient-Gated Scoring Logic

`patient_gated_score = barneo_score * disease_gate * presentation_gate * antigen_gate * immune_context_gate * model_confidence_gate`

PAAD high-priority context is resected/MRD/low-burden disease. THCA high-priority research context is ATC, PDTC, progressive radioiodine-refractory DTC, or high-risk recurrence. Unknown patient metadata are marked as uncertainty rather than inferred.

## Limitations

BAR-Neo currently reflects available public/local benchmark labels and sparse patient metadata. It is a research triage and reliability layer, not a clinical treatment selector.
