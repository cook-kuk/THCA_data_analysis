# CROSS-Neo Kaggle Winner Playbook

Date: 2026-05-10

Scope: public Kaggle winning/top-solution writeups that are closest to our cancer-vaccine / neoantigen setting: RNA vaccine stability, RNA folding/reactivity, molecular binding, single-cell perturbation, MoA, protein function, and cancer pathology. This is not a literal archive of every Kaggle competition. It is a filtered import list for `CROSS-Neo`, `BAR-Neo`, TCR extension, hard-decoy repair, and patient-gated vaccine triage.

## Bottom line

The useful Kaggle lesson is not "make one bigger model." The recurring winner pattern is:

1. build a validation contract that mimics the hidden/private distribution,
2. use multiple complementary experts,
3. select ensemble weights by OOF quality plus prediction diversity,
4. add leak/overlap kill-switches,
5. use pseudo-labeling or self-supervised pretraining only where it does not contaminate claims,
6. postprocess predictions to obey the biological hierarchy,
7. keep a manual-review/abstention path for rows the model should not claim.

Our current system already has several of these pieces: source/HLA/near-peptide stress splits, BAR-Neo-BMA, BAR-Neo-X abstention, public-overlap audit, hard decoys, TCR diagnostic branch, and patient gates. The highest-yield import is to tighten them into Kaggle-style locked experiments rather than starting a new model family.

## P0 imports

### 1. Ensemble diversity weighting, not just "best mean AUPRC"

Kaggle MoA winner logic: blend strong but low-correlated models; use OOF to optimize weights and avoid redundant models. CROSS-Neo translation:

- Add an ensemble-correlation penalty to BAR-Neo-BMA weights.
- Keep the existing reviewer-safe score, but add `pairwise_oof_corr`, `incremental_topk_gain`, and `stress_min_floor`.
- Downweight methods that are individually strong but collapse in the same failure slice as the current winner.

Why it matters here: `source_balanced_*` models look strong on overall AUPRC but source-heldout deltas are often negative. Diversity-aware weights should favor experts that rescue different slices: HLA stress, low-prevalence TESLA, NEPdb, Korean HLA, TCR-supported rows.

### 2. Validity DAG / conditional caps

CAFA top solutions exploited label hierarchy: child probabilities must be consistent with parent terms. CROSS-Neo should do the same biologically:

`patient_useful_vaccine_candidate <= disease_context_gate <= immune_context_gate <= TCR_recognition_gate <= presentation_gate <= antigen_gate`

Concrete cap examples:

- If `hla_loh` or `hla_expression` is unknown, cap clinical/patient score and route to `research_triage_only`.
- If public overlap is unresolved, allow "caveated comparator" but block "clean claim."
- If presentation is weak, TCR evidence cannot rescue a candidate into a claim-safe vaccine hit.
- If exact TCR evidence exists but pMHC presentation is weak, route to structure/manual review rather than promote.

This is the most direct way to convert Kaggle hierarchy tricks into reviewer-safe translational logic.

### 3. Label-noise audit, not blind relabeling

PANDA cancer pathology winners used OOF prediction disagreement to detect noisy labels. CROSS-Neo translation:

- Generate `label_noise_audit_queue.tsv` from OOF residuals under exact peptide-HLA, near-peptide, HLA-heldout, and source-heldout splits.
- Mark rows as `trusted_positive`, `trusted_negative`, `ambiguous_assay`, `source_conflict`, or `manual_review`.
- Do not auto-relabel in manuscript-facing benchmarks. Train an auxiliary "denoised training variant" only for sensitivity analysis.

Why this matters here: public neoantigen labels mix assay systems, pathogens/cancer contexts, and positive-unlabeled behavior. A PANDA-style audit is safer than pretending the labels are clean.

### 4. Hard-negative miner v2

Our hard-decoy repair is already in the right Kaggle direction. The v2 should split decoys by failure mechanism:

- `anchor_preserved_decoy`: same length and anchor-like residues, scrambled non-anchor positions.
- `near_positive_decoy`: high peptide similarity but different HLA/source label.
- `hla_confuser_decoy`: same peptide, related HLA supertype where presentation should differ.
- `tcr_shuffle_decoy`: same peptide-HLA, shuffled TCR for pMTnet/TEPCAM branch.
- `source_shift_decoy`: CEDAR-like positives against NEPdb/TESLA-like negatives.

Train these as bounded auxiliary signals, not as a standalone "motif predictor." The current report already shows why standalone decoys are brittle.

### 5. Leakage kill-switch as a first-class model feature

Ribonanza/OpenVaccine competitions emphasized true private/future testing and test leakage removal. CROSS-Neo should make overlap state part of every model card:

- exact peptide-HLA overlap,
- near-peptide cluster overlap,
- source overlap,
- HLA supertype overlap,
- public pretrained tool training-corpus overlap,
- TCR-resource overlap.

Allowed behavior:

- Use overlap flags for abstention, model selection, and caveat labels.
- Do not use unresolved public predictor scores as clean features.

### 6. Weighted sampling by evidence quality

Ribonanza winners/top models used signal-quality handling instead of treating every profile equally. CROSS-Neo equivalent:

- sample/weight by assay quality, source reliability, HLA allele support, cancer-context evidence, patient metadata completeness, and label confidence;
- compare hard filtering vs soft weights under source/HLA/near-peptide locked splits;
- report low-quality rows as abstain/manual review rather than hiding them.

## P1 imports

### 7. Self-supervised peptide/pMHC pretraining

BELKA winner pattern: simple transformer plus task-relevant pretraining beat many complex alternatives. CROSS-Neo translation:

- Masked amino-acid modeling on peptide strings and HLA pseudo-sequences.
- Descriptor prediction head: length, anchor class, hydrophobicity profile, charge, instability index, cleavage/TAP-like descriptors.
- Contrastive positive pairs: same peptide across equivalent HLA formatting; same HLA across pseudo-sequence aliases.
- Strict rule: no benchmark labels and no public comparator score labels in pretraining.

This is a practical way to get representation learning without importing public-tool leakage.

### 8. Pairwise attention bias from pMHC structure

OpenVaccine/Ribonanza winners repeatedly used 1D sequence plus 2D pairwise matrices as attention bias: BPP, distance, adjacency, contact-style features. CROSS-Neo equivalent:

- peptide position x HLA pocket/contact matrix,
- peptide position x peptide position distance/anchor matrix,
- TCR CDR3 position x peptide position contact prior when exact or modelable TCR exists,
- MD-derived stability/contact occupancy as a diagnostic attention-bias feature, not a primary label.

Short-term: export these as low-dimensional tabular features for RF/LR/MoE.
Long-term: sequence/pair two-track model for pMHC/TCR-pMHC.

### 9. Length-generalization guard

RNA competitions had train/test length mismatch, and winners used relative/dynamic position encodings or length-aware validation. CROSS-Neo translation:

- separate peptide length buckets: 8, 9, 10, 11, 12+;
- force length-heldout and long-peptide stress splits;
- avoid absolute positional features that overfit 9-mers;
- evaluate class-II-like/long class-I candidates separately.

### 10. Patient/context embeddings

Open Problems single-cell perturbation winners used cell-type and molecule text embeddings plus target encoding/one-hot features. CROSS-Neo equivalent:

- disease context embedding: THCA, PAAD, MRD/adjuvant/metastatic, fusion/BRAF/KRAS state;
- immune context embedding: TLS/IFNG/cytolytic/HLA expression;
- antigen context embedding: expression, clonality, VAF, driver/passenger;
- missingness bits as real features, not silent zeroes.

This should feed patient-gated score caps, not replace pMHC evidence.

## P2 / optional imports

### 11. DeepInsight-style feature map

MoA winners added CNN diversity by mapping correlated gene/cell features into a 2D image. CROSS-Neo equivalent could map expert scores and evidence axes into a small "candidate evidence image" and train a CNN/TabNet auxiliary model. This is lower priority because reviewer interpretability and small-n source stress matter more than squeezing leaderboard points.

### 12. External "top-3 blend" pseudo-labels

Ribonanza released top-model blends and used pseudo-labels post-competition. CROSS-Neo can use OOF soft labels on unlabeled candidate bags, but only for discovery triage. Do not train benchmark claims on test-set pseudo-labels.

## Explicit no-go list

- No leaderboard probing logic in manuscript-facing benchmarks.
- No test-set pseudo-labels for clean claims.
- No unresolved public predictor score as a clean model feature.
- No automatic label rewrite without a sensitivity branch and audit trail.
- No quantum/headline claim from QK branch; keep it bounded fallback/fusion.
- No TCR-recognition claim from synthetic-decoy public benchmarks alone.

## Files in this package

- `kaggle_winner_tricks_to_cross_neo.tsv`: detailed trick-to-action table.
- `implementation_queue.tsv`: practical P0/P1/P2 experiment queue.
- `MERGED_WITH_SEUNGHO_IDEAS_KR.md`: Korean merged strategy that treats the existing CROSS-Neo/BAR-Neo/TCR/MD/patient-gate ideas as the main architecture and Kaggle tricks as reinforcements.
- `seungho_idea_x_kaggle_matrix.tsv`: row-level crosswalk from existing ideas to Kaggle imports and next actions.
- `source_manifest.tsv`: source URLs and what was extracted from each.
