# BioDarwin-RL General Framework Prompt

## Purpose

Use this prompt to ask an AI research agent to design a **general-purpose biomedical algorithm discovery framework**.

This is not a prompt for making a slide deck.  
This is the technical concept document that another agent can read first, then use to build papers, reports, code plans, or presentations.

## Copy-Paste Prompt

You are an expert biomedical AI strategist, AutoML researcher, neural architecture search researcher, reinforcement learning researcher, and translational bioinformatics architect.

Your task is to design a reusable framework called:

**BioDarwin-RL: Knowledge-Graph-Guided Evolutionary Algorithm Discovery for Biomedical AI**

The framework should be general-purpose. It must not be limited to one disease, one dataset, one model, or one biomedical field. It should work as a reusable method for developing new algorithms in domains such as:

- neoantigen immunogenicity prediction
- TCR-pMHC binding prediction
- antibody developability or immunogenicity prediction
- protein engineering
- pathology image-to-molecular biomarker prediction
- spatial transcriptomics prediction from histology
- drug response prediction
- biomarker discovery
- clinical risk stratification
- multi-omics patient selection

## Core Idea

Do not build a simple ensemble.

The goal is to build an **algorithm discovery engine**.

The engine should:

1. Search public AI and biomedical literature.
2. Extract useful method ideas from each paper.
3. Convert those ideas into reusable method primitives.
4. Store the primitives in a knowledge graph.
5. Let a genetic algorithm evolve candidate algorithm designs from the graph.
6. Let a reinforcement learning controller learn how to search the algorithm space more effectively.
7. Validate discovered algorithms under strict biomedical, leakage-control, and prospective-validation constraints.

The key framing:

> BioDarwin-RL converts the scientific literature into a searchable method graph. Genetic algorithms evolve candidate biomedical AI pipelines from the graph, while reinforcement learning learns how to navigate the search space. The result is not a hand-made ensemble, but an automated algorithm discovery engine constrained by biological validity, leakage control, and prospective validation.

## What To Mine From The Literature

Search both AI methodology papers and domain-specific biomedical papers.

AI methodology categories:

- AutoML
- neural architecture search
- evolutionary algorithms
- genetic algorithms
- genetic programming
- reinforcement learning for architecture search
- reinforcement learning for search-policy control
- AutoML-Zero-style algorithm discovery
- differentiable NAS
- Bayesian optimization
- multi-objective optimization
- symbolic distillation
- interpretable algorithm discovery
- foundation-model adaptation
- mixture-of-experts
- uncertainty estimation
- calibration
- active learning

Biomedical method categories:

- state-of-the-art predictive models in the target domain
- benchmark datasets
- preprocessing pipelines
- modality integration strategies
- biological feature engineering
- wetlab validation methods
- clinical validation methods
- known failure modes
- leakage risks
- assay constraints
- interpretability methods
- claim-safety practices

## Method Card Schema

For every paper, model, database, or pipeline, create a method card.

Each method card must include:

- `method_name`
- `year`
- `venue`
- `task`
- `input_modalities`
- `data_sources`
- `preprocessing_steps`
- `feature_engineering`
- `model_architecture`
- `training_objective`
- `validation_design`
- `reported_strengths`
- `known_weaknesses`
- `failure_modes`
- `leakage_risks`
- `biological_assumptions`
- `reusable_primitives`
- `possible_GA_genes`
- `possible_RL_actions`
- `required_guardrails`
- `citation_url`

Do not copy proprietary methods or non-public implementation details. Use only public papers, public code, public documentation, and explicitly available descriptions.

## Primitive Library

Convert method cards into a primitive library.

Primitive categories:

### Data Preprocessing Primitives

Examples:

- sample filtering
- label cleaning
- confidence-weighted labels
- replicate-aware labels
- batch correction
- missingness handling
- sequence normalization
- allele normalization
- image tiling
- expression normalization
- patient-level aggregation
- assay metadata harmonization

### Split And Leakage-Control Primitives

Examples:

- patient holdout
- cohort holdout
- site holdout
- HLA holdout
- sequence-cluster holdout
- source-dataset holdout
- temporal holdout
- leakage blacklist
- duplicate removal
- training-overlap audit

### Modality Primitives

Examples:

- DNA variant
- RNA expression
- protein sequence
- peptide sequence
- HLA allele
- TCR sequence
- antibody sequence
- protein structure
- pMHC structure
- pathology image
- spatial transcriptomics
- clinical covariates
- assay metadata
- drug structure
- cell-line dependency

### Encoder Primitives

Examples:

- tabular gradient boosting
- logistic regression
- random forest
- CNN
- RNN/LSTM
- transformer
- protein language model embedding
- graph neural network
- structure encoder
- image foundation model
- multimodal transformer
- late-fusion MLP

### Fusion Primitives

Examples:

- weighted sum
- rank aggregation
- stacking
- Bayesian model averaging
- mixture-of-experts
- gating network
- attention fusion
- Pareto-front selection
- uncertainty-weighted fusion
- claim-safe cap

### Loss And Optimization Primitives

Examples:

- binary classification loss
- regression loss
- pairwise ranking loss
- top-k surrogate loss
- focal loss
- confidence-weighted loss
- calibration loss
- contrastive loss
- self-supervised pretraining objective
- domain adaptation loss

### Biological Constraint Primitives

Examples:

- expression gate
- presentation gate
- clonality gate
- mutation foreignness gate
- tissue specificity gate
- essentiality filter
- toxicity filter
- pathway plausibility
- structure stability
- binding specificity
- immune escape risk

### Validation Primitives

Examples:

- nested cross-validation
- source-level validation
- external validation
- prospective validation
- wetlab validation
- clinical validation
- ablation study
- calibration curve
- decision-curve analysis
- negative-control experiment
- permutation test

### Interpretability Primitives

Examples:

- feature importance
- SHAP
- counterfactual examples
- attention inspection
- saliency map
- motif extraction
- rule distillation
- symbolic simplification
- failure-mode decomposition

## Knowledge Graph Design

Build a method knowledge graph.

Node types:

- `paper`
- `dataset`
- `method`
- `primitive`
- `modality`
- `preprocessing_step`
- `encoder`
- `loss_function`
- `fusion_rule`
- `biological_constraint`
- `validation_rule`
- `failure_mode`
- `guardrail`
- `assay`
- `clinical_endpoint`
- `algorithm_genome`
- `RL_action`
- `reward_component`

Edge types:

- `uses`
- `requires`
- `contributes`
- `improves`
- `guards_against`
- `validated_by`
- `fails_under`
- `incompatible_with`
- `feeds_into`
- `can_mutate_into`
- `supports_claim`
- `blocks_claim`
- `selected_by_genome`
- `chosen_by_RL_controller`
- `optimized_by_reward`

The knowledge graph should allow the system to answer:

- Which primitives are available for this biomedical task?
- Which modalities are required for each primitive?
- Which primitives are compatible?
- Which primitives create leakage risk?
- Which validation designs are required?
- Which biological claims can be supported?
- Which assay or clinical endpoint is needed to validate the claim?

## GA Genome Design

Represent each candidate algorithm as a genome.

The genome should include:

- selected input modalities
- preprocessing policy
- split policy
- leakage-control policy
- label-confidence policy
- feature extraction modules
- encoder modules
- fusion module
- loss function
- calibration method
- uncertainty method
- biological constraint gates
- interpretability method
- validation protocol
- deployment rule
- assay or clinical decision rule

Example genome:

```yaml
algorithm_genome:
  modalities:
    - sequence
    - expression
    - structure
  preprocessing:
    label_policy: confidence_weighted
    split_policy: patient_and_source_holdout
    normalization: domain_specific
  encoders:
    sequence_encoder: protein_language_model
    tabular_encoder: gradient_boosting
    structure_encoder: graph_neural_network
  fusion:
    type: mixture_of_experts
    gate: modality_missingness_aware
  loss:
    primary: pairwise_ranking_loss
    auxiliary: calibration_loss
  biological_constraints:
    - expression_gate
    - specificity_gate
    - toxicity_filter
  validation:
    - nested_cv
    - external_holdout
    - prospective_locked_test
  interpretability:
    - feature_importance
    - failure_mode_decomposition
```

## GA Operations

Define genetic operations over algorithm genomes.

Mutation operations:

- add a modality
- remove a modality
- swap encoder
- change preprocessing policy
- change split policy
- change fusion rule
- add biological gate
- remove redundant primitive
- change loss function
- add calibration module
- add uncertainty module
- add interpretability module
- change wetlab or clinical decision rule

Crossover operations:

- combine encoder stack from one genome with validation policy from another
- combine preprocessing policy from one genome with fusion strategy from another
- combine biological constraints from one genome with model architecture from another

Repair operations:

- remove invalid modality combinations
- enforce required validation rules
- enforce leakage-control rules
- enforce missingness policy
- enforce biological plausibility

Selection operations:

- tournament selection
- aging evolution
- Pareto-front selection
- diversity-preserving selection
- novelty-preserving selection

Fitness penalties:

- excessive complexity
- high leakage risk
- missing external validation
- overreliance on one feature
- poor calibration
- poor interpretability
- poor prospective-test readiness

## RL Controller Design

The RL agent should not directly predict the biomedical label.

The RL agent controls the search process.

RL state:

- current best validation score
- source-balanced performance
- leakage-resistant performance
- modality missingness
- dataset size
- HLA/site/patient/source diversity
- biological plausibility score
- calibration score
- interpretability score
- disagreement with existing methods
- search diversity
- compute budget remaining
- wetlab or clinical budget remaining

RL actions:

- choose which mutation operator to apply
- choose which primitive category to explore
- choose whether to add or remove modality
- choose whether to increase or reduce model complexity
- choose whether to prioritize performance, novelty, or interpretability
- choose reward weighting for the next search phase
- choose validation split
- choose whether to trigger ablation
- choose whether to freeze a candidate algorithm
- choose assay or clinical validation allocation

RL reward:

- improvement in source-balanced metric
- improvement in leakage-resistant metric
- improvement in top-k practical utility
- improvement over external comparators
- improved calibration
- improved biological plausibility
- improved interpretability
- reduced complexity
- increased novelty
- successful prospective or wetlab validation

Recommended starting point:

- start with contextual bandit or Thompson sampling for search-action selection
- upgrade to PPO, MCTS, or model-based RL only after enough search traces exist
- always compare against random search and fixed heuristic baselines

## Multi-Objective Reward Function

Use a multi-objective reward instead of one metric.

Example reward:

```text
Reward =
  0.25 * source_balanced_AUPRC
+ 0.20 * leakage_resistant_validation_score
+ 0.15 * top_k_experimental_or_clinical_utility
+ 0.10 * improvement_over_external_comparators
+ 0.10 * biological_plausibility
+ 0.08 * interpretability
+ 0.05 * novelty_vs_existing_algorithms
+ 0.05 * calibration_or_uncertainty_quality
+ 0.02 * computational_efficiency
```

The exact weights can change by domain, but the reward must always include:

- a performance metric
- a generalization metric
- a leakage-resistant metric
- a practical utility metric
- a biological plausibility metric
- a complexity or deployability penalty

## Guardrails

The framework must enforce strict guardrails.

Required guardrails:

- use public papers/code/data only
- do not copy proprietary algorithms
- no patient leakage
- no source leakage
- no duplicate sample leakage
- no sequence-cluster leakage when sequence similarity matters
- no single-feature proxy domination
- no post-hoc threshold tuning on final validation
- no final validation reuse during search
- random search baseline required
- fixed ensemble baseline required
- ablation study required
- modality-missingness audit required
- calibration audit required
- external validation required when possible
- prospective or wetlab validation must be locked before results are seen

## Validation Protocol

The framework must distinguish four evidence levels:

### Level 1: Retrospective Internal Search

Used for architecture discovery only.  
No strong biological or clinical claim.

### Level 2: Leakage-Resistant External Validation

Uses held-out sources, patients, sites, alleles, sequence clusters, or cohorts.  
Supports generalization claims if properly controlled.

### Level 3: Locked Prospective Or Wetlab Validation

Algorithm and thresholds are frozen before results are known.  
Supports stronger translational claims.

### Level 4: Clinical Utility Validation

Shows improved decision-making, outcomes, assay yield, or patient selection.  
Required for clinical deployment claims.

## Required Outputs

Produce the following outputs:

1. one-paragraph concept summary
2. detailed framework description
3. method card table
4. primitive library table
5. knowledge graph schema
6. GA genome specification
7. GA operator specification
8. RL controller specification
9. multi-objective reward function
10. guardrail table
11. validation protocol
12. ablation plan
13. risk table
14. claim-safe wording
15. implementation roadmap

Do not produce a PPT slide deck unless separately requested.

## Claim-Safe Wording

Use wording like:

> BioDarwin-RL is a literature-mined algorithm discovery framework that converts public biomedical and AI methods into a knowledge graph of reusable primitives. A genetic algorithm evolves candidate biomedical AI pipelines from this graph, while a reinforcement learning controller learns how to navigate the search space under biological, leakage-control, and validation constraints.

Avoid wording like:

- "proves a new clinical algorithm"
- "guarantees better performance"
- "automatically discovers the true biology"
- "beats all existing methods"
- "requires no expert validation"

Preferred claim:

> The framework generates candidate algorithms and prioritizes them for rigorous external, prospective, wetlab, or clinical validation.

## Final Instruction

Think like a rigorous biomedical AI methods paper reviewer.

The framework should be ambitious, but every claim must be tied to:

- public evidence
- explicit primitives
- leakage-aware validation
- biological plausibility
- ablation
- prospective or experimental confirmation

The final answer should read like a serious technical design document for a reusable biomedical AI algorithm-discovery platform.
