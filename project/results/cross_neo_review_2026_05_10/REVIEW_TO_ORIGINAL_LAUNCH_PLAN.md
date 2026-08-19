# Two-Paper Launch Plan

## North Star

Do not publish CROSS-Neo as "one more neoantigen predictor." Publish the Review first to define the field's missing standard, then publish CROSS-Neo as the first implementation of that standard.

## Paper 1: Review / Perspective

### Title
From neoantigen discovery to vaccine-ready prioritization

### Main Job
Make the field accept the premise:
- vaccine candidate selection is a top-k ranking problem;
- AUROC is insufficient;
- public benchmark overlap must be explicitly audited;
- HLA/source/study shifts are not edge cases;
- incomplete negatives require PU-aware interpretation;
- OOD-aware abstention is clinically useful, not a weakness.

### Winning Contribution
A practical evaluation contract for neoantigen prioritization.

### Do Not Do
- Do not center CROSS-Neo.
- Do not attack public predictors.
- Do not claim external validation.
- Do not include quantum claims.

## Paper 2: Original CROSS-Neo

### Title
CROSS-Neo: contamination-controlled pan-allele neoantigen prioritization for vaccine candidate triage

### Main Job
Show one working implementation of the Review's evaluation contract:
- pan-allele counterfactual mutant-WT branch;
- hard-decoy and PU-aware pressure;
- explicit public-comparator separation;
- train-only retrieval and leakage audit;
- top-k, AUPRC and enrichment as primary outcomes;
- OOD/source-aware abstention;
- business/product track separated from clean scientific claim.

### Clean Claim
Internal locked evidence that CROSS-Neo improves practical top-k prioritization over strong public predictor-assisted baselines under defined split conditions.

### Forbidden Claim
External validation, universal superiority, or quantum advantage.

## Timeline

### Week 0
- Finish Review synopsis package.
- Send Nature Reviews Clinical Oncology presubmission enquiry.
- Send Nature Cancer Perspective enquiry if no conflict with journal policy.
- Build graphical abstract and figure mockups.

### Week 1
- Expand Review skeleton to 3,500-4,500 words.
- Fill Table 1 clinical vaccine signal table.
- Fill Table 2 predictor taxonomy.
- Prepare benchmark failure-mode figure.

### Week 2
- If editor interest: tailor manuscript to requested format.
- If no response: prepare Trends in Cancer and JITC versions.
- In parallel freeze CROSS-Neo original paper claim boundary and final benchmark tables.

### Week 3-4
- Submit Review/Perspective.
- Build original CROSS-Neo preprint draft around the evaluation contract.
- Separate product demo white paper from clean manuscript.

## Business Use

The Review creates category leadership:

> We are not merely selling another predictor; we are defining the evaluation standard for vaccine-ready neoantigen prioritization.

The product demo then becomes:

> A working pan-allele implementation of the standard, already outperforming public predictor-only baselines on an internal retrospective triage set.

## Best Possible Outcome

1. Nature Reviews/Nature Cancer engages with the Review topic.
2. CROSS-Neo original is framed as an implementation paper rather than a model paper.
3. Product demo benefits from category authority before the algorithm paper is accepted.
