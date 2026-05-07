# Paper 3 ICI Track-A+ Strengthening Addendum

Date: 2026-05-06  
Scope: design-only addendum; no execution; no new data acquisition

## 1. Executive Verdict

Paper 3 should be strengthened as an **ICI vulnerability / readiness atlas**, not as an ICI response predictor.

The paper's strongest defensible claim is that dedifferentiated or dark-matter thyroid cancer can be immunogenomically stratified into states with different levels of ICI readiness, immune exclusion, antigen-presentation competence, and TLS-associated immune organization.

Track A remains the frozen foundation. Track B remains blocked until the required protected-access and external ICI datasets are available. This addendum does not authorize Track B execution.

## 2. Strengthened Central Hypothesis

Dedifferentiated / dark-matter thyroid cancer is not immunologically uniform. Instead, it separates into clinically relevant immune ecotypes:

- **ICI-ready ecotype**: IFNG/TIS-high, MHC-I/MHC-II competent, cytolytic, and potentially TLS/CXCL13 enriched.
- **Immune-excluded ecotype**: inflammatory or stromal context without effective cytotoxic infiltration.
- **HLA-lost ecotype**: immune pressure or dedifferentiation accompanied by loss of antigen presentation.
- **TLS-rich ecotype**: CXCL13/TLS/B-cell-associated immune organization that may indicate immune engagement even without direct response data.

The paper should argue for **prioritization and biological readiness**, not therapeutic response prediction.

## 3. Track-A+ Modules

### A. Bulk Immune Ecotype

Purpose: define immune-state groups from bulk tumor transcriptomes.

Core readouts:

- IFNG/TIS activation
- Cytolytic activity
- MHC-I and MHC-II presentation
- T-cell, B-cell, macrophage, myeloid, and stromal/exclusion features
- Dedifferentiation-context interaction

Deliverable: thyroid-specific immune ecotype map anchored to lineage-silencing status.

### B. scRNA Immune-State Atlas

Purpose: resolve whether bulk immune signals reflect true immune-cell states, malignant epithelial dedifferentiation, or stromal contamination.

Core readouts:

- CXCL13-positive exhausted/activated CD8 states
- inflammatory TAM states
- M2-like TAM states
- Treg enrichment
- MDSC-like myeloid states
- dedifferentiated epithelial compartment

Deliverable: cell-state-level interpretation of bulk immune readiness and exclusion.

### C. HLA Intactness + Neoantigen Presentability

Purpose: distinguish inflamed tumors that can present antigen from inflamed tumors with impaired presentation.

Core readouts:

- HLA-I expression
- HLA-II expression
- B2M and antigen-processing genes
- HLA LOH or other HLA disruption when protected-access data permit
- mutation-derived neoantigen presentability when WES/BAM/dbGaP gates permit

Deliverable: antigen-presentation architecture layered onto immune ecotypes.

### D. DIAL Audit of Pan-Cancer ICI Signatures

Purpose: test whether established pan-cancer ICI-associated signatures retain their biological directionality in thyroid cancer.

Verdict classes:

- **PASS**: expected immune/readiness direction is preserved in thyroid context.
- **FLIP**: thyroid association reverses relative to expected ICI biology.
- **COLLAPSE**: signature loses dynamic range, specificity, or interpretable association.
- **AMBIGUOUS**: mixed or unstable behavior; not eligible for score inclusion.

Only PASS signatures can enter the thyroid ICI-readiness score.

### E. Integrated ICI-Readiness Score

Purpose: convert module-level evidence into a transparent prioritization score.

The score should be interpretable, modular, and audit-friendly. It should not be presented as a clinical response model.

## 4. New Signature Candidates

Candidate signatures for Track-A+ development:

- IFNG/TIS
- CYTOLYTIC
- MHC1_CORE
- MHC2_CORE
- TLS_CABRITA
- CXCL13_AXIS
- Treg
- M2_TAM
- MDSC_like
- scCXCL13_CD8
- scTAM_INFLAM
- scDEDIFF_EPI

Each candidate should be assigned to one of three roles:

- readiness-positive signal
- resistance/exclusion signal
- context/modifier signal

## 5. Strengthened DIAL Audit

The DIAL audit should become a formal gatekeeper for importing pan-cancer ICI signatures into thyroid cancer.

Required signature-level fields:

- source signature
- intended biological meaning
- expected ICI-associated direction
- thyroid association with lineage silencing
- thyroid association with immune ecotype
- PASS / FLIP / COLLAPSE / AMBIGUOUS verdict
- eligibility for integrated score

Rule: **Only PASS signatures enter the thyroid readiness score.**

Signatures marked FLIP, COLLAPSE, or AMBIGUOUS may be discussed as biology, but they cannot contribute to the final readiness score.

## 6. Integrated Score Proposal

Proposed conceptual score:

```text
ICI_READINESS =
  IFNG/TIS
  + MHC1
  + MHC2
  + TLS/CXCL13
  + neoantigen presentability
  - HLA_LOH
  - M2_TAM
  - MDSC
  - Treg
  + dedifferentiation-context interaction
```

Interpretation:

- high score: ICI-ready or immunogenomically prioritized tumor state
- intermediate score: immune-active but incomplete readiness
- low score: immune-cold, HLA-impaired, or suppressive/excluded state

The dedifferentiation-context interaction is essential. It prevents the score from becoming a generic pan-cancer immune score and anchors the model to thyroid dark-matter biology.

## 7. Claim Guard

Allowed claims:

- ICI vulnerability
- ICI readiness
- immunogenomic prioritization
- immune ecotype atlas
- antigen-presentation-aware thyroid immune stratification

Forbidden claims:

- ICI response predictor
- patient treatment recommendation
- thyroid ICI response validation
- clinically actionable treatment selection
- prospective clinical utility

Preferred language:

- "prioritizes tumors for future ICI-focused evaluation"
- "defines immunogenomic readiness states"
- "identifies antigen-presentation-aware immune ecotypes"

Avoid:

- "predicts ICI response"
- "selects patients for ICI"
- "validates ICI benefit in thyroid cancer"

## 8. Figure Upgrades

### F1. Cohort Landscape

Show available Track A cohorts, lineage-silencing states, immune-score distributions, and access gates for unavailable Track B assets.

### F2. Bulk Immune Ecotypes

Cluster tumors by IFNG/TIS, cytolytic, MHC, TLS/CXCL13, suppressive myeloid, Treg, and dedifferentiation features.

### F3. scRNA Immune States

Map immune and epithelial cell states, emphasizing scCXCL13_CD8, scTAM_INFLAM, M2-like TAM, Treg, MDSC-like, and scDEDIFF_EPI programs.

### F4. HLA/Neoantigen Architecture

Layer HLA intactness, antigen-processing expression, HLA LOH where available, and neoantigen presentability where permitted.

### F5. DIAL Audit

Display pan-cancer ICI signatures with PASS / FLIP / COLLAPSE / AMBIGUOUS verdicts and eligibility for thyroid score inclusion.

### F6. Integrated Readiness Atlas

Combine immune ecotype, HLA/neoantigen status, DIAL-pass signatures, and dedifferentiation interaction into ICI-readiness groups.

## 9. Dataset Access Gates

Track-A+ execution requires explicit future access decisions. No access or execution is authorized by this addendum.

Required gates:

- TCGA-THCA WES/BAM/dbGaP access for HLA LOH and neoantigen presentability.
- scRNA cohort availability and permission for immune-state atlas construction.
- At least 3 pan-cancer ICI cohorts for DIAL signature auditing.

Until these gates are cleared, Paper 3 remains a Track A frozen concept with Track-A+ design only.

## 10. Why This Is Paper 3, Not Paper 1/2/4

Paper 1 defines the driver-orthogonal thyroid-lineage differentiation axis and establishes the 8-gene compact RAI-lineage readout. It should not be expanded into an ICI atlas.

Paper 2 and Paper 4 remain untouched by this addendum.

Paper 3 is the appropriate home because it asks how lineage-silenced thyroid tumors differ in immune vulnerability, antigen presentation, TLS organization, and ICI readiness. It builds on Paper 1's state definition but does not modify Paper 1's manuscript claims.

