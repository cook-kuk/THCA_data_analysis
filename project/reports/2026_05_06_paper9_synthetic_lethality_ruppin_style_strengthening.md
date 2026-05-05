# Paper 9 Synthetic Lethality Ruppin-Style Strengthening Addendum

Date: 2026-05-06  
Scope: design-only addendum; no execution; no new data acquisition

## 1. Executive Verdict

Paper 9 should be treated as a future therapeutic-vulnerability paper that follows Paper 1.

Paper 1 defines the lineage-silenced thyroid cancer state. Paper 9 asks whether that state creates targetable dependencies that are not captured by canonical driver mutations such as BRAF, RAS, or TERT.

No execution should occur now. No DepMap, PRISM, GDSC, CTRP, CCLE, or other new data should be downloaded for this addendum.

## 2. Ruppin-Style Conceptual Template

Paper 9 can use a Ruppin-style conceptual template without copying text, figures, or claims.

Relevant conceptual components:

- **ISLE-like synthetic lethality**: infer candidate synthetic-lethal relationships by connecting genetic or state-level lesions to selective vulnerabilities.
- **SELECT-like transcriptome-guided therapy prediction**: use tumor transcriptome state to prioritize therapies or dependencies beyond directly actionable mutations.
- **Tumor transcriptome beyond actionable mutations**: treat the lineage-silenced expression state as a therapeutic vulnerability context, not merely as a diagnostic marker.

The adaptation must be thyroid-specific and state-specific. It should not imply that ISLE or SELECT results have already been reproduced in thyroid cancer.

## 3. Thyroid-Specific Adaptation

Input lesion/state:

- RAI-lineage silenced state
- DM1-like dedifferentiated state
- low thyroid-lineage differentiation / compact 8-gene RAI-lineage readout

Target:

- synthetic-lethal or synthetic-rescue dependencies of lineage-silenced thyroid cancer
- druggable vulnerabilities that emerge from the loss of thyroid differentiation programs

The central design principle is to move from mutation-first oncology to state-first vulnerability discovery.

## 4. Main Hypothesis

Lineage-silenced thyroid cancers acquire targetable dependencies that are not captured by BRAF/RAS/TERT status.

These dependencies may include metabolic stress, inflammatory signaling, epigenetic repression, DNA-damage response, kinase rewiring, and lineage-state-associated drug sensitivity.

## 5. Candidate Target Classes

Candidate classes for future evaluation:

- NAMPT / NAD salvage
- STAT3-JAK / IL6 axis
- DNMT / epigenetic repression
- LYN-SRC-FYN
- KCNN4
- ATR-CHEK2
- MYC-GLS-LDHA
- OSMR/IL6R
- TROP2 as non-SL ADC vulnerability only

TROP2 should remain demoted from title/main claim and must not be framed as a synthetic-lethal target unless future evidence directly supports that claim.

## 6. Data Needed Later

Future execution would require:

- DepMap CRISPR Chronos
- CCLE expression
- PRISM/GDSC/CTRP drug response
- TCGA-THCA expression/mutation/clinical
- GSE76039 + GPL570 external cohorts
- thyroid cell lines/organoids

These are future requirements only. This addendum does not authorize downloading or analyzing these resources now.

## 7. Analysis Plan

Future execution plan:

1. Define DM1-like lineage-silenced state.
2. Map candidate target expression and dependency in thyroid-relevant models.
3. Infer synthetic-lethal and synthetic-rescue pairs.
4. Rank druggable vulnerabilities by state specificity, dependency strength, and therapeutic tractability.
5. Validate in thyroid cell lines and organoids.

Key design constraints:

- separate state-associated expression from functional dependency
- separate drug sensitivity from synthetic lethality
- distinguish thyroid-specific vulnerabilities from pan-cancer proliferation artifacts
- preserve Paper 1's role as state definition rather than therapy claim

## 8. Figure Plan

### F1. Lineage-Silenced State

Show Paper 1-derived state definition and the DM1-like / RAI-lineage-silenced input used for vulnerability discovery.

### F2. Dependency Landscape

Display candidate target dependencies across thyroid-relevant cell lines or organoids, benchmarked against lineage-silencing status.

### F3. SL/SR Network

Map inferred synthetic-lethal and synthetic-rescue relationships connected to lineage silencing, candidate targets, and druggable pathways.

### F4. Druggability Map

Prioritize targets by drug availability, assayability, selectivity, and translational feasibility.

### F5. Validation in Cell Lines / Drug Screens

Show future perturbation or drug-screen validation, clearly separated from computational nomination.

### F6. Patient-Level Vulnerability Report

Prototype a research-use report that summarizes lineage state, candidate vulnerabilities, and evidence tiers without making treatment recommendations.

## 9. Business/IP Angle

Potential product concept:

- Thyroid synthetic-lethality vulnerability engine
- research-use report first
- state-first target nomination for lineage-silenced thyroid cancer
- future trial-matching layer only after validation

Boundary:

- no clinical treatment recommendation yet
- no patient-level therapy selection
- no clinical utility claim

## 10. What Not To Claim

Do not claim:

- proven drug response
- patient treatment selection
- TROP2 synthetic lethality
- clinical utility
- validated therapy prediction
- treatment recommendation

Allowed language:

- "future vulnerability discovery framework"
- "candidate dependency nomination"
- "state-guided therapeutic hypothesis"
- "research-use prioritization"

## 11. Why This Is Paper 9, Not Paper 1/3

Paper 1 defines the driver-orthogonal thyroid-lineage differentiation axis and the compact 8-gene RAI-lineage readout.

Paper 3 defines immune vulnerability, antigen-presentation architecture, TLS biology, and ICI readiness.

Paper 9 defines cell-intrinsic targetable dependencies, synthetic-lethal or synthetic-rescue hypotheses, druggability, and future functional validation.

The three papers are connected but should not collapse into one claim space.

## 12. 12-Week Future Execution Plan After Marathon

### Weeks 1-2: State Lock

Freeze the DM1-like / RAI-lineage-silenced state definition inherited from Paper 1. Define inclusion rules for tumors, cell lines, and organoids.

### Weeks 3-4: Data Acquisition and Harmonization

Acquire approved dependency, expression, drug-response, and thyroid cohort resources. Harmonize identifiers, lineage-state scores, mutations, and target classes.

### Weeks 5-6: Dependency Mapping

Map lineage-state associations with CRISPR dependency and candidate target expression. Flag proliferation-linked artifacts and low-confidence models.

### Weeks 7-8: SL/SR Inference

Infer candidate synthetic-lethal and synthetic-rescue pairs. Prioritize state-specific dependencies over broad cancer-essential targets.

### Weeks 9-10: Drug Response Integration

Integrate PRISM/GDSC/CTRP-style drug-response evidence where available. Separate expression, dependency, and drug sensitivity evidence tiers.

### Weeks 11-12: Validation Blueprint and Manuscript Skeleton

Design cell-line/organoid validation assays, finalize figure shells, and draft the Paper 9 manuscript architecture as a future therapeutic-vulnerability paper.

