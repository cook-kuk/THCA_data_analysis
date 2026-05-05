# Paper 9 Sprint 2 Metabolic Vulnerability Plan

## Reframed Goal
Paper 9 is framed as a metabolic vulnerability prioritization map for lineage-silenced thyroid cancer, not as validated synthetic lethality.

## Current Strongest Claim
Lineage-silenced thyroid cancer may exhibit glutamine-axis vulnerability candidates requiring experimental validation.

## Sprint 2 Modules
- Glutamine-axis genes: GLS, GLUD1, GLUD2, SLC1A5, SLC7A5, ASNS, GOT1, GOT2, GPT2.
- MYC glutamine-addiction module.
- LDHA/glycolysis comparator.
- OXPHOS comparator.
- NAMPT/NAD salvage comparator.

## Analysis Tasks
- Use existing DepMap/CCLE CRISPR and expression matrices only.
- Score individual genes and modules for dependency strength and expression.
- Estimate pan-cancer and cancer-type-centered associations with the lineage-silenced score.
- Audit common-essential, nonessential, tissue-context, MYC/LDHA/NAMPT artifact risks.
- Summarize processed PRISM/GDSC/CTRP drug evidence when metadata matches are available.

## Current Prioritization Snapshot
- GLUD1: PRIORITIZE_FOR_WETLAB (MEDIUM); claim=allowed_candidate_hypothesis.
- SLC1A5: PRIORITIZE_FOR_WETLAB (MEDIUM); claim=allowed_candidate_hypothesis.
- GLS: KEEP_AS_CANDIDATE (MEDIUM_MINUS); claim=allowed_candidate_hypothesis.
- MYC_glutamine_addiction_module: KEEP_AS_AXIS_OR_SECONDARY (HYPOTHESIS_SUPPORT_ARTIFACT_GUARDED); claim=weak_artifact_guard_required.
- LDHA_glycolysis_comparator: COMPARATOR_ONLY (LOW_ARTIFACT_GUARDED); claim=weak_comparator_only.
- NAMPT_NAD_salvage_comparator: COMPARATOR_ONLY (LOW_ARTIFACT_GUARDED); claim=weak_comparator_only.

## Decision Boundary
- Allowed: candidate metabolic vulnerability and research-use prioritization.
- Weak: thyroid-specific claims, because strict thyroid CRISPR overlap remains small.
- Forbidden: validated synthetic lethality, patient selection, clinical treatment recommendation.
