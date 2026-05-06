# Paper 9 Drug Directionality Audit

Local directionality rules:
- AUC matrices: higher raw AUC generally indicates less sensitivity/resistance. Paper 9 Sprint 2 uses direction-normalized sensitivity for association summaries.
- LFC matrices: lower LFC indicates stronger depletion/sensitivity. Direction must not be inferred from treatment metadata alone.
- Metadata-only hits are not interpretable as drug response.

Local result:
- GDSC1 BPTES is available for GLS, but the lineage-sensitivity association is weak/discrepant in the processed summary.
- PRISM 24Q2 treatment metadata did not yield processed SLC1A5, GLUD1, or GLS inhibitor matches in the prior Paper 9 pass.
- SLC1A5/GLUD1 pharmacologic validation remains wet-lab feasibility, not public drug-response support.
