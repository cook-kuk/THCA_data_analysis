# Grand Challenge Manifesto

Personalized cancer vaccines require a short list of candidates. That list is increasingly shaped by computational models, but the field lacks a shared standard for deciding whether a ranked list is trustworthy.

The central issue is not whether one more model can improve AUROC. The central issue is whether any model can show that its shortlist is not driven by public-data overlap, near-peptide memorization, HLA shortcuts, source-specific assay bias, or ambiguous negative labels.

NEO-PRIOR defines the minimum evidence that vaccine-ready neoantigen AI should report. NeoBench-Vax makes the standard operational through locked splits, public-overlap audits, model cards and leaderboard tracks. CROSS-Neo is the first implementation, evaluated conservatively under this framework.

The clinical promise is not a magic oracle. The practical promise is auditable prioritization: fewer candidates, clearer evidence tags, better top-k enrichment, and explicit abstention when the evidence is unsafe.
