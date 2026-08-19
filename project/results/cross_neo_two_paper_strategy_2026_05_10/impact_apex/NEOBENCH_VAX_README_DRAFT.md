# NeoBench-Vax

NeoBench-Vax evaluates whether neoantigen AI systems can produce high-precision vaccine candidate shortlists while controlling public overlap, near-peptide memorization, source/study shift, HLA shortcuts, ambiguous negatives, calibration and OOD risk.

## Tracks
- Clean Comparator: public predictor scores are not training features.
- Product-Assisted: public predictor scores allowed if disclosed.
- Source-Stress: source/study/prevalence shift.
- OOD/Abstention: safe non-decision behavior.

## Primary Metrics
AUPRC, top5/top10 precision, enrichment over prevalence, recall@k and abstention coverage versus precision.
