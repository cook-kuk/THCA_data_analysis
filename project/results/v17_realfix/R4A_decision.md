# R4-A · REAL FIX decision matrix

Generated: 2026-04-27 12:50

## Result class: SCENARIO A (8-gene VALID, framing fix)

| Criterion | Threshold | Result | Pass? |
|-----------|-----------|--------|-------|
| Leak-free CV AUC > 0.85 (any variant) | 0.85 | R1-A 0.962, R1-B 0.925, R1-D 0.957 | YES |
| ΔAUC vs BRAF > +0.05 | +0.05 | +0.111 ~ +0.130 | YES |
| 95% CI excludes 0 (ΔAUC) | exclude 0 | yes for R1-A,B,D | YES |
| External validation (real ground truth) | AUC > 0.7 | GSE76039 0.935 | YES |
| Cluster robustness (ARI vs orig) | > 0.4 | R1-A 0.642, R1-D 0.680 | YES |

→ **Scenario A confirmed**: 8-gene panel is a valid contribution. The original AUC 0.954 was high-but-real; not circular artifact. The ΔAUC +0.132 over BRAF is sustained at +0.111 ~ +0.130 across leak-free cluster definitions. Original cluster definition included 8-gene, but the biology these 8 genes encode is genuine (TIERA67 - 8 still recovers very similar partitioning; ARI 0.642).

## What needs to change in the manuscript

1. **Reframe the 8-gene panel as a parsimonious classifier of the leak-free cluster**, not as predicting an "independent" DM1/DM2. Acknowledge that DM1/DM2 was originally defined from a 67-gene panel including the 8-gene set, and demonstrate that re-definition without the 8-gene set produces concordant clusters (ARI 0.64 - 0.68) with sustained AUC and ΔAUC.

2. **Replace 4-cohort meta-analysis (proxy labels)** with the GSE76039 ATC/PDTC histology test (real ground truth) as the primary external validation. Defer the 4-cohort meta to a per-cohort proxy-label table in supplement, with explicit note that these labels are 8-gene-derived.

3. **DCA and subgroup forest** were computed against proxy labels — re-compute against R1-A leak-free labels (1-day revision-round work).

4. **Add GSE151180 RAI-refractory vs RAI-avid as the most direct clinical ground-truth test** — defer to revision round (probe→symbol mapping needed for the GPL21575 microarray data).

5. **TERT 4-group survival, Hot/Cold composite, DM1/DM2 vs BRAF/RAS orthogonality** — UNCHANGED (independent of 8-gene panel).

## Venue recommendation

→ **Stay with npj Precision Oncology** (Scenario A). With honest reframing, P estimate:
- Honest external: 50-60% (was 45-55% pre-fix; +5pt for transparency about leak-free re-validation)
- Self: 75-80% (was 78-82%; -2pt for slightly lower headline number 0.962 vs 0.954)

Korean cohort outreach (Bundang SNUH) becomes ESSENTIAL for revision round — not just nice-to-have.

## Decision Tree

- IF Korean cohort answers within 1-2 weeks → wait, integrate, submit revision-round-ready version
- IF Korean cohort answers later → submit v6 (this scenario A version) with explicit "Korean cohort in active outreach" commitment
- IF Korean cohort says no → submit v6 anyway, with GSE76039 + GSE151180 (revision-round) + MSK-IMPACT prevalence as the external triangle

## What is genuinely strong (post-fix)

1. 8-gene CV AUC 0.96 against leak-free cluster (NOT circular)
2. 8-gene AUC 0.93 against zero-overlap functional cluster (NOT possibly circular)
3. ΔAUC +0.11 ~ +0.13 over BRAF baseline, consistent across cluster definitions
4. GSE76039 AUC 0.94 with histology ground truth
5. TERT 4-group survival p = 3.8e-5 (independent)
6. Hot/Cold Cohen's d = 1.68 (independent of 8-gene)
7. MSK-IMPACT prevalence cross-cohort (independent)
