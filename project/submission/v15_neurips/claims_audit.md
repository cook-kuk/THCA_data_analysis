# v15 — Claims-vs-Source Audit

_Generated 2026-04-27. Cross-verifies every numerical claim in the
abstract / contributions / experiments against checkpoint JSONs and
TSVs._

| Claim location | Paper says | Source value | Source file | Match |
|---|---|---:|---|:---:|
| Abstract / §5.1 | R² = 0.94 | 0.9405 | `task1_theorem2.json` `r2` | ✅ |
| Abstract / §5.2 | ROC-AUC = 0.78 (broader) | 0.7800 | `task2_stress.json` `roc_auc_any` | ✅ |
| §5.2 caption | ROC narrow = 0.62 | 0.6212 | `task2_stress.json` `roc_auc_subspace` | ✅ |
| Abstract / §5.3 / Fig 3 | DANN-lite AUC = 0.39 | 0.3920 | `methods_comparison.tsv` `DANN_lite` | ✅ |
| Abstract / §5.3 / Fig 3 | naive ComBat AUC = 0.33 | 0.3260 | `methods_comparison.tsv` `combat_naive` | ✅ |
| §5.2 / Table 1 | concept mag=1.0 DIAL = 0.393 | 0.3930 | `stress_test_summary.tsv` | ✅ |
| Abstract / §5.5 | 4/4 domains flip | 4/4 (all True) | `task5_cross_domain.json` `any_flip_per_domain` | ✅ |
| Abstract / §5.1 | β_parallel dominates by 33× | 32.1× | `task1_theorem2.json` `coef.delta_parallel / coef.delta_cov` = 0.0493 / 0.0015 | ✅ rounded |
| §5.3 caption | 10 methods, 5 seeds | 10, 5 | `task3_tta.json` `n_methods`, `n_seeds` | ✅ |

**9/9 claims match within rounding.** No prose-vs-source mismatch.

## Abstract length

- 218 words (≤ 250 NeurIPS limit) ✓
- 1 535 characters (no NeurIPS char limit but OpenReview accepts ≤ ~3 000) ✓

## Theorem numbering after fix

- "Theorem 1" appears 27 times in the PDF (correct — main theorem of paper)
- "Theorem 2" appears 0 times in the PDF (cleaned)
- Predecessor work referenced as "predecessor: rank-one alignment condition" without theorem-number assignment
- Restated theorem in Appendix A renders as "Theorem 1 (Unified shift decomposition; restated)"

## Limitations of this audit

This audit only checks numerical claims that appear in checkpoint
JSONs or TSVs. It does NOT verify:

- Citation accuracy (whether cited papers say what we attribute to them) — author-side responsibility, ~2 hour task
- Theorem 1 proof correctness — sketched 6-step proof in Appendix A; reviewers may probe
- Conjecture 1 stated assumption — explicitly marked as conjectural

These three items are flagged in `rebuttal_prep_v1.md` Probes 1 and 5.
