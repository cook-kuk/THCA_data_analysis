# Paper 2 HLA — DQB1*02:01 zero-cell QC audit report
**Date:** 2026-05-06
**Author:** Seungho Cook
**Scope:** Korean PTC pool n = 874 (3 subcohorts: K2 n=235, Lee2024 n=630, GSE286332_PTC n=9–18)
**Question:** Is DQB1*02:01 = 0 carriers in this pool a real biological depletion or a technical artifact of the arcasHLA imputation pipeline?

**Final intermediate verdict (this audit):** **Technically suspicious — likely arcasHLA subtype-resolution artifact.** Locus-level callability fails (28% missing) and the DQB1*02 *family* is present at the expected Korean ~13% rate, but split among DQB1*02:02 and ~40 rare/suspect subtypes instead of the canonical DQB1*02:01. The zero-cell depletion signal **does not pass the locus-level callability + resolution checks** required to remain a high-priority candidate.

This audit does *not* claim a real association either way. It recommends suspending the DQB1*02:01-specific depletion claim until orthogonal NGS-typer validation is performed.

---

## 1. Inputs

- Per-sample 4-digit calls: `project/results/d4p1_panasian_meta/korean_PTC_pool_n908.tsv` (n=874 rows, 12 4-digit columns + cohort)
- K2 raw arcasHLA: `project/results/v17_korean/arcasHLA/K2_arcasHLA_genotypes.tsv`
- GSE286332 raw arcasHLA: `project/results/d4p1_panasian_meta/GSE286332_arcasHLA_genotypes.tsv`
- Per-sample arcasHLA JSON: `project/results/v17_korean/arcasHLA_GSE286332/*.genotype.json` and `*.genes.json`

---

## 2. Per-locus callability (callable = both alleles present at 4-digit)

| Locus | Callable both | Total | Callable % |
|---|---|---|---|
| A | 874 | 874 | **100.0%** |
| B | 874 | 874 | **100.0%** |
| C | 873 | 874 | **99.9%** |
| DRB1 | 861 | 874 | **98.5%** |
| DPB1 | 779 | 874 | **89.1%** |
| **DQB1** | **631** | **874** | **72.2%** ⚠ |

**Finding:** DQB1 has the **lowest callability by ≈17–28 percentage points** vs all other loci. ~243 samples (27.8%) lack DQB1 calls entirely. This is the first technical red flag: a true biological zero-rate would not preferentially co-occur with elevated locus-level missingness.

See `F21_DQB1_locus_callability_barplot.png` for visual.

---

## 3. DQB1 4-digit allele distribution (within 631 callable samples)

- Unique 4-digit DQB1 alleles observed: **67**
- Unique 2-digit DQB1 families observed: **5** (DQB1*02, *03, *04, *05, *06 — biologically expected)
- Top-1 allele share: 10.3% (no single allele dominates → diversity is plausible *within* the called set)

---

## 4. DQB1*02 family detail

| 4-digit subtype | Observations (allele copies) | Plausibility |
|---|---|---|
| **DQB1*02:02** | **53** | common Korean variant |
| DQB1*02:276 | 14 | ⚠ extremely high subtype number — typing-suspect |
| DQB1*02:253 | 10 | ⚠ high subtype number — typing-suspect |
| DQB1*02:225Q | 4 | ⚠ Q = "questionable expression" suffix; rare |
| DQB1*02:199 | 3 | ⚠ rare/suspect |
| DQB1*02:214N | 1 | ⚠ N = "null allele" suffix; non-expressed |
| DQB1*02:272 | 1 | ⚠ rare/suspect |
| DQB1*02:215 | 1 | ⚠ rare/suspect |
| DQB1*02:274 | 1 | ⚠ rare/suspect |
| DQB1*02:205 | 1 | ⚠ rare/suspect |
| **DQB1*02:01** | **0** | **target — completely absent** |
| **TOTAL DQB1*02 family** | **89** | — |

**Pool-wide DQB1*02 family carrier count = 84** (sum of `carriers_DQB1*02_any` across subcohorts).
**Pool-wide DQB1*02 family carrier rate = 84 / 631 callable = 13.3%** — *in the published Korean range* (typically 13–17%).

**This is the smoking gun.** The DQB1*02 family is present at the expected Korean rate, but the canonical DQB1*02:01 (the dominant *02 subtype in published Korean HLA literature, e.g., In JW 2015) is **completely absent**, while large counts (28 of 89 *02-family observations = 31%) are spread across rare or suspect subtype-number assignments (*02:253, *02:276, *02:225Q, *02:199, *02:214N, etc.).

In published Korean HLA cohorts (e.g., In JW 2015 Ann Lab Med, n=613×2), DQB1*02:01 is the dominant *02 subtype at ~2.1% allele frequency. Expected count in n=631 callable samples (1262 alleles) at 2.1% AF: **~26 alleles, ~25 carriers** under HWE-like assumptions. **Observed: 0.**

Probability of zero observed carriers if true rate is 2.1% is ≈ exp(−25) ≈ 10⁻¹¹ — biologically near-impossible *unless* the typing pipeline is misclassifying reads.

---

## 5. Subcohort consistency

| Subcohort | n_total | DQB1 callable | DQB1*02:01 carriers | DQB1*02 family any carriers | Family-level rate |
|---|---|---|---|---|---|
| K2 | 235 | 163 (69.4%) | 0 | 18 | 11.0% |
| Lee2024 | 630 | 464 (73.7%) | 0 | 66 | 14.2% |
| GSE286332_PTC | 9–18 | 4 (44.4%) | 0 | 0 | 0% (n too small to interpret) |
| **Pool** | **874** | **631 (72.2%)** | **0** | **84** | **13.3%** |

**Finding:** The DQB1*02:01 absence is **consistent across all three subcohorts** at the 4-digit level — but the DQB1*02 family carrier rate is also consistent with Korean published baselines at the 2-digit level (K2 11%, Lee2024 14%, pool 13.3%).

A real biological *02:01 deficiency would show as either (a) elevated *02:02 alone or (b) reduced overall *02 family. We observe neither. We observe a **redistribution within the *02 family from canonical *02:01/*02:02 toward 40 rare/suspect subtypes**.

K2 raw arcasHLA file (n=63 sub-sample with full arcasHLA output) has DQB1 callable in 100% of those 63 — so K2 callability degradation in the pool comes from samples where DQB1 was added from a non-arcasHLA call source (typing collapse / external manifest). This itself is suspicious: it means the K2 235 cohort merges arcasHLA calls with a separate pipeline whose DQB1 resolution behavior differs.

---

## 6. Cross-locus diversity comparison

| Locus | n observations | Unique 4-digit | Unique 2-digit families | Top-1 share |
|---|---|---|---|---|
| A | 1748 | 64 | 13 | 21.2% |
| B | 1748 | 86 | 26 | 9.7% |
| C | 1746 | 89 | 11 | 12.9% |
| DRB1 | 1722 | 36 | 13 | 9.6% |
| **DQB1** | **1262** | **67** | **5** | **10.3%** |
| DPB1 | 1558 | 41 | 39 | 36.2% |

DQB1's 5 2-digit families is biologically expected (DQB1 has the smallest 2-digit family count among classical loci). Top-1 share is in the normal range. **Diversity at 2-digit level is normal**; the anomaly is the 4-digit subtype assignment within DQB1*02.

---

## 7. Rare/impossible allele patterns

**Observations with subtype number ≥ 50 (potential typing-suspect):** 140 across 40 unique alleles in DQB1.

Top such alleles: DQB1*03:483 (35 obs), DQB1*02:276 (14), DQB1*02:253 (10), DQB1*03:567 (9), DQB1*04:93 (6), DQB1*05:359 (5), DQB1*06:335 (5), DQB1*03:94 (5), DQB1*05:193 (4), DQB1*03:117 (4)…

Subtype numbers in the 100s–500s are uncommon in published Korean cohorts (most carriers fall into single-/double-digit subtypes). 140/1262 ≈ **11.1% of DQB1 allele observations are in the high-subtype-number tail** — consistent with arcasHLA assigning ambiguous-coverage reads to rare allele numbers rather than collapsing to the canonical 4-digit code.

By comparison, K2's separate arcasHLA file (the source where the full pipeline ran with all loci called) shows DQB1 alleles like DQB1*02:02:01, DQB1*05:01:01, DQB1*03:02:01 at the *6-digit* level — well-resolved canonical alleles. **The pool table merges multiple typing pipelines; the high-subtype anomaly is concentrated in samples that did not pass through the full arcasHLA + canonical-collapse step.**

---

## 8. arcasHLA per-sample QC (GSE286332 = 18 samples)

- 18 `*.genotype.json` files present, 18 `*.genes.json` files present.
- 12 of 18 (66.7%) have a non-empty DQB1 entry in `genotype.json`.
- 6 of 18 (33.3%) have an *empty* DQB1 list — i.e., arcasHLA failed to call DQB1.
- For the 12 successful DQB1 calls: distributions are dominated by DQB1*03 (multiple subtypes), with one DQB1*02:02 observation but **zero DQB1*02:01**.

This is consistent with the pool finding: arcasHLA fails DQB1 in a meaningful fraction of samples, and on the called fraction, *02:01 is replaced by *02:02 + non-canonical subtype numbers.

---

## 9. Verdict — biological vs technical

| Audit check | Pass? | Implication |
|---|---|---|
| DQB1 locus callability ≥ 95% | ❌ 72.2% | technical issue likely |
| DQB1*02 family rate matches published Korean ~13–17% | ✅ 13.3% | family is present at expected level |
| Canonical DQB1*02:01 dominates within *02 family | ❌ absent (vs *02:02 + rare *02:xxx) | resolution-collapse artifact |
| Rare-subtype tail < 5% of locus observations | ❌ 11.1% have subtype ≥ 50 | typing-pipeline ambiguity |
| Cross-pipeline consistency (K2 raw arcasHLA fully resolved) | ⚠ partial | pool merges arcasHLA + non-arcasHLA |
| Subtype absence consistent across subcohorts | ✅ all 3 zero | systematic, not random |

**Verdict: technically suspicious — most consistent with arcasHLA / typing-pipeline subtype-collapse artifact, not a biological depletion.**

The zero-cell signal at DQB1*02:01 should **not** be reported as a Korean PTC depletion candidate until orthogonal NGS-typer validation (e.g., HLA-LA, OptiType, T1K) confirms whether the canonical *02:01 vs the rare-tail *02:xxx assignment difference reflects biology or pipeline behavior.

---

## 10. Final allowed conclusion

> "DQB1*02:01 remains a high-priority depletion candidate **only if** locus-level callability and resolution checks do not indicate technical failure. In the current Korean PTC pool, locus-level callability (72.2%) and 4-digit subtype distribution within the DQB1*02 family (canonical *02:01 absent, rare-tail *02:xxx subtypes inflated) **do indicate technical issues**, so DQB1*02:01 zero-cell depletion is **suspended pending orthogonal NGS-typer validation**. The DQB1*02 family-level carrier rate (13.3%) is consistent with published Korean baselines, supporting the technical-artifact interpretation."

This conclusion supersedes any prior text in `paper2_hla_professor_packet_kr.md` and `paper2_hla_claim_boundary.md` that listed DQB1*02:01 as a depletion candidate surviving strict harmonization. Those documents should be read with this audit's caveat applied.

---

## 11. Recommended next actions (validation only)

1. **Re-type DQB1 with an orthogonal NGS HLA caller** on at least 50 samples (sampled across K2, Lee2024). Compare 4-digit assignments to arcasHLA pool calls. If *02:xxx high-subtype tail collapses to *02:01 in orthogonal typer, the artifact is confirmed.
2. **Check arcasHLA reference IMGT/HLA database version** used for the pool. Newer IMGT releases include many rare *02:xxx subtypes that may be absent from the typical In JW 2015 published Korean cohorts (which used older databases) — this version-skew alone could explain the rare-tail inflation.
3. **Audit per-sample DQB1 call confidence** (arcasHLA `genes.json` reports per-allele likelihoods). Samples whose DQB1 calls fall in the rare *02:xxx tail likely have low confidence; flag these as "DQB1 ambiguous" rather than as carriers/non-carriers of *02:01.
4. **Re-state Paper 2 manuscript framing** to remove any DQB1*02:01-specific claim. Remaining valid claims (after this audit): family-level 13.3% carrier rate is in the Korean published range. Nothing more.
5. **Do NOT** connect this audit's findings to AITD/Hashimoto/Graves' biology. The zero-cell signal is no longer interpretable as biology.

---

## 12. Files produced

- `project/results/p2_pillar1_forest_v2/paper2_dqb1_zero_cell_qc_audit.tsv` — consolidated per-locus + DQB1 specific + subcohort metrics
- `project/results/p2_pillar1_forest_v2/_paper2_dqb1_audit_locus_callability.tsv`
- `project/results/p2_pillar1_forest_v2/_paper2_dqb1_audit_subcohort.tsv`
- `project/results/p2_pillar1_forest_v2/_paper2_dqb1_audit_diversity.tsv`
- `project/results/p2_pillar1_forest_v2/_paper2_dqb1_audit_K2_missing_samples.tsv` *(empty — K2 arcasHLA file has 100% DQB1 calls; missingness is in the merged pool)*
- `project/papers_hub_2026_05_04/assets/paper2_hla/F20_DQB1_zero_cell_qc_matrix.png` (3-panel matrix)
- `project/papers_hub_2026_05_04/assets/paper2_hla/F21_DQB1_locus_callability_barplot.png` (callability + diversity overlay)
- `project/notebooks_or_scripts/paper2_dqb1_zero_cell_qc_audit.py` (this audit script)

---

**End of QC audit. No association claim is implied. The DQB1*02:01-specific depletion signal is suspended pending orthogonal NGS-typer validation.**
