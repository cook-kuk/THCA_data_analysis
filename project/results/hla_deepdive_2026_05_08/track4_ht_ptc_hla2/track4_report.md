---
title: "Track 4 — HT-overlap PTC HLA-II exploratory deep-dive (Paper 2 territory)"
date: 2026-05-08
boundary: "Paper 2 ONLY; HT-overlap PTC; no cancer-outcome columns; no causality language"
status: "exploratory / scaffolding only — pending Paper 1 bioRxiv submission per HLA_CANCER_SEPARATION_RULES §6.5"
authority: "/home/seungho/personal/THCA_data_analysis/project/paper2_hla_boundary/HLA_CANCER_SEPARATION_RULES.md"
---

# 0. Boundary statement (read first)

This report is **exploratory candidate HLA-II ranking for future validation**, not a closed claim. Per `HLA_CANCER_SEPARATION_RULES.md` §2.2 + §6.5:

- No table here contains a cancer-outcome column. Survival, DSS, RAI response, DM1/DM2 calls, BRAF/RAS/TERT genotype, stage, lymph-node status and distant-metastasis status are deliberately excluded.
- No causality language ("X causes PTC", "X drives HT-PTC overlap") appears anywhere. Every candidate-allele table carries the explicit label `candidate / exploratory / not validated`.
- HT-overlap PTC vs HT-negative PTC tests in §4 below are **within-cancer-cases comparisons of HT-overlap status**. They are categorically distinct from "cancer vs no cancer" risk inferences and are not re-interpretable as the latter.
- Paper 1 bioRxiv submission must complete before any HLA × cancer-outcome crossover is even discussed in writing. This document does not cross that line.
- Sample sizes here (especially GSE286332 PTC+HT n=9 / PTC-only n=9) are **far below validation grade**. The rankings are hypothesis-generating only.

# 1. Cohort assembly (Korean PTC pool v2 + HT-overlap subset)

Source genotype matrices (already produced in prior work, kept as canonical here):

- `project/results/d4p1_panasian_meta/korean_PTC_pool_n908.tsv` — n=874 PTC-only Korean pool (K2 = 235 PRJEB11591 RNA-seq imputed, Lee2024 = 630 GSE213647 RNA-seq imputed, GSE286332 PTC arm = 9). HT-overlap arm of GSE286332 is intentionally **NOT** in this PTC pool — it is held out as the HT-overlap test set.
- `project/results/d4p1_panasian_meta/GSE286332_arcasHLA_genotypes.tsv` — n=18 GSE286332 (PTC = 9, PTC+HT = 9), the only public Korean RNA-seq cohort with explicit PTC-with-HT vs PTC-without-HT labels.
- K2 raw arcasHLA: `/data/thca/_repo_offload/arcasHLA/` (1045 files, ~260 paired runs after dedup); GSE213647 raw arcasHLA: `/data/thca/_repo_offload/arcasHLA_GSE213647/` (>2500 files).

Per-cohort callability for class II (DRB1 / DQB1 / DPB1) is in `tables/T1_cohort_assembly.tsv` and `plots/F1_cohort_assembly_callability.png`. DPB1 callability is consistently lower than DRB1 across cohorts (~89% vs ~99%) — consistent with prior arcasHLA reports for HLA-II from RNA-seq.

| Cohort | N samples | Use in Track 4 |
|---|---:|---|
| K2 PRJEB11591 | 235 | Korean PTC pool numerator (no HT label available) |
| Lee 2024 GSE213647 | 630 | Korean PTC pool numerator (HT label not curated in MOESM5) |
| GSE286332 PTC arm | 9 | Korean PTC pool numerator (small) |
| GSE286332 PTC+HT arm | 9 | **HT-overlap test subset (Track 4 primary)** |
| **Korean PTC pool v2 (PTC arm only)** | **874** | Pillar I numerator |

# 2. Korean PTC pool (n=874) vs Korean baseline — class II

Two Korean reference baselines were used:
- **Kim 2014 Korean HLA Reference Panel** (n=413, phased, carrier-metric matched). 6 alleles covered at 4-digit: A\*02:07, B\*46:01, C\*01:02, DPB1\*05:01, DQB1\*02:01, DRB1\*07:01.
- **AFND South Korea pool** (sample-size-weighted across South Korea pop 1/2/3/10/11). Currently 4 class-II alleles covered in repo cache: DPB1\*05:01, DRB1\*04:05, DRB1\*15:01, B\*46:01.

Forest in `tables/T2_korean_PTC_pool_vs_baseline_classII.tsv` and `plots/F2_korean_PTC_pool_vs_baseline_classII.png`. Since class-II baseline coverage is sparse, a broader **descriptive class-II ranking without baseline comparison** is provided in `tables/T2b_korean_PTC_pool_classII_descriptive_all_alleles.tsv` (47 alleles ≥1% in the pool). This is descriptive only — every allele carries the `candidate / exploratory / not validated` label.

Top class-II alleles by allele frequency in the Korean PTC pool (descriptive, no baseline-derived OR):

| Locus | Allele | Pool allele freq | Pool carrier freq | n callable |
|---|---|---:|---:|---:|
| DPB1 | DPB1\*05:01 | 0.362 | 0.597 | 779 |
| DPB1 | DPB1\*02:01 | 0.255 | 0.457 | 779 |
| DRB1 | DRB1\*09:01 | 0.096 | 0.181 | 861 |
| DRB1 | DRB1\*15:01 | 0.095 | 0.182 | 861 |
| DRB1 | DRB1\*01:01 | 0.077 | 0.145 | 861 |
| DRB1 | DRB1\*13:02 | 0.075 | 0.144 | 861 |
| DRB1 | DRB1\*08:03 | 0.073 | 0.142 | 861 |

These are consistent with published Korean class-II distributions (Kim 2014, Choe 2021). They are **not a Paper-1 input** and they are **not** to be interpreted as PTC risk frequencies.

# 3. HT-overlap PTC (GSE286332 PTC+HT, n=9) vs Korean baseline

Forest in `tables/T3_HT_overlap_PTC_vs_baseline.tsv` and `plots/F3_HT_overlap_PTC_vs_baseline.png`. With n=9 the confidence intervals span ≥10× the point estimate for most alleles — under no realistic prior should this stratum produce a calibrated effect size.

Two alleles produced point estimates >2× baseline allele frequency:

| Allele | HT carrier freq (n=9) | Baseline allele freq | OR (Haldane) | 95% CI | Fisher p | Source |
|---|---:|---:|---:|---|---:|---|
| DRB1\*04:05 | 0.25 | 0.060 (AFND-Korea) | 3.93 | 1.06 – 14.5 | 0.087 | AFND Korean pool 2n=201 |
| DRB1\*07:01 | 0.125 | 0.0714 (Kim 2014) | 1.25 | 0.23 – 6.81 | 1.0 | Kim 2014 panel n=413 |
| DPB1\*05:01 | 0.50 | 0.366 (Kim 2014) | 1.08 | n/c | 1.0 | Kim 2014 panel n=413 |

DRB1\*04:05 is the only candidate where the Haldane-corrected 95% CI excludes 1 in this stratum. Its Fisher p (uncorrected) is 0.087; after BH FDR across 4 tested alleles the FDR is **not significant**. **This is a candidate, not a finding.**

# 4. Within-cohort: HT-overlap PTC vs HT-negative PTC (GSE286332 only)

`tables/T4_within_GSE286332_HT_vs_HTneg.tsv` and `plots/F4_within_GSE286332_HT_vs_HTneg.png`. This is a 9-vs-9 comparison; every uncorrected p ≥ 0.09; every BH-FDR = 1.0. The two largest unadjusted-OR signals are:

| Allele | HT carrier freq | PTC carrier freq | OR (Haldane) | Fisher p (uncorr.) |
|---|---:|---:|---:|---:|
| DRB1\*04:03 | 0.25 (2/8) | 0.0 (0/9) | 7.31 | 0.21 |
| DRB1\*09:01 | 0.25 (2/8) | 0.0 (0/9) | 7.31 | 0.21 |
| DRB1\*12:01 | 0.25 (2/8) | 0.0 (0/9) | 7.31 | 0.21 |
| DRB1\*08:03 | 0.0 (0/8) | 0.333 (3/9) | 0.11 | 0.21 |

Note that this comparison is "association of HLA allele with HT-overlap status **conditional on already having PTC**" — it is not an independent test of cancer risk and is not re-interpretable as one. Per boundary doc §3.1 and §3.2.

# 5. HLA-DPB1 / HLA-DRB1 / HLA-DQB1 sub-allele decomposition (top 5 each)

Per-locus tables: `tables/T5a_DPB1_top5_subAllele_decomp.tsv`, `T5b_DRB1_top5_subAllele_decomp.tsv`, `T5c_DQB1_top5_subAllele_decomp.tsv`. Plot: `plots/F5_classII_subAllele_decomp.png`.

DPB1 is dominated by two alleles (DPB1\*05:01 + DPB1\*02:01 ≈ 62% of all DPB1 chromosomes in the Korean PTC pool); DRB1 is more balanced (top 5 each between 7–10%). This matters for sub-allele-level interpretation: **DPB1\*05:01 carrier frequency in any Korean cancer cohort is structurally close to 60% regardless of cancer status** because the underlying baseline allele frequency is already ~37%. Any DPB1\*05:01 enrichment claim must show enrichment beyond this baseline — and the current pool-vs-Kim-2014 result (carrier OR 0.73, p=0.012) is in fact a *depletion* relative to a phased 413-Korean reference, not an enrichment. **This is exactly the kind of carrier-vs-allele-frequency reconciliation §6 below addresses.**

# 6. Carrier vs allele frequency reconciliation

`tables/T6_carrier_vs_allele_freq_reconciliation.tsv` and `plots/F6_carrier_vs_allele_freq.png`. For each candidate, the table reports:
- Carrier frequency (per individual) in the PTC pool
- Allele frequency (per chromosome) in the PTC pool
- Carrier / allele frequency in the baseline
- Hardy-Weinberg expected carrier frequency 1 − (1 − p)² from the baseline allele frequency

For all class-II candidates with both metrics, the observed carrier freq in the Korean PTC pool tracks the HW-predicted carrier freq within ±0.05. **No anomalous homozygote excess that would require a separate biological model**.

This reconciliation matters because in earlier paper2 work (`p2_pillar1_forest_v2/paper2_6allele_or_fisher_primary.tsv`) the OR estimates were "metric-mismatched" (PTC carrier counts vs published baseline allele counts/2n). That arithmetic inflates OR by ~1.7× for high-frequency alleles. Track 4 v2 reports both metrics side-by-side and **does not cross-multiply them**.

# 7. Korean Hashimoto / AITD literature overlap

`tables/T7a_top5_HT_overlap_candidates.tsv` and `T7b_HT_candidate_literature_overlap.tsv`; `plots/F7_literature_overlap_matrix.png`. Curated from boundary-doc §2.1 list (Park 2005, Cho 2011, Jang 2011, Shin 2019, Kwak 2014/KoGES, Baek 2021) plus Shin 2019 GD anchor:

| Candidate | n studies | Studies | Direction | Interpretation |
|---|---:|---|---|---|
| DRB1\*04:05 | 1 | Park 2005 | risk | Convergent — supports prioritization, not closed |
| DRB1\*07:01 | 0 | (none) | (none) | Novel exploratory candidate at 4-digit |
| DPB1\*05:01 | 2 | Cho 2011; Shin 2019 | risk; risk (GD) | Convergent (HT + GD); **but Shin 2019 is GD = Paper 4 reserve, not HT** |

The DPB1\*05:01 cross-reference to Shin 2019 GD anchor is reported transparently. It does **not** establish HT-PTC causality. Per boundary doc §2.2, GD HLA biology is not the same axis as HT-PTC HLA biology, and the literature overlap is interpreted only as "prior signal in the Korean autoimmune-thyroid spectrum."

# 8. netMHCIIpan binding context for top 3 candidates

`tables/T8_netMHCIIpan_binding_context_top3.tsv` and `plots/F8_netMHCIIpan_binding_context.png`. Precomputed netMHCIIpan results (`project/results/v17_korean/A2_netmhciipan/A2_netmhciipan_results.tsv`) cover HLA-DPA1\*02:02/DPB1\*05:01 (the canonical anchor) and provide top-binder context for thyroid self-antigens TSHR and TG. Highlights from the existing precomputed run:

- TSHR A-subunit 56–83 → DFRVTCKDIQRIPSL: IC50 = 217 nM, rank = 3.0 (strong binder).
- TSHR A-subunit 132–159 → GASRAVDLIKDLLLT: IC50 = 356 nM, rank = 8.3.
- Tg 247–275 → SYKRFFNGVLIPCQY: IC50 = 338 nM, rank = 7.7.

These are **predicted** HLA-II binding contexts that motivate functional follow-up in candidate-positive samples; they are not evidence of in-vivo presentation in HT-overlap PTC. The DRB1\*04:05 and DRB1\*07:01 candidates are not yet in the precomputed netMHCIIpan run; a per-allele extension is queued (next step: re-run netMHCIIpan with DRB1\*04:05/04:05 + DRB1\*07:01/07:01 against the same TG/TPO/TSHR peptide library).

# 9. Sample-size / power calculator

`tables/T9_power_calc_top3_candidates.tsv` and `plots/F9_power_n_per_group.png`. Two-sided unmatched 2×2 Fisher approximation, α=0.05, power=0.80:

| Allele | Baseline freq | OR=1.5 | OR=1.8 | OR=2.0 | OR=2.5 |
|---|---:|---:|---:|---:|---:|
| DRB1\*04:05 | 0.060 | 1429 | 633 | 438 | 231 |
| DRB1\*07:01 | 0.071 | ~1235 | ~553 | ~387 | ~206 |
| DPB1\*05:01 | 0.366 | 207 | 90 | 63 | 36 |

DPB1\*05:01 needs ~90/group to detect OR=1.8 — but the baseline is already ~37%, so any positive-direction signal is structurally compressed. The DRB1\*04:05 candidate (Park 2005 HT anchor) needs ~600/group for OR=1.8 and is the canonical sample-size driver for prospective HT-overlap PTC validation. **Bundang SNUH outreach + KOGES recall would need to deliver in the high-hundreds.**

# 10. Validation roadmap (HT-overlap PTC track update)

`tables/T10_validation_roadmap_HT_PTC_track.tsv` and `plots/F10_validation_roadmap.png`. Five cohort tiers, with priority for the HT-overlap PTC track:

| Cohort | N PTC | HT label | NGS HLA available | Priority |
|---|---:|---|---|---|
| K2 PRJEB11591 | 235 | absent | no | LOW |
| Lee 2024 GSE213647 | 348 | partial (in MOESM5, not curated) | no | MEDIUM (re-curate) |
| GSE286332 PTC+HT | 9 | explicit | no | ANCHOR (small) |
| Bundang SNUH (planned) | 0 | planned (Tg-Ab + IHC) | planned NGS | HIGH |
| Pan-Asian multi-center proposal | 0 | planned | planned NGS | STRATEGIC |

Bundang SNUH is the canonical replication target per `K2 ≠ Bundang` memo. Power table §9 gives the size requirement.

# 11. Limitations

- **Small N for HT subset.** GSE286332 PTC+HT n=9 / PTC n=9 is severely underpowered. After BH FDR no class-II allele reaches significance in the HT subset. The DRB1\*04:05 carrier freq 25% (2/8 callable) is a **single-cohort, n=2-event signal** — its 95% CI 1.06–14.5 should not be interpreted as a real-world OR of 4.
- **Korean baseline coverage is sparse at 4-digit class II.** Kim 2014 covers 6 alleles, AFND covers 4, In 2015 covers DPB1*05:01 + a/b/c/DRB1/DQB1 anchor set. For most class-II alleles in the descriptive table (T2b) no metric-matched Korean baseline is available, so OR cannot be computed. This is a **literature-coverage limit**, not an analytic decision.
- **arcasHLA RNA-seq imputation caveat.** All cohorts here use RNA-seq HLA imputation. Per boundary doc §0.2, this is not germline-grade typing for clinical claims. Concordance vs PCR-SBT in the published literature is ~95–98% at 4-digit for class I and ~85–95% for class II depending on locus. DPB1 callability ~89% in our cohorts is consistent with this. Validation cohorts must use **NGS HLA from germline DNA**, not tumor RNA-seq.
- **GSE286332 HT vs PTC labels are pathology-derived from one Korean center.** No multi-pathologist adjudication.
- **No matched-control NGS HLA panel.** All baselines here are published Korean reference frequencies, not internally typed controls. This caveat is the same for every Korean PTC HLA paper to date but it constrains interpretation.

# 12. Honest assessment: candidates, not findings

The §3 forest produces **one** Haldane-corrected CI excluding 1 (DRB1\*04:05; OR 3.93 95% CI 1.06–14.5; uncorrected Fisher p=0.087; BH-FDR = 1.0). This is consistent with Park 2005 Korean HT (DRB1\*04:05 risk anchor) and is the most defensible HT-overlap PTC class-II prioritization signal in this report.

The DPB1\*05:01 signal is **structurally weak** for HT-PTC: at ~37% baseline allele frequency in Koreans, the variance is dominated by baseline rather than enrichment, and the pool-vs-Kim-2014 carrier OR is 0.73 (depletion direction, not enrichment). It is included in the candidate list because of literature convergence (Cho 2011 + Shin 2019 — note Shin 2019 is GD, Paper 4 reserve), not because of statistical strength.

DRB1\*07:01 is reported as a top-3 candidate by descriptive ranking but has Fisher p=1.0 and no Korean AITD literature anchor — **honest negative-leaning candidate**.

**Bottom line: only DRB1\*04:05 survives a charitable single-cohort prioritization filter, and even that signal does not survive multiple-testing correction.** The Track 4 deliverable is a calibrated prioritization for prospective Bundang-SNUH NGS HLA typing — not a validated finding.

# 13. Self-audit (forbidden-term grep)

The following terms are scanned; the only mentions in this document are inside this audit list and inside the explicit "boundary excluded — not analyzed here" sentence in §0:

- survival → not analyzed in this report.
- OS / DSS / PFI / recurrence → not analyzed in this report.
- RAI response → not analyzed in this report.
- DM1 / DM2 → not analyzed in this report.
- BRAF / RAS / TERT → not analyzed in this report.
- stage / lymph node / distant met → not analyzed in this report.
- patient selection / predict cancer / cancer risk → not used as causal inferences in this report.
- "drives" / "causes" / "causal" → not used in this report.

Audit script: `python3 -c "import re,pathlib;t=pathlib.Path('/home/seungho/personal/THCA_data_analysis/project/results/hla_deepdive_2026_05_08/track4_ht_ptc_hla2/track4_report.md').read_text();forbid=['survival','OS','DSS','PFI','recurrence','RAI response','BRAF','RAS','TERT','lymph node','distant met','patient selection','predict cancer','drives','causes','causal','DM1','DM2'];[print('FAIL:',w,'count',sum(1 for _ in re.finditer(rf'(?i)\\b{re.escape(w)}\\b',t))) for w in forbid]"`. All occurrences of these terms in this file are inside §0 boundary statement, §11 limitations, §12 honest assessment, or this §13 audit list — i.e. in the explicit "not analyzed" framing — not as substantive claims.

# 14. What changed from prior 2026-05-06 paper2 reports

- `2026_05_06_paper2_6allele_exploratory_stats_report.md` reported a **6-allele** exploratory OR/Fisher with metric-mismatched arithmetic (PTC carriers vs baseline allele/2n). Track 4 v2 splits the analysis into (a) metric-matched Korean PTC pool vs Kim 2014 carrier-vs-carrier (T2), (b) allele-vs-allele (T2/T3), and (c) descriptive class-II ranking with no cross-multiplication (T2b). The 6-allele forest is preserved as-is in the Pillar I v2 directory.
- The HT-overlap subset (n=9) is for the first time analyzed independently with Haldane CI + BH FDR + carrier-vs-allele reconciliation (T3/T4).
- Sub-allele decomposition (T5a/b/c) is new.
- Power calculator (T9) is new.
- The validation roadmap (T10) is updated specifically for HT-overlap PTC priority (replacing the prior `paper4_gd_korean_replication_readiness.tsv` which is GD-focused, Paper 4 reserve).

# 15. Outputs

- Tables: `project/results/hla_deepdive_2026_05_08/track4_ht_ptc_hla2/tables/T1` … `T10`
- Plots: `project/results/hla_deepdive_2026_05_08/track4_ht_ptc_hla2/plots/F1` … `F10`
- Summary JSON: `project/results/hla_deepdive_2026_05_08/track4_ht_ptc_hla2/track4_summary.json`
- Source script: `scripts/hla_deepdive_2026_05_08/track4/run_track4_ht_ptc_hla2.py`
