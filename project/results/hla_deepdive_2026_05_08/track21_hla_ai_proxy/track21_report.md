# Track 21 — TCGA-THCA HLA-I locus allelic-imbalance read-ratio proxy x DM1

**Track:** Multi-track HLA deep-dive (2026-05-08), Track 21
**Cohort:** TCGA-THCA Primary Tumor (n=505 with DM1 score), Solid Tissue Normal (n=59)
**Caption boilerplate (must appear on every figure):** "HLA-I locus allelic imbalance read-ratio proxy — transcriptomic signal, not allele genotype. Cancer-cohort allele genotyping is out of scope per separation rules."

## 0. Boundary statement (read first)

Per `project/paper2_hla_boundary/HLA_CANCER_SEPARATION_RULES.md` Section 1.2, HLA allele typing on TCGA-THCA (or any cancer cohort) is **forbidden** in Paper 1 territory. This track therefore analyses HLA-A / HLA-B / HLA-C strictly as **per-locus gene-expression read counts** and inter-locus *ratios*. The "allelic imbalance" of the title is **locus-level** (HLA-A vs HLA-B vs HLA-C), not allele-level. We did **not** infer or report any HLA allele genotype, allele frequency, or per-allele read ratio. The narrative below was self-audited for forbidden language (no DRB1/DQB1/etc. alleles named, no carrier-frequency claims, no susceptibility claims, no clinical-risk claims). Section 7 documents the wall-stop pivot away from allele-resolution AI.

## 1. Data

- **Expression matrix:** TCGA pancan log2(norm_count+1) at `project/data/raw/TCGA_pancan/pancan_geneExp.gz`, restricted to thyroid carcinoma primary tumors and solid-tissue normals (572 thyroid samples annotated; 564 had complete HLA-I module + 8-panel expression).
- **DM1 score:** 8-panel anti-thyroid-differentiation module (TG, TPO, TSHR, SLC5A5, FOXE1, PAX8, NKX2-1, DIO1) z-scored against the primary-tumor reference; pancan DM1 (`project/results/paper11_pancancer/pancan_dm1_scored.tsv`) merged where available.
- **Driver labels:** `project/results/dark_matter_phase1/tcga_dm_master_with_pfi.tsv` (BRAF / RAS / Fusion / TripleNeg).
- **Final usable n:** 505 primary tumors with DM1; 59 normals with locus expression.

## 2. Per-locus HLA-I expression x DM1 (Deliverable 1, 2)

Spearman ρ in TCGA-THCA primary tumors (Table T02; Fig F01 forest, F02 scatter):

| Locus | n | rho | 95% CI |
|---|---|---|---|
| HLA-A | 505 | **+0.399** | +0.323, +0.469 |
| HLA-B | 505 | **+0.478** | +0.408, +0.542 |
| HLA-C | 505 | **+0.430** | +0.357, +0.498 |
| B2M | 505 | **+0.409** | +0.334, +0.479 |

**All three classical class-I loci are concordant** (rho range across A/B/C = 0.079, all p<1e-20). HLA-B has the strongest anti-correlation with thyroid differentiation, but the three loci move together. There is **no per-locus signal of one class-I gene preferentially silenced or up-regulated as DM1 rises** at the cohort mean level. This is the expected pattern for a coordinated transcriptomic immune-context shift (consistent with the IFN-γ / NLRC5 / CIITA-driven HLA-I induction program), and *inconsistent* with one allele being differentially silenced.

This is the appropriate first-pass test: if one locus had been preferentially silenced under DM1 immune escape, we would have expected a divergent rho (one near 0 or negative while others positive). We see no such divergence at the locus level.

## 3. HLA-A:B:C ratio + tumor-vs-normal (Deliverable 3)

Per-sample log2 ratios across the three classical loci (Table T03; Fig F03):

| Ratio | Median tumor | Median normal | MW p (tumor vs normal) |
|---|---|---|---|
| log2(A/B) | −0.12 | +0.13 | **5.2e-4** |
| log2(A/C) | +0.59 | +0.70 | 0.40 (NS) |
| log2(B/C) | +0.73 | +0.59 | **7.4e-3** |

In normal thyroid, A and B are roughly equal (median ratio ~0.13); in tumors, B is slightly elevated relative to A (median ratio −0.12). Both significant shifts (A/B and B/C) are in the same direction: **HLA-B is mildly *over*-expressed in tumor relative to normal, not silenced**. This is consistent with IFN-γ-driven HLA-I induction acting strongest on B. The ratio shift is a small *mean* shift, not a heavy-tailed silencing pattern; per-sample IQR is similar in tumor and normal.

## 4. Locus dominance label x DM1 (Deliverables 4, 5)

Per-sample dominance label (z-score of one locus exceeds the median of the other two by >1.0 SD; mirror rule for "silenced") was assigned within primary tumors. Distribution (Table T04):

- balanced: 482 (95.4%)
- A_dominant: 12 (2.4%)
- C_dominant: 4
- C_silenced: 3, B_silenced: 2, A_silenced: 1
- B_dominant: 1

**95% of tumors are "balanced" at the locus level.** The chi-square test for DM1 tertile × locus-dominance label gives chi2=18.50, dof=12, **p=0.101 (not significant)** (Table T05b; Fig F04). With only 23/505 imbalanced samples, the test is also under-powered; even a permissive read of the unadjusted p=0.10 does not support enrichment of any specific locus-imbalance pattern in DM1-high tumors. **Practical conclusion: DM1-high tumors do not preferentially silence one class-I locus** at the resolution this transcriptomic proxy can see.

## 5. B2M re-confirmation (Deliverable 6)

Track 5 reported ρ(DM1, B2M) = +0.42 (n=527, HLA-I module data). This track recomputes on n=505 primary tumors with the same DM1-from-8-panel definition: **ρ = +0.409, p = 8.4e-22** (Table T06). Match within ±0.01. Track 5 finding is **fully reproduced**. B2M is the canonical class-I LOH-target gene; its tracking with DM1 confirms the broader Paper 1 signal (DM1 ↔ class-I induction) extends to the light chain, which is consistent with a coordinated antigen-presentation program rather than allele-specific loss.

## 6. Class-I module heterogeneity (Deliverable 7)

Per-sample coefficient of variation of the 19-gene HLA-I module (linear-space CV) (Table T07; Fig F05):

- Median CV in tumor = 1.455 (n=505)
- Median CV in normal = 1.401 (n=59)
- MW p (tumor vs normal) = **0.43 (NS)**
- Spearman ρ(DM1, CV) = +0.073, p = 0.10 (NS)

**Tumor module CV is not elevated relative to normal, and CV does not track DM1.** If DM1-high tumors had highly heterogeneous class-I expression (one or two genes dropping out), CV would have risen. It does not. This is again consistent with a coordinated up-regulation pattern, not selective silencing.

## 7. Allelic imbalance attempt — wall-stop pivot (Deliverable 8)

**Wall reached. Pivot executed.** Two reasons:

1. **Boundary forbids it.** Allele-resolution AI from RNA-seq requires either (a) haplotype-phased BAMs over HLA-A/B/C exons (LOHHLA, ASE-pipelines, or a phasing pipeline) or (b) per-allele typing followed by SNP read-ratio. **Both routes begin with HLA allele typing on cancer-cohort BAMs**, which is forbidden in Paper 1 per HLA_CANCER_SEPARATION_RULES.md Section 1.2. The "AI proxy" framing of the brief was honored as locus-level imbalance, not allele-level.

2. **No data for it locally.** `find project/data/external -name '*.bam'` returns 0 hits. No `optitype/`, `lohhla/`, or `arcasHLA` outputs for TCGA-THCA exist in the repo (the existing `arcasHLA_GSE213647` and Korean K2 outputs are *non-cancer / Paper 2-3 territory* and are out of scope here).

What was done instead (sections 2–6, deliverables 1–7): per-locus expression z-scores, inter-locus log2 ratios, locus-dominance labels, DM1-tertile chi-square crosstab, B2M re-confirmation, and HLA-I module CV. These are **gene-expression-only** quantities; no allele typing was performed.

## 8. Thorsson 2018 cross-reference (Deliverable 9)

The local Thorsson cache at `project/results/v17p35/tables/FIX2_thorsson_subtype.tsv` is a placeholder ("MISSING" rows; the v17p35 sprint did not pull the per-sample HLA expression / neoantigen burden columns from cBioPortal or the Synapse table). Cross-validation with Thorsson per-sample HLA expression / SNV neoantigen burden columns is therefore **deferred** to a future track. We did not invoke a fresh Thorsson fetch in this track because (a) brief states "IF accessible" and (b) the 2026-05-08 sprint window precludes a network pull.

## 9. Limitations

- **Proxy, not direct LOH.** Locus-level ratios from bulk RNA-seq cannot separate (a) allele-specific loss of one HLA-A/B/C allele on the LOH'd haplotype, (b) transcriptional silencing without DNA loss, (c) IFN-γ-driven coordinated up-regulation, or (d) tumor-cell-intrinsic vs immune-infiltrate-driven HLA-I expression.
- **Bulk RNA hides cell-type heterogeneity.** A DM1-high tumor with HLA-B-dominant infiltrating immune cells could read identically to a tumor-cell HLA-B silencing — bulk cannot tell. Single-cell analysis would be needed to resolve, but is also subject to the same allele-typing boundary.
- **Locus dominance under-powered.** Only 23/505 primary tumors fall into any imbalance class; the chi-square is not adequately powered to detect modest enrichment.
- **No DNA layer.** WES/WGS-based LOHHLA-style true allelic-loss calls were not run (boundary).
- **Reference frame.** DM1 score is an anti-thyroid-differentiation module; it is partially confounded with tumor purity (Track 5's purity-decoupling figure is the relevant control). Track 5 already showed the HLA-I module rho persists after purity-proxy partialling.

## 10. Implication for Paper 1 (residualization sanity)

The Track 5 conclusion — DM1 ↔ HLA-I module is a coordinated transcriptomic shift, not allele-specific silencing — **survives finer-grained inspection at the locus level**:

1. Per-locus rho is concordant across HLA-A / -B / -C (range 0.079, all strongly +).
2. Tumor-vs-normal locus ratios shift in the *induction* direction (HLA-B mildly over-expressed in tumor), not the silencing direction.
3. 95% of tumors are locus-balanced; no DM1-tertile enrichment for any specific locus dominance/silencing.
4. Module CV is not elevated in tumor and does not track DM1.
5. B2M, the canonical class-I LOH-target gene, mirrors HLA-A/B/C (ρ=+0.41 ≈ Track 5's +0.42), confirming the program is class-I-coordinated, not allele-selective.

**For Paper 1 residualization:** when HLA-I gene-expression module is used as a covariate (Track 5 Cox forest), it acts as a single coordinated transcriptomic axis. There is no evidence at the locus level of a hidden allelic-imbalance variable that would split that axis or change the residualization. The Paper 1 framing — "HLA-I module = coordinated transcriptomic immune-context proxy, not allele genotype, deeper HLA = Paper 2 territory" — is supported.

## Files

- Script: `/home/seungho/personal/THCA_data_analysis/scripts/hla_deepdive_2026_05_08/track21/run_track21.py`
- Results: `/home/seungho/personal/THCA_data_analysis/project/results/hla_deepdive_2026_05_08/track21_hla_ai_proxy/`
- Figures (6): F01–F06 in `figs/`
- Tables (9): T01–T08 + T05b in `tables/`
- Summary: `track21_summary.json`

## Self-audit checklist (boundary)

- [x] No HLA allele names typed (DRB1/DQB1/DPB1 etc.).
- [x] No allele-typing pipeline invoked on cancer cohort.
- [x] No allele frequency / carrier frequency claim.
- [x] No HLA-based susceptibility / risk / prognosis / patient-selection claim.
- [x] No HT / GD / AITD HLA claim.
- [x] All figures carry the caption boilerplate.
- [x] Wall-stop pivot documented openly when allele-resolution AI was infeasible (Section 7).
- [x] Paper-1 residualization framing only (Section 10).
