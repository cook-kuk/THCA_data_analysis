# Track 27 — HLA-G + HLA-E + HLA-F non-classical class-I axis × DM1

**Date:** 2026-05-08
**Scope:** Cancer-related (Paper 1 territory + Paper 11 hook). Gene-expression module analysis only.

---

## 0. Boundary

> HLA gene-expression module — not allele genotype.
> Cancer-cohort allele genotyping is out of scope per
> `project/paper2_hla_boundary/HLA_CANCER_SEPARATION_RULES.md` (Sections 1.1, 1.2, 3.1).

All correlations and quadrants in this track use **bulk mRNA expression of HLA loci**, not patient HLA allele typing. We do not claim allele-level association, autoimmune susceptibility, or cancer risk from these data; this is a transcriptomic immune-context module analysis.

---

## 1. Data sources & inputs

| Source | Use | Path |
|---|---|---|
| TCGA pancan expression | per-sample HLA-A/B/C/E/F/G + ligands | `project/data/raw/TCGA_pancan/pancan_geneExp.gz` |
| TCGA pancan phenotype | filter to Primary Tumor | `project/data/raw/TCGA_pancan/phenotype.tsv.gz` |
| Pancan DM1 score | DM1 portable score per sample | `project/results/paper11_pancancer/pancan_dm1_scored.tsv` |
| TCGA-THCA DM master | driver class, PFI, DM1_use | `project/results/dark_matter_phase1/tcga_dm_master_with_pfi.tsv` |
| Track 5 per-sample | reuse DM1_use, driver_simple, purity proxy | `track5_dm1_hla1_module/tables/T01_per_sample_module_scores.tsv` |
| Phase C ICI cohorts | HLA-G × ICI response (4 cohorts) | `paper11_pancancer/phase_C_ICI/processed/{IMvigor210,GSE176307,riaz_GSE91061,MGH_GSE115821}/` |
| 8-gene methylation | global hypomethylation proxy (caveat) | `audit_2026_04_30/round5/r5_2_sample_methylation_8gene.tsv` |

- N TCGA-THCA primary tumors with full module: **n=505 (505 expression × 527 merged with DM master)**.
- All 12 target genes (HLA-A/B/C/E/F/G, LILRB1/2, KIR2DL4, KLRC1/2, KLRD1) present in TCGA pancan and in all 4 ICI cohorts (IMvigor210, GSE176307, riaz_GSE91061, MGH_GSE115821 — note Riaz/MGH lacked KIR2DL4 in some splits but KLRC family complete).

---

## 2. Classical vs non-classical pairwise (TCGA-THCA)

Pairwise Spearman ρ across HLA-A/B/C × HLA-E/F/G (`F01`, `T02`):

| | HLA-E | HLA-F | HLA-G |
|---|---|---|---|
| HLA-A | high | high | mid |
| HLA-B | high | high | mid |
| HLA-C | high | high | mid |

Non-classical loci are correlated with classical, but **HLA-G is the most decoupled** of the three — exactly the divergence we expected (HLA-G is methylation- and tissue-restricted, not constitutively co-regulated with HLA-A/B/C).

Per-locus expression distribution (`F02`): HLA-B has the highest median; HLA-E is high (housekeeping-like); HLA-G has wide dynamic range (most variable across samples), consistent with selective derepression.

---

## 3. DM1 × classical vs non-classical (TCGA-THCA, n=527)

`F03`, `T03_dm1_classical_nonclassical_corr.tsv`:

| Gene | Class | Spearman ρ | 95 % CI | p |
|---|---|---|---|---|
| HLA-A | classical | +0.404 | 0.330–0.473 | 4.3e-22 |
| HLA-B | classical | +0.484 | 0.416–0.547 | 2.5e-32 |
| HLA-C | classical | +0.433 | 0.361–0.500 | 1.6e-25 |
| HLA-E | nonclassical | +0.241 | 0.159–0.320 | 2.1e-08 |
| HLA-F | nonclassical | +0.336 | 0.258–0.410 | 2.2e-15 |
| **HLA-G** | nonclassical | **+0.428** | **0.356–0.496** | **6.4e-25** |

**Key finding:** HLA-G ρ(DM1) = **+0.428** in TCGA-THCA, which is *as strong as classical HLA-A/B/C* (+0.40 to +0.48) and **stronger than HLA-E** (+0.24) and **HLA-F** (+0.34). This is a meaningful upgrade over the per-gene split in Track 5 (which had reported ρ≈+0.43); we now confirm the magnitude with full Fisher CIs and side-by-side classical comparator.

DM1-high THCA tumors thus express *all six* HLA-I loci more — including the tolerogenic HLA-G — not selectively suppressing classical HLA-I.

---

## 4. Decoupling score (z(non-classical) − z(mean A/B/C))

`F04`, `F05`, `T04_decoupling_dm1.tsv`:

| Non-classical gene | n | ρ(decoupling, DM1) | p |
|---|---|---|---|
| HLA-G | 527 | **+0.006** | 0.887 |
| HLA-E | 527 | −0.387 | 3.0e-20 |
| HLA-F | 527 | −0.290 | 1.0e-11 |

**Headline:** **HLA-G tracks classical HLA-I almost perfectly across DM1 levels** — the decoupling score has ρ≈0 with DM1. In contrast, HLA-E and HLA-F decouple **negatively** with DM1 (DM1-high tumors have non-classical-low / classical-high relative to baseline). So in TCGA-THCA, DM1-high tumors are *not* selectively up-regulating HLA-G beyond what classical HLA-I co-regulation would predict — the rise of HLA-G is **proportional** to the general HLA-I induction. This argues against a simple "HLA-G escape" model in DM1; the HLA-G high signal is a co-component of the inflammatory module documented in Track 5/6.

---

## 5. HLA-G ligand axis (LILRB1/ILT2, LILRB2/ILT4, KIR2DL4)

`F06` (top row), `T05_hlaG_ligand_corr.tsv`:

| Ligand | ρ(HLA-G, ligand) | p | ρ(DM1, ligand) | p |
|---|---|---|---|---|
| LILRB1 (ILT2) | +0.381 | 1.1e-19 | +0.142 | 1.1e-3 |
| LILRB2 (ILT4) | +0.447 | 2.6e-27 | +0.192 | 8.6e-6 |
| KIR2DL4 | +0.097 | 0.026 | −0.146 | 8.0e-4 |

LILRB1/2 (ILT2/4 — main inhibitory HLA-G receptors on myeloid cells) co-express with HLA-G strongly (ρ≈+0.38 to +0.45). DM1 also tracks LILRB2 weakly positive (+0.19), consistent with the myeloid-suppressive Track-B-lite signal. **The tolerogenic HLA-G → ILT2/4 axis is biochemically plausible at the bulk-tumor level in DM1-high THCA**, though ρ alone doesn't prove engagement.

KIR2DL4 (NK receptor) is essentially decoupled from HLA-G in this cohort — KIR2DL4 expression in bulk RNA-seq is sparse and may reflect NK fraction rather than the receptor pathway.

---

## 6. HLA-E ligand axis (KLRC1/NKG2A, KLRC2/NKG2C, KLRD1/CD94)

`F06` (bottom row), `T06_hlaE_ligand_corr.tsv`:

| Ligand | ρ(HLA-E, ligand) | p | ρ(DM1, ligand) | p |
|---|---|---|---|---|
| KLRC1 (NKG2A, inhibitory) | +0.574 | 1.9e-47 | +0.338 | 1.4e-15 |
| KLRC2 (NKG2C, activating) | +0.580 | 1.2e-48 | +0.302 | 1.4e-12 |
| KLRD1 (CD94) | +0.644 | 5.4e-63 | +0.078 | 0.075 |

HLA-E co-expresses with both inhibitory (NKG2A) and activating (NKG2C) ligand subunits at near-identical strength. DM1 tracks both NKG2A and NKG2C similarly (+0.30 to +0.34) — i.e., **DM1-high tumors recruit NK cells but with a balanced inhibitory/activating receptor profile, not a clean NKG2A-dominant suppression**. CD94 (the obligate dimerization partner) is the strongest correlate of HLA-E expression but only weakly tied to DM1, consistent with CD94 being broadly expressed across NK/T subsets.

---

## 7. DM1-high × HLA-G-high quadrant (TCGA-THCA)

`F07` (KM curves), `T07_quadrant_summary_THCA.tsv`, `T08_survival_THCA_dm1_hlaG_PFI.tsv`:

Quadrant | n | median HLA-G | median DM1 | median classical_z | median nonclassical_z
---|---|---|---|---|---
DM1lo_HLAGlo | 178 | 6.02 | −0.94 | −0.84 | −0.78
DM1lo_HLAGhi | 86 | 8.58 | −0.78 | +0.22 | +0.36
DM1hi_HLAGlo | 86 | 6.74 | −0.22 | +0.01 | −0.37
DM1hi_HLAGhi | **177** | 9.37 | +0.61 | +0.54 | +0.54

The DM1hi_HLAGhi quadrant (n=177, ~34 % of cohort) is a coherent group: high DM1, high HLA-G, high classical HLA-I, high non-classical HLA-I — i.e., **the inflamed dark-matter tumor**.

PFI Cox (n=482, 50 events):
- DM1 alone: HR=1.36 per σ, p=0.010 — DM1 is prognostic, replicates Round 3.
- HLA-G alone: HR=1.04, p=0.76 — HLA-G alone has no prognostic signal in THCA.
- DM1 + HLA-G + interaction: DM1 HR=1.40 p=0.007; HLA-G HR=0.95 p=0.66; **interaction HR=1.09 p=0.51** — no significant DM1×HLA-G interaction on PFI in THCA.

Bottom line: in TCGA-THCA the HLA-G-high status by itself does **not** add prognostic value over DM1, and there is no DM1×HLA-G synergy. The HLA-G axis here is a *passenger* of DM1-driven HLA-I induction at the population level. THCA event count (50) is small; we caveat this is power-limited and the pancan setting (Sec. 8) is the relevant generalization test.

---

## 8. Pan-cancer HLA-G distribution (`F08`, `F09`, `T09`, `T10`)

**Top 10 lineages by median HLA-G expression** (n=10,593 primary tumors, 33 lineages):

| Rank | Lineage | n | Median HLA-G |
|---|---|---|---|
| 1 | kidney clear cell carcinoma | 533 | 9.06 |
| 2 | diffuse large B-cell lymphoma | 48 | 8.14 |
| 3 | kidney papillary cell carcinoma | 290 | 7.81 |
| 4 | **thyroid carcinoma** | 505 | 7.73 |
| 5 | head & neck squamous | 520 | 7.40 |
| 6 | skin cutaneous melanoma | 103 | 7.37 |
| 7 | pancreatic adenocarcinoma | 178 | 7.35 |
| 8 | lung adenocarcinoma | 515 | 6.98 |
| 9 | uveal melanoma | 80 | 6.86 |
| 10 | stomach adenocarcinoma | 415 | 6.74 |

Sanity: choriocarcinoma is not a TCGA lineage (so the trophoblast positive control isn't in this cohort), but **clear cell renal, DLBCL, and head & neck** are well-known HLA-G-high cancers — pan-cancer ranking is biologically sensible. THCA is **#4 / 33** in absolute HLA-G expression.

**Per-lineage ρ(DM1, HLA-G)** (`T10_pancan_per_lineage_dm1_nonclassical_corr.tsv`):

| Rank | Lineage | n | ρ(DM1, HLA-G) | p |
|---|---|---|---|---|
| **1** | **thyroid carcinoma** | 505 | **+0.438** | 4.0e-25 |
| 2 | testicular germ cell | 150 | +0.279 | 5.6e-4 |
| 3 | thymoma | 120 | +0.248 | 6.3e-3 |
| 4 | brain LGG | 516 | +0.190 | 1.3e-5 |
| 5 | GBM | 160 | +0.187 | 0.018 |
| 6 | pheochromocytoma & paraganglioma | 179 | +0.179 | 0.016 |

**Thyroid is the #1 lineage** for DM1↔HLA-G coupling pan-cancer — twice the magnitude of the next lineage (testis germ-cell, an immune-privileged site). This is the strongest DM1-locus link we have observed in the non-classical class-I axis.

---

## 9. ICI cohorts: HLA-G vs response (4 phase-C cohorts)

`F10`, `T11_ici_HLAG_response.tsv`:

| Cohort | n_R | n_NR | median HLA-G (R) | median HLA-G (NR) | MW p (HLA-G) | Cohen's d | MW p (decoupling) | d_decoupling |
|---|---|---|---|---|---|---|---|---|
| IMvigor210 (urothelial, anti-PD-L1) | 68 | 230 | 1.71 | 1.45 | 0.123 | **+0.25** | 0.71 | +0.06 |
| GSE176307 (urothelial, atezo) | 11 | 50 | −1.22 | −1.00 | 0.086 | **−0.42** | 0.052 | −0.55 |
| riaz_GSE91061 (mel, anti-PD-1) | 10 | 39 | +0.10 | +0.25 | 0.147 | **−0.41** | 0.28 | −0.42 |
| MGH_GSE115821 | n.a. | n.a. | — | — | — | — | — | — |

(MGH cohort dropped: pre-treatment R/NR ≥3 not satisfied after timepoint filter.)

**Result:** mixed-direction across cohorts (IMvigor +0.25 favors R; both Riaz and GSE176307 ≈−0.4 favor NR), no cohort reaches FDR significance; pooled signal is weak and direction-inconsistent. The decoupling score in GSE176307 trends d=−0.55 (p=0.05) — HLA-G *down*-regulated relative to classical HLA-I in responders — which would be consistent with a mild "less tolerogenic" pattern in ICI-responsive urothelial tumors but is not robust.

**Interpretation: no robust HLA-G→ICI signal in our 3 evaluable cohorts** — Reviewer Q reserve only, not a primary claim.

Per-sample ICI table at `T11b_ici_per_sample_nonclassical.tsv` (all 4 cohorts × HLA loci × ligands).

---

## 10. Methylation × HLA-G

`F11`, `T12_meth_per_gene.tsv`, `T13_meth_dm1_hlaG_summary.tsv`:

Available proxy: **mean β across 8 thyroid-differentiation genes** (`mean_8g_beta`, n=523 with HLA-G overlap). HLA-G-specific HM450 probe data is **not currently in the repo** — true HLA-G promoter methylation requires fetching the HM450 matrix and the HLA-G CpG island.

| Comparison | n | ρ | p |
|---|---|---|---|
| DM1 vs mean_8g_beta | 523 | **+0.491** | 5.1e-33 |
| DM1 vs HLA-G | 523 | **+0.431** | 5.1e-25 |
| mean_8g_beta vs HLA-G | 523 | **+0.562** | 7.3e-45 |
| HLA-G vs DM1 partial | mean_8g_beta | 523 | +0.221 | 3.3e-7 |

Per-gene methylation correlation (mean_8g_beta vs each HLA gene): **all six HLA loci correlate +0.55 to +0.68 with mean_8g_beta** — meaning higher 8-gene panel methylation (DM1 phenotype) co-occurs with higher HLA-I module expression broadly. This is **not a HLA-G-specific epigenetic signature**; it's the DM1↔inflammation co-direction expressed through both methylation and expression layers.

After partialling out mean_8g_beta, DM1↔HLA-G ρ drops from +0.43 to +0.22 (still p=3.3e-7), meaning **about half of the DM1↔HLA-G signal is co-explained by the DM1-typical methylation phenotype**. The other half is residual DM1-specific signal.

**Important caveat:** mean_8g_beta is hyper-methylation of *thyroid-differentiation gene promoters* in DM1, not global hypomethylation. We do not have HLA-G CpG-island data; this "methylation × HLA-G" finding is therefore a co-occurrence with the DM1 methylation phenotype, not a mechanism statement on HLA-G promoter regulation. The Round 4 memory (DM1 mean_8g_beta d=−1.75 p=3.3e-39; sign convention: in DM1, panel-gene promoters are *hyper*-methylated) is the right framing.

---

## 11. Driver-stratified HLA-G

`F12`, `T14_driver_stratified_HLAG.tsv`, `T15_driver_pairwise_HLAG.tsv`:

| Driver | n | Median HLA-G | Median decoupling | ρ(DM1, HLA-G) | p |
|---|---|---|---|---|---|
| BRAF | 294 | **8.43** | +0.097 | +0.049 | 0.40 |
| RAS | 54 | 5.53 | −0.114 | +0.138 | 0.32 |
| Triple-negative (mut-neg, fusion-neg) | 168 | 6.88 | **−0.212** | **+0.405** | 5.4e-8 |

Pairwise (`T15`):
- TripleNeg vs BRAF: MW p=1.5e-16, d=−0.88 (BRAF >> TripleNeg in HLA-G).
- TripleNeg vs RAS: MW p=3.1e-8, d=+0.91 (TripleNeg >> RAS).

**Counter-intuitive but interesting result:** HLA-G is **highest in BRAF**-mutant THCA (median 8.4) and **lowest in RAS** (median 5.5). The triple-negative compartment is intermediate. **But the only driver class where DM1↔HLA-G correlation is significant is triple-negative** (ρ=+0.405, p=5e-8). In BRAF tumors HLA-G is uniformly high regardless of DM1 (BRAF→HLA-G is driver-driven, not DM1-driven); in triple-negative tumors HLA-G tracks DM1 — i.e., **the DM1-driven HLA-G signal lives in the triple-negative compartment**, exactly the compartment Paper 1's molecular dark matter targets.

This supports the Paper-1 framing: in BRAF/RAS-negative (triple-neg) THCA, the dark-matter axis is also where the tolerogenic HLA-G axis activates.

---

## 12. Limitations

1. **Bulk RNA-seq HLA isoform-blind.** HLA-G has ~7 alternatively spliced isoforms (HLA-G1–G7); soluble HLA-G5/G6 vs membrane HLA-G1 have different immune-evasion roles. Our HLA-G value sums all transcripts and cannot distinguish.
2. **Bulk = mixture.** HLA-G expression in bulk tumor can come from tumor cells, infiltrating myeloid cells, or trophoblast-like cells; we cannot resolve cell-of-origin without scRNA. (Track 12 has scRNA HLA-I context that could be re-interrogated.)
3. **Correlative immune-evasion claim.** ρ(HLA-G, LILRB2)≈+0.45 demonstrates module co-expression, **not engagement** of the tolerogenic axis. Functional confirmation requires cell-cell ligand-receptor scoring on scRNA / spatial.
4. **HLA-G-specific methylation absent.** We use 8-gene panel β as a proxy; this conflates DM1-phenotype methylation with HLA-G promoter regulation. The Track-26/Round-4 pipeline could be re-run on HLA-G CpG islands when full HM450 is available.
5. **ICI cohorts underpowered for non-classical HLA.** Only 3/4 cohorts had pre-treatment R/NR ≥3; effect sizes ±0.3 to ±0.4 in either direction; pooled meta would be at chance. Treat as Reviewer-Q reserve.
6. **No allele-level claim.** Per HLA_CANCER_SEPARATION_RULES.md, we cannot infer patient HLA-G genotype (the *G\*01:01:01* vs *G\*01:04* etc. distinctions or the well-known +14bp/-14bp 3' UTR insertion polymorphism) from bulk tumor RNA — that requires germline-source HLA typing not in scope.
7. **THCA event-deficient.** PFI = 50 events on n=482 — power to detect HR=1.1 effects of HLA-G is low. Pancan thyroid+other lineages (Track 10) is the right next test.

---

## 13. Implications

### For Paper 1
- **HLA-G ρ(DM1) = +0.428 in TCGA-THCA** is a clean addition to the inflammatory-context narrative: DM1-high tumors up-regulate the **tolerogenic** HLA-G in lockstep with classical HLA-A/B/C — i.e., the hot-but-tolerant phenotype. This complements Track 5 (classical HLA-I module with DM1) by showing the *tolerogenic* arm rises in parallel.
- **Driver-stratified result is paper-relevant:** in the BRAF/RAS-negative compartment, DM1↔HLA-G ρ=+0.405 (p=5e-8). This is the compartment Paper 1 frames as the dark-matter sub-stratifier, and the HLA-G coupling is exactly there.
- **The decoupling-score result** (HLA-G tracks classical HLA-I; HLA-E and HLA-F decouple negatively) means we *cannot* claim "HLA-G is selectively up-regulated as immune-escape in DM1." The honest framing is: DM1-high tumors are inflamed and **all six HLA-I loci including HLA-G rise together**.

### For Paper 11 (pan-cancer)
- **Thyroid is #1 lineage** for DM1↔HLA-G correlation pan-cancer (rho +0.438, ~2× the next lineage). Adds a novel non-classical-HLA result to the pan-cancer DM1 axis manuscript.
- HLA-G median expression ranks THCA #4/33 — alongside RCC, DLBCL, papillary kidney — also consistent with the DM1 / inflammatory dark-matter cluster pattern in those cancers.

### For Paper 3 (ICI vulnerability) — Track-B-lite reserve
- **HLA-G x ICI response is direction-inconsistent across 4 cohorts** (IMvigor +d, urothelial GSE176307 −d, melanoma Riaz −d). Not a primary claim. Reviewer Q reserve only.

---

## Outputs

- Results dir: `/home/seungho/personal/THCA_data_analysis/project/results/hla_deepdive_2026_05_08/track27_hla_g_nonclassical/`
- Figures (12): `figs/F01_classical_vs_nonclassical_pairs.{png,pdf}` … `figs/F12_driver_stratified_HLAG.{png,pdf}`
- Tables (15): `tables/T01_thca_per_sample.tsv` … `tables/T15_driver_pairwise_HLAG.tsv`
- Summary JSON: `track27_summary.json`
- Script: `/home/seungho/personal/THCA_data_analysis/scripts/hla_deepdive_2026_05_08/track27/run_track27.py`
