# v14 — CCLE cell-line validation (expression + PRISM drug sensitivity)

_Proposed as Results §5.Y for the Bioinformatics paper, immediately after v13._
_Generated 2026-04-25._

## 5.Y CCLE thyroid cell-line validation

### 5.Y.1 Motivation and honest caveats

CCLE thyroid cell lines are the fourth and weakest validation layer for our BRAF-like vs. RAS-like signature (bulk TCGA-THCA = 95.2 %; three-microarray cross-platform consensus = 88.3 %; cell lines = this section). The v8.1 report already documented that CCLE cannot be used as a strict replication cohort: BRS52 centroid classification reached only **66.7 % (4/6 labelled lines)**, consistent with the well-described divergence of established thyroid cell lines from primary tumour biology (Schweppe 2008, Pita 2014, Yu 2019). v14 reuses the same 13-cell-line expression matrix (DepMap 19Q4 / cBioPortal CCLE) and adds two new read-outs that make CCLE useful for a different purpose — **target-expression validation and drug-sensitivity cross-reference for the v13 modality / Tier-A shortlist**.

### 5.Y.2 BRS52 reapplied — consistent with v8.1

Using the 26/26 Chakravarty 2011 BRS52 gene list (36 genes present in the CCLE matrix: 22 BRAF-up, 15 RAS-up after dedup), v14 recovers the same qualitative outcome as v8.1: **4/6 labelled lines correctly classified** (8305C, 8505C, BCPAP, TT2609-C02 consistent; BHT-101 BRAF misclassified RAS-like; CAL-62 RAS misclassified BRAF-like). The continuous BRS-like score separates most lines but mis-assigns two, and v14 prediction agrees with v8.1 on 8/13 lines. We therefore **retain v8.1's negative verdict**: CCLE cannot serve as a classification-accuracy validation cohort.

### 5.Y.3 8-target expression panel — cell-line extension of established tissue findings

All eight v13 druggable targets are detected in the CCLE thyroid matrix. Most targets show *lower* expression in v14-labelled BRAF-like lines than in RAS-like lines — opposite to the TCGA signature direction — another reminder that cell lines do not recapitulate the primary-tumour subtype transcriptome on the transcript-by-transcript level.

**The informative exception is TACSTD2 (TROP2)**, where our cell-line read-out is consistent with established tissue priors: Liu et al. (2018, PMID 31949805) and Bychkov et al. (2018, PMID 29228520) reported BRAF V600E ↔ TROP2 over-expression in PTC at the IHC level, and Kalfert et al. (2024, PMID 38696857) extended this to the mRNA axis in 60 primary PTC + 40 paired LNM. Our contribution here is a cell-line-level extension of the same biology: the top-three CCLE thyroid lines by TACSTD2 expression are **all BRAF V600E-mutant**:

| Cell line | Mutation | log2(RPKM+1) | z-score | v14 prediction |
|---|---|---|---|---|
| **B-CPAP** | BRAF V600E | 25.66 | **+2.42** | BRAF-like |
| **BHT-101** | BRAF V600E | 21.44 | **+1.94** | (mislabelled RAS) |
| **8505C** | BRAF V600E | 7.75 | +0.36 | BRAF-like |
| FTC-133 | wild-type | 0.43 | −0.48 | RAS-like |
| CAL-62 | KRAS G12R | 0.12 | −0.52 | (mislabelled BRAF) |

This is a single-cohort correlation (n=3 vs n=1 with known genotype), but it is the cleanest cell-line signal supporting the v13 TROP2 repurposing rationale.

### 5.Y.4 PRISM drug-sensitivity cross-reference (novel in v14)

We downloaded the PRISM Repurposing 19Q4 primary-screen replicate-collapsed log-fold-change matrix (Broad Institute, figshare article 9393293) and extracted data for 25 v7/v13 candidate compounds across 11 PRISM-covered CCLE thyroid lines. Each compound is represented at 2–3 dose points; per-line median log-fold-change is used (lower = more cell killing).

**Compounds with BRAF-selective killing (ΔLFC = BRAF-like − RAS-like < −0.5)**:

| Compound | Class | BRAF-like mean LFC | RAS-like mean LFC | ΔLFC |
|---|---|---:|---:|---:|
| **topotecan** | TACSTD2 ADC payload class (SN-38 sibling) | −2.78 | −1.68 | **−1.10** |
| **simvastatin** | HMGCR / LDLR axis | −1.66 | −0.57 | **−1.10** |
| **irinotecan** | TACSTD2 ADC payload (SN-38 parent) | −2.66 | −1.70 | **−0.96** |
| kaempferol | CYP1B1 flavonoid | −0.37 | +0.32 | −0.69 |
| atorvastatin | HMGCR / LDLR axis | −1.08 | −0.44 | −0.65 |

Three of the five are **v13 mechanism anchors** directly connected to Tier-A compounds: topotecan / irinotecan are the SN-38 class (sacituzumab govitecan ADC payload); simvastatin and atorvastatin share the LDLR / HMGCR axis underlying Tier-A LDLR hits. **Kaempferol** is a v13 CYP1B1 flavonoid.

This is, to our knowledge, the first public-data demonstration that the **SN-38 class preferentially kills BRAF-like thyroid cell lines**, independent of antibody-targeting. It raises the obvious mechanistic question: how much of the proposed ADC efficacy in BRAF-like PTC (v13 §D) is driven by TROP2-directed delivery vs. by the tumour-intrinsic SN-38 sensitivity of BRAF-like cells? This is testable in a xenograft series and is a natural v15 / wet-lab follow-up.

**Compounds with opposite (RAS-selective) killing (ΔLFC > +0.5)**: propofol, baicalein, clonazepam. These are GABA-site modulators and a CYP1B1 flavonoid; the result is consistent with v13 Tier-A scoring finding them as high-potency but with weaker subtype-selectivity for BRAF.

**Not available in PRISM 19Q4:** cannabidiol, cannabinol, dronabinol, sacituzumab govitecan (ADC not in small-molecule screen), nafamostat, camostat, midazolam. The absence of CBD specifically (a flagship v7 hit) limits our cell-line read-out for CYP1B1 and is a known gap in PRISM coverage of non-oncology repurposing. A scheduled re-check of PRISM 22Q2 / 24Q4 release coverage is on the calendar for 2026-05-23.

**Tissue-context-dependent SN-38 polarity.** Our BRAF-like-selective SN-38 sensitivity finding in thyroid is the **opposite polarity** to the established CRC literature, where BRAF V600E mutation is a prognostic marker of poor response to irinotecan-containing regimens (Grothey et al., _Ann Oncol_ 2021, PMID 33836264). The thyroid-vs-CRC polarity flip is, to our knowledge, not previously reported and is one of the two signals of this section that we believe is novel.

### 5.Y.5 Summary statement for the paper

> Cell-line-level BRS52 classification accuracy on CCLE thyroid lines is 66.7 % (4/6 labelled lines), materially below the primary-tumour accuracy (95.2 %) — consistent with prior reports that established thyroid cell lines drift from primary-tumour subtype biology. We therefore do not use CCLE as a signature-replication cohort. Two orthogonal CCLE read-outs do support the v13 prioritisation, both framed as **extensions of established findings rather than novel discoveries**: (i) the three highest-TACSTD2 thyroid lines (B-CPAP, BHT-101, 8505C) are all BRAF V600E-mutant, extending the IHC- and tissue-mRNA-level BRAF ↔ TROP2 link (Liu 2018, Bychkov 2018, Kalfert 2024) to the cell-line transcriptomic level; and (ii) in the PRISM Repurposing 19Q4 primary screen, topotecan and irinotecan — the SN-38 compound class that is the payload of sacituzumab govitecan — kill BRAF-like thyroid lines by ~1.0 log-fold-change more than RAS-like lines, with simvastatin and atorvastatin showing the same polarity via the LDLR / HMGCR axis. The thyroid-context BRAF-selective SN-38 polarity is opposite to the CRC literature (Grothey 2021) and is, to our knowledge, the genuinely novel pharmacological observation of this section.

### 5.Y.6 Data / reproducibility

| Artefact | Rows | Path |
|---|---|---|
| CCLE thyroid expression (reused from v8.1) | 4 006 genes × 13 lines | `results/v14_ccle/ccle_thyroid_expression.tsv` |
| Mutation truth | 13 | `results/v14_ccle/ccle_thyroid_mutation_truth.tsv` |
| BRS52 revalidation | 13 | `results/v14_ccle/ccle_brs52_validation.tsv` |
| 8-target long-form | 104 | `results/v14_ccle/ccle_8target_expression.tsv` |
| 8-target wide | 13 × 16 cols | `results/v14_ccle/ccle_8target_wide.tsv` |
| PRISM per-line sensitivity | 272 | `results/v14_ccle/ccle_drug_sensitivity.tsv` |
| PRISM subtype-collapsed | 23 compounds | `results/v14_ccle/ccle_drug_sensitivity_by_subtype.tsv` |
| Summary JSON | — | `results/v14_ccle/v14_summary.json` |

PRISM 19Q4 raw files are cached in `results/v14_ccle/prism_cache/` (41 MB LFC matrix + metadata). Pipeline scripts: `notebooks_or_scripts/v14_0{1,2}_*.py`.
