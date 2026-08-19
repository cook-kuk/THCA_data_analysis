# Paper E · Synthesis manuscript draft

## Title
**A unified leakage-aware framework for neoantigen-vaccine prioritization, ICI-vulnerability prediction, and synthetic-lethality discovery across pancreatic and thyroid cancer**

Alternative title: *"From dark thyroid cancer to KRAS-driven pancreatic cancer: an integrated systematic-evaluation framework for personalized immunotherapy"*

---

## Why this synthesis paper

The Lumenix project has produced five interlocking lines of work. Individually each is publishable; together they form a coherent framework for cold-tumor immunotherapy that no individual paper captures. This synthesis paper:

1. Frames neoantigen-vaccine prediction (Paper 5, this session) and ICI-vulnerability axis discovery (Paper 3) as two sides of the same problem — both ask "which immunologically-cold cancer subtype can we re-sensitize?"
2. Demonstrates the framework on two opposite-extreme cancers — KRAS-driven pancreatic (PAAD, ~40% G12D) and driver-poor thyroid (THCA, lowest TMB).
3. Identifies synthetic-lethality vulnerabilities for combo strategies (vaccine + ICI; BRAF inhibitor + vaccine + ICI).
4. Provides a leakage-aware, HLA-aware, source-bias-audited benchmark resource that any future immunotherapy prediction paper should adopt.

---

## Abstract

**Background.** Cold tumors — those with low tumor mutational burden (TMB) and poor immune-checkpoint-inhibitor (ICI) response — represent a major unmet need in cancer immunotherapy. Two complementary computational strategies are emerging: (i) neoantigen-vaccine prediction, which seeks to prime de-novo T-cell responses against shared driver-mutation peptides, and (ii) ICI-vulnerability axis discovery, which identifies immune-context features that predict ICI response without relying on TMB alone. Both rely on cross-cohort generalization that has been under-validated.

**Methods.** We built a leakage-aware 13-source post-TESLA neoantigen-vaccine benchmark (n=119,237 cleaned master records, n=84,549 labeled, leakage-free n=19,817) with explicit `test_set_safety` annotation, and systematically evaluated peptide-HLA representations under five strict-generalization split designs. We re-analyzed two contrasting cancer types: pancreatic adenocarcinoma (PAAD, KRAS-driven, low TMB) and thyroid carcinoma (THCA, driver-poor, lowest TMB), both targets of an integrated 5-paper Lumenix pipeline.

**Results.** Conventional random k-fold CV substantially overestimates strict generalization (source-only shortcut classifier reaches AUROC 0.93 random-CV but 0.50 LOSO). Among representations, ESM2-150M protein-language-model embeddings + Random Forest achieves the highest mean LOSO AUROC (0.75 ± 0.16, n=3 sources), with NEPdb-out lifted from 0.40 (biophys baseline) to 0.86. MHCflurry binding-affinity features alone reach 0.92 on TESLA-out — task-aligned binding features dominate PLM for some sources. PLM scaling saturates at 150M parameters. We re-prioritize KRAS G12D vaccine candidates for PAAD/COAD (top: GADGVGKSAL, computational score 0.69 under strict training) and KRAS G12C for LUAD (smoking signature). For THCA, the per-allele HLA-A*24:02 model achieves AUROC 0.74 — the highest of any individual HLA — directly supporting Korean-cohort prioritization. We synthesize this with prior Lumenix findings on the ICI vulnerability axis of dark thyroid cancer (anaplastic vs differentiated, d=+2.53 myeloid / -2.51 thyroid_diff), showing that ICI-resistance and neoantigen-prediction are linked by shared immunogenomic context.

**Conclusions.** A leakage-aware, HLA-stratified neoantigen-vaccine benchmark is a necessary substrate for cold-tumor immunotherapy prediction. Representation choice (ESM2-150M, mhcflurry binding features) matters more than architecture under strict generalization. KRAS G12D in PAAD and BRAF V600E + Hashimoto-PTC in THCA are computationally prioritized synthetic-lethality combination targets for vaccine + ICI strategies, pending wet-lab validation.

---

## Three pillars of the framework

### Pillar 1 · Leakage-aware benchmark + representation comparison (Paper 5 = this session)
- 13 sources · 119,237 records · leakage-free 19,817
- Random CV inflated 0.85-0.96 → LOSO 0.40-0.50 — source memorization quantified
- ESM2-150M sweet spot · mhcflurry 0.92 on TESLA-out
- HLA-stratified analysis attenuates pooled biophysical attributions (Simpson's paradox)
- Reusable resource for any cancer immunotherapy prediction paper

### Pillar 2 · ICI vulnerability axis (Paper 3 frozen)
- ATC vs DTC: d=+2.53 myeloid / -2.51 thyroid_diff (4/4 sign across 4 thyroid cohorts)
- RAI-refractory module on same axis as ATC immunogenomic phenotype
- HT-PTC subgroup: TLS d=+1.96, IGHV clonality d>0.5, AICDA up
- Track B-lite: pooled n=415 across 8 cohorts; ATC d=+2.53 myeloid / -2.51 thyroid_diff (4/4 sign)
- Provides immune-context features that predict ICI response without TMB

### Pillar 3 · Synthetic lethality / combination targeting
- **PAAD**: KRAS G12D vaccine + ICI + (potential KRAS G12D inhibitor sotorasib-like)
- **LUAD**: KRAS G12C vaccine + sotorasib + ICI (smoking signature)
- **THCA**: BRAF V600E vaccine + dabrafenib + anti-PD-1 in HT-PTC subgroup
- Combination strategy because monotherapy ICI has <5% ORR in PAAD and ~0% in DTC

---

## Two complementary cancer cases

### Case 1 · Pancreatic adenocarcinoma (PAAD)

| Property | Value |
|---|---|
| Median TMB | ~3 mut/Mb |
| KRAS-mutant fraction | >90% |
| Dominant mutation | G12D ~45% |
| ICI monotherapy ORR | <5% |
| Top vaccine candidate | KRAS G12D · GADGVGKSAL (HLA-A*11:01 restricted) |
| TCGA-PAAD G12D Cox HR | 2.17 (p=0.002, n=185) |
| Pan-cancer meta G12D | n=3,113, HR=1.41 |
| Korean A*11:01 freq | 9.0% (vs A*02:01 11.3%) |
| World A*11:01 freq | 9.6% |

**Synthetic lethality framing**: KRAS G12D mutation is the single largest oncogenic driver of PAAD; it is also the worst-survival molecular subgroup (HR=2.17). Combining a G12D peptide vaccine (re-prioritized score 0.69 strict-LR) with anti-PD-1 ICI creates the priming step that ICI alone cannot deliver in PAAD's immune-excluded TME.

### Case 2 · Thyroid carcinoma (THCA)

| Property | Value |
|---|---|
| Median TMB | ~0.5 mut/Mb (lowest of major solids) |
| BRAF-mutant fraction (PTC) | ~60% |
| ICI monotherapy ORR (DTC) | ~0% |
| ATC d (myeloid) | +2.53 vs DTC (4/4 sign) |
| ATC d (thyroid_diff) | -2.51 vs DTC (4/4 sign) |
| HT-PTC TLS d | +1.96 |
| Per-allele A*24:02 AUROC (this work) | 0.74 ± 0.03 |
| Korean A*24:02 master under-rep | 10× |

**Synthetic lethality framing**: dark thyroid cancer is paradoxical — TMB is essentially zero, yet the HT-PTC subgroup is immune-active and the ATC subtype shows myeloid-dominant immune infiltration. <b>BRAF V600E vaccine + dabrafenib (BRAF inhibitor) + anti-PD-1</b> in the HT-PTC subset specifically exploits both the existing HT-driven immune priming and the BRAF-mutant tumor dependency.

### Why PAAD and THCA together?

PAAD and THCA are at opposite extremes of the TMB axis (PAAD ~3, THCA ~0.5) but face the same clinical problem: poor ICI response. The unified framework asks: in either extreme, can a leakage-aware benchmark + representation-aware predictor + ICI-vulnerability axis identify a synthetic-lethality combination? The answer is conditional on representation choice and strict-evaluation discipline. We provide both.

---

## Six headline results

1. **Conventional random CV overestimates neoantigen prediction by 0.4 AUROC**: source_only shortcut classifier reaches AUROC 0.926 random-CV but 0.500 LOSO.
2. **ESM2-150M is the LOSO sweet spot**: 0.75 ± 0.16, with NEPdb-out 0.86 (biophys baseline 0.40).
3. **MHCflurry binding-affinity (2-d) wins on some sources**: AUROC 0.92 on TESLA-out — task-aligned binding features beat 1280-d PLM there.
4. **PLM scaling saturates at 150M**: 35M → 150M → 650M gives 0.74 → 0.75 → 0.71. Larger PLM not the answer.
5. **HLA-A*24:02 per-allele AUROC = 0.74**, the highest of any single HLA. Korean cohort is the missing dataset.
6. **KRAS G12D for PAAD/COAD vs G12C for LUAD** is a cancer-type-specific prioritization with verifiable population-frequency × score expected values.

---

## Sections

### 1. Introduction
- Cold-tumor immunotherapy unmet need
- TMB-based prediction limits
- Two computational strategies: vaccine prediction + ICI vulnerability
- Lumenix integrated framework rationale

### 2. Results
- 2.1 Leakage-aware benchmark construction (Paper 5)
- 2.2 Random CV overestimates strict generalization
- 2.3 Representation comparison: ESM2 vs mhcflurry vs biophys
- 2.4 PLM scaling saturates at 150M
- 2.5 HLA-stratified analysis (per-allele A*24:02 = 0.74)
- 2.6 Case 1 — PAAD KRAS G12D candidate prioritization
- 2.7 Case 2 — THCA ICI vulnerability axis (Paper 3)
- 2.8 Synthetic-lethality combination framing
- 2.9 Korean cohort gap and Paper 4 backlog

### 3. Discussion
- Why representation matters more than architecture for cross-source generalization
- Why TMB-blind ICI vulnerability axis is necessary for dark cancers
- Combination therapy strategy under strict-generalization framework
- Field needs standard splits and leakage audit

### 4. Limitations
- No wet-lab validation in this work
- Population coverage uses simplified independence
- ICI vulnerability axis is bulk RNA-seq based
- Synthetic-lethality is speculative without preclinical models
- Korean cohort backlog not yet ingested

### 5. Methods (full reproducibility)
- Datasets, cleaning, leakage audit
- Feature engineering, models, splits
- Bootstrap CI, Wilcoxon paired tests
- HLA stratification and shortcut tests
- Candidate prioritization with autoimmunity safety

### 6. Data + Code availability
- Public hub: 40.82.129.113/papers_hub_2026_05_04/
- Frozen sha256 manifest in `manuscript/data_snapshots/MANIFEST.tsv`

---

## Why this synthesis paper is publishable

| Strength | Evidence |
|---|---|
| Novel framework | Leakage-aware + HLA-stratified + ICI-vulnerability axis combined |
| Two contrasting cancer cases | PAAD (KRAS-driven, low-TMB) + THCA (driver-poor, lowest-TMB) |
| Cross-paper integration | References 5 Lumenix papers (P1 DM1, P2 HT-isolated, P3 ICI, P4 K2-backlog, P5 benchmark) |
| Honest evaluation | Random CV vs LOSO, Simpson's paradox, no wet-lab overclaim |
| Reusable resource | 13-source benchmark + 6 RCSB pMHC + ESMFold structures + 4 ESM2 caches |
| Synthetic-lethality angle | Combo strategy concrete (KRAS G12D + ICI; BRAF + dabrafenib + ICI) |

## Target venues

| Tier | Venue | Fit | Risk |
|---|---|---|---|
| **A** | **Nature Communications** | Cross-cancer immunotherapy framework | Medium (no wet-lab) |
| **A** | **Cell Reports Medicine** | Translational synthetic-lethality angle | Medium |
| **A-** | **Briefings in Bioinformatics** | Strong methodology + representation comparison | Low ★ |
| **A-** | **Cancer Immunology Research** | Synthetic-lethality + ICI vulnerability angle | Low |
| **B+** | **NAR Database / Bioinformatics App Note** | Resource-paper companion | Very low |

**Recommendation**: Submit Paper 5 (leakage-aware benchmark) first to **BIB** (1-2 weeks), then this synthesis Paper E to **Cancer Immunology Research** or **Cell Reports Medicine** as a follow-up that builds on the BIB paper.

---

## Two-track manuscript strategy

```
[NOW]            Paper 5 → BIB submit (verified ESM2-150M, mhcflurry, KRAS candidate)
[+1 week]        Paper B (resource companion) → Bioinformatics App Note
[+4 weeks]       Paper E (this synthesis) → Cancer Immunology Research / Cell Rep Med
[parallel track] Paper 1 DM1 → npj submission (already ship-ready)
[parallel track] Paper 2 HT-isolated → ongoing
[parallel track] Paper 3 ICI Track A FROZEN → submit when wet-lab not gating
[parallel track] Paper 4 Korean GD HLA → backlog, gated on K2/Bundang FFPE
```

Total realistic output: **5-7 papers** across the next 6-12 months from the integrated Lumenix framework.

---

**Status**: synthesis manuscript outline ready; full draft requires 3-4 days of polish + Phase 3 ICI cohort details.
**Decision pending**: Seungho approval for Paper E to be added to active submission queue.
