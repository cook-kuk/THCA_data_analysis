# Advance paper feasibility map after Paper 1 freeze

**Date:** 2026-05-06  
**Status:** strategy / feasibility memo only. No new analysis. No new data download. No Paper 1 edit.  
**Paper 1 fixed identity:** “A thyroid-lineage differentiation axis stratifies thyroid cancer beyond canonical driver mutations.”  
**Important boundary:** Paper 1 is waiting/packaging mode. The ideas below are for separate advance-paper trajectories, not Paper 1 expansion.

---

## 1. Executive verdict

Paper 1 should not absorb more driver/fusion/metabolism/germline material. The strongest next move is to separate the candidate ideas into executable tracks:

| Rank | Candidate | Verdict | Why |
|---:|---|---|---|
| 1 | Paper 2 HLA Pillar I v2 | **Most executable now** | Uses published Korean baseline / AITD sources; no raw data; already in progress. |
| 2 | BRAF/RAS-negative molecular taxonomy | **Potentially viable, but not Paper 1** | Existing dark-matter analyses found DICER1/EIF1AX/PPM1D enrichment and single-cell FVPTC/cPTC gradient, but survival claim failed. |
| 3 | Multi-driver BRAF/RAS/TERT/RET/NTRK map | **Useful support layer, weak as standalone** | Already has TERT recovery and driver maps; novelty is moderate unless external driver/fusion cohort is added. |
| 4 | Synthetic-lethal / metabolism vulnerability map | **Good future Paper 9, not now** | Conceptually strong after Paper 1 is published, but requires DepMap/CCLE/PRISM and careful SL framing. |
| 5 | Korean NRG1 germline x somatic | **Blocked** | No matched Korean germline + somatic data in project. Needs KoGES/Bundang/dbGaP access. |
| 6 | Ethnicity/age/methylation comparison | **High 삽질 risk now** | Requires controlled Korean genomic/epigenomic data and causal environmental interpretation. |

---

## 2. Existing work already found in repo

| Topic | Existing files | What they already show |
|---|---|---|
| Dark matter viability | `project/results/dark_matter_phase1/viability_dark_matter.md` | TCGA BRAF/RAS-negative dark matter = 137/482 = 28.4%; DICER1/EIF1AX/PPM1D enrichment in one cluster; OS underpowered. |
| Single-cell dark matter follow-up | `project/results/dark_matter_phase1/viability_dark_matter_step4_5.md` | 8-gene score behaves as a continuous FVPTC/cPTC differentiation gradient; single-cell survival/prognostic framing not viable. |
| Full phase2 dark matter dashboard | `project/results/dark_matter_phase2/web/v2.html` | Large explanatory dashboard already exists; useful as archive, not current Paper 1 identity. |
| Driver class phenotype | `project/results/dark_matter_phase2/p2d_driver_class_phenotype.tsv` | Driver class table exists; true driver-negative and DICER1/EIF1AX classes are separated. |
| v17 driver/fusion tables | `project/results/v17/tables/driver_landscape_v17_summary.tsv`, `project/results/v17/tables/dark_matter_fusion_overlay.tsv` | Driver/fusion overlay exists; rare fusion counts are small. |
| TERT recovery | `project/results/v17_tert_recovery/` | TERT crosstab and external BRAF/TERT work exist. |
| Korean NRG1 plan | `project/results/audit_2026_04_29/nrg1_plan.md` | Explicitly blocked on germline data. |
| Paper 9 metabolism / SL plan | `project/papers_hub_2026_05_04/paper9_plan.html` | Good future translational plan, but marked post-marathon / after Paper 1. |
| Paper 2 HLA | `project/reports/2026_05_06_hla_more_data_registry.md`, `project/reports/2026_05_06_paper2_korean_baseline_extraction_report.md` | Current most executable advance track. |

---

## 3. Candidate-by-candidate feasibility

### 3.1 BRAF/RAS-negative “Dark Matter” 신규 driver / target paper

**User idea:** BRAF V600E and RAS hotspot negative samples, then search fusion genes, lncRNA/miRNA networks, DICER1/EIF1AX/non-BRAF/RAS cluster, RAI resistance / recurrence.

**Existing evidence:**

- TCGA BRAF/RAS-negative compartment already tested: **137/482 = 28.4%**.
- Within dark matter, DICER1/EIF1AX/PPM1D enrichment existed in the older frame:
  - DM2: 6/55 = 10.9%
  - DM1: 1/81 = 1.2%
  - Fisher p = 0.0175 in the archived analysis.
- Single-cell follow-up reframed the result:
  - 8-gene score tracks a **continuous FVPTC/cPTC differentiation gradient**, not a clean binary tumor subtype.
  - Within-tumor heterogeneity was stronger than between-sample variation in the tested sc dataset.
- Survival/prognostic claim failed or was endpoint-limited:
  - TCGA-THCA OS events are too sparse.
  - Within-DM outcome claim should not be headline.

**Feasibility verdict:** **Medium-high for molecular taxonomy; low for prognostic/therapy paper.**

**Safe version:**

> A BRAF/RAS-negative thyroid cancer molecular-taxonomy paper centered on DICER1/EIF1AX/PPM1D enrichment and FVPTC/cPTC transcriptional architecture.

**Unsafe version:**

> New therapeutic target discovery, fusion-driven dark matter, RAI resistance proven, recurrence validated.

**Why it is not Paper 1 now:** Paper 1 has already moved to a driver-orthogonal lineage axis. Reopening dark-matter/fusion framing would weaken the clean Paper 1 identity.

---

### 3.2 BRAF/RAS + TERT / high-risk multi-mutation signature

**User idea:** BRAF only, RAS only, BRAF+TERT, RAS+TERT, other high-risk combinations, triple-negative; compare prognosis, RAI genes, TDS.

**Existing evidence:**

- `v17_tert_recovery` already contains TERT integration and BRAF/TERT tables.
- `dark_matter_phase2/p2d_driver_class_phenotype.tsv` already separates:
  - Class1_BRAF_V600E: n=291
  - Class2_RAS_hotspot: n=54
  - Class4_DICER1_EIF1AX: n=7
  - Class5_TERT_only: n=5
  - Class6_True_driver_neg: n=125
- TERT-only and co-mutation groups are small in TCGA-THCA.
- This is useful for claim boundary and supplementary stratification, but it is not a very novel standalone story unless an external cohort with reliable TERT/fusion calls is added.

**Feasibility verdict:** **Medium as support paper / low-medium as standalone.**

**Best use:**

- Supplementary driver map for lineage-axis paper family.
- A methods/resource table for multi-driver risk stratification.
- Not a main “advance paper” unless Korean or external mutation-rich cohort is obtained.

**Main risk:** BRAF/RAS/TERT stratification is heavily published and may feel incremental.

---

### 3.3 Korean NRG1 germline risk variant x somatic profile

**User idea:** Korean germline NRG1 risk variant plus somatic BRAF/RAS/TERT/RET pattern; outcome, age, multifocality.

**Existing evidence:**

- There is already a dedicated plan: `project/results/audit_2026_04_29/nrg1_plan.md`.
- That plan explicitly says:
  - no Korean germline genotype data in project;
  - no matched Korean germline + somatic data;
  - analysis blocked until KoGES / Bundang / dbGaP-like controlled access.

**Feasibility verdict:** **Scientifically interesting but blocked now.**

**Can do now:**

- Literature registry.
- Data-access checklist.
- SNP/haplotype extraction plan.

**Cannot do now:**

- Germline x somatic association.
- Korean-specific risk modeling.
- NRG1 effect on age/multifocality/outcome.

**Decision:** Defer until real germline data arrives.

---

### 3.4 RET/ALK/NTRK/fusion full MAPK driver map

**User idea:** Full driver classification including kinase fusions and targetable alterations.

**Existing evidence:**

- Fusion overlay tables exist in v17 and audit folders.
- Counts are very small in TCGA:
  - RET_fusion: 1
  - NTRK_fusion: 1
  - ALK_fusion: 1
  - RET: 2
- Prior Paper 1 audit explicitly warned that fusion claims should not be headline because calls are incomplete / underpowered.

**Feasibility verdict:** **Good as a registry/resource layer; weak as a discovery analysis unless a proper fusion cohort is obtained.**

**Safe use:**

- “Full driver map” descriptive table.
- Candidate cohort registry.
- Targetable-driver landscape appendix.

**Unsafe use:**

- Fusion-driven subtype claim.
- Reflex fusion testing recommendation.
- RET/NTRK/ALK therapy-response implication without clinical response data.

---

### 3.5 Ethnicity / age / Korean vs Western genomic comparison

**User idea:** TCGA vs Korean cohorts; why Korean BRAF rate is high; methylation/environment signatures.

**Existing evidence:**

- Korean K2 / Yoo / GSE213647 related outputs exist, but not enough for controlled ancestry-causal genomic inference.
- NRG1 plan is blocked on germline.
- Methylation is deliberately excluded from Paper 1 strategy.

**Feasibility verdict:** **High 삽질 risk now.**

**Why risky:**

- “Why Koreans have more BRAF” needs strong germline/environment/ascertainment design.
- TCGA vs Korean differences can be confounded by histology mix, age, sequencing panel, sample selection, hospital referral bias.
- Methylation/environment causal claims are very hard without designed data.

**Possible safe version:**

- Descriptive registry of Korean thyroid cancer genomic cohorts.
- No causal ancestry/environment claim.

---

### 3.6 Metabolic reprogramming classification

**User idea:** OXPHOS vs glycolytic thyroid cancer states, maybe prognosis beyond BRAF.

**Existing evidence:**

- Paper 9 plan already exists and includes metabolic vulnerability candidates:
  - GLS
  - LDHA
  - MYC
  - OXPHOS/glycolysis compensation
  - NAMPT/NAD salvage
- Paper 1 mechanism figures already show pathway-level support, but not enough for a standalone metabolism paper.

**Feasibility verdict:** **Good future translational track, not immediate.**

**Best framing:** Paper 9 synthetic-lethal / vulnerability map after Paper 1 is accepted or at least stable.

**Why not now:**

- Requires CCLE/DepMap/PRISM/GDSC.
- Needs thyroid cell-line caveat handling.
- If started now, it will distract from Paper 1 and Paper 2.

---

## 4. Recommended advance-paper queue

### Immediate now

1. **Paper 2 HLA Pillar I v2**
   - Finish Korean baseline source/provenance extraction.
   - No forest until source table is stable.
   - This is the cleanest executable next step.

2. **Freeze Paper 1**
   - No more data.
   - Continue only web/writing/figure packaging.

### Next, after Paper 2 source table

3. **Dark-matter molecular taxonomy salvage memo**
   - Use existing `dark_matter_phase1/2` outputs.
   - Do not run new data first.
   - Decide whether it becomes Paper 1b, Paper 7 case study, or archive.

### Future only

4. **Paper 9 synthetic-lethal / metabolism**
   - Start only when Paper 1 is in print or stable enough to cite internally.

5. **NRG1 Korean germline x somatic**
   - Start only after germline data access.

---

## 5. Practical next command options

### Best next command

```text
extract Paper2 Korean baseline tables only
```

Reason: executable now, low-risk, no Paper 1 conflict, no raw data.

### If user wants advance-paper strategy instead

```text
write Dark Matter molecular taxonomy salvage memo from existing outputs only
```

Allowed scope: read existing `dark_matter_phase1/2` and `v17` outputs, no new GEO, no fusion overclaim, no survival headline.

### If user wants future translational track

```text
freeze Paper9 synthetic-lethal roadmap
```

Allowed scope: roadmap only, no DepMap/CCLE/PRISM execution yet.

---

## 6. Final decision

The listed ideas are not equally ready. The most executable next analysis is **Paper 2 HLA baseline extraction**. The most scientifically interesting archived advance-paper candidate is **BRAF/RAS-negative molecular taxonomy**, but it must be reframed away from survival/fusion/therapy claims. NRG1 and ethnicity-causal work are blocked by data access. Metabolism is a good Paper 9 direction, not a current sprint.

No Paper 1 expansion recommended.
