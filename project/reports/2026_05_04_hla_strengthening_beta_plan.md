---
title: "HLA strengthening β plan — Tier 1 + Tier 2 spec (HT-isolated, no Paper 4 / Paper 3 / large new analysis)"
date: 2026-05-04
selected_option: β (from `2026_05_04_HLA_strengthening_situation.md` §5)
parent_audit: project/reports/2026_05_04_HLA_strengthening_situation.md
status: SPEC/PLAN ONLY — no analysis execution, no file edits, no manuscript prose, no voice-protected writing
scope_lock: Paper 2 (Hashimoto-overlap PTC) HLA work strengthening only
forbidden: Paper 4 GD analysis · Chu 2018 GD reactivation · Bundang Graves · cookHLA SNP integration · new large data acquisition · pathology image analysis · voice-protected Discussion prose · Paper 3 touch
---

# HLA strengthening β plan

---

## 1. Executive decision

**β selected**: Tier 1 (immediate execution path) + Tier 2 (spec only).

**Goal**: Strengthen Paper 2 HLA work without entering Paper 4/GD territory and without launching new large analyses.

**Allowed scope**:
- A2 Lee 2014 baseline frequency lookup PLAN + table extraction PREP
- A3 AFND South Korea baseline pool composition DECISION MEMO
- A1 Pillar I v2 6-allele forest EXECUTION PLAN
- B5 / B3 / B4 / C4 SPEC ONLY (no execution)

**Hard forbidden**: as stated in frontmatter.

**Marathon mode**: 5/4–6/13 active. No execution within this plan unless explicitly approved per §7.

---

## 2. Tier 1 — immediate plan

### 2.1 A2 — Lee 2014 baseline frequency lookup PLAN

**Objective**: Extract Korean general population HLA allele baseline frequencies for Pillar I v2 6-allele forest, specifically the 2 alleles currently missing from `p2_pillar1_forest_v2/` (C\*01:02, DQB1\*02:01).

**Source paper (verify exact citation at execution)**:
- Likely candidate: Lee MN et al. *Tissue Antigens* 2014 (KOTRY Korean donor cohort)
- Verification required (do NOT hardcode citation — avoid Krishnamoorthy/Chen misattribution repeat)
- Verification method: PubMed search "Lee Tissue Antigens 2014 Korean HLA donor" + cross-check against AFND South Korea entries

**Target fields per allele**:

| Field | Type | Purpose |
|---|---|---|
| `allele_4digit` | str | e.g., C\*01:02, DQB1\*02:01 |
| `carrier_count` | int | n carriers in Lee 2014 cohort |
| `n_total` | int | total Lee 2014 sample size for that allele class |
| `carrier_frequency` | float | proportion |
| `wilson_lower_95` | float | Wilson 95% CI lower |
| `wilson_upper_95` | float | Wilson 95% CI upper |
| `table_reference` | str | Lee 2014 Table N location |
| `doi_or_pmid` | str | citation anchor |
| `extraction_date` | str | 2026-05-04 |
| `extraction_method` | str | "manual_pdf" / "paperclip" / "AFND_cross_check" |

**Priority alleles**:

| Allele | Status in p2_pillar1_forest_v2 | Action |
|---|---|---|
| **C\*01:02** | ⏸ MISSING (Lee 2014 fill needed) | Extract |
| **DQB1\*02:01** | ⏸ MISSING (Lee 2014 fill needed) | Extract |
| DPB1\*05:01 | ✅ in v2 (verify against Lee 2014) | Cross-check only |
| A\*02:07 | ✅ in v2 (verify against Lee 2014) | Cross-check only |
| B\*46:01 | ✅ in v2 (verify against Lee 2014) | Cross-check only |
| DRB1\*07:01 | ✅ in v2 (verify against Lee 2014) | Cross-check only |

**Output schema (proposed)**:
- File: `project/results/p2_pillar1_forest_v2/lee2014_baseline_lookup.tsv` (NEW, on execution)
- Format: long-form, one row per allele
- 10 columns per the table above

**Execution decision**: AWAIT advisor approval (§5 Q1, Q2) + user `run A2 lookup only` command (§7).

### 2.2 A3 — AFND South Korea baseline pool composition DECISION MEMO

**Objective**: Decide Pillar I v2 §8 G2 — what sources compose the "AFND South Korea baseline pool"?

**Decision options** (for advisor):

| Option | Composition | Pros | Cons |
|---|---|---|---|
| **G2-1** | AFND South Korea raw entries only | Pure AFND, single-source provenance | Variable per allele class, AFND may have low-n entries |
| **G2-2** | AFND + KOTRY donor (Lee 2014) merged | Larger n, well-established Korean reference | Two-source heterogeneity to disclose |
| **G2-3** | AFND + KOTRY + KBP (Korean Bone Marrow Donor Program) | Largest n, broad Korean coverage | Three-source heterogeneity, reference selection bias |
| **G2-4** | AFND + KOTRY + KBP + cookHLA Korean reference panel | Largest n + cookHLA continuity | cookHLA panel uses SNP imputation (different platform); cross-platform mixing |
| **G2-5** | Lee 2014 only (single-source) | Cleanest attribution | Smaller n, paper-of-record dependence |
| **G2-6** | AFND only (single-source) | Pure database | Same as G2-1 |
| **G2-7** | Advisor specifies other | TBD | TBD |

**Claude recommendation (advisory only)**: **G2-2** (AFND + KOTRY/Lee 2014 merged) as middle ground — adequate n, two-source disclosure manageable, KOTRY = recognized Korean reference.

**Decision memo output**:
- File: `project/results/p2_pillar1_forest_v2/G2_baseline_pool_decision.md` (NEW, on advisor sign-off)
- Content: locked composition + n per source + heterogeneity disclosure + sensitivity plan

**Execution decision**: AWAIT advisor approval (§5 Q1, Q4).

### 2.3 A1 — Pillar I v2 6-allele forest EXECUTION PLAN

**Objective**: Extend existing 4-allele forest in `p2_pillar1_forest_v2/` to 6-allele forest by integrating A2 + A3 outputs.

**Required inputs**:

| Input | Source | Status |
|---|---|---|
| Korean PTC pool allele freq (4 alleles) | `p2_pillar1_forest_v2/korean_ptc_allele_freq.tsv` | ✅ exists |
| Korean PTC pool allele freq (C\*01:02, DQB1\*02:01) | Re-extract from arcasHLA outputs (no new imputation, existing data) | ⏸ pending |
| Korean baseline allele freq (4 alleles) | `p2_pillar1_forest_v2/korean_baseline_allele_freq.tsv` | ✅ exists (verify per G2 composition) |
| Korean baseline allele freq (C\*01:02, DQB1\*02:01) | A2 Lee 2014 lookup output | ⏸ pending |
| Statistical method | random-effects DerSimonian-Laird per allele OR | locked spec |

**Required outputs** (extending existing):

| Output | Format | Notes |
|---|---|---|
| `korean_ptc_allele_freq_v2_6alleles.tsv` | TSV | 6 rows, replaces current 4-allele |
| `korean_baseline_allele_freq_v2_6alleles.tsv` | TSV | 6 rows |
| `korean_PTC_vs_korean_baseline_forest_v2_6alleles.tsv` | TSV | Per-allele OR + 95% CI + Fisher p + BH FDR |
| `korean_PTC_vs_korean_baseline_forest_v2_6alleles.json` | JSON | Metadata + provenance |
| `forest_paper2_HT_only_6alleles.pdf` | PDF | Cell Press style 6-allele forest |
| `forest_paper2_HT_only_6alleles.png` | PNG | 300 dpi raster |
| `forest_v2_4allele_vs_6allele_comparison.md` | MD | Sensitivity contrast (4 vs 6); per-allele OR stability |
| `discussion_paragraph_v2_6alleles.md` | MD | TEMPLATE only with placeholders, NOT prose generation |

**Required script**:
- File: `project/notebooks_or_scripts/v17_paper2_pillar1_forest_v2_6alleles.py` (NEW on execution)
- Reuses existing v2 4-allele logic, extends to 6 alleles
- NO new data acquisition (A2 lookup is the only data input)

**Decision gates that must clear before A1 execution**:

| Gate | Status |
|---|---|
| G1 advisor approval of v2 spec | ⏸ pending |
| G2 AFND baseline pool composition (A3) | ⏸ pending |
| G3 Lee 2014 fill | ✅ YES locked 5/4 |
| G4 Harbin Korean fallback exclusion (sensitivity) | ⏸ pending |
| G5 Figure 1 re-render with 6-allele | ⏸ pending (post-execution) |
| G6 manuscript paragraph framing | ⏸ pending (post-execution) |

**Execution decision**: AWAIT G1 + G2 + G4 + user `prepare A1 forest script only` or `run full A1 forest` command (§7).

---

## 3. Tier 2 — spec only (no execution)

### 3.1 B5 — Antigen presentation gene score SPEC

**Objective**: Define a Paper 2 antigen presentation module score capturing HLA class II + class I + processing machinery activation in Hashimoto-overlap PTC.

**Gene set (proposed)**:

| Category | Genes | Rationale |
|---|---|---|
| HLA class II (primary) | HLA-DRA, HLA-DRB1, HLA-DRB3/4/5, HLA-DPA1, HLA-DPB1, HLA-DQA1, HLA-DQB1 | Paper 2 mechanism core (HLA-II–mediated) |
| HLA class II master regulator | CIITA | Drives HLA-II expression |
| Invariant chain | CD74 | HLA-II antigen loading chaperone |
| HLA class I (supportive) | HLA-A, HLA-B, HLA-C, B2M | Antigen presentation context |
| Antigen processing | TAP1, TAP2, ERAP1, ERAP2, CALR, CANX | Peptide loading machinery |
| Proteasome (immunoproteasome) | PSMB8, PSMB9, PSMB10 | IFN-γ-induced peptide processing |

**Method**:
- Per-sample z-score (within cohort) of each gene
- Module score = mean z across genes (class II module separate from class I)
- Optional: GSVA / ssGSEA per sample

**Cohort application**:

| Cohort | Application | Comparison |
|---|---|---|
| TCGA-THCA | Module score per sample | Paper 2 substrate cross-reference |
| Korean PTC pool n=874 | Module score per sample | Paper 2 main |
| GSE286332 | Module score per sample | PTC+HT vs PTC differential (already in P3 work) |
| Pu 2021 GSE184362 | Per-thyrocyte module score | sc validation (overlaps B3) |

**Output spec**:
- `project/results/p2_b5_antigen_presentation/per_sample_scores.tsv`
- `project/results/p2_b5_antigen_presentation/group_comparison.tsv` (DM1 vs DM2 + PTC+HT vs PTC)
- `project/results/p2_b5_antigen_presentation/B5_summary.json`

**Forbidden in B5**: Paper 4 GD context · Paper 3 ICI / HLA LOH / neoantigen · pathology image · voice-protected prose.

**Execution decision**: SPEC ONLY in this plan. Execution requires advisor approval + user explicit `execute B5` command.

### 3.2 B3 — sc HLA-II thyrocyte cluster analysis SPEC

**Objective**: Quantify HLA-II expression specifically in thyrocyte clusters from Pu 2021 + Lu 2023 single-cell datasets, in Hashimoto-overlap context.

**Datasets**:
- Pu 2021 GSE184362 (n=6 Fudan PTC patients) — already loaded for Paper 1 sc validation
- Lu 2023 GSE193581 (n=23) — already loaded for Paper 1 sc validation

**Analysis spec**:
1. Subset cells: thyrocyte markers (KRT8 + KRT19 + EPCAM positive); exclude immune (PTPRC+), stromal (PDGFRA+, COL1A1+)
2. Per-thyrocyte HLA-II module score (gene set from B5)
3. Per-patient pseudo-bulk: tumor thyrocyte mean HLA-II vs adjacent normal thyrocyte mean HLA-II
4. Compare to bulk DM1 score (Paper 1 cross-reference) — but DO NOT deep-dive DM1 mechanism (Paper 1 territory)
5. Hashimoto-overlap context: if patient-level Hashimoto status available in Pu 2021 metadata, sub-stratify

**Forbidden in B3**:
- ICI prediction (Paper 3)
- HLA LOH (Paper 3)
- Cancer driver mechanism deep dive (Paper 1)
- pathology image
- GD/Graves' context (Paper 4)

**Output spec**:
- `project/results/p2_b3_sc_hla2/per_thyrocyte_hla2_score.tsv`
- `project/results/p2_b3_sc_hla2/per_patient_pseudobulk_hla2.tsv`
- `project/results/p2_b3_sc_hla2/B3_summary.json`

**Execution decision**: SPEC ONLY. Execution requires user explicit command.

### 3.3 B4 — BCR/TLS × HLA-II crosstab SPEC

**Objective**: Cross-tabulate B cell receptor (BCR) clonality + tertiary lymphoid structure (TLS) score + HLA-II module activation per sample, in Hashimoto-overlap PTC context.

**Inputs (existing)**:

| Input | Source |
|---|---|
| BCR per-sample diversity | `project/results/d5p6_bcr_repertoire/per_sample_diversity.tsv` |
| TLS 12-gene score (Cabrita 2020) | `project/results/d5p6_bcr_repertoire/tls_score_per_sample.tsv` |
| HLA-II module score | B5 output (or recompute from existing v17_hla data) |

**Analysis spec**:
1. Per-sample matrix: BCR clonality × TLS score × HLA-II module
2. Spearman correlation matrix
3. Multivariate analysis: HLA-II as predictor of BCR clonality + TLS, controlling for sample n
4. Cohort: GSE286332 PTC+HT vs PTC (n=18) — primary
5. Optional: extend to Korean PTC pool n=874 if HLA-II module available per sample

**Forbidden in B4**:
- ICI response prediction (Paper 3)
- Cross-paper TLS-pan-cancer framing (Paper 3 territory)
- GD/Graves' context (Paper 4)
- voice-protected Discussion prose

**Output spec**:
- `project/results/p2_b4_bcr_tls_hla2/correlation_matrix.tsv`
- `project/results/p2_b4_bcr_tls_hla2/multivariate_summary.json`
- `project/results/p2_b4_bcr_tls_hla2/B4_summary.json`

**Execution decision**: SPEC ONLY. Execution requires user explicit command.

### 3.4 C4 — TCGA HLA imputation consistency audit SPEC

**Objective**: QC audit of TCGA-THCA HLA imputation consistency across methods (arcasHLA RNA-seq vs OptiType RNA-seq vs TCGA-released calls). NOT a methods paper analysis; QC for Paper 2 cohort confidence.

**Methods to compare**:
- arcasHLA from RNA-seq (Orenbuch et al.)
- OptiType from RNA-seq (Szolek et al.)
- TCGA Pan-Cancer Immune Project HLA calls (released)

**NOT included** (per "no cookHLA SNP integration" forbidden):
- ❌ cookHLA SNP imputation (requires SNP data acquisition — forbidden)
- ❌ HIBAG SNP imputation (same forbidden category)

**Concordance metrics**:
- Per-sample, per-gene 4-digit allele agreement rate
- Aggregate concordance per HLA gene (A, B, C, DRB1, DQB1, DPB1)
- Sensitivity by RNA-seq depth

**Cohort**: TCGA-THCA (n=504) only. Korean cohorts excluded from this QC (single-method arcasHLA only available).

**Forbidden in C4**:
- New SNP data download
- cookHLA / HIBAG SNP path
- Korean cohort cross-method comparison (single-method only)
- Paper 4 GD context

**Output spec**:
- `project/results/p2_c4_tcga_hla_qc/per_sample_concordance.tsv`
- `project/results/p2_c4_tcga_hla_qc/per_gene_concordance_summary.tsv`
- `project/results/p2_c4_tcga_hla_qc/C4_summary.json`

**Execution decision**: SPEC ONLY. Execution requires user explicit command.

---

## 4. Cross-paper boundary

### Paper 2 ALLOWED (this plan's scope)
- Hashimoto-overlap PTC mechanism
- HLA-II antigen presentation
- Korean PTC vs Korean baseline forest
- AFND South Korea baseline
- HT-specific tissue context
- BCR/TLS within Paper 2 scope
- mediation model
- destructive lymphocytic thyroiditis
- thyroid-lineage dedifferentiation

### Paper 4 FORBIDDEN (HARD STOP)
- Graves, GD, TSAb, TSI
- thyroid-stimulating antibody, hyperthyroidism, thyrotoxicosis
- exophthalmos, Graves' ophthalmopathy, TED, thyroid eye disease
- Bundang Graves cohort
- Chu 2018 GD reactivation (any direct comparison)

### Paper 3 FORBIDDEN (HARD STOP)
- ICI response prediction
- HLA LOH
- neoantigen
- DIAL audit
- pan-cancer ICI
- C5AR1 anti-PD-1 synergy
- Paper 3 ICI Track A bundle (chmod 444 FROZEN — do not touch)

### Paper 1 LIMITED (Paper 2 cross-reference only, no deep dive)
- DM1/DM2 cluster identity OK as substrate label
- BRAF/RAS/TERT/fusion mechanism deep dive forbidden
- TROP2, sacituzumab, H&E-DM1, image-DM1, TCGA WSI forbidden

### Forbidden categories outside paper boundaries
- cookHLA SNP integration (data acquisition forbidden)
- new large data acquisition (any direction)
- pathology image analysis (Paper 1 H&E-DM1 dropped; Paper 2 image scope undefined; HARD STOP this plan)
- voice-protected Discussion prose (always)

---

## 5. Advisor questions (Yu professor)

To bring to next Yu meeting:

1. **Q1 (G1)** — Is Pillar I v2 spec (Korean PTC vs Korean baseline) acceptable as Paper 2 Pillar I framing, replacing the deprecated v1 (vs Chu 2018 GD)?
2. **Q2 (G3 confirm)** — Should Lee 2014 fill be used to extend the forest from 4 alleles to 6 alleles? (G3 was YES locked 5/4 by user; this is advisor confirm)
3. **Q3 (G5)** — If 6-allele forest is generated, should it replace the existing 4-allele v2 in Paper 2 figures and Methods? (Or keep 4-allele as primary, 6-allele as supplementary?)
4. **Q4 (G2)** — AFND South Korea baseline pool composition: G2-1 / G2-2 / G2-3 / G2-4 / G2-5 / G2-6 (see §2.2 table). Recommendation G2-2 (AFND + KOTRY/Lee 2014). Advisor decides.
5. **Q5 (G4)** — Harbin Korean fallback entries in AFND (if exist): exclude in primary forest, retain in sensitivity, or include in primary?
6. **Q6 (Tier 2 prioritization)** — Which Tier 2 spec to advance first if execution capacity available: B5 (antigen presentation), B3 (sc HLA-II thyrocytes), B4 (BCR/TLS × HLA-II), C4 (TCGA HLA QC)?
7. **Q7 (timing)** — Pillar I v2 6-allele forest execution: now (after Q1-Q5 sign-off) or after Paper 1 bioRxiv?

**Advisor packet structure** (deferred until user confirms with `prepare advisor packet` command): 1-page bullets per Q1-Q7, no prose generation.

---

## 6. Execution gates

### What can run NOW (within this plan, no further approval needed)
- ✅ This spec/plan document (already written)

### What requires advisor (Yu professor) approval before execution
- A1 6-allele forest execution → Q1 + Q2 + Q4 + Q5 sign-off
- A2 Lee 2014 lookup execution → Q2 confirm + lookup approach approval
- A3 G2 baseline pool decision lock → Q4 sign-off
- B5/B3/B4/C4 spec → execution → Q6 prioritization

### What requires explicit user command (after advisor where applicable)
- `run A2 lookup only` → execute Lee 2014 baseline lookup
- `prepare A1 forest script only` → write Python script for 6-allele forest, NO execution
- `run full A1 forest` → execute 6-allele forest pipeline (after Q1+Q4 sign-off)
- `execute B5` / `execute B3` / `execute B4` / `execute C4` → Tier 2 spec → execution
- `prepare advisor packet` → 1-page bullet packet for Yu meeting

### What is DEFERRED (always)
- Voice-protected Discussion prose (本人 keyboard, always)
- Manuscript paragraph fill (post-execution + user voice)
- Figure 1 re-render decision (G5, post-execution)
- Paper 4 work (Bundang outreach + Chu 2018 reactivation — gating not met)
- Paper 3 work (Track A FROZEN; Track B BLOCKED)
- cookHLA SNP integration (forbidden in this plan)
- Pathology image analysis (forbidden in this plan)

---

## 7. Next command options

본인이 다음 중 specific 1개 (또는 조합) 명시:

- **`run A2 lookup only`** — execute Lee 2014 baseline frequency lookup for C\*01:02 + DQB1\*02:01 + cross-check 4 existing alleles. Output: `lee2014_baseline_lookup.tsv`. Requires advisor Q2 confirm OR user override.
- **`prepare A1 forest script only`** — write `v17_paper2_pillar1_forest_v2_6alleles.py` Python script for 6-allele forest, but DO NOT execute. Awaits all advisor gates (Q1+Q2+Q4+Q5).
- **`stop and bring to Yu meeting`** — freeze all execution. Package §5 Q1-Q7 advisor questions into 1-page bullet brief. No further action until Yu meeting + sign-off returns.
- **`advance B5 spec only`** / **`advance B3 spec only`** / **`advance B4 spec only`** / **`advance C4 spec only`** — extend the spec for one Tier 2 item with deeper detail (still spec only, no execution).
- **`prepare advisor packet`** — write `_paper2_HLA_advisor_packet_2026_05_04.md` with Q1-Q7 + 1-page bullets, ready for Yu meeting.
- **`other`** (specify)

**No analysis execution unless explicitly approved per command above.**

---

# END β PLAN — Marathon-compliant scaffolding only. Awaiting next command.
