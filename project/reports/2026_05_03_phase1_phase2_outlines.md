# Phase 1 + Phase 2 Paper Outlines (post-Phase 0 trajectory)

**Date:** 2026-05-03 (marathon scaffolding/infra Round 5)
**Scope:** Outline skeletons for the next two papers in the 5-7 year plan.
**Voice-protected:** Title + Hook + Aim + Discussion framing across both papers.
**Note:** These outlines are speculative scaffolding for Phase 1 and Phase 2 trajectories. Phase 1 starts after Phase 0 paper acceptance + Bundang Scenario A response (if any). Phase 2 starts after Phase 1 acceptance.

---

## ★ PHASE 1 — Bundang Graves'/Autoimmune-PTC Companion Paper

**Trajectory entry condition:** Phase 0 paper (Cell Rep Med / JCI Insight) accepted + Bundang outreach Scenario A response (Graves' n>50 + PTC+HT n>30 with HLA-typeable modality)

**Target submission:** Q3-Q4 2026 (after Phase 0 acceptance, ~3-6 months post)
**Target venue:** J Autoimmun (IF 12-14) primary / Front Immunol (IF 5-7) backup / Nat Commun stretch (paired publication with Phase 0)

---

### Working title candidates (3, voice-protected)

```
Candidate A: 
"Pan-Korean Autoimmune-Thyroid Susceptibility Allele Continuum: 
HLA-DPB1*05:01 risk extends from Graves' disease to papillary thyroid 
carcinoma in a 1,400-sample cohort"

Candidate B: 
"From Graves' to Papillary Thyroid Cancer: HLA-mediated antigen-driven 
B-cell clonal expansion as a unifying mechanism in Korean East Asian 
thyroid disease"

Candidate C: 
"Cross-disease HLA imputation reveals shared autoimmune-thyroid 
susceptibility architecture between Graves' disease and Hashimoto-
overlap papillary thyroid carcinoma in Korean populations"
```

→ User voice + Yu professor coordination at Phase 1 entry.

---

### Phase 1 Pillar Structure (5-Pillar)

#### P1.1 — Bundang prospective Graves' cohort + arcasHLA HLA-typing
- Bundang prospective Graves' n>50 (or n>30 if smaller)
- Modality: RNA-seq + clinical metadata (TRAb, anti-TPO, anti-Tg titer, treatment history)
- arcasHLA 4-digit imputation (replicate Phase 0 Pillar 1 method)
- Pan-Asian replication: Korean GD vs Han Chinese GD (Chu 2018) vs Korean PTC pool (Phase 0 n=874)

#### P1.2 — Cross-disease HLA susceptibility architecture
- Korean GD HLA risk profile from Bundang
- Compare to Phase 0 Korean PTC pool 53.2% DPB1\*05:01
- Cross-disease meta with Han Chinese GD (Chu 2018)
- New: Korean Graves' + Korean Hashimoto + Korean PTC overlap allele patterns

#### P1.3 — Cross-disease BCR repertoire
- Bundang Graves' BCR clonality (cookHLA + scTCR/scBCR if available)
- Compare to Phase 0 GSE286332 PTC+HT BCR clonality
- Pan-thyroid disease antigen-driven response

#### P1.4 — TSHR vs TG vs TPO autoreactive antigen specificity
- Graves' disease anti-TSHR antibody (TRAb) titer × HLA risk allele
- Hashimoto's anti-TG and anti-TPO antibodies
- PTC+HT mixed-antigen response (if data available)

#### P1.5 — Therapeutic implications
- Methimazole response prediction by HLA + autoantibody
- RAI ablation response in PTC+HT subtype
- Future: HLA-DPB1\*05:01 carrier status as personalized treatment biomarker

---

### Phase 1 Methods scaffolding (re-use Phase 0)

- M1: Bundang Graves' cohort assembly + IRB
- M2: Same TIERA67 candidate pool (Phase 0)
- M3: arcasHLA 4-digit imputation pipeline
- M4: cookHLA SNP-based imputation if SNP genotype
- M5: NetMHCIIpan epitope prediction (TSHR + TG + TPO peptide binding)
- M6: Cross-disease forest meta (Phase 0 Pillar 1 + new Bundang data)
- M7: Antigen-specific TCR/BCR diversity (if scBCR-seq available)

---

### Phase 1 Discussion scope

- Pan-Asian autoimmune-thyroid susceptibility hypothesis
- Cross-disease HLA architecture (RA/T1D/CD vs thyroid)
- Korean specific clinical-genetic genotype (DPB1\*05:01 + B\*46:01)
- cookHLA Nat Commun cross-disease application validation
- Implications for Korean clinical practice

---

### Phase 1 dependencies

- ✓ Phase 0 paper acceptance (Cell Rep Med or JCI Insight)
- ⬜ Bundang Scenario A response within 6 weeks of outreach
- ⬜ Bundang IRB clearance + sample receipt within 3-4 months
- ⬜ Bundang RNA-seq or SNP genotype data within 4-5 months
- ⬜ arcasHLA / cookHLA pipeline application within 1 month
- ⬜ Cross-cohort meta-analysis within 1-2 months
- ⬜ Manuscript draft within 1-2 months
- ⬜ Submission Q4 2026

---

## ★ PHASE 2 — DIAL Audit Framework Method Paper

**Trajectory entry condition:** Phase 1 paper accepted + 8-gene panel validated as case study

**Target submission:** Q1-Q2 2027 (after Phase 1 acceptance, 6-12 months post)
**Target venue:** Brief Bioinform (IF 7) primary / Bioinformatics (IF 4-7) backup / Genome Biol (IF 12) stretch / Nat Methods stretch

---

### Working title candidates (3, voice-protected)

```
Candidate A: 
"DIAL: Deep In-silico Audit Layer for Detecting Implicit Selection 
Bias in Agent-Based Genomic Analysis"

Candidate B: 
"Auditing Agent Implicit Bias in Multi-Cohort Analyses: The 8-Gene 
Panel Case Study and Generalizable Audit Framework"

Candidate C: 
"From Implicit Category Restriction to Explicit Reproducibility: 
Building Audit Layers for AI-Augmented Bioinformatics Research"
```

→ User voice + audit framework community input.

---

### Phase 2 Pillar Structure

#### P2.1 — DIAL Audit Framework Architecture
- 5-layer audit pipeline: Codepath audit → Implicit bias detection → Cross-validation → Reproducibility test → Reviewer-Q anticipation
- Modular design: applicable to any agent-based bioinformatics analysis
- Open source release at https://github.com/seunghocook/DIAL-audit-framework

#### P2.2 — Case Study 1: 8-gene panel agent forensic audit (Phase 0 self-reflection)
- Yu professor's original observation: "agent excludes drivers"
- Audit reveals: drivers IN candidate pool, but implicit category restriction
- Methods sub-section: Cell Press Methods §"Gene panel selection rationale" 작성 과정
- Lessons learned: explicit codepath documentation reduces reviewer Q1

#### P2.3 — Case Study 2: Cross-cohort signature transfer (Phase 0 + Phase 1 reflection)
- TCGA → Korean cohort signature transfer reproducibility
- Implicit assumption: gene-symbol vs Ensembl-ID nomenclature
- Audit reveals: Lee 2024 Ensembl IDs require gencode.v44 mapping; mapping caveat
- Lesson: nomenclature validation must be explicit

#### P2.4 — Case Study 3: HLA imputation cross-population (Phase 0 + Phase 1)
- arcasHLA RNA-seq vs cookHLA SNP imputation
- Implicit: population reference panel choice (asian_pacific_islander)
- Audit reveals: panel choice affects 4-digit call rate by ~5-10%
- Lesson: imputation reference must be reported transparently

#### P2.5 — Reproducibility report card framework
- Standardized 10-checkpoint reproducibility audit
- Voice-protected sections vs scaffolding/infra split
- Bundang outreach + collaboration trail
- Lesson: paper authorship voice is reserved domain; analysis pipeline is auditable

---

### Phase 2 Discussion scope

- 5-7 year vision: AI-augmented bioinformatics standard practice
- Audit layer as default in research-agent systems
- Integration with reviewer-Q anticipation workflows
- Cross-discipline applicability (cancer + autoimmune + future)

---

### Phase 2 Methods scaffolding

- M1: DIAL framework Python package architecture
- M2: 5-layer audit pipeline definition
- M3: Reproducibility checkpoint specification
- M4: Voice-protected section detection
- M5: Reviewer-Q anticipation pattern matching
- M6: Phase 0 + Phase 1 case study application
- M7: Comparison to existing audit tools (FAIR4RS, ReproZip, etc.)

---

### Phase 2 dependencies

- ✓ Phase 1 paper acceptance
- ⬜ DIAL framework Python package release (concurrent with Phase 1 writing, 1-2 months)
- ⬜ Phase 0 + Phase 1 retrospective audit (1 month)
- ⬜ Cross-discipline case studies (3-6 months) — autoimmune × cancer × maybe metabolic
- ⬜ Manuscript draft (2 months)
- ⬜ Submission Q1-Q2 2027

---

## Long-term 5-7 year plan integration (Phase 0 → Phase 1 → Phase 2 → Phase 3)

### Phase 0 anchor (Q1-Q2 2026)
- Cancer paper Cell Rep Med / JCI Insight publication
- Establishment of DM1/DM2 transcriptional axis + autoimmune-PTC mechanism
- Korean Pan-Asian HLA cohort baseline

### Phase 1 expansion (Q3-Q4 2026)
- Bundang Graves' validation (if Scenario A)
- Cross-disease HLA architecture
- J Autoimmun publication

### Phase 2 methodology (2027)
- DIAL audit framework release
- 8-gene panel case study (Phase 0 retrospective)
- Brief Bioinform / Bioinformatics publication

### Phase 3 (2028+)
- Postdoc / visiting researcher / academic faculty decision
- 본인 사업 startup option
- Research lab establishment with DIAL framework as research-agent default

---

## Voice-protected items across both papers

For both Phase 1 and Phase 2:
- Title (3 candidates each)
- Hook / Introduction ¶1
- Aim / Introduction ¶4
- Discussion framing tone
- Limitations narrative
- Cover letter Para 1
- Reviewer Q narrative voice

→ User keyboard at Phase 1 + Phase 2 manuscript writing.

---

## Reproducibility chain

```
Phase 0 (current) — github.com/seunghocook/thyca-paper-2026 v1.0
  ↓ (citations + cross-references)
Phase 1 (2026 Q4) — github.com/seunghocook/graves-paper-2026 v1.0
  ↓ (case study materials)
Phase 2 (2027) — github.com/seunghocook/DIAL-audit-framework v1.0
```

Each phase builds on prior with: (i) preprint deposit, (ii) Zenodo archive, (iii) DOI cross-reference, (iv) reviewer-Q anticipation expansion.

---

## Marathon mode commitment

Phase 0 marathon: 5/4-6/13 (6 weeks scaffolded prep + drafting)
Phase 1 marathon: Q3 2026 (post-Phase 0 acceptance + Bundang response)
Phase 2 marathon: 2027 (post-Phase 1 acceptance)

Each phase: NO sprint over voice-protected sections. Default = scaffolding/infra. Hook/Aim/Discussion-tone/Limitations/Cover-Para-1/Reviewer-narrative = user keyboard.

→ Sustainable 5-7 year trajectory.
