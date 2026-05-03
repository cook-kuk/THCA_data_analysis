# Paper 3 ICI — HLA & Neoantigen Feasibility

**Working title:** HLA loss, neoantigen architecture, and immune ecotypes define ICI vulnerability in molecularly dark thyroid cancer
**Track:** A — design only. No HLA typing runs, no LOHHLA execution, no NetMHCpan invocation. arcasHLA results from `v17_arcasHLA_korean_k2` are referenced read-only and **not** re-run.
**Author:** Seungho Cook
**Status date:** 2026-05-04

**Claim guard:** HLA loss + neoantigen load are reported as *immunogenomic features*, not as standalone "ICI response signals". Combined into the integrated readiness score (`paper3_ici_signature_registry.md` §9) only after DIAL audit.

---

## 1. Per-cohort feasibility matrix

| Cohort | RNA-seq | WES / WGS | HLA Class I | HLA Class II | LOHHLA-feasible | Neoantigen-feasible | Notes |
|---|---|---|---|---|---|---|---|
| TCGA-THCA | yes | yes (BAM via GDC) | OptiType / Polysolver from BAM | arcasHLA from RNA BAM | yes (paired tumor-normal BAM) | yes (MAF + RNA expression filter) | Primary neoantigen + HLA LOH cohort. |
| GSE76039 PDTC/ATC (Landa 2016) | yes | yes (Landa 2016 published WES) | OptiType from RNA BAM (no germline normal) | arcasHLA | partial — paired normal availability per-sample needed | yes if MAF + paired normal | Verify normal availability; without it, somatic calls less reliable. |
| GSE33630 / GSE65144 / GSE54958 / GSE53157 / GSE29265 | microarray | no | no | no | no | no | Microarray — no HLA, no neoantigen. Excluded from Module C. |
| Yoo SK 2019 Korean ATC | yes (RNA-seq) | likely yes (matched WES) | arcasHLA | arcasHLA | conditional on access mode | yes | EGA-controlled candidate; verify at Track B kickoff. |
| PRJKA210106 Korean PTC n=282 | yes (RNA-seq) | unclear | arcasHLA | arcasHLA | likely no (no paired normal) | no (no MAF) | HLA frequency anchor only; cross-paper boundary care. |
| PRJEB11591 Korean PTC n=260 | yes | no | arcasHLA already done | arcasHLA already done | no | no | Read-only reference per `v17_arcasHLA_korean_k2`. **Do not re-run.** |
| Pan-cancer ICI cohorts (DIAL audit set) | yes | partial | per-cohort published HLA | per-cohort | per-cohort | per-cohort | Used as DIAL anchors, not as thyroid HLA evidence. |

**Effective Module C cohort = TCGA-THCA + GSE76039 + (Yoo 2019 if accessible).**

---

## 2. Tool selection (locked at design)

| Step | Primary | Sensitivity / fallback | Why |
|---|---|---|---|
| HLA Class I from WES/WGS BAM | OptiType (paired-end exome) | Polysolver | OptiType is high-accuracy on Class I; Polysolver published for cancer cohorts. |
| HLA Class II from WES | xHLA / HLA*LA | — | WES Class II coverage is uneven; require ≥30× over HLA-DRB1, HLA-DPB1, HLA-DQB1 exons. |
| HLA from RNA-seq | arcasHLA | seq2HLA | arcasHLA already validated by us on PRJEB11591 (`v17_arcasHLA_korean_k2`). Default for cohorts without WES. |
| HLA SNP-imputed 4-digit | cookHLA (Cook et al. *Nat Commun*) | — | Reserved for cohorts with SNP arrays. **Cookhla scope is Paper 4 (GD HLA backlog).** Do not introduce here unless absolutely necessary. |
| HLA LOH | LOHHLA | HLAthena (rare) | Requires paired tumor-normal BAM, allele-specific. |
| Somatic mutation calling | Mutect2 (single-sample mode if no normal); paired Mutect2 if normal exists | Strelka2 | Single-sample mode flagged; PoN required. |
| Variant filtering | dbSNP common-variant filter + gnomAD AF<1e-4 | — | |
| Neoantigen prediction (MHC-I) | NetMHCpan 4.1 | MHCflurry 2.0 | Standard. |
| Neoantigen prediction (MHC-II) | NetMHCIIpan 4.0 | — | Class II is noisier; report binders ≤500 nM with annotation that values are less reliable. |
| Neoantigen RNA-expression filter | TPM ≥1 in tumor RNA-seq for the source gene | — | Required to avoid silent-allele false positives. |
| Self-similarity filter | NetMHCpan against UniProt human reference; remove peptides identical to self | — | |
| Pipeline orchestration | nf-core/epitopeprediction or pVACseq | — | Pick at Wk 1 of Track B based on cohort BAM availability. |

---

## 3. Computational scope (estimate, NOT executed)

Per cohort, single-machine cost estimate (Azure NC-series or local HPC):

| Module | TCGA-THCA (~500) | GSE76039 (~37) | Yoo 2019 (~? aggressive) |
|---|---|---|---|
| HLA typing (RNA-seq + WES) | ~30 CPU-hr | ~3 CPU-hr | ~5 CPU-hr |
| LOHHLA (paired tumor-normal) | ~50 CPU-hr | ~5 CPU-hr | ~5 CPU-hr |
| Mutect2 (paired) | ~200 CPU-hr if BAMs already aligned | ~20 CPU-hr | ~30 CPU-hr |
| NetMHCpan + RNA-expression filter | ~10 CPU-hr | ~1 CPU-hr | ~2 CPU-hr |

**Total Module C estimate:** ~300–400 CPU-hr on TCGA-THCA dominated by Mutect2 if re-running from BAM. If re-using GDC MC3 MAF, Mutect2 cost drops to near-zero and Module C becomes ~50 CPU-hr.

**Decision for Track B:** Reuse GDC MC3 MAF for TCGA somatic calls. Run LOHHLA fresh because MC3 does not provide HLA LOH calls. Run OptiType/Polysolver fresh because per-sample HLA calls may not be packaged.

---

## 4. Headline Module C analyses (defined now, executed later)

1. **HLA typing & frequency** — per cohort, allele frequencies at 4-digit Class I/II. Compare to Korean baseline (PRJEB11591) and TCGA-THCA pan-Asian/Caucasian breakdown.
2. **HLA LOH frequency** — TCGA-THCA paired BAM via LOHHLA. Stratify by BRAF / RAS / fusion / dark-matter / TDS tertile.
3. **Neoantigen load** — per-sample SNV-derived neoantigen count (Class I, Class II), expression-filtered, self-filtered.
4. **Clonal vs subclonal neoantigens** — using PyClone-VI or sciClone CCF estimates if available.
5. **Driver-derived neoantigens** — restrict to BRAF V600E, RAS hotspots, TERT promoter mutations (non-coding so neoantigen unlikely; reported separately), TP53 hotspots — per known immunogenicity priors.
6. **HLA-intact + neoantigen-positive subgroup** — defined as no HLA LOH × neoantigen load ≥ cohort median. Cross-tabulate with bulk immune ecotype (Module A) and dedifferentiation tertile.
7. **PTC vs PDTC vs ATC neoantigen burden** — expected: ATC > PDTC > PTC in TMB; test whether HLA presentation (LOH frequency) co-evolves.

---

## 5. Risks & mitigations

- **R1: TCGA-THCA paired-normal BAMs are slow to retrieve.** Mitigation — request GDC dbGaP-controlled access at Track B Wk 1; in the interim, design and dry-run pipeline against open-access TCGA tutorial BAM.
- **R2: GSE76039 lacks paired normals for some samples.** Mitigation — drop those from LOHHLA + neoantigen; report HLA typing + somatic-MAF-based neoantigen only with caveat.
- **R3: Yoo 2019 may be EGA-controlled.** Mitigation — apply at Wk 1; if not granted by Wk 8, frame as planned validation, drop from primary Module C.
- **R4: Class II prediction is unreliable.** Mitigation — report Class II as descriptive, highlight Class I in headline; Class II analyses framed as exploratory.
- **R5: TMB in PTC is low.** Many PTC samples will have <5 SNV-derived neoantigens. Statistical handling — model as Poisson; report rate, not naive count comparisons.
- **R6: HLA LOH calls noisy in low-purity samples.** Mitigation — purity ≥ 0.4 inclusion gate; sensitivity at 0.3 reported.
- **R7: cookHLA is the user's first-author tool but scope is Paper 4.** Mitigation — do not introduce cookHLA here. Class I from OptiType/Polysolver/arcasHLA only.

---

## 6. Cross-paper boundary (binding)

- **Paper 1 (DM1 molecular)** — Paper 1 does not include HLA typing or neoantigen prediction. Module C is fully owned by Paper 3.
- **Paper 2 (H&E-DM1 / HT-overlap PTC)** — Paper 2's Pillar 1 forest is HLA *susceptibility* allele frequency, Korean PTC vs Korean baseline. Different question (susceptibility ≠ tumor HLA loss). Paper 3 Module C analyzes *somatic* HLA loss in tumors and *neoantigen* presentation; explicitly distinct from Paper 2 Pillar 1.
- **Paper 4 (GD HLA backlog)** — cookHLA SNP imputation belongs to Paper 4. Paper 3 does not run cookHLA.

---

## 7. Outputs (filenames to be produced in Track B Module C)

- `project/results/paper3_ici/hla/hla_class1_calls.tsv`
- `project/results/paper3_ici/hla/hla_class2_calls.tsv`
- `project/results/paper3_ici/hla/lohhla_calls.tsv`
- `project/results/paper3_ici/neo/neoantigens_class1.tsv`
- `project/results/paper3_ici/neo/neoantigens_class2.tsv`
- `project/results/paper3_ici/integrated/hla_intact_neoag_pos_subgroup.tsv`

---

## 8. Module C → ICI vulnerability score linkage

- `HLA_INTACTNESS = 1 − HLA_LOH_indicator × allele_lost_fraction`
- `NEO_LOAD_LOG = log10(neoantigen_count_class1 + 1)` (expression-filtered)
- `NEO_PRESENTABLE = NEO_LOAD_LOG × HLA_INTACTNESS`

These three feed Module E's integrated score with weight calibration deferred to Track B.

---

Track A completed. No marathon violation.
