# Paper 3 ICI — 12-Week Execution Plan

**Working title:** HLA loss, neoantigen architecture, and immune ecotypes define ICI vulnerability in molecularly dark thyroid cancer
**Track:** A (this document) — design only. Plan below describes Track B execution; Track B starts only after entry gates G1–G6 close.
**Author:** Seungho Cook
**Status date:** 2026-05-04

**Claim guard:** All weeks below assume "ICI vulnerability / readiness / immunogenomic prioritization" framing. If thyroid ICI-treated raw RNA-seq is acquired during the plan, framing may upgrade — only after re-decision with user, never silently.

---

## 0. Entry preconditions (G1–G6)

Plan does not start until:
- G1 Paper 1 bioRxiv submitted (target 2026-06-13)
- G2 Paper 2 Task A/B/C closed
- G3 TCGA dbGaP access decided
- G4 ≥3 of 6 scRNA cohorts accessible
- G5 ≥3 of 5 Tier 1/2 pan-cancer ICI cohorts accessible
- G6 User explicit "Track B 시작"

If G1/G2 not closed, the plan re-baselines its Wk 1 day to the post-ship date.

---

## 1. Week-by-week schedule

### Wk 1 — Dataset acquisition & ETL scaffolding

- Verify all `to_verify` accessions in `paper3_ici_dataset_registry.md` §1–§4. Produce `paper3_ici_dataset_registry_VERIFIED.tsv`.
- Submit dbGaP / EGA applications in parallel (Liu melanoma, Yoo 2019 if EGA, TCGA paired-BAM controlled if needed).
- Stand up `project/results/paper3_ici/` directory tree per registry / atlas / hla / neo / dial / integrated outputs.
- Reuse Paper 1 ETL outputs as read-only inputs (driver / fusion calls, BRS, TDS).
- Lock signature gene lists to actually-measurable subset per cohort. Produce `paper3_ici_signature_registry_VERIFIED.tsv`.

Deliverables: verified registries, directory scaffolding, dbGaP/EGA applications submitted.

### Wk 2 — Bulk preprocessing & cohort harmonization

- Re-pull / re-normalize bulk cohorts to a common log2(TPM+1) baseline (RNA-seq) and RMA log-intensity (microarray). Reuse Paper 1 pipelines where applicable; do not re-run if no scope difference.
- ComBat / ComBat-seq batch correction per cohort family with subtype as protected covariate; LODO sensitivity per `v52_lodo_finding`.
- QC pass: read count, mapping rate, dropout, gene detection.
- Compute all bulk signatures from `paper3_ici_signature_registry.md` §1–§7 on harmonized matrix; lock `paper3_ici_bulk_signature_scores.parquet`.

Deliverables: harmonized bulk matrix, QC tables, all-signature score matrix.

### Wk 3 — Bulk immune ecotype discovery (Module A part 1)

- NMF rank selection (K = 3..8) on signature-z-score matrix; cophenetic correlation, dispersion. Pick K.
- Consensus clustering as alternative; compare ARI with NMF.
- Stability: 100 bootstrap subsamples; per-sample assignment confidence.
- Headline figure 2A/B drafts.

Deliverables: ecotype labels, stability metrics, fig 2A/B drafts.

### Wk 4 — Bulk ecotype × clinical/molecular axes (Module A part 2)

- Cross-tabulate ecotype × subtype (PTC/PDTC/ATC), driver (BRAF/RAS/fusion/dark-matter), TDS tertile.
- Survival in aggressive subset (TCGA stage III/IV + GSE76039 PDTC/ATC) with Cox.
- Sensitivity to scoring method (ssGSEA vs Singscore) and to NMF seed.
- Lock `paper3_ici_bulk_immune_ecotype_summary.tsv`.
- Fig 2C/D/E/F drafts.

Deliverables: ecotype summary table, fig 2 first complete draft.

### Wk 5 — scRNA atlas integration (Module B part 1)

- Per-cohort QC, doublet removal, normalization.
- HVG selection, scVI training, scANVI fine-tune, Harmony cross-check.
- kBET / LISI / iLISI metrics; integration acceptance per `paper3_ici_scRNA_reference_atlas_summary.md` §3.
- Coarse Level 1 cluster + annotation (celltypist seeded).

Deliverables: integrated AnnData, integration QC, Level 1 annotations.

### Wk 6 — scRNA sub-state mapping & signature lock (Module B part 2)

- Level 2 sub-state annotation per cell type.
- DE per sub-state; lock `sc_signatures_locked.json` per registry §8.
- inferCNV / CopyKAT for malignant epithelial calling.
- Per-patient state frequency table (`state_frequency_per_patient.tsv`).
- Bulk projection of sc-derived signatures via ssGSEA + CIBERSORTx cross-check.
- Fig 3 drafts.

Deliverables: locked sc signatures, per-patient state frequencies, fig 3 first complete draft.

### Wk 7 — HLA typing (Module C part 1)

- TCGA-THCA WES BAM via GDC controlled access — OptiType + Polysolver Class I; xHLA / HLA*LA Class II.
- arcasHLA on RNA-seq for cohorts without WES.
- Lock `hla_class1_calls.tsv`, `hla_class2_calls.tsv`.
- Allele frequency comparison vs Korean baseline (PRJEB11591) and TCGA pan-Asian/Caucasian breakdown.

Deliverables: locked HLA call tables, allele-frequency tables.

### Wk 8 — HLA LOH + Neoantigen prediction (Module C part 2)

- LOHHLA on TCGA-THCA paired tumor-normal BAM with purity ≥0.4 inclusion.
- Mutect2 paired re-run only if MC3 MAF inadequate; otherwise reuse MC3.
- NetMHCpan 4.1 Class I + NetMHCIIpan 4.0 Class II prediction with RNA-expression filter (TPM ≥1 in tumor).
- Self-similarity filter against UniProt human reference.
- Lock `lohhla_calls.tsv`, `neoantigens_class1.tsv`, `neoantigens_class2.tsv`.
- HLA-intactness and neoantigen-presented composite per `paper3_ici_HLA_neoantigen_feasibility.md` §8.
- Fig 4 drafts.

Deliverables: HLA LOH + neoantigen tables, HLA-intact + neo+ subgroup roster, fig 4 first complete draft.

### Wk 9 — DIAL audit (Module D)

- Pan-cancer ICI cohort signature scoring on the same registry.
- Logistic regression (response endpoint) and TIS-regression (T-cell-inflamed endpoint) per cohort.
- Bootstrap × 1000, permutation × 1000 for sign + effect size CI.
- Tissue-transfer check on TCGA aggressive subset + GSE76039 (TIS endpoint).
- Lock `dial_verdict.tsv`.
- Fig 5 drafts.

Deliverables: DIAL verdict per signature, fig 5 first complete draft.

### Wk 10 — Integrated readiness score (Module E)

- Define `ICI_VULN_DIAL` using DIAL-passing signatures only (default = equal weight on survivors).
- Sensitivity weightings (a/b/c/d) per signature registry §9.
- Stratify score by subtype, driver, TDS, ecotype, HLA-intact + neo+ subgroup.
- Pan-cancer transfer validation — predict response in IMvigor210 / Hugo / Riaz with the *thyroid-adjusted* score (sanity).
- Spatial overlay (read-only from Paper 1/2 spatial pipeline).
- Lock `integrated_readiness_score.tsv`.
- Fig 6 drafts.

Deliverables: readiness score per sample, fig 6 first complete draft.

### Wk 11 — Manuscript drafting (figures locked)

- Figure captions finalized; supplementary figures S1–S14 drafted.
- Methods written from M-tier scaffolding (mirroring Paper 1's `methods_M1_M11_scaffold.md` style).
- Results sections following figure order: F1 cohort → F2 ecotype → F3 sc atlas → F4 HLA/neo → F5 DIAL → F6 integrated score.
- Discussion drafted by user (voice-protected; per `v17_sprint_vs_marathon_violation.md` no Claude generation in voice-protected sections; Claude provides scaffolding).
- Limitations section explicit on thyroid ICI raw RNA-seq gap.
- References built on shared `2026_05_03_references.bib` extended with ICI literature.

Deliverables: manuscript v1 with locked figures + methods + results; user-drafted discussion outline.

### Wk 12 — Internal review & revision

- User-led discussion completion.
- Internal QA: claim guard sweep — flag any "response predictor" language; rewrite to "vulnerability / readiness".
- Cross-paper boundary sweep — flag any prose that overlaps Paper 1 (DM1 cluster) or Paper 2 (PTC+HT BCR/TLS) or Paper 4 (GD HLA forest).
- Cohort and accession verification refresh — every accession in dataset registry confirmed final.
- Reviewer Q anticipation list (mirroring Paper 1's `2026_05_03_reviewer_QA_consolidated.md`).
- Cover letter draft (mirroring Paper 1's `2026_05_03_cover_letter_assembly.md`).

Deliverables: manuscript v2 ship-ready, reviewer Q list, cover letter draft.

---

## 2. Critical-path & parallelism

Critical path: Wk 1 (data) → Wk 2 (preprocess) → Wk 3–4 (Module A) → Wk 7–8 (Module C) → Wk 9 (Module D) → Wk 10 (Module E) → Wk 11 (manuscript).

Parallel branches:
- Wk 5–6 (Module B scRNA) runs alongside Wk 3–4 / Wk 7–8 because it does not block Module A or C.
- DIAL pan-cancer ICI ETL can begin Wk 5 in parallel with Module B.
- Manuscript scaffolding (Methods / supp tables) can begin Wk 8 in parallel with Module D.

---

## 3. Compute budget (rough estimate)

| Module | Compute | Notes |
|---|---|---|
| Bulk preprocess + signatures | low (~50 CPU-hr) | mostly normalization. |
| NMF + bootstrap | low | |
| scRNA atlas integration | medium (1–2 GPU-day on A100; or CPU-only ~2–3 days) | scVI on 100–300K cells. |
| LOHHLA + Mutect2 | medium-high (~200–400 CPU-hr) | biggest variable. |
| NetMHCpan / NetMHCIIpan | low-medium | per-allele × per-peptide. |
| DIAL audit | low | bootstrap dominates. |

Total estimate: ~1500–2500 CPU-hr + ~2 GPU-days. Single Azure NC-series burst over 1–2 weeks suffices; the `v17_arcasHLA_korean_k2` $4.80 burst pattern generalizes.

---

## 4. Risk & contingency triggers

| Trigger | Contingency |
|---|---|
| dbGaP TCGA paired BAM not granted by Wk 4 | Drop LOHHLA from Module C; report HLA call frequencies + MAF-based neoantigen only; mark fig 4B as future-work. |
| <3 scRNA cohorts accessible | Run Module B on best single cohort; fig 3 narrows; sc signature lock proceeds with caveat. |
| <3 ICI reference cohorts accessible | DIAL audit narrows; fig 5 reduced; if <2, defer DIAL to follow-up paper, drop fig 5, retain Module E with unaudited signatures + explicit caveat. |
| MC3 MAF disagreement with Paper 1 driver calls | Halt Module C until reconciled with Paper 1 ETL. |
| Paper 1 ship slips past late August 2026 | Push Track B start; do NOT begin compute before Paper 1 ship. |
| User wants to claim ICI response prediction in thyroid mid-plan | Hard stop; re-decide framing; do not silently proceed. |

---

## 5. Manuscript venue ladder (target ladder, decided at Wk 10)

- **Reach** — *Nature Cancer* / *Cancer Cell* — only if DIAL produces a thyroid-specific failure-mode story AND HLA-intact + neo+ + dedifferentiated subgroup is well-defined AND spatial cohort validates TLS/CXCL13 niche.
- **Default base** — *Nature Communications* / *Cell Reports Medicine* / *JCI Insight* — the registered design supports this tier even with degraded scope.
- **Safe** — *npj Precision Oncology* / *Genome Medicine* — fallback if multiple G-gates fail.

Avoid Frontiers/MDPI per `v17_graves_pivot`. ASCO / SITC abstract spin-off optional after submission.

---

## 6. Deliverables checklist (Track B end of Wk 12)

- [ ] `paper3_ici_dataset_registry_VERIFIED.tsv`
- [ ] `paper3_ici_signature_registry_VERIFIED.tsv`
- [ ] `paper3_ici_bulk_immune_ecotype_summary.tsv`
- [ ] `scrna_atlas/atlas_obs.parquet`, `atlas_markers.tsv`, `sc_signatures_locked.json`, `state_frequency_per_patient.tsv`
- [ ] `hla/hla_class1_calls.tsv`, `hla_class2_calls.tsv`, `lohhla_calls.tsv`
- [ ] `neo/neoantigens_class1.tsv`, `neoantigens_class2.tsv`
- [ ] `dial/dial_verdict.tsv`, `audit_report.md`
- [ ] `integrated/integrated_readiness_score.tsv`, `hla_intact_neoag_pos_subgroup.tsv`
- [ ] Manuscript v2 + supp + cover letter draft

---

## 7. Track A → Track B handoff

This 12-week plan, the registries, the atlas plan, the HLA/neoantigen feasibility, the DIAL plan, the figure plan, and the go/no-go verdict are the complete design freeze. Track B begins by re-reading these documents, executing G1–G6 closure, and starting Wk 1.

---

Track A completed. No marathon violation.
