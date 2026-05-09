# Submission readiness report

**Date**: 2026-05-07
**Status**: DRAFT after Phase 1-23 of "what matters" reframing.
**Recommendation**: **B+ → A-** depending on next round (see below).

---

## Headline grade per target journal class

| Target | Grade | Rationale |
|---|---|---|
| **NAR Database / Resource paper** | **A-** | Strongest fit. 13-source integration, leakage-aware benchmark, dbPepNeo rescue, public hub, reproducible pipeline. Needs final figure polish + supplementary tables. |
| **Briefings in Bioinformatics** | **A-** | Excellent fit for systematic-benchmark + leakage-audit framing. The "what actually matters" angle directly maps to BIB's audience. |
| **Cell Reports Methods / Patterns** | **B+** | Possible. Methods angle (strict-eval framework) is solid. ESM2 LOSO partial-fix needs scaling for higher impact. |
| **npj Digital Medicine** | **B** | Possible companion if KRAS prioritization + safety flags are framed as digital-twin-style design tool. Risky as primary venue. |
| **Cell Reports Medicine** | **B-** | Risky without wet-lab validation. The candidate prioritization is honest but still computational. |
| **Nature Cancer / Nat Mach Intell** | **C** | Not ready. Without GPU-scale ESM2-650M + full task-aligned features (NetMHCpan/MHCflurry under strict splits) + at least one wet-lab pilot, the novelty bar is too high. |

---

## What is ready (✅)

- 13-source integrated benchmark with `test_set_safety` annotation
- Leakage-aware audit (TRAINING_OVERLAP=61,656; HELD_OUT=15,765; EXTERNAL_TEST=1,197)
- Strict evaluation under random / peptide-GroupKFold / LOSO
- ESM2-35M LOSO benchmark (partial cross-source rescue)
- HLA confounding analysis (Simpson's paradox flag)
- Shortcut test suite (source-only / HLA-only / label-shuffle within-source)
- Candidate re-prioritization under strict-trained model
- Self-similarity autoimmunity safety flagging
- 6 RCSB pMHC + TCR-pMHC reference structures
- 290 VDJdb cross-annotation hits
- Reproducible scripts (`scripts/what_matters_*`)
- Frozen data snapshot manifest (sha256)
- Reviewer-2 attack list with 20+ pre-empted criticisms
- Manuscript draft (`manuscript/draft.md`) with HONEST framing
- Manuscript figures 1-6, 8 (Fig 7 architecture ablation deferred)

## What is risky (⚠)

- ESM2-35M LOSO improvement is modest (~+0.07 over biophys); reviewers may say "not transformative"
- Random-CV vs LOSO drop is the headline; depends on reviewer accepting that framing
- HLA confounding finding is partially sample-size dependent (stronger on full corpus, attenuated on strict subset)
- Architecture ablation limited to LR/RF (no MLP/attention/adversarial)
- No NetMHCpan/MHCflurry/PRIME baseline under strict splits
- Population coverage uses independence assumption
- KRAS candidates remain top but absolute scores recalibrated downward

## What must be fixed before submission (🔧)

1. Run NetMHCpan or MHCflurry under our strict splits — even a single tool would dramatically strengthen R6/R7. Wrappers in `scripts/` are the next code task.
2. ESM2-150M or ESM2-650M with GPU access (Azure burst ~$5-10) to test PLM-scaling as a primary experimental contribution.
3. Architecture ablation expansion: at minimum add a 2-layer MLP and a class-balanced LR; preferably add a source-adversarial head.
4. Tighten Figure 7 (interpretability triangulation) — currently deferred. Run permutation importance + ablation overlap on the strict-LR model.
5. Population coverage redo with IEDB population coverage tool.
6. Self-similarity safety filter expand to BLAST-against-human-proteome.
7. One or two figures need clearer captions / legend / pLDDT-color scale.

## What claims are allowed (✅)

- "Leakage-aware post-TESLA neoantigen vaccine benchmark"
- "Source-bias and HLA-composition confounding audit"
- "Strict-generalization benchmark under LOSO / peptide-GroupKFold"
- "Representation comparison: biophys vs PLM vs HLA shortcut"
- "Candidate prioritization with autoimmunity safety flags"
- "Pre-clinical validation roadmap"
- "Reproducible pipeline + public hub"
- "External TCR / VDJdb annotation layer"
- "Structural reference annotation layer (RCSB + ESMFold)"

## What claims must be removed (❌)

(See `manuscript/overclaim_audit.md` for full list of 128 flags found in earlier markdown — must be cleaned up before final draft.)

- "Validated vaccine"
- "Clinical efficacy"
- "Patient selection"
- "Predictor of response"
- "Confirmed immunogenic"
- "Causal"
- "Structure proves binding"
- "AlphaFold confirms"
- "Safe vaccine"
- "Universal vaccine"
- "Production-grade"
- "Clinical-grade"
- "Patient-ready"

---

## Figure checklist

- [x] Figure 1 — Project overview
- [x] Figure 2 — Source × label composition
- [x] Figure 3 — Evaluation design changes the conclusion
- [x] Figure 4 — Representation benchmark under LOSO
- [x] Figure 5 — Pooled vs per-HLA stratified feature effects
- [x] Figure 6 — Shortcut tests
- [ ] Figure 7 — Interpretability triangulation **DEFERRED**
- [x] Figure 8 — Candidate prioritization under strict model
- [x] Extended Data 1 — TCR / VDJdb annotation (provided as table)
- [x] Extended Data 2 — pMHC / TCR-pMHC structural visualization (live HTML viewer)
- [x] Extended Data 3 — Self-similarity safety flags (in candidate table)
- [ ] Extended Data 4 — Architecture ablation **DEFERRED**
- [ ] Extended Data 5 — PLM scaling comparison **PARTIAL** (35M only; 150M/650M deferred)

## Table checklist

- [x] Table 1 — Benchmark composition (`manuscript/tables/benchmark_composition.tsv`)
- [x] Table 2 — Candidate prioritization (`manuscript/tables/candidate_prioritization.tsv`)
- [x] Table S1 — Source × label matrix (`source_label_matrix.tsv`)
- [x] Table S2 — HLA × label matrix (`hla_label_matrix.tsv`)
- [x] Table S3 — Cross-source conflicts (`conflict_labels.tsv`)
- [x] Table S5 — Representation benchmark all configs (`representation_benchmark_summary.tsv`)
- [x] Table S6 — Feature confounding (`feature_confounding.tsv`)
- [x] Table S7 — Shortcut tests (`shortcut_tests.tsv`)
- [x] Table S8 — Claim audit (`claim_audit_table.tsv`)
- [x] Table S9 — Artifact inventory (`artifact_inventory.tsv`)

## Reproducibility checklist

- [x] All scripts in `scripts/` directory
- [x] Sha256 manifest of key cached files
- [x] Run log with environment info
- [x] YAML config file
- [x] Run instructions in `manuscript/reproducibility_audit.md`

## Data availability checklist

- [x] Master corpus on public hub
- [x] All result TSVs / JSONs on public hub
- [x] Six RCSB pMHC PDBs available for download
- [x] ESMFold structures cached
- [x] ESM2-35M embeddings cached
- [ ] Code repository (e.g. GitHub) — pending publication

## Code availability checklist

- [x] Run scripts present
- [x] Configuration file
- [x] Reproducibility audit
- [ ] Public GitHub mirror — TODO
- [ ] LICENSE file — TODO
- [ ] CITATION.cff — TODO

---

## Immediate next-action recommendation for Seungho

**Option A — submit now to NAR Database / BIB** (estimated 2 days of polish + 1 week peer review).
- Pros: Acceptance very likely (resource paper fit). Locks in priority on the leakage-aware benchmark claim. dbPepNeo rescue is a citable contribution.
- Cons: ESM2 LOSO improvement is modest in absolute terms; reviewer may push for scaling.

**Option B — run ESM2-650M (GPU $5-10 Azure burst) + NetMHCpan baseline under strict splits + 2-layer MLP architecture ablation, then submit to Cell Reports Methods or BIB** (estimated 1 week).
- Pros: Stronger ML novelty. PLM-scaling story matters. NetMHCpan baseline closes the "but you didn't compare to standard tool" gap.
- Cons: Marginal extra time. GPU burst cost.

**Option C — pilot wet-lab tetramer assay on 5 KRAS G12D candidates (collaborator + ~$5k + 6-12 weeks) and submit to Cell Reports Medicine or npj Digital Medicine.**
- Pros: Highest impact. Validation of even one candidate substantially strengthens the story.
- Cons: Time + cost + collaborator dependency.

**Recommendation**: **Option B**. ESM2-650M + NetMHCpan baseline + manuscript polish in 1 week → submit to Briefings in Bioinformatics. Hold Option C as follow-up companion paper.
