<!-- Paste 전체 to Codex / 또는 § 1 + § 2 + 해당 paper § 3 block 만 paste -->
---
title: "Codex prompt — THCA project 8-paper manuscript writing"
date: 2026-05-04
author_handoff: Seungho Cook (kukshomr@gmail.com)
codex_target: "GPT-5 / Claude / OpenAI Codex / Gemini — agnostic"
purpose: "8 papers 전체 context + per-paper writing assignment + SCOPE rules + verification"
---

# Codex Prompt — THCA Project 8-Paper Writing Assignment

> **Codex 사용법**
> (a) **Whole-portfolio context 가 필요할 때** → § 0 + § 1 + § 2 + § 3 (모든 paper) + § 4 + § 5 paste
> (b) **단일 paper 만 쓸 때** → § 0 + § 1 + § 2 + § 3 의 해당 paper block + § 4 + § 5 paste
> (c) Codex 가 file system access 가능 → § 6 의 file path 그대로 read 권장
> (d) Codex 가 file system access 불가 → § 7 의 inline fact bundle 사용

---

## § 0 — Codex, read this prompt completely BEFORE writing anything

You are writing manuscripts for an 8-paper portfolio in thyroid cancer / autoimmune-thyroid genomics. The author (Seungho Cook) has built up substantial scaffolding (briefs, dashboards, audit reports, source data) which YOU MUST reference. Do NOT invent facts. Do NOT cross-contaminate paper SCOPEs. Verify against source files before output.

**Critical SCOPE rules** (violation = manuscript rejection):
- Paper 2 ≠ Paper 4 (Hashimoto vs Graves)
- Paper 1 ≠ Paper 3 (driver-mutation cancer vs ICI vulnerability)
- Paper 2 ≠ Paper 1 (autoimmune layer vs cancer mechanism)
- Citation: Chu X et al. 2018 (NOT "Chen 2018")
- Author 본인 voice preserved where author has marked sections (Hook / Aim / Disc 3.1 / Limit / Cover Para 1 / Reviewer Q9) — Codex writes drafts, author edits

---

## § 1 — Project context (THCA portfolio overview)

### 1.1 Project meta

- **Domain**: Thyroid cancer (papillary thyroid carcinoma, PTC) + autoimmune thyroid disease genomics
- **Author**: Seungho Cook (1st across all 8 papers); Senior author: Yu Hyeong-won
- **Phase**: Marathon mode 2026-05-04 to 2026-06-13 (6-week manuscript writing window)
- **Hard deadlines**: NeurIPS 2026 = May 4 abstract / May 6 full paper (AOE)
- **Cell Press / npj / Bioinf / Genome Med = rolling submission**

### 1.2 8 papers at a glance

| ID | Paper | Target venue | Status | Codex writing scope |
|---|---|---|---|---|
| **0** | npj 8-gene RAI biomarker | npj Precision Oncology | ship-ready, outreach pending | minor edits only (manuscript_v6 done) |
| **1** | Cancer paper v8 (DM1 fusion + epigenetic dark matter) | Cell Rep Med 1순위 → JCI Insight + Nat Commun reach | drafting W1-W6 | full manuscript writing (00-09 sections) |
| **2** | Hashimoto-overlap PTC (HT-only) | Cell Rep Med / JCI Insight | brief done, manuscript TBD | full manuscript after Yu meeting |
| **3** | ICI vulnerability dark thyroid cancer | TBD | FROZEN bundle (chmod 444) | NO TOUCH — reserve only |
| **4** | Korean GD HLA / Pan-Asian (backlog) | TBD | gated 4/4 | NO WRITING until Bundang n>50 + Yu approval |
| **B** | Bioinformatics OUP — DIAL framework methods | Bioinformatics OUP | cover letter done, manuscript TBD | full manuscript writing |
| **G** | Genome Medicine — DM1/DM2 multi-omics axis | Genome Medicine | cover letter done, manuscript TBD | full manuscript writing |
| **N** | NeurIPS DIAL theorem | NeurIPS 2026 main | 12 pages compiled, GO verdict | revisions only (full paper compiled) |

### 1.3 Cross-paper boundary (CRITICAL — Codex must respect)

| Paper | Allowed terms | Forbidden terms (= immediate STOP + rewrite) |
|---|---|---|
| **Paper 0** | 8-gene panel, BRAF V600E, RAI, TCGA, MSK, K2 pilot | DM1/DM2 mechanism deep dive (Paper 1), Hashimoto-overlap (Paper 2), DIAL (Paper B/N) |
| **Paper 1** | DM1, DM2, fusion (RET/NTRK/ALK/BRAF), epigenetic, methylation, dark matter, 8-gene, TCGA, MSK | Graves'/GD/TSAb, ICI vulnerability (Paper 3), DIAL methods deep dive (Paper B/N) |
| **Paper 2** | Hashimoto thyroiditis (HT), HT-overlap PTC, HLA-II, TLS, BCR clonality, AICDA, mediation, Korean baseline | **Graves'/GD/TSAb/TSI/hyperthyroidism/exophthalmos/TED, Pan-Asian autoimmune-thyroid, Chu 2018 GD comparator (Paper 4 reserve), autoimmune-PTC** |
| **Paper 3** | ICI vulnerability, dark thyroid cancer | (FROZEN — no writing) |
| **Paper 4** | Korean Graves' disease, Pan-Asian autoimmune-thyroid, Chu 2018 forest, cookHLA SNP-based | (gated — no writing until Paper 1 bioRxiv + Paper 2 A/B/C + Bundang n>50 + Yu approval) |
| **Paper B** | DIAL framework, leak detection, batch correction, ComBat-seq, gene-ID recovery, GSEA, reproducibility | Clinical translation deep dive (Paper 0/1), DIAL theorem proofs (Paper N) |
| **Paper G** | DM1/DM2 axis, multi-omics integration, PRISM drug screen, pan-cancer transfer, immune evasion (PD-L1/IDO1/HLA/CTLA-4), single-cell | DIAL methodology deep dive (Paper B), clinical biomarker primary framing (Paper 0) |
| **Paper N** | DIAL theorem, subspace-aligned conditional shift, batch correction theory, 4-domain synthetic | Biomedical translation (Paper B/0), thyroid cancer biology (Paper 1/2) |

### 1.4 Citation rule (every paper)

- ✅ **Chu X et al. 2018** *J Med Genet* 55(10):685–692, doi:10.1136/jmedgenet-2017-105146 (PMC 6161647)
- ❌ NEVER write "Chen 2018" — this is a corrected misattribution

### 1.5 Author voice preferences

- 본인 = methods-driven + transparent disclosure + multi-cohort validation 강조
- Avoid: speculative mechanism claims, single-cohort overclaim, predatory venue (Frontiers / MDPI 회피)
- Prefer: honest limitations 8-row table, 5-Pillar evidence framework, paper-defining figure 표시 + WHY 박스
- Prior tools 강조 (cookHLA Nat Commun first-author, agentic-research framework)

---

## § 2 — Master file index (Codex reads BEFORE writing)

### 2.1 Top-level dashboards (read first)

```
project/submission/papers_overview.html        # 6 papers master dashboard
project/three_papers_index.html                # 3-trajectory landing
project/manuscript_p2_brief/SITUATION_OVERVIEW.md  # Paper 2 current state
project/manuscript_v8/dashboard.html            # Paper 1 manuscript navigation
```

### 2.2 Per-paper anchor docs

```
# Paper 0 (npj)
project/submission/npj/README.md
project/submission/npj/manuscript_v6.md         # current submission text
project/submission/npj/reviewer_faq.html        # 20+ Q&A
project/submission/npj/easy_explainer.html      # plain-Korean explanation
project/submission/npj/SUBMISSION_CHECKLIST.md

# Paper 1 (Cancer v8)
project/manuscript_v8/README.md
project/manuscript_v8/00_title_candidates.md    # 3 후보 + Hybrid Cand 1 ★
project/manuscript_v8/01_abstract.md            # Cell Press 153 words v2
project/manuscript_v8/02_outline.md             # Mermaid + 5/4-6/13 schedule
project/manuscript_v8/03_introduction.md
project/manuscript_v8/04_intro_1_1_hook.md      # author voice anchor
project/manuscript_v8/04_results.md
project/manuscript_v8/05_figure_captions.md     # Figure 1-8 + Suppl S1-S9
project/manuscript_v8/06_discussion.md          # 3.1-3.4 + Limitations
project/manuscript_v8/07_star_methods.md
project/manuscript_v8/08_cover_letter.md
project/manuscript_v8/09_reviewer_qa.md         # 12 Q&A pre-empt
project/manuscript_v8/13_supplementary_tables.md
project/manuscript_v8/03_intro_references.bib   # 17 cites
project/reports/yu_review_dashboard/index.html  # 5-pillar visual

# Paper 2 (Hashimoto-overlap PTC)
project/manuscript_p2_brief/paper2_brief.html       # ★ canonical story brief (38 pages)
project/manuscript_p2_brief/p2_advisor_discussion.html  # extended kitchen-sink
project/manuscript_p2_brief/SITUATION_OVERVIEW.md   # current state dashboard
project/manuscript_p2_brief/YU_MEETING_REVIEW_PACKET_2026_05_04.md
project/manuscript_p2_brief/YU_MEETING_TALK_CARD_2026_05_04.md
project/manuscript_p2_brief/COHORT_ACCESS_GUIDE.md
project/manuscript_p2_brief/DATA_SOURCES_INDEX.md   # 26/26 paths verified
project/results/p2_pillar1_forest_v2/               # ★ ACTIVE (HT context)
project/results/p3_gse286332/                       # Pillar II
project/results/d4p2_tcga_hashimoto_signature/      # Pillar III
project/results/d5p6_bcr_repertoire/                # Pillar IV
project/results/d6p7_dm1_subcluster/                # Pillar V
project/results/d3p5_pdm1_gradient/                 # § 10 Mediation

# Paper 3 (ICI — FROZEN)
project/reports/paper3_ici/PAPER3_ICI_TRACK_A_BUNDLE.md  # chmod 444, NO TOUCH

# Paper 4 (Korean GD — backlog gated)
# NO files to write — gated until 4/4 conditions

# Paper B (Bioinformatics OUP — DIAL framework)
project/submission/bioinformatics/cover_letter.md
project/submission/bioinformatics/README.md         # ★ NEW context

# Paper G (Genome Medicine — DM1/DM2 multi-omics)
project/submission/genome_med/cover_letter.md
project/submission/genome_med/README.md             # ★ NEW context

# Paper N (NeurIPS — DIAL theorem)
project/submission/v15_neurips/v15_neurips_SUBMIT_FINAL.pdf  # 12 pages compiled
project/submission/v15_neurips/v15_NEURIPS_DASHBOARD.html
project/submission/v15_neurips/audit_report.md
project/submission/v15_neurips/claims_audit.md      # 9/9 verified
project/submission/v15_neurips/rebuttal_prep_v1.md
project/submission/v15_neurips/camera_ready_prep_plan.md
```

### 2.3 Cross-paper memory (auto-context)

```
~/.claude/projects/-home-seungho-personal-THCA-data-analysis/memory/MEMORY.md  # 24 entries index
- v18_paper2_HT_isolated.md      # Paper 2 SCOPE canonical
- v17_paper2_pillar1_forest_strong.md  # SUPERSEDED (Chu 2018 framing)
- v17_graves_pivot.md            # SUPERSEDED (2026-04-29 v0)
- v17_npj_ship_status.md         # Paper 0 status
- v17_marathon_mode_post_pillar1.md  # 5/4-6/13 marathon
- v17_sprint_vs_marathon_violation.md  # voice-protected sections rule
- v17_K2_vs_bundang_distinction.md  # K2 ≠ Bundang
- v17_korean_k2_calibration.md   # 8-gene mini-index TPM 10-100× inflation
- paper_numbering_2026_05_04.md  # canonical 4-paper portfolio
- v19_paper3_ici_track_a.md      # Paper 3 FROZEN
- v19_paper4_GD_backlog.md       # Paper 4 gating
```

---

## § 3 — Per-paper writing blocks

### 3.0 Paper 0 — npj 8-gene RAI biomarker (MINOR EDITS ONLY)

```yaml
status: ship-ready 2026-04-27, outreach pending
target: npj Precision Oncology
title: "An 8-gene RAI-responsiveness biomarker outperforms BRAF V600E status in papillary thyroid carcinoma"
codex_scope: minor edits only (manuscript_v6 done; do NOT rewrite)
cohort: TCGA n=504 + MSK n=117 + K2 PRJEB11591 Korean pilot n=9
headline: ΔAUC = +0.132 vs BRAF V600E alone, 5-fold CV AUC 0.954
survival: lifelines TERT⁺ Cox HR=4.33 (univariate joint 4-group), HR=0.95 multivariate
read_first: project/submission/npj/manuscript_v6.md, reviewer_faq.html, SUBMISSION_CHECKLIST.md
forbidden: DM1/DM2 fusion+epigenetic mechanism deep dive (Paper 1), Hashimoto layer (Paper 2)
output: only minor revisions if author asks (e.g., word reduction, citation update)
```

### 3.1 Paper 1 — Cancer paper v8 (Dark Matter)

```yaml
status: drafting W1-W6 (5/4-6/13)
target: Cell Reports Medicine 1순위 → JCI Insight + Nat Commun reach → npj Prec Onc fallback
title: "An 8-gene panel reveals fusion-driven, epigenetically silenced dark matter in BRAF/RAS-negative thyroid cancer"
codex_scope: full manuscript writing (00-09 sections)
sections_to_write:
  - 00 title (3 candidates exist; Hybrid Cand 1 ★)
  - 01 abstract (Cell Press structured 153 words target; v2 exists, can polish)
  - 03 introduction 1.1-1.4 (715 w; Hook 1.1 = author voice anchor in 04_intro_1_1_hook.md)
  - 04 results 5 sub-results (3,250 w target; 대안 A ordering)
  - 05 figure captions Fig 1-8 + Suppl S1-S9
  - 06 discussion 3.1-3.4 + Limitations (1,380 w; 3.1 = author voice; Limitations = author voice)
  - 07 STAR Methods 5 sub-sections
  - 08 cover letter (Para 1 = author voice; rest = Codex)
  - 09 reviewer Q&A 12 questions (Q9 = author voice; Q1-Q8, Q10-Q12 = Codex)
voice_protected_for_author:
  - 04_intro_1_1_hook.md (Hook only)
  - 06_discussion 3.1 paragraph
  - 06_discussion Limitations
  - 08_cover_letter Paragraph 1
  - 09_reviewer_qa Q9
voice_codex_writes:
  - All other sections
cohort: TCGA-THCA n=504 + MSK-IMPACT n=117 + Korean cohorts n=874 + cBioPortal SV n=542 + HM450 n=503 + sc Pu 2021 (n=7) + Lu 2023 (n=23)
3_layer_DM1_mechanism:
  L1_genetic: 76.8% tyrosine kinase fusion+ (RET/NTRK/ALK/BRAF), OR 7.41 vs DM2
  L2_heterogeneity: fusion+ young/well-diff vs fusion- older/immune-hot/Hashimoto-overlap
  L3_epigenetic: TPO promoter d=2.30, DIO1 d=1.24, TSHR d=1.20 (fusion-independent)
key_numbers:
  - DM1 captures 81.8% of TCGA RET-fusion-positive
  - Pooled OS HR 2.53 [1.31, 4.89] (TCGA + MSK meta, I²=0%)
  - 8-gene panel ΔAUC=+0.132 vs BRAF V600E
  - mean 8-gene methylation β 0.385 (DM1) vs 0.253 (DM2)
read_first: project/manuscript_v8/README.md, 02_outline.md, 01_abstract.md, results/audit_2026_04_30/FINAL_COMPREHENSIVE_SUMMARY_v9.md
forbidden: Hashimoto-overlap PTC mechanism deep dive (Paper 2), Graves'/GD (Paper 4), DIAL methodology (Paper B/N), ICI vulnerability (Paper 3)
output_format: Cell Press structured, ~5,500 words main + supplementary
```

### 3.2 Paper 2 — Hashimoto-overlap PTC (HT-ONLY)

```yaml
status: brief done, manuscript TBD (after Yu meeting)
target: Cell Reports Medicine / JCI Insight
title_proposal: "A Hashimoto-overlap molecular layer of papillary thyroid carcinoma — Korean PTC HLA susceptibility, in-cohort mechanism, and HLA-II–mediated dedifferentiation"
codex_scope: full manuscript writing AFTER Yu meeting (not before)
SCOPE_RULES_HARD:
  - Hashimoto thyroiditis (HT) overlap mechanism ONLY
  - NEVER mention: Graves, GD, TSAb, TSI, hyperthyroidism, thyrotoxicosis, exophthalmos, Graves' ophthalmopathy, TED, thyroid eye disease
  - NEVER use: "autoimmune-PTC", "autoimmune-thyroid axis", "Pan-Asian autoimmune-thyroid susceptibility"
  - Use instead: "Hashimoto-overlap PTC", "Hashimoto-thyroid overlap axis", "Pan-Asian thyroid HLA susceptibility (HT context)"
  - Pillar I comparator = Korean baseline AFND South Korea pool (NOT Chu 2018 GD)
  - Chu 2018 GD = ONE Discussion line cross-ref to Paper 4 only
  - HLA substrate = HT-overlap explanation only; GD-specific 해석 = Paper 4
5_pillars:
  - I: Korean PTC HLA pool n=874 vs Korean baseline (HT context). DPB1*05:01 OR 1.96 [1.60, 2.41] vs Korean baseline, p=1.1e−10
  - II: GSE286332 PTC vs PTC+HT in-cohort. 8-gene Cohen d=−1.60, HLA-II d=+3.65, IFN-γ FDR=2e−4, 10,380 DEGs
  - III: TCGA + Korean cross-cohort. DM2 Hashimoto OR=0.20, p=6.4e−10. 18-30% prevalence robust
  - IV: BCR clonal + TLS antigen-driven. TLS Cabrita d=+1.96, IGHV clonality d=+2.09, AICDA up
  - V: DM1 sub-B = NBNR cluster. 96% mutation-negative; Korean transfer 47-53%
synthesis: PTC+HT → HLA-II → P(DM1) ↓; HLA-II 140% mediation, generic immune NS
cohort: TCGA n=500 + GSE286332 (Lim 2025 Dongguk n=18) + GSE213647 (Lee 2024 n=632) + K2 PRJEB11591 n=235 typeable
read_first: project/manuscript_p2_brief/paper2_brief.html, SITUATION_OVERVIEW.md, YU_MEETING_REVIEW_PACKET_2026_05_04.md
read_data: project/results/p2_pillar1_forest_v2/PILLAR1_FOREST_V2_SUMMARY.md (★ active), discussion_paragraph_v2.md (paste-ready)
forbidden: All GD-specific terms; Chu 2018 GD as primary comparator; "autoimmune-PTC" framing
voice_protected_for_author: Hook, Aim, Disc 3.1, Limitations, Cover Para 1, Q9 (same as Paper 1)
output_format: Cell Press structured, ~5,500 words target
```

### 3.3 Paper 3 — ICI vulnerability (FROZEN)

```yaml
status: FROZEN bundle (chmod 444 at project/reports/paper3_ici/PAPER3_ICI_TRACK_A_BUNDLE.md)
codex_scope: NO WRITING — reserve only
note: Track B BLOCKED until Paper 1 bioRxiv + Paper 2 A/B/C + explicit "Track B 시작" by author
forbidden: ANY content generation for Paper 3
```

### 3.4 Paper 4 — Korean GD HLA (backlog gated)

```yaml
status: gated until 4/4 conditions: (a) Paper 1 bioRxiv (b) Paper 2 A/B/C complete (c) Bundang Graves' n>50 (d) Yu approval
codex_scope: NO WRITING until all 4 met
note: Chu 2018 GD forest material reserved here; Paper 2 cross-references in ONE Discussion line
forbidden: ANY content generation for Paper 4 until ungated
```

### 3.B Paper B — Bioinformatics OUP (DIAL framework methods)

```yaml
status: cover letter done 2026-04-27, manuscript TBD
target: Bioinformatics (Oxford University Press)
title: "Recovered external validation and proper GSEA define DM1/DM2 as an actionable thyroid cancer transcriptional axis"
codex_scope: full manuscript writing (methods voice)
triple_deliverable:
  - DIAL framework — leak-detection safeguard for batch-correction protocols (v5.1 → v5.2 self-audit)
  - 5-cohort validation pipeline — DM1/DM2 axis recovery
  - Proper preranked GSEA — pathway-name normalisation, 100% sign-agreement with prior proxy GSEA
cohort: TCGA-THCA + GSE27155 + GSE33630 + GSE29265 + GSE76039
reproducibility:
  - random_state=42 deterministic seeds
  - sklearn 1.8.0, gseapy 1.13, scanpy 1.12.1 version-pinned
  - per-fold ComBat-seq batch correction
  - notebooks_or_scripts/v17p35_*.py
clinical_relevance: 8-gene RAI panel ΔAUC=+0.132 vs BRAF V600E (downstream application, NOT main framing)
suggested_reviewers: Aedin Culhane (U Limerick), Lior Pachter (Caltech), Casey Greene (UC Boulder)
read_first: project/submission/bioinformatics/cover_letter.md, README.md
forbidden: Clinical translation primary framing (= Paper 0 territory), DIAL theorem proofs (= Paper N territory), DM1 mechanism deep dive (= Paper 1 territory)
output_format: Bioinformatics OUP methods voice, ~6,000 words
```

### 3.G Paper G — Genome Medicine (DM1/DM2 multi-omics axis)

```yaml
status: cover letter done 2026-04-27, manuscript TBD
target: Genome Medicine
title: "A driver-orthogonal transcriptomic axis (DM1/DM2) stratifies thyroid cancer immune state and 8-gene RAI responsiveness"
codex_scope: full manuscript writing (mechanistic + multi-omics integration voice)
key_findings:
  - 513 TCGA-THCA primary tumor unsupervised → DM1/DM2 axis
  - Differentiation continuum: TPO/DIO1/FOXE1 high → DM2; low → DM1
  - Immune continuum: cold OxPhos → DM2; hot inflammatory/IFN-γ → DM1
  - Driver-orthogonal: 17/335 mutation-carrying tumor violate canonical BRAF/RAS map
  - 8-gene panel 5-fold CV AUC 0.954, ΔAUC=+0.132 vs BRAF V600E
6_layer_multi_omics:
  - bulk RNA × 5 cohorts (TCGA + GSE27155 + GSE33630 + GSE29265 + GSE76039)
  - sc-RNA 66K cells / 7 patient (GSE184362)
  - GSEA Hallmark v2024.1.Hs (51 pathways, FDR<1e−15 inflammatory in DM1)
  - immune-evasion (PD-L1, IDO1, HLA, CTLA-4) bulk + sc
  - PRISM Repurposing 19Q4 (4,517 compounds; MEK + HMGCR DM1-selectivity)
  - pan-cancer transfer LUAD/COAD/LGG/SKCM (CDKN2A, FOSL1, ETV4, DUSP6)
honest_disclosure:
  - 2/5 robust external transfer rate
  - n=5 vs n=5 PRISM constraint (FDR-grade single-drug discovery limited)
  - wet-lab validation absent
suggested_reviewers: Aleix Prat (Hospital Clínic Barcelona), Charles Perou (UNC), Yuri Nikiforov (Pittsburgh)
read_first: project/submission/genome_med/cover_letter.md, README.md
forbidden: DIAL methodology deep dive (= Paper B territory), clinical biomarker primary framing (= Paper 0 territory), Paper 1 fusion+epigenetic primary mechanism
output_format: Genome Medicine mechanistic depth voice, ~6,000 words
```

### 3.N Paper N — NeurIPS DIAL theorem (REVISIONS ONLY)

```yaml
status: 12 pages compiled, GO verdict, 2 manual steps remaining (OpenReview profile + anonymous_code.zip URL)
target: NeurIPS 2026 main track, double-blind
title: "DIAL: A post-hoc diagnostic for subspace-aligned conditional shift under linear batch correction"
codex_scope: revisions only (full paper compiled, claims 9/9 verified)
hard_deadline: May 4 abstract / May 6 full paper AOE
content_summary:
  - DIAL framework — 1-number AUC-based diagnostic for batch correction class signal leak
  - Theorem 2 stays (8/8 todos resolved)
  - 3 lemmas + 1 proposition + 1 conjecture
  - 5 experiments (synthetic 4-domain reproduction + THCA/ComBat)
  - 52 bibliography entries
fallback_venues: TMLR (rolling, no deadline) → ICLR 2027 → AISTATS 2027
read_first: project/submission/v15_neurips/v15_neurips_SUBMIT_FINAL.pdf, audit_report.md, claims_audit.md
forbidden: ANY content generation that conflicts with anonymisation (double-blind)
output: rebuttal responses (after reviewer comments arrive Aug 2026), camera-ready additions (if accepted Sep 2026)
```

---

## § 4 — Output requirements

### 4.1 Per-venue style

| Venue | Voice | Word count | Structure |
|---|---|---|---|
| Cell Reports Medicine (Paper 1, Paper 2) | Mechanistic depth + clinical relevance + multi-cohort validation | ~5,500 main + suppl | Cell Press structured Abstract / 4-section Intro / 5 sub-Results / 4-section Discussion / STAR Methods |
| Genome Medicine (Paper G) | Mechanistic + multi-omics integration | ~6,000 | Standard journal sections + extensive supplementary |
| Bioinformatics OUP (Paper B) | Methods + reproducibility | ~6,000 | Methods-paper structure: Introduction / Methods (heavy) / Results / Discussion |
| npj Precision Oncology (Paper 0) | Clinical biomarker translation | ~5,000 (existing v6) | npj structured (Background / Methods / Results / Discussion) |
| NeurIPS 2026 (Paper N) | ML theory + empirical | 8 pages main + 4 appendix | NeurIPS LaTeX class, double-blind |

### 4.2 Output naming convention (Paper 1/2/B/G drafts)

```
00_title_candidates.md
01_abstract.md
02_outline.md          # Mermaid + bullet schedule
03_introduction.md
04_results.md
05_figure_captions.md
06_discussion.md
07_methods.md          # or STAR Methods for Cell Press
08_cover_letter.md
09_reviewer_qa.md
10_full_manuscript_compiled.md
13_supplementary_tables.md
```

### 4.3 Always include in output

- Word count at top of each file
- Source citations from the manuscript's reference.bib (do NOT invent citations)
- Cross-paper boundary check at bottom of each file (forbidden terms grep result)

---

## § 5 — Self-verification checklist (Codex runs BEFORE output)

After writing each draft, verify the following before returning to author:

### 5.1 Citation correctness

- [ ] All "Chen 2018" → "Chu X et al. 2018" (J Med Genet 55:685–692)
- [ ] No invented citations (every cite traceable to bib file or established literature)
- [ ] DOI format: `10.xxxx/yyy.zzz`
- [ ] PMC IDs use format `PMC1234567`

### 5.2 Cross-paper boundary

For Paper 2 specifically:
- [ ] grep -i "graves" output = 0 (or only Discussion 한 줄 cross-ref to Paper 4)
- [ ] grep "TSAb\|TSI\|hyperthyroid" output = 0
- [ ] grep "autoimmune-PTC\|autoimmune-thyroid axis" output = 0
- [ ] grep "Pan-Asian autoimmune-thyroid" output = 0
- [ ] Pillar I comparator = Korean baseline (NOT Chu 2018 GD)

For Paper 1 specifically:
- [ ] No deep dive into Hashimoto-overlap PTC mechanism (= Paper 2)
- [ ] No DIAL methodology details (= Paper B/N)

For Paper B specifically:
- [ ] Methods voice (NOT clinical translation primary)
- [ ] No DIAL theorem proofs (= Paper N)

For Paper G specifically:
- [ ] Multi-omics integration framing
- [ ] No DIAL methodology deep dive (= Paper B)

### 5.3 Voice protection

- [ ] Sections marked "voice_protected_for_author" 는 작성 X 또는 placeholder 만 (e.g., `[AUTHOR HOOK — voice insertion point]`)
- [ ] Sections "voice_codex_writes" = full Codex draft

### 5.4 Numerical accuracy

- [ ] Cohort n values match source (TCGA n=500 or 504; MSK n=117; K2 n=235 typeable; etc.)
- [ ] Effect sizes (Cohen d, OR, HR) match source TSV/JSON
- [ ] p-values match (1.1e−10, 6.4e−10, etc.)
- [ ] Reference cohort sizes (Chu 2018: 1,468 GD / 1,490 ctrl)

### 5.5 Honest disclosure

- [ ] Limitations section addresses known weaknesses (cohort size, retrospective, single-cohort, etc.)
- [ ] Avoid overclaim language ("this proves", "definitively shows" → "supports", "consistent with")
- [ ] Brand expressions preserved (cookHLA Nat Commun first-author 강조)

---

## § 6 — File path quick reference (Codex with file system access)

```bash
# Project root
PROJ=/home/seungho/personal/THCA_data_analysis

# 8 paper anchors
PAPER0=$PROJ/project/submission/npj
PAPER1=$PROJ/project/manuscript_v8
PAPER2=$PROJ/project/manuscript_p2_brief
PAPER3=$PROJ/project/reports/paper3_ici  # FROZEN
PAPER4=  # gated, no dir yet
PAPERB=$PROJ/project/submission/bioinformatics
PAPERG=$PROJ/project/submission/genome_med
PAPERN=$PROJ/project/submission/v15_neurips

# Master indexes
$PROJ/project/submission/papers_overview.html
$PROJ/project/three_papers_index.html

# Source data dirs (Paper 1 + 2 share)
$PROJ/project/results/p2_pillar1_forest_v2/    # Paper 2 Pillar I active
$PROJ/project/results/p3_gse286332/            # Paper 2 Pillar II
$PROJ/project/results/d4p2_tcga_hashimoto_signature/  # Paper 2 Pillar III
$PROJ/project/results/d5p6_bcr_repertoire/     # Paper 2 Pillar IV
$PROJ/project/results/d6p7_dm1_subcluster/     # Paper 2 Pillar V (+ Paper 1 DM1 sub-B)
$PROJ/project/results/d3p5_pdm1_gradient/      # Paper 2 Mediation
$PROJ/project/results/audit_2026_04_30/        # Paper 1 v9 comprehensive summary

# Memory
$HOME/.claude/projects/-home-seungho-personal-THCA-data-analysis/memory/
  - MEMORY.md                          # 24-entry index
  - v18_paper2_HT_isolated.md          # Paper 2 SCOPE canonical
  - paper_numbering_2026_05_04.md      # canonical 4-paper portfolio
```

---

## § 7 — Inline fact bundle (Codex without file system access)

If Codex cannot read files, use this fact bundle directly:

### 7.1 Cohort facts

```
TCGA-THCA: n=504 primary tumor (some analyses 500-513), Pan-Cancer multi-omics
MSK-IMPACT thyroid 2016 (Landa Cell): n=117 advanced disease, mutation-focused
K2 / PRJEB11591 (Yoo 2016 SNU-GMI): n=235 typeable Korean PTC (260 manifest)
Lee 2024 / GSE213647: n=632 Korean PTC validation
GSE286332 (Lim 2025 Dongguk Univ Macrogen Illumina NovaSeq X): n=18 (9 PTC + 9 PTC+HT)
Chu 2018 (PMC 6161647) Han Chinese GD: n=1,468 GD vs n=1,490 ctrl, 4-digit SNP2HLA
Pu 2021 sc (GSE184362): n=7 PTC patients, 66K cells
Lu 2023 sc (GSE193581): n=23 samples
PRISM Repurposing 19Q4: 4,517 compounds cell-line drug screen
```

### 7.2 Paper 1 key numbers

```
DM1 fusion+ rate: 76.8% (RET/NTRK/ALK/BRAF), OR 7.41 vs DM2
Pooled OS HR: 2.53 [1.31, 4.89] (TCGA + MSK meta), I²=0%
8-gene panel 5-fold CV AUC: 0.954
ΔAUC vs BRAF V600E alone: +0.132
DM1 captures 81.8% of TCGA RET-fusion-positive
TPO promoter methylation (DM1 vs DM2): Cohen d=2.30
DIO1 promoter methylation: Cohen d=1.24
TSHR promoter methylation: Cohen d=1.20
mean 8-gene methylation β: DM1 0.385 vs DM2 0.253 (fusion-independent)
Driver_anchor cluster ARI: 0 (random; vs DM1/DM2)
TIERA67 ARI: 0.90 ≈ pan-genome top-5000 ARI 0.92
8-gene alone ARI: 0.49
BRAF mRNA Cohen d (V600E vs WT): -0.04 (NS)
```

### 7.3 Paper 2 key numbers (HT-only, v2 Korean baseline)

```
Korean PTC HLA pool n=874 = K2 235 + Lee 630 + GSE286332-PTC 9
Korean baseline = AFND South Korea pool (sample-weighted)
DPB1*05:01: PTC 53.2% vs baseline 36.7%, OR 1.96 [1.60, 2.41], p=1.1e-10
A*02:07: PTC 8.1% vs baseline 3.4%, OR 2.54 [1.90, 3.40], p=3.1e-10
B*46:01: PTC 10.3% vs baseline 4.7%, OR 2.33 [1.80, 3.01], p=1.2e-10
DRB1*07:01: PTC 11.4% vs baseline 5.0%, OR 2.46 [1.26, 4.79], p=0.0084 (direction discord caveat)
C*01:02 / DQB1*02:01: AFND South Korea unavailable (Lee 2014 Korean ref TBD)

GSE286332 PTC vs PTC+HT (n=18):
- 10,380 DEGs (PyDESeq2 padj<0.05)
- 8-gene RAI Cohen d = -1.602 (MW p=4.1e-4)
- HLA-II module Cohen d = +3.65 (MW p=4.1e-4) ★ saturated
- HLA-I module Cohen d = +2.34
- Hallmark Allograft Rejection NES=+2.12 (FDR=0)
- KEGG Type I Diabetes NES=+1.92 (FDR=0)
- Hallmark IFN-γ NES=+1.80 (FDR=2e-4)

TCGA cross-cohort (n=500):
- DM1 Hashimoto-like 10.7% vs DM2 37.5%, Fisher OR=0.20, p=6.4e-10
- Stromal+immune residualized OR=0.29, p=8e-9 (specific axis, not generic immune)
- Hashimoto-like prevalence 18-30% across 4 thresholds

Korean GSE213647 (n=632):
- Hashimoto-like prevalence 22.8% (GMM) / 28.2% (Otsu)

BCR/TLS (Pillar IV):
- TLS Cabrita 2020 12-gene Cohen d = +1.96
- IGHV clonality Cohen d = +2.09 (p=1.1e-3)
- IGKV clonality d = +1.41 (p=0.034)
- AICDA up

§ 10 Mediation (PTC+HT → mediator → P_DM1):
- HLA-II 140% mediated (p=0.023, 95% CI [-0.31, -0.03])
- 8-gene RAI 63% partial (p=0.002)
- HLA-I 87% (p=0.042)
- Generic immune 88% NS (p=0.120) ← confounder, not mediator

DM1 sub-B (n=56): 96% mutation-negative (NBNR), Korean transfer 47-53%
```

### 7.4 Paper B / G / N key numbers

```
Paper B (Bioinformatics):
- DIAL framework: 1-number AUC-based leak diagnostic
- v5.1 → v5.2 self-audit: ΔDIAL=-0.49 caught
- 5 cohort validation: 100% sign-agreement with prior proxy GSEA on overlapping pathways
- Reproducibility: random_state=42, sklearn 1.8.0, gseapy 1.13, scanpy 1.12.1, per-fold ComBat-seq
- Clinical relevance: 8-gene ΔAUC=+0.132 (downstream)

Paper G (Genome Medicine):
- 513 TCGA-THCA primary tumor unsupervised → DM1/DM2 axis
- 17 of 335 mutation-carrying tumor violate canonical BRAF/RAS map
- 8-gene panel 5-fold CV AUC 0.954, ΔAUC=+0.132 vs BRAF V600E
- GSEA Hallmark v2024.1.Hs: 51 pathways, FDR<1e-15 inflammatory in DM1
- PRISM 19Q4: 4,517 compounds, MEK + HMGCR DM1-selectivity (mechanism-class)
- Pan-cancer transfer (LUAD/COAD/LGG/SKCM): conserved markers CDKN2A, FOSL1, ETV4, DUSP6
- 2/5 robust external transfer rate
- n=5 vs n=5 PRISM constraint

Paper N (NeurIPS):
- DIAL theorem (Theorem 2 stays — 8/8 todos resolved)
- 3 lemmas + 1 proposition + 1 conjecture
- 5 experiments: synthetic 4-domain + THCA/ComBat
- R²=0.9405 (Theorem 2 fit)
- ROC-AUC=0.78 (broader stress test)
- naive ComBat AUC=0.33, DANN-lite AUC=0.39 (baselines)
- 4/4 domains flip in cross-domain test
- 12 pages (8 main + 4 appendix), 52 bibliography entries
```

### 7.5 Universal forbidden terms (per Paper 2)

```
Graves, GD, TSAb, TSI, thyroid-stimulating antibody
hyperthyroidism, thyrotoxicosis, exophthalmos
Graves' ophthalmopathy, TED, thyroid eye disease
Bundang Graves' (cohort direct mention as Paper 2 data)
autoimmune-PTC, autoimmune PTC
autoimmune-thyroid axis
Pan-Asian autoimmune-thyroid susceptibility
"Chen 2018" (always Chu 2018)
```

---

## § 8 — Codex output handoff to author

After writing, return to Seungho Cook (kukshomr@gmail.com) with:

1. **Output files** in conventional naming (00_title.md ~ 09_reviewer_qa.md)
2. **Self-verification report** — what was verified per § 5 checklist
3. **Voice insertion points** — list of `[AUTHOR ... VOICE INSERTION POINT]` placeholders
4. **Open questions** — any ambiguity that needs author judgment
5. **Citation list** — all cites used + their source bib entries
6. **Cross-paper boundary check log** — grep results for forbidden terms

---

## § 9 — Tool calls Codex should make (if agentic)

```
1. read_file(<§ 2 anchor docs>)         # Establish context
2. read_file(<§ 3 per-paper read_first>) # Paper-specific deep context
3. write_file(00_title_candidates.md)
4. write_file(01_abstract.md)
... (per § 4.2 naming)
5. grep_check(forbidden_terms)           # § 5.2 self-verification
6. report_to_author(deliverables + verification_log)
```

---

*Generated 2026-05-04. Author: Seungho Cook. Use this prompt with any LLM (Codex / GPT / Claude / Gemini). For single-paper writing, paste § 0 + § 1 + § 2 + § 3 (해당 paper) + § 4 + § 5 + § 7 (해당 paper) only.*

*Last fact-check date: 2026-05-04. Verify against source files for any analysis-dated facts before publishing.*
