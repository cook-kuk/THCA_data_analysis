---
title: "Paper 2 HLA β plan — Advisor packet (Yu professor meeting, 1 page)"
date: 2026-05-04
audience: Yu Hyeong-won (advisor) — also web-served for co-author / external advisor read
parent_plan: project/reports/2026_05_04_hla_strengthening_beta_plan.md
status: SPEC/PLAN ONLY — no execution performed during packet preparation
scope_lock: Paper 2 (Hashimoto-overlap PTC) HLA. Paper 3 + Paper 4 + voice-protected prose hard-blocked.
web_hub: papers_hub_2026_05_04 (40.82.129.113 local)
last_updated: 2026-05-04
---

# Paper 2 HLA β plan — Advisor packet

> **Status badges**: 🟢 v2 framing locked · 🟡 6-allele forest awaiting Q1+Q2+Q4+Q5 · 🔴 Paper 4 GD work BLOCKED · 🔒 Paper 3 ICI Track A FROZEN

## Quick TOC

1. [One-line Paper 2 HLA scope](#1-one-line-paper-2-hla-scope)
2. [Why v2 matters](#2-why-v2-matters)
3. [Current Pillar I v2 status](#3-current-pillar-i-v2-status)
4. [Advisor decisions Q1–Q7](#4-advisor-decisions-q1q7)
5. [Recommended defaults](#5-recommended-defaults-claudes-pre-fill-advisor-confirms-or-overrides)
6. [Will be executed only after approval](#6-will-be-executed-only-after-approval)
7. [What remains BLOCKED](#7-what-remains-blocked-non-negotiable)
8. [Related documents (web hub)](#8-related-documents-web-hub)
9. [Decision timeline](#9-decision-timeline)
10. [Glossary (for web readers)](#10-glossary-for-web-readers)
11. [How to act on this packet](#11-how-to-act-on-this-packet)

---

## 1. One-line Paper 2 HLA scope

Paper 2 HLA work = **Korean Hashimoto-overlap PTC HLA susceptibility (Pillar I) + HLA-II–mediated dedifferentiation mechanism (Pillars II/III)**, with all GD / Graves' / Bundang work moved to Paper 4 backlog.

---

## 2. Why v2 matters

- **v1 conflated HT and GD**: Korean PTC pool (HT-substrate, ~20-30% Hashimoto comorbidity) compared directly to Chu 2018 Han Chinese **GD** cohort. Two distinct disease entities → cross-disease mismatch → reviewer red flag.
- **v2 fixes by within-population framing**: Korean PTC pool **vs Korean general population baseline** (AFND South Korea pool ± Lee 2014 KOTRY donor reference).
- Aligns Paper 2 substrate (HT) with comparator substrate (Korean general population, which includes general thyroid health states).
- GD comparison moves to Paper 4 (backlog, gated).

---

## 3. Current Pillar I v2 status

| Item | Status |
|---|---|
| v1 forest meta (Korean PTC vs Chu GD) | ❌ DEPRECATED 5/4 |
| v2 4-allele forest (Korean PTC vs Korean baseline) | ✅ Already executed (`p2_pillar1_forest_v2/`) |
| Alleles in v2: DPB1\*05:01, A\*02:07, B\*46:01, DRB1\*07:01 | ✅ done |
| 2 additional alleles (C\*01:02, DQB1\*02:01) | ⏸ awaiting Lee 2014 fill |
| Cohort | Korean PTC pool n=874 (K2 235 + Lee 2024 630 + GSE286332-PTC 9) |
| Terminology cleanup ("autoimmune-PTC" → "Hashimoto-overlap PTC") | ✅ Paper 1 manuscript core (00-09) clean |
| Paper 4 GD framing | ✅ Moved to backlog (`v19_paper4_GD_backlog`, gated 4/4) |
| Paper 3 ICI Track A | 🔒 FROZEN (chmod 444), separate trajectory |

---

## 4. Advisor decisions Q1–Q7

| # | Question | Decision needed |
|---|---|---|
| **Q1** | Approve Pillar I v2 framing (Korean PTC vs Korean baseline) as Paper 2 Pillar I? | YES / NO / modify |
| **Q2** | Approve Lee 2014 fill for C\*01:02 + DQB1\*02:01 (extends 4-allele → 6-allele forest)? | YES / NO / defer |
| **Q3** | If 6-allele forest generated, primary in figures/Methods OR supplementary? | primary / supplementary / both panels |
| **Q4** | AFND South Korea baseline pool composition: AFND only / AFND + KOTRY (Lee 2014) / AFND + KOTRY + KBP / Lee 2014 only / other? | pick one |
| **Q5** | Harbin Korean fallback (if AFND has these entries): primary / sensitivity / exclude? | pick one |
| **Q6** | Tier 2 priority — B5 (antigen presentation gene score) vs B3 (sc HLA-II thyrocytes) vs B4 (BCR/TLS × HLA-II) vs C4 (TCGA HLA QC)? | rank order |
| **Q7** | Execute Pillar I v2 6-allele forest now (after Q1–Q5 sign-off) OR defer until after Paper 1 bioRxiv? | now / after Paper 1 |

---

## 5. Recommended defaults (Claude's pre-fill, advisor confirms or overrides)

| # | Recommended | Rationale |
|---|---|---|
| Q1 | ✅ **Approve v2** | Within-population framing fixes v1 cross-disease conflation; reviewer-defensible |
| Q2 | ✅ **YES Lee 2014 fill** | 6-allele forest > 4-allele coverage; Lee 2014 is established Korean HLA reference |
| Q3 | **6-allele primary IF Lee 2014 extraction clean; otherwise 4-allele primary + 6-allele supplementary** | Quality-conditional default |
| Q4 | **AFND + KOTRY (Lee 2014)** (G2-2) | Two-source pool: adequate n, KOTRY = recognized Korean reference, heterogeneity disclosure manageable |
| Q5 | **Harbin fallback in sensitivity panel only, not primary** (assuming Lee 2014 covers primary alleles) | Cleaner Korean baseline; Harbin = Korean-Manchurian, ancestry note |
| Q6 | **B5 first** (antigen presentation gene score) | Builds on existing v17_hla data; lowest new-execution cost; Paper 2 mechanism core (HLA-II axis activation) |
| Q7 | **After Paper 1 Hook/Aim user-voice cleanup; before Paper 1 bioRxiv submission** | Paper 1 Aim must remain HT-isolated; v2 forest result useful as Paper 2 Pillar I lock before Paper 1 ships |

---

## 6. Will be executed ONLY after approval

| Action | Trigger |
|---|---|
| A2 Lee 2014 baseline lookup | Q2 YES + user `run A2 lookup only` |
| A3 G2 baseline pool decision lock | Q4 sign-off |
| A1 6-allele forest pipeline (script + execution + outputs) | Q1 + Q2 + Q4 + Q5 + user `prepare A1 forest script only` or `run full A1 forest` |
| Tier 2 spec → execution (B5 / B3 / B4 / C4) | Q6 prioritization + user `execute B5` (etc.) |
| Figure 1 re-render | post-A1 execution + G5 sign-off |
| Manuscript paragraph fill | post-A1 execution + 본인 voice (always voice-protected) |

---

## 7. What remains BLOCKED (non-negotiable)

| Item | Reason |
|---|---|
| Paper 4 GD analysis (Chu 2018 reactivation, Bundang Graves' work) | Paper 4 backlog, 4/4 gating not met (Paper 1 bioRxiv + Paper 2 A/B/C + Bundang n>50 + Yu 합의) |
| Paper 3 ICI Track A bundle | chmod 444 FROZEN |
| Paper 3 ICI Track B (HLA LOH, neoantigen, DIAL audit, pan-cancer ICI) | BLOCKED until Paper 1 bioRxiv + Paper 2 A/B/C + explicit Track B 시작 |
| cookHLA SNP integration | Requires SNP data acquisition (currently absent for K2/Lee); large new data forbidden in β plan |
| Pathology image multimodal layer (HLA-II IHC, image-DM1, TCGA WSI) | Forbidden in β plan; Paper 1 strategy "H&E-DM1 dropped"; Paper 2 image scope undefined |
| Voice-protected prose (Paper 2 Hook / Aim / mechanism story / Limitations / Cover letter Para 1) | Always 본인 keyboard, 본인 voice |
| Manuscript edit beyond terminology cleanup | Marathon mode: spec/scaffolding only |
| New large data download | Marathon constraint |

---

---

## 8. Related documents (web hub)

For web readers wanting the full chain of evidence behind each section:

| Document | Purpose | Path |
|---|---|---|
| **β plan (full spec)** | Tier 1 + Tier 2 detailed spec for all options | `project/reports/2026_05_04_hla_strengthening_beta_plan.md` |
| HLA strengthening situation audit | Cross-paper situational + 24 options + scoring | `project/reports/2026_05_04_HLA_strengthening_situation.md` |
| HT vs GD critical audit | Conflation 1/2/3 catch — why v2 was needed | `project/manuscript_v8/_audit_HT_vs_GD_2026_05_02.md` |
| Pillar I v2 forest spec | Original spec defining 8 decision gates | `project/manuscript_v8/PAPER2_PILLAR1_V2_FOREST_META_SPEC_2026_05_04.md` |
| Lee 2014 fill decision lock (G3) | YES locked 5/4 | `project/manuscript_v8/_LEE_2014_FILL_DECISION_LOCK_2026_05_04.md` |
| Terminology correction plan | "autoimmune-PTC" → "Hashimoto-overlap PTC" sed plan | `project/manuscript_v8/PAPER2_TERMINOLOGY_CORRECTION_PLAN_2026_05_04.md` |
| Terminology scan results | 523 hits across 55 files (90 actionable initially → 85 after user cleanup) | `project/manuscript_v8/_terminology_scan_hits.tsv` |
| Existing v2 4-allele forest | Already executed | `project/results/p2_pillar1_forest_v2/PILLAR1_FOREST_V2_SUMMARY.md` |
| Paper 2 master view | Full session origin + scope lock | `project/manuscript_v8/_PAPER2_MASTER_VIEW.md` |
| Introduction v1 strategy audit | Section-by-section verdict (KEEP/MODIFY/MAJOR REFRAME/REWRITE) | `project/reports/2026_05_04_intro_v1_strategy_audit.md` |

**Memory entries (cross-reference)**:
- `v18_paper2_HT_isolated.md` — Paper 2 HT-isolation 5/4 advisor decision
- `v17_paper2_pillar1_forest_strong.md` — v1 history (SUPERSEDED)
- `v19_paper4_GD_backlog.md` — Paper 4 GD backlog gating
- `v19_paper3_ici_track_a.md` — Paper 3 ICI Track A FROZEN
- `paper_numbering_2026_05_04.md` — canonical Paper 1-4 numbering
- `v17_marathon_mode_post_pillar1.md` — marathon mode rules
- `v17_sprint_vs_marathon_violation.md` — voice-protected sections + sprint trigger

---

## 9. Decision timeline

| Date | Event |
|---|---|
| 2026-04-29 | Initial 6-audit-session sweep, 5-Pillar paper structure emerges |
| 2026-04-30 | Marathon mode locked (`v17_marathon_mode_post_pillar1`) — Pillar 1 STRONG = analysis end, 6-week manuscript writing |
| 2026-05-02 | ★ User HT vs GD critical catch → audit doc → conflation 1/2/3 identified |
| 2026-05-03 | Pillar I v1 (Korean PTC vs Chu GD direct comparison) recorded in memory; later DEPRECATED |
| 2026-05-04 (Yu meeting) | Paper renumbering: Paper 3 = ICI vulnerability NEW; Paper 4 = GD backlog; Paper 2 HT-isolated activated |
| 2026-05-04 | Pillar I v2 4-allele forest already executed in `p2_pillar1_forest_v2/` |
| 2026-05-04 | β plan written: Tier 1 + Tier 2 spec |
| 2026-05-04 | G3 Lee 2014 fill = YES locked |
| 2026-05-04 | This advisor packet prepared |
| **PENDING** | Yu meeting — Q1-Q7 decisions |
| PENDING | A1 6-allele forest execution (post Q1+Q2+Q4+Q5) |
| PENDING | Tier 2 spec → execution (per Q6 prioritization) |
| PENDING | Paper 1 bioRxiv submission (W6 marathon) |
| BLOCKED | Paper 4 GD work (4/4 gating) |
| BLOCKED | Paper 3 ICI Track B |

---

## 10. Glossary (for web readers)

| Term | Definition |
|---|---|
| **HT** | Hashimoto's thyroiditis — destructive lymphocytic infiltration of thyroid → hypothyroidism. Paper 2 substrate. |
| **GD** | Graves' disease — TSH receptor stimulating antibody → hyperthyroidism. **Paper 4 backlog**, NOT Paper 2. |
| **PTC** | Papillary thyroid carcinoma. |
| **PTC+HT** | PTC with concurrent Hashimoto's thyroiditis (~20-30% of Korean PTC). |
| **Pillar I** | Paper 2 first pillar = HLA susceptibility analysis. Other pillars (II/III) = mechanism (HLA-II–mediated dedifferentiation, mediation analysis). |
| **AFND** | Allele Frequency Net Database — public reference for HLA frequencies by population. |
| **KOTRY** | Korean Organ Transplantation Registry — Korean donor cohort, source for Lee 2014 baseline. |
| **KBP** | Korean Bone Marrow Donor Program — alternative Korean reference. |
| **arcasHLA** | RNA-seq–based HLA imputation tool (Orenbuch et al.). Used for Korean PTC pool typing (no SNP data needed). |
| **cookHLA** | SNP-based HLA imputation tool (Cook et al. *Nat Commun* 2021, first author = user). Requires SNP data; **deferred for Korean PTC cohorts** (no SNP data available). |
| **DPB1\*05:01** | East-Asian thyroid autoimmunity allele. Shared HT + GD risk. Korean PTC pool 53.2%, Korean baseline ~38-42% (Lee 2014). |
| **DM1 / DM2** | Paper 1 8-gene panel cluster (DM1 = differentiation-low / fusion-enriched; DM2 = differentiation-high). Paper 2 uses DM cluster as substrate label only, no deep mechanism. |
| **Marathon mode** | User-locked 5/4-6/13 6-week manuscript writing rule: no sprint generation, no new analysis except paper-blocking, voice-protected sections always 본인 keyboard. |
| **Voice-protected sections** | Hook / Aim / Discussion mechanism story / Limitations / Cover letter Para 1 / Reviewer Q9 — always 본인 keyboard, never Claude prose. |
| **Q1-Q7** | The 7 advisor decisions in §4 of this packet. |

---

## 11. How to act on this packet

### For Yu professor (advisor)
1. Read §1-3 (scope + why v2 + status) — ~2 minutes
2. Read §4 Q1-Q7 — ~3 minutes
3. Compare against §5 recommended defaults — ~2 minutes
4. Mark each Q1-Q7 as approve / modify / reject
5. Return to user (Seungho) with decisions

### For user (Seungho) post-Yu-meeting
1. Lock advisor decisions in `_LEE_2014_FILL_DECISION_LOCK_2026_05_04.md` style memos per Q
2. Issue execution command per §7 of β plan: `run A2 lookup only` / `prepare A1 forest script only` / etc.
3. Update `MEMORY.md` with G2/G4/G5/G6 status changes
4. Re-run terminology scanner if framing changes ripple to manuscript

### For external readers (papers hub web)
1. §1 + §2 = scope orientation
2. §10 glossary if HLA / Korean cohort terminology unfamiliar
3. §8 related documents for full evidence chain
4. §9 decision timeline for project history context

---

**Advisor packet ready. No execution performed.**
