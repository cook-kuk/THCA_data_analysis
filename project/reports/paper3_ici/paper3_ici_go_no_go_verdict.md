# Paper 3 ICI — Go / No-Go Verdict

**Working title:** HLA loss, neoantigen architecture, and immune ecotypes define ICI vulnerability in molecularly dark thyroid cancer
**Track:** A — design only. Verdict below is **conditional GO** for Track B execution at the post-marathon entry point.
**Author:** Seungho Cook
**Status date:** 2026-05-04

**Claim guard:** Below verdict assumes Paper 3 is framed as ICI vulnerability / readiness / immunogenomic prioritization. If the user later wants to claim "ICI response prediction" without thyroid ICI-treated raw RNA-seq, verdict drops to NO-GO and must be re-decided.

---

## 1. Module-level verdicts

| Module | Verdict | Rationale | Critical dependency |
|---|---|---|---|
| **A — Bulk immune ecotype** | **GO** | TCGA-THCA + GSE76039 + microarray cohorts give ≥600 samples after pooling; ssGSEA + NMF mature; thyroid ecotypes literature-supported. | None blocking. |
| **B — scRNA atlas** | **CONDITIONAL GO** | 6 candidate cohorts, total ~100–300K cells. PDTC/ATC scRNA <5K — atlas viable for state annotation, not for de novo dedifferentiation trajectory. Mitigation framed in atlas plan §1. | At least 3 of 6 cohorts must be accessible at Wk 5. |
| **C — HLA & neoantigen** | **GO (degraded scope)** | TCGA-THCA paired BAM + GDC MC3 MAF available; LOHHLA + NetMHCpan well-established; GSE76039 paired-normal access uncertain per sample; Yoo 2019 may be EGA-controlled. | TCGA-THCA dbGaP controlled-access for paired BAM by Wk 1 of Track B. |
| **D — DIAL audit** | **GO** | Pan-cancer ICI cohorts mostly accessible (IMvigor210 R-package; Hugo / Riaz / Gide GEO; Liu / Cho controlled). Statistical machinery already in `v18_agentic_research framework`. | At least 3 of 5 Tier 1/2 melanoma + UC cohorts accessible by Wk 7. |
| **E — Integrated readiness score** | **GO conditional on D** | Score is well-defined; weighting calibration deferred. | DIAL audit must produce ≥4 PASS signatures for the score to be meaningful. |

**Overall verdict — CONDITIONAL GO.** Paper 3 is feasible. The conditions are dataset-access events that map to the 12-week plan's Wk 1 / Wk 5 / Wk 7 gates, not analytical uncertainties.

---

## 2. Hard gates (must close before any Track B compute spend)

| Gate | Condition | Action if not met |
|---|---|---|
| **G1** | Paper 1 bioRxiv submitted | Track B blocked. Marathon mode preserved per `v17_marathon_mode_post_pillar1`. |
| **G2** | Paper 2 Task A/B/C closure | Track B blocked. |
| **G3** | TCGA dbGaP access (or confirmation that GDC open access is sufficient for our scope) | Module C drops to "RNA-seq HLA only via arcasHLA + MAF-based neoantigen using MC3 SNV calls without LOHHLA". Headline figure 4B (HLA LOH) downgraded to "future-work" if G3 fails. |
| **G4** | At least 3 of 6 scRNA cohorts accessible | Module B drops to "single-cohort sub-state map" with reduced replication; figure 3 panels narrow. |
| **G5** | At least 3 of 5 Tier 1/2 pan-cancer ICI cohorts accessible | Module D narrows; if <3, DIAL is reported as plan-only and Module E uses unaudited signatures with explicit caveat. |
| **G6** | User explicit "Track B 시작" — not "고고" / "다 해줘" | Per `v17_sprint_vs_marathon_violation.md`. |

---

## 3. Risks that DO NOT trigger NO-GO

- Yoo 2019 inaccessible — replace with TCGA aggressive subset; degrade Asian-validation claim, do not block.
- Cho NSCLC / Kim gastric inaccessible — drop, retain Tier 1 melanoma + UC anchors, do not block.
- HCA thyroid reference unavailable — use within-study normal cells, do not block.
- Class II HLA calls noisy — frame as exploratory, do not block.

---

## 4. Risks that WOULD trigger NO-GO (kill switches)

| Trigger | Reason |
|---|---|
| K1 — User attempts to claim "ICI response predictor in thyroid" without thyroid ICI-treated raw RNA-seq | Mis-claim invalidates paper. Hard stop. |
| K2 — Paper 1 ship slips past late August 2026 such that marathon mode is still active | Track B must remain blocked; no execution until ship. |
| K3 — DIAL audit shows zero PASS signatures | Module E has no inputs. Pivot Paper 3 to a feasibility-only / negative-result paper, but this is a different paper. |
| K4 — Paper-boundary contamination — Paper 3 prose pulls Paper 2's PTC+HT TLS / AICDA / BCR claims as Paper 3's own findings | Reviewer scope-rejection certain. Hard stop until rewritten. |
| K5 — TCGA-THCA driver-call disagreement vs Paper 1 ETL produces inconsistent dark-matter rosters | Recompute against Paper 1 source-of-truth; if cannot reconcile, halt Module C/E until resolved. |
| K6 — User pivots priority again (e.g., back to GD or to a new Paper 5) | Re-decide stack; do not silently continue. |

---

## 5. Decision dependencies on Paper 1 / Paper 2

- Paper 1 ship (target 2026-06-13) frees marathon constraint per `v17_marathon_mode_post_pillar1`.
- Paper 2 Task A/B/C closure clarifies which PTC+HT analyses are Paper 2 and which are off-limits to Paper 3.
- Yu professor agreement (per cross-paper boundary discipline in `v19_paper3_ici_track_a.md`) before Paper 3 enters compute.

---

## 6. Headline summary

Paper 3 (ICI vulnerability in molecularly dark thyroid cancer) is **feasible in design** with **CONDITIONAL GO** for Track B. The five modules are mature; the dependencies are dataset-access events and the Paper 1/2 gating events, not analytical uncertainties.

The only framing in which the verdict reverses to NO-GO is over-claiming "ICI response predictor in thyroid" without thyroid ICI-treated raw RNA-seq. Track A's claim guard makes this explicit and binds the Track B prose.

**Recommendation:** lock Track A artifacts as the paper's design freeze; do not begin Track B compute until G1–G6 close.

---

Track A completed. No marathon violation.
