# Paper 3 ICI — Decision Brief (2026-05-04)

**Purpose:** Single-screen judgment document. Read this offline, decide one of three paths, return with the chosen letter.
**Author:** Seungho Cook
**Date:** 2026-05-04

---

## 1. Status snapshot

| Item | State |
|---|---|
| Paper 3 working title | HLA loss, neoantigen architecture, and immune ecotypes define ICI vulnerability in molecularly dark thyroid cancer |
| Track A (design) | COMPLETED + FROZEN (chmod 444) |
| Track B (execution) | BLOCKED |
| Allowed claim | "ICI vulnerability" / "ICI-readiness" / "immunogenomic prioritization" |
| Forbidden claim | "ICI response predictor" (no thyroid ICI raw RNA-seq exists) |
| Marathon mode | ACTIVE (5/4–6/13, Paper 1 + Paper 2 writing) |
| Paper 1 ship target | bioRxiv 2026-06-13 |
| Paper 2 status | Task A/B/C in progress |
| Paper 4 (Korean GD HLA) | backlog, gated (4/4 entry conditions) |

**Track B unlock requires all 3 of:** (i) Paper 1 bioRxiv submitted, (ii) Paper 2 Task A/B/C closed, (iii) explicit user phrase "Paper 3 Track B 시작" — "고고" / "다 해줘" / "faster" / "추가 분석 더 해줘" do NOT count (per `v17_sprint_vs_marathon_violation.md`).

---

## 2. What is already on disk

**Frozen design (read-only):**
- `reports/paper3_ici/PAPER3_ICI_TRACK_A_BUNDLE.md` — 72 KB, 1236 lines, 8 deliverables consolidated
- `reports/paper3_ici/paper3_ici_*.md` — 8 individual design docs (chmod 444)
- `reports/paper3_ici/paper3_ici_track_a.tar.gz` — 51 KB archive (sha256 `290349f6...8936aa`)
- `reports/PAPER3_ICI_TRACK_A_FREEZE_2026_05_04.md` — root-level freeze notice (chmod 444)

**Editable todo (kept outside freeze):**
- `reports/paper3_ici/todo/paper3_ici_track_b_kickoff_checklist.md`
- `reports/paper3_ici/todo/paper3_ici_controlled_access_checklist.md`
- `reports/paper3_ici/todo/paper3_ici_dataset_accession_verification_checklist.md`

**Memory anchors:** `v19_paper3_ici_track_a.md` (Track A FROZEN), `v19_paper4_GD_backlog.md` (GD renumbered), `paper_numbering_2026_05_04.md` (canonical numbering), `MEMORY.md` index updated.

---

## 3. Decision tree — pick one letter

### (A) — Unlock Paper 3 Track B now

**Action:** Begin Wk 1 of 12-week execution plan immediately.
**Cost:** Halts Paper 1 bioRxiv writing (target 2026-06-13 slips). Halts Paper 2 Task A/B/C. Compute spend ~1500–2500 CPU-hr + 2 GPU-day starts.
**Gate compliance:** FAILS — Paper 1 not submitted, Paper 2 not closed.
**Verdict:** Recommend NO. Choosing (A) breaks marathon and the user's own freeze decision. Only choose if you have a concrete reason that overrides the marathon (e.g., Yu professor demanded it today).

### (B) — Add more Track-A-compatible scaffolding

**Action:** Generate one or more of: ① reviewer-Q anticipation list ② methods M1–Mxx scaffolding ③ statistical analysis plan (SAP) ④ pre-registration template ⑤ cover-letter scaffolding ⑥ code-repo skeleton README ⑦ ethics / data-use statement boilerplate.
**Cost:** Paper 1/2 writing time displaced (~2–6 hr per item depending on depth). No compute spend. Track B remains blocked.
**Gate compliance:** Track B still blocked. Marathon DISPLACED but not violated in the strict (no-compute) sense.
**Verdict:** Acceptable IF the chosen scaffolding genuinely accelerates eventual Paper 3 ship without harming Paper 1/2. Specify which item (①–⑦) when picking.

### (C) — Return to Paper 1/2 marathon (default)

**Action:** Stop touching Paper 3. Resume Paper 1 manuscript / Paper 2 Pillar I v2 forest + H&E-DM1 projection + TCGA validation work.
**Cost:** None to Paper 3 (already frozen). Paper 1 stays on 2026-06-13 ship target.
**Gate compliance:** Full marathon discipline preserved.
**Verdict:** RECOMMENDED. This is the default per the freeze decision and per `v17_marathon_mode_post_pillar1.md`.

---

## 4. Tradeoff summary (1 line each)

- (A) Fast Paper 3 progress, slow Paper 1 ship, freeze decision reversed.
- (B) Modest Paper 3 progress, modest Paper 1/2 cost, freeze decision honored, scope creep risk.
- (C) Zero Paper 3 progress, full Paper 1/2 throughput, freeze decision honored, recommended.

---

## 5. Recommended action

**(C). Return to Paper 1/2 marathon.** Paper 3 design is locked and ship-ready for Track B execution at the post-marathon entry point. Adding more design now does not make Paper 3 ship faster after Paper 1/2 close; it only delays Paper 1.

If you choose (B), pick exactly one of items ①–⑦ and time-box it (≤4 hours). If you choose (A), state the override reason explicitly so the marathon-violation memory can record the exception.

---

## 6. To return with your decision

Reply with one of:
- `(A) — reason: <why override marathon>`
- `(B) — item: <①–⑦>, timebox: <hours>`
- `(C)` (default; no further input needed)

---

Paper 3 Track A frozen. Return to Paper 1/2 marathon.
