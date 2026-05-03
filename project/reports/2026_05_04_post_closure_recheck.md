# Paper 1 — Post-Closure Recheck (2026-05-04)

**Author:** Seungho Cook
**Date:** 2026-05-04 (post-crash recovery, 3-group commit, marathon scaffolding/infra session)
**Mode:** Read-only verification of HM closure (`PAPER1_MARATHON_AUDIT_BUNDLE_2026_05_04.md` + `2026_05_04_paper1_HM_closure_report.md`).
**Scope:** 9-point verification. No manuscript file modified.

---

## 0. Commit context

| # | Hash | Files | +Lines | Message |
|---|---|---:|---:|---|
| 1 | `c55d582` | 15 | +2,962 | freeze: paper3 ici track a design bundle |
| 2 | `ef023dd` | 15 | +2,058 | fix: close paper1 marathon audit high medium issues |
| 3 | `25d693c` | 205 | +44,635 | chore: archive spatial and pathology feasibility artifacts |

Branch `main`. Working tree post-commits: 2 stragglers (see §3.10).

---

## 1. Check results — 9-point verification

### 1.1 ✓ PASS — No remaining "96%" mut-neg literal in Paper 1 manuscript files

```
$ grep -n "96%" {OUTLINE, R-PRO, CAP, QA, SUPP}.md
(no matches)
```

The literal "96%" (sub-B mutation-negative claim) is fully eliminated from manuscript files. "95% CI" / "95% Wilson" remain — these are confidence-interval references, completely different semantics; not flagged. ✓

### 1.2 ✓ PASS — sub-B "53/56 (94.6%)" propagated correctly

| File | Line | Excerpt |
|---|---|---|
| `figure_captions_all.md` | 153 | sub-B: 1 BRAF+, 2 RAS+, **53 mutation-negative (53/56, 94.6%)** |
| `figure_captions_all.md` | 190 | **53 mutation-negative (53/56, 94.6%)** |
| `results_R1_R5_prose.md` | 73 | Sub-B contained only 1 BRAF+ and 2 RAS+ of 56 samples **(53/56, 94.6% mutation-negative)** |
| `manuscript_v8_OUTLINE.md` | 170 | sub-B (n=56; **53/56 mut-neg = 94.6%**) |
| `paper_supp_tables_draft.md` | 149 | sub-B \| 56 \| 1 (2%) \| 2 (4%) \| **53 (94%)** \| ... |

SUPP T8 retains rounded "94%" by closure decision (count cells preserved). ✓

### 1.3 ✓ PASS — sub-A 69% / 61% denominator disambiguated

| File | Line | Excerpt |
|---|---|---|
| `manuscript_v8_OUTLINE.md` | 170 | sub-A (n=84; **51/74 mut-tested RAS+ = 69%, 51/84 = 61% of total**) |
| `results_R1_R5_prose.md` | 73 | **51/74 mutation-tested sub-A samples were RAS-positive (69%; 51/84, 61% of total sub-A)** |
| `figure_captions_all.md` | 152 | sub-A: 51 RAS+ **(69% of mutation-tested; 61% of total n=84)** |
| `figure_captions_all.md` | 188 | sub-A (n=84) 51 RAS+ **(69% of mutation-tested; 61% of total n=84)** |
| `reviewer_QA_consolidated.md` | 160 | sub-A: 51/74 mutation-tested RAS+, 1/74 BRAF+ **(69% of mutation-tested; 51/84 = 61% of total sub-A; FVPTC core)** |

Both denominators (74 mut-tested, 84 total) explicit in every cell. ✓

### 1.4 ✓ PASS — Figure 6 + Figure 7 callouts present in R-PRO

| Section | R-PRO line | Callout |
|---|---|---|
| R5b (Lee 2024 replication) | 63 | "(**Figure 6**; Suppl Figure S3; Suppl Table S7)" |
| R5d (DM1 sub-A/B) | 73 | "(Figure 5D; **Figure 7**)" |

M2 closure verified. ✓

### 1.5 ✓ PASS — S1 cohort_assembly.tsv exists with expected schema

```
$ head -1 project/results/p2_pillar1_forest/cohort_assembly.tsv
cohort_id  source_study  n_total  n_valid_HLA  modality  population
DM_cluster_avail  HLA_avail  use_case  accession

$ wc -l project/results/p2_pillar1_forest/cohort_assembly.tsv
6  (header + 5 cohort rows)
```

10 columns × 5 cohorts. H2 build verified. ✓

### 1.6 ✓ PASS — S3 GSE286332 DEGs XLSX

| Property | Value |
|---|---|
| Path | `project/results/p3_gse286332/SuppTable_S3_GSE286332_DEGs.xlsx` |
| Size | 2,236,160 bytes (~2.1 MB) |
| Sheets | `['S3_DEG_PTC+HT_vs_PTC']` |
| Sheet rows | 29,673 (header + 29,672 DEG rows) |

Matches HM_closure_report claim (29,672 rows, padj-sorted). ✓

### 1.7 ✓ PASS — S8 DM1 subcluster XLSX

| Property | Value |
|---|---|
| Path | `project/results/d6p7_dm1_subcluster/SuppTable_S8_DM1_subcluster.xlsx` |
| Size | 3,088,130 bytes (~2.9 MB) |
| Sheets | `['S8a_score_profile', 'S8b_clinical', 'S8c_mut_x_hashi', 'S8d_subBvA_DEGs']` |
| Sheet rows | S8a=8, S8b=3, S8c=3, S8d=51,712 |

All 4 sheets present. S8d row count (51,711 + header) matches expected DEG list size. ✓

### 1.8 ✓ PASS — Paper 3 frozen files unchanged after Group 1 commit

```
$ git log --oneline -- project/reports/paper3_ici/PAPER3_ICI_TRACK_A_BUNDLE.md
c55d582 freeze: paper3 ici track a design bundle  (only commit)

$ ls -l project/reports/paper3_ici/PAPER3_ICI_TRACK_A_BUNDLE.md
-r--r--r-- ... 71684 bytes  (chmod 444 preserved)
```

Same for `PAPER1_MARATHON_AUDIT_BUNDLE` (chmod 444), `PAPER3_ICI_TRACK_A_FREEZE_2026_05_04.md` (chmod 444), `PAPER3_ICI_DECISION_BRIEF_2026_05_04.md` (chmod 444). No mutation. ✓

### 1.9 ✓ PASS — Voice-protected sections untouched

The following sections were **not modified** by either HM closure or recheck (verified by source-of-truth in HM_closure_report.md "Voice-protected sections untouched" attestation, plus this recheck performing zero edits):

- Hook ¶1
- Aim ¶4
- Discussion §3.1 (Krishnamoorthy 2025 → Landa 2016 cite save pending user voice draft)
- Discussion §3.4 Limitations
- Cover letter ¶1
- Reviewer Q9

Per `v17_sprint_vs_marathon_violation` memory rule. ✓

### 1.10 △ NOTE — Working tree stragglers

```
$ git status --short
?? project/.runpod_pods.json
?? project/notebooks_or_scripts/runpod_pod_B_WSI_pathology.sh
```

- `project/.runpod_pods.json` — RunPod runtime state JSON. Likely pid/endpoint state. Recommend gitignore rather than commit.
- `project/notebooks_or_scripts/runpod_pod_B_WSI_pathology.sh` — sibling of runpod_pod_C/D scripts already committed in Group 3. Likely missed because file appeared after staging. Safe to add as a follow-up commit.

Neither blocks marathon work.

---

## 2. Closure status by issue (re-stated post-recheck)

### 🔴 HIGH — both CLOSED

| ID | State | Verification |
|---|---|---|
| H1 sub-B "96%" inconsistency | ✓ closed | §1.1 + §1.2 |
| H2 S1 cohort_assembly.tsv "build TBD" | ✓ closed | §1.5 |

### 🟡 MEDIUM — 4/6 CLOSED, 2 deferred

| ID | State | Verification / Note |
|---|---|---|
| M1 sub-A 69% vs 61% denominator | ✓ closed | §1.3 |
| M2 F6/F7 R-PRO callouts | ✓ closed | §1.4 |
| **M3 Lee2024 GSE213647 bib author placeholder** | ⏳ deferred | Requires GEO/PubMed lookup. Marathon-displacement-safe — handle in next scaffolding session. ~15 min. |
| **M4 Lim2025 GSE286332 bib title verification** | ⏳ deferred | Verify against PMID 41113708. ~5 min. |
| M5 S3 DEG XLSX | ✓ closed | §1.6 |
| M6 S8d sub-A/B DEG XLSX | ✓ closed | §1.7 |

### 🟢 LOW — 9 deferred to W6 cosmetic batch (6/8–6/13)

| ID | Issue | Action |
|---|---|---|
| L1 | TIERA67 / pan-genome ARI rounding (0.90 vs 0.903) | unify or accept rounding |
| L2 | TCGA Hashi top-30% Fisher p (6e-10 vs 6.4e-10) | unify precision |
| L3 | BRAF mRNA Cohen d (-0.04 vs -0.044) | unify to -0.044 |
| L4 | references.bib header "~22 papers" — actual 29 | update header |
| L5 | ★ Unicode in Chu2018JMG note may break BibTeX | replace with ASCII |
| L6 | Pan2025NatComm pages format unusual | verify Nat Comm convention |
| L7 | 4 bib (Tuttle2019, Yi2016, Chen2024, Kim2014) not yet cited | defer to Discussion drafting |
| L8 | "★ exceptional" verb in F2D caption — borderline subjective | user voice polish |
| L9 | Possible Krishnamoorthy 2025 bib gap | decide at §3.1 voice draft |

---

## 3. Voice-protected work (NOT scheduled by Claude)

Per `v17_sprint_vs_marathon_violation`, these are owned by the user and Claude must not pre-draft. Tracked here only as work-remaining inventory:

- Hook ¶1
- Aim ¶4
- Discussion §3.1 (with Landa 2016 cite save per `v17_landa2016_cite_save` memory: 3-layer reverse-causality argument blocking 8-gene 5/8 ATC overlap misattribution to Krishnamoorthy 2025)
- Discussion §3.4 Limitations
- Cover letter ¶1
- Reviewer Q9

Suggested user slot: **6/4–6/13** (W5–W6 of marathon).

---

## 4. Next recommended action

Two safe options. Neither violates marathon mode.

### Option A — `bib-m3m4` (low-risk scaffolding, ~20 min)

1. Look up `Lee2024 GSE213647` author list via GEO landing page or PMID search.
2. Update `references.bib` Lee2024 entry — replace "Y and others" placeholder with full author list.
3. Look up `Lim2025 GSE286332` against PMID 41113708 — verify title matches the PTC-vs-PTC+HT framing used in our manuscript. If mismatch, update bib title or add note.
4. Single-file edit (`references.bib`); diff-friendly; no manuscript file touched.

### Option B — `voice-hook` (user keyboard, voice-protected)

User opens `manuscript_v8_OUTLINE.md` and writes Hook ¶1 directly. Claude does not pre-draft. Earliest start of voice-protected work.

### Other queued items (not recommended next, but on the radar)

- `low-batch` — L1–L9 cosmetic batch. Originally scheduled W6 (6/8–6/13). Pulling earlier is OK if voice-protected work is already underway, but premature now.
- `audit-recheck-figs` — Visual inspection of S8 XLSX sheet contents (open in Excel) to confirm column ordering matches manuscript description. Optional sanity check.
- Followup commit for stragglers (§1.10): add `project/.runpod_pods.json` to .gitignore + commit `runpod_pod_B_WSI_pathology.sh`. <2 min, can batch with `bib-m3m4`.

---

## 5. Marathon discipline attestation

- Voice-protected sections untouched ✓
- No new analysis ✓
- No new data download ✓
- Paper 3 design bundle untouched (chmod 444 preserved) ✓
- Paper 4 not touched ✓
- All Claude work this session = scaffolding / infra / audit / commit ✓

---

Recheck closed. Bundle integrity confirmed. Ready for next user-decision turn.
