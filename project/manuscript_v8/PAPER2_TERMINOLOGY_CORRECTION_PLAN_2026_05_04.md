---
title: "Paper 2 terminology correction PLAN (sed-level dry-run, no execution)"
date: 2026-05-04
status: PLAN ONLY — no scan execution, no edits. Awaiting user approval per category.
scope: manuscript_v8/*.md + memory/*.md (Paper 2 relevant). Paper 1/3/4 file lists for boundary clarity, NOT for editing in this plan.
purpose: HT vs GD conflation cleanup at terminology level (sed-level word/phrase replacement) before any Paper 2 outline/prose generation.
---

# Paper 2 terminology correction plan

---

## 1. Terms to replace

| Search term (case-insensitive, word-boundary) | Replacement | Rationale |
|---|---|---|
| `autoimmune-PTC` | `Hashimoto-overlap PTC` | Specifies HT (the project substrate), not generic autoimmune (which lumps HT + GD) |
| `autoimmune PTC` (with space, no hyphen) | `Hashimoto-overlap PTC` | Same as above; covers both punctuations |
| `autoimmune-thyroid continuum` | `Hashimoto-overlap context` (or **REMOVE** per sentence; manual review) | Continuum framing implies HT ≡ GD which is incorrect |
| `autoimmune thyroid continuum` | Same as above | |
| `Pan-Asian autoimmune-thyroid` | **REMOVE** or move to Paper 4 reserve | Pan-Asian + autoimmune-thyroid framing was tied to Chu 2018 GD comparison (deprecated) |
| `Pan-Asian autoimmune` (general) | **REMOVE** or replace with `East-Asian Hashimoto-overlap` per sentence context | Manual review per occurrence |

**Replacement rule**: Per category, present diff preview to user before execution. NO automatic apply.

---

## 2. Forbidden terms

These should produce **scan hits** (not auto-replacements). User decides per occurrence whether to remove, rephrase, or move to a different paper's file.

### 2.1 Paper 4 territory (GD)

```
Graves
GD                                  (whole-word; avoid matching "GDT" etc.)
TSAb
TSI
thyroid-stimulating antibody
hyperthyroidism
thyrotoxicosis
exophthalmos
Graves' ophthalmopathy
Graves ophthalmopathy
TED                                  (whole-word; avoid matching "TEDx" etc.)
thyroid eye disease
Bundang Graves cohort
Bundang Graves' cohort
```

### 2.2 Paper 1 territory (cancer driver mechanism deep)

```
TROP2
sacituzumab
H&E-DM1                              (escape & for grep)
image-DM1
TCGA WSI
DM1 cluster mechanism deep
BRAF/RAS driver mechanism deep
TERT promoter kinetics
kinase fusion mechanism
dark matter subtype mechanism
```

### 2.3 Paper 3 territory (ICI / immune)

```
ICI response prediction
HLA LOH
neoantigen
DIAL audit
pan-cancer ICI
C5AR1 anti-PD-1 synergy
```

**Scan output format**:
```
[file_path]:[line_no]:[matching_text]:[paper_territory]:[suggested_action]
```

**Suggested actions** (suggested only, user decides):
- `REMOVE` — delete the sentence/clause
- `REPHRASE` — replace with HT-only language
- `MOVE_TO_PAPER_X` — copy to Paper X file, remove from Paper 2
- `KEEP_AS_BACKGROUND` — leave as 1-line background reference (e.g., DPB1\*05:01 GD context citation per Pillar I v2 spec §7)

---

## 3. Allowed terms

These should NOT be flagged. Whitelist:

```
Hashimoto-overlap PTC
Hashimoto's thyroiditis
HT (in HT context with disambiguation, not "GD/HT" combo)
destructive lymphocytic thyroiditis
antigen presentation
HLA-II
HLA class II
thyroid-lineage dedifferentiation
thyroid lineage dedifferentiation
RAI/thyroid differentiation loss
RAI machinery
BCR/TLS                               (within Paper 2 scope only — see §6)
B cell receptor                       (within Paper 2 scope only)
tertiary lymphoid structure           (within Paper 2 scope only)
mediation model
Korean PTC vs Korean baseline
AFND South Korea baseline
HT-specific tissue context
Hashimoto-like signature
```

**Boundary note**: BCR/TLS terms are allowed only when used to describe the Hashimoto-overlap antigen-driven B cell response within Paper 2 mechanism (i.e., not when describing pan-cancer TLS-ICI response, which is Paper 3 territory).

---

## 4. Dry-run scan plan

### 4.1 Target files

**Primary scope (Paper 2 terminology cleanup)**:
```
project/manuscript_v8/00_title_candidates.md
project/manuscript_v8/01_abstract.md
project/manuscript_v8/02_outline.md
project/manuscript_v8/03_introduction.md
project/manuscript_v8/04_results.md
project/manuscript_v8/05_figure_captions.md
project/manuscript_v8/06_discussion.md
project/manuscript_v8/07_star_methods.md
project/manuscript_v8/08_cover_letter.md
project/manuscript_v8/09_reviewer_qa.md
project/manuscript_v8/13_supplementary_tables.md
project/manuscript_v8/README.md
project/manuscript_v8/_audit_HT_vs_GD_2026_05_02.md           (audit doc itself; review-only, NO edit since it's the catch source)
```

**Secondary scope (memory cleanup)**:
```
~/.claude/projects/-home-seungho-personal-THCA-data-analysis/memory/MEMORY.md
~/.claude/projects/-home-seungho-personal-THCA-data-analysis/memory/v17_paper2_pillar1_forest_strong.md   (PRIMARY conflation source — heavy edit needed)
~/.claude/projects/-home-seungho-personal-THCA-data-analysis/memory/v17_2026_04_30_pivot.md
~/.claude/projects/-home-seungho-personal-THCA-data-analysis/memory/v17_arcasHLA_korean_k2.md
~/.claude/projects/-home-seungho-personal-THCA-data-analysis/memory/v17_D4P2_tcga_hashimoto_generalization.md
~/.claude/projects/-home-seungho-personal-THCA-data-analysis/memory/v17_D5P6_BCR_clonal_TLS.md
~/.claude/projects/-home-seungho-personal-THCA-data-analysis/memory/v17_D6P7_dm1_subB_NBNR.md
~/.claude/projects/-home-seungho-personal-THCA-data-analysis/memory/v17_gse286332_strong_go.md
~/.claude/projects/-home-seungho-personal-THCA-data-analysis/memory/v18_paper2_HT_isolated.md            (already clean, verify)
```

**Excluded scope (DO NOT EDIT in this plan)**:
- `_archive_2026_05_03_sprint_draft.md` — archive, never modify
- `04_intro_1_1_hook.md` — empty W1 voice-first file, never auto-edit
- `_PAPER2_MASTER_VIEW.md`, `_session_origin_paper2_isolated.md`, `_paper2_isolated_resume_response.md` — meta-documents about this audit, scan but no edit
- `paper3_ici/PAPER3_ICI_TRACK_A_BUNDLE.md` — chmod 444 FROZEN, Paper 3 territory
- Memory files for Paper 1/3/4 specific scope (e.g., `v17_landa2016_cite_save.md`, `v19_paper3_ici_track_a.md`, `v19_paper4_GD_backlog.md`) — keep paper-specific content, don't merge

### 4.2 Grep patterns

For each search term in §1 + §2:
```bash
# Case-insensitive, line-numbered, with surrounding context (3 lines before/after)
grep -rn -i -B 3 -A 3 -E '\b<TERM>\b' <target_file_or_directory>
```

For replacement candidates (§1):
```bash
# Show before/after diff per occurrence (sed -n preview, NOT in-place edit)
sed -n 's|<old>|<new>|p' <file>
```

For each forbidden term (§2):
```bash
# Output to scan_hits.tsv with path, line_no, matching_text, paper_territory, suggested_action
grep -rn -i -E '\b<term>\b' <target_files> | awk -v territory='<territory>' '{...}'
```

**Output**: Single consolidated TSV `_terminology_scan_hits.tsv` (located at `project/manuscript_v8/`). NOT auto-applied.

### 4.3 Expected false positives

| Pattern | Likely false positive | Mitigation |
|---|---|---|
| `\bGD\b` | "GD" in cite key (e.g., `Chu2018GD` BibTeX), or "GDP", "PhD", etc. (whole-word boundary should handle most) | Manual review per hit |
| `\bTED\b` | "TEDx", "tested", "wanted" (lowercase) — should be excluded by case-insensitive whole-word | Manual review |
| `\bGraves\b` | "Graves" as a person's name (unlikely in our domain; very rare) | Manual review |
| `autoimmune` (general) | Used outside HT context — could be in BCR/TLS context, immune escape context, etc. — manual review per occurrence | Per-occurrence review |
| `Pan-Asian` | If used outside autoimmune-thyroid context (e.g., "Pan-Asian PTC cohort generalizability") — could be retained | Manual review |
| `TROP2`, `sacituzumab` | If forming part of Paper 1 forward-implication 1-line discussion (Paper 1 has its own scope) | Skip if in Paper 1 file |

### 4.4 Approval flow

1. Dry-run scan execution (user-approved)
2. Output `_terminology_scan_hits.tsv` written, NO file edits
3. User reviews TSV, approves/rejects per category (NOT per file)
4. After category approval, sed-level execution per approved category, with diff preview per file before commit
5. Final file changes summary written to `_terminology_correction_log.md`

---

## 5. Edit policy

### 5.1 No voice-protected prose modification

Voice-protected sections (per `v17_marathon_mode_post_pillar1.md` + `v17_sprint_vs_marathon_violation.md`):

**Paper 1**:
- Hook 첫 줄 (Section 1.1)
- Aim / preview paragraph (Section 1.4)
- Discussion 3.1 first paragraph (mechanism story)
- Limitations (Section 3.4)
- Cover letter Para 1
- Reviewer Q9

**Paper 2** (TBD, likely similar):
- Hook 첫 줄
- Aim / preview
- HT mechanism story (Discussion equivalent)
- Limitations
- Cover letter Para 1

**Edit rule for voice-protected sections**:
- Forbidden terms in voice-protected sections → flag in scan output, **DO NOT auto-replace**
- Suggested action = `MANUAL_REVIEW` (user keyboard required)
- If section is empty (e.g., `04_intro_1_1_hook.md`), skip entirely

### 5.2 No automatic rewrite

- No prose paragraph rewrites
- No sentence-level rephrasing without user approval
- sed-level word/phrase replacements only

### 5.3 sed-level only after user approval

- Each category (§1 replacements + §2.1 GD forbidden + §2.2 Paper 1 forbidden + §2.3 Paper 3 forbidden) requires separate user approval
- Approval = "execute category X" (e.g., "execute §1 only" or "execute §1 + §2.1, hold §2.2 + §2.3")
- Per-file diff preview before commit (`git diff` style or `diff -u` output)

### 5.4 Reversibility

- All edits performed via `git`-tracked workflow
- Pre-edit commit checkpoint: `git add . && git commit -m "checkpoint before terminology correction"`
- Per-category commit: `git commit -m "Paper 2 terminology: §1 replacements applied"`
- Rollback path: `git revert <commit>`

---

## 6. Cross-paper boundary

| Paper | Reserved terms (Paper 2 should NOT use) | If found in non-Paper 2 file | Action |
|---|---|---|---|
| **Paper 1** | DM1/DM2 cluster mechanism deep, BRAF/RAS driver mechanism deep, TERT promoter kinetics, kinase fusion mechanism, TROP2, sacituzumab, H&E-DM1, image-DM1, TCGA WSI | If found in Paper 1 manuscript file (e.g., `04_results.md` 2.3) | KEEP — Paper 1 territory |
| **Paper 3** | ICI response prediction, HLA LOH, neoantigen, DIAL audit, pan-cancer ICI, C5AR1 anti-PD-1 synergy | If found in `paper3_ici/PAPER3_ICI_TRACK_A_BUNDLE.md` (chmod 444 FROZEN) | DO NOT TOUCH |
| **Paper 4** | Graves, GD, TSAb, TSI, hyperthyroidism, thyrotoxicosis, exophthalmos, Graves' ophthalmopathy, TED, thyroid eye disease, Bundang Graves cohort | If found in (future) Paper 4 files | KEEP — Paper 4 territory |

**Cross-paper boundary rule**: Terminology cleanup applies only to **Paper 2 manuscript files + Paper 2-relevant memory entries**. Paper 1/3/4 files are scanned for boundary verification but NOT edited.

**Special handling**:
- `MEMORY.md` index file: contains references to all papers. Edit only the Paper 2 specific lines if terminology issues found there. Other paper lines untouched.
- `_audit_HT_vs_GD_2026_05_02.md`: audit document discussing HT vs GD by name — `Graves`, `GD`, `TSAb` etc. appear deliberately in audit context. Whitelist this file (review-only, no edit).
- `_PAPER2_MASTER_VIEW.md`, `_session_origin_paper2_isolated.md`, `_paper2_isolated_resume_response.md`: meta-documents — same whitelist (review-only).

---

# END PLAN — execution requires user approval per category (§1, §2.1, §2.2, §2.3 separately).
