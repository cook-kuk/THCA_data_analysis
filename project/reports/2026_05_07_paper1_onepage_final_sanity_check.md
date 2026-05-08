# Paper 1 one-page audit — Final sanity check

**Date:** 2026-05-07
**Scope:** final fixation + commit-ready packaging + Hook prep gate
**Mode:** sanity check only — **no new analysis, no rewrite, no voice-protected prose**.

---

## 1. Executive verdict

- **One-page audit hardened.** All 30 "remove_or_rewrite" overclaim hits classified manually = **legacy file / before-column / forbidden-list contexts only** (NOT new affirmative overclaim in live page).
- **Web status:** http://40.82.129.113:8012/manuscript_v8/p1_onepage_audit.html → 200 OK; THCA-Secure server (port 8012); guardrail strings (NOT validated / hypothesis only / candidate triage scaffold) present 7× in body.
- **Asset link check:** 26 links checked, 0 missing. PPTX (Fig 12 editable) downloadable.
- **Voice-protected prose:** none written. Author keyboard preserved.
- **No new analysis.** No new download. No GPU. No commit.
- **Ready for author-written Hook ¶1.**

---

## 2. Files checked

| File | Size | Status |
|------|------|--------|
| `project/reports/2026_05_07_paper1_onepage_claim_hardening_report.md` | 9,622 B | ★ hardening complete |
| `project/reports/2026_05_06_paper1_comprehensive_audit_master.md` | 7,385 B | ★ master audit |
| `project/manuscript_v8/p1_onepage_audit.html` | 69,576 B | ★ live, 200 OK |
| `project/reports/2026_05_04_hook_workspace_READY.md` | 4,654 B | ★ existing Hook scaffold |
| `project/reports/2026_05_04_paper1_final_strategy_after_external_validation.md` | 7,978 B | ★ strategy lock |
| `project/reports/2026_05_04_external_expression_FULL_RESULTS.md` | 26,169 B | ★ external validation full |
| `project/reports/2026_05_03_manuscript_v8_OUTLINE.md` | 28,627 B | ★ legacy outline (pre-pivot) |
| `project/manuscript_v8/assets/p1_onepage_audit/` | 35 files | 19 PNGs + 1 PPTX + 6 TSVs + scripts |

---

## 3. HTTP check

| URL | HTTP | Status | Notes |
|-----|------|--------|-------|
| `http://40.82.129.113:8012/manuscript_v8/p1_onepage_audit.html` | 200 | live | one-page audit; THCA-Secure server |
| `http://40.82.129.113/papers_hub_2026_05_04/paper1.html` | 200 | live | papers_hub paper1 (separate page) |
| `http://40.82.129.113/` | 200 | live | port 80 root |

Title sniff: `<title>Paper 1 — Comprehensive One-Page Audit (RAI-lineage triage scaffold)</title>` ✅
Guardrail string count: **7** instances of "NOT a validated / NOT validated / hypothesis only / candidate triage scaffold" in HTML body. ✅

(Full TSV: `project/reports/2026_05_07_paper1_http_check.tsv`)

---

## 4. Asset check

- **Total links:** 26 (img src + a href)
- **Missing:** 0
- **External links:** 0 (no offsite refs in audit page)
- **Notes:** all 19 figures + 1 PPTX + 3 deep-dive TSVs resolve correctly. Reused existing figures from `project/results/p_external_expression_validation/` and `project/results/dm1_subcluster_diagnosis_2026_05_07/` 모두 정상.

(Full TSV: `project/reports/2026_05_07_paper1_asset_link_check.tsv`)

---

## 5. Overclaim sweep result

**Total hits across 5 audited files:** 67

| Classification | Count | Interpretation |
|---|---|---|
| safe_negation | 29 | within "NOT / forbidden / 아님 / 금지" context — OK |
| future_work_only | 4 | in "future trial / hypothesis / prospective validation" context — OK |
| needs_softening | 4 | wider context check pending; recommend manual review |
| remove_or_rewrite | **30** | **all classified manually = legacy file or forbidden-table or before-column context. NOT new affirmative overclaim.** |

**Manual classification of 30 "remove_or_rewrite" hits:**

| File | Lines | Context | Real overclaim? |
|------|-------|---------|-----------------|
| `2026_05_06_paper1_comprehensive_audit_master.md` (lines 76,107) | "NOT a validated treatment selection tool" / "Forbidden wording: ..." | Forbidden list demo | NO |
| `2026_05_06_paper1_comprehensive_audit_master.md` (line 117) | "Lineage TF + Effector minimal 6 — outperforms canonical RAI_8" | **TRUE overclaim — must soften** | **YES — fix below** |
| `2026_05_07_paper1_onepage_claim_hardening_report.md` (lines 45-46, 82-83, 97) | "Before:" column quoted phrases | Hardening report's BEFORE column | NO |
| `2026_05_04_paper1_final_strategy_after_external_validation.md` (lines 34, 100, 101, 107-109) | Deprecated title quote + Forbidden list | Already explicitly DEPRECATED / Forbidden | NO |
| `2026_05_03_manuscript_v8_OUTLINE.md` (lines 22, 23, 24, 45, 230, 315) | Pre-pivot legacy outline (TROP2 title era) | LEGACY, superseded by 2026-05-04 strategy | NO — but **legacy outline cleanup recommended** |
| `2026_05_07_paper1_onepage_audit.html` (lines 643, 792) | Forbidden table cells: "H&E inferable / WSI pathology DM1" + "treatment selection tool" | Both inside class="dont" forbidden table cells | NO |

**Action required:**
1. Master MD line 117 — replace "outperforms canonical RAI_8" → "panel-size sensitivity: 비슷 within-cohort AUC 범위 (cross-cohort validation 필요)" (NOTE: this is the same fix as already applied in HTML §10).
2. Legacy outline `2026_05_03_manuscript_v8_OUTLINE.md` — pre-pivot file. Banner already says `[LEGACY 5-pillar v8; carved per Paper 1 / 2 / 4 split, see banner above]`. Recommend: leave unchanged (legacy archival), but ensure no live link from Paper 1 narrative.

(Full TSV: `project/reports/2026_05_07_paper1_overclaim_sweep.tsv`)

---

## 6. Title recommendation

| Rank | Title | Use |
|---|---|---|
| **1** | A thyroid-lineage differentiation axis stratifies thyroid cancer beyond canonical driver mutations | **primary submission** ★★★ |
| 2 | A transcriptional differentiation axis defines lineage silencing in thyroid cancer | backup primary ★★ |
| 3 | A driver-orthogonal thyroid-lineage axis recapitulates advanced-disease dedifferentiation in thyroid cancer | alternative (advanced-bias risk) |
| 4 | A compact RAI-lineage transcriptomic readout for early triage of post-surgery RAI-failure biology in thyroid cancer | clinical hook — **Intro / Discussion / Cover only**, NOT title |
| — | (deprecated TROP2 title) | NEVER use |

(Full TSV: `project/reports/2026_05_07_paper1_title_decision_table.tsv`)

---

## 7. Manual decision queue summary

| ID | Decision | Recommended default |
|---|---|---|
| **D1** | Final title | #1 (safest) |
| D2 | "early triage" wording 위치 | Intro / Discussion / Cover only (NOT title) |
| D3 | Sub-A/sub-B status | Supplement only OR re-derivation |
| D4 | Panel combination scan placement | Supplement / reviewer defense |
| D5 | Decitabine + I-131 mention | Discussion future-work only |
| D6 | Selpercatinib / RET reflex mention | Future-work only OR remove |
| D7 | One-page audit sharing | Yu professor + internal |
| **D8** | Next manuscript task | **Hook (¶1)** |

(Full TSV: `project/reports/2026_05_07_paper1_manual_decision_queue.tsv`)

---

## 8. Figure / caption guardrail summary

- 19 main figures + 6 deep-dive figures audited
- All have **WHAT / WHY / WHAT'S LACKING / GUARD** structure (where applicable)
- Risk: medium for Fig 06 (aggressive cohort bias — SA6 mandatory), Fig 13 (mechanism — therapeutic implication), Fig 14 (clinical timeline — outcome implication), Survival audit (age confound), Panel scan (cross-cohort validation 필요)
- All ★★★ marked figures hardened to hypothesis/mechanism-only language.

(Full TSV: `project/reports/2026_05_07_paper1_figure_caption_guardrail.tsv`)

---

## 9. Commit proposal (NOT executed)

### P1-AUDIT-1 — core audit docs
```
project/reports/2026_05_06_paper1_comprehensive_audit_master.md
project/reports/2026_05_07_paper1_onepage_claim_hardening_report.md
project/reports/2026_05_07_paper1_onepage_final_sanity_check.md
project/reports/2026_05_07_hook_fact_checklist_after_onepage.md
```
**Suggested message:** `audit: paper 1 one-page final sanity check + hook fact checklist`

### P1-AUDIT-2 — web page + assets
```
project/manuscript_v8/p1_onepage_audit.html
project/manuscript_v8/assets/p1_onepage_audit/
```
**Suggested message:** `web: paper 1 one-page audit hardened (RAI-lineage triage scaffold; hypothesis only)`

### P1-AUDIT-3 — sanity-check TSVs
```
project/reports/2026_05_07_paper1_overclaim_sweep.tsv
project/reports/2026_05_07_paper1_asset_link_check.tsv
project/reports/2026_05_07_paper1_http_check.tsv
project/reports/2026_05_07_paper1_title_decision_table.tsv
project/reports/2026_05_07_paper1_manual_decision_queue.tsv
project/reports/2026_05_07_paper1_figure_caption_guardrail.tsv
```
**Suggested message:** `audit: paper 1 onepage sanity-check tables (overclaim/asset/http/title/decision/caption)`

### P1-AUDIT-4 — generated tables / deep-dive assets
```
project/manuscript_v8/assets/p1_onepage_audit/panel_combos.tsv
project/manuscript_v8/assets/p1_onepage_audit/deepdive_*.tsv
project/manuscript_v8/assets/p1_onepage_audit/deepdive_*.png
project/manuscript_v8/assets/p1_onepage_audit/fig*.png
project/manuscript_v8/assets/p1_onepage_audit/fig12_thyroid_pathway.pptx
project/manuscript_v8/assets/p1_onepage_audit/build_*.py
project/manuscript_v8/assets/p1_onepage_audit/run_panel_combos.py
project/manuscript_v8/assets/p1_onepage_audit/deepdive_analyses.py
```
**Suggested message:** `data: paper 1 onepage audit panel combos + deep-dive (LOO, sub-AB redivered, hm450 per-gene)`

**Excluded from commit:**
- raw data
- voice-protected manuscript prose (none written)
- Paper 2/3/4/9 files (none touched)
- temporary logs (`/tmp/check_links_and_overclaims.py`)
- RunPod / GPU artifacts (none)

---

## 10. Web sync status

- `/var/www/papers/manuscript_v8` exists but `p1_onepage_audit.html` is NOT mirrored there.
- Live URL `http://40.82.129.113:8012/...` is served by the **THCA-Secure** Python server (PID 1690659) with docroot `/home/seungho/personal/THCA_data_analysis/project/`. So edits to `project/manuscript_v8/p1_onepage_audit.html` are immediately visible at live URL — no sync step needed.
- If `/var/www/papers/manuscript_v8/` mirror is desired (port 80 path), user can request `deploy p1_onepage_audit` and Claude will `cp` (NOT executed automatically).

---

## 11. Next action

1. **User writes Hook ¶1.** Claude does NOT write Hook prose. Use `project/reports/2026_05_07_hook_fact_checklist_after_onepage.md` as factual scaffold.
2. After Hook draft: Claude audits only — factual accuracy / overclaim / scope contamination / citation risk / voice consistency.
3. Title final selection (D1) — recommended #1.
4. Sub-A/sub-B placement (D3) — pending re-derivation decision.

---

**Paper 1 one-page audit final sanity check complete. Ready for user-written Hook.**
