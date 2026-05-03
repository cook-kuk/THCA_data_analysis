# Marathon snapshot — after image-DM1 NO-GO

**Date:** 2026-05-04 · **Owner:** Seungho Cook
**Marathon window:** 2026-05-04 → 2026-06-13 (6 weeks). bioRxiv 6/13 = critical milestone.
**Trigger event:** Closure battery final NO-GO on image-DM1 (`2026_05_04_image_dm1_final_nogo_decision.md`). ~4–6 weeks of would-be analysis time freed.

---

## 1. Per-paper state

| paper | state | active scope | image-DM1 impact | next concrete action |
|---|---|---|---|---|
| **Paper 1** — 8-gene RAI / DM1 dark-matter sub-stratifier (cancer) | **active**, marathon | DM1 molecular subtype (RNA-seq), TERT recovery v2 (HR=4.33), 8-gene driver-excluded framing, ST supplementary (GSE250521) | **none** (image-DM1 was not Paper 1's contribution) | bib-m3m4 if not done → otherwise voice-hook (author keyboard) → manuscript_v8 / submission/npj prose |
| **Paper 2** — Hashimoto-overlap PTC mechanism | **active**, marathon | Pillar I v2 forest (Korean PTC vs AFND-South-Korea) DONE; GSE286332 PTC+HT signature DONE; TCGA Hashimoto-like generalization DONE; BCR/TLS DONE; DM1 sub-B = NBNR DONE; **Pillar II/III scaffolding TBD** | image-DM1 / pathology-AI angle **dropped**, brief audit pending (`2026_05_04_paper2_post_image_dm1_nogo_status.md`) | Pillar II/III definition + scaffolding; brief audit (delete deprecated wording) |
| **Paper 3** — Korean GD HLA / pan-Asian | **frozen** | gated on Paper 1 bioRxiv + Paper 2 A/B/C + Bundang n>50 + Yu 4/4 consensus | none (image-DM1 was never Paper 3's territory) | **no touch** until gate (`v19_paper3_GD_gating` memory) |
| **Paper 4** — backlog (TBD) | **backlog** | not yet defined | none | parked until at least Paper 1 bioRxiv |

## 2. Marathon Task ledger

### Paper 1 marathon (2026-05-04 → 2026-06-13)

- **Voice-protected (author keyboard only):** Hook · Aim · Discussion 3.1 · Limitations · Cover Letter Para 1 · Reviewer Q9
- **Done (5/4 work):** HM closure (`2026_05_04_paper1_HM_closure_report.md`); bib citation audit; methods prose gap report; numerical consistency audit; supp tables completeness; xref audit; post-closure recheck
- **Pending bib-m3m4:** `2026_05_04_bib_m3m4_straggler_cleanup.md` exists — verify completion or pick up
- **Next:** if bib-m3m4 not done → finish that. Otherwise → voice-hook (author) or supp/methods scaffolding (Claude OK)
- **Venue ladder:** Sci Rep base; Cell Rep Med / JCI Insight reach (post-Yu split re-evaluation)
- **bioRxiv target:** **2026-06-13**
- **Open risks:** 4 outreach drafts NOT sent (`v17_npj_ship_status` memory); bioRxiv slip cascades to Paper 3

### Paper 2 marathon (2026-05-04 → 2026-06-13)

- **Voice-protected (author keyboard only):** Hook · Aim · Discussion 3.1 · Limitations · Cover Para 1 · Q9 (per `STATUS_PAPER2_2026_05_02.md`)
- **Done:** Pillar I v2 forest (5/4 commit); terminology audit Task B (5/4 commit); v1 → v2 deprecation banner; GSE286332/TCGA-Hashimoto/BCR/DM1-sub-B analyses
- **Now newly freed (4–6 weeks):** image-DM1 path closed; capacity available for Pillar II/III definition + supplementary build-out
- **Pending immediate:** brief/manuscript audit per `2026_05_04_paper2_post_image_dm1_nogo_status.md` §4 (author grep + delete-or-reframe deprecated image-DM1 wording)
- **Pending mid-term:** Pillar II/III scoping with Yu (next meeting)
- **Venue:** TBD (candidates: JCI Insight, Cell Rep Med — HT-overlap mechanism + TLS + IGHV + TCGA generalization package)

### Paper 3 frozen — `v19_paper3_GD_gating`

- Gate to entry (4/4 consensus required):
  1. Paper 1 bioRxiv submitted
  2. Paper 2 Pillar I + II + III committed
  3. Bundang FFPE cohort n > 50 received
  4. Yu advisor green light
- Until all 4: **no touch.** No HT/PTC/8-gene wording. Forbidden words list active.

### Paper 4 — backlog

- Not yet scoped. Parked until ≥ Paper 1 bioRxiv.

## 3. Cross-paper guards

| guard | status |
|---|---|
| Paper 1 cancer narrative isolated from HT/Pillar I | active |
| Paper 2 HT-only, no GD/TSAb/exophthalmos | active (`STATUS_PAPER2_2026_05_02.md` forbidden list) |
| Paper 2 image-DM1 / pathology-AI deprecated wording | NEW (this memo + `2026_05_04_paper2_post_image_dm1_nogo_status.md`) |
| Paper 3 frozen, no touch | active |
| No new analyses except paper-blocking (marathon mode) | active |

## 4. What was archived (image-DM1 closure outputs)

All preserved on disk and committed where lightweight; **not deleted** (provides reproducibility + future re-entry baseline if §5 conditions of decision memo are ever met):

```
project/results/03_pathology_poc/
  closure_battery_metrics.tsv             # 26 rows (A/B/D/E/F)
  closure_battery_predictions.tsv.gz      # per-tile preds
  closure_battery_summary.png             # bar plot
  loso_metrics_resnet50.tsv               # 224 + 448 LOSO
  loso_predictions_resnet50.tsv.gz        # per-tile preds
  negative_controls_summary.tsv           # 30 random + housekeeping (resid)
  negative_controls_raw_summary.tsv       # 50 random + housekeeping + RAI real (raw)
  pred_vs_obs_resnet50.png                # 224 scatter
  pred_vs_obs_resnet50_per_slide.png      # 16-slide bar
  tile_metadata*.tsv.gz                   # 8,521 tiles 224/448/672

project/reports/
  pathology_dm1_phaseA_cpu_verdict_2026_05_04.md
  pathology_dm1_closure_battery_2026_05_04.md
  2026_05_04_image_dm1_final_nogo_decision.md           # decision authority
  2026_05_04_paper2_post_image_dm1_nogo_status.md       # P2 scope cleanup
  2026_05_04_marathon_state_after_image_dm1_nogo.md     # this file

repo root:
  CLOSURE_BATTERY_2026_05_04.md           # top-level Korean summary
  PHASE_A_NOGO_REPORT_2026_05_04.md       # top-level Phase A summary
  PATHOLOGY_DM1_FEASIBILITY_2026_05_03.md # initial plan (now superseded)
```

Embeddings (npz files) and large H&E tiles remain in `.gitignore` (regenerable from raw via `project/src/05_pathology_poc/`).

## 5. Recommended next action

**Order of priority (Claude-safe defaults; voice-protected sections still author-only):**

1. **Verify bib-m3m4 completion** for Paper 1 — if `2026_05_04_bib_m3m4_straggler_cleanup.md` indicates incomplete, finish that block.
2. **Otherwise**: Paper 1 voice-protected sections (Hook / Aim / Discussion 3.1 / Limitations / Cover Para 1 / Q9) — **author keyboard**, not Claude.
3. **Otherwise**: Paper 2 brief audit per `2026_05_04_paper2_post_image_dm1_nogo_status.md` §4 — `grep` for deprecated image-DM1 wording, delete or reframe (author-led).
4. **Otherwise**: Paper 2 Pillar II/III scaffolding scoping (await Yu meeting input).
5. **Optional / safe-to-defer**: Paper 1 supplementary methods scaffolding for the closure battery negative-feasibility entry (Claude OK, not voice-protected).

**Do NOT do during this marathon:**
- new image-DM1 analyses
- foundation model retries
- TCGA WSI download
- Azure/RunPod GPU re-spinup for image-DM1
- Paper 3 touches
- new dataset downloads
- voice-protected prose generation by Claude

## 6. Cross-references

- `STATUS_PAPER1_2026_05_02.md`, `STATUS_PAPER2_2026_05_02.md`, `STATUS_PAPER3_2026_05_02.md`
- `2026_05_04_image_dm1_final_nogo_decision.md`
- `2026_05_04_paper2_post_image_dm1_nogo_status.md`
- `CLOSURE_BATTERY_2026_05_04.md`
- Memory pointers: `v17_marathon_mode_post_pillar1`, `v17_sprint_vs_marathon_violation`, `v18_paper2_HT_isolated`, `v19_paper3_GD_gating`
