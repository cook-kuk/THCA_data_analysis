# AGENTS.md — Codex / agent configuration for THCA data analysis repo

This file is the entry point for Codex (and any other AI coding agent that reads `AGENTS.md`). Claude Code reads `CLAUDE.md` at the same path; both files share the disk-layout / RunPod / marathon rules below. **Read CLAUDE.md as well** — it has the bind-mount layout details that matter for any process writing large outputs.

## Identity + handoff

- Project: Paper 1 (DM1 dark matter / 8-gene panel / fusion-driven thyroid dedifferentiation) primary; Paper 2 (image-DM1) and Paper 3 (ICI vulnerability) tracked.
- Author: Seungho Cook (kukshomr@gmail.com) — Co-PI own first-author target = Nature Communications-tier or above; below-NC Q1 (npj / JCI Insight) is fallback reserve.
- **Latest handoff doc**: [`project/HANDOFF_2026_05_09.md`](project/HANDOFF_2026_05_09.md) — read this first to understand what was just completed and what's open.
- **Recent commit chain (paper9-perturbation-extension-20260506)**:
  - `0bdc243` v14 cross-cohort forest + Q14 + reviewer defense dashboard
  - `982b05e` paper1 Fig 8 mechanism dossier (v6→v13 forensic, dossier pattern)
  - `e190064` paper1 Fig 8 mechanism wiring — Q13 / Supp SX_v13 / Fig 8 cross-ref
  - `9f85fe4` deconv v13 — TDS-16 × MAPK panel (analysis + figure + 7 TSVs)

## Marathon discipline (CRITICAL — applies to all agents)

We are in **marathon mode** 2026-05-04 → 2026-06-13: 6-week manuscript writing window. New analyses must be **paper-blocking only**; default mode for any agent is **scaffolding / infra / data-block / figure** — NOT prose generation.

**Voice-protected sections — author keyboard only, no AI prose generation:**

| Section | File anchor |
|---|---|
| Hook | `project/manuscript_v8/03_introduction.md` + `04_intro_1_1_hook.md` |
| Aim | §Aim within Intro |
| Discussion §3.1 (mechanism interpretation) | `project/manuscript_v8/06_discussion.md` §3.1 |
| Limitations §3.4 | `06_discussion.md` §3.4 |
| Cover Letter Paragraph 1 | `08_cover_letter.md` ¶1 |
| Q9 (mechanism of silencing) | `09_reviewer_qa.md` Q9 — "motivates rather than confirms" boundary |

**Korean enthusiasm phrases** (`고고`, `faster`, `다 해줘`, `대박 고고고`) do NOT override voice-protection. They authorize scaffolding/infra/data work, not voice-protected sprint generate.

**OK to do** (scaffolding/infra):
- New analysis pipelines + plotting scripts under `project/results/p_<topic>_YYYY_MM_DD/`
- Update `02_outline.md`, `05_figure_captions.md`, Q-numbered factual blocks in `09_reviewer_qa.md` (Q1–Q8 / Q10–Q14 only; Q9 is voice-protected)
- HTML dossier pages under `project/papers_hub_2026_05_04/` with auto-deploy to `/var/www/papers/papers_hub_2026_05_04/`
- SUMMARY.md updates for analysis result dirs
- Memory writes (see Memory section below)

**NOT OK** without explicit authorization:
- Generate Hook / Aim / §3.1 / §3.4 / Cover Para 1 / Q9 prose
- Edit voice-protected sections beyond mechanical scaffolding (e.g. inserting `(Figure SX_v14)` cross-references is OK; rewriting sentences is not)

References: `~/.claude/projects/-home-seungho-personal-THCA-data-analysis/memory/v17_marathon_mode_post_pillar1.md` + `v17_sprint_vs_marathon_violation.md`.

## Disk layout (IMPORTANT — bind mounts; mirrors CLAUDE.md)

Root disk is small (123G); a 512G Premium SSD is mounted at `/data`. All large outputs MUST land on `/data`:

```
project/results  ← bind mount /data/thca/repo_results   (xfs on /dev/sdb)
project/data     ← bind mount /data/thca/repo_data
```

You can write to `project/results/...` or `project/data/...` like normal — storage lands on `/data`, paths look identical to git. Verify: `findmnt project/results` should show `/dev/sdb[/thca/repo_results]`.

Pre-existing offload paths (do not re-create on root):
- `project/results/v17_korean/arcasHLA{,_GSE213647}` → `/data/thca/_repo_offload/...`

`/tmp` is on root disk (123G). For STAR / arcasHLA / GNU parallel, use `--tmpdir /data/thca/_tmp --compress`.

## Auto-deploy convention (papers_hub)

- HTML pages live in `project/papers_hub_2026_05_04/`
- Always `cp` to `/var/www/papers/papers_hub_2026_05_04/` on update — bind-mounted live URL serves from there.
- Memory: `~/.claude/projects/-home-seungho-personal-THCA-data-analysis/memory/papers_hub_deploy_target.md`.
- Standard deploy command:
  ```bash
  cp project/papers_hub_2026_05_04/<page>.html /var/www/papers/papers_hub_2026_05_04/
  ```
- For new figures: also mirror to `project/papers_hub_2026_05_04/assets/paper1/` AND `/var/www/papers/papers_hub_2026_05_04/assets/paper1/`.

## Dossier page pattern (validated)

User said "대박" 2026-05-08 about `cancer_vaccine_full_dossier.html`. Default for multi-component project review:

1. Hero (kicker + title + lead + 6–8 stat strip + crumbs)
2. Sticky TOC (left sidebar, sub-anchors, scrollable)
3. TL;DR + headline metric table — every value cites the CSV/JSON it came from
4. Numbered sections (`<span class="num">NN</span>`) — every claim has data table or fig
5. **Honest caveats co-located with results** — Simpson's paradox, LOSO failure, weakness boxes inline (not a separate page)
6. Decision matrix (~9 rows, "Disposition" column gold-highlighted)
7. All-pages card grid linking siblings
8. Sources + paths section
9. Self-contained CSS, dark theme, JetBrains Mono / Cormorant Garamond

References: `paper1_fig8_mechanism_dossier.html` + `paper1_reviewer_defense_dashboard.html` + `cancer_vaccine_full_dossier.html`.

## RunPod conventions

- API key + SSH key: `~/.runpod/config.toml` and `~/.runpod/ssh/RunPod-Key-Go`
- Active pod (2026-05-07): `thca-spark-dm-a6000-v4` (id `uvp9i2r9s6l85y`), RTX A6000 48GB, $0.33/hr.
  - SSH: `ssh -i ~/.runpod/ssh/RunPod-Key-Go -p 20788 root@135.84.176.142` (port/IP via API if changed).
- 50GB `/runpod-volume` is pod-specific persistent (not migratable). `/workspace`, `/tmp`, `/root` are ephemeral.
- See `CLAUDE.md` for full RunPod section.

## Memory (cross-session context)

Claude Code maintains a persistent memory at `~/.claude/projects/-home-seungho-personal-THCA-data-analysis/memory/`. Codex does not share this directory by default, but:

- **Read-only reference**: every memory file is plain Markdown with frontmatter; Codex can `cat` them for project context. The `MEMORY.md` index is the entry point.
- **Most relevant memories for current work** (Paper 1 Fig 8 mechanism layer):
  - `deconv_v5_v12_findings_2026_05_08.md`
  - `deconv_v13_tds16_panel_mapk_2026_05_09.md`
  - `deconv_v14_cross_cohort_forest_2026_05_09.md`
  - `v17_marathon_mode_post_pillar1.md` + `v17_sprint_vs_marathon_violation.md`
  - `feedback_dossier_pattern.md` + `user_style_figure_heavy_briefs.md`
  - `paper_numbering_2026_05_04.md` (Paper 1/2/3/4 canonical numbering)
- **If Codex discovers something new** worth remembering across sessions: write a Markdown memo to that directory with the same frontmatter format (`name`, `description`, `type` ∈ {user, feedback, project, reference}), and add a line to `MEMORY.md` index.

## Repo conventions

- Branch: `paper9-perturbation-extension-20260506` (current). Main branch is `main`.
- Commit message style: lowercase concise headline, then dash-bulleted details, ending with `Co-Authored-By: <agent name>`.
- Never use `--no-verify`, `--no-gpg-sign`, `git reset --hard`, `git push --force` without explicit authorization.
- Pre-commit hooks may run; if they fail, fix the underlying issue and create a new commit (don't `--amend`).
- For data-only commits, prefer staging specific files over `git add -A`.

## Suggested next steps (paper-blocking only)

See `project/HANDOFF_2026_05_09.md` §"Open / suggested next moves" for the current candidate list. As of 2026-05-09:
- v15 candidate A: K2 mini-index Panel-only z-score overlay (Korean cohort generalizability lock; 8-gene only)
- v15 candidate B: spatial Lu 2023 GSE250521 — does MAPK output × Panel-8 anti-correlate spatially within tumors?
- v15 candidate C: DepMap PRISM MAPK-inhibitor sensitivity overlay on DM1-high cell lines (memory `paper11_pancancer_2026_05_08`: PRISM top 15 = ALL MAPK inhibitors)
- Methods dossier: reviewer reproducibility checklist (gene lists, formulas, cohort filters, edge cases)
