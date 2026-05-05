---
title: "Paper 2 and Paper 4 HLA professor-facing pages"
date: 2026-05-06
status: "COMPLETE"
---

# Executive verdict

Paper 2 and Paper 4 HLA are now both represented as professor-facing web pages with explicit scientific flow, red caveat boxes, figure-by-figure interpretation, tables, limitations, and next-step plans.

Paper 2 focuses on Korean PTC / Hashimoto-overlap exploratory HLA signal.

Paper 4 focuses on GD/Graves HLA. The executable Paper 4 analysis currently available locally is the Chu 2018 Han Chinese GD anchor forest. Full Pan-Asian pooled meta-analysis still requires source-level extraction from Liao 2022, Shin 2019, Park 2005, Chen 2011, Ueda 2014, and historical sources.

# URLs

- `http://40.82.129.113/paper2-hla`
- `http://40.82.129.113/paper4-hla`

# Paper 2 content

- Existing F13/F14/F15 exploratory OR/Fisher figure explanation retained.
- Added explicit Paper 2 vs Paper 4 red-box section.
- Added explanation that shared HLA allele names do not make the two analyses interchangeable.
- Added link from Paper 2 page to Paper 4 page.

# Paper 4 content

Created `project/paper4-hla/index.html`.

Sections include:

- Executive verdict.
- Critical red box: registry is not meta-analysis.
- Why this analysis is being done.
- Chu 2018 GD anchor forest.
- Source extractability matrix.
- Paper 2/Paper 4 bridge-not-input figure.
- Validation gap ladder.
- Primary Chu 2018 anchor table.
- Broader Paper 4 source table.
- What is still missing.
- Allowed vs forbidden interpretation.
- Korean professor-facing summary.
- Next action plan.

# Paper 4 assets created

- `project/papers_hub_2026_05_04/assets/paper4_hla/P4_F11_chu2018_gd_anchor_forest.png`
- `project/papers_hub_2026_05_04/assets/paper4_hla/P4_F12_source_extractability_matrix.png`
- `project/papers_hub_2026_05_04/assets/paper4_hla/P4_F13_paper2_paper4_bridge_not_input.png`
- `project/papers_hub_2026_05_04/assets/paper4_hla/P4_F14_validation_gap_ladder.png`

# Paper 4 result table

- `project/results/paper4_gd_hla/paper4_chu2018_gd_anchor_forest.tsv`
- `project/results/paper4_gd_hla/paper4_source_extractability.tsv`

# Script created

- `project/notebooks_or_scripts/paper4_gd_hla_professor_assets.py`

# Deployment

Synced to:

- `/var/www/papers/paper2-hla/index.html`
- `/var/www/papers/paper4-hla/index.html`
- `/var/www/papers/papers_hub_2026_05_04/assets/paper2_hla/`
- `/var/www/papers/papers_hub_2026_05_04/assets/paper4_hla/`

Services:

- `thyroid-http.service`: active.
- `papers-web.service`: active.

# HTTP checks

- `http://40.82.129.113/paper2-hla` -> 301 to `/paper2-hla/`, content verified.
- `http://40.82.129.113/paper4-hla` -> 301 to `/paper4-hla/`, content verified.
- `P4_F11_chu2018_gd_anchor_forest.png` -> HTTP 200.
- `P4_F14_validation_gap_ladder.png` -> HTTP 200.

# Boundary

No Paper 3 work was performed. No manuscript prose was edited. Bundang validation was not executed because no local Bundang GD phenotype table was available for analysis in this workspace. Full Paper 4 random-effects Pan-Asian meta-analysis remains the next step after source-level extraction.
