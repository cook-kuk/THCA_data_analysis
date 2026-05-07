---
title: "Paper 2 HLA professor-facing exploratory figure page"
date: 2026-05-06
status: "COMPLETE"
---

# Files created or changed

- `project/paper2-hla/index.html`
- `project/scripts/serve_secure.py`
- `project/papers_hub_2026_05_04/index.html`
- `project/three_papers_index.html`

# Page content

Created a professor-facing web report for:

`Paper 2 HLA Pillar I: 6-Allele Exploratory OR/Fisher Forest`

The page includes:

- Hero section with red exploratory-only warning badge.
- Executive verdict card.
- Three summary cards: what the figure shows, why it matters, what it does not prove.
- Strong red metric-mismatch caveat section before the forest plot.
- Source/provenance flow diagram.
- Figure F13 forest explanation.
- Result callouts for DPB1*05:01, A*02:07/B*46:01, DQB1*02:01, C*01:02/DRB1*07:01.
- Primary exploratory statistics table.
- Figure F14 and F15 sections.
- Paper 2 story importance cards.
- Missing validation list.
- Allowed vs forbidden interpretation table.
- Korean professor-facing summary block.
- Next action section.

# Deployment

Copied the route to:

- `project/paper2-hla/index.html`
- `/var/www/papers/paper2-hla/index.html`

Also synced Paper 2 HLA assets into:

- `/var/www/papers/papers_hub_2026_05_04/assets/paper2_hla/`

# Server status

- `thyroid-http.service` restarted and active on `:8012`.
- `papers-web.service` started and active on port `80`.
- The old manual Streamlit process on port `80` was stopped; the port `8001` Streamlit process remains active.

# HTTP check

- `http://40.82.129.113/paper2-hla` -> HTTP 301 to `/paper2-hla/`, then HTTP 200.
- `http://40.82.129.113/papers_hub_2026_05_04/paper1.html` -> HTTP 200.
- `http://40.82.129.113/papers_hub_2026_05_04/paper2_hla.html` -> HTTP 200.
- `http://40.82.129.113/papers_hub_2026_05_04/assets/paper2_hla/F13_exploratory_or_fisher_forest.png` -> HTTP 200.

# Boundary

This page explains the already executed Paper 2 exploratory OR/Fisher/forest.
It does not execute Paper 4 GD meta-analysis, Bundang integration, Paper 3 work, or manuscript prose editing.
