---
title: "Paper 4 GD HLA web integration full result"
date: 2026-05-06
status: "COMPLETE - backlog registry web page only"
---

# Files created

- `project/papers_hub_2026_05_04/paper4_gd_hla.html`
- `project/papers_hub_2026_05_04/assets/paper4_hla/`

# Figures created

1. `P4_F1_gd_scope_boundary.png`
2. `P4_F2_pan_asian_source_map_or_timeline.png`
3. `P4_F3_gd_evidence_registry_matrix.png`
4. `P4_F4_paper4_gate_flow.png`
5. `P4_F5_source_priority_bands.png`
6. `P4_F6_key_allele_theme_map.png`
7. `P4_F7_execution_block_dashboard.png`
8. `P4_F8_future_extraction_schema_flow.png`
9. `P4_F9_paper2_paper4_separation_wall.png`
10. `P4_F10_future_work_ladder.png`

# Key page sections

- One-line verdict: Paper 4 is Pan-Asian GD HLA backlog.
- Why this is not Paper 2 plus reader glossary.
- Paper 1-style logic flow, analysis story, and evidence ladder.
- Source registry.
- Key allele themes.
- What can be done later.
- What is blocked now.
- Figure panel P4-F1 to P4-F10 with captions and click-to-expand lightbox.
- Expanded registry planning tables.
- Source priority, allele-theme, and execution-block dashboards.
- Safe wording vs unsafe wording.
- Claim boundary table.

# Claim boundary

- Registry: yes.
- Meta-analysis: no.
- Forest: no.
- Bundang integration: no.
- Paper 2 input: no.
- Future Paper 4: yes, after gates.

# Deployment status

Local project page exists and assets render. Page was expanded to match the richer Paper 1-style content pattern: reader glossary, logic flow, analysis story, evidence ladder, safe/unsafe wording, more registry tables, deeper boundary explanation, and clickable image enlargement.

Copied to `/var/www/papers/papers_hub_2026_05_04/`, but public port `:80` currently returns the Streamlit app shell rather than the static THCA secure server. The THCA static server on `:8012` serves the new page correctly.

# HTTP check

- `http://40.82.129.113:8012/papers_hub_2026_05_04/paper4_gd_hla.html` -> HTTP 200, content length 22450.
- `http://40.82.129.113:8012/papers_hub_2026_05_04/assets/paper4_hla/P4_F1_gd_scope_boundary.png` -> HTTP 200.
- `http://40.82.129.113/papers_hub_2026_05_04/paper4_gd_hla.html` -> HTTP 200 but served Streamlit shell, not static HTML.

# Commit proposal

HLA-W2 Paper 4 web:

- `project/papers_hub_2026_05_04/paper4_gd_hla.html`
- `project/papers_hub_2026_05_04/assets/paper4_hla/`
- `project/reports/2026_05_06_paper4_gd_hla_web_integration_full_result.md`

HLA-W3 indexes:

- `project/papers_hub_2026_05_04/index.html`
- `project/three_papers_index.html`

# Next action

Freeze Paper 4 registry until Paper 4 gates clear, or return to Paper 1 Hook.
