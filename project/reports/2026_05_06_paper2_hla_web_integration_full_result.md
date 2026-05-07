---
title: "Paper 2 HLA web integration full result"
date: 2026-05-06
status: "COMPLETE - explanatory web page plus Paper 2 exploratory statistics"
---

# Files created

- `project/papers_hub_2026_05_04/paper2_hla.html`
- `project/papers_hub_2026_05_04/assets/paper2_hla/`
- `project/results/p2_pillar1_forest_v2/paper2_6allele_or_fisher_primary.tsv`
- `project/results/p2_pillar1_forest_v2/paper2_6allele_or_fisher_primary.json`
- `project/results/p2_pillar1_forest_v2/paper2_6allele_exploratory_or_fisher_forest.png`
- `project/results/p2_pillar1_forest_v2/paper2_6allele_exploratory_or_fisher_forest.pdf`
- `project/reports/2026_05_06_paper2_6allele_exploratory_stats_report.md`

# Figures created

1. `F1_paper2_hla_scope_boundary.png`
2. `F2_six_allele_readiness_matrix.png`
3. `F3_source_coverage_heatmap.png`
4. `F4_metric_mismatch_diagram.png`
5. `F5_descriptive_values_not_for_forest.png`
6. `F6_decision_gate_flow.png`
7. `F7_allele_source_role_ladder.png`
8. `F8_table_extraction_status_counts.png`
9. `F9_metric_compatibility_checklist.png`
10. `F10_claim_boundary_dashboard.png`
11. `F11_source_hierarchy_board.png`
12. `F12_no_statistics_signoff.png`
13. `F13_exploratory_or_fisher_forest.png`
14. `F14_exploratory_stats_table.png`
15. `F15_exploratory_metric_inputs.png`

# Key page sections

- One-line verdict: Paper 2 HLA is Hashimoto-overlap PTC only, with exploratory statistics now executed after approval.
- Non-specialist HLA background plus reader glossary.
- Paper 1-style logic flow, analysis story, and evidence ladder.
- Paper 2 vs Paper 4 boundary.
- 6-allele source/provenance summary.
- Exploratory OR/Fisher statistics table.
- Primary baseline source per allele.
- Metric compatibility warning.
- Figure panel F1-F15 with captions and click-to-expand lightbox.
- Expanded source/provenance tables.
- Descriptive-only readiness and status summaries.
- Safe wording vs unsafe wording.
- Claim boundary table.
- Yu decision checklist Q1-Q5.
- Next action options.

# Claim boundary

- Korean baseline source table prepared: yes.
- 6-allele exploratory OR/Fisher/forest: executed after approval.
- C*01:02 fill: source available from In 2015.
- DQB1*02:01 fill: source available from In 2015.
- DPB1*05:01: current v2 AFND South Korea primary; Jung 2023 cross-check.
- Harbin fallback: sensitivity only.
- OR/Fisher/forest: executed as metric-mismatched exploratory statistics.
- Paper 4 meta-analysis: not executed.
- Paper 4 GD: blocked in Paper 2.

# Deployment status

Local project page exists and assets render. Page was expanded to match the richer Paper 1-style content pattern: reader glossary, logic flow, analysis story, evidence ladder, figure interpretation notes, safe/unsafe wording, more tables, and clickable image enlargement.

Copied to `/var/www/papers/papers_hub_2026_05_04/`, but public port `:80` currently returns the Streamlit app shell rather than the static THCA secure server. The THCA static server on `:8012` serves the new page correctly.

# HTTP check

- `http://40.82.129.113:8012/papers_hub_2026_05_04/paper2_hla.html` -> HTTP 200, content length 30029.
- `http://40.82.129.113:8012/papers_hub_2026_05_04/assets/paper2_hla/F1_paper2_hla_scope_boundary.png` -> HTTP 200.
- `http://40.82.129.113:8012/papers_hub_2026_05_04/assets/paper2_hla/F13_exploratory_or_fisher_forest.png` -> HTTP 200.
- `http://40.82.129.113/papers_hub_2026_05_04/paper2_hla.html` -> HTTP 200 but served Streamlit shell, not static HTML.

# Commit proposal

HLA-W1 Paper 2 web:

- `project/papers_hub_2026_05_04/paper2_hla.html`
- `project/papers_hub_2026_05_04/assets/paper2_hla/`
- `project/reports/2026_05_06_paper2_hla_web_integration_full_result.md`
- `project/notebooks_or_scripts/paper2_hla_6allele_exploratory_stats.py`
- `project/results/p2_pillar1_forest_v2/paper2_6allele_or_fisher_primary.tsv`
- `project/results/p2_pillar1_forest_v2/paper2_6allele_or_fisher_primary.json`
- `project/results/p2_pillar1_forest_v2/paper2_6allele_exploratory_or_fisher_forest.png`
- `project/results/p2_pillar1_forest_v2/paper2_6allele_exploratory_or_fisher_forest.pdf`
- `project/reports/2026_05_06_paper2_6allele_exploratory_stats_report.md`

HLA-W3 indexes:

- `project/papers_hub_2026_05_04/index.html`
- `project/three_papers_index.html`

# Next action

Decide whether the metric-mismatched exploratory OR/Fisher/forest belongs in main, supplement, or advisor-only material. Paper 4 remains separated.
