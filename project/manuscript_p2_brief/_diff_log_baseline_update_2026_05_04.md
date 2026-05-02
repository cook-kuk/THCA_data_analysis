---
title: Paper 2 brief — Pillar I baseline update diff log
date: 2026-05-04
session: Paper 2 isolated session (Korean general-pop baseline + Taiwan caveat)
target_file: p2_advisor_discussion.html
backup: p2_advisor_discussion_v1_pre_baseline_update.html (pre-update snapshot)
---

# Paper 2 Pillar I baseline update · diff log 2026-05-04

## 6 changes applied

| # | Spec | Action | Location | Status |
|---|---|---|---|---|
| 1 | Pillar I forest 4-arm (Korean general-pop baseline 추가) | fig2 Plotly script 수정 — 4번째 trace `Korean general-pop (n=680, AFND)` 추가, color muted gold (#B8893C), error_x [33.2, 40.4], DPB1*05:01 only (other 5 alleles null) | ~line 1241 (script) | ✅ |
| 2 | Pillar I 본문 inline sentence 추가 | "Importantly, Korean PTC pool 53.2% also exceeds an independent Korean general-population baseline ..." sentence inserted into "① Pillar I (HLA frequency)" paragraph | line 506 | ✅ |
| 3 | Limitations table (6) + (7) 추가 | Row 6 = Korean baseline source (AFND 3 studies, n=680, I²=0%) · Row 7 = Taiwan caveat (I²=94%, τ²=0.045, n_studies=5, total_n=704) | line 670–674 | ✅ |
| 4 | New Figure 1C — Pan-Asian baseline panel | New `<h3>3.1.2c</h3>` + figure block + `fig1c` Plotly trace · 4 East-Asian populations · Taiwan sensitivity-only x-marker | between §3.1.2 forest + §3.1.3 pooled OR | ✅ |
| 5 | Executive Summary table description | DPB1*05:01 stat-note updated — "두 독립 baseline (Korean general-pop 36.7%, Chu Chinese GD 44.0%) 모두 대비 enrichment — population-stratification artifact 가 아닌 진정한 PTC-specific signal" | line 414 | ✅ |
| 6 | Glossary entry (AFND + I²/τ²) | New amber callout at top of §3.1 — defines AFND, I², τ², RE pooled | line 555 | ✅ |

## NOT changed (per spec)

- Title / hero / TOC headline (Korean baseline 추가 강조 X — 자연스럽게 Pillar I 안에 통합)
- Status badge `STRONG` (이미 격상됨, 변경 X)
- Venue ladder (Cell Rep Med / JCI Insight / Nat Commun) 변경 X
- Pillar II / III / IV / V section 변경 X
- 색상 팔레트 (burgundy / gold / sage) 유지
- 5-Pillar 구조 유지 (Pillar 6+ 추가 X)
- 다른 figure 추가 X (Figure 1C 만)

## Validation checklist

- ✅ forbidden term grep on NEW edits: 0
  - 4 pre-existing hits in P2 (cookHLA SNP-based, DRB1-DQA1-DQB1 trio at lines 591, 670, 674) — these are pre-existing P2 content, not introduced by this session
- ✅ Korean baseline statistical statement: 36.7% [33.2-40.4%], Wald z = 6.58, p = 4.74 × 10⁻¹¹
- ✅ Taiwan caveat: I²=94%, τ²=0.045, n_studies=5, total_n=704
- ✅ AFND 3 studies, total n=680
- ✅ Figure 1 4-arm forest 색상 일관성 (burgundy / gold / sage)
- ✅ Figure 1C Pan-Asian panel Taiwan sensitivity-only label (x-marker, muted gray)
- ✅ Print-friendly 유지 (@media print rules unchanged)
- ✅ Paper 1 / Paper 3 영역 단어 등장 0건 (in NEW edits)
- ✅ Status badge `STRONG` 유지
- ✅ Venue ladder 변경 X
- ✅ Pillar II/III/IV/V section 변경 X

## Cross-paper boundary integrity

- Paper 1 territory terms (BRAF V600E, RET fusion, DM1, 8-gene panel mechanism, LIBRETTO-001, selpercatinib, TIERA67) NOT introduced in this update
- Paper 3 territory terms (Graves' disease mechanism, TSAb, hyperthyroidism) NOT introduced in this update
- Paper 3 boundary phrases (cookHLA SNP-based, DRB1-DQA1-DQB1 trio) — 4 pre-existing hits in P2 reflect cookHLA cross-disease bridge mention (line 591) + limitations (line 670, 674) — these were present BEFORE this session's update and are within Paper 2's own historical narrative scope

## Files

- `p2_advisor_discussion.html` — updated in-place
- `p2_advisor_discussion_v1_pre_baseline_update.html` — pre-update backup
- `_diff_log_baseline_update_2026_05_04.md` — this file
