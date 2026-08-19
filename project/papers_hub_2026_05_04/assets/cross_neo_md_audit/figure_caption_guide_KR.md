# CROSS-Neo Figure 설명 가이드

## 전체 메시지

CROSS-Neo는 cheap DL/uncertainty로 후보를 먼저 줄이고, TCR/구조/MD는 마지막 검증/해석 계층으로만 쓴다.

## Fig. A / fig_md24 — End-to-end DL-first recognition-aware funnel

- 보여주는 것: 649 candidates are narrowed by cheap DL, uncertainty, TCR, structure, and MD gates.
- 핵심 메시지: The pipeline turns a broad labeled candidate set into a tiny, auditable wetlab/MD shortlist.
- 쓰는 위치: Opening overview figure
- 주의 문장: Enrichment and prioritization, not proof of immunogenicity.

## Fig. B / fig_md25 — Data usage and evidence lineage matrix

- 보여주는 것: Which data artifacts feed each model, uncertainty, TCR, structure, MD, optimizer, and UI layer.
- 핵심 메시지: Every claim has a traceable source table and a stated boundary.
- 쓰는 위치: Methods/reproducibility slide
- 주의 문장: Lineage map does not validate model performance by itself.

## Fig. C / fig_md26 — Representative candidate evidence board

- 보여주는 것: GADGVGKSAL and HMTEVVRHC across main DL, Bayesian uncertainty, TCR branch, Baker fallback, and completed MD.
- 핵심 메시지: HMTEVVRHC now has the strongest completed MD signal, while GADGVGKSAL remains the more balanced DL+TCR+MD candidate with controls queued.
- 쓰는 위치: Candidate selection/case study slide
- 주의 문장: Candidate evidence does not establish clinical utility.

## Table 1 — Stage-by-stage label composition

- 보여주는 것: Counts of positive/negative rows after each gate.
- 핵심 메시지: Positive fraction rises from 24.8% in the full set to 100% in MD escalation/wetlab shortlist under current settings.
- 쓰는 위치: Funnel performance table
- 주의 문장: Current-label enrichment only; external holdout needed for final claims.

## Table 2 — Threshold operating presets

- 보여주는 것: TP/TN/FP/FN, precision, recall, and F1 for optimized slider presets.
- 핵심 메시지: A no-false-positive preset found 13 labeled positives with 0 false positives in the current table.
- 쓰는 위치: Decision threshold slide
- 주의 문장: Preset is optimized on current labels and must be validated externally.

## Table 3 — Data usage manifest

- 보여주는 것: Where each dataset/result is used and what claims it can support.
- 핵심 메시지: The pipeline is auditable from input label table to final UI.
- 쓰는 위치: Supplement/methods table
- 주의 문장: Documentation table, not performance evidence.

