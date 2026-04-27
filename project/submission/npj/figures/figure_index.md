# Figure index — npj submission

## Main figures (composite 4-panel)

- **Figure 6 ★ npj headline** — `Fig6.html` (4-panel: ROC + feature importance + cell-line scatter + PDTC reframe)

## Component figures (individual panels)

- `AMP4_rai_decision_roc.html` — 8-gene LogReg / RF ROC (CV AUC 0.954 / 0.975), BRAF baseline 0.822
- `AMP4_feature_importance.html` — RF feature importance (TPO 0.27, DIO1 0.21, FOXE1 0.14)
- `FIX1_celline_dm_scatter.html` — CCLE 13 lines, BRAF 4/4 in DM1 cluster
- `FIX1_drug_volcano_v2.html` — PRISM 19Q4 mechanism-class enrichment

## Deferred to next sprint

- Figure 1–5, Figure 7 multi-panel composites (require Plotly subplot composition for: driver landscape, marker heatmap, trajectory, orthogonality, hot/cold)
- Supplementary Figure S1–S15

## Notes

- All HTMLs use Plotly CDN, are interactive, and are platform-portable.
- PNG / PDF export requires `kaleido` package; deferred to next sprint.
- DM1 colour code: orange (#FF7F0E); DM2 colour code: blue (#1F77B4).
