# NeoICI Nature Cancer synthesis

## What this package now shows

- NeoICI is not the strongest universal scorer yet; `KG_GA_evolved` still leads the local validation-like split at AUPRC 0.978.
- `NeoICI_GA_RL_combo_v1` is strong but lower at AUPRC 0.863.
- Industrial locked still favors `BioDarwin_mode_bank_v4` at AUPRC 0.628.

## Combo-replay evidence

- The open combo GEO replay retained 9 paired samples across GSE222011 and GSE255830.
- Median response-like index was 0.033; median top-clone fraction was 0.030.

## Proxy comparison

- NeoICI-style proxy Spearman vs response-like index: 0.817.
- NeoPrecis-style proxy Spearman vs response-like index: 0.417.
- NeoICI-style proxy Spearman vs top-clone fraction: 0.233.
- NeoPrecis-style proxy Spearman vs top-clone fraction: -0.150.

## Write-up rule

- Use NeoICI as the combo-therapy readiness layer.
- Do not claim universal superiority over NeoPrecis or any published benchmark.
- If a figure is needed for the manuscript, the synthesis overview should be the top-level entry point, with the proxy comparison as the supporting panel.
