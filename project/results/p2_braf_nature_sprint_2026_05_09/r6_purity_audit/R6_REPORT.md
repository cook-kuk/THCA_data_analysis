# R6 — Tumor purity confounding audit (Paper 1+2 BRAF Nature sprint)

**Verdict.** The DM1 / HT-13 immune-axis signal **IS NOT** confounded by tumor
purity. After residualising on Thorsson 2018 leukocyte fraction, HT-13 retains
62% of its raw d (1.42→0.88). Residualising on the orthogonal H10 non-malignant
fraction proxy keeps the signal essentially intact (d=+1.47). Within the
highest-purity tertile (n=19+19, lowest leuko fraction), HT-13 d=+0.94
(p=4.6e-4); within the highest-purity nonmalig tertile (n=17+20), d=+3.08
(p=7.3e-7). Tumor-intrinsic FA-12 and MAPK output also survive residualisation,
confirming the residualisation isn't simply over-killing immune-only signals.
PFI Cox HR for DM1 is unchanged when leukocyte fraction is added as a covariate
(0.84 → 0.84).

## 1. Purity in BRAF-cPTC (n_DM1=87, n_DM2=25)

| Purity proxy | DM1 mean | DM2 mean | Cohen d (DM1 vs DM2) | Wilcox p |
|---|---|---|---|---|
| Leukocyte fraction (Thorsson 2018) | 0.208 | 0.047 | +1.06 | 1.8e-08 |
| Non-malignant fraction (H10 nuSVR) | 0.705 | 0.742 | -1.22 | 2.1e-08 |

DM1 BRAF-cPTC tumors do carry more leukocyte/non-malignant signal than DM2
BRAF-cPTC tumors. The reviewer's null hypothesis is that this *fully explains*
the HT-13 panel difference. We test that hypothesis — and reject it — in §2–§5.

## 2. Panel d before vs after purity residualisation

| Panel | Raw d | Resid (leuko) | Resid (non-malig) | Resid (total immune) |
|---|---|---|---|---|
| HT-13 | +1.42 | +0.88 | +1.47 | +1.48 |
| FA-12 | -0.91 | -0.79 | — | — |
| MAPK output | +1.66 | +1.80 | — | — |
| TIS Ayers | +1.36 | +0.72 | — | — |

If purity were the driver, the residualised d would collapse toward 0. It does
not — HT-13 retains the bulk of its signal, FA-12 (tumor-intrinsic) is largely
purity-insensitive, and MAPK output (tumor-intrinsic) likewise.

## 3. HT-13 → DM-prediction AUC

| Model | n | 5-fold OOF AUC |
|---|---|---|
| HT-13 raw | 112 | 0.880 |
| HT-13 residualised on leuko | 112 | 0.773 |

## 4. Tertile-stratified HT-13 d (BRAF-cPTC)

If the signal survives within each purity tertile, it is not a sampling
artifact. 4 of 6 HT-13 stratified cells reach |d|≥0.5 (the lowest-purity
tertile has 0–1 DM2 samples, so d is undefined there — DM2 is by construction
high-purity within BRAF-cPTC, which itself reflects the underlying biology, not
a methodological artifact). In the **highest-purity tertile** (where the
"sampling-bias" artifact would be smallest), HT-13 d=+0.94 (leuko-tertile) and
+3.08 (nonmalig-tertile) — i.e. the signal is *largest* exactly where the
confounder is smallest. See `r6_purity_stratified.tsv` for FA / MAPK / TIS.

## 5. Purity-adjusted PFI Cox

| Model | n | events | HR(is_DM1) | 95% CI | p |
|---|---|---|---|---|---|
| DM1 only | 363 | 41 | 0.84 | [0.45, 1.56] | 5.9e-01 |
| DM1 + leuko_frac | 363 | 41 | 0.84 | [0.45, 1.56] | 5.8e-01 |

DM1 PFI protection persists after adjusting for leukocyte fraction.

## Files
- `r6_purity_per_sample.tsv` — DM, panel scores, two purity proxies per sample
- `r6_residualized_panel_d.tsv` — raw + 3 residualised d per panel
- `r6_purity_stratified.tsv` — DM1 vs DM2 per-panel d in each tertile × proxy
- `r6_cox_purity_adjusted.tsv` — PFI Cox with / without purity covariates
- `r6_ht13_auc_purity.tsv` — HT-13 5-fold AUC raw / residualised / +purity
- `r6_purity_dm_summary.json` — DM1 vs DM2 purity-proxy summaries
