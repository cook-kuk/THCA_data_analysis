# H30 — Lee 2024 (GSE213647 n=632) comprehensive external replication

**Verdict: Lee 2024 replicates 4 of 6 testable findings; 4 of 10 (driver-stratified survival, TERT interaction) are NOT TESTABLE — GSE213647 carries no per-sample driver / TERT / follow-up annotation in either GEO SOFT (`genotype: NA`, `treatment: NA`) or supplementary MOESM5 (only fixation, library, dx_year, age, sex).**

## Cohort and panels

- 632 samples; 369 tumors (PTC=348, PDTC=5, ATC=16) + 263 normals
- 100 % gene coverage across HT-13, FA-12, MAPK-9, 8-gene differentiation
- DM-surrogate (GMM(2) on tumor 8-gene): DM1=353 / DM2=16, BIC bimodal (Δ=−84.8). DM2 minority concentrates the 16 ATC + 5 PDTC dediff tumors, matching Korean overdiagnosis paradigm.

## Replication of TCGA targets

| Test | TCGA target | Lee 2024 | Verdict |
|---|---|---|---|
| AUC HT-13 vs Hashimoto-GMM | 0.948 (H3) | **0.929** | REPLICATES |
| MAPK-9 d (DM1 − DM2) | +1.81 (H2) | **+1.00, p=0.002** | REPLICATES_SIGN |
| FA-12 d (DM1 − DM2) | −0.64 (H2) | **−1.22, p=0.002** | REPLICATES_SIGN (stronger) |
| Pearson r(HT × FA) tumors | −0.16 (H2) | **−0.006, p=0.91** | REPLICATES_ORTHOGONAL |
| HT-13 d (DM1 − DM2) | +1.58 (H2) | −0.02, p=0.83 | FLIP/NULL |
| AUC HT+FA combo vs DM1 | 0.974 (BRAF-cPTC) | 0.652 | WEAK |

The HT-13 effect-size flat is **not** a contradiction of H3: the same panel hits AUC=0.929 against Hashimoto-GMM in the same cohort. Flat d arises because the n=16 DM-surrogate-DM2 stratum is dominated by ATC/PDTC dediff tumors whose HT distribution overlaps PTC; HT and dediff are orthogonal axes in Lee 2024 (HT-vs-ATC/PDTC AUC=0.39, FA-vs-ATC/PDTC AUC=0.79). The combo-AUC weakness shares the n=16 origin and partially recovers against Hashimoto-GMM (combo AUC=0.728).

## Convergent two-axis signal

- FA axis tracks dediff: AUC(FA, ATC/PDTC vs PTC) = 0.79
- HT axis tracks Hashimoto immune infiltration: AUC(HT, Hashimoto-GMM) = 0.929
- HT × FA tumor r=−0.006 — full orthogonality preserved (TCGA r=−0.16)
- FA × MAPK r=−0.42 — additional Korean-specific MAPK-output / FA inverse coupling

## Not testable

PFI HR=0.21 (DM1 vs DM2 in BRAF-cPTC), OS HR=0.07, DSS HR=8.78, DM1×TERT PFI HR=6.94 — all NOT TESTABLE in Lee 2024. Driver, TERT, follow-up absent. These stay anchored to TCGA H6/H8; Korean survival replication needs the SNUBH retrospective audit (out of sprint scope).

## Files

- `h30_lee2024_panel_scores.tsv` — 632 × HT/FA/MAPK/p8 + DM-surrogate
- `h30_lee2024_replication_table.tsv` — 21 rows with TCGA target + verdict
- `h30_metrics.json` — summary + annotation-gap audit
