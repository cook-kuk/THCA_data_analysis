# R10 — Anti-PD-1 vs anti-PD-L1 regimen stratification

**DM1-adaptive responder phenotype is regimen-specific to anti-PD-L1 in this 5-cohort pool — but regimen and tumor type are perfectly confounded, so the asymmetry cannot yet be attributed to regimen alone.**

The R5 headline (DM1-adaptive HT-13×HLA-I high → response, n=448 pooled OR=1.67 [1.06, 2.62] p=0.026) is **driven by the anti-PD-L1 stratum** (UC; atezolizumab/avelumab/durvalumab; n=320, OR=1.61 [0.95, 2.71] p=0.077). The anti-PD-1 stratum (mel + UC pembro/nivo; n=123) returns OR=1.35 [0.55, 3.34] p=0.51 (3/4 sign-concordant). Cox OS shows the same asymmetry: in anti-PD-L1 every panel passes (HT-13 HR=0.84 p=0.017; TIS HR=0.77 p=3.6e-4; DM1_inflam HR=0.80 p=0.002, n=320, cohort-stratified); anti-PD-1 returns null HRs (all p>0.4, n=39 with OS, low events).

## Regimen × tumor-type confound — must be honest
PD-L1 stratum = 100% urothelial; PD-1 stratum = 88% melanoma + 12% UC pembro/nivo. The asymmetry could be **regimen-specific** (anti-PD-L1 axis preferentially senses HT-13/HLA-I-driven adaptive resistance), **tumor-type-specific** (UC's TMB makes the axis decisive while mel's baseline immunogenicity saturates it), or **power** (n=123 underpowers a 1.5× OR). Path forward: melanoma anti-PD-L1 (rare; IMspire170 if RNA released); UC anti-PD-1 (CheckMate 275 nivolumab-UC, KEYNOTE-052 pembro-UC).

## 2×2 HLA × HT partition (mechanism cross-check)
Anti-PD-L1: DM1-adaptive 30.1% (37/123) > DM2-classical 21.1% > DM2-adaptive 16.2% > DM1-classical 13.5% — monotone gradient, diagonal OR=1.60. Anti-PD-1: DM1-adaptive 26.5% > DM2-classical 20.0%, but DM2-adaptive 41.7% (n=12, unstable). H19 mechanism (DM1=adaptive=responder) is **directionally preserved in both regimens**; only PD-L1 reaches significance. Consistent with power-loss not mechanistic reversal.

## PD-L1 IHC concordance (IMvigor210, n=298) — most important honest finding
Ventana SP142 IC2+ vs IC0/IC1 OR=2.57 (p=0.0012) — the FDA companion biomarker beats every R5 panel as a univariate responder predictor. HT-13 (median-split) OR=1.35 (p=0.33). After joint logistic adjustment, **HT-13 collapses to OR=1.01 (p=0.96) while IC2+ retains OR=2.57 (p=0.0016)**. HT-13 does NOT add independent responder information beyond PD-L1 IHC in UC. TIS shows the same pattern.

## Implications for Paper 3
- The value-add of HT-13 must be argued **on the OS axis** (anti-PD-L1 UC: TIS HR=0.77 p=3.6e-4; HT-13 HR=0.84 p=0.017; DM1_inflam HR=0.80 p=0.002), not as a standalone responder marker.
- DM1-adaptive composite OR=1.61 is hypothesis-generating; needs replication in melanoma anti-PD-L1 + UC anti-PD-1 to break the confound.
- Limitations must state: HT-13 does not beat IC2+ for response; its independent value is on survival.

Files: `r10_regimen_stratified_meta.tsv`, `r10_per_cohort_within_regimen.tsv`, `r10_regimen_cox_os.tsv`, `r10_2x2_partition.tsv`, `r10_pdl1_ihc_concordance.tsv`, `r10_pdl1_ihc_HT_crosstab.tsv`, `r10_headline.json`.
