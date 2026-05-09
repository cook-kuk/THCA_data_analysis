# H13 — causal mediation: MAPK → 8-gene methylation → HT immune axis

**Verdict:** *Stratum-dependent.* The chain is **forward (MAPK → 8-gene methylation → HT)** in **RAS_like / FVPTC** (almost full mediation, prop = 0.86, 95 % CI [0.43, 2.64]). In **BRAF_like / cPTC** the same forward chain fails (a-path n.s., total effect is *negative* −0.155 and acts directly, not through methylation). Across pooled TCGA-THCA the directional model with best AIC is **MAPK → HT → 8-gene methylation** (AIC 1039 vs 1070 forward / 1402 reverse), i.e. methylation behaves as a *downstream* readout of the MAPK + immune axes rather than a mediator between them.

## 1. Three-layer cohort

| stratum | n with all three layers |
|---|---|
| BRAF_like + cPTC | 363 |
| RAS_like + FVPTC | 91 |
| All TCGA-THCA primary tumors | 509 |

X = MAPK panel z-mean (DUSP4/5/6, SPRY2/4, ETV4/5, PHLDA1, CCND1).
M = 8-gene HM450 mean β (DIO1, FOXE1, NKX2-1, PAX8, SLC5A5, TG, TPO, TSHR).
Y = HT-13 panel z-mean (HLA-DR/DP/DQ, CD79A/B, MS4A1, AICDA, CXCL13, CCR6, IFNG).
Standardised, OLS Baron–Kenny + 5 000-bootstrap CI.

## 2. Mediation results (key columns)

| stratum | direction | total c | indirect a·b | prop mediated [95 % CI] | AIC |
|---|---|---|---|---|---|
| BRAF_like_cPTC | forward MAPK→M→HT | **−0.155 (p=.003)** | +0.031 (n.s.) | −0.20 [−3.55, +0.33] | 791 |
| BRAF_like_cPTC | reverse HT→M→MAPK | −0.155 | +0.185 | −1.20 [−6.77, −0.27] | 1011 |
| BRAF_like_cPTC | alt MAPK→HT→M | +0.045 (n.s.) | −0.108 | −2.39 [−23.2, +18.8] | 799 |
| RAS_like_FVPTC | **forward MAPK→M→HT** | **+0.241 (p=.021)** | **+0.206** | **+0.86 [+0.43, +2.64]** | 190 |
| RAS_like_FVPTC | reverse HT→M→MAPK | +0.241 | +0.170 | +0.70 [−0.23, +4.42] | 256 |
| RAS_like_FVPTC | alt MAPK→HT→M | +0.281 (p=.007) | +0.173 | +0.61 [+0.14, +1.05] | **188** |
| All_TCGA_THCA | forward MAPK→M→HT | +0.066 (n.s.) | +0.186 | +2.83 [−16.8, +19.9] | 1070 |
| All_TCGA_THCA | alt MAPK→HT→M | +0.250 (p<10⁻⁴) | +0.046 | +0.19 [−0.18, +0.37] | **1039** |

Direction interpretation (`a` is X→M, `b` is M→Y | X):

* **RAS_like / FVPTC** — a = +0.28 (p=.007), b = +0.73 (p<10⁻³⁰), direct c′ collapses to 0.035 (p=.64). 86 % of the MAPK→HT effect runs through 8-gene methylation. The chain is real here.
* **BRAF_like / cPTC** — a = +0.045 (p=.39): MAPK output is **decoupled** from 8-gene methylation within this stratum (everything is already silenced). The total MAPK→HT effect is *negative* and stays negative when M is added (c′ = −0.186, p<10⁻⁴): more MAPK output → less HT signal, *directly*. The "alt" model (MAPK→HT→M, AIC 799 vs forward 791) shows HT itself partially mediates MAPK's residual effect on methylation (indirect = −0.108, CI [−0.21, −0.01]).
* **Pooled** — best AIC is alt (1039), worst is reverse (1402). Methylation looks like a *downstream consequence* of both MAPK and immune state, not the bottleneck connecting them.

## 3. DM-polarity audit (mean 8-gene β, master `dm` label)

| stratum | DM1 | DM2 | not_DM | DM1−DM2 (Cohen d) |
|---|---|---|---|---|
| BRAF_like_cPTC | 0.390 (n=87) | 0.256 (n=25) | — | +1.50 |
| RAS_like_FVPTC | 0.362 (n=18) | 0.263 (n=39) | — | +1.21 |
| All TCGA-THCA | 0.385 (n=109) | 0.258 (n=68) | 0.355 (n=332) | +1.34 |

**Apparent contradiction with `dm1_round4` memo (claimed DM1 < DM2, d=−1.75) is a label-flip artifact.** Round4 used `n_DM1=138, n_DM2=358` — i.e. it pooled `not_DM` (n=334) into the comparator AND assigned the LOW-β / high-expression cluster as "DM1". Master / paper-1 standard (`r9_1_master_with_rai.tsv`, n=110/69/334) defines DM1 as the **HIGH-β / 8-gene-silenced** dark-matter group, consistent with the BRAF multimodal report (DM1 β = 0.384, DM2 β = 0.258). All downstream paper-1 narrative should use the master convention; round4's per-gene table needs a sign flip when ported into manuscript figures. Histograms saved to `h13_polarity_histograms.png` confirm bimodality with DM1 right-shifted in every stratum.

## 4. Implications for paper 1 + 2

1. The "MAPK output → 8-gene silencing → HT immune axis" chain is a **RAS-like / FVPTC story**, not a BRAF-like / cPTC one. Within cPTC, methylation is at-floor regardless of MAPK level (no a-path), so the immunity correlate is direct from MAPK.
2. The FVPTC mediation result is *consistent* with deconv v9–v12: MAPK output drives panel hypermethylation, which then sets the immune permissiveness. The mediating role of methylation in this stratum is the cleanest causal evidence we have to date.
3. The polarity audit closes the round4-vs-master sign-flip; future panels must label DM1 = high-β.
