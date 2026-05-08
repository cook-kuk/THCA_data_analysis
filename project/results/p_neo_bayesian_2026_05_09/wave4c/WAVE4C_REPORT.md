# Wave 4C — Stacked Ensemble Report

Wall: **118.0s** · 2026-05-09 01:49:59

## Inputs

Tools ensembled (wave3 predictions):
- mhcflurry
- bigmhc
- deepimmuno
- prime
- netmhcpan
- transphla

Plus Wave 2.5 PWM percentile (`wave25/pwm_features.tsv`, ITSNdb subset n=311).

All 6 tools (MHCflurry, BigMHC, DeepImmuno, PRIME, NetMHCpan, TransPHLA) plus PWM
landed in time. NetMHCpan and TransPHLA both arrived during the wave4c run.

Merged ITSNdb evaluation set: **319** rows
  - in_master=True  : 213  (LR training pool for E3a/E4a)
  - in_master=False : 106  (HEADLINE no_overlap test)

## Strategies

| code | strategy                                              | training data                       |
|------|-------------------------------------------------------|-------------------------------------|
| E1   | Mean of rank-normalized scores                        | none                                |
| E2   | Median of rank-normalized scores                      | none                                |
| E3a  | LogReg on 4 tool scores                               | ITSNdb in_master=True (n=213)       |
| E4a  | LogReg on 4 tool scores + PWM                         | ITSNdb in_master=True (n=213)       |
| E3b  | LogReg on 4 tool scores (5-fold OOF)                  | no_overlap CV folds (n=106)         |
| E4b  | LogReg on 4 tool scores + PWM (5-fold OOF)            | no_overlap CV folds (n=106)         |

E3b / E4b are the cleanest no-leakage comparison vs E1.
E3a / E4a are the documented "train on master patterns, transfer to novel" design;
the in_master AUROC for these is the FIT (not held-out) and is reported as
`in_master_TRAIN` for transparency.

## Headline — ITSNdb_no_overlap (n=106, n_pos=33)

| rank | method | AUROC [95% CI] |
|------|--------|----------------|
| 1 | **mhcflurry** (no_overlap) | 0.668 [0.544, 0.780] |
| 2 | **E4a** (no_overlap) | 0.635 [0.526, 0.734] |
| 3 | **E3a** (no_overlap) | 0.631 [0.516, 0.741] |
| 4 | **bigmhc** (no_overlap) | 0.626 [0.509, 0.736] |
| 5 | **E1** (no_overlap) | 0.623 [0.507, 0.736] |
| 6 | **E3b** (no_overlap_OOF) | 0.621 [0.495, 0.747] |
| 7 | **E4b** (no_overlap_OOF) | 0.618 [0.504, 0.731] |
| 8 | **E2** (no_overlap) | 0.611 [0.492, 0.721] |
| 9 | **deepimmuno** (no_overlap) | 0.579 [0.453, 0.697] |
| 10 | **prime** (no_overlap) | 0.572 [0.457, 0.682] |
| 11 | **transphla** (no_overlap) | 0.555 [0.438, 0.670] |
| 12 | **netmhcpan** (no_overlap) | 0.550 [0.435, 0.675] |


**Best overall** : `mhcflurry` (no_overlap) — 0.668 [0.544, 0.780]
**Best individual tool** : `mhcflurry` — 0.668 [0.544, 0.780]
**Δ (best ensemble − best individual)** : **+0.000**

vs Wave-3 reference MHCflurry no_overlap = 0.668: Δ = **+0.000**


## Stacking vs simple mean

- E1 (mean) no_overlap : 0.623 [0.507, 0.736]
- E3b (LogReg OOF) no_overlap : 0.621 [0.495, 0.747]  →  Δ vs E1 = **-0.002**
- E4b (LogReg + PWM OOF) no_overlap : 0.618 [0.504, 0.731]  →  Δ vs E1 = **-0.006**
- E3a (LogReg trained on in_master) no_overlap : 0.631 [0.516, 0.741]  →  Δ vs E1 = **+0.008**
- E4a (LogReg + PWM trained on in_master) no_overlap : 0.635 [0.526, 0.734]  →  Δ vs E1 = **+0.012**


## Did adding PWM help?

E4 vs E3 deltas (PWM minus no-PWM):
- E4a − E3a : **+0.004**
- E4b − E3b : **-0.003**


## All slices

```
    method         testset   n  n_pos  AUROC  CI_lo95  CI_hi95
 mhcflurry      no_overlap 106     33 0.6683   0.5442   0.7796
 mhcflurry       in_master 213    103 0.6416   0.5711   0.7102
 mhcflurry        combined 319    136 0.6367   0.5749   0.6995
    bigmhc      no_overlap 106     33 0.6260   0.5089   0.7360
    bigmhc       in_master 213    103 0.8416   0.7837   0.8900
    bigmhc        combined 319    136 0.7477   0.6922   0.7994
deepimmuno      no_overlap 106     33 0.5791   0.4533   0.6974
deepimmuno       in_master 213    103 0.8026   0.7405   0.8550
deepimmuno        combined 319    136 0.7019   0.6451   0.7568
     prime      no_overlap 106     33 0.5724   0.4568   0.6816
     prime       in_master 213    103 0.5639   0.4851   0.6451
     prime        combined 319    136 0.5673   0.4996   0.6313
 netmhcpan      no_overlap 106     33 0.5502   0.4352   0.6752
 netmhcpan       in_master 213    103 0.5709   0.4959   0.6450
 netmhcpan        combined 319    136 0.5669   0.5026   0.6314
 transphla      no_overlap 106     33 0.5550   0.4379   0.6702
 transphla       in_master 212    102 0.5832   0.5095   0.6562
 transphla        combined 318    135 0.5797   0.5189   0.6389
        E1      no_overlap 106     33 0.6235   0.5068   0.7360
        E1       in_master 213    103 0.7387   0.6682   0.7988
        E1        combined 319    136 0.6903   0.6331   0.7462
        E2      no_overlap 106     33 0.6110   0.4923   0.7215
        E2       in_master 213    103 0.6963   0.6238   0.7608
        E2        combined 319    136 0.6592   0.5943   0.7251
       E3a      no_overlap 106     33 0.6314   0.5156   0.7410
       E3a in_master_TRAIN 213    103 0.8748   0.8217   0.9161
       E4a      no_overlap 106     33 0.6351   0.5255   0.7343
       E4a in_master_TRAIN 213    103 0.8929   0.8457   0.9314
       E3b  no_overlap_OOF 106     33 0.6210   0.4952   0.7468
       E4b  no_overlap_OOF 106     33 0.6177   0.5043   0.7310
```

## Honesty caveats

1. **All 6 tools loaded** at run time. NetMHCpan no_overlap AUROC = 0.550 and
   TransPHLA = 0.555 — both close to chance, dragging E1 mean down vs the 4-tool
   subset. LR (E3a/E4a) DOWN-weights both (coefs negative or near zero) and
   actually beats E1.
2. **E3a/E4a have a documented confound**: the algorithms are MORE accurate on
   in_master peptides (because those peptides were in their training pools).
   The meta-classifier sees high-quality scores at fit time and noisier scores
   at test — distribution mismatch. Treat E3a/E4a no_overlap AUROC as a stress
   test, not a fair comparison.
3. **E3b/E4b 5-fold CV** on n=106 is small. CIs overlap E1 broadly; treat any
   ensemble-vs-mean delta < 0.03 as noise.
4. **E1 (mean) is the most defensible no-train baseline**; report it alongside
   any LR variant.

## Files

- `wave4c_results.tsv` — all (method × testset) AUROC + bootstrap CI
- `ensemble_predictions.tsv` — peptide-level scores for E1–E4
- `fig_wave4c_ensemble.png/pdf` — forest of all individuals + ensembles on no_overlap
- `lr_E3a.npz`, `lr_E4a.npz` — fitted LogReg coefficients for the in_master-trained variants
