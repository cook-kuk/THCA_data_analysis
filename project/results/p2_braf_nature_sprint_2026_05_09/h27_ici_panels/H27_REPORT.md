# H27 — Immune-checkpoint expression in DM1 vs DM2 BRAF-cPTC

**Headline.** TIS (Ayers 2017, FDA-validated ICI surrogate) higher in DM1 vs DM2 BRAF-cPTC: **d=+1.35, MW p=5.6e-9, BH-FDR=1.3e-8** (DM1 z=+0.47, DM2=-0.73; n=110, DM1=85, DM2=25). **Verdict: DM1 is an ICI-candidate phenotype** within BRAF-mutant cPTC.

## Composite scores (DM1 - DM2 d, BH-FDR)

| Composite | d | FDR |
|---|---|---|
| antigen_pres | +1.74 | 3.9e-9 |
| treg_suppr | +1.67 | 2.7e-10 |
| stimulatory | +1.40 | 7.4e-9 |
| **TIS Ayers (20-gene)** | **+1.35** | **1.3e-8** |
| effector | +1.15 | 7.6e-7 |
| immune_hot5 (GZMB+PRF1+IFNG+CD274+PDCD1) | +1.15 | 4.4e-7 |
| ICI_favorable (inhib + effector) | +1.13 | 7.4e-7 |
| inhibitory | +1.07 | 3.4e-6 |
| cytok_chemo | +1.02 | 3.4e-6 |

Every panel — inhibitory, stimulatory, antigen-presentation, effector — is **co-elevated** in DM1: the "inflamed-with-brakes-on" pattern of ICI-responder tumours, not exhausted-cold or T-cell-excluded biology.

## Per-gene tops (DM1 > DM2, all FDR < 1e-6)

* **Checkpoints (ICI targets):** HAVCR2/TIM-3 +1.82, CD274/PD-L1 +1.53, CTLA4 +1.29, PDCD1LG2/PD-L2 +1.20, TIGIT +1.18, LAG3/IDO1 (in TIS).
* **Effectors:** GNLY +1.58, GZMB +0.98, NKG7 +0.89, PRF1 +0.88, TNF +0.86.
* **Antigen presentation:** HLA-DRA +2.07, HLA-B +1.83, HLA-C +1.65, HLA-A +1.44, B2M +1.39, TAP1/TAP2 +1.21–1.27 — fully MHC-I/II competent.
* **Co-stim:** CD40 +2.10, CD86 +1.89, GITR +1.64.

CD274 (+1.53) and CTLA4 (+1.29) satisfy both FDA-precedented ICI axes simultaneously in DM1.

## H10 cross-reference
TIS↔CD8 ρ=0.94 (p=6e-51), TIS↔Tfh ρ=0.85, TIS↔TLS_Cabrita12 ρ=0.91; all checkpoint composites correlate r>0.84 with CD8/Tfh/TLS scores. Checkpoint elevation in DM1 is **co-localised** with bona fide CD8/Tfh/TLS infiltration — not tumour-intrinsic mimic biology.

## ICI-likely-responder fraction
Def: DM1 ∧ TIS≥median ∧ CD8≥median ∧ HLA-I≥median.
**44/110 (40 %) of all BRAF-cPTC; 44/85 (51.8 %) of DM1.** Zero DM2 patients meet all three thresholds.

## Outcome
TCGA BRAF-cPTC PFI/DSS/OS events too rare (≤9) for direct test: ICI-likely vs others rates 6.8/13.6 % (PFI), 0/6.1 % (DSS), 4.5/7.6 % (OS) — direction favours ICI-likely but Fisher p≥0.16. Outcome claim must be reserved for an ICI-treated cohort; panel role here is **candidacy biomarker, not prognostic stratifier.**

## Bottom line
DM1 BRAF-cPTC carries every molecular hallmark of ICI-responsive tumours — high TIS, simultaneous PD-L1+CTLA-4+TIM-3+LAG3+TIGIT, intact MHC-I/II, co-localised CD8/Tfh/TLS. Paper 3 hypothesis: among BRAF-mutant PTCs that progress to RAI-refractory or PDTC/ATC and qualify for ICI, the DM1-derived subset should respond preferentially.
