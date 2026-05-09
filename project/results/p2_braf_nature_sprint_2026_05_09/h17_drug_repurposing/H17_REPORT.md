# H17 — connectivity-score drug repurposing (BRAF_like / cPTC, N=41)

DM1 n=26, DM2 n=15. 8-gene panel excluded. Signature export = top 200↑(DM1) + 200↓(DM2-up) at FDR<0.05; connectivity score uses the **full** DE table (50,112 genes) ∩ drug KB. Score = weighted Σ sign(d_DM1>DM2)·impact_drug, range [-1, +1]; **+1 = drug perfectly reverses DM2 toward DM1**.

## Top 3 drug candidates (REVERSE DM2 → DM1)

| rank | drug | score | n∩ | MoA | MAPK-class confirmation |
|---:|---|---:|---:|---|---|
| 1 | **IFNG-proxy** | +1.000 | 28 | IFN-γ signaling activator (IFN-γ / TLR3 / TLR9 agonist) | **not** a MAPK inhibitor — operates orthogonal IFN-γ axis |
| 2 | **TLR3-agonist** | +1.000 | 11 | TLR3 agonist (poly-I:C-like) | **not** a MAPK inhibitor — operates orthogonal IFN-γ axis |
| 3 | **anti-PD1** | +0.813 | 10 | PD-1 immune checkpoint blockade | **not** a MAPK inhibitor — operates orthogonal IFN-γ axis |

## All H17 drugs (sorted reverser → opposing)

| drug | score | n∩ | FDR<.05 hits | MoA | verdict |
|---|---:|---:|---:|---|---|
| IFNG-proxy | +1.000 | 28 | 27 | IFN-γ signaling activator (IFN-γ / TLR3 / TLR9 agonist) | reverser |
| TLR3-agonist | +1.000 | 11 | 10 | TLR3 agonist (poly-I:C-like) | reverser |
| anti-PD1 | +0.813 | 10 | 10 | PD-1 immune checkpoint blockade | reverser |
| decitabine | +0.189 | 19 | 18 | DNMT inhibitor | reverser |
| azacitidine | +0.148 | 14 | 14 | DNMT inhibitor | reverser |
| vorinostat | -0.051 | 22 | 21 | pan-HDAC inhibitor | neutral |
| JQ1 | -0.149 | 16 | 13 | BET bromodomain inhibitor (MYC-axis) | opposing |
| dabrafenib | -0.212 | 29 | 23 | BRAF V600E inhibitor | opposing |
| vemurafenib | -0.212 | 32 | 26 | BRAF V600E inhibitor | opposing |
| romidepsin | -0.251 | 14 | 13 | class-I HDAC inhibitor | opposing |
| OTX015 | -0.276 | 12 | 10 | BET bromodomain inhibitor (birabresib) | opposing |
| selumetinib | -0.292 | 32 | 26 | MEK1/2 inhibitor | opposing |
| encorafenib | -0.300 | 24 | 18 | BRAF V600E inhibitor (3rd gen) | opposing |
| cobimetinib | -0.318 | 21 | 16 | MEK1/2 inhibitor | opposing |
| ulixertinib | -0.334 | 22 | 16 | ERK1/2 inhibitor | opposing |
| trametinib | -0.356 | 30 | 24 | MEK1/2 inhibitor | opposing |
| GSK2879552 | -0.369 | 6 | 6 | LSD1 inhibitor (epigenetic re-diff) | opposing |
| AZD-0364 | -0.404 | 27 | 21 | ERK1/2 inhibitor (MAPK-pathway blockade) | opposing |
| tanespimycin | -0.549 | 11 | 9 | HSP90 inhibitor (clears BRAF V600E) | opposing |
| sorafenib | -0.778 | 15 | 13 | VEGFR/RAF/PDGFR multikinase | opposing |
| lenvatinib | -0.782 | 15 | 13 | VEGFR/FGFR/RET/KIT multikinase | opposing |
| FK866 | -1.000 | 7 | 6 | NAMPT inhibitor (NAD+ depletion) | opposing |
| everolimus | -1.000 | 8 | 7 | mTORC1 inhibitor | opposing |
| cabozantinib | -1.000 | 11 | 10 | MET/VEGFR/RET multikinase | opposing |

Verdict counts: **5 reverser** (score>0.10), **18 opposing** (score<-0.10). MAPK inhibitors are OPPOSING, not reversers — see below.

## Module direction within BRAF-cPTC (DM1 minus DM2)

| module | n genes | mean Cohen d | FDR<.05 frac | interpretation |
|---|---:|---:|---:|---|
| TDS_residual | 8/8 | -1.02 | 50% | DM2>>DM1 |
| IFNG_response | 17/18 | +1.86 | 88% | DM1>>DM2 |
| MAPK_output | 11/11 | +1.31 | 45% | DM1>>DM2 |
| Bcell_TLS | 10/10 | +0.64 | 50% | DM1>>DM2 |
| MYC_axis | 7/7 | +0.92 | 43% | DM1>>DM2 |
| ECM_aggressive | 8/8 | +2.58 | 100% | DM1>>DM2 |

**Interpretation.** Within BRAF-cPTC, DM1 = immunogenic + MAPK-active + ECM-aggressive arm; DM2 = immune-cold + thyroid-residue arm. Per H6 the DM2 arm carries PFI HR=5.91. High-risk DM2 phenotype is **immune-cold despite preserved differentiation** — risk reverses by inducing inflammation, NOT by adding re-differentiation.

## H17 × H9 PRISM cross-reference (MAPK class)

| drug | PRISM d (kills DM1-high) | H17 score (reverses DM2→DM1) | verdict |
|---|---:|---:|---|
| encorafenib | -0.60 | -0.300 | opposing |
| AZD-0364 | -0.09 | -0.404 | opposing |

Both axes negative: BRAFi/ERKi selectively kill DM1-high (MAPK-active) lines (PRISM viability) but do NOT reverse the DM2 signature toward DM1 (H17 connectivity). MAPK inhibitor will not rescue DM2 high-risk arm; ICI / IFN-γ axis is the DM2-pertinent intervention.

## Druggable DEG targets (top-30 DM1↑ + top-30 DM2↑)

- **TACSTD2/Trop-2** (DM2-up): sacituzumab govitecan ADC — clinically active TNBC/UC; thyroid basket
- **NAMPT** (DepMap DM1-low d=-0.45): FK866, KPT-9274 — synthetic-lethal in DM2
- **MET** (DM2-up): crizotinib / capmatinib / tepotinib
- **MYC** (DM1-up proliferation feedback): BETi (JQ1, OTX015, BMS-986158); MYCi975
- **STING1** (DM1-up): ADU-S100 / MIW815 — adjuvant to anti-PD1
- **CDK4/6** via CCND1: palbociclib / ribociclib / abemaciclib

## Bottom line

**DM2 BRAF-cPTC reversal modality = IFN-γ axis activation, not MAPK blockade.** Top reversers (IFN-γ proxy +1.0, TLR3 agonist +1.0, anti-PD1 +0.81) push the IFN-γ-response module that DM1 has and DM2 lacks. DNMT inhibitors (decitabine +0.19, azacitidine +0.15) are second-tier via MHC-II + AICDA induction. BRAFi/MEKi/ERKi score -0.21 to -0.40 (opposing) — they shut off DUSP4/5/PLAU/FN1 which are DM1-up here. H6 (DM2 PFI HR=5.91) calls for immunotherapy, not MAPK targeting, in this stratum.

**Caveats.** (1) Drug→gene KB curated from published pharmacology, NOT live LINCS L1000 (clue.io auth-walled, no local L1000 dump). `h17_dm2_signature.tsv` upload-ready for external L1000 run. (2) N=41 cohort: top |d| stable (ICAM1 d=+4.63); tail fold-fragile. (3) Score = sign-of-effect, not potency — PRISM (H9) orthogonal viability axis.

## Files

- `h17_dm2_signature.tsv` — top 200↑(DM1) + 200↓(DM2-up), FDR<0.05, ranked by |d| (clue.io upload-ready)
- `h17_connectivity_top50.tsv` — full drug × score table (30 drugs, sign + intersection)
- `h17_targetable_dm2_genes.tsv` — per-gene druggability for top 30+30 DEGs
- `h17_module_directionality.tsv` — TDS / IFNG / MAPK-output / TLS / MYC / ECM module direction
- `h17_h9_prism_crossref.tsv` — H17 × H9 PRISM cross-reference
- `run_h17_drug_repurposing.py` — full pipeline
