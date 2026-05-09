# H19 — TMB / neoantigen / immune-evasion in DM1 vs DM2 BRAF-cPTC

**Stratum**: TCGA-THCA, BRAF × cPTC. DM1 n=87, DM2 n=25 (true-DM2; `not_DM` band excluded). Sources: track29 T01 + Thorsson 2018 + TCGA z-score RNA + K2 arcasHLA.

## Verdict — DM1 immune-hot via preserved HLA + IFN-γ axis, NOT via TMB

DM1 BRAF-cPTC has **lower** mutation burden than DM2 (d = −0.58, p = 1e-3) yet a **dramatically preserved** antigen-presentation apparatus. Mechanism is **HLA preservation + adaptive IFN-γ resistance**, not neoantigen excess and not allelic priming.

## Layers

### (1) TMB — DM1 is colder, not hotter
| feature | d (DM1−DM2) | p |
|---|---:|---:|
| Thorsson nonsil/Mb | **−0.58** | 1.1e-3 |
| n_nonsyn_SNV | **−0.57** | 6.3e-4 |
| n_binding_pMHC | **−0.72** | 7.7e-4 |

DM1 has fewer mutations and fewer raw HLA-I binders. TMB-driven recognition ruled out.

### (2) HLA / antigen-presentation — preserved in DM1, lost in DM2
All d > 0, all p < 1e-7: HLA-I composite **+1.73**; HLA-II **+1.97** (DPB1 +1.83, DRB1 +2.04); B2M +1.36, TAP1/TAP2 +1.21; STAT1 +1.63, JAK1/2 +0.47–0.60 (IFN-γ pathway intact); CD8A +0.86, GZMA +1.21, CXCL9/10/13 +0.78–0.99 (cytotoxic infiltrate).

### (3) Immune-escape concentrated in DM2
Escape = HLA-I-z<−0.5 OR B2M-z<−0.5:

| group | n | HLA-I-low | B2M-low | either |
|---|---:|---:|---:|---:|
| DM1 | 87 | 18% | 23% | **24%** |
| DM2 | 25 | 84% | 72% | **84%** |

Fisher p (either-low) = **8.7e-8**. DM2 = classical immune-escape; DM1 = presentation-competent.

### (4) Adaptive resistance, not exclusion
DM1 elevated for CTLA4 (d=+1.28), PDCD1 (+0.98), FOXP3 (+1.33), CD274 (+1.52), IFNG (+0.75) — Wherry adaptive-resistance pattern, "recognized but suppressed". Composite **TMB × HLA-I = +1.43 (p=4.6e-7)** despite lower raw TMB: presentation rescues the lower input.

### (5) Hashimoto-allele priming — INDETERMINATE
K2 arcasHLA n=12 DM1 vs n=223 DM2: HLA-C\*14:03 (DM1 21% vs DM2 5.6%, raw p=0.014) and HLA-B\*44:03 (21% vs 7%, p=0.029) trend DM1-enriched, but **BH-FDR = 1.0**. Underpowered. TCGA arcasHLA not run.

## Caveats
- DM2 n=25 small (mitigated: every HLA gene d>1, Wilcoxon ≤1e-7).
- "B2M-low" is mRNA proxy; LOH/mutation-level escape needs MAF re-parse (473 .maf.gz at `/data/thca/data_raw/gdc/TCGA-THCA/mutation/`).
- K2 HLA-allele test underpowered; TCGA arcasHLA proposed as A6000 burst.

## Outputs
`h19_tmb_results.tsv`, `h19_hla_preservation.tsv`, `h19_immune_escape_score.tsv`, `h19_immune_escape_fisher.json`, `h19_hla_alleles_dm1.tsv`, `h19_all_features_combined.tsv`, `h19_braf_cptc_merged.tsv`, `h19_headline.json`.

**For Paper 1+2 BRAF framing**: DM1's heat is preservation + IFN-γ adaptive resistance, not neoantigen excess. Direct Paper-3 ICI-vulnerability link — DM2 BRAF-cPTC = escape-driven, predicted refractory; DM1 BRAF-cPTC = rational ICI candidate.
