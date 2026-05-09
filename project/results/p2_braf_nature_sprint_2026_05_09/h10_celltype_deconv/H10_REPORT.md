# H10 — Cell-type deconvolution within BRAF-cPTC (DM1 vs DM2)

**Bottom line.** DM1 is **NOT a B-cell-only phenotype**. Within BRAF-cPTC the
DM1↑ compartment is a **coordinated multi-lineage immune infiltrate** — M2
macrophage, DC, Treg, CD4-Tfh, CD8 effector, naïve-B, plasma, M1 — all up
with d≈+0.6 to +1.8. Strongest: **M2 macrophage (d=+1.83)** and **DC
(d=+1.74)**, not B cells. Memory-B is flat (d=+0.30, NS). The HT-13 panel
works because it samples this coordinated program (HLA-II + B-axis +
IFNG + chemokines), not as a B-cell detector.

## Cohort & methods

- **Master**: n=110 (DM1=85, DM2=25), `molecular_subtype="BRAF_like"
  ∧ histology_subtype="cPTC" ∧ dm∈{DM1,DM2}` (`r9_1_master_with_rai.tsv`).
- **CLAM cut**: n=41 (DM1=26, DM2=15) — same stratum ∩ Phase-2 CLAM ids.
  Replicates with equal/stronger d (M2 +2.32, DC +1.82, Treg +1.70, TLS +1.09).

Three layers on z-scored TCGA-THCA bulk RNA: (1) nu-SVR @ 8-class reusing
`p_deconv_2026_05_08/fractions_nu_SVR.tsv` (Lu 2023 pseudobulk; only v2
method that resolves T cell out of NNLS sparse-collapse); (2) LM22-grade
mean-z gene-set scores per subtype (markers in `run_h10.py`); (3) Cabrita
2020 12-gene TLS (12/12 present). Cohen's d (DM1−DM2), MW p, BH-FDR.

## Top hits (BRAF-cPTC n=110, all FDR < 1e-4)

| Compartment           | Layer    | Cohen's d | BH-FDR  |
|-----------------------|----------|-----------|---------|
| Macrophage M2         | gene-set | **+1.83** | 8.5e-10 |
| Dendritic cells       | gene-set | **+1.74** | 1.1e-10 |
| Epithelial (purity↓)  | nu-SVR   | −1.49     | 8.3e-08 |
| Treg (FOXP3+)         | gene-set | **+1.43** | 9.0e-10 |
| Naïve B               | gene-set | +1.26     | 1.2e-09 |
| Malignant (purity↑)   | nu-SVR   | +1.22     | 8.3e-08 |
| CD4 Tfh (CXCL13+)     | gene-set | +1.04     | 2.9e-07 |
| CD8 T (effector)      | gene-set | +1.02     | 1.0e-05 |
| Myeloid (broad)       | nu-SVR   | +1.00     | 7.4e-06 |
| TLS Cabrita-12        | TLS sig  | +0.96     | 6.4e-07 |
| Macrophage M1         | gene-set | +0.95     | 1.0e-04 |

Plasma d=+0.57 (FDR 0.014); Memory-B d=+0.30 NS; broad B-cell (nu-SVR)
d=+0.54 (FDR 0.07). nu-SVR T-cell d=−0.90 reflects fractional dilution by
myeloid expansion, not absolute loss — gene-set CD8 d=+1.02 confirms
absolute increase.

## Interpretation

HT-13 signal is dominated by **myeloid + DC + Treg + Tfh**; B-cell genes
(CD79A/B, MS4A1) ride the Naïve-B / TLS axis, not memory-B expansion.
Immunophenotype = antigen-presenting, TLS-organising, Treg-rich
Hashimoto-overlap PTC — coordinated immune niche, not B-cell-only.

**Files.** `h10_celltype_d.tsv`, `h10_per_sample_scores.tsv`.
