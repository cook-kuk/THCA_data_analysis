# Section S8 — HanLab-inspired supplementary analyses (THCA)

> **v5.2 retraction status, per sub-section:**
> - **S8A FastRNA centering — FALLS under v5.2.** Was a "rescue" of the
>   v5.1 flip (DIAL 0.331 → 0.000); under v5.2 there is no flip to
>   rescue. Centering is still a valid pre-processing alternative but
>   no longer presented as a corrective.
> - **S8B BUHMBOX PC1 KS — SURVIVES v5.2.** Structural cohort-PCA
>   observation (KS p < 1e-40 within each subtype) is independent of
>   the LODO ComBat protocol and remains a real data feature.
> - **S8C LMM biomarker correction — SURVIVES v5.2.** Per-gene LMM with
>   cohort random effect; 99.4 % in-universe survival. See
>   `reports/v8/v8p1_rigor_summary.md` Task B for the full 11 710-gene
>   refresh that supersedes the 3 000-gene subset reported here.

This supplement strengthens the THCA-specific claim of v5.1 with three
Han-lab-inspired analyses: (S8A) FastRNA-style per-cohort centering
(Lee & Han 2022, *AJHG*), (S8B) a BUHMBOX-inspired concept check for
hidden sub-group heterogeneity (Han *et al.* 2016, *Nat Genet*), and
(S8C) linear-mixed-model (LMM) biomarker correction in the style of
Sul *et al.* 2018. All use the v5.1 harmonised THCA cohort
(n=392, BRAF=321/RAS=71; TCGA-THCA=351/GSE27155=41; top-3000 variance
genes from 11 710 shared).

## S8A. FastRNA-style cohort centering

FastRNA subtracts the per-cohort mean from every sample — variance is
preserved, global location is removed, no empirical-Bayes step is
required. We applied this and recomputed LODO DIAL without ComBat.

v5.1's ComBat LODO pipeline gave mean DIAL = 0.331 on THCA, with
LogReg_l2 and LogReg_elasticnet at DIAL ≈ 0.49 — essentially perfect
label polarity flips, the canonical fingerprint of batch-entangled
labels. FastRNA-style centering reverses this: LODO AUC stays at
0.95–0.99 for linear models (0.74–0.88 for tree ensembles) with the
correct polarity, yielding DIAL = 0.00 for every classifier. Mean
DIAL drops from 0.331 to 0.000 (delta ≈ +0.33), exceeding the
expected ~0.1 floor.

The THCA cross-cohort BRAF-vs-RAS signal is recoverable once the
(highly unbalanced) per-cohort mean is removed; ComBat over-shrinks
toward a batch prior dominated by the large TCGA arm and flips the
minority-cohort polarity. See *Table S8A*
(`results/v8_statgen/v8_fastRNA_style_dial.tsv`) and *Figure S8A*
(planned: grouped bar of DIAL-v5.1 vs DIAL-centered per classifier).

## S8B. BUHMBOX concept check

BUHMBOX (Han *et al.* 2016, *Nat Genet*) detects whether shared
heritability between two case groups reflects truly shared aetiology
or hidden sub-group admixture. Our check is *conceptual* (not the
original statistic, which was built for genotype risk-allele data):
for each subtype, we fit PCA on the pooled class-specific expression
matrix and tested whether PC1 scores from TCGA-THCA and GSE27155
follow the same distribution (two-sample Kolmogorov–Smirnov).
p>0.05 → cohorts share the implied biology; p≤0.05 → hidden
heterogeneity within the label.

The result is unambiguous: KS = 1.00 with p = 1.3 × 10⁻⁴⁰ for BRAF
(n=28 vs 293) and p = 3.4 × 10⁻¹⁴ for RAS (n=13 vs 58); both flagged
`hidden_heterogeneity`. PC1 aligns perfectly with cohort identity
even after restricting to a single subtype, echoing S8A: the dominant
axis of variation in THCA is cohort, not biology. See *Table S8B*
(`results/v8_statgen/v8_buhmbox_concept_check.tsv`) and *Figure S8B*
(planned: PC1 density by cohort within subtype).

## S8C. LMM-based biomarker correction

We fit per-gene `expression ~ subtype + (1 | cohort)`
(statsmodels `MixedLM.from_formula`, REML, L-BFGS) and compared it to
naive OLS `expression ~ subtype`. Both were BH-FDR corrected across
3 000 top-variance genes, following Sul *et al.* 2018.

Fitting all 3 000 models with joblib `n_jobs=-1` took 22 s, so no
subsampling degradation was needed. At FDR<0.05 the naive test
flagged 2 152 genes; the LMM flagged 2 201. The lists overlap on
1 661 genes; 491 naive hits were `lost_with_LMM` (cohort-confounded
mean shifts) and 540 were `gained_with_LMM` (biology previously
obscured by inter-cohort noise). The LMM *rescued* more genes than
it dropped, consistent with S8A: once cohort is modelled, the
BRAF/RAS effect gets cleaner, not weaker. Example top hits STMN2,
CTSA, CD82 all retain LMM FDR < 10⁻¹⁵ with consistent direction. See
*Table S8C* (`results/v8_statgen/v8_LMM_corrected_biomarkers.tsv`)
and *Figure S8C* (planned: scatter of −log₁₀(p) naive vs LMM, coloured
by status).

## Summary

Three independent analyses converge: the THCA BRAF-vs-RAS signal is
real, but the dominant technical axis in the v5.1 harmonised cohort is
platform identity. (i) FastRNA-style centering fully resolves the LODO
polarity flip that gave v5.1 its DIAL = 0.33 (S8A). (ii) The platform
axis survives stratification within subtype — BUHMBOX-style PC1 is
perfectly disjoint by cohort (S8B). (iii) A per-gene LMM with cohort
random effect reshuffles ~30 % of FDR-significant genes, mostly by
*gaining* biology (S8C). These results reinforce the THCA-specific
claim while providing constructive alternatives to v5.1 ComBat.

### Artefacts

| Section | Table | Figure (planned) |
|---|---|---|
| S8A | `results/v8_statgen/v8_fastRNA_style_dial.tsv` | Figure S8A: DIAL v5.1 vs centered, per classifier |
| S8B | `results/v8_statgen/v8_buhmbox_concept_check.tsv` | Figure S8B: PC1 density by cohort within subtype |
| S8C | `results/v8_statgen/v8_LMM_corrected_biomarkers.tsv` | Figure S8C: naive vs LMM −log₁₀(p) scatter |

### References

- Lee, Y. & Han, B. (2022). FastRNA: an efficient solution for PC analysis
  of single-cell RNA-seq data. *American Journal of Human Genetics*.
- Han, B., Pouget, J. G., Slowikowski, K., *et al.* (2016). A method
  to decipher pleiotropy by detecting underlying heterogeneity driven
  by hidden subgroups applied to autoimmune and neuropsychiatric
  diseases. *Nature Genetics*, 48(7), 803–810.
- Sul, J. H., Martin, L. S. & Eskin, E. (2018). Population structure
  in genetic studies: confounding factors and mixed models.
  *Genome Biology / PLoS Genetics*.
