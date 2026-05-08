# Post-QC paper decision - 2026-05-08

## Bottom line

The paper is still viable, but the headline must be narrowed.

Do **not** lead with "CNV-high territories drive TROP2." After QC and tissue-composition residualization, TROP2/TACSTD2 reverses direction, so that claim is not robust.

The defensible high-impact claim is:

> Driver-negative thyroid cancer contains a recurrent copy-number residual genomic class, and spatial transcriptomics reveals coherent CNV-like RNA territories whose most robust downstream ecosystem signal is an immune interface, especially TAM/cytotoxic/TLS programs. H&E/UNI can partly read downstream phenotypes, but not the CNV territory itself.

## What strengthened today

- Gene-order inferCNV-lite finished across **57,997 spots / 16 slides**.
- Cancer-slide median KNN coherence z for inferCNV-lite T00 territories: **57.9**.
- Cancer-slide median bin-bootstrap ARI: **0.859**.
- Gene-bin territory concordance with the prior arm-level territory:
  - pooled cancer ARI **0.706**
  - high-territory Jaccard **0.801**
- This supports the spatial territory figure as a sensitivity layer.

## What weakened today

- Both arm-level and gene-bin expression-CNV high territories are strongly coupled to sequencing depth:
  - arm-level slide-median high-low log-count delta **2.01**
  - arm-level slide-median high-low n-gene delta **2.26e3**
  - gene-bin slide-median high-low log-count delta **2.04**
  - gene-bin slide-median high-low n-gene delta **2.31e3**
- After slide + QC + composition residualization:
  - arm-level TROP2_raw beta **-0.382**, r **-0.0477**
  - arm-level TACSTD2_z beta **-0.360**, r **-0.0875**
  - gene-bin TROP2_raw beta **-0.354**, r **-0.0434**
  - gene-bin TACSTD2_z beta **-0.370**, r **-0.0879**

This means TROP2/TACSTD2 should be moved to exploratory/multimodal phenotype, not causal CNV-territory biology.

## What survives best

After slide + QC + composition residualization:

- arm-level Macrophage_TAM_z beta **0.176**, r **0.0742**
- arm-level Tcell_cytotoxic_z beta **0.190**, r **0.0855**
- arm-level TLS_B_z beta **0.0400**, r **0.0205**
- gene-bin Macrophage_TAM_z beta **0.203**, r **0.0839**
- gene-bin Tcell_cytotoxic_z beta **0.208**, r **0.0915**
- gene-bin TLS_B_z beta **0.0498**, r **0.0251**

These are modest but consistent. The immune-interface claim is now the safest spatial biology layer.

## Pathology/multimodal role

Pathology should not be framed as a CNV replacement:

- UNI image-only to ST-CNV high territory AUROC **0.534**.
- UNI image-only to inferCNV-lite T00-high AUROC **0.506**.
- UNI image-only to continuous inferCNV-lite aneuploidy score rho **0.465** and bin SD rho **0.740**, but binary high-territory prediction remains poor.

Use pathology as a boundary figure: H&E/UNI can read tissue phenotype and some downstream ecosystem signal, but not robustly recover molecular CNV territory.

## Revised figure spine

1. **Bulk genomic anchor**: TCGA/cBio Class6-DM2 CNV residual class, DM2 35.4% vs DM1 1.3%, OR 41.1, p=1.6e-7.
2. **Robustness**: histology, purity, stage, BRAF/RAS strata.
3. **Spatial territory map**: arm-level expression-CNV territories plus gene-order inferCNV-lite sensitivity.
4. **QC-aware ecosystem layer**: immune/TAM/cytotoxic/TLS enrichment survives residualization; TROP2 is QC-sensitive.
5. **Multimodal boundary**: H&E/UNI reads ecosystem/phenotype better than CNV territory.
6. **Validation ask**: allele-specific CalicoST/infercnvpy and IF/IHC for TAM/T cell/TLS rather than TROP2-first.

## Decision

Proceed, but with the title softened:

> Copy-number residual thyroid cancers exhibit spatially organized immune ecosystems

Avoid:

> Copy-number territories drive TROP2 niches

The former is defensible with the current data. The latter is likely to fail reviewer scrutiny.

## One-hour reviewer stress update

I added a reviewer-style stress battery after the first QC memo. This used slide-level residual effects, condition-stratified residual effects, leave-one-slide influence, spatial block permutation, and QC/composition predictability of the territory labels.

New stress-test result:

- QC/composition predicts territory labels very strongly within slide:
  - arm-level territory median apparent AUROC **0.970**
  - gene-bin inferCNV-lite territory median apparent AUROC **0.985**
- Leave-one-slide-out QC/composition predictability is weaker but still non-trivial:
  - arm-level AUROC **0.601**
  - gene-bin AUROC **0.621**

Reviewer-safe scorecard:

| Outcome | Arm-level territory | Gene-bin inferCNV-lite territory | Decision |
|---|---|---|---|
| Macrophage_TAM_z | survives_strong | survives_strong | Main spatial ecosystem claim |
| Tcell_cytotoxic_z | survives_moderate | survives_moderate | Secondary immune-interface claim |
| TLS_B_z | survives_moderate | survives_moderate | Secondary immune-interface claim |
| TROP2_raw | fails_or_reverses | directionally_inconsistent | Exploratory only |
| TACSTD2_z | directionally_inconsistent | directionally_inconsistent | Exploratory only |

Key numbers:

- Macrophage_TAM_z:
  - arm-level median slide beta **0.107**, positive slides **10/12**, conditions **3/3**, spatial block p **0.00249**
  - gene-bin median slide beta **0.0916**, positive slides **10/12**, conditions **3/3**, spatial block p **0.00249**
- Tcell_cytotoxic_z:
  - arm-level median slide beta **0.0469**, positive slides **9/12**, conditions **2/3**, spatial block p **0.00249**
  - gene-bin median slide beta **0.0643**, positive slides **9/12**, conditions **2/3**, spatial block p **0.00249**
- TLS_B_z:
  - arm-level median slide beta **0.0394**, positive slides **8/12**, conditions **2/3**, spatial block p **0.0698**
  - gene-bin median slide beta **0.0800**, positive slides **10/12**, conditions **2/3**, spatial block p **0.0299**

Updated final title direction:

> Copy-number residual thyroid cancers exhibit spatially organized myeloid and lymphoid immune ecosystems

Updated claim hierarchy:

1. Bulk genomic claim: true driver-negative THCA contains a recurrent CNV residual class.
2. Spatial method claim: arm-level and gene-order expression-CNV territories are coherent, but QC-sensitive.
3. Biology claim: the robust downstream spatial signal is TAM/cytotoxic/TLS immune ecology, not TROP2.
4. Multimodal claim: H&E/UNI reads downstream phenotype better than CNV territory; use as boundary, not replacement.

This is now tighter and more reviewer-resistant than the earlier TROP2-forward version.
