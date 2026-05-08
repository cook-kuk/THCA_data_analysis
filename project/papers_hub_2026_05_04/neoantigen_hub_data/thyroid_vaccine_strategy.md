# Thyroid vaccine strategy · permitted vs forbidden claims

Generated: 2026-05-07

## Forbidden claims (DO NOT make)

1. **❌ Validated clinical efficacy** — no published clinical neoantigen vaccine has demonstrated efficacy in thyroid carcinoma as of 2026-05-07.
2. **❌ HLA-binding alone implies vaccine target** — TESLA (Wells 2020 Cell) showed binding ≠ immunogenic. Require T-cell assay (LEVEL_4) at minimum, patient-matched (LEVEL_5) preferred.
3. **❌ Allele frequency = carrier frequency** — DPB1*05:01 frequency ≈ 56% in Korean does NOT mean ~56% of patients have homozygous coverage. Phenotype frequency calculations are required.
4. **❌ Melanoma TMB-high logic transferred to thyroid** — PTC TMB ≈ 0.4 mut/Mb, ATC ≈ 6.1; melanoma ≈ 13.5. Different antigen quantity regime.
5. **❌ "Personalized neoantigen vaccine" without per-patient WES + HLA + RNA + immunogenicity prediction pipeline** — claims must include the full computational chain (and verifier).

## Permitted claims (with evidence ladder)

| Claim level | Required evidence | Example claim |
|---|---|---|
| **Hypothesis generation** | LEVEL_0 + biological rationale | "BRAF V600E is a candidate target for HLA-DQ-restricted CD4 vaccine in BRAFm PTC" |
| **Exploratory prioritization** | LEVEL_2 (MS) + cohort enrichment | "Cancer-testis antigens MAGE-A/MAGE-C1 are upregulated in ATC and warrant immunopeptidomics validation" |
| **Validation roadmap** | LEVEL_4+ (T-cell assay) on the candidate | "We propose ELISpot + tetramer + autologous-organoid cytotoxicity for top 5 candidates" |
| **Companion biomarker design** | LEVEL_5 (patient-matched T-cell) | "Patient HLA-A*11:01 + KRAS G12V detected → lymph-node-targeted amphiphile vaccine candidate per Pant 2024" |
| **Low-TMB feasibility discussion** | TMB measured + neoantigen quality scored | "Even at low TMB, ATC carries 6 mut/Mb plus CT antigens; off-the-shelf shared neoantigens (BRAF, RET) preferred over personalized" |
| **Antigen presentation atlas** | bulk RNA + sc + spatial + HLA typing | "PTC+HT subset has HLA-II d=+3.65 and TLS biology; vaccine-permissive niche identified" |

## Evidence-aligned roadmap

### Step 1 (now)
- HLA-typed thyroid cohort (TCGA-THCA + K2 PRJEB11591 + GSE286332) → arcasHLA + AFND prior + supertype mapping.
- Bulk RNA + WES → driver landscape (BRAF/RAS/TERT/RET/fusions) + Moffitt-equivalent subtype.
- Module scoring (HLA-I/II, IFN-γ, TLS, myeloid, CAF) → readiness map.

### Step 2 (3-6 months)
- Apply BNT122 / ELI-002 / KEYNOTE-942 design templates.
- Identify HLA-A*02:01 / A*11:01 / B*07:02-restricted candidate peptides for top-3 thyroid hotspots.
- ELISpot + tetramer + IFN-γ release on patient PBMCs from the K2 cohort.

### Step 3 (6-18 months)
- Phase I trial design for highest-priority candidate × HLA pairing.
- Companion biomarker: HLA + KRAS/BRAF allele genotype (matches Pant 2024 ELI-002 inclusion criteria).

### Step 4 (18+ months)
- Phase II adjuvant trial.
- Long-term TCR persistence tracking (per Sethna 2025 Nature).

## Comparison to PDAC/RCC/HCC clinical vaccine evidence

- **PDAC** (Rojas 2023 BNT122): high private-neoantigen yield, mFOLFIRINOX + atezo backbone, 8/16 vaccine-responders.
  - Translatable to thyroid: no — TMB too low for personalized pool.
  - Translatable: validation pipeline (TCR persistence tracking).
- **Melanoma** (Weber 2024 KEYNOTE-942): HR 0.561 RFS, mRNA-4157 + pembro.
  - Translatable: ICI-naïve setting + adjuvant timing.
  - Not transferable: TMB-driven antigen quantity rationale.
- **RCC** (Braun 2025 Nature): peptide vaccine + ipi.
  - Translatable: ICI-resistant cold-tumor framework.
  - Cautious: thyroid HLA-II axis differs from RCC.
- **HCC** (Yarchoan 2024 Nat Med): DNA vaccine + checkpoint.
  - Translatable: minimal-residual-disease setting after surgery.

## Bottom line

**Thyroid neoantigen vaccine is at TRL 1-2 (basic principles → technology concept).**
The data hub + benchmark we built supports hypothesis generation and validation roadmap. It does not support clinical efficacy claims. Future progress depends on:

1. Patient-matched T-cell assay data on thyroid mutations (we have ZERO entries in our master at LEVEL_5 from thyroid).
2. Korean institutional cohort (Yonsei drop-in template ready).
3. ICGC / dbGaP controlled-access thyroid datasets.
4. Wet-lab partnership for ELISpot + tetramer + organoid co-culture.
