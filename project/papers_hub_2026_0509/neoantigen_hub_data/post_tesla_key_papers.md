# Post-TESLA key papers · neoantigen immunogenicity field map

Generated: 2026-05-07

## The trajectory: binding → presentation → recognition → clinical response

### 2017 — Binding-only era foundation
- **Łuksza et al. 2017 Nature** — neoantigen recognition fitness model (foreignness × dissimilarity-to-self)
- **Balachandran et al. 2017 Nature** — long-survivor PDAC neoantigen quality framework (R × D)
- **Sahin et al. 2017 Nature** — first personalized melanoma vaccine (peptide IVAC)

### 2020 — TESLA: the canonical benchmark
- **Wells et al. 2020 Cell** — TESLA consortium · 608 pMHCs tested · 37 patient-matched T-cell-bound · independent cohort 310 epitopes
- Demonstrated: HLA binding alone fails; presentation + recognition + tumor context all matter

### 2020–2023 — Predictor revolution
- **Reynisson et al. 2020 NAR** — NetMHCpan-4.1 (binding + presentation)
- **O'Donnell et al. 2020 Cell Systems** — MHCflurry-2.0
- **Schmidt et al. 2021 Cell Systems** — PRIME (binding + immunogenicity)
- **Sarkizova et al. 2020 Nat Biotech** — HLAthena (95k MS-validated peptides)
- **Albert et al. 2023 Bioinformatics** — BigMHC (presentation + transfer learning to immunogenicity, MANAFEST)

### 2023 — Harmonization era
- **Müller et al. 2023 Immunity** — NeoRanking · harmonized WES/RNA-seq benchmark · ML methods comparison
- **Borch et al. 2024 Frontiers in Immunology** — IMPROVE · 17,500 tested epitopes · 467 T-cell-recognized · 70 patients

### 2023–2025 — Clinical vaccine breakthroughs
- **Rojas et al. 2023 Nature** — BNT122 (autogene cevumeran) personalized mRNA in PDAC · 8/16 responders · RFS HR ≈ 0.08
- **Sethna et al. 2025 Nature** — BNT122 long-FU · vaccine-induced TCRs persist ≥3 yr
- **Weber et al. 2024 Lancet** — KEYNOTE-942 mRNA-4157/V940 + pembro adjuvant melanoma · RFS HR 0.561
- **Pant et al. 2024 Nat Med** — ELI-002 mKRAS amphiphile in PDAC + CRC · off-the-shelf cassette anchor
- **Braun et al. 2025 Nature** — RCC neoantigen vaccine
- **Yarchoan et al. 2024 Nat Med** — HCC DNA neoantigen vaccine

### 2024–2026 — Database era
- **CEDAR** (Koşaloğlu-Yalçın et al. 2023 NAR) — curated cancer epitope DB (5000+ records ingested in our hub)
- **NEPdb** — positive + negative validated neoepitopes
- **dbPepNeo / dbPepNeo2.0** — LC/MC/HC confidence levels
- **TSNAdb v2.0** — SNV/INDEL/fusion-derived
- **Neodb** — Val-Neo + Driver-Neo subsets
- **NeoTCR** — neoantigen-reactive TCR catalog

### 2026 — TESLA-aware benchmarking
- **Berbís et al. 2026 Nat Med (SPARK)** — agentic AI for cancer pathology — architecture template our hub adopts.

## Field shift summary

> "From binding-only prediction (NetMHCpan era, 2018–2020) to presentation + TCR recognition + tumor context + clinical vaccine response (TESLA + BigMHC + PRIME + Müller + clinical vaccine triplet, 2020–2025)."

## Thyroid cancer relevance

- **No clinical neoantigen vaccine evidence in thyroid cancer to date.**
- **Low-TMB context**: thyroid carcinoma TMB ≈ 0.4–6.1 mut/Mb (PTC ≈ 0.4, ATC ≈ 6.1). The PDAC/melanoma high-TMB neoantigen logic does NOT directly transfer.
- **Available targets**:
  - BRAF V600E (PTC ~50%, ATC ~30%) — HLA-DQ-restricted CD4 evidence (Veatch 2018 JCI); HLA-A*02:01 class-I MS validation NOT robust.
  - RET fusions (PTC); NTRK fusions; ALK fusions — fusion junction neoantigens.
  - TERT promoter mutations — NOT classical neoantigens (non-coding).
  - Cancer-testis antigens (MAGE-A 62%, MAGE-C1 57%, NY-ESO-1 14% in ATC).
- **Strategy**: treat thyroid as exploratory; use PDAC G12D / KRAS G12 cassette designs (ELI-002) and BNT122 framework as design templates, not as validated transferable interventions.

## Citations

All papers above are anchored in our `data_sources.yaml` with PMID + DOI where available.
