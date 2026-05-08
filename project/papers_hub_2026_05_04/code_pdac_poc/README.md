# PDAC PoC code (Lumenix Cancer Vaccine Agent)

Live run on TCGA-PAAD via cBioPortal API. Run order:

```bash
python 01_fetch_tcga_paad.py      # cBioPortal: clinical + mutations + mRNA z-scores
python 02_pdac_registry.py        # 37-cohort registry (Nature/Cell/Science/GEO/ICGC/CPTAC/EGA/PRIDE)
python 03_moffitt_classifier.py   # Moffitt 2015 basal/classical + KM
python 04_neoantigen_quality.py   # Balachandran 2017 R × D framework
python 05_pdac_immune_readiness.py# 9 PDAC TME modules + Cohen's d
python 06_vaccine_target_priority.py  # 11-target ranking + Korean HLA bias audit
python 07_agent_orchestrator.py "Korean PDAC + KRAS G12V + HLA-A*11:01"
python 08_figures.py              # 5 PNG/SVG figures
python 09_consolidate.py          # single JSON for HTML consumption
```

Live run today on TCGA-PAAD:
- 184 samples · 424 driver mutations · 177×165 mRNA z-score matrix
- 94 basal-like / 83 classical (53.1% basal — within Moffitt 2015 expected 30–70%)
- 110 KRAS-G12 hotspot mutations across 121 samples (92.4% of KRAS)
- Top NeoQ gene = KRAS (PDAC expected, verifier PASS)
- 11-target vaccine panel ranked; top-1 = KRAS G12D (priority 0.237)
- Korean off-the-shelf cassette population coverage = 14.93% (verifier WARN — recommend personalized track)

Architecture: SPARK-inspired (Berbís et al. Nat Med 2026 s41591-026-04357-y).
Components: Planner / ToolRouter / 8 tools (T1–T8) / Verifier / Memory.
