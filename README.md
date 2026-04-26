# THCA Multi-Omics Dashboard — Code Archive

Self-contained archive of the THCA (papillary thyroid carcinoma) multi-omics analysis
pipeline and its interactive dashboard. Built on top of TCGA-THCA + 4 external GEO
cohorts (GSE27155, GSE76039, GSE97466, GSE126698, GSE213647).

## What's inside

| Path | Purpose |
|---|---|
| `project/notebooks_or_scripts/` | Main pipeline scripts (preprocessing, ML, quantum, biomarker, pathway) |
| `project/scripts/` | Dashboard build, serve, injectors, UX polish |
| `project/metadata/` | Sample master, gene panels (TDS16, TierA67, BRS71), dataset manifest |
| `project/qc/` | Per-modality QC metrics |
| `project/reports/html/` | The deployable static dashboard (22 pages + assets) |
| `project/results/` | ML, biomarker, drug-discovery, pathway result tables and figures |
| `project/state/` | Runtime PID files |
| `project/thyroid-http.service` | Hardened systemd unit for the serve_secure.py wrapper |
| `node_modules/` | Playwright (for `ui_bug_hunt.js`) |

## What's excluded (can be regenerated)

- `data_raw/` — 6.1 GB, raw TCGA STAR counts + MAFs + GEO SOFT/family matrices. Re-download via `run_thca_pipeline.py`.
- `data_processed/` — 2.1 GB, log2-transformed matrices. Rebuild via `run_thca_pipeline.py` → `rerun_v2.py`.
- `.venv/` — Python 3.12 virtual environment. Recreate with the pip list below.
- `logs/ui_screenshots/` — screenshot artifacts.

## Requirements

Python 3.12 venv with:

```
pandas numpy scipy scikit-learn matplotlib plotly statsmodels shap requests
torch --index-url https://download.pytorch.org/whl/cpu
pennylane pennylane-lightning qiskit qiskit-machine-learning qiskit-aer
qiskit-optimization qiskit-algorithms dimod dwave-neal xgboost
```

Node 20+ with `playwright` (installed under `node_modules/`).

## Run the dashboard

```bash
sudo cp project/thyroid-http.service /etc/systemd/system/thyroid-http.service
sudo systemctl daemon-reload
sudo systemctl enable --now thyroid-http
# Dashboard: http://<host>:8012/reports/html/index.html
```

Or one-off:

```bash
python3 project/scripts/serve_secure.py
```

## Rebuild pipeline (requires raw data)

```bash
cd project
.venv/bin/python notebooks_or_scripts/run_thca_pipeline.py
.venv/bin/python notebooks_or_scripts/rerun_v2.py
.venv/bin/python notebooks_or_scripts/biomarker_analysis.py
.venv/bin/python notebooks_or_scripts/ml_dl_quantum.py
.venv/bin/python notebooks_or_scripts/ml_quantum_full.py
.venv/bin/python notebooks_or_scripts/pathway_immune_meth.py
.venv/bin/python scripts/build_html_reports.py
```

## Key findings (summarized in the dashboard)

- Classical LogReg on `TierA67_clean` (54 genes, no driver-gene leakage) reaches AUC = 1.00 on TCGA internal CV and 0.97–0.89 on external GEO cohorts.
- Deep Learning MLP ties internal, overfits externally.
- Among 12 quantum / quantum-inspired algorithms, none beat classical. QUBO-based feature selection (QAOA or dwave-neal) picks a 14-gene sub-panel that slightly improves external AUC (0.97 → 0.98, 0.89 → 0.94).
- 2,773 novel-validated biomarkers, top-5: TACSTD2 · PLEKHA6 · CYP1B1 · TMPRSS4 · LDLR — all replicated in both external cohorts.
- Biomarker-to-drug pipeline: 8 targets → 15 compounds via ChEMBL, 1 validated (CYP1B1) + 2 emerging + 5 novel.
- BRAF-like: Proteases + Inflammation + Cell-cycle + MAPK. RAS-like retains TDS differentiation signal.

## License / caveats

Computational triage only; not medical or investment advice. Label provenance for
GSE27155 and GSE126698 is histology-proxy, not mutation-verified. See
`reports/html/pages/14_caveats.html`.
