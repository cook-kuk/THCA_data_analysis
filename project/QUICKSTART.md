# THCA Dashboard — Quickstart

End-to-end walkthrough to reproduce the THCA multi-omics dashboard on a
fresh machine. Two paths: **(A) native Python 3.12 venv** or
**(B) Docker**.

---

## A. Native (Python 3.12 venv)

### 1. Clone
```bash
git clone <your-fork-or-mirror>.git THCA_data_analysis
cd THCA_data_analysis
```

### 2. Create venv and install pinned deps
```bash
python3.12 -m venv project/.venv
source project/.venv/bin/activate
pip install --upgrade pip
pip install --index-url https://download.pytorch.org/whl/cpu torch==2.11.0+cpu
pip install -r requirements.txt
```

> The CPU-only torch wheel is ~200 MB; the total install is ~3 GB.

### 3. (Optional) Drop raw data into `project/data_raw/`

The dashboard ships with all *processed* tables under `project/results/`,
so you can skip this step if you only want to view the dashboard.
If you want to **rerun the full pipeline** you'll need:

- `tcga_thca.tar.gz` from GDC (PanCanAtlas THCA)
- `GSE27155_series_matrix.txt.gz` from GEO
- `GSE126698_series_matrix.txt.gz` from GEO

Place them under `project/data_raw/` matching the layout documented in
`project/notebooks_or_scripts/run_thca_pipeline.py`.

### 4. Run the pipeline (optional — outputs are already in `project/results/`)
```bash
cd project
.venv/bin/python notebooks_or_scripts/run_thca_pipeline.py
```

### 5. Run meta-analysis + generate power/forest figures
```bash
.venv/bin/python notebooks_or_scripts/meta_power.py
```

### 6. Start the secure server
```bash
.venv/bin/python scripts/serve_secure.py
# listens on 0.0.0.0:8012 by default (override with PORT env var)
```

### 7. Open in browser
- Home: <http://localhost:8012/reports/html/index.html>
- Meta-analysis / Power / Reproducibility: <http://localhost:8012/reports/html/pages/20_meta_power.html>

---

## B. Docker

### 1. Build
```bash
docker build -t thca-dashboard .
```

### 2. Run
```bash
docker run --rm -p 8012:8012 thca-dashboard
# or with user-supplied raw/processed data:
docker compose up -d
```

### 3. Browse
Same URLs as path A.

---

## Troubleshooting

| Symptom | Fix |
|---|---|
| `torch` wheel not found | install separately first: `pip install --index-url https://download.pytorch.org/whl/cpu torch==2.11.0+cpu` |
| `qiskit` import errors | ensure `qiskit-aer==0.17.2` matches `qiskit==2.4.0` — both are pinned in `requirements.txt` |
| Plotly figures blank | your browser is blocking the vendored `plotly-2.35.2.min.js`; check devtools network tab |
| Port 8012 in use | `PORT=9000 python scripts/serve_secure.py` |

---

## File layout reference
```
THCA_data_analysis/
├── Dockerfile
├── docker-compose.yml
├── requirements.txt
├── QUICKSTART.md              ← you are here
├── REPRODUCIBILITY.md
└── project/
    ├── data_raw/              ← user supplies
    ├── data_processed/        ← auto-generated
    ├── metadata/
    ├── notebooks_or_scripts/
    │   └── meta_power.py      ← Phase 1 addition
    ├── scripts/
    │   └── serve_secure.py    ← Phase 3 entrypoint
    ├── results/
    │   └── tables/
    │       └── meta_analysis.tsv    ← Phase 1 output
    └── reports/html/
        ├── pages/20_meta_power.html ← Phase 4 page
        ├── figs_interactive/
        │   ├── meta_forest_top20.html
        │   └── meta_funnel.html
        └── assets/
            ├── css/power-calculator.css
            └── js/power-calculator.js
```
