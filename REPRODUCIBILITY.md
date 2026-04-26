# THCA Dashboard — Reproducibility Notes

This document captures the methodological + environment choices that
determine whether another team can exactly reproduce the numbers in the
dashboard.

---

## 1. Environment

| Layer | Version |
|---|---|
| Python | **3.12.3** |
| numpy | 2.4.4 |
| pandas | 3.0.2 |
| scipy | 1.17.1 |
| scikit-learn | 1.8.0 |
| statsmodels | 0.14.6 |
| xgboost | 3.2.0 |
| shap | 0.51.0 |
| torch | 2.11.0+cpu |
| plotly | 6.7.0 |
| qiskit | 2.4.0 |
| pennylane | 0.44.1 |
| umap-learn | 0.5.12 |

All pins are in `requirements.txt`. The CPU-only torch wheel is pulled
from `https://download.pytorch.org/whl/cpu` because the default PyPI
wheel is CUDA-linked.

## 2. Random seeds

| Component | Seed | Location |
|---|---|---|
| `numpy.random.default_rng` (pipeline) | `42` | `run_thca_pipeline.py` |
| `sklearn` splits | `random_state=42` | ML sections |
| `pennylane.devices("default.qubit")` | Deterministic | `ml_quantum_full.py` |
| `qiskit-aer` shots | `shots=1024, seed_simulator=42` | quantum VQC / QAOA |
| Meta-analysis funnel subsample | `random_state=42` | `meta_power.py` |

## 3. Deterministic vs non-deterministic

| Stage | Deterministic? | Notes |
|---|---|---|
| Differential expression (t-test / Cohen's d) | **Yes** | Closed-form stats |
| Meta-analysis (FE, DL random-effect, Fisher) | **Yes** | Closed-form |
| FDR (Benjamini–Hochberg) | **Yes** | Rank-based |
| Power / sample-size calculator | **Yes** | Analytic formula |
| ML training (XGBoost, MLP) | Seeded but **not bit-exact** | Threading & BLAS kernels introduce ±1e-6 drift |
| UMAP / t-SNE projections | Seeded but **not bit-exact** | Numba JIT order can differ across CPUs |
| SHAP TreeExplainer | **Yes** on same model | Sample order in swarm plot is seeded |
| Quantum VQC (statevector) | **Yes** | ideal simulator |
| Quantum VQC (shot-based) | Seeded (`seed_simulator=42`) | Identical shots → identical outcome |
| QAOA | Seeded | Classical optimiser may still hit different local minima at very tight tolerances |

## 4. Data provenance

| Source | Samples | Accession |
|---|---|---|
| TCGA-THCA | 503 RNA-seq (BRAF_like=392, RAS_like=111) | GDC PanCanAtlas |
| GSE27155 | 72 microarrays | GEO |
| GSE126698 | 12 microarrays | GEO |

## 5. Meta-analysis specifics (`meta_power.py`)

- **Fixed-effect (FE)**: inverse-variance weighted mean of log2FC.
- **Random-effect (RE)**: DerSimonian–Laird τ² estimator applied on top of FE.
- **Heterogeneity**: Cochran's Q, I² = max(0, (Q − df)/Q).
- **Standard errors**:
  - TCGA: derived from log2FC + Cohen's d + n_braf/n_ras.
  - GSE27155 / GSE126698: derived from log2FC + reported two-sided p-value via `SE = |log2FC| / z_{p/2}`.
  Both approximations ignore heteroscedasticity between the two BRAF/RAS subgroups and are therefore *approximate* confidence intervals. The point estimate is exact; the CI width should be read as an upper bound.
- **Fisher's combined p**: `-2 Σ log(p_i) ~ χ²_{2k}`.
- **Concordance flag**: `True` iff sign(log2FC) is identical across **all** cohorts where the gene is present.
- **FDR**: Benjamini–Hochberg on RE p across all genes with ≥2 cohorts.

## 6. Power calculator specifics

- Two-sample t-test approximated by z-test (large n): `n = ((z_{α/2} + z_β)² × (σ1² + σ2²)) / d²`.
- σ=1 because d is already standardised (Cohen's d).
- Cohort per-group n approximated as total/2 (TCGA: 333→166, GSE27155: 72→36, GSE126698: 12→6). For TCGA this under-estimates power for the majority class (BRAF_like=392) and slightly over-estimates for the minority class (RAS_like=111) — the harmonic mean is closer to reality; we use total/2 for pedagogical symmetry.
- Achievable power: `Φ(sqrt(n d² / 2) − z_{α/2})`.

## 7. What we do NOT guarantee

- Identical rendered pixel output across Plotly versions > 6.7.0.
- Identical embedding coordinates if you rebuild on a different CPU architecture (UMAP numba output depends on SIMD path).
- Identical floating-point p-values past the 10th decimal place across BLAS backends (OpenBLAS vs MKL).

## 8. Repro check-list

1. `python3.12 -m venv .venv && pip install -r requirements.txt`
2. `project/.venv/bin/python project/notebooks_or_scripts/meta_power.py`
3. Compare `project/results/tables/meta_analysis.tsv` row count to the in-repo snapshot (should be ≥ 25,000 rows, ≥ 10,000 concordant).
4. Compare top-10 3-cohort concordant genes by `meta_q` — the gene names should match regardless of BLAS backend.
