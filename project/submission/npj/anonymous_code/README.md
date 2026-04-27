# Anonymous code archive — npj submission supplementary

This archive contains all pipeline scripts for the manuscript.
Author identifiers have been removed for double-blind review.

## Layout
- `code/v17p35_*.py` — main analysis pipeline (Phase A → SYNTH-1).
- `code/v17_FINAL_*.py` — final ship sprint (manuscript v2 + Fig8 + audit).

## Replication
1. `python -m venv .venv && source .venv/bin/activate`
2. `pip install plotly pandas numpy scipy scikit-learn 'kaleido==0.2.1'`
3. Run `code/v17p35_SYNTH1_Fig1.py` … `code/v17_FINAL_F1B_fig8.py` to regenerate all figures.
4. Outputs deterministic with `random_state=42` / `np.random.seed(42)`.

## License
MIT.

## Data
Raw TCGA + GEO data described in Methods §4. cBioPortal `thca_tcga_pub` mirror used for TERT recovery.
