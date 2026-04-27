# v15 — DIAL: Direction-Invariant AUC Leakage (anonymous review release)

This is the supplementary code release for v15 NeurIPS 2026 submission
"DIAL: A post-hoc diagnostic for subspace-aligned conditional shift
under linear batch correction". Author identity withheld per
double-blind review.

## Contents

```
scripts/
  v15_theorem2_verify.py    §5.1 numerical validation of Theorem 2
  v15_stress_test.py        §5.2 4×5×5×20 synthetic shift grid
  v15_tta_benchmark.py      §5.3 ten-method TTA benchmark
  v15_scaling.py            §5.4 eight-encoder scaling sweep
  v15_cross_domain.py       §5.5 four-domain reproducibility
checkpoints/
  task[1-5]_*.json          headline metric + seed list per experiment
results/
  theory_validation/        §5.1 long-form TSV (630 rows)
  synthetic_stress/         §5.2 grid + summary TSVs
  tta_benchmark/            §5.3 method comparison TSVs
  foundation_model_scaling/ §5.4 scaling-law TSV
  cross_domain_validation/  §5.5 4-domain TSV
```

## Reproduce

```bash
python -m venv .venv && source .venv/bin/activate
pip install numpy pandas scikit-learn matplotlib plotly
for s in scripts/v15_*.py; do python $s; done
```

Wall time on a single 8-core CPU, 32 GB RAM, no GPU: ~13 minutes
total across all five experiments. Checkpoints reproduce the headline
numbers in the paper (Theorem 2 R²=0.94, ROC=0.78, etc.).

## License

MIT (see LICENSE).
