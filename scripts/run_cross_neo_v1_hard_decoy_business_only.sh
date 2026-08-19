#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/.."

export OMP_NUM_THREADS="${OMP_NUM_THREADS:-4}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-4}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-4}"

python scripts/train_cross_neo_v1_hard_decoy_focal.py
python scripts/evaluate_cross_neo_v1_contract.py
python scripts/write_cross_neo_v1_decision_report.py
python scripts/build_cross_neo_v1_business_package.py

python - <<'PY'
from pathlib import Path
import pandas as pd

out = Path("project/results/cross_neo_v1")
hard = pd.read_csv(out / "hard_decoy_focal_metrics.tsv", sep="\t")
comp = pd.read_csv(out / "v1_model_comparison.tsv", sep="\t")
print("=== hard-decoy/focal summary ===")
print(hard.sort_values(["split_name", "AUPRC"], ascending=[True, False]).to_string(index=False))
hla = comp[(comp["family"] == "v1_hard_decoy_focal") & (comp["split_name"] == "hla_stratified_group_5fold") & (comp["subset"] == "all_rows")]
if len(hla):
    print("")
    print("HLA-stratified hard-decoy rows:")
    print(hla.sort_values("AUPRC", ascending=False)[["model", "AUPRC", "AUROC", "top10_precision", "enrichment_at_10"]].to_string(index=False))
print("")
print(f"business package: {out / 'business_package'}")
PY
