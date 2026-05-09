#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/../../.."

python project/scripts/cross_neo_v2/00_build_registry.py
python project/scripts/cross_neo_v2/01_overlap_audit.py
python project/scripts/cross_neo_v2/02_make_strict_splits.py
python project/scripts/cross_neo_v2/03_make_features.py
python project/scripts/cross_neo_v2/04_train_models.py
python project/scripts/cross_neo_v2/05_evaluate_locked.py
python project/scripts/cross_neo_v2/06_source_collapse_rescue.py
python project/scripts/cross_neo_v2/07_case_audit.py
python project/scripts/cross_neo_v2/08_make_figures.py
python project/scripts/cross_neo_v2/09_write_manuscript_packages.py
python project/scripts/cross_neo_v2/10_write_final_report.py

python - <<'PY'
from pathlib import Path
import pandas as pd

out = Path("project/results/cross_neo_v2_sota_sprint_2026_05_09")
print("\n=== CROSS-Neo v2 outputs ===")
print(out / "CROSS_Neo_v2_SOTA_sprint_decision_report.md")
print(out / "CROSS_Neo_v2_all_results_summary.xlsx")
cmp = out / "metrics/all_model_all_split_metrics.tsv"
if cmp.exists():
    df = pd.read_csv(cmp, sep="\t")
    cols = ["split_name", "model_name", "n", "n_pos", "prevalence", "AUPRC", "AUROC", "top10_precision", "top20_precision", "claim_status"]
    print("\n=== top 20 all-model rows by split/AUPRC ===")
    print(df.sort_values(["split_name", "AUPRC", "top10_precision"], ascending=[True, False, False])[cols].head(20).to_string(index=False))
rep = out / "CROSS_Neo_v2_SOTA_sprint_decision_report.md"
if rep.exists():
    txt = rep.read_text().splitlines()
    print("\n=== decision report head ===")
    print("\n".join(txt[:60]))
PY
