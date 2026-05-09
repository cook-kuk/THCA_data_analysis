#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/../../.."

python project/scripts/cross_neo_v2/11_make_real_esm2_features.py
python project/scripts/cross_neo_v2/04_train_models.py
python project/scripts/cross_neo_v2/05_evaluate_locked.py
python project/scripts/cross_neo_v2/06_source_collapse_rescue.py
python project/scripts/cross_neo_v2/07_case_audit.py
python project/scripts/cross_neo_v2/08_make_figures.py
python project/scripts/cross_neo_v2/09_write_manuscript_packages.py
python project/scripts/cross_neo_v2/10_write_final_report.py
