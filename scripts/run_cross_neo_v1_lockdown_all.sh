#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/.."

export OMP_NUM_THREADS="${OMP_NUM_THREADS:-8}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-8}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-8}"

python scripts/lock_cross_neo_v1_anchor.py
python scripts/test_cross_neo_v1_foldsafe_fusion.py
python scripts/analyze_cross_neo_v1_qk_rescue_harm.py
python scripts/diagnose_cross_neo_v1_source_collapse.py
python scripts/test_cross_neo_v1_source_topk_rescue.py
python scripts/build_cross_neo_v1_overlap_manifest.py
python scripts/evaluate_cross_neo_v1_lockdown.py
python scripts/make_cross_neo_v1_lockdown_figures.py
python scripts/write_cross_neo_v1_lockdown_report.py
python scripts/build_cross_neo_v1_reviewer_audit_pack.py

python - <<'PY'
from pathlib import Path
import pandas as pd

out = Path("project/results/cross_neo_v1_lockdown")
report = out / "CROSS_Neo_v1_lockdown_decision_report.md"
comparison = pd.read_csv(out / "v1_lockdown_comparison.tsv", sep="\t")
print("")
print("=== CROSS-Neo v1 lockdown final ===")
print(f"decision report: {report}")
print(f"reviewer audit pack: {out / 'CROSS_Neo_v1_reviewer_audit_pack.md'}")
print("")
print("top 20 rows of v1_lockdown_comparison.tsv:")
cols = [
    "split_name",
    "method",
    "method_family",
    "n",
    "n_pos",
    "prevalence",
    "AUPRC",
    "AUROC",
    "top5_precision",
    "top10_precision",
    "top20_precision",
    "enrichment_at_10",
]
print(comparison.sort_values(["AUPRC", "top10_precision"], ascending=False)[cols].head(20).to_string(index=False))
print("")
print("source_collapse_report.md summary:")
txt = (out / "source_collapse_report.md").read_text().splitlines()
for line in txt[:45]:
    print(line)
print("")
print("foldsafe_fusion_report.md summary:")
txt = (out / "foldsafe_fusion_report.md").read_text().splitlines()
for line in txt[:45]:
    print(line)
print("")
print("QK rescue/harm summary:")
qk = pd.read_csv(out / "qk_rescue_harm_summary.tsv", sep="\t")
print(qk.sort_values(["split_name", "comparator", "event"]).head(30).to_string(index=False))
decision = report.read_text().split("Decision: **", 1)[1].split("**", 1)[0]
print("")
print(f"final KEEP/HOLD/KILL: {decision}")
PY
