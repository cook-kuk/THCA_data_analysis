#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/.."

export OMP_NUM_THREADS="${OMP_NUM_THREADS:-8}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-8}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-8}"

mkdir -p project/results/cross_neo_v1/figures

python scripts/analyze_cross_neo_v0_branch_errors.py
python scripts/diagnose_cross_all_failure.py
python scripts/train_cross_neo_v1_gated_moe.py
python scripts/train_cross_neo_v1_source_bias_corrected.py
python scripts/train_cross_neo_v1_pu_ranking.py
python scripts/train_cross_neo_v1_decoy_focal.py
python scripts/build_cross_neo_v1_public_overlap_audit.py
python scripts/evaluate_cross_neo_v1_contract.py
python scripts/make_cross_neo_v1_figures.py
python scripts/write_cross_neo_v1_decision_report.py

python - <<'PY'
from pathlib import Path
import pandas as pd

out = Path("project/results/cross_neo_v1")
comp = pd.read_csv(out / "v1_model_comparison.tsv", sep="\t")
v1 = comp[(comp["family"].str.startswith("v1")) & (comp["subset"] == "all_rows")].copy()
hla = v1[v1["split_name"] == "hla_stratified_group_5fold"].sort_values("AUPRC", ascending=False)
top10 = v1.sort_values("top10_precision", ascending=False)
clean = comp[(comp["subset"] == "retrieval_clean_only") & (comp["family"].str.startswith("v1"))].sort_values("AUPRC", ascending=False)
source = v1[v1["split_name"].str.startswith("source_heldout", na=False)].sort_values("AUPRC", ascending=False)
fixed = comp[(comp["model"] == "v0_fixed_late_fusion_C_QK_no_anchor") & (comp["split_name"] == "hla_stratified_group_5fold") & (comp["subset"] == "all_rows")]
best = hla.iloc[0] if len(hla) else None
beats_fixed = bool(best is not None and len(fixed) and best["AUPRC"] > fixed.iloc[0]["AUPRC"])
decision = (out / "CROSS_Neo_v1_decision_report.md").read_text().split("Executive Verdict: **", 1)[1].split("**", 1)[0]
manifest = pd.read_csv(out / "public_overlap_download_manifest_needed.tsv", sep="\t") if (out / "public_overlap_download_manifest_needed.tsv").exists() else pd.DataFrame()

print("")
print("=== CROSS-Neo v1 summary ===")
if best is not None:
    print(f"best v1 model by HLA-stratified AUPRC: {best['family']} / {best['model']} AUPRC={best['AUPRC']:.3f} AUROC={best['AUROC']:.3f} top10={best['top10_precision']:.3f}")
if len(top10):
    r = top10.iloc[0]
    print(f"best v1 model by top10 precision: {r['family']} / {r['model']} split={r['split_name']} top10={r['top10_precision']:.3f} AUPRC={r['AUPRC']:.3f}")
if len(clean):
    r = clean.iloc[0]
    print(f"best retrieval-clean-only model: {r['family']} / {r['model']} split={r['split_name']} AUPRC={r['AUPRC']:.3f} top10={r['top10_precision']:.3f}")
if len(source):
    r = source.iloc[0]
    print(f"best source-heldout model: {r['family']} / {r['model']} split={r['split_name']} AUPRC={r['AUPRC']:.3f} top10={r['top10_precision']:.3f}")
print(f"v1 beats v0 fixed late fusion on HLA-stratified AUPRC: {beats_fixed}")
print(f"final KEEP/HOLD/KILL: {decision}")
if len(manifest):
    print("next missing data/downloads:")
    print(manifest[["dataset", "status", "local_target_path"]].to_string(index=False))
else:
    print("next missing data/downloads: none recorded")
PY
