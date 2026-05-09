#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/.."

export OMP_NUM_THREADS="${OMP_NUM_THREADS:-8}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-8}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-8}"

python scripts/build_cross_neo_v0_master_table.py
python scripts/build_cross_neo_v0_train_only_retrieval.py
python scripts/build_cross_neo_v0_public_overlap_audit.py
python scripts/build_cross_neo_v0_counterfactual_embeddings.py
python scripts/build_cross_neo_v0_structure_geometry.py
python scripts/build_cross_neo_v0_quantum_fallback.py
python scripts/train_cross_neo_v0.py
python scripts/build_cross_neo_v0_late_fusion.py
python scripts/evaluate_cross_neo_v0_source_heldout.py
python scripts/calibrate_cross_neo_v0_ood.py
python scripts/evaluate_cross_neo_v0_contract.py
python scripts/make_cross_neo_v0_figures.py
python scripts/write_cross_neo_v0_decision_report.py

python - <<'PY'
from pathlib import Path
import pandas as pd
out = Path("project/results/cross_neo_v0")
m = pd.read_csv(out / "metrics_by_split.tsv", sep="\t")
late = pd.read_csv(out / "late_fusion_oof_metrics.tsv", sep="\t")
source = pd.read_csv(out / "source_heldout_metrics.tsv", sep="\t")
q = pd.read_csv(out / "qk_fallback_metrics.tsv", sep="\t")
eligible = m[m["n"] >= 50].copy()
if eligible.empty:
    eligible = m.copy()
raw_best_auprc = m.sort_values("AUPRC", ascending=False).iloc[0]
best_auprc = eligible.sort_values("AUPRC", ascending=False).iloc[0]
best_late = late[late["status"].eq("prespecified_equal_weight")].sort_values("AUPRC", ascending=False).iloc[0]
best_top10 = eligible.sort_values(["top10_precision", "AUPRC"], ascending=False).iloc[0]
internal = m[m["split_name"] == "repeated_stratified_5x5_internal"].sort_values("AUPRC", ascending=False).iloc[0]
strict = m[m["split_name"] == "near_peptide_cluster_holdout"].sort_values("AUPRC", ascending=False)
strict = strict.iloc[0] if len(strict) else internal
retr = pd.read_csv(out / "retrieval_leakage_audit.tsv", sep="\t")
hla = eligible[eligible["split_name"].isin(["hla_stratified_group_5fold", "hla_supertype_heldout"])].sort_values("AUPRC", ascending=False)
hla_row = hla.iloc[0] if len(hla) else None
decision = (out / "CROSS_Neo_v0_decision_report.md").read_text().split("Decision: **", 1)[1].split("**", 1)[0]
print("")
print("=== CROSS-Neo v0 summary ===")
print(f"raw best model by AUPRC: {raw_best_auprc['feature_group']} / {raw_best_auprc['model']} / {raw_best_auprc['split_name']} n={int(raw_best_auprc['n'])} AUPRC={raw_best_auprc['AUPRC']:.3f} AUROC={raw_best_auprc['AUROC']:.3f}")
print(f"eligible best model by AUPRC (n>=50): {best_auprc['feature_group']} / {best_auprc['model']} / {best_auprc['split_name']} AUPRC={best_auprc['AUPRC']:.3f} AUROC={best_auprc['AUROC']:.3f}")
print(f"best prespecified late fusion: {best_late['fusion']} / {best_late['split_name']} AUPRC={best_late['AUPRC']:.3f} AUROC={best_late['AUROC']:.3f} top10={best_late['top10_precision']:.3f}")
print(f"best model by top10 precision: {best_top10['feature_group']} / {best_top10['model']} / {best_top10['split_name']} top10={best_top10['top10_precision']:.3f}")
print(f"internal repeated result: {internal['feature_group']} / {internal['model']} AUPRC={internal['AUPRC']:.3f} AUROC={internal['AUROC']:.3f} top10={internal['top10_precision']:.3f}")
print(f"strict no-reference result: {strict['feature_group']} / {strict['model']} / {strict['split_name']} AUPRC={strict['AUPRC']:.3f} AUROC={strict['AUROC']:.3f} top10={strict['top10_precision']:.3f}")
print("retrieval-contaminated result counts:")
print(retr.to_string(index=False))
if hla_row is not None:
    print(f"HLA-heldout result: {hla_row['feature_group']} / {hla_row['model']} split={hla_row['split_name']} AUPRC={hla_row['AUPRC']:.3f} AUROC={hla_row['AUROC']:.3f}")
else:
    print("HLA-heldout result: not available")
src = source[source["method"].eq("sourceheld_prespecified_late_fusion_w0.5")].sort_values("AUPRC", ascending=False)
if len(src):
    print("source-heldout prespecified late fusion:")
    print(src[["heldout_study","n","n_pos","AUPRC","AUROC","top10_precision"]].to_string(index=False))
print(f"KEEP/HOLD/KILL: {decision}")
PY
