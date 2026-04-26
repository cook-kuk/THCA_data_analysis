#!/usr/bin/env python3
"""v8.1 quantum finish — run the LUAD and COAD QSVC + VQC cells that
v8 S5 skipped under the 22-min wall-clock budget.

Imports the helper functions from v8_quantum_dial so the runs are
bit-for-bit identical (same PCA dims, same downsample seed, same
ZZFeatureMap / VQC ansatz). Only runs QSVC and VQC; SVC_RBF baselines
are already in v8_quantum_comparison.tsv.

Appends to results/v8_statgen/v8_quantum_comparison.tsv, de-duplicating
by (cancer, classifier_family) so re-runs overwrite the skipped rows.
"""
from __future__ import annotations
import sys, time
from pathlib import Path
import numpy as np
import pandas as pd

PROJECT = Path("/opt/thyroid-dash/project")
RES = PROJECT / "results" / "v8_statgen"
LOGS = PROJECT / "logs"
sys.path.insert(0, str(PROJECT / "notebooks_or_scripts"))

from v8_quantum_dial import (
    load_cancer, stratified_downsample,
    preprocess_to_pcs, combat_preserve_pcs,
    lodo_dial_for_family, N_CAP, N_QUBITS, RANDOM_SEED,
)

LOG = LOGS / "v8p1_quantum_finish.log"
_log_f = open(LOG, "w")
def log(msg: str) -> None:
    stamp = time.strftime("[%Y-%m-%d %H:%M:%S]")
    print(f"{stamp} {msg}", flush=True)
    _log_f.write(f"{stamp} {msg}\n"); _log_f.flush()

log("v8.1 quantum finish — running LUAD + COAD QSVC/VQC (no budget cap)")
log("=" * 60)

# Early import smoke-test
try:
    import qiskit_machine_learning  # noqa
    qsvc_ok = True
except Exception as e:
    log(f"qiskit-machine-learning import failed: {e}")
    qsvc_ok = False
try:
    import pennylane  # noqa
    vqc_ok = True
except Exception as e:
    log(f"pennylane import failed: {e}")
    vqc_ok = False

TARGET_CANCERS = ["LUAD", "COAD"]
TARGET_FAMILIES = ["QSVC", "VQC"]

new_rows = []
t_global = time.time()

for cancer in TARGET_CANCERS:
    log(f"--- {cancer} ---")
    X, Y, B = load_cancer(cancer)
    Xs, Ys, Bs = stratified_downsample(X, Y, B, N_CAP, RANDOM_SEED)
    n_used = Xs.shape[0]
    log(f"  downsampled n={n_used}; running preprocess_to_pcs + combat_preserve_pcs")
    X_pre = preprocess_to_pcs(Xs, N_QUBITS, RANDOM_SEED)
    X_post = combat_preserve_pcs(X_pre, Ys, Bs)
    classes_sorted = sorted(np.unique(Ys).tolist())
    Y_bin = (Ys == classes_sorted[0]).astype(int)

    for family in TARGET_FAMILIES:
        if family == "QSVC" and not qsvc_ok:
            log(f"  [{cancer}/{family}] SKIP (no qiskit-ml)")
            continue
        if family == "VQC" and not vqc_ok:
            log(f"  [{cancer}/{family}] SKIP (no pennylane)")
            continue
        elapsed = time.time() - t_global
        log(f"  [{cancer}/{family}] start  (global elapsed {elapsed:.0f}s)")
        try:
            ap, ao, dv, it, dt = lodo_dial_for_family(
                family, X_pre, X_post, Y_bin, Bs, RANDOM_SEED
            )
        except Exception as e:
            log(f"  [{cancer}/{family}] FATAL {type(e).__name__}: {e}")
            ap, ao, dv, it, dt = (
                float("nan"), float("nan"), float("nan"),
                f"error:{type(e).__name__}", 0.0,
            )
        log(f"  [{cancer}/{family}] auc_pre={ap:.3f} auc_post={ao:.3f} "
            f"dial={dv:.3f} {it}  ({dt:.1f}s)")
        new_rows.append(dict(
            cancer=cancer, classifier_family=family,
            auc_pre=ap, auc_post=ao, dial=dv,
            interpretation=it, n_used=n_used, seconds=round(dt, 2),
        ))

# Merge with existing TSV, de-dup by (cancer, family)
tsv = RES / "v8_quantum_comparison.tsv"
existing = pd.read_csv(tsv, sep="\t")
log(f"[merge] existing rows: {len(existing)}")
new_df = pd.DataFrame(new_rows)
log(f"[merge] new rows: {len(new_df)}")

if len(new_df) == 0:
    log("[merge] no new rows; skipping write")
else:
    combined = pd.concat([existing, new_df], ignore_index=True)
    # For (cancer, classifier_family), keep last (newer rows preferred)
    combined = combined.drop_duplicates(subset=["cancer", "classifier_family"],
                                         keep="last")
    # Restore canonical order: THCA, SKCM, LGG, LUAD, COAD × SVC_RBF, QSVC, VQC
    order_c = {"THCA":0, "SKCM":1, "LGG":2, "LUAD":3, "COAD":4}
    order_f = {"SVC_RBF":0, "QSVC":1, "VQC":2}
    combined["_oc"] = combined["cancer"].map(order_c)
    combined["_of"] = combined["classifier_family"].map(order_f)
    combined = combined.sort_values(["_oc", "_of"]).drop(columns=["_oc", "_of"])
    combined.to_csv(tsv, sep="\t", index=False, float_format="%.4f")
    log(f"[merge] wrote {tsv} with {len(combined)} rows")

dt_total = time.time() - t_global
log(f"v8.1 quantum finish complete ({dt_total:.1f}s wall)")
_log_f.close()
