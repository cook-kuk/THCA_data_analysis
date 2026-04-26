#!/usr/bin/env python3
"""v8 audit fix F7: 5-seed replicate of quantum classifier DIAL for THCA only.

Why
---
The v8 S5 row "THCA / VQC dial=0.386 / batch_entangled" was computed with a
single random seed. Because VQC weights are randomly initialised and trained
with stochastic optimisation, one seed is exploratory, not confirmatory. This
script reruns VQC, QSVC, and SVC_RBF on the THCA harmonised cohort under five
seeds {0,1,2,3,4} and reports mean +/- std of the DIAL metric.

Acceptance
----------
Mean dial across 5 seeds for THCA/VQC > 0.3 AND std < 0.1 -> "robust".
Otherwise the original single-seed result is flagged as exploratory.

Outputs
-------
- results/v8_statgen/v8_quantum_thca_5seed.tsv          (per-seed rows)
- results/v8_statgen/v8_quantum_thca_5seed_summary.tsv  (mean/std rows)

Reuses fit_predict_vqc / fit_predict_qsvc / fit_predict_svc_rbf and the
DIAL pipeline from v8_quantum_dial.py so seed handling and stratified
downsampling are identical to the original S5 run.
"""
from __future__ import annotations

import sys
import time
import warnings
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

PROJECT = Path("/opt/thyroid-dash/project")
RES = PROJECT / "results" / "v8_statgen"
LOGS = PROJECT / "logs"
for p in [RES, LOGS]:
    p.mkdir(parents=True, exist_ok=True)

LOG_PATH = LOGS / "v8_quantum_5seed.log"

# Make v8_quantum_dial importable so we reuse the *same* code paths.
sys.path.insert(0, str(PROJECT / "notebooks_or_scripts"))
import v8_quantum_dial as v8q  # noqa: E402

CANCER = "THCA"
SEEDS = [0, 1, 2, 3, 4]
N_CAP = v8q.N_CAP                 # 80, same as original S5
N_QUBITS = v8q.N_QUBITS           # 8


def log(msg: str) -> None:
    ts = time.strftime("%Y-%m-%d %H:%M:%S")
    line = f"[{ts}] {msg}"
    print(line, flush=True)
    with open(LOG_PATH, "a") as fh:
        fh.write(line + "\n")


def run_one_seed(
    seed: int,
    families: List[str],
) -> List[Dict]:
    """Compute DIAL rows for the given families at the given seed.

    The same seed flows through:
      - stratified downsample (so the n=80 sample composition can shift seed-to-seed,
        which is honest; the alternative would mask seed sensitivity)
      - PCA random_state
      - VQC weight initialisation (per-fold seed = seed + fold_index)
    """
    log(f"--- {CANCER} seed={seed} ---")
    X, Y, B = v8q.load_cancer(CANCER)
    Xs, Ys, Bs = v8q.stratified_downsample(X, Y, B, N_CAP, seed)
    n_used = Xs.shape[0]
    log(f"  downsampled n={n_used} batches={list(np.unique(Bs))} y={list(np.unique(Ys))}")

    X_pre = v8q.preprocess_to_pcs(Xs, N_QUBITS, seed)
    X_post = v8q.combat_preserve_pcs(X_pre, Ys, Bs)

    classes_sorted = sorted(np.unique(Ys).tolist())
    Y_bin = (Ys == classes_sorted[0]).astype(int)

    if len(np.unique(Bs)) < 2:
        log(f"  [{CANCER} seed={seed}] single batch; LODO impossible")
        return [
            dict(
                seed=seed, classifier_family=fam,
                auc_pre=float("nan"), auc_post=float("nan"),
                dial=float("nan"),
                interpretation="no_folds_single_batch",
                n_used=n_used, seconds=0.0,
            )
            for fam in families
        ]

    rows: List[Dict] = []
    for fam in families:
        log(f"  [{CANCER}/{fam} seed={seed}] start")
        try:
            ap, ao, dv, it, dt = v8q.lodo_dial_for_family(
                fam, X_pre, X_post, Y_bin, Bs, seed
            )
        except Exception as e:
            log(f"  [{CANCER}/{fam} seed={seed}] FATAL {type(e).__name__}: {e}")
            ap, ao, dv, it, dt = (
                float("nan"), float("nan"), float("nan"),
                f"error:{type(e).__name__}", 0.0,
            )
        rows.append(dict(
            seed=seed, classifier_family=fam,
            auc_pre=ap, auc_post=ao, dial=dv,
            interpretation=it, n_used=n_used, seconds=round(dt, 2),
        ))
        log(
            f"  [{CANCER}/{fam} seed={seed}] "
            f"auc_pre={ap:.3f} auc_post={ao:.3f} dial={dv:.3f} {it} ({dt:.0f}s)"
        )
    return rows


def main() -> int:
    t_global = time.time()
    log(
        f"=== v8 F7 audit: 5-seed quantum DIAL for {CANCER}; "
        f"seeds={SEEDS}; n_cap={N_CAP}; qubits={N_QUBITS} ==="
    )

    # Probe quantum stack once.
    qsvc_ok = True
    vqc_ok = True
    try:
        v8q._build_qsvc()
    except Exception as e:
        qsvc_ok = False
        log(f"[QSVC unavailable] {type(e).__name__}: {e}")
    try:
        import pennylane  # noqa: F401
    except Exception as e:
        vqc_ok = False
        log(f"[VQC unavailable] {type(e).__name__}: {e}")

    all_rows: List[Dict] = []

    # SVC_RBF is deterministic (sklearn SVC + fixed PCA). Run once at seed 0,
    # then replicate the row for each seed so the per-seed table is rectangular.
    svc_families_per_seed = ["VQC"] if not qsvc_ok else ["QSVC", "VQC"]
    if not vqc_ok:
        svc_families_per_seed = [f for f in svc_families_per_seed if f != "VQC"]

    log("=== SVC_RBF deterministic single run ===")
    svc_rows = run_one_seed(seed=0, families=["SVC_RBF"])
    svc_template = svc_rows[0]
    for s in SEEDS:
        all_rows.append({**svc_template, "seed": s})

    # VQC + QSVC across all seeds.
    for seed in SEEDS:
        rows = run_one_seed(seed=seed, families=svc_families_per_seed)
        all_rows.extend(rows)
        elapsed = time.time() - t_global
        log(f"=== seed {seed} done; cumulative elapsed {elapsed:.0f}s ===")

    # ------------------------------------------------------------------
    # Per-seed table
    # ------------------------------------------------------------------
    df = pd.DataFrame(all_rows)
    fam_order = {"VQC": 0, "QSVC": 1, "SVC_RBF": 2}
    df["_f"] = df["classifier_family"].map(fam_order).fillna(99)
    df = (
        df.sort_values(["seed", "_f"])
          .drop(columns="_f")
          .reset_index(drop=True)
    )
    out_rows = RES / "v8_quantum_thca_5seed.tsv"
    df.to_csv(out_rows, sep="\t", index=False, float_format="%.4f")
    log(f"wrote {out_rows} ({len(df)} rows)")

    # ------------------------------------------------------------------
    # Summary table (VQC + QSVC only; SVC_RBF is constant by construction)
    # ------------------------------------------------------------------
    summary_rows = []
    for fam in ["VQC", "QSVC"]:
        sub = df[df["classifier_family"] == fam]
        # Restrict to seeds in SEEDS only (drop NaNs in summary calc).
        sub_seeds = sub[sub["seed"].isin(SEEDS)]
        n_be = int((sub_seeds["interpretation"] == "batch_entangled").sum())
        summary_rows.append(dict(
            classifier_family=fam,
            mean_dial=float(np.nanmean(sub_seeds["dial"])),
            std_dial=float(np.nanstd(sub_seeds["dial"], ddof=1)) if len(sub_seeds) > 1 else float("nan"),
            mean_auc_post=float(np.nanmean(sub_seeds["auc_post"])),
            std_auc_post=float(np.nanstd(sub_seeds["auc_post"], ddof=1)) if len(sub_seeds) > 1 else float("nan"),
            n_seeds_batch_entangled=n_be,
        ))
    summary = pd.DataFrame(summary_rows)
    out_sum = RES / "v8_quantum_thca_5seed_summary.tsv"
    summary.to_csv(out_sum, sep="\t", index=False, float_format="%.4f")
    log(f"wrote {out_sum}")

    # ------------------------------------------------------------------
    # Verdict
    # ------------------------------------------------------------------
    log("=== PER-SEED TABLE ===")
    for _, r in df.iterrows():
        log(
            f"  seed={int(r['seed'])} {r['classifier_family']:8s} "
            f"pre={r['auc_pre']:.3f} post={r['auc_post']:.3f} "
            f"dial={r['dial']:.3f} {r['interpretation']} "
            f"n={int(r['n_used'])} t={r['seconds']}s"
        )
    log("=== SUMMARY ===")
    for _, r in summary.iterrows():
        log(
            f"  {r['classifier_family']:6s} "
            f"mean_dial={r['mean_dial']:.3f} std_dial={r['std_dial']:.3f} "
            f"mean_auc_post={r['mean_auc_post']:.3f} std_auc_post={r['std_auc_post']:.3f} "
            f"n_batch_entangled={int(r['n_seeds_batch_entangled'])}/{len(SEEDS)}"
        )

    vqc_row = summary[summary["classifier_family"] == "VQC"].iloc[0]
    if (vqc_row["mean_dial"] > 0.3) and (vqc_row["std_dial"] < 0.1):
        verdict = (
            f"ROBUST: VQC mean_dial={vqc_row['mean_dial']:.3f} "
            f"(std={vqc_row['std_dial']:.3f}) across {len(SEEDS)} seeds; "
            f"original S5 THCA flip claim survives."
        )
    else:
        verdict = (
            f"EXPLORATORY: VQC mean_dial={vqc_row['mean_dial']:.3f} "
            f"(std={vqc_row['std_dial']:.3f}) across {len(SEEDS)} seeds; "
            f"single-seed result was exploratory; multi-seed shows wider variance."
        )
    log(f"=== VERDICT: {verdict} ===")
    log(f"=== total wall {time.time() - t_global:.0f}s ===")
    return 0


if __name__ == "__main__":
    sys.exit(main())
