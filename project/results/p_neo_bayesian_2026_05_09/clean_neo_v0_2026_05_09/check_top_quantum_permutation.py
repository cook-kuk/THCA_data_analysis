#!/usr/bin/env python3
"""Permutation sanity check for the top strict-set quantum-kernel fallback."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

import run_fallback_multiprocess_sweep as sweep


ROOT = Path(__file__).resolve().parent


def score_for_labels(x: np.ndarray, y: np.ndarray, gamma: float = 1.0) -> tuple[float, float]:
    pred = sweep.quantum_oof_predict(x, y, gamma)
    return float(roc_auc_score(y, pred)), float(average_precision_score(y, pred))


def main() -> None:
    df = sweep.load_data()
    features = tuple([*[c for c in sweep.QUANTUM_SCORES if c in df], *sweep.STRUCTURE_FEATS, "tcr_motif_score"])
    x = df[list(features)].to_numpy(float)
    y = df["label"].to_numpy(int)
    obs_auc, obs_auprc = score_for_labels(x, y)
    rng = np.random.default_rng(20260509)
    rows = [{"perm": "observed", "auroc": obs_auc, "auprc": obs_auprc}]
    n_perm = 100
    for i in range(n_perm):
        yp = rng.permutation(y)
        auc, auprc = score_for_labels(x, yp)
        rows.append({"perm": i, "auroc": auc, "auprc": auprc})
    out = pd.DataFrame(rows)
    out.to_csv(ROOT / "top_quantum_permutation.tsv", sep="\t", index=False)
    null = out[out["perm"] != "observed"]
    p_auc = (1 + (null["auroc"] >= obs_auc).sum()) / (1 + n_perm)
    p_auprc = (1 + (null["auprc"] >= obs_auprc).sum()) / (1 + n_perm)
    summary = f"""# Top Quantum Permutation Check

Candidate: `quantum_kernel_no_anchor_gamma1.0`

Features: `{','.join(features)}`

Observed AUROC: {obs_auc:.3f}
Observed AUPRC: {obs_auprc:.3f}

Permutation n={n_perm}

- AUROC null mean: {null['auroc'].mean():.3f}; 95% range: {null['auroc'].quantile(0.025):.3f}-{null['auroc'].quantile(0.975):.3f}; empirical p={p_auc:.3f}
- AUPRC null mean: {null['auprc'].mean():.3f}; 95% range: {null['auprc'].quantile(0.025):.3f}-{null['auprc'].quantile(0.975):.3f}; empirical p={p_auprc:.3f}

Interpretation: permutation checks whether this small strict-set internal CV signal is obviously label-random. It does not replace an external leakage-controlled benchmark.
"""
    (ROOT / "TOP_QUANTUM_PERMUTATION_CHECK.md").write_text(summary)
    print(summary)


if __name__ == "__main__":
    main()
