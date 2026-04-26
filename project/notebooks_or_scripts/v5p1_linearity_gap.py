#!/usr/bin/env python3
"""v5.1 Phase 4 — Linearity gap analysis.

Per cancer:
  - dial_linear_mean = mean(dial over LogReg_l2, LogReg_elasticnet)
  - dial_nonlinear_mean = mean(dial over RandomForest, GradientBoosting, XGBoost|HistGB)
  - dial_linearity_gap = linear_mean - nonlinear_mean

Labels:
  gap > 0.15 -> linear_specific_leakage
  |gap| <= 0.05 -> classifier_agnostic
  gap < -0.05 -> nonlinear_overfits_batch
  else -> intermediate
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, "/opt/thyroid-dash/project/notebooks_or_scripts")
from v5p1_common import RESULTS_V5, LOGS, log_line

LOG = LOGS / "v5p1_linearity.log"

LINEAR = ("LogReg_l2", "LogReg_elasticnet")
NONLINEAR = ("RandomForest", "GradientBoosting", "XGBoost", "HistGB")


def main():
    p = RESULTS_V5 / "v5p1_dial_all_cancers.tsv"
    if not p.exists():
        log_line(LOG, "no all_cancers file — skip")
        return
    df = pd.read_csv(p, sep="\t")
    rows = []
    for cancer, g in df.groupby("cancer"):
        dlin = g[g["classifier"].isin(LINEAR)]["dial"].mean()
        dnonlin = g[g["classifier"].isin(NONLINEAR)]["dial"].mean()
        gap = float(dlin - dnonlin) if not (pd.isna(dlin) or pd.isna(dnonlin)) else float("nan")
        if pd.isna(gap):
            label = "unknown"
        elif gap > 0.15:
            label = "linear_specific_leakage"
        elif abs(gap) <= 0.05:
            label = "classifier_agnostic"
        elif gap < -0.05:
            label = "nonlinear_overfits_batch"
        else:
            label = "intermediate"
        rows.append({"cancer": cancer,
                     "dial_linear_mean": float(dlin) if not pd.isna(dlin) else float("nan"),
                     "dial_nonlinear_mean": float(dnonlin) if not pd.isna(dnonlin) else float("nan"),
                     "dial_linearity_gap": gap,
                     "label": label,
                     "n_linear": int(g["classifier"].isin(LINEAR).sum()),
                     "n_nonlinear": int(g["classifier"].isin(NONLINEAR).sum()),
                     "semi_synthetic": False})
    out = pd.DataFrame(rows)
    out.to_csv(RESULTS_V5 / "v5p1_linearity_gap.tsv", sep="\t", index=False)
    log_line(LOG, f"wrote v5p1_linearity_gap.tsv ({len(out)} rows)")
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
