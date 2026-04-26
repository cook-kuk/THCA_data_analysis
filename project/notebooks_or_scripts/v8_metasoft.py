"""
v8 Task 2: Random-effects meta-analysis of per-cancer DIAL values.

For each classifier (5 classifiers x 5 cancers), treat DIAL as an effect size.
- Compute Hanley-McNeil SE(AUC_post) and approximate SE(DIAL) ~ SE(AUC_post).
- Random-effects meta using DerSimonian-Laird tau^2 estimator.
- Compute Cochran's Q, p-value (chi^2, k-1 df), I^2.
- Han 2011 m-values via Gaussian mixture posterior with prior 0.5.

Inputs:
  /opt/thyroid-dash/project/results/v5/v5p1_dial_all_cancers.tsv
  /opt/thyroid-dash/project/results/v5/v5p1_harmonization.tsv

Outputs:
  /opt/thyroid-dash/project/results/v8_statgen/v8_metasoft_results.tsv
  /opt/thyroid-dash/project/results/v8_statgen/v8_metasoft_forest_data.tsv

  When --bootstrap N is supplied (audit F6), per-cancer SE comes from a
  classical bootstrap (B=N resamples) over the 5 classifier DIAL values
  for that cancer instead of Hanley-McNeil SE(AUC_post). This is a coarse
  cell-level bootstrap; a fully defensible version would resample the
  per-fold (Y_bin, scores) tensors but those are not persisted to disk.
  The bootstrap-mode tables land at:
    v8_metasoft_results_bootstrap.tsv
    v8_metasoft_forest_data_bootstrap.tsv

Run:
  python -u /opt/thyroid-dash/project/notebooks_or_scripts/v8_metasoft.py
  python -u /opt/thyroid-dash/project/notebooks_or_scripts/v8_metasoft.py --bootstrap 1000
"""

from __future__ import annotations

import argparse
import math
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

PROJECT = Path("/opt/thyroid-dash/project")
DIAL_TSV = PROJECT / "results" / "v5" / "v5p1_dial_all_cancers.tsv"
HARM_TSV = PROJECT / "results" / "v5" / "v5p1_harmonization.tsv"
OUT_DIR = PROJECT / "results" / "v8_statgen"
OUT_RESULTS = OUT_DIR / "v8_metasoft_results.tsv"
OUT_RESULTS_BOOT = OUT_DIR / "v8_metasoft_results_bootstrap.tsv"
OUT_FOREST = OUT_DIR / "v8_metasoft_forest_data.tsv"
OUT_FOREST_BOOT = OUT_DIR / "v8_metasoft_forest_data_bootstrap.tsv"
LOG_DIR = PROJECT / "logs"

CANCERS = ["THCA", "SKCM", "LGG", "LUAD", "COAD"]
CLASSIFIERS = [
    "LogReg_l2",
    "RandomForest",
    "XGBoost",
    "LogReg_elasticnet",
    "GradientBoosting",
]


def hanley_mcneil_se(auc: float, n1: int, n2: int) -> float:
    """Hanley-McNeil standard error for an AUC.

    AUC is folded to [0.5, 1] first because the statistic is symmetric under
    label flip (a fundamental ambiguity for unlabelled/adversarial batches).
    """
    a = float(auc)
    if a < 0.5:
        a = 1.0 - a
    a = min(max(a, 0.5 + 1e-9), 1.0 - 1e-9)
    q1 = a / (2.0 - a)
    q2 = 2.0 * a * a / (1.0 + a)
    num = a * (1.0 - a) + (n1 - 1) * (q1 - a * a) + (n2 - 1) * (q2 - a * a)
    denom = n1 * n2
    if num <= 0 or denom <= 0:
        return float("nan")
    return math.sqrt(num / denom)


def dersimonian_laird(theta: np.ndarray, se: np.ndarray) -> dict:
    """DerSimonian-Laird random-effects meta-analysis.

    Returns a dict with mean_effect, se_mean, tau2, Q, pvalue_Q, I2,
    and the final RE weights.
    """
    theta = np.asarray(theta, dtype=float)
    se = np.asarray(se, dtype=float)
    k = theta.size
    if k < 2:
        raise ValueError("Need at least 2 studies for meta-analysis")

    w_fe = 1.0 / (se * se)
    theta_fe = np.sum(w_fe * theta) / np.sum(w_fe)

    q_stat = float(np.sum(w_fe * (theta - theta_fe) ** 2))
    df = k - 1
    c = float(np.sum(w_fe) - np.sum(w_fe * w_fe) / np.sum(w_fe))
    tau2 = max(0.0, (q_stat - df) / c) if c > 0 else 0.0

    w_re = 1.0 / (se * se + tau2)
    mean_effect = float(np.sum(w_re * theta) / np.sum(w_re))
    se_mean = float(math.sqrt(1.0 / np.sum(w_re)))

    pvalue_q = float(1.0 - stats.chi2.cdf(q_stat, df)) if df > 0 else float("nan")
    i2 = float(max(0.0, (q_stat - df) / q_stat)) if q_stat > 0 else 0.0

    return {
        "mean_effect": mean_effect,
        "se_mean": se_mean,
        "tau2": float(tau2),
        "Q": q_stat,
        "pvalue_Q": pvalue_q,
        "I2": i2,
        "w_re": w_re,
        "theta_fe": float(theta_fe),
    }


def han_m_values(
    theta: np.ndarray,
    se: np.ndarray,
    mean_effect: float,
    tau2: float,
    prior: float = 0.5,
) -> np.ndarray:
    """Han 2011 m-values (Gaussian mixture analog).

    For each study i:
      P(data_i | effect)    = N(theta_i; mean_effect, sqrt(se_i^2 + tau2))
      P(data_i | no effect) = N(theta_i; 0,           se_i)
    m_i = posterior P(effect present at study i | data_i) under prior 0.5.
    """
    theta = np.asarray(theta, dtype=float)
    se = np.asarray(se, dtype=float)

    sd_eff = np.sqrt(se * se + tau2)
    ll_eff = stats.norm.logpdf(theta, loc=mean_effect, scale=sd_eff)
    ll_null = stats.norm.logpdf(theta, loc=0.0, scale=se)

    log_prior_eff = math.log(prior)
    log_prior_null = math.log(1.0 - prior)

    log_num = log_prior_eff + ll_eff
    log_den = np.logaddexp(log_prior_eff + ll_eff, log_prior_null + ll_null)
    return np.exp(log_num - log_den)


def build_class_counts(harm: pd.DataFrame) -> dict:
    out = {}
    for _, row in harm.iterrows():
        out[str(row["cancer"])] = (int(row["n_class_A"]), int(row["n_class_B"]))
    return out


def bootstrap_se_per_cancer(dial: pd.DataFrame, B: int,
                            seed: int = 20260424) -> dict:
    """Per-cancer bootstrap SE for DIAL.

    For each cancer, we resample (with replacement) the 5 classifier DIAL
    values B times, compute the mean of each resample, and report the
    standard deviation of those means as SE_boot(cancer).

    This is a coarse cell-level bootstrap (5 classifier resampling, NOT a
    full LODO score resampling) because the underlying per-fold
    (Y_bin, scores) tensors were not persisted.
    """
    rng = np.random.default_rng(seed)
    out = {}
    for cancer, sub in dial.groupby("cancer"):
        vals = sub["dial"].astype(float).to_numpy()
        n = vals.size
        if n < 2:
            out[str(cancer)] = float("nan")
            continue
        # Resample indices [B, n], take mean over axis=1 -> [B] sample means.
        idx = rng.integers(0, n, size=(B, n))
        boot_means = vals[idx].mean(axis=1)
        # Use ddof=1 for an unbiased SD.
        out[str(cancer)] = float(np.std(boot_means, ddof=1))
    return out


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="v8 METASOFT-style meta-analysis.")
    p.add_argument(
        "--bootstrap",
        type=int,
        default=0,
        help=("If >0, replace Hanley-McNeil SE with per-cancer bootstrap SE "
              "estimated from the 5 classifier DIAL values (B resamples). "
              "Default 0 keeps the v8 original Hanley-McNeil behaviour."),
    )
    return p.parse_args()


def main() -> int:
    args = parse_args()
    use_bootstrap = args.bootstrap > 0
    B = int(args.bootstrap)
    print(f"[v8-metasoft] reading {DIAL_TSV}", flush=True)
    dial = pd.read_csv(DIAL_TSV, sep="\t")
    print(f"[v8-metasoft] reading {HARM_TSV}", flush=True)
    harm = pd.read_csv(HARM_TSV, sep="\t")

    class_counts = build_class_counts(harm)
    for c in CANCERS:
        if c not in class_counts:
            print(f"[v8-metasoft] WARN: no class counts for {c}", flush=True)

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    boot_se = None
    if use_bootstrap:
        boot_se = bootstrap_se_per_cancer(dial, B=B)
        print(f"[v8-metasoft] bootstrap mode ON (B={B}); per-cancer SE_boot:",
              flush=True)
        for cancer in CANCERS:
            print(f"    SE_boot[{cancer}] = {boot_se.get(cancer, float('nan')):.6f}",
                  flush=True)

    forest_rows = []
    results_rows = []

    for clf in CLASSIFIERS:
        sub = dial[dial["classifier"] == clf].copy()
        sub = sub.set_index("cancer")
        missing = [c for c in CANCERS if c not in sub.index]
        if missing:
            print(
                f"[v8-metasoft] WARN: classifier={clf} missing cancers={missing}",
                flush=True,
            )
            continue

        thetas = []
        ses = []
        for cancer in CANCERS:
            row = sub.loc[cancer]
            auc_post = float(row["auc_post"])
            dial_val = float(row["dial"])
            if use_bootstrap:
                se_use = boot_se.get(cancer, float("nan"))
            else:
                n1, n2 = class_counts[cancer]
                se_use = hanley_mcneil_se(auc_post, n1, n2)
            # Guard against degenerate zeros (perfect separation); use a small floor.
            if not math.isfinite(se_use) or se_use < 1e-4:
                se_use = 1e-4
            thetas.append(dial_val)
            ses.append(se_use)

        thetas = np.array(thetas, dtype=float)
        ses = np.array(ses, dtype=float)

        meta = dersimonian_laird(thetas, ses)
        m_vals = han_m_values(
            thetas, ses, meta["mean_effect"], meta["tau2"], prior=0.5
        )

        print(
            f"[v8-metasoft] {clf}: mean={meta['mean_effect']:.4f} "
            f"tau2={meta['tau2']:.4f} Q={meta['Q']:.2f} "
            f"p={meta['pvalue_Q']:.3g} I2={meta['I2']*100:.1f}%",
            flush=True,
        )
        for cancer, m in zip(CANCERS, m_vals):
            print(f"    m[{cancer}] = {m:.4f}", flush=True)

        row_out = {
            "classifier": clf,
            "mean_effect": meta["mean_effect"],
            "tau2": meta["tau2"],
            "Q": meta["Q"],
            "pvalue_Q": meta["pvalue_Q"],
            "I2": meta["I2"],
        }
        for cancer, m in zip(CANCERS, m_vals):
            row_out[f"m_{cancer}"] = float(m)
        results_rows.append(row_out)

        # Forest data: per-study 95% CI using Wald on theta with its own SE.
        for cancer, theta_i, se_i, m_i in zip(CANCERS, thetas, ses, m_vals):
            lo = theta_i - 1.96 * se_i
            hi = theta_i + 1.96 * se_i
            forest_rows.append(
                {
                    "classifier": clf,
                    "cancer": cancer,
                    "theta": float(theta_i),
                    "se": float(se_i),
                    "lower": float(lo),
                    "upper": float(hi),
                    "m_value": float(m_i),
                }
            )

    if not results_rows:
        print("[v8-metasoft] ERROR: no classifier results produced", flush=True)
        return 2

    col_order = [
        "classifier",
        "mean_effect",
        "tau2",
        "Q",
        "pvalue_Q",
        "I2",
        "m_THCA",
        "m_SKCM",
        "m_LGG",
        "m_LUAD",
        "m_COAD",
    ]
    df_res = pd.DataFrame(results_rows)[col_order]
    out_results_path = OUT_RESULTS_BOOT if use_bootstrap else OUT_RESULTS
    out_forest_path = OUT_FOREST_BOOT if use_bootstrap else OUT_FOREST
    df_res.to_csv(out_results_path, sep="\t", index=False, float_format="%.6g")
    print(f"[v8-metasoft] wrote {out_results_path}", flush=True)

    df_forest = pd.DataFrame(forest_rows)[
        ["classifier", "cancer", "theta", "se", "lower", "upper", "m_value"]
    ]
    df_forest.to_csv(out_forest_path, sep="\t", index=False, float_format="%.6g")
    print(f"[v8-metasoft] wrote {out_forest_path}", flush=True)

    # Echo the 5-row table so it shows up in the log.
    print("\n[v8-metasoft] results table:", flush=True)
    with pd.option_context("display.max_columns", None, "display.width", 200):
        print(df_res.to_string(index=False), flush=True)

    return 0


if __name__ == "__main__":
    sys.exit(main())
