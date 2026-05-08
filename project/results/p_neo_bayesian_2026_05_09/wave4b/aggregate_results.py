"""Aggregate M5/M6/M7/M8 + Wave 1 baseline into wave4b_results.tsv with bootstrap CIs,
then plot fig_wave4b_uplift.png/pdf.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.metrics import roc_auc_score

ROOT = Path(__file__).resolve().parent
WAVE1 = ROOT.parent

OUT_TSV = ROOT / "wave4b_results.tsv"
OUT_FIG_PNG = ROOT / "fig_wave4b_uplift.png"
OUT_FIG_PDF = ROOT / "fig_wave4b_uplift.pdf"

WAVE1_BASELINE_NO_OVERLAP = 0.4106514084507042  # from results_summary.tsv
MHCFLURRY_NO_OVERLAP = 0.668                    # 6-algorithm forest baseline
N_BOOT = 1000
RNG_SEED = 0


def auroc_with_ci(y, p, n_boot=N_BOOT, seed=RNG_SEED):
    y = np.asarray(y).astype(int); p = np.asarray(p, dtype=float)
    if len(set(y)) < 2 or len(y) < 5:
        return float("nan"), float("nan"), float("nan")
    auc = roc_auc_score(y, p)
    rng = np.random.default_rng(seed)
    boots = []
    n = len(y)
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        if len(set(y[idx])) < 2:
            continue
        boots.append(roc_auc_score(y[idx], p[idx]))
    if not boots:
        return float(auc), float("nan"), float("nan")
    lo, hi = np.percentile(boots, [2.5, 97.5])
    return float(auc), float(lo), float(hi)


def collect():
    rows = []

    # ----- Wave 1 baseline (reference) ----------------------------------------
    base = pd.read_csv(WAVE1 / "predictions_itsndb.tsv", sep="\t")
    in_master = base["in_master"].astype(bool).values
    for testset, mask in [
        ("ITSNdb_main", base["split"] == "ext_itsndb_main"),
        ("ITSNdb_Val",  base["split"] == "ext_itsndb_val"),
        ("ITSNdb_combined", np.ones(len(base), dtype=bool)),
        ("ITSNdb_in_master",  in_master),
        ("ITSNdb_no_overlap", ~in_master),
    ]:
        sub = base[mask]
        if len(sub) < 5 or sub["label"].nunique() < 2:
            continue
        auc, lo, hi = auroc_with_ci(sub["label"], sub["pred_mean"])
        rows.append({
            "method": "Wave1_baseline",
            "testset": testset,
            "in_master": "TRUE" if testset.endswith("in_master") else ("FALSE" if testset.endswith("no_overlap") else "BOTH"),
            "n": len(sub),
            "AUROC": auc, "CI_lo95": lo, "CI_hi95": hi,
            "delta_vs_baseline": 0.0,
        })

    # ----- M5 kNN -------------------------------------------------------------
    p_knn = ROOT / "predictions_knn.tsv"
    if p_knn.exists():
        df = pd.read_csv(p_knn, sep="\t")
        in_master = df["in_master"].astype(bool).values
        for testset, mask in [
            ("ITSNdb_main", df["split"] == "ext_itsndb_main"),
            ("ITSNdb_Val",  df["split"] == "ext_itsndb_val"),
            ("ITSNdb_combined", np.ones(len(df), dtype=bool)),
            ("ITSNdb_in_master",  in_master),
            ("ITSNdb_no_overlap", ~in_master),
        ]:
            sub = df[mask]
            if len(sub) < 5 or sub["label"].nunique() < 2:
                continue
            auc, lo, hi = auroc_with_ci(sub["label"], sub["pred_knn"])
            base_row = next((r for r in rows if r["method"] == "Wave1_baseline" and r["testset"] == testset), None)
            base_auc = base_row["AUROC"] if base_row else float("nan")
            rows.append({
                "method": "M5_kNN",
                "testset": testset,
                "in_master": "TRUE" if testset.endswith("in_master") else ("FALSE" if testset.endswith("no_overlap") else "BOTH"),
                "n": len(sub),
                "AUROC": auc, "CI_lo95": lo, "CI_hi95": hi,
                "delta_vs_baseline": auc - base_auc if not np.isnan(base_auc) else float("nan"),
            })

    # ----- M6 TENT ------------------------------------------------------------
    p_tent = ROOT / "predictions_tent.tsv"
    if p_tent.exists():
        df = pd.read_csv(p_tent, sep="\t")
        in_master = df["in_master"].astype(bool).values
        for testset, mask in [
            ("ITSNdb_main", df["split"] == "ext_itsndb_main"),
            ("ITSNdb_Val",  df["split"] == "ext_itsndb_val"),
            ("ITSNdb_combined", np.ones(len(df), dtype=bool)),
            ("ITSNdb_in_master",  in_master),
            ("ITSNdb_no_overlap", ~in_master),
        ]:
            sub = df[mask]
            if len(sub) < 5 or sub["label"].nunique() < 2:
                continue
            # Report TENT (not the surrogate baseline) under the M6 row
            auc, lo, hi = auroc_with_ci(sub["label"], sub["pred_tent"])
            base_row = next((r for r in rows if r["method"] == "Wave1_baseline" and r["testset"] == testset), None)
            base_auc = base_row["AUROC"] if base_row else float("nan")
            rows.append({
                "method": "M6_TENT",
                "testset": testset,
                "in_master": "TRUE" if testset.endswith("in_master") else ("FALSE" if testset.endswith("no_overlap") else "BOTH"),
                "n": len(sub),
                "AUROC": auc, "CI_lo95": lo, "CI_hi95": hi,
                "delta_vs_baseline": auc - base_auc if not np.isnan(base_auc) else float("nan"),
            })

    # ----- M7 conformal -------------------------------------------------------
    p_conf = ROOT / "predictions_conformal.tsv"
    if p_conf.exists():
        df = pd.read_csv(p_conf, sep="\t")
        in_master = df["in_master"].astype(bool).values
        for testset, mask in [
            ("ITSNdb_main", df["split"] == "ext_itsndb_main"),
            ("ITSNdb_Val",  df["split"] == "ext_itsndb_val"),
            ("ITSNdb_combined", np.ones(len(df), dtype=bool)),
            ("ITSNdb_in_master",  in_master),
            ("ITSNdb_no_overlap", ~in_master),
        ]:
            sub = df[mask]
            if len(sub) < 5 or sub["label"].nunique() < 2:
                continue
            auc, lo, hi = auroc_with_ci(sub["label"], sub["pred_p"])
            base_row = next((r for r in rows if r["method"] == "Wave1_baseline" and r["testset"] == testset), None)
            base_auc = base_row["AUROC"] if base_row else float("nan")
            rows.append({
                "method": "M7_conformal",
                "testset": testset,
                "in_master": "TRUE" if testset.endswith("in_master") else ("FALSE" if testset.endswith("no_overlap") else "BOTH"),
                "n": len(sub),
                "AUROC": auc, "CI_lo95": lo, "CI_hi95": hi,
                "delta_vs_baseline": auc - base_auc if not np.isnan(base_auc) else float("nan"),
                "coverage": float(sub["covers_truth"].mean()),
                "mean_set_size": float(sub["set_size"].mean()),
            })

    # ----- M8 WiSE-FT (skipped) -----------------------------------------------
    p_w = ROOT / "predictions_wisefit.tsv"
    if p_w.exists():
        df = pd.read_csv(p_w, sep="\t")
        # M8 is α=0 placeholder = Wave 1 baseline numerically. Tag as such.
        in_master = df["in_master"].astype(bool).values
        for testset, mask in [
            ("ITSNdb_main", df["split"] == "ext_itsndb_main"),
            ("ITSNdb_Val",  df["split"] == "ext_itsndb_val"),
            ("ITSNdb_combined", np.ones(len(df), dtype=bool)),
            ("ITSNdb_in_master",  in_master),
            ("ITSNdb_no_overlap", ~in_master),
        ]:
            sub = df[mask]
            if len(sub) < 5 or sub["label"].nunique() < 2:
                continue
            auc, lo, hi = auroc_with_ci(sub["label"], sub["pred_wisefit"])
            base_row = next((r for r in rows if r["method"] == "Wave1_baseline" and r["testset"] == testset), None)
            base_auc = base_row["AUROC"] if base_row else float("nan")
            rows.append({
                "method": "M8_WiSE-FT_skipped",
                "testset": testset,
                "in_master": "TRUE" if testset.endswith("in_master") else ("FALSE" if testset.endswith("no_overlap") else "BOTH"),
                "n": len(sub),
                "AUROC": auc, "CI_lo95": lo, "CI_hi95": hi,
                "delta_vs_baseline": auc - base_auc if not np.isnan(base_auc) else 0.0,
            })
    return rows


def plot_uplift(df_results):
    sub = df_results[df_results["testset"] == "ITSNdb_no_overlap"].copy()
    methods = sub["method"].tolist()
    aucs = sub["AUROC"].astype(float).values
    lo = sub["CI_lo95"].astype(float).values
    hi = sub["CI_hi95"].astype(float).values
    err_lo = aucs - lo
    err_hi = hi - aucs

    fig, ax = plt.subplots(figsize=(8.5, 4.5))
    x = np.arange(len(methods))
    colors = ["#888888" if "Wave1" in m or "skipped" in m else "#1f77b4" for m in methods]
    bars = ax.bar(x, aucs, yerr=[err_lo, err_hi], capsize=4, color=colors, edgecolor="black", linewidth=0.7)
    ax.axhline(WAVE1_BASELINE_NO_OVERLAP, color="#444444", linestyle="--", linewidth=0.9, label=f"Wave 1 baseline ({WAVE1_BASELINE_NO_OVERLAP:.3f})")
    ax.axhline(MHCFLURRY_NO_OVERLAP, color="#d62728", linestyle="--", linewidth=0.9, label=f"MHCflurry ({MHCFLURRY_NO_OVERLAP:.3f})")
    ax.set_xticks(x); ax.set_xticklabels(methods, rotation=20, ha="right")
    ax.set_ylim(0.0, 1.0)
    ax.set_ylabel("AUROC (ITSNdb_no_overlap)")
    ax.set_title("Wave 4B inference-time methods — AUROC on ITSNdb no_overlap")
    ax.grid(axis="y", linestyle=":", linewidth=0.5, alpha=0.6)
    for xi, a in zip(x, aucs):
        ax.text(xi, a + 0.02, f"{a:.3f}", ha="center", va="bottom", fontsize=9)
    ax.legend(loc="upper right", fontsize=9)
    fig.tight_layout()
    fig.savefig(OUT_FIG_PNG, dpi=150)
    fig.savefig(OUT_FIG_PDF)
    plt.close(fig)


def main():
    rows = collect()
    df = pd.DataFrame(rows)
    df.to_csv(OUT_TSV, sep="\t", index=False)
    print(f"saved {len(df)} rows → {OUT_TSV}")
    print(df[df["testset"] == "ITSNdb_no_overlap"][["method", "n", "AUROC", "CI_lo95", "CI_hi95", "delta_vs_baseline"]].to_string(index=False))
    plot_uplift(df)
    print(f"saved {OUT_FIG_PNG}\nsaved {OUT_FIG_PDF}")


if __name__ == "__main__":
    main()
