"""Build paper-grade figures comparing Bayesian ESM2 vs RF baseline.

Inputs (all read from same dir):
  results_summary.tsv
  predictions_in_domain.tsv
  predictions_itsndb.tsv
  predictions_venus.tsv
  per_allele_loso_bayesian.tsv

External (from RF baseline directory):
  cancer_vaccine_robustness_2026_05_09/itsndb_results.tsv
  cancer_vaccine_robustness_2026_05_09/per_allele_loso.tsv
  cancer_vaccine_robustness_2026_05_09/venusvaccine_results.tsv
  cancer_vaccine_robustness_2026_05_09/source_balanced_loso.tsv

Outputs:
  fig_auroc_comparison.png/pdf
  fig_calibration.png/pdf
  fig_uncertainty_ood.png/pdf
  fig_per_allele_forest.png/pdf
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("Agg")
from sklearn.metrics import roc_auc_score
from sklearn.calibration import calibration_curve

ROOT = Path(__file__).parent
RF_DIR = Path("/home/seungho/personal/THCA_data_analysis/project/results/cancer_vaccine_robustness_2026_05_09")


def f1_auroc_comparison():
    """Bar chart: Bayesian vs RF on each external bench."""
    bayes = pd.read_csv(ROOT / "results_summary.tsv", sep="\t")
    rf_itsn = pd.read_csv(RF_DIR / "itsndb_results.tsv", sep="\t")
    rf_venus = pd.read_csv(RF_DIR / "venusvaccine_results.tsv", sep="\t")
    rf_loso = pd.read_csv(RF_DIR / "source_balanced_loso.tsv", sep="\t")

    points = []
    # ITSNdb subsets
    for sub in ["ITSNdb_main", "ITSNdb_Val", "ITSNdb_combined", "ITSNdb_no_overlap"]:
        rf_v = rf_itsn[rf_itsn["subset"] == sub]["AUROC"].values
        b_v = bayes[(bayes["eval"] == "ITSNdb") & (bayes.get("subset") == sub)]["AUROC"].values
        if len(rf_v) and len(b_v):
            points.append({"bench": sub, "RF": rf_v[0], "Bayes": b_v[0]})
    # VenusVaccine top10_mean
    for split in ["test", "valid"]:
        rf_v = rf_venus[(rf_venus["split"] == split) & (rf_venus["aggregator"] == "top10_mean")]["AUROC"].values
        b_v = bayes[(bayes["eval"] == "VenusVaccine") & (bayes.get("split") == f"ext_venus_{split}") & (bayes.get("aggregator") == "top10_mean")]["AUROC"].values
        if len(rf_v) and len(b_v):
            points.append({"bench": f"Venus_{split}_top10", "RF": rf_v[0], "Bayes": b_v[0]})
    # Cross-source LOSO mean
    if len(rf_loso) and "balanced" in rf_loso.columns:
        rf_mean = rf_loso["balanced"].mean()
    else:
        rf_mean = float("nan")
    b_loso = bayes[bayes["eval"] == "cross_source_loso"]["AUROC"].mean()
    points.append({"bench": "LOSO_mean", "RF": rf_mean, "Bayes": b_loso})
    # In-domain 5-fold
    rf_indom = 0.854  # from ROBUSTNESS_REPORT (within-source 5-fold AUROC)
    b_indom = bayes[(bayes["eval"] == "in_domain_5fold") & (bayes["model"] == "DeepEnsemble")]["AUROC"].mean()
    points.append({"bench": "InDomain_5fold", "RF": rf_indom, "Bayes": b_indom})

    pf = pd.DataFrame(points)
    pf.to_csv(ROOT / "auroc_comparison.tsv", sep="\t", index=False)

    fig, ax = plt.subplots(figsize=(11, 5))
    x = np.arange(len(pf)); w = 0.4
    ax.bar(x - w/2, pf["RF"], w, label="RF biophys+HLA (baseline)", color="#888")
    ax.bar(x + w/2, pf["Bayes"], w, label="Bayesian ESM2-150M", color="#1f77b4")
    ax.set_xticks(x); ax.set_xticklabels(pf["bench"].tolist(), rotation=30, ha="right")
    ax.axhline(0.5, ls="--", color="r", lw=1, label="chance")
    ax.set_ylabel("AUROC"); ax.set_ylim(0.3, 1.0)
    ax.legend(loc="lower right")
    ax.set_title("External benchmarks: Bayesian ESM2 vs RF biophys")
    plt.tight_layout()
    plt.savefig(ROOT / "fig_auroc_comparison.png", dpi=150)
    plt.savefig(ROOT / "fig_auroc_comparison.pdf")
    plt.close()
    print(f"saved fig_auroc_comparison")


def f2_calibration():
    pred = pd.read_csv(ROOT / "predictions_in_domain.tsv", sep="\t")
    if "p_ens" not in pred.columns or len(pred) < 30:
        print("calibration: insufficient preds, skipping"); return
    fig, ax = plt.subplots(figsize=(6, 6))
    for col, name in [("p_bayes", "Bayesian (1 seed, MC=30)"), ("p_ens", "Deep Ensemble (5 seeds × MC=30)")]:
        prob_true, prob_pred = calibration_curve(pred["y"], pred[col], n_bins=10, strategy="quantile")
        ax.plot(prob_pred, prob_true, "o-", label=name)
    ax.plot([0, 1], [0, 1], "k--", lw=1, label="ideal")
    ax.set_xlabel("Predicted probability"); ax.set_ylabel("Observed positive rate")
    ax.set_title("Calibration — in-domain 5-fold (ESM2 Bayesian)")
    ax.legend()
    plt.tight_layout()
    plt.savefig(ROOT / "fig_calibration.png", dpi=150)
    plt.savefig(ROOT / "fig_calibration.pdf")
    plt.close()
    print(f"saved fig_calibration")


def f3_uncertainty_ood():
    """Predictive entropy density: in-domain vs ITSNdb_no_overlap (proxy OOD)."""
    pred_in = pd.read_csv(ROOT / "predictions_in_domain.tsv", sep="\t")
    if not (ROOT / "predictions_itsndb.tsv").exists():
        print("OOD: no ITSNdb preds"); return
    pred_it = pd.read_csv(ROOT / "predictions_itsndb.tsv", sep="\t")
    # entropy from predictive mean
    def H(p):
        eps = 1e-9; p = np.clip(p, eps, 1-eps)
        return -(p*np.log(p) + (1-p)*np.log(1-p))
    h_in = H(pred_in["p_ens"].values)
    h_it_full = H(pred_it["pred_mean"].values)
    h_it_clean = H(pred_it[~pred_it["in_master"]]["pred_mean"].values) if "in_master" in pred_it.columns else h_it_full

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    ax = axes[0]
    bins = np.linspace(0, np.log(2), 30)
    ax.hist(h_in, bins=bins, alpha=0.5, label=f"In-domain n={len(h_in)}", density=True, color="#1f77b4")
    ax.hist(h_it_clean, bins=bins, alpha=0.5, label=f"ITSNdb_no_overlap n={len(h_it_clean)}", density=True, color="#ff7f0e")
    ax.set_xlabel("Predictive entropy"); ax.set_ylabel("Density")
    ax.set_title("Uncertainty: in-domain vs OOD (ITSNdb no_overlap)")
    ax.legend()

    # Predictive std (epistemic) density — Bayesian/ensemble std
    ax = axes[1]
    s_in = pred_in["p_std"].values
    s_it_clean = pred_it[~pred_it["in_master"]]["pred_std"].values if "in_master" in pred_it.columns else pred_it["pred_std"].values
    ax.hist(s_in, bins=30, alpha=0.5, label=f"In-domain n={len(s_in)}", density=True, color="#1f77b4")
    ax.hist(s_it_clean, bins=30, alpha=0.5, label=f"ITSNdb_no_overlap n={len(s_it_clean)}", density=True, color="#ff7f0e")
    ax.set_xlabel("Predictive std (epistemic)"); ax.set_ylabel("Density")
    ax.set_title("Predictive std")
    ax.legend()
    plt.tight_layout()
    plt.savefig(ROOT / "fig_uncertainty_ood.png", dpi=150)
    plt.savefig(ROOT / "fig_uncertainty_ood.pdf")
    plt.close()

    # OOD-AUROC: classify in-domain (0) vs OOD (1) by entropy
    y_ood = np.array([0]*len(h_in) + [1]*len(h_it_clean))
    s_ood = np.concatenate([h_in, h_it_clean])
    try:
        auc = roc_auc_score(y_ood, s_ood)
    except ValueError:
        auc = float("nan")
    print(f"OOD-AUROC by entropy (in-domain vs ITSNdb_no_overlap): {auc:.3f}")
    # Also try by predictive std
    s2 = np.concatenate([s_in, s_it_clean])
    try:
        auc2 = roc_auc_score(y_ood, s2)
    except ValueError:
        auc2 = float("nan")
    print(f"OOD-AUROC by predictive std: {auc2:.3f}")
    # Save
    (ROOT / "ood_auroc.json").write_text(f'{{"by_entropy": {auc:.4f}, "by_pred_std": {auc2:.4f}}}\n')
    print("saved fig_uncertainty_ood")


def f4_per_allele_forest():
    if not (ROOT / "per_allele_loso_bayesian.tsv").exists():
        print("per-allele: no preds"); return
    bayes = pd.read_csv(ROOT / "per_allele_loso_bayesian.tsv", sep="\t")
    rf = pd.read_csv(RF_DIR / "per_allele_loso.tsv", sep="\t")
    df = bayes.merge(rf[["allele", "AUROC", "n_test"]],
                     on="allele", suffixes=("_bayes", "_rf"))
    df = df.sort_values("AUROC_bayes")
    fig, ax = plt.subplots(figsize=(7, max(3, 0.4 * len(df) + 2)))
    y = np.arange(len(df))
    ax.scatter(df["AUROC_rf"], y, label="RF biophys+HLA", color="#888", s=70, marker="s")
    ax.scatter(df["AUROC_bayes"], y, label="Bayesian ESM2", color="#1f77b4", s=70, marker="o")
    for i, (rfa, ba) in enumerate(zip(df["AUROC_rf"], df["AUROC_bayes"])):
        ax.plot([rfa, ba], [i, i], color="#bbb", lw=1, zorder=0)
    ax.set_yticks(y); ax.set_yticklabels(df["allele"].tolist())
    ax.axvline(0.5, ls="--", color="r", lw=1)
    ax.set_xlim(0.3, 1.0); ax.set_xlabel("Per-allele LOSO AUROC")
    ax.set_title("Per-allele LOSO — Bayesian ESM2 vs RF baseline")
    ax.legend(loc="lower right")
    plt.tight_layout()
    plt.savefig(ROOT / "fig_per_allele_forest.png", dpi=150)
    plt.savefig(ROOT / "fig_per_allele_forest.pdf")
    plt.close()
    print(f"saved fig_per_allele_forest")


def main():
    f1_auroc_comparison()
    f2_calibration()
    f3_uncertainty_ood()
    f4_per_allele_forest()


if __name__ == "__main__":
    main()
