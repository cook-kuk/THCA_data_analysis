"""M7 — Split conformal prediction over Wave 1 Bayesian outputs (Angelopoulos & Bates 2023).

Inputs:
  Wave 1 in-domain predictions: predictions_in_domain.tsv  (y, p_bayes, p_ens, p_std)
    → use as calibration
  Wave 1 ITSNdb predictions:    predictions_itsndb.tsv     (label, pred_mean, pred_std, in_master, split)
    → use as test

For each test point produce a 90 % prediction set (α=0.10) via the standard
binary-classification non-conformity score s(x, y) = 1 − p̂_y.

Reports:
  - Coverage rate (fraction of points whose set contains the true label)
  - Mean set size
  - AUROC of mean-prediction (sanity)
  - Selective AUROC: confidence-thresholded subset (set size = 1)
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

ROOT = Path(__file__).resolve().parent
WAVE1 = ROOT.parent

OUT_PRED = ROOT / "predictions_conformal.tsv"
OUT_LOG = ROOT / "m7_conformal.log"
ALPHA = 0.10


def conformity(y, p):
    """Non-conformity = 1 - p̂_y for binary class."""
    p_y = np.where(y == 1, p, 1 - p)
    return 1.0 - p_y


def prediction_set(p, q_hat):
    """For test prob p (P(y=1)), label y is in set iff 1 - p̂_y ≤ q_hat."""
    in_y1 = (1.0 - p) <= q_hat       # include label 1?
    in_y0 = (1.0 - (1.0 - p)) <= q_hat  # include label 0?
    return in_y0, in_y1


def main():
    log = open(OUT_LOG, "w")
    def L(*a):
        s = " ".join(str(x) for x in a)
        print(s); log.write(s + "\n"); log.flush()

    L("loading Wave 1 predictions…")
    cal = pd.read_csv(WAVE1 / "predictions_in_domain.tsv", sep="\t")
    test = pd.read_csv(WAVE1 / "predictions_itsndb.tsv", sep="\t")
    L(f"  calibration n={len(cal)} pos={int(cal['y'].sum())}")
    L(f"  itsndb test  n={len(test)} pos={int(test['label'].sum())}")

    # Use deep-ensemble prediction (p_ens) for calibration; fall back to p_bayes
    if "p_ens" in cal.columns:
        cal_p = cal["p_ens"].values
    else:
        cal_p = cal["p_bayes"].values
    cal_y = cal["y"].values.astype(int)
    test_p = test["pred_mean"].values
    test_y = test["label"].values.astype(int)

    # ----- Split conformal calibration -----------------------------------------
    rng = np.random.default_rng(0)
    n_cal_total = len(cal)
    # we don't actually need a 20% holdout because cal is already the held-out
    # 5-fold predictions on train pool — already disjoint from training data
    # (these were predicted on each fold's test split). Use them all.
    s_cal = conformity(cal_y, cal_p)
    n = len(s_cal)
    # finite-sample-corrected quantile: q = ceil((n+1)(1-α))/n
    q_level = np.ceil((n + 1) * (1 - ALPHA)) / n
    q_level = min(q_level, 1.0)
    q_hat = float(np.quantile(s_cal, q_level, method="higher"))
    L(f"\n  α={ALPHA}  q_level={q_level:.4f}  q_hat={q_hat:.4f}  (n_cal={n})")

    # ----- Test set: prediction sets ----------------------------------------------
    in_y0, in_y1 = prediction_set(test_p, q_hat)
    set_sizes = in_y0.astype(int) + in_y1.astype(int)
    contains_truth = np.where(test_y == 1, in_y1, in_y0)
    coverage = float(contains_truth.mean())
    mean_size = float(set_sizes.mean())
    L(f"\n  ALL ITSNdb (n={len(test)}):")
    L(f"    coverage={coverage:.4f}  (target ≥ {1-ALPHA:.2f})")
    L(f"    mean set size={mean_size:.3f}")
    L(f"    set-size dist: 0={int((set_sizes==0).sum())}  1={int((set_sizes==1).sum())}  2={int((set_sizes==2).sum())}")
    try:
        auc_all = roc_auc_score(test_y, test_p)
    except ValueError:
        auc_all = float("nan")
    L(f"    AUROC (mean-prediction)={auc_all:.4f}")

    # selective AUROC on points with set size = 1
    sel_mask = (set_sizes == 1)
    if sel_mask.sum() >= 5 and len(set(test_y[sel_mask])) == 2:
        auc_sel = float(roc_auc_score(test_y[sel_mask], test_p[sel_mask]))
    else:
        auc_sel = float("nan")
    L(f"    SELECTIVE AUROC (set size = 1): n_kept={int(sel_mask.sum())} retained={sel_mask.mean():.3f} AUROC={auc_sel:.4f}")

    # ----- Per-subset reporting --------------------------------------------------
    L("\n  per-subset:")
    in_master = test["in_master"].astype(bool).values
    rows = []
    for name, mask in [
        ("ITSNdb_main", test["split"] == "ext_itsndb_main"),
        ("ITSNdb_Val",  test["split"] == "ext_itsndb_val"),
        ("ITSNdb_combined", np.ones(len(test), dtype=bool)),
        ("ITSNdb_in_master", in_master),
        ("ITSNdb_no_overlap", ~in_master),
    ]:
        m_idx = np.where(mask)[0]
        if len(m_idx) < 5 or len(set(test_y[m_idx])) < 2:
            continue
        sub_y = test_y[m_idx]
        sub_p = test_p[m_idx]
        sub_size = set_sizes[m_idx]
        sub_cover = contains_truth[m_idx]
        sel = (sub_size == 1)
        if sel.sum() >= 5 and len(set(sub_y[sel])) == 2:
            auc_sel = float(roc_auc_score(sub_y[sel], sub_p[sel]))
        else:
            auc_sel = float("nan")
        try:
            auc_full = float(roc_auc_score(sub_y, sub_p))
        except ValueError:
            auc_full = float("nan")
        L(f"    {name:>22s}: n={len(m_idx):>4d} cov={sub_cover.mean():.3f} mean_size={sub_size.mean():.3f} AUROC={auc_full:.4f} sel_AUROC={auc_sel:.4f} (n_sel={int(sel.sum())})")
        rows.append({
            "subset": name, "n": len(m_idx), "coverage": sub_cover.mean(),
            "mean_set_size": sub_size.mean(), "AUROC": auc_full,
            "selective_AUROC": auc_sel, "n_selective": int(sel.sum()),
            "selective_retain_rate": float(sel.mean()),
        })

    # ----- Save predictions table ------------------------------------------------
    out = test.copy()
    out["pred_p"] = test_p
    out["set_in_y0"] = in_y0.astype(int)
    out["set_in_y1"] = in_y1.astype(int)
    out["set_size"] = set_sizes
    out["covers_truth"] = contains_truth.astype(int)
    out["q_hat"] = q_hat
    out["alpha"] = ALPHA
    out.to_csv(OUT_PRED, sep="\t", index=False)
    L(f"\n  saved → {OUT_PRED}")
    pd.DataFrame(rows).to_csv(ROOT / "m7_conformal_subsets.tsv", sep="\t", index=False)
    with open(ROOT / "m7_conformal_summary.json", "w") as f:
        json.dump({
            "alpha": ALPHA, "q_hat": q_hat, "n_cal": n,
            "coverage_all": coverage, "mean_set_size_all": mean_size,
        }, f, indent=2)

    log.close()
    print("\nM7 done.")


if __name__ == "__main__":
    main()
