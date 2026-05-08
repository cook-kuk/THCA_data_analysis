"""Source-balanced LOSO — does inverse-frequency reweighting improve generalization?

Compares two regimes:
  A. Naive LOSO (no sample weights, raw concatenation)
  B. Source-balanced LOSO (sample_weight = 1 / n_source per row, normalized)

For each held-out source: train on remaining sources (with/without balanced
weights), score the held-out source, log AUROC.

Tested with biophys + HLA-onehot features (no ESM2 → external compatible).
"""
from __future__ import annotations
import sys, json
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from _common import (
    load_master_benchmark, filter_loso_eligible, build_hla_onehot_factory,
    build_features, fit_rf, metrics, cap_per_source, RESULTS,
)


def loso_run(df, balanced=False, cap=5000):
    """Run LOSO with or without source-balanced sample weights."""
    df = cap_per_source(df, cap=cap)
    rows = []
    sources = sorted(df["source"].unique())
    for s_held in sources:
        train = df[df["source"] != s_held].copy()
        test = df[df["source"] == s_held].copy()
        if len(test) < 30 or test["label"].nunique() < 2:
            continue
        if train["label"].sum() < 30 or (1 - train["label"]).sum() < 30:
            continue

        encode_hla, _, _ = build_hla_onehot_factory(train["HLA_norm"].tolist())
        Xtr = build_features(train, encode_hla)
        ytr = train["label"].values
        Xte = build_features(test, encode_hla)
        yte = test["label"].values

        sw = None
        if balanced:
            n_per_src = train["source"].value_counts().to_dict()
            w = train["source"].map(lambda x: 1.0 / n_per_src[x]).values
            sw = w * len(train) / w.sum()
        clf = fit_rf(Xtr, ytr, sample_weight=sw, n_estimators=200, max_depth=10)
        proba = clf.predict_proba(Xte)[:, 1]
        m = metrics(yte, proba)
        m["heldout_source"] = s_held
        m["n_train"] = len(train)
        m["n_test"] = len(test)
        m["mode"] = "source_balanced" if balanced else "naive"
        rows.append(m)
        print(f"  [{m['mode']:16s}] heldout={s_held:42s} n_test={len(test):>5d}  AUROC={m['AUROC']:.3f}")
    return rows


def main():
    print("=== Source-balanced LOSO benchmark ===")
    master = load_master_benchmark()
    df = filter_loso_eligible(master,
                              exclude_safety=("TRAINING_OVERLAP", "DEMO_ONLY", "PREDICTED_ONLY"),
                              min_n_per_source=50,
                              require_both_classes=True)
    print(f"LOSO-eligible: n={len(df):,}, sources={sorted(df['source'].unique())}")

    print("\n--- mode A: naive LOSO ---")
    rows_naive = loso_run(df, balanced=False)
    print("\n--- mode B: source-balanced LOSO ---")
    rows_bal = loso_run(df, balanced=True)

    all_rows = rows_naive + rows_bal
    out = pd.DataFrame(all_rows)
    out.to_csv(RESULTS / "source_balanced_loso.tsv", sep="\t", index=False)
    print(f"\nsaved → {RESULTS / 'source_balanced_loso.tsv'}")

    # Summary uplift
    if rows_naive and rows_bal:
        n_mean = float(np.nanmean([r["AUROC"] for r in rows_naive if r["AUROC"] is not None]))
        b_mean = float(np.nanmean([r["AUROC"] for r in rows_bal if r["AUROC"] is not None]))
        delta = b_mean - n_mean
        print(f"\nMean LOSO AUROC: naive={n_mean:.3f}  balanced={b_mean:.3f}  Δ={delta:+.3f}")

        summary = {
            "naive_mean_AUROC": n_mean,
            "balanced_mean_AUROC": b_mean,
            "delta": delta,
            "n_folds": len(rows_naive),
            "honest_note": ("Δ may be 0 or negative because biophys+HLA-onehot is dominated by HLA "
                           "feature; per-source distribution shift in HLA usage is the bottleneck. "
                           "Reweighting helps when source bias is in label prevalence; here it is "
                           "primarily in feature coverage."),
        }
        (RESULTS / "source_balanced_summary.json").write_text(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
