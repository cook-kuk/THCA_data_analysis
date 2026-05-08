"""Per-allele LOSO — hold out each major HLA allele, train on rest, score.

Reveals which alleles the model generalizes to vs which collapse. The
dominant alleles (A*02:01, A*11:01, A*24:02, B*07:02) are well-represented
in training; minor alleles are not.
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


def main():
    print("=== Per-allele LOSO benchmark ===")
    master = load_master_benchmark()
    df = filter_loso_eligible(master,
                              exclude_safety=("TRAINING_OVERLAP", "DEMO_ONLY", "PREDICTED_ONLY"),
                              min_n_per_source=50,
                              require_both_classes=True)
    df = cap_per_source(df, cap=5000)
    print(f"LOSO-eligible: n={len(df):,}")

    # Allele frequency table
    allele_counts = df.groupby("HLA_norm")["label"].agg(["count", "sum", "mean"])
    allele_counts.columns = ["n", "n_pos", "pos_rate"]
    allele_counts = allele_counts.sort_values("n", ascending=False)
    print(f"\nTop alleles in pool:\n{allele_counts.head(15)}")

    # Cumulative coverage to ≥80%
    allele_counts["cum_frac"] = allele_counts["n"].cumsum() / allele_counts["n"].sum()
    eligible = allele_counts[(allele_counts["n"] >= 50) &
                             (allele_counts["pos_rate"] > 0.05) &
                             (allele_counts["pos_rate"] < 0.95)]
    print(f"\n>=50 examples + both classes (5–95% pos): {len(eligible)} alleles")

    rows = []
    for allele in eligible.index:
        train = df[df["HLA_norm"] != allele].copy()
        test = df[df["HLA_norm"] == allele].copy()
        if len(test) < 30 or test["label"].nunique() < 2:
            continue
        if train["label"].sum() < 30 or (1 - train["label"]).sum() < 30:
            continue

        # Source-balanced sample weights
        n_per_src = train["source"].value_counts().to_dict()
        w = train["source"].map(lambda x: 1.0 / n_per_src[x]).values
        sw = w * len(train) / w.sum()

        encode_hla, _, _ = build_hla_onehot_factory(train["HLA_norm"].tolist())
        Xtr = build_features(train, encode_hla)
        ytr = train["label"].values
        Xte = build_features(test, encode_hla)
        yte = test["label"].values

        clf = fit_rf(Xtr, ytr, sample_weight=sw, n_estimators=200, max_depth=10)
        proba = clf.predict_proba(Xte)[:, 1]
        m = metrics(yte, proba)
        m["allele"] = allele
        m["n_test"] = len(test)
        m["n_pos_test"] = int(yte.sum())
        m["pos_rate_test"] = float(yte.mean())
        rows.append(m)
        print(f"  {allele:18s} n_test={len(test):>4d}  pos={yte.sum():>3d}  AUROC={m['AUROC']:.3f}  AUPRC={m['AUPRC']:.3f}")

    out = pd.DataFrame(rows).sort_values("AUROC", ascending=False)
    out.to_csv(RESULTS / "per_allele_loso.tsv", sep="\t", index=False)
    print(f"\nsaved → {RESULTS / 'per_allele_loso.tsv'}")

    # Top / bottom 3
    valid = out[out["AUROC"].notna()]
    if len(valid) >= 3:
        print(f"\nTop 3 alleles (best generalization):")
        print(valid.head(3)[["allele", "n_test", "n_pos_test", "AUROC"]])
        print(f"\nBottom 3 alleles (worst generalization):")
        print(valid.tail(3)[["allele", "n_test", "n_pos_test", "AUROC"]])

    summary = {
        "n_alleles_tested": len(out),
        "mean_AUROC": float(np.nanmean(out["AUROC"].values)) if len(out) else None,
        "median_AUROC": float(np.nanmedian(out["AUROC"].values)) if len(out) else None,
        "min_AUROC": float(np.nanmin(out["AUROC"].values)) if len(out) else None,
        "max_AUROC": float(np.nanmax(out["AUROC"].values)) if len(out) else None,
        "top3_alleles": valid.head(3)["allele"].tolist() if len(valid) >= 3 else [],
        "bottom3_alleles": valid.tail(3)["allele"].tolist() if len(valid) >= 3 else [],
    }
    (RESULTS / "per_allele_summary.json").write_text(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
