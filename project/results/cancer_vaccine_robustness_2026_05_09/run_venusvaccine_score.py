"""VenusVaccine TumorBinary (ICLR 2025) — protein-level antigen benchmark.

VenusVaccine antigens are FULL PROTEINS (mean length ~500 AAs), not peptides.
Strategy: enumerate 8-11mer windows, score each window with the source-balanced
RF, aggregate to protein-level via max + mean, then compute protein-level AUROC.

This is a stress-test of generalization: if the peptide-level model captures
real immunogenicity signal, top-scoring windows should concentrate in
"protective" antigens (label=1).

Output: venusvaccine_results.tsv
"""
from __future__ import annotations
import sys, json
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from _common import (
    load_master_benchmark, filter_loso_eligible, build_hla_onehot_factory,
    build_features, biophys, fit_rf, metrics, normalize_hla, cap_per_source,
    EXT, RESULTS,
)

# Most common HLA alleles in master cohort (use as evaluation panel)
EVAL_HLAS = [
    "HLA-A*02:01", "HLA-A*01:01", "HLA-A*11:01", "HLA-A*03:01",
    "HLA-A*24:02", "HLA-B*07:02", "HLA-B*08:01", "HLA-B*15:01",
    "HLA-B*44:02", "HLA-B*35:01"
]


def enumerate_windows(seq, lengths=(9, 10, 11)):
    """Yield (window_seq, start_pos, length)."""
    s = "".join(c for c in str(seq).upper() if c.isalpha())
    for L in lengths:
        for i in range(0, len(s) - L + 1):
            yield s[i:i + L], i, L


def main():
    print("=== VenusVaccine TumorBinary external benchmark ===")
    # 1. Load tumor antigen test set
    paths = [EXT / "VenusVaccine/TumorBinary/test.csv",
             EXT / "VenusVaccine/TumorBinary/valid.csv",
             EXT / "VenusVaccine/TumorBinary/train.csv"]
    parts = []
    for p in paths:
        if p.exists():
            d = pd.read_csv(p)
            d["split"] = p.stem
            parts.append(d)
    venus = pd.concat(parts, ignore_index=True)
    print(f"VenusVaccine tumor n={len(venus)} (cols: {venus.columns.tolist()})")
    print(f"  label dist: {venus['label'].value_counts().to_dict()}")
    print(f"  seq length: median={venus['aa_seq'].str.len().median():.0f}, "
          f"min={venus['aa_seq'].str.len().min()}, max={venus['aa_seq'].str.len().max()}")

    # Use test split for primary evaluation (78 antigens)
    test = venus[venus["split"] == "test"].reset_index(drop=True)
    valid = venus[venus["split"] == "valid"].reset_index(drop=True)
    print(f"  test n={len(test)} (pos={test['label'].sum()}), valid n={len(valid)} (pos={valid['label'].sum()})")

    # 2. Train pool same as ITSNdb
    print("\nLoading master + training source-balanced RF...")
    master = load_master_benchmark()
    train_pool = filter_loso_eligible(master,
                                      exclude_safety=("TRAINING_OVERLAP", "DEMO_ONLY", "PREDICTED_ONLY"),
                                      min_n_per_source=50,
                                      require_both_classes=True)
    train_pool = cap_per_source(train_pool, cap=5000)
    print(f"  train pool n={len(train_pool):,}")

    n_per_src = train_pool["source"].value_counts().to_dict()
    train_pool["sw"] = train_pool["source"].map(lambda s: 1.0 / n_per_src[s])
    train_pool["sw"] = train_pool["sw"] * len(train_pool) / train_pool["sw"].sum()

    encode_hla, n_hla_dim, _ = build_hla_onehot_factory(train_pool["HLA_norm"].tolist())
    Xtr = build_features(train_pool, encode_hla)
    ytr = train_pool["label"].values
    sw = train_pool["sw"].values
    clf = fit_rf(Xtr, ytr, sample_weight=sw, n_estimators=200, max_depth=10)

    # 3. For each protein in venus, enumerate 9/10/11-mers x EVAL_HLAS panel,
    #    score, aggregate.
    print("\nEnumerating peptides + scoring...")
    aggregated = []
    for split_name, dset in [("test", test), ("valid", valid)]:
        for ridx, row in dset.iterrows():
            seq = row["aa_seq"]
            label = int(row["label"])
            name = row["name"]
            # All 9-11mer windows
            wins = list(enumerate_windows(seq, lengths=(9, 10, 11)))
            if not wins:
                continue
            # Cap at 2000 windows per protein for speed (proteins ~500 aa generate ~1500 windows)
            if len(wins) > 2500:
                wins = wins[:2500]
            wseqs = [w[0] for w in wins]

            # Score against panel of common HLAs (max over alleles for each window)
            top_window_scores = []
            allele_max = {}
            for hla in EVAL_HLAS:
                rep = pd.DataFrame({"peptide": wseqs, "HLA_norm": [hla] * len(wseqs)})
                Xv = build_features(rep, encode_hla)
                proba = clf.predict_proba(Xv)[:, 1]
                allele_max[hla] = float(proba.max())
                top_window_scores.append(proba)
            # Best across HLAs for each window
            top_window_scores = np.stack(top_window_scores).max(axis=0)

            aggregated.append({
                "name": name,
                "split": split_name,
                "label": label,
                "n_windows": len(wseqs),
                "max_score": float(top_window_scores.max()),
                "mean_score": float(top_window_scores.mean()),
                "top10_mean": float(np.sort(top_window_scores)[-10:].mean() if len(top_window_scores) >= 10 else top_window_scores.mean()),
                "top1pct_mean": float(np.sort(top_window_scores)[-max(1, len(top_window_scores)//100):].mean()),
                **{f"hla_{h.split('*')[1].replace(':','_')}_max": v for h, v in allele_max.items()},
            })
        print(f"  scored {split_name}: {len(dset)} antigens")

    df_agg = pd.DataFrame(aggregated)
    df_agg.to_csv(RESULTS / "venusvaccine_scored.tsv", sep="\t", index=False)

    # 4. Aggregation-level AUROC
    rows = []
    for split_name in ["test", "valid", "test+valid"]:
        if split_name == "test+valid":
            sub = df_agg
        else:
            sub = df_agg[df_agg["split"] == split_name]
        if len(sub) < 5 or sub["label"].nunique() < 2:
            continue
        for agg in ["max_score", "mean_score", "top10_mean", "top1pct_mean"]:
            m = metrics(sub["label"].values, sub[agg].values)
            m["split"] = split_name
            m["aggregator"] = agg
            m["n_pos"] = int(sub["label"].sum())
            m["n_total"] = len(sub)
            rows.append(m)
            print(f"  {split_name:12s} agg={agg:14s} n={len(sub)} (pos={int(sub['label'].sum())}) "
                  f"AUROC={m['AUROC']:.3f} AUPRC={m['AUPRC']:.3f}")

    out = RESULTS / "venusvaccine_results.tsv"
    pd.DataFrame(rows).to_csv(out, sep="\t", index=False)
    print(f"\nsaved → {out}")
    print(f"saved → {RESULTS / 'venusvaccine_scored.tsv'}")

    composition = {
        "venusvaccine_n_total": int(len(venus)),
        "venusvaccine_test_n": int(len(test)),
        "venusvaccine_test_pos": int(test["label"].sum()),
        "n_eval_hlas": len(EVAL_HLAS),
        "lengths": [9, 10, 11],
        "aggregation": "max over HLA panel; per-protein {max, mean, top10_mean, top1pct_mean}",
        "caveat": "VenusVaccine antigens are FULL proteins, not 8-11mers — peptide-level model is being stressed at protein level. Compatibility is imperfect; per-window negatives are ambiguous (no per-window label).",
    }
    (RESULTS / "venusvaccine_composition.json").write_text(json.dumps(composition, indent=2))


if __name__ == "__main__":
    main()
