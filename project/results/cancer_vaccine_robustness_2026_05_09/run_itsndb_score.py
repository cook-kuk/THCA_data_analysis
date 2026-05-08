"""ITSNdb (Fasoulis 2024) — independent tumor-specific neoantigen DB external benchmark.

Trains source-balanced RF on master benchmark_clean.tsv (excluding TRAINING_OVERLAP),
scores ITSNdb peptide×HLA pairs, reports AUROC + overlap audit.

Output: itsndb_results.tsv
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


def load_itsndb():
    """Load main ITSNdb (200 SNV peptides) and Val_dataset (120 peptides)."""
    main_csv = EXT / "ITSNdb/data/ITSNdb.csv"
    val_csv = EXT / "ITSNdb/data/Val_dataset.csv"

    main = pd.read_csv(main_csv)
    main_clean = pd.DataFrame({
        "peptide": main["Neoantigen"].astype(str).str.upper(),
        "HLA_norm": main["HLA"].apply(normalize_hla),
        "label": (main["NeoType"].astype(str).str.lower() == "positive").astype(int),
        "subset": "ITSNdb_main",
    })

    val = pd.read_csv(val_csv)
    # Sample column has Pos/Neg as label proxy
    val_clean = pd.DataFrame({
        "peptide": val["Neoantigen"].astype(str).str.upper(),
        "HLA_norm": val["HLA"].apply(normalize_hla),
        "label": (val["Sample"].astype(str).str.lower().str.startswith("pos")).astype(int),
        "subset": "ITSNdb_Val",
    })

    df = pd.concat([main_clean, val_clean], ignore_index=True)
    df = df[df["HLA_norm"].notna() & df["peptide"].str.len().between(8, 15)].reset_index(drop=True)
    return df


def main():
    print("=== ITSNdb external benchmark ===")
    # 1. Load ITSNdb
    itsn = load_itsndb()
    print(f"ITSNdb loaded: n={len(itsn)} (main={ (itsn['subset']=='ITSNdb_main').sum()}, "
          f"Val={(itsn['subset']=='ITSNdb_Val').sum()})")
    print(f"  label distribution: {itsn['label'].value_counts().to_dict()}")
    print(f"  unique HLAs: {itsn['HLA_norm'].nunique()}")

    # 2. Load master, filter, train source-balanced RF
    print("\nLoading master benchmark...")
    master = load_master_benchmark()
    print(f"  master n={len(master):,}")

    # Overlap audit (peptide × HLA)
    master_keys = set(zip(master["peptide"], master["HLA_norm"]))
    itsn["in_master"] = [(p, h) in master_keys for p, h in zip(itsn["peptide"], itsn["HLA_norm"])]
    n_overlap = itsn["in_master"].sum()
    print(f"\n  ITSNdb peptide×HLA overlap with master: {n_overlap}/{len(itsn)}")
    if n_overlap > 0:
        print("  Flagging overlapping rows (will exclude from external AUROC)")

    # 3. Train pool: master excluding TRAINING_OVERLAP / DEMO_ONLY / PREDICTED_ONLY
    train_pool = filter_loso_eligible(master,
                                      exclude_safety=("TRAINING_OVERLAP", "DEMO_ONLY", "PREDICTED_ONLY"),
                                      min_n_per_source=50,
                                      require_both_classes=True)
    print(f"  train pool after LOSO-eligible filter: n={len(train_pool):,}")
    print(f"  sources in train pool: {train_pool['source'].value_counts().to_dict()}")

    # Cap each source at 5000 (CEDAR + NEPdb dominate)
    train_pool = cap_per_source(train_pool, cap=5000)
    print(f"  after per-source cap=5000: n={len(train_pool):,}")

    # 4. Source-balanced sample weights: each source contributes equally
    n_per_src = train_pool["source"].value_counts().to_dict()
    train_pool["sw"] = train_pool["source"].map(lambda s: 1.0 / n_per_src[s])
    train_pool["sw"] = train_pool["sw"] * len(train_pool) / train_pool["sw"].sum()  # normalize

    # 5. Build features (biophys + HLA one-hot, master HLA vocabulary)
    encode_hla, n_hla_dim, _ = build_hla_onehot_factory(train_pool["HLA_norm"].tolist())
    Xtr = build_features(train_pool, encode_hla)
    ytr = train_pool["label"].values
    sw = train_pool["sw"].values

    print(f"\n  Xtr shape: {Xtr.shape}, ytr pos={ytr.sum()}/{len(ytr)} ({ytr.mean():.3f})")
    print("  fitting source-balanced RF (n=200 trees)...")
    clf = fit_rf(Xtr, ytr, sample_weight=sw, n_estimators=200, max_depth=10)

    # 6. Score ITSNdb
    Xte = build_features(itsn.assign(HLA_norm=itsn["HLA_norm"]), encode_hla)
    yte = itsn["label"].values
    proba = clf.predict_proba(Xte)[:, 1]
    itsn["pred_proba"] = proba

    # 7. Per-subset metrics
    rows = []
    for sub in ["ITSNdb_main", "ITSNdb_Val", "ITSNdb_combined", "ITSNdb_no_overlap"]:
        if sub == "ITSNdb_combined":
            mask = np.ones(len(itsn), dtype=bool)
        elif sub == "ITSNdb_no_overlap":
            mask = ~itsn["in_master"].values
        else:
            mask = itsn["subset"].values == sub
        if mask.sum() < 5 or len(set(yte[mask])) < 2:
            print(f"  {sub}: skipping (n={mask.sum()}, classes={set(yte[mask])})")
            continue
        m = metrics(yte[mask], proba[mask])
        m["subset"] = sub
        m["n_pos"] = int(yte[mask].sum())
        m["n_total"] = int(mask.sum())
        rows.append(m)
        print(f"  {sub}: n={mask.sum()} (pos={int(yte[mask].sum())}) AUROC={m['AUROC']:.3f} AUPRC={m['AUPRC']:.3f}")

    out = RESULTS / "itsndb_results.tsv"
    pd.DataFrame(rows).to_csv(out, sep="\t", index=False)
    print(f"\nsaved → {out}")

    # Save scored peptides (for downstream inspection)
    itsn.to_csv(RESULTS / "itsndb_scored.tsv", sep="\t", index=False)
    print(f"saved → {RESULTS / 'itsndb_scored.tsv'}")

    # Save train-pool composition for transparency
    composition = {
        "train_pool_n": int(len(train_pool)),
        "train_pool_per_source": train_pool["source"].value_counts().to_dict(),
        "train_pool_pos_rate": float(ytr.mean()),
        "feature_dim": int(Xtr.shape[1]),
        "model": "RandomForestClassifier(n=200,depth=10,balanced)",
        "sample_weight": "1/n_source normalized",
        "itsndb_n": int(len(itsn)),
        "itsndb_pos_rate": float(yte.mean()),
        "n_overlap_with_master": int(n_overlap),
    }
    (RESULTS / "itsndb_train_composition.json").write_text(json.dumps(composition, indent=2))
    print(f"saved → {RESULTS / 'itsndb_train_composition.json'}")


if __name__ == "__main__":
    main()
