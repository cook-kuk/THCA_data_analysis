"""Wave 5A — Step 1: Build distillation pool.

Sources combined into `distill_pool.tsv`:
  1. Master TSV (benchmark_clean.tsv at /data/neoantigen_vaccine_hub/...) —
     all rows with HLA normalization OK and 8-11mer peptide. Primary IEDB+CEDAR+
     NEPdb+TESLA+Neodb+... ~27.7k pep×HLA pairs.
  2. Wave 1 train pool (n=2396) is already a subset of (1) — kept implicitly.
  3. ITSNdb peptides EXCLUDED (we evaluate on this; would leak).

Output columns:
  peptide | hla | source | label_hard (or NaN) | in_master_train (bool flag)

Filter:
  - peptide length 8..11 (MHCflurry domain)
  - HLA in MHCflurry-supported list
  - dedupe (peptide, hla) — keep first occurrence (priority by source order)
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, "/home/seungho/personal/THCA_data_analysis/project/results/cancer_vaccine_robustness_2026_05_09")
from _common import load_master_benchmark, normalize_hla

OUT = Path("/home/seungho/personal/THCA_data_analysis/project/results/p_neo_bayesian_2026_05_09/wave5a")
BUNDLE = Path("/home/seungho/personal/THCA_data_analysis/project/results/p_neo_bayesian_2026_05_09/bundle.tsv")
ITSN_EXT = Path("/home/seungho/personal/THCA_data_analysis/project/data/external_benchmarks/ITSNdb/data")


def load_itsndb_peptides():
    main = pd.read_csv(ITSN_EXT / "ITSNdb.csv")
    val = pd.read_csv(ITSN_EXT / "Val_dataset.csv")
    peps = set(main["Neoantigen"].astype(str).str.upper().tolist())
    peps |= set(val["Neoantigen"].astype(str).str.upper().tolist())
    return peps


def get_supported_alleles():
    from mhcflurry import Class1PresentationPredictor
    predictor = Class1PresentationPredictor.load()
    return set(predictor.supported_alleles)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    print("=== Wave 5A — Build distill pool ===")

    print("\n[1/4] Loading master benchmark...")
    master = load_master_benchmark()
    print(f"  master raw: {len(master):,}")

    print("\n[2/4] Filter: 8-11mer + HLA-norm OK")
    master = master[master["HLA_norm"].notna()]
    master = master[master["peptide"].str.len().between(8, 11)].copy()
    print(f"  after length+HLA filter: {len(master):,}")

    print("\n[3/4] Filter: MHCflurry-supported alleles")
    supported = get_supported_alleles()
    master = master[master["HLA_norm"].isin(supported)].copy()
    print(f"  after MHCflurry-allele filter: {len(master):,}")

    print("\n[4/4] Exclude ITSNdb peptides (eval set)")
    itsn_peps = load_itsndb_peptides()
    master = master[~master["peptide"].isin(itsn_peps)].copy()
    print(f"  after ITSNdb exclusion: {len(master):,}")

    # Dedupe peptide × HLA (keep first; sort first so deterministic)
    master = master.sort_values(["source", "peptide", "HLA_norm"]).reset_index(drop=True)
    before = len(master)
    master = master.drop_duplicates(subset=["peptide", "HLA_norm"], keep="first").reset_index(drop=True)
    print(f"  after pep×HLA dedupe: {len(master):,} (removed {before-len(master):,})")

    # Mark which rows are in Wave 1 train pool (n≈2396) for hard-label use
    bundle = pd.read_csv(BUNDLE, sep="\t")
    train = bundle[bundle["split"] == "train"].copy()
    train_keys = set(zip(train["peptide"].astype(str).str.upper(),
                         train["HLA_norm"].astype(str)))
    master["in_train_pool"] = [
        (p, h) in train_keys for p, h in zip(master["peptide"], master["HLA_norm"])
    ]
    print(f"\n  rows flagged as Wave 1 train pool: {master['in_train_pool'].sum():,}")
    # NB: master rows have hard labels (the original `label` col). For rows in
    # train_pool we'll keep them; for rows OUTSIDE train_pool we still have
    # labels but they may be lower-quality / weakly-positive — we'll still
    # surface them but mark which ones we trust for HARD-label loss.

    # Output
    out = pd.DataFrame({
        "peptide": master["peptide"].astype(str),
        "hla": master["HLA_norm"].astype(str),
        "source": master["source"].astype(str),
        "label_hard": master["label"].astype(int),
        "in_train_pool": master["in_train_pool"].astype(bool),
        "safety": master.get("safety", pd.Series([""]*len(master))).astype(str),
    })
    out_tsv = OUT / "distill_pool.tsv"
    out.to_csv(out_tsv, sep="\t", index=False)
    print(f"\nsaved: {out_tsv}  rows={len(out):,}")

    # Stats
    print("\n--- pool stats ---")
    print(f"  unique peptides: {out['peptide'].nunique():,}")
    print(f"  unique HLAs:     {out['hla'].nunique()}")
    print(f"  unique pep×HLA:  {out[['peptide','hla']].drop_duplicates().shape[0]:,}")
    print(f"  source breakdown:")
    for s, n in out["source"].value_counts().items():
        print(f"    {s:30s} {n:>7,}")
    print(f"  label hard pos rate: {out['label_hard'].mean():.3f}")
    print(f"  train_pool rows:     {out['in_train_pool'].sum():,}")


if __name__ == "__main__":
    main()
