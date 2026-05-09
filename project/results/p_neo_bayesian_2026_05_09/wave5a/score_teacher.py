"""Wave 5A — Step 2: Score distillation pool with MHCflurry teacher.

Outputs `distill_teacher_scores.tsv`:
  peptide | hla | mhcflurry_presentation | mhcflurry_processing | mhcflurry_affinity

Notes:
  - We use Class1PresentationPredictor.predict() with sample_names → 1:1 alignment
    of peptides with alleles.
  - affinity is reported in nM; we ALSO compute a transformed affinity_score in
    [0,1] consistent with MHCflurry's convention: 1 - log(aff)/log(50000), clipped.
  - processing_score in [0,1].
  - presentation_score in [0,1] (target soft label for distillation).
"""
from __future__ import annotations
import sys, time
from pathlib import Path
import numpy as np
import pandas as pd

OUT = Path("/home/seungho/personal/THCA_data_analysis/project/results/p_neo_bayesian_2026_05_09/wave5a")
POOL_TSV = OUT / "distill_pool.tsv"


def affinity_to_score(aff_nm):
    """Map affinity (nM) to score in [0,1] following MHCflurry convention."""
    a = np.asarray(aff_nm, dtype=np.float64)
    a = np.clip(a, 1.0, 50000.0)
    return np.clip(1.0 - np.log(a) / np.log(50000.0), 0.0, 1.0)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    print("=== Wave 5A — Score distill pool with MHCflurry ===")

    pool = pd.read_csv(POOL_TSV, sep="\t")
    print(f"pool: {len(pool):,} pairs")

    from mhcflurry import Class1PresentationPredictor
    predictor = Class1PresentationPredictor.load()

    # Group by allele — fastest path is single-allele predict (no per-row sample
    # overhead). Each allele call internally batches the peptides.
    results = []
    t0 = time.time()
    n_alleles = pool["hla"].nunique()
    for i, (allele, g) in enumerate(pool.groupby("hla", sort=True)):
        peps = g["peptide"].tolist()
        # Dedupe peptides per allele for predict speed; merge back below
        u_peps = list(dict.fromkeys(peps))
        try:
            df = predictor.predict(peptides=u_peps, alleles=[allele], verbose=0)
        except Exception as e:
            print(f"  [warn] {allele} predict failed: {e}")
            continue
        df = df[["peptide", "affinity", "processing_score",
                 "presentation_score"]].copy()
        df["hla"] = allele
        results.append(df[["peptide", "hla", "presentation_score",
                           "processing_score", "affinity"]])
        if (i + 1) % 10 == 0 or i == 0 or i == n_alleles - 1:
            elapsed = time.time() - t0
            done = sum(len(r) for r in results)
            print(f"  [{i+1}/{n_alleles}] {allele:>15s} "
                  f"n_peps={len(u_peps):>5,} cum_done={done:,}  elapsed={elapsed:.1f}s",
                  flush=True)

    teacher = pd.concat(results, ignore_index=True)
    # add affinity_score (transformed)
    teacher["affinity_score"] = affinity_to_score(teacher["affinity"].values)

    # rename to consistent column names
    teacher = teacher.rename(columns={
        "presentation_score": "mhcflurry_presentation",
        "processing_score": "mhcflurry_processing",
        "affinity": "mhcflurry_affinity",
        "affinity_score": "mhcflurry_aff_score",
    })

    # Merge with pool (some rows may have failed → NaN)
    merged = pool.merge(teacher, on=["peptide", "hla"], how="left")
    print(f"\n  merged: {len(merged):,}  (NaN presentation: "
          f"{merged['mhcflurry_presentation'].isna().sum():,})")

    out_tsv = OUT / "distill_teacher_scores.tsv"
    merged.to_csv(out_tsv, sep="\t", index=False)
    print(f"saved: {out_tsv}")

    # Sanity: teacher AUROC on labeled rows in pool (label_hard != NaN, both classes)
    have_pres = merged["mhcflurry_presentation"].notna()
    sub = merged[have_pres & merged["label_hard"].notna()].copy()
    sub["label_hard"] = sub["label_hard"].astype(int)
    if sub["label_hard"].nunique() == 2 and len(sub) > 100:
        from sklearn.metrics import roc_auc_score
        auc_pres = roc_auc_score(sub["label_hard"].values,
                                 sub["mhcflurry_presentation"].values)
        auc_proc = roc_auc_score(sub["label_hard"].values,
                                 sub["mhcflurry_processing"].values)
        auc_aff = roc_auc_score(sub["label_hard"].values,
                                sub["mhcflurry_aff_score"].values)
        print(f"\n  teacher sanity AUROC on pool's hard labels (n={len(sub):,}, "
              f"pos_rate={sub['label_hard'].mean():.3f}):")
        print(f"    presentation_score : {auc_pres:.3f}")
        print(f"    processing_score   : {auc_proc:.3f}")
        print(f"    aff_score          : {auc_aff:.3f}")

    print(f"\ntotal time: {(time.time()-t0):.1f}s")


if __name__ == "__main__":
    main()
