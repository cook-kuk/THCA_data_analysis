#!/usr/bin/env python3
"""Compute leakage-stratified AUROC for each algorithm × test set.

Outputs:
  algorithm_sweep_results.tsv — long-form table
  algorithm_sweep_pivot.tsv  — wide pivot for forest plot
"""
import json
import numpy as np
import pandas as pd
from pathlib import Path
from sklearn.metrics import roc_auc_score

ROOT = Path(__file__).resolve().parent
RNG = np.random.default_rng(42)
N_BOOT = 1000

scoring = pd.read_csv(ROOT / "scoring_set_itsndb.tsv", sep="\t")
print(f"scoring set n={len(scoring)}")

# Load each algorithm's scored output
mhcf = pd.read_csv(ROOT / "mhcflurry_scored.tsv", sep="\t")
bigmhc = pd.read_csv(ROOT / "bigmhc_scored.tsv", sep="\t")
bigmhc_el_raw = pd.read_csv(ROOT / "bigmhc_el_output.csv")
prime = pd.read_csv(ROOT / "prime_scored.tsv", sep="\t") if (ROOT / "prime_scored.tsv").exists() else None

# Merge BigMHC EL into a common dataframe (it shares peptide+HLA)
bigmhc_el_renamed = bigmhc_el_raw.rename(columns={"mhc": "HLA_norm", "pep": "peptide"})[
    ["peptide", "HLA_norm", "BigMHC_EL"]
]

# Build a master per-peptide score table
master = scoring.copy()
master = master.merge(
    mhcf[["peptide", "HLA_norm", "mhcflurry_presentation", "mhcflurry_affinity"]],
    on=["peptide", "HLA_norm"], how="left"
)
master = master.merge(
    bigmhc[["peptide", "HLA_norm", "bigmhc_BigMHC_IM"]],
    on=["peptide", "HLA_norm"], how="left"
).rename(columns={"bigmhc_BigMHC_IM": "bigmhc_im"})
master = master.merge(bigmhc_el_renamed, on=["peptide", "HLA_norm"], how="left")
master = master.rename(columns={"BigMHC_EL": "bigmhc_el"})
if prime is not None:
    master = master.merge(
        prime[["peptide", "HLA_norm", "prime_score", "prime_rank"]],
        on=["peptide", "HLA_norm"], how="left"
    )

master.to_csv(ROOT / "scored_master.tsv", sep="\t", index=False)
print(f"master shape: {master.shape}, columns: {list(master.columns)}")

# Bring in our own ESM2-Bayesian + biophys+RF predictions if available
EXT_DIR = ROOT.parent
our_pred = EXT_DIR / "predictions_itsndb.tsv"
ours = None
if our_pred.exists():
    ours = pd.read_csv(our_pred, sep="\t")
    print(f"loaded our wave1 predictions: {ours.shape}, columns: {list(ours.columns)}")
    keep = [c for c in ["peptide", "HLA_norm", "pred_mean", "pred_std"] if c in ours.columns]
    master = master.merge(ours[keep], on=["peptide", "HLA_norm"], how="left", suffixes=("", "_w1"))

# Wave 2 structure_lr (biophys/structure + Wave1 ensemble) predictions
w2_pred = EXT_DIR / "wave2" / "structure_lr_predictions.tsv"
if w2_pred.exists():
    w2 = pd.read_csv(w2_pred, sep="\t")
    print(f"loaded wave2 structure_lr: {w2.shape}, columns: {list(w2.columns)}")
    keep = [c for c in ["peptide", "HLA_norm", "p_struct_lr", "p_wave1_plus_struct_avg50"] if c in w2.columns]
    master = master.merge(w2[keep], on=["peptide", "HLA_norm"], how="left")

# Wave 2 VQC predictions (quantum, exploratory)
vqc_pred = EXT_DIR / "wave2" / "vqc_predictions.tsv"
if vqc_pred.exists():
    v = pd.read_csv(vqc_pred, sep="\t")
    keep = [c for c in ["peptide", "HLA_norm", "p_vqc", "p_lr_8d_classical"] if c in v.columns]
    master = master.merge(v[keep], on=["peptide", "HLA_norm"], how="left")

# Algorithm spec: name -> (column, higher_is_positive_for_label_1)
# label 1 = immunogenic. presentation_score / IM / score: high=positive.
ALG_SPECS = {
    "MHCflurry_presentation": ("mhcflurry_presentation", True),
    "MHCflurry_affinity_neg": ("mhcflurry_affinity", False),  # affinity nM lower=binder; flip
    "BigMHC_IM":              ("bigmhc_im", True),
    "BigMHC_EL":              ("bigmhc_el", True),
}
if "prime_score" in master.columns:
    ALG_SPECS["PRIME_score"] = ("prime_score", True)
    ALG_SPECS["PRIME_rank_neg"] = ("prime_rank", False)

# Add our models if columns are present
for col, name, sign_pos in [
    ("pred_mean",                  "Ours_ESM2_Bayesian",      True),
    ("p_struct_lr",                "Ours_BiophysStructLR",    True),
    ("p_wave1_plus_struct_avg50",  "Ours_Wave1_plus_Struct",  True),
    ("p_vqc",                      "Ours_VQC",                True),
    ("p_lr_8d_classical",          "Ours_LR8d",               True),
]:
    if col in master.columns:
        ALG_SPECS[name] = (col, sign_pos)

print(f"\nAlgorithms to evaluate: {list(ALG_SPECS.keys())}")

def boot_auroc(y, s, n=N_BOOT, rng=RNG):
    aurocs = []
    n_samp = len(y)
    if len(np.unique(y)) < 2:
        return np.nan, np.nan, np.nan
    base = roc_auc_score(y, s)
    for _ in range(n):
        idx = rng.integers(0, n_samp, n_samp)
        ys, ss = y[idx], s[idx]
        if len(np.unique(ys)) < 2:
            continue
        aurocs.append(roc_auc_score(ys, ss))
    if not aurocs:
        return base, np.nan, np.nan
    return base, float(np.percentile(aurocs, 2.5)), float(np.percentile(aurocs, 97.5))

results = []
for alg, (col, sign_pos) in ALG_SPECS.items():
    if col not in master.columns:
        continue
    for testset in [("ext_itsndb_main",), ("ext_itsndb_val",), ("ext_itsndb_main", "ext_itsndb_val")]:
        ts_name = "+".join(testset)
        sub_all = master[master["split"].isin(testset)]
        for in_master_val in [True, False, None]:
            if in_master_val is None:
                sub = sub_all
                im_lbl = "all"
            else:
                sub = sub_all[sub_all["in_master"] == in_master_val]
                im_lbl = str(in_master_val)
            valid = sub[col].notna()
            sub_v = sub[valid]
            n = len(sub_v)
            n_pos = int((sub_v["label"] == 1).sum())
            n_neg = n - n_pos
            n_skipped = int((~valid).sum())
            if n < 10 or n_pos < 1 or n_neg < 1:
                results.append({
                    "algorithm": alg, "testset": ts_name, "in_master": im_lbl,
                    "n": n, "n_pos": n_pos, "n_neg": n_neg,
                    "AUROC": np.nan, "AUROC_lo95": np.nan, "AUROC_hi95": np.nan,
                    "n_skipped_unsupported_hla": n_skipped,
                    "note": "insufficient class balance" if n>0 else "no data"
                })
                continue
            y = sub_v["label"].values.astype(int)
            s = sub_v[col].values.astype(float)
            if not sign_pos:
                s = -s
            auc, lo, hi = boot_auroc(y, s)
            results.append({
                "algorithm": alg, "testset": ts_name, "in_master": im_lbl,
                "n": n, "n_pos": n_pos, "n_neg": n_neg,
                "AUROC": auc, "AUROC_lo95": lo, "AUROC_hi95": hi,
                "n_skipped_unsupported_hla": n_skipped,
                "note": ""
            })

res = pd.DataFrame(results)
res.to_csv(ROOT / "algorithm_sweep_results.tsv", sep="\t", index=False)
print(f"\nsaved {ROOT / 'algorithm_sweep_results.tsv'}")
print(res.to_string(index=False))

# Pivot for forest plot: rows = (algorithm), cols = (testset, in_master)
pivot = res[res["testset"] == "ext_itsndb_main+ext_itsndb_val"].pivot_table(
    index="algorithm", columns="in_master", values=["AUROC", "AUROC_lo95", "AUROC_hi95"]
)
pivot.to_csv(ROOT / "algorithm_sweep_pivot.tsv", sep="\t")
print("\nPivot (combined ITSNdb):")
print(pivot)
