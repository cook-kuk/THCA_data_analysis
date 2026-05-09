#!/usr/bin/env python3
"""Wave 11 — import wave3 ITSNdb predictions as wave11 test bundles.

Existing predictions (already scored) live in:
  wave3_mhcflurry/predictions.tsv      score_mhcflurry
  wave3_bigmhc/predictions.tsv         score_bigmhc_im
  wave3_deepimmuno/predictions.tsv     score_deepimmuno
  wave3_transphla/predictions.tsv      score_transphla
  wave3_netmhcpan/predictions.tsv      score_netmhcpan
  wave3_prime/predictions.tsv          score_prime

We re-emit them as wave11_predictions/<algo>__itsndb.tsv with the wave11 schema:
  peptide | hla | label | source | length | n_overlap_with_other_sources | <score_col>

Source split: ITSNdb_main + ITSNdb_val combined into one bundle "itsndb".
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path("/home/seungho/personal/THCA_data_analysis/project/results/p_neo_bayesian_2026_05_09")
WAVE11 = ROOT / "wave11"
OUT = WAVE11 / "wave11_predictions"

WAVE3_FILES = {
    "MHCflurry": ("wave3_mhcflurry/predictions.tsv", "score_mhcflurry"),
    "BigMHC_IM": ("wave3_bigmhc/predictions.tsv", "score_bigmhc_im"),
    "DeepImmuno": ("wave3_deepimmuno/predictions.tsv", "score_deepimmuno"),
    "TransPHLA": ("wave3_transphla/predictions.tsv", "score_transphla"),
    "NetMHCpan_4.1": ("wave3_netmhcpan/predictions.tsv", "score_netmhcpan"),
    "PRIME": ("wave3_prime/predictions.tsv", "score_prime"),
}

# Also add an ITSNdb test bundle to test_bundles for record-keeping
itsndb_bundle = None
for algo, (fp, score_col) in WAVE3_FILES.items():
    src_path = ROOT / fp
    if not src_path.exists():
        print(f"[skip] {algo}: {fp} missing"); continue
    df = pd.read_csv(src_path, sep="\t")
    df = df.dropna(subset=[score_col])
    if len(df) == 0: continue
    df["length"] = df["peptide"].str.len()
    df["source"] = "ITSNdb"
    df["n_overlap_with_other_sources"] = 0  # ITSNdb is independent
    out = df[["peptide", "hla", "label", "source", "length",
              "n_overlap_with_other_sources", score_col]].copy()
    out.to_csv(OUT / f"{algo}__itsndb.tsv", sep="\t", index=False)
    print(f"[wave3->wave11] {algo}: {len(out)} ITSNdb preds emitted", flush=True)
    if itsndb_bundle is None:
        # Save a representative bundle file (use peptide+hla+label, no scores)
        itsndb_bundle = df[["peptide", "hla", "label", "source", "length",
                            "n_overlap_with_other_sources"]].drop_duplicates(["peptide", "hla"])
        itsndb_bundle.to_csv(WAVE11 / "test_bundles" / "itsndb.tsv", sep="\t", index=False)
        print(f"[bundle] saved itsndb.tsv with {len(itsndb_bundle)} rows", flush=True)

# Also import our wave1+wave2 predictions if present (Bayesian + Structure_LR + GP_quantum)
WAVE12_FILES = {
    "ESM2_Bayesian": ("wave1/predictions_itsndb.tsv", None),  # auto-detect
    "Structure_LR":  ("wave2/structure_lr_predictions.tsv", "p_struct_lr"),
    "GP_quantum":    ("wave2/qk_svm_predictions.tsv", None),
    "VQC":           ("wave2/vqc_predictions.tsv", None),
}

for algo, (fp, sc) in WAVE12_FILES.items():
    src_path = ROOT / fp
    if not src_path.exists(): continue
    df = pd.read_csv(src_path, sep="\t")
    if "label" not in df.columns: continue
    if sc is None:
        # Find first numeric col not in {peptide,hla,label,...}
        cands = [c for c in df.columns if c.startswith(("p_", "pred_", "prob_", "score"))
                 and df[c].dtype in [np.float32, np.float64]]
        if not cands: continue
        sc = cands[0]
    if sc not in df.columns: continue
    if "HLA_norm" in df.columns:
        df = df.rename(columns={"HLA_norm": "hla"})
    df = df.dropna(subset=[sc])
    if len(df) == 0: continue
    df["length"] = df["peptide"].str.len()
    df["source"] = "ITSNdb"
    df["n_overlap_with_other_sources"] = 0
    out = df[["peptide", "hla", "label", "source", "length",
              "n_overlap_with_other_sources", sc]].rename(columns={sc: f"score_{algo.lower()}"})
    out.to_csv(OUT / f"{algo}__itsndb.tsv", sep="\t", index=False)
    print(f"[wave1/2->wave11] {algo}: {len(out)} ITSNdb preds emitted from {fp} (col={sc})", flush=True)

print("[done] wave3+wave1/2 import complete")
