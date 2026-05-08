#!/usr/bin/env python3
"""
Decisive split-stress retraining for the RAS-like AUC=1.000 audit.

Strategy A — Multi-seed shuffled KFold(5)  (re-run original recipe with
            different random_state to test if seed=42 was lucky)
Strategy B — StratifiedKFold(3) on (label × molecular_subtype) so each
            fold contains ~1 RAS-like positive + ~4 RAS-like negatives
Strategy C — Leave-one-out within the 16 RAS-like slides
            (each held out once, trained on 58 others)

Outputs to analysis_supp/audit_ras_auc100/:
  retrain_strategyA_multi_seed.tsv   (overall + RAS-like AUC × seed)
  retrain_strategyB_stratified.tsv
  retrain_strategyC_loo_RASlike.tsv
  RETRAIN_SUMMARY.json
"""
from __future__ import annotations
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import KFold, StratifiedKFold

import sys
sys.path.insert(0, str(Path(__file__).parent))
from audit_clam_train_lib import train_one_fold, load_features  # noqa: E402

ROOT = Path("/home/seungho/personal/THCA_data_analysis/project/results/"
            "p2_image_dm1_v2_foundation_clam_2026_05_07")
META = Path("/home/seungho/personal/THCA_data_analysis/project/metadata/"
            "sample_master_v3.tsv")
FEAT_DIR = ROOT / "phase2_tcga_clam" / "features"
OUT = ROOT / "analysis_supp" / "audit_ras_auc100"
OUT.mkdir(parents=True, exist_ok=True)

DEVICE = "cuda" if torch.cuda.is_available() else "cpu"
EPOCHS = 30


def load_dataset():
    man = pd.read_csv(ROOT / "phase2_tcga_clam" / "slide_manifest.tsv", sep="\t")
    man = man[man["dm"].isin(["DM1", "DM2"])].copy()
    man["label"] = man["dm"].map({"DM2": 0, "DM1": 1})
    # only keep slides whose features exist
    man["has_feat"] = man["file_id"].apply(
        lambda f: (FEAT_DIR / f"{f}.pt").exists())
    man = man[man["has_feat"]].reset_index(drop=True)

    # subgroup labels
    meta = pd.read_csv(META, sep="\t", low_memory=False)
    meta = meta[meta["dataset"] == "TCGA-THCA"].copy()
    meta["case_short"] = meta["sample_id"].str.extract(r"^(TCGA-[^-]+-[^-]+)")
    meta_case = (meta.sort_values("sample_id")
                     .drop_duplicates("case_short", keep="first")
                     [["case_short", "histology_subtype", "molecular_subtype"]])
    man = man.merge(meta_case, left_on="submitter_id",
                    right_on="case_short", how="left")
    print(f"[load] {len(man)} slides with features and DM1/DM2 label")
    print(man["molecular_subtype"].value_counts(dropna=False).to_string())

    print(f"[load] loading 1024-dim UNI features ...")
    bags = load_features(FEAT_DIR, man["file_id"].tolist())
    labels = man["label"].values
    return man, bags, labels


def compute_subgroup_auc(man: pd.DataFrame, oof_prob: np.ndarray,
                         group_col: str, group_val) -> dict:
    mask = man[group_col].values == group_val
    y = man["label"].values[mask]
    p = oof_prob[mask]
    if len(set(y)) < 2:
        return {"n": int(mask.sum()), "auc": np.nan, "n_pos": int(y.sum()),
                "n_neg": int(len(y) - y.sum())}
    return {"n": int(mask.sum()), "auc": float(roc_auc_score(y, p)),
            "n_pos": int(y.sum()), "n_neg": int(len(y) - y.sum())}


def strategy_a_multi_seed(man, bags, labels, seeds=(42, 1, 7, 13, 99, 2026)):
    """Original KFold(5, shuffle, random_state=seed) with multiple seeds."""
    rows = []
    for seed in seeds:
        t0 = time.time()
        kf = KFold(n_splits=5, shuffle=True, random_state=seed)
        oof = np.full(len(labels), np.nan)
        for fold_i, (tr, va) in enumerate(kf.split(bags)):
            tr_bags = [bags[i] for i in tr]
            tr_lab = [int(labels[i]) for i in tr]
            va_bags = [bags[i] for i in va]
            va_lab = [int(labels[i]) for i in va]
            _, probs = train_one_fold(tr_bags, tr_lab, va_bags, va_lab,
                                       epochs=EPOCHS, device=DEVICE, seed=seed)
            for j, idx in enumerate(va):
                oof[idx] = probs[j]
        overall = roc_auc_score(labels, oof)
        ras = compute_subgroup_auc(man, oof, "molecular_subtype", "RAS_like")
        braf = compute_subgroup_auc(man, oof, "molecular_subtype", "BRAF_like")
        elapsed = time.time() - t0
        rows.append({
            "strategy": "A_multi_seed",
            "seed": seed,
            "overall_auc": float(overall),
            "ras_like_auc": ras["auc"],
            "ras_like_n": ras["n"], "ras_like_n_pos": ras["n_pos"],
            "braf_like_auc": braf["auc"],
            "braf_like_n": braf["n"], "braf_like_n_pos": braf["n_pos"],
            "seconds": round(elapsed, 1),
        })
        print(f"  seed={seed}: overall={overall:.3f}  RAS={ras['auc']:.3f}  "
              f"BRAF={braf['auc']:.3f}  ({elapsed:.0f}s)")
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "retrain_strategyA_multi_seed.tsv", sep="\t", index=False)
    return df


def strategy_b_stratified(man, bags, labels):
    """StratifiedKFold on (label × molecular_subtype) — n_splits=3
    because (DM1, RAS_like) has only n=3 samples."""
    strat_key = (man["label"].astype(str) + "_"
                 + man["molecular_subtype"].astype(str)).values
    print(f"  strat_key counts: "
          f"{pd.Series(strat_key).value_counts().to_dict()}")
    skf = StratifiedKFold(n_splits=3, shuffle=True, random_state=42)
    oof = np.full(len(labels), np.nan)
    fold_aucs = []
    for fold_i, (tr, va) in enumerate(skf.split(bags, strat_key)):
        tr_bags = [bags[i] for i in tr]
        tr_lab = [int(labels[i]) for i in tr]
        va_bags = [bags[i] for i in va]
        va_lab = [int(labels[i]) for i in va]
        auc, probs = train_one_fold(tr_bags, tr_lab, va_bags, va_lab,
                                     epochs=EPOCHS, device=DEVICE, seed=42)
        for j, idx in enumerate(va):
            oof[idx] = probs[j]
        n_ras_pos = int(((man.iloc[va]["molecular_subtype"] == "RAS_like")
                         & (np.array(va_lab) == 1)).sum())
        n_ras_neg = int(((man.iloc[va]["molecular_subtype"] == "RAS_like")
                         & (np.array(va_lab) == 0)).sum())
        fold_aucs.append({"fold": fold_i + 1, "val_auc": auc,
                          "n_val": len(va), "n_ras_pos_val": n_ras_pos,
                          "n_ras_neg_val": n_ras_neg})
        print(f"  fold {fold_i+1}: val_auc={auc:.3f}  "
              f"RAS pos={n_ras_pos} neg={n_ras_neg} val_n={len(va)}")
    overall = roc_auc_score(labels, oof)
    ras = compute_subgroup_auc(man, oof, "molecular_subtype", "RAS_like")
    braf = compute_subgroup_auc(man, oof, "molecular_subtype", "BRAF_like")
    out = {
        "strategy": "B_stratified_n3",
        "overall_auc": float(overall),
        "ras_like_auc": ras["auc"],
        "braf_like_auc": braf["auc"],
        "fold_breakdown": fold_aucs,
    }
    pd.DataFrame([out]).to_csv(OUT / "retrain_strategyB_stratified.tsv",
                                sep="\t", index=False)
    print(f"  STRATIFIED: overall={overall:.3f}  RAS={ras['auc']:.3f}  "
          f"BRAF={braf['auc']:.3f}")
    return out


def strategy_c_loo_raslike(man, bags, labels):
    """For each of the 16 RAS-like slides, retrain on the other 58 slides
    and predict the held-out one. The other RAS-like slides are still in
    training, but the *same* slide is never seen. Then compute AUC over
    the 16 stringent OOF predictions."""
    ras_idx = np.where(man["molecular_subtype"].values == "RAS_like")[0]
    print(f"  LOO over {len(ras_idx)} RAS-like slides")
    preds = []
    for k, hold_idx in enumerate(ras_idx):
        tr = np.array([i for i in range(len(bags)) if i != hold_idx])
        tr_bags = [bags[i] for i in tr]
        tr_lab = [int(labels[i]) for i in tr]
        va_bags = [bags[hold_idx]]
        va_lab = [int(labels[hold_idx])]
        # AUC undefined for single sample → take final-epoch prob
        _, probs = train_one_fold(tr_bags, tr_lab, va_bags, va_lab,
                                   epochs=EPOCHS, device=DEVICE, seed=42)
        preds.append({
            "slide": man.iloc[hold_idx]["file_id"],
            "submitter_id": man.iloc[hold_idx]["submitter_id"],
            "label": va_lab[0],
            "prob_DM1_loo": float(probs[0]),
            "histology_subtype": man.iloc[hold_idx]["histology_subtype"],
        })
        print(f"  [{k+1}/{len(ras_idx)}] {preds[-1]['slide'][:8]}  "
              f"label={va_lab[0]}  prob={probs[0]:.3f}")
    df = pd.DataFrame(preds)
    df.to_csv(OUT / "retrain_strategyC_loo_RASlike.tsv",
              sep="\t", index=False)
    if df["label"].nunique() > 1:
        loo_auc = float(roc_auc_score(df["label"], df["prob_DM1_loo"]))
    else:
        loo_auc = np.nan
    out = {"strategy": "C_loo_raslike", "n": len(df),
           "n_pos": int(df["label"].sum()),
           "n_neg": int(len(df) - df["label"].sum()),
           "loo_auc": loo_auc}
    print(f"  LOO RAS-like AUC = {loo_auc:.3f}")
    return out


def main():
    man, bags, labels = load_dataset()
    summary = {}
    print("\n[Strategy A] Multi-seed KFold(5) — testing seed sensitivity")
    df_a = strategy_a_multi_seed(man, bags, labels)
    summary["strategy_A"] = df_a.to_dict(orient="records")

    print("\n[Strategy B] StratifiedKFold(3) on label × molecular_subtype")
    summary["strategy_B"] = strategy_b_stratified(man, bags, labels)

    print("\n[Strategy C] Leave-one-out within 16 RAS-like slides")
    summary["strategy_C"] = strategy_c_loo_raslike(man, bags, labels)

    with (OUT / "RETRAIN_SUMMARY.json").open("w") as f:
        json.dump(summary, f, indent=2, default=str)
    print(f"\n[done] → {OUT}/RETRAIN_SUMMARY.json")


if __name__ == "__main__":
    main()
