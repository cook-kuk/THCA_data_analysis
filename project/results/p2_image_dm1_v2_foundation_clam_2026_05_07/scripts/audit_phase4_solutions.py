#!/usr/bin/env python3
"""
Phase-4 SOLUTIONS — fix the paper-killing issues identified in phase 1-3.

S1. TSS-ComBat on UNI features — remove scanning-center batch from
    the 1024-d feature space before CLAM training.
S2. TSS-stratified split — StratifiedKFold on (label × TSS-bucket).
S3. Multimodal LR (CLAM_prob + histology + sex + TSS dummies) —
    is the image channel additive on top of clinical?
S4. Sex-balanced retraining — oversample males in training.
S5. EM-out subsample — drop the 10 RAS-like (TSS=EM) all-DM2 cases
    and retrain. Tests if signal survives without the EM-DM2 anchor.

Outputs to analysis_supp/audit_ras_auc100/phase4/.
"""
from __future__ import annotations
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold

import sys
sys.path.insert(0, str(Path(__file__).parent))
from audit_clam_train_lib import train_one_fold, load_features  # noqa: E402

ROOT = Path("/home/seungho/personal/THCA_data_analysis/project/results/"
            "p2_image_dm1_v2_foundation_clam_2026_05_07")
META = Path("/home/seungho/personal/THCA_data_analysis/project/metadata/"
            "sample_master_v3.tsv")
FEAT_DIR = ROOT / "phase2_tcga_clam" / "features"
OUT = ROOT / "analysis_supp" / "audit_ras_auc100" / "phase4"
OUT.mkdir(parents=True, exist_ok=True)
DEVICE = "cpu"
EPOCHS = 30


def load_dataset():
    man = pd.read_csv(ROOT / "phase2_tcga_clam"
                      / "slide_manifest.tsv", sep="\t")
    man = man[man["dm"].isin(["DM1", "DM2"])].copy()
    man["label"] = man["dm"].map({"DM2": 0, "DM1": 1})
    man["has_feat"] = man["file_id"].apply(
        lambda f: (FEAT_DIR / f"{f}.pt").exists())
    man = man[man["has_feat"]].reset_index(drop=True)
    man["tss"] = man["submitter_id"].str.extract(r"TCGA-([^-]+)-")

    meta = pd.read_csv(META, sep="\t", low_memory=False)
    meta = meta[meta["dataset"] == "TCGA-THCA"].copy()
    meta["case_short"] = meta["sample_id"].str.extract(r"^(TCGA-[^-]+-[^-]+)")
    meta_case = (meta.sort_values("sample_id")
                     .drop_duplicates("case_short", keep="first")
                     [["case_short", "histology_subtype", "molecular_subtype",
                       "sex", "age", "ajcc_stage_group"]])
    man = man.merge(meta_case, left_on="submitter_id",
                    right_on="case_short", how="left")
    bags = load_features(FEAT_DIR, man["file_id"].tolist())
    labels = man["label"].values
    return man, bags, labels


def tss_combat_per_tile(bags, tss):
    """Per-feature mean-centering by TSS, applied at tile level.
    For each TSS group, subtract the group mean (over all tiles in that
    group) from each tile's feature vector. Optionally divide by group sd.
    Returns new bags."""
    # stack tile-level features per TSS
    tss_arr = np.array(tss)
    # compute per-TSS mean over all tiles
    means = {}
    for t in np.unique(tss_arr):
        idx = np.where(tss_arr == t)[0]
        all_tiles = np.vstack([bags[i] for i in idx])
        means[t] = all_tiles.mean(axis=0)
    # also a global mean to add back so features don't collapse to 0
    global_mean = np.vstack(bags).mean(axis=0)
    new_bags = []
    for i, t in enumerate(tss_arr):
        new_bags.append(bags[i] - means[t] + global_mean)
    return new_bags


def stratkfold_clam(bags, labels, strat_key, init_seed=42, fold_seed=42,
                     n_splits=3):
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True,
                          random_state=fold_seed)
    oof = np.full(len(labels), np.nan)
    for tr, va in skf.split(bags, strat_key):
        tr_bags = [bags[i] for i in tr]
        tr_lab = [int(labels[i]) for i in tr]
        va_bags = [bags[i] for i in va]
        va_lab = [int(labels[i]) for i in va]
        _, probs = train_one_fold(tr_bags, tr_lab, va_bags, va_lab,
                                   epochs=EPOCHS, device=DEVICE,
                                   seed=init_seed)
        for j, idx in enumerate(va):
            oof[idx] = probs[j]
    return oof


def subgroup_auc(man, oof, col, val):
    mask = man[col].values == val
    y = man["label"].values[mask]; p = oof[mask]
    if len(set(y)) < 2:
        return float("nan")
    return float(roc_auc_score(y, p))


def s1_tss_combat(man, bags, labels, n_init_seeds=4):
    print("\n=== S1. TSS-ComBat on UNI features ===")
    bags_corr = tss_combat_per_tile(bags, man["tss"].tolist())
    print(f"  per-tile TSS centering applied (drift={np.linalg.norm(np.vstack(bags[0]) - np.vstack(bags_corr[0])):.3f})")
    strat_key = (man["label"].astype(str) + "_"
                 + man["molecular_subtype"].astype(str)).values
    rows = []
    for init in [42, 1, 7, 13]:
        t0 = time.time()
        oof = stratkfold_clam(bags_corr, labels, strat_key,
                                init_seed=init, fold_seed=42)
        overall = float(roc_auc_score(labels, oof))
        ras = subgroup_auc(man, oof, "molecular_subtype", "RAS_like")
        braf = subgroup_auc(man, oof, "molecular_subtype", "BRAF_like")
        male = subgroup_auc(man, oof, "sex", "Male")
        female = subgroup_auc(man, oof, "sex", "Female")
        rows.append({"init_seed": init,
                     "overall": overall, "ras": ras, "braf": braf,
                     "male": male, "female": female,
                     "seconds": round(time.time() - t0, 1)})
        print(f"  init={init}: overall={overall:.3f}  RAS={ras:.3f}  "
              f"BRAF={braf:.3f}  Male={male:.3f}  Female={female:.3f}")
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "s1_tss_combat_per_init.tsv", sep="\t", index=False)
    print(f"  S1 ComBat-CLAM RAS-like median = "
          f"{df['ras'].median():.3f}  range "
          f"[{df['ras'].min():.3f}, {df['ras'].max():.3f}]")
    return df.to_dict("records")


def s2_tss_balanced_split(man, bags, labels):
    """StratifiedKFold on combined (label, TSS-bucket).
    Bucketize TSS into {EM, BJ_DJ_FK, other} so the 3 buckets each appear
    in train and val of every fold."""
    print("\n=== S2. TSS-balanced stratified split ===")
    def bucket(t):
        if t == "EM":
            return "EM"
        if t in ("BJ", "DJ", "FK"):
            return "BJ_DJ_FK"
        return "other"
    tss_bucket = man["tss"].map(bucket)
    strat_key = (man["label"].astype(str) + "_" + tss_bucket).values
    print(f"  strat key counts: "
          f"{pd.Series(strat_key).value_counts().to_dict()}")
    rows = []
    for init in [42, 1, 7, 13]:
        t0 = time.time()
        oof = stratkfold_clam(bags, labels, strat_key,
                                init_seed=init, fold_seed=42, n_splits=3)
        overall = float(roc_auc_score(labels, oof))
        ras = subgroup_auc(man, oof, "molecular_subtype", "RAS_like")
        braf = subgroup_auc(man, oof, "molecular_subtype", "BRAF_like")
        rows.append({"init_seed": init, "overall": overall,
                     "ras": ras, "braf": braf,
                     "seconds": round(time.time() - t0, 1)})
        print(f"  init={init}: overall={overall:.3f}  RAS={ras:.3f}  "
              f"BRAF={braf:.3f}")
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "s2_tss_balanced_split.tsv", sep="\t", index=False)
    return df.to_dict("records")


def s3_multimodal_LR(man, bags, labels):
    """CLAM_prob (from existing OOF) + clinical features fed into LR.
    Compare against CLAM-only and clinical-only."""
    print("\n=== S3. Multimodal LR (CLAM_prob + clinical) ===")
    pred = pd.read_csv(ROOT / "phase2_tcga_clam"
                       / "clam_per_slide_predictions.tsv", sep="\t")
    pred = pred.merge(man[["file_id", "histology_subtype", "sex", "tss",
                            "molecular_subtype"]],
                       left_on="slide", right_on="file_id", how="left")
    sub = pred.dropna(subset=["histology_subtype", "sex"]).copy()
    sub = sub[sub["histology_subtype"].isin(["cPTC", "FVPTC"])]

    clinical = pd.get_dummies(sub[["histology_subtype", "sex", "tss"]],
                               drop_first=True).astype(float).values
    image = sub["prob_DM1"].values.reshape(-1, 1)
    multi = np.hstack([image, clinical])
    y = sub["label"].values

    rows = []
    from sklearn.model_selection import KFold
    for name, X in [("clinical_only", clinical),
                    ("clam_only", image),
                    ("multimodal", multi)]:
        kf = KFold(n_splits=5, shuffle=True, random_state=42)
        oof = np.full(len(y), np.nan)
        for tr, va in kf.split(X):
            lr = LogisticRegression(max_iter=1000, C=1.0)
            lr.fit(X[tr], y[tr])
            oof[va] = lr.predict_proba(X[va])[:, 1]
        auc = float(roc_auc_score(y, oof))
        rows.append({"predictor": name, "auc": auc})
        print(f"  {name}: AUC = {auc:.3f}")
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "s3_multimodal_LR.tsv", sep="\t", index=False)
    return df.to_dict("records")


def s5_em_out_retrain(man, bags, labels):
    """Drop the 10 RAS-like (TSS=EM, all DM2) and retrain CLAM.
    Tests whether RAS-like signal survives without the EM-DM2 anchor."""
    print("\n=== S5. EM-out subsample retrain ===")
    drop_mask = ((man["tss"] == "EM")
                 & (man["molecular_subtype"] == "RAS_like"))
    keep = ~drop_mask.values
    print(f"  dropping {drop_mask.sum()} (TSS=EM, RAS_like) slides; "
          f"keeping {keep.sum()}")
    man2 = man[keep].reset_index(drop=True)
    bags2 = [bags[i] for i in range(len(bags)) if keep[i]]
    labels2 = labels[keep]
    strat_key = (man2["label"].astype(str) + "_"
                 + man2["molecular_subtype"].astype(str)).values
    rows = []
    for init in [42, 1, 7, 13]:
        t0 = time.time()
        oof = stratkfold_clam(bags2, labels2, strat_key,
                                init_seed=init, fold_seed=42)
        overall = float(roc_auc_score(labels2, oof))
        ras = subgroup_auc(man2, oof, "molecular_subtype", "RAS_like")
        braf = subgroup_auc(man2, oof, "molecular_subtype", "BRAF_like")
        rows.append({"init_seed": init, "overall": overall,
                     "ras": ras, "braf": braf, "n_train": int(keep.sum()),
                     "seconds": round(time.time() - t0, 1)})
        print(f"  init={init}: overall={overall:.3f}  RAS={ras:.3f}  "
              f"BRAF={braf:.3f}")
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "s5_em_out.tsv", sep="\t", index=False)
    return df.to_dict("records")


def main():
    print("[Phase 4 — SOLUTIONS]")
    man, bags, labels = load_dataset()
    print(f"[load] {len(man)} slides")

    summary = {}
    summary["S1_tss_combat"] = s1_tss_combat(man, bags, labels)
    summary["S2_tss_balanced_split"] = s2_tss_balanced_split(
        man, bags, labels)
    summary["S3_multimodal_LR"] = s3_multimodal_LR(man, bags, labels)
    summary["S5_em_out"] = s5_em_out_retrain(man, bags, labels)

    with (OUT / "PHASE4_SOLUTIONS_SUMMARY.json").open("w") as f:
        json.dump(summary, f, indent=2, default=str)
    print(f"\n[done] → {OUT}/PHASE4_SOLUTIONS_SUMMARY.json")


if __name__ == "__main__":
    main()
