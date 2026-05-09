#!/usr/bin/env python3
"""
Phase-3 decisive retraining audit.

E6. Init-only variance — fix StratifiedKFold(3, random_state=42), vary
    model init seed (6 seeds). Isolates stochastic init noise from split.
E7. Label-shuffle null — globally shuffle DM1/DM2 labels (3 shuffles),
    retrain via StratifiedKFold(3). Should give RAS-like AUC ≈ 0.5.
E8. Leave-one-TSS-out — for each TSS, hold out all that TSS's slides,
    train on the rest. If model relies on TSS batch, RAS-like AUC
    collapses for the TSS that contains DM1 vs DM2 imbalance.
"""
from __future__ import annotations
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch
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
OUT = ROOT / "analysis_supp" / "audit_ras_auc100"
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
                     [["case_short", "histology_subtype", "molecular_subtype"]])
    man = man.merge(meta_case, left_on="submitter_id",
                    right_on="case_short", how="left")
    print(f"[load] {len(man)} slides")
    bags = load_features(FEAT_DIR, man["file_id"].tolist())
    labels = man["label"].values
    return man, bags, labels


def stratkfold_eval(bags, labels, strat_key, init_seed=42, fold_seed=42):
    skf = StratifiedKFold(n_splits=3, shuffle=True, random_state=fold_seed)
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
        return np.nan
    return float(roc_auc_score(y, p))


def e6_init_variance(man, bags, labels):
    print("\n=== E6. Init-only variance (StratifiedKFold split fixed) ===")
    strat_key = (man["label"].astype(str) + "_"
                 + man["molecular_subtype"].astype(str)).values
    rows = []
    for init_seed in [42, 1, 7, 13, 99, 2026]:
        t0 = time.time()
        oof = stratkfold_eval(bags, labels, strat_key,
                               init_seed=init_seed, fold_seed=42)
        overall = float(roc_auc_score(labels, oof))
        ras = subgroup_auc(man, oof, "molecular_subtype", "RAS_like")
        braf = subgroup_auc(man, oof, "molecular_subtype", "BRAF_like")
        elapsed = time.time() - t0
        rows.append({"init_seed": init_seed,
                     "overall_auc": overall,
                     "ras_like_auc": ras,
                     "braf_like_auc": braf,
                     "seconds": round(elapsed, 1)})
        print(f"  init_seed={init_seed}: overall={overall:.3f}  "
              f"RAS={ras:.3f}  BRAF={braf:.3f}  ({elapsed:.0f}s)")
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "e6_init_variance.tsv", sep="\t", index=False)
    print(f"  RAS-like median={df['ras_like_auc'].median():.3f}  "
          f"range=[{df['ras_like_auc'].min():.3f}, "
          f"{df['ras_like_auc'].max():.3f}]")
    return df.to_dict("records")


def e7_label_shuffle_null(man, bags, labels):
    """Train CLAM with shuffled labels (global shuffle, preserve marginal).
    Should give RAS-like AUC ≈ 0.5 across reps. If it doesn't, the
    architecture itself extracts class-correlated structure even from
    randomized labels — i.e. memorisation."""
    print("\n=== E7. Label-shuffle null (3 reps) ===")
    rng = np.random.default_rng(2026)
    rows = []
    for rep in range(3):
        shuf_labels = rng.permutation(labels)
        strat_key = (pd.Series(shuf_labels).astype(str) + "_"
                     + man["molecular_subtype"].astype(str)).values
        oof = stratkfold_eval(bags, shuf_labels, strat_key,
                               init_seed=42, fold_seed=42)
        overall = float(roc_auc_score(shuf_labels, oof))
        # subgroup AUC computed against shuffled labels
        ras_mask = man["molecular_subtype"].values == "RAS_like"
        if len(set(shuf_labels[ras_mask])) > 1:
            ras = float(roc_auc_score(shuf_labels[ras_mask], oof[ras_mask]))
        else:
            ras = np.nan
        rows.append({"rep": rep + 1, "overall_shuf_auc": overall,
                     "ras_like_shuf_auc": ras})
        print(f"  rep {rep+1}: overall={overall:.3f}  RAS={ras:.3f}")
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "e7_label_shuffle_null.tsv", sep="\t", index=False)
    return df.to_dict("records")


def e8_leave_one_tss_out(man, bags, labels):
    """For each TSS group, hold all its slides out, train on the rest.
    Compute test-fold AUC and aggregate-OOF AUC. Reveals TSS batch effect."""
    print("\n=== E8. Leave-one-TSS-out cross-validation ===")
    tss_list = sorted(man["tss"].unique())
    oof = np.full(len(labels), np.nan)
    rows = []
    for tss in tss_list:
        va_mask = man["tss"].values == tss
        va_idx = np.where(va_mask)[0]
        tr_idx = np.where(~va_mask)[0]
        tr_bags = [bags[i] for i in tr_idx]
        tr_lab = [int(labels[i]) for i in tr_idx]
        va_bags = [bags[i] for i in va_idx]
        va_lab = [int(labels[i]) for i in va_idx]
        t0 = time.time()
        _, probs = train_one_fold(tr_bags, tr_lab, va_bags, va_lab,
                                   epochs=EPOCHS, device=DEVICE, seed=42)
        for j, idx in enumerate(va_idx):
            oof[idx] = probs[j]
        if len(set(va_lab)) > 1:
            auc = float(roc_auc_score(va_lab, probs))
        else:
            auc = np.nan
        ras_n = int(((man.iloc[va_idx]["molecular_subtype"] == "RAS_like")
                    ).sum())
        ras_pos = int(((man.iloc[va_idx]["molecular_subtype"] == "RAS_like")
                       & (np.array(va_lab) == 1)).sum())
        rows.append({"tss_held_out": tss, "n_test": len(va_idx),
                     "n_pos": int(sum(va_lab)),
                     "n_ras_test": ras_n, "n_ras_pos_test": ras_pos,
                     "tss_holdout_auc": auc,
                     "seconds": round(time.time() - t0, 1)})
        print(f"  TSS={tss}: n={len(va_idx)} ({sum(va_lab)}+) "
              f"RAS={ras_n}({ras_pos}+)  auc={auc}")
    overall = float(roc_auc_score(labels, oof))
    ras_mask = man["molecular_subtype"].values == "RAS_like"
    ras_auc = (float(roc_auc_score(labels[ras_mask], oof[ras_mask]))
               if len(set(labels[ras_mask])) > 1 else np.nan)
    print(f"  Pooled LOTO overall AUC = {overall:.3f}")
    print(f"  Pooled LOTO RAS-like AUC = {ras_auc:.3f}")
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "e8_leave_one_tss_out.tsv", sep="\t", index=False)
    return {"per_tss": df.to_dict("records"),
            "pooled_overall_auc": overall,
            "pooled_ras_like_auc": ras_auc}


def main():
    man, bags, labels = load_dataset()
    summary = {}

    print("\n[E6] Starting (~6 init seeds × 3 folds, ~5-7 min)")
    summary["E6_init_variance"] = e6_init_variance(man, bags, labels)

    print("\n[E7] Starting (~3 reps × 3 folds, ~3-4 min)")
    summary["E7_label_shuffle_null"] = e7_label_shuffle_null(
        man, bags, labels)

    print("\n[E8] Starting (~13 TSS LOTO folds, ~6-8 min)")
    summary["E8_leave_one_tss_out"] = e8_leave_one_tss_out(
        man, bags, labels)

    with (OUT / "PHASE3_AUDIT_SUMMARY.json").open("w") as f:
        json.dump(summary, f, indent=2, default=str)
    print(f"\n[done] → {OUT}/PHASE3_AUDIT_SUMMARY.json")


if __name__ == "__main__":
    main()
