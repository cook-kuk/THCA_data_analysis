#!/usr/bin/env python3
"""
Audit: RAS-like AUC=1.000 — is it a data-split / shortcut artifact?

Checks (no model retraining yet):
  C1. RAS-like vs FVPTC subgroup overlap (same 16 slides?)
  C2. prob_DM1 distribution within RAS-like / FVPTC (raw separation)
  C3. Permutation null AUC distribution (label shuffle within subgroup)
  C4. Histology-only / molecular-only LogReg baseline (shortcut probe)
  C5. Per-fold composition — how many RAS-like positives in each test fold

Subsequent script (audit_stratified_retrain.py) re-runs CLAM with
StratifiedGroupKFold on (molecular_subtype × label) to test whether
random-shuffle KFold artificially inflated the RAS-like AUC.
"""
from __future__ import annotations
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold

ROOT = Path("/home/seungho/personal/THCA_data_analysis/project/results/"
            "p2_image_dm1_v2_foundation_clam_2026_05_07")
META = Path("/home/seungho/personal/THCA_data_analysis/project/metadata/"
            "sample_master_v3.tsv")
OUT = ROOT / "analysis_supp" / "audit_ras_auc100"
OUT.mkdir(parents=True, exist_ok=True)

RNG = np.random.default_rng(42)
N_PERM = 10000
N_BOOT = 1000


def load_predictions_with_subgroups() -> pd.DataFrame:
    preds = pd.read_csv(ROOT / "phase2_tcga_clam"
                        / "clam_per_slide_predictions.tsv", sep="\t")
    man = pd.read_csv(ROOT / "phase2_tcga_clam"
                      / "slide_manifest.tsv", sep="\t")
    df = preds.merge(
        man[["file_id", "case_id", "submitter_id"]],
        left_on="slide", right_on="file_id", how="left"
    )

    meta = pd.read_csv(META, sep="\t", low_memory=False)
    meta = meta[meta["dataset"] == "TCGA-THCA"].copy()
    # sample_id like TCGA-DJ-A2Q6-01A → reduce to TCGA-DJ-A2Q6 to match submitter_id
    meta["case_short"] = meta["sample_id"].str.extract(r"^(TCGA-[^-]+-[^-]+)")
    # one row per case (some cases appear with different sample suffixes; first row OK
    # because histology/molecular subtype is case-level for THCA)
    meta_case = (meta.sort_values("sample_id")
                     .drop_duplicates("case_short", keep="first")
                     [["case_short", "histology_subtype", "molecular_subtype"]])
    df = df.merge(meta_case, left_on="submitter_id",
                  right_on="case_short", how="left")
    return df


def c1_subgroup_overlap(df: pd.DataFrame) -> dict:
    ras = df[df["molecular_subtype"] == "RAS_like"]
    fv = df[df["histology_subtype"] == "FVPTC"]
    inter = ras["slide"].isin(fv["slide"]).sum()
    union = pd.concat([ras["slide"], fv["slide"]]).drop_duplicates()
    out = {
        "n_RAS_like": len(ras),
        "n_FVPTC": len(fv),
        "n_intersection": int(inter),
        "n_union": int(len(union)),
        "RAS_like_AND_FVPTC": int(inter),
        "RAS_like_NOT_FVPTC": int(len(ras) - inter),
        "FVPTC_NOT_RAS_like": int(len(fv) - inter),
    }
    print(f"\n[C1] RAS_like ∩ FVPTC = {inter}/{len(ras)} (RAS) "
          f"= {inter}/{len(fv)} (FVPTC); union={len(union)}")
    return out


def c2_prob_dist(df: pd.DataFrame) -> dict:
    rows = []
    for grp_col, grp_val in [("molecular_subtype", "RAS_like"),
                             ("molecular_subtype", "BRAF_like"),
                             ("histology_subtype", "FVPTC"),
                             ("histology_subtype", "cPTC")]:
        sub = df[df[grp_col] == grp_val].copy()
        if len(sub) == 0:
            continue
        for lbl, name in [(1, "DM1"), (0, "DM2")]:
            ssub = sub[sub["label"] == lbl]
            if len(ssub) == 0:
                continue
            rows.append({
                "group": f"{grp_col}={grp_val}",
                "label": name,
                "n": len(ssub),
                "min_prob_DM1": float(ssub["prob_DM1"].min()),
                "max_prob_DM1": float(ssub["prob_DM1"].max()),
                "mean_prob_DM1": float(ssub["prob_DM1"].mean()),
                "median_prob_DM1": float(ssub["prob_DM1"].median()),
            })
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "c2_prob_distribution.tsv", sep="\t", index=False)
    print(f"\n[C2] prob_DM1 by subgroup × label written")
    print(out.to_string(index=False))

    # also dump raw 16-slide RAS-like predictions
    ras = df[df["molecular_subtype"] == "RAS_like"].sort_values("prob_DM1")
    ras_out = ras[["slide", "submitter_id", "histology_subtype",
                   "label", "prob_DM1", "fold"]].copy()
    ras_out.to_csv(OUT / "c2_RAS_like_16_predictions.tsv",
                   sep="\t", index=False)
    return {"n_rows": len(out)}


def c3_permutation_null(df: pd.DataFrame) -> dict:
    """Within each subgroup (RAS_like / FVPTC), shuffle labels N_PERM times
    and compute AUC. Compare observed AUC against null distribution.

    Also do BOOTSTRAP CI (resample slides with replacement) — this is
    the *correct* uncertainty bound for the observed AUC.
    """
    results = {}
    for grp_col, grp_val in [("molecular_subtype", "RAS_like"),
                             ("histology_subtype", "FVPTC"),
                             ("molecular_subtype", "BRAF_like"),
                             ("histology_subtype", "cPTC")]:
        sub = df[df[grp_col] == grp_val].copy()
        if len(sub) < 2 or sub["label"].nunique() < 2:
            continue
        y = sub["label"].values
        p = sub["prob_DM1"].values
        observed = roc_auc_score(y, p)

        # PERMUTATION: shuffle labels (preserve marginal counts)
        null_aucs = np.empty(N_PERM)
        for i in range(N_PERM):
            y_perm = RNG.permutation(y)
            try:
                null_aucs[i] = roc_auc_score(y_perm, p)
            except Exception:
                null_aucs[i] = 0.5
        p_two = float(np.mean(np.abs(null_aucs - 0.5)
                              >= abs(observed - 0.5)))
        p_one = float(np.mean(null_aucs >= observed))

        # BOOTSTRAP: resample slides with replacement
        boot_aucs = []
        n = len(y)
        for _ in range(N_BOOT):
            idx = RNG.integers(0, n, size=n)
            yb = y[idx]; pb = p[idx]
            if len(np.unique(yb)) < 2:
                continue
            boot_aucs.append(roc_auc_score(yb, pb))
        boot_aucs = np.array(boot_aucs) if boot_aucs else np.array([np.nan])
        ci_lo = float(np.percentile(boot_aucs, 2.5))
        ci_hi = float(np.percentile(boot_aucs, 97.5))

        results[f"{grp_col}={grp_val}"] = {
            "n": int(len(y)),
            "n_pos": int(y.sum()),
            "n_neg": int(len(y) - y.sum()),
            "observed_auc": float(observed),
            "permutation_p_two_sided": p_two,
            "permutation_p_one_sided": p_one,
            "null_auc_mean": float(null_aucs.mean()),
            "null_auc_p99": float(np.percentile(null_aucs, 99)),
            "null_auc_max": float(null_aucs.max()),
            "bootstrap_ci95_lo": ci_lo,
            "bootstrap_ci95_hi": ci_hi,
            "n_perms": int(N_PERM),
            "n_boots": int(len(boot_aucs)),
        }
    with (OUT / "c3_permutation_test.json").open("w") as f:
        json.dump(results, f, indent=2)
    print(f"\n[C3] permutation null + bootstrap CI written")
    for k, v in results.items():
        print(f"  {k}: AUC={v['observed_auc']:.3f}  "
              f"perm-p={v['permutation_p_two_sided']:.4f}  "
              f"boot95=[{v['bootstrap_ci95_lo']:.3f},{v['bootstrap_ci95_hi']:.3f}]")
    return results


def c4_shortcut_baselines(df: pd.DataFrame) -> dict:
    """Can we predict prob_DM1 (or DM1 label) from histology / molecular
    subtype ALONE? If yes, the CLAM model's RAS-like advantage may be
    a histology shortcut, not a true image signal.
    """
    sub = df.dropna(subset=["histology_subtype", "molecular_subtype"]).copy()
    sub = sub[sub["histology_subtype"].isin(["FVPTC", "cPTC"])
              & sub["molecular_subtype"].isin(["RAS_like", "BRAF_like"])]

    out = {}
    # C4a: histology-only LogReg → label
    X = (sub["histology_subtype"] == "FVPTC").astype(int).values.reshape(-1, 1)
    y = sub["label"].values
    lr = LogisticRegression().fit(X, y)
    out["c4a_histology_only_LR_AUC"] = float(
        roc_auc_score(y, lr.predict_proba(X)[:, 1]))

    # C4b: molecular-only LogReg → label
    X = (sub["molecular_subtype"] == "RAS_like").astype(int).values.reshape(-1, 1)
    lr = LogisticRegression().fit(X, y)
    out["c4b_molecular_only_LR_AUC"] = float(
        roc_auc_score(y, lr.predict_proba(X)[:, 1]))

    # C4c: how well does histology PREDICT prob_DM1 (not label)? If R^2 high,
    # CLAM has internalized histology as a feature.
    p_dm1 = sub["prob_DM1"].values
    is_fv = (sub["histology_subtype"] == "FVPTC").astype(float).values
    out["c4c_corr_histology_predDM1"] = float(np.corrcoef(is_fv, p_dm1)[0, 1])

    # C4d: cross-tab of histology × label and within-group AUC chance
    ct = pd.crosstab(sub["histology_subtype"], sub["label"],
                     margins=True, margins_name="all")
    out["c4d_crosstab_histology_label"] = ct.to_dict()

    # C4e: combinatorial chance — within RAS-like (3 pos, 13 neg),
    # P(perfect ranking | random) = 1 / C(16, 3)
    from math import comb
    out["c4e_random_perfect_RASlike_prob"] = 1.0 / comb(16, 3)

    with (OUT / "c4_shortcut_baselines.json").open("w") as f:
        json.dump(out, f, indent=2, default=str)
    print(f"\n[C4] shortcut baselines:")
    for k, v in out.items():
        if isinstance(v, float):
            print(f"  {k}: {v:.4f}")
    return out


def c5_per_fold(df: pd.DataFrame) -> dict:
    """How many RAS-like positives ended up in each test fold?"""
    ras = df[df["molecular_subtype"] == "RAS_like"]
    fold_compo = (ras.groupby("fold")["label"]
                     .agg(n="count", n_pos="sum")
                     .reset_index())
    fold_compo["n_neg"] = fold_compo["n"] - fold_compo["n_pos"]
    fold_compo.to_csv(OUT / "c5_RAS_like_fold_composition.tsv",
                      sep="\t", index=False)
    print(f"\n[C5] RAS-like fold composition:")
    print(fold_compo.to_string(index=False))
    return fold_compo.to_dict("records")


def main():
    df = load_predictions_with_subgroups()
    df.to_csv(OUT / "merged_preds_subgroups.tsv", sep="\t", index=False)
    print(f"[load] {len(df)} OOF predictions with subgroup labels")
    print(df["molecular_subtype"].value_counts(dropna=False).to_string())
    print(df["histology_subtype"].value_counts(dropna=False).to_string())

    summary = {}
    summary["C1_subgroup_overlap"] = c1_subgroup_overlap(df)
    summary["C2_prob_dist"] = c2_prob_dist(df)
    summary["C3_permutation"] = c3_permutation_null(df)
    summary["C4_shortcut"] = c4_shortcut_baselines(df)
    summary["C5_per_fold"] = c5_per_fold(df)

    with (OUT / "AUDIT_SUMMARY.json").open("w") as f:
        json.dump(summary, f, indent=2, default=str)
    print(f"\n[done] AUDIT_SUMMARY.json written → {OUT}")


if __name__ == "__main__":
    main()
