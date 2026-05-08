"""Wave 8 — Evaluate each new feature block standalone, then merge.

For each block (TCR / self-similarity):
  1. Train logistic regression on TRAIN POOL (split=='train', n=2396) using only that block.
  2. Score ITSNdb_no_overlap (split in ['ext_itsndb_main','ext_itsndb_val'] AND in_master==False).
  3. AUROC + 95% CI bootstrap (1000 resamples).

Also run on full ITSNdb (combined) and ITSNdb_in_master for context.

Produces:
  - wave8_combined_features.tsv  (peptide | hla | tcr_* | self_*)
  - wave8_for_wave7.tsv          (same; ready to merge with wave7 features)
  - wave8_standalone_results.tsv (per-block AUROCs)
  - fig_wave8_modality_uplift.png/pdf
"""
from __future__ import annotations
import os, sys, time, json
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.metrics import roc_auc_score
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path("/home/seungho/personal/THCA_data_analysis/project/results/p_neo_bayesian_2026_05_09")
WAVE8 = ROOT / "wave8"
BUNDLE = ROOT / "bundle.tsv"
TCR = WAVE8 / "tcr_features.tsv"
SELFSIM = WAVE8 / "self_similarity_features.tsv"

COMBINED = WAVE8 / "wave8_combined_features.tsv"
FOR_WAVE7 = WAVE8 / "wave8_for_wave7.tsv"
RESULTS = WAVE8 / "wave8_standalone_results.tsv"
FIG_PNG = WAVE8 / "fig_wave8_modality_uplift.png"
FIG_PDF = WAVE8 / "fig_wave8_modality_uplift.pdf"
META = WAVE8 / "wave8_eval_meta.json"

RNG = np.random.default_rng(42)


def boot_ci(y, p, n=1000, seed=42):
    rng = np.random.default_rng(seed)
    n_obs = len(y)
    aurocs = []
    y = np.asarray(y); p = np.asarray(p)
    for _ in range(n):
        idx = rng.integers(0, n_obs, size=n_obs)
        yy = y[idx]; pp = p[idx]
        if yy.sum() == 0 or yy.sum() == n_obs:
            continue
        try:
            aurocs.append(roc_auc_score(yy, pp))
        except Exception:
            continue
    aurocs = np.asarray(aurocs)
    if len(aurocs) < 50:
        return float("nan"), float("nan"), float("nan")
    return float(np.mean(aurocs)), float(np.percentile(aurocs, 2.5)), float(np.percentile(aurocs, 97.5))


def main():
    t0 = time.time()
    print("[load] bundle...", flush=True)
    df = pd.read_csv(BUNDLE, sep="\t", low_memory=False)
    df = df[df["peptide"].apply(lambda s: isinstance(s, str))].copy()
    df["HLA_norm"] = df["HLA_norm"].fillna("").astype(str)

    print("[load] track8A TCR...", flush=True)
    tcr = pd.read_csv(TCR, sep="\t").drop_duplicates(["peptide", "hla"])
    print(f"[load] tcr rows={len(tcr)} (unique pep,hla) cols={list(tcr.columns)}", flush=True)

    print("[load] track8C self-similarity...", flush=True)
    sf = pd.read_csv(SELFSIM, sep="\t").drop_duplicates(["peptide", "hla"])
    print(f"[load] selfsim rows={len(sf)} (unique pep,hla) cols={list(sf.columns)}", flush=True)

    # Merge into bundle
    base = df[["peptide", "HLA_norm", "label", "source", "split", "in_master"]].copy()
    base["hla"] = base["HLA_norm"]
    merged = base.merge(tcr, on=["peptide", "hla"], how="left").merge(
        sf, on=["peptide", "hla"], how="left"
    )

    # Fill NaNs (any peptide that fell off either side)
    fill_cols = [c for c in merged.columns if c.startswith("tcr_") or c.startswith("self_")]
    for c in fill_cols:
        if merged[c].isna().any():
            med = merged[c].median(skipna=True)
            n_nan = merged[c].isna().sum()
            merged[c] = merged[c].fillna(med)
            print(f"[fillna] {c}: {n_nan} -> median={med:.4f}", flush=True)

    print(f"[merged] {len(merged)} rows, feature cols: {fill_cols}", flush=True)

    # Save merged feature matrix
    out_combined = merged[
        ["peptide", "hla", "label", "source", "split", "in_master"] + fill_cols
    ].copy()
    out_combined.to_csv(COMBINED, sep="\t", index=False)
    out_combined.to_csv(FOR_WAVE7, sep="\t", index=False)
    print(f"[write] {COMBINED}", flush=True)

    # Define test splits
    train_mask = merged["split"] == "train"
    itsndb_all_mask = merged["split"].isin(["ext_itsndb_main", "ext_itsndb_val"])
    no_overlap_mask = itsndb_all_mask & (merged["in_master"].astype(str).isin(["False", "false", "FALSE"]))
    in_master_mask = itsndb_all_mask & ~no_overlap_mask

    splits = {
        "ITSNdb_no_overlap": merged[no_overlap_mask].copy(),
        "ITSNdb_in_master": merged[in_master_mask].copy(),
        "ITSNdb_combined": merged[itsndb_all_mask].copy(),
    }
    for name, dat in splits.items():
        print(f"[split] {name}: n={len(dat)} pos={int((dat['label']==1).sum())}", flush=True)

    # Define feature blocks
    blocks = {
        "TCR_only": [c for c in fill_cols if c.startswith("tcr_")],
        "SelfSim_only": [c for c in fill_cols if c.startswith("self_")],
        "TCR+SelfSim": fill_cols,
    }

    # Eval
    train_df = merged[train_mask].copy()
    print(f"[train] n={len(train_df)} pos={int((train_df['label']==1).sum())}", flush=True)

    rows = []
    coef_log = {}
    for block_name, cols in blocks.items():
        Xtr = train_df[cols].to_numpy(dtype=float)
        ytr = train_df["label"].astype(int).to_numpy()
        if Xtr.shape[1] == 0 or len(np.unique(ytr)) < 2:
            continue
        clf = Pipeline([
            ("scaler", StandardScaler()),
            ("lr", LogisticRegression(max_iter=2000, C=1.0, solver="lbfgs")),
        ])
        clf.fit(Xtr, ytr)
        # Log coefficient signs
        coefs = dict(zip(cols, clf.named_steps["lr"].coef_.ravel().tolist()))
        coef_log[block_name] = coefs

        for split_name, dat in splits.items():
            X = dat[cols].to_numpy(dtype=float)
            y = dat["label"].astype(int).to_numpy()
            if len(y) == 0 or len(np.unique(y)) < 2:
                continue
            p = clf.predict_proba(X)[:, 1]
            auroc = roc_auc_score(y, p)
            mean_b, lo_b, hi_b = boot_ci(y, p, n=1000, seed=42)
            rows.append({
                "feature_block": block_name,
                "n_features": len(cols),
                "test_set": split_name,
                "n": int(len(y)),
                "n_pos": int(y.sum()),
                "AUROC": float(auroc),
                "boot_mean": float(mean_b),
                "CI_lo95": float(lo_b),
                "CI_hi95": float(hi_b),
            })
            print(f"[eval] {block_name} | {split_name}: AUROC={auroc:.4f} (95% CI {lo_b:.3f}-{hi_b:.3f})", flush=True)

    res_df = pd.DataFrame(rows)
    res_df.to_csv(RESULTS, sep="\t", index=False)
    print(f"[write] {RESULTS}", flush=True)

    # Figure: per modality bar with CI on no_overlap
    no_o = res_df[res_df["test_set"] == "ITSNdb_no_overlap"].copy()
    fig, ax = plt.subplots(figsize=(7, 4.5))
    xs = np.arange(len(no_o))
    bars = ax.bar(xs, no_o["AUROC"], color="#4c72b0")
    yerr_lo = no_o["AUROC"] - no_o["CI_lo95"]
    yerr_hi = no_o["CI_hi95"] - no_o["AUROC"]
    ax.errorbar(xs, no_o["AUROC"], yerr=[yerr_lo, yerr_hi], fmt="none", ecolor="#222", capsize=4, lw=1.2)
    ax.axhline(0.5, color="#999", lw=1, ls="--")
    ax.axhline(0.55, color="#c44", lw=1, ls=":", label="0.55 threshold")
    ax.set_xticks(xs)
    ax.set_xticklabels(no_o["feature_block"], rotation=15)
    ax.set_ylabel("AUROC (ITSNdb_no_overlap)")
    ax.set_title("Wave 8 — Standalone modality AUROC (TRAIN→ITSNdb_no_overlap)")
    ax.set_ylim(0.3, 0.85)
    for b, v in zip(bars, no_o["AUROC"]):
        ax.text(b.get_x() + b.get_width()/2, v + 0.01, f"{v:.3f}", ha="center", fontsize=9)
    ax.legend(loc="upper right", fontsize=8)
    plt.tight_layout()
    plt.savefig(FIG_PNG, dpi=140)
    plt.savefig(FIG_PDF)
    plt.close()
    print(f"[fig] {FIG_PNG}", flush=True)

    META.write_text(json.dumps({
        "n_combined_rows": int(len(out_combined)),
        "feature_blocks": {k: v for k, v in blocks.items()},
        "coefficients": coef_log,
        "elapsed_sec": time.time() - t0,
    }, indent=2))
    print(f"[done] elapsed {time.time()-t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()
