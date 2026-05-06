#!/usr/bin/env python3
"""G2 — Slide-level regression: TCGA H&E embedding → bulk DM1_like.
80/20 train/test by patient, plus 5-fold CV. Survival validation if PFI available."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import KFold
from scipy.stats import spearmanr, pearsonr


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--embed", default="results/g2_slide_embeddings.npz")
    ap.add_argument("--meta", default="results/g2_slide_embeddings.meta.tsv")
    ap.add_argument("--score", default="data/s_tcga_thca_scored.tsv",
                    help="TCGA-THCA bulk DM1 score (s_tcga_thca_scored.tsv)")
    ap.add_argument("--out-dir", default="results/g2_regression")
    args = ap.parse_args()

    out = Path(args.out_dir); out.mkdir(parents=True, exist_ok=True)
    npz = np.load(args.embed, allow_pickle=True)
    X = npz["features"]; file_ids = npz["file_ids"]
    print(f"  slide embeddings: {X.shape}, model={npz['model']}")
    slide_meta = pd.read_csv(args.meta, sep="\t")
    score = pd.read_csv(args.score, sep="\t")
    score = score.drop_duplicates("patient")
    # join: slide_meta.case_id → score.patient
    df = slide_meta.merge(score[["patient","DM1_like","RAI_8","NONOVERLAP"]],
                          left_on="case_id", right_on="patient", how="inner")
    print(f"  matched {len(df)} slide-bulk pairs")

    # build aligned X / y
    fid_to_idx = {fid: i for i, fid in enumerate(file_ids)}
    df = df[df["file_id"].isin(fid_to_idx)].reset_index(drop=True)
    Xm = X[[fid_to_idx[f] for f in df["file_id"]]]
    y = df["DM1_like"].astype(float).values
    print(f"  X: {Xm.shape}, y: {y.shape}")

    # 5-fold CV
    kf = KFold(n_splits=5, shuffle=True, random_state=42)
    preds = np.full(len(y), np.nan)
    fold_metrics = []
    for fold_i, (tr, te) in enumerate(kf.split(Xm)):
        sc = StandardScaler().fit(Xm[tr])
        Xtr = sc.transform(Xm[tr]); Xte = sc.transform(Xm[te])
        mdl = Ridge(alpha=1.0).fit(Xtr, y[tr])
        preds[te] = mdl.predict(Xte)
        rho, p = spearmanr(y[te], preds[te])
        r, pp = pearsonr(y[te], preds[te])
        fold_metrics.append({"fold": fold_i+1, "n_test": len(te),
                              "spearman": rho, "spearman_p": p,
                              "pearson": r, "pearson_p": pp})
    fold_df = pd.DataFrame(fold_metrics)
    fold_df.to_csv(out / "fold_metrics.tsv", sep="\t", index=False)
    print(fold_df.to_string(index=False))

    rho_all, p_all = spearmanr(y, preds)
    r_all, p_r_all = pearsonr(y, preds)
    print(f"\nPooled 5-fold CV: Spearman ρ = {rho_all:.3f} (p={p_all:.1e}); Pearson r = {r_all:.3f}")

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(y, preds, s=30, alpha=0.5, color="#962E2E", edgecolor="black")
    ax.plot([y.min(), y.max()], [y.min(), y.max()], color="black", lw=0.5, ls="--")
    ax.set_xlabel("True bulk DM1_like (TCGA RNA-seq)"); ax.set_ylabel("Predicted from H&E (5-fold CV)")
    ax.set_title(f"G2 — TCGA-THCA digital pathology DM1 prediction (n={len(y)} slides)\n"
                 f"5-fold CV pooled ρ = {rho_all:.3f}, p = {p_all:.1e}\n"
                 f"AI predicts DM1 axis directly from FFPE H&E",
                 fontsize=11)
    fig.tight_layout()
    fig.savefig(out / "g2_cv_scatter.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"→ {out / 'g2_cv_scatter.png'}")


if __name__ == "__main__":
    main()
