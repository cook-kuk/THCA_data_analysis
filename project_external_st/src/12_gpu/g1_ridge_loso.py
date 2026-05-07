#!/usr/bin/env python3
"""G1 — Ridge regression LOSO on H&E tile embeddings → DM1_like_resid.
Output: per-tile predictions + per-fold metrics + summary plot."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler
from scipy.stats import spearmanr, pearsonr
from sklearn.metrics import roc_auc_score


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--embed", default="results/embeddings.npz")
    ap.add_argument("--meta", default="data/all_tile_metadata.tsv.gz")
    ap.add_argument("--target", default="DM1_like_score_resid",
                    help="column name in metadata to predict")
    ap.add_argument("--out-dir", default="results/g1_ridge")
    ap.add_argument("--alpha", type=float, default=1.0)
    args = ap.parse_args()

    out = Path(args.out_dir); out.mkdir(parents=True, exist_ok=True)
    npz = np.load(args.embed, allow_pickle=True)
    X = npz["features"]
    spot_ids = npz["spot_id"]
    sample_ids = npz["sample_id"]
    print(f"  embeddings: {X.shape}, model={npz['model']}")

    meta = pd.read_csv(args.meta, sep="\t").set_index("spot_id")
    # align meta with embed order
    meta_aligned = meta.reindex(spot_ids)
    if args.target not in meta_aligned.columns:
        print(f"target {args.target} not in metadata; available numeric:")
        print(meta_aligned.select_dtypes(include="number").columns.tolist()[:30])
        return
    y = meta_aligned[args.target].astype(float).values
    valid = ~np.isnan(y)
    print(f"  valid target values: {valid.sum()}/{len(y)}")
    X = X[valid]; y = y[valid]; spot_ids = spot_ids[valid]; sample_ids = sample_ids[valid]
    sample_ids_arr = sample_ids
    unique_samples = np.unique(sample_ids_arr)
    print(f"  unique slides: {len(unique_samples)}")

    fold_metrics = []
    preds = np.full(len(y), np.nan)
    for held in unique_samples:
        train_mask = sample_ids_arr != held
        test_mask = sample_ids_arr == held
        if train_mask.sum() < 10 or test_mask.sum() < 5: continue
        scaler = StandardScaler().fit(X[train_mask])
        Xtr = scaler.transform(X[train_mask]); Xte = scaler.transform(X[test_mask])
        mdl = Ridge(alpha=args.alpha).fit(Xtr, y[train_mask])
        pred = mdl.predict(Xte)
        preds[test_mask] = pred
        rho, p = spearmanr(y[test_mask], pred)
        r, pp = pearsonr(y[test_mask], pred)
        fold_metrics.append({"slide": held, "n_test": int(test_mask.sum()),
                              "spearman": rho, "spearman_p": p,
                              "pearson": r, "pearson_p": pp})
    fold_df = pd.DataFrame(fold_metrics)
    fold_df.to_csv(out / "fold_metrics.tsv", sep="\t", index=False)

    # Pooled metric
    valid2 = ~np.isnan(preds)
    rho_all, p_all = spearmanr(y[valid2], preds[valid2])
    r_all, p_all_r = pearsonr(y[valid2], preds[valid2])
    # AUROC top 25 vs bottom 25
    q25, q75 = np.percentile(y[valid2], [25, 75])
    binary = np.full_like(y[valid2], -1, dtype=int)
    binary[y[valid2] <= q25] = 0
    binary[y[valid2] >= q75] = 1
    valid3 = binary >= 0
    if valid3.sum() > 10:
        try:
            auc = roc_auc_score(binary[valid3], preds[valid2][valid3])
        except Exception: auc = np.nan
    else:
        auc = np.nan

    summary = {
        "n_tiles": int(valid2.sum()), "n_slides": int(len(unique_samples)),
        "model": str(npz["model"]),
        "pooled_spearman": rho_all, "pooled_spearman_p": p_all,
        "pooled_pearson": r_all, "pooled_pearson_p": p_all_r,
        "auroc_top25_vs_bot25": auc,
        "mean_per_slide_spearman": fold_df["spearman"].mean() if len(fold_df) else np.nan,
        "median_per_slide_spearman": fold_df["spearman"].median() if len(fold_df) else np.nan,
    }
    pd.Series(summary).to_csv(out / "summary.tsv", sep="\t")
    print("\n=== G1 LOSO Ridge result ===")
    for k, v in summary.items(): print(f"  {k}: {v}")

    # Figure
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    ax = axes[0]
    ax.scatter(y[valid2], preds[valid2], s=5, alpha=0.4, color="#962E2E")
    ax.plot([y.min(), y.max()], [y.min(), y.max()], color="black", lw=0.5, ls="--")
    ax.set_xlabel(f"True {args.target}"); ax.set_ylabel("LOSO predicted")
    ax.set_title(f"G1.A — All-tile LOSO scatter\n"
                 f"pooled Spearman ρ = {rho_all:.3f} (p = {p_all:.1e}), AUROC = {auc:.3f}",
                 fontsize=11)
    ax = axes[1]
    fold_df_sorted = fold_df.sort_values("spearman")
    ax.barh(fold_df_sorted["slide"], fold_df_sorted["spearman"],
            color=["#962E2E" if x>0 else "#3C6B4F" for x in fold_df_sorted["spearman"]],
            edgecolor="black")
    ax.axvline(0, color="grey", lw=0.5, ls=":")
    ax.tick_params(axis="y", labelsize=7)
    ax.set_xlabel("Per-slide Spearman ρ")
    ax.set_title(f"G1.B — Per-slide LOSO ρ\nmean={fold_df['spearman'].mean():.3f}",
                 fontsize=11)
    ax = axes[2]
    if valid3.sum() > 10:
        from sklearn.metrics import roc_curve
        fpr, tpr, _ = roc_curve(binary[valid3], preds[valid2][valid3])
        ax.plot(fpr, tpr, color="#962E2E", lw=2)
        ax.plot([0,1], [0,1], color="grey", lw=0.5, ls=":")
        ax.set_xlabel("FPR"); ax.set_ylabel("TPR")
        ax.set_title(f"G1.C — ROC top25 vs bot25\nAUROC = {auc:.3f}", fontsize=11)
    fig.tight_layout()
    fig.savefig(out / "g1_loso_summary.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"→ {out / 'g1_loso_summary.png'}")


if __name__ == "__main__":
    main()
