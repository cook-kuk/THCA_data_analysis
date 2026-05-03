#!/usr/bin/env python3
"""Leave-one-slide-out Ridge / ElasticNet regression on tile embeddings.

Targets: DM1_like_score_resid (primary), RAI_8_score_resid, TDS_like_score_resid

For each fold:
  - StandardScaler fit on train embeddings only
  - Ridge / ElasticNet fit on train
  - predict held-out slide
  - Spearman, Pearson, R2, MAE
  - DM1_high (top quartile in TRAIN) AUROC on test
  - sensitivity: DM1_high in GLOBAL top quartile

Outputs:
  loso_metrics_resnet50.tsv         per (target, model, tile_size, fold)
  loso_predictions_resnet50.tsv.gz  per tile predictions
"""
from __future__ import annotations
import argparse, logging, os, time
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, pearsonr
from sklearn.linear_model import Ridge, ElasticNet
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import r2_score, mean_absolute_error, roc_auc_score

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger("loso")

ROOT = Path(os.environ.get("THCA_ROOT", "/home/seungho/personal/THCA_data_analysis"))
META = ROOT / "project/results/03_pathology_poc/tile_metadata_resid.tsv.gz"
TARGETS = ["DM1_like_score_resid", "RAI_8_score_resid", "TDS_like_score_resid"]


def safe_metric(fn, *args, default=np.nan):
    try: return float(fn(*args))
    except Exception: return default


def loso_run(meta: pd.DataFrame, X: np.ndarray, target: str,
             model_name: str, alpha: float = 1.0, l1: float = 0.5):
    rows = []
    preds = np.full(len(meta), np.nan)
    samples = sorted(meta.sample_id.unique())
    for held in samples:
        tr = meta.sample_id != held
        te = ~tr
        y_tr = meta.loc[tr, target].values.astype(float)
        y_te = meta.loc[te, target].values.astype(float)
        keep_tr = ~np.isnan(y_tr)
        keep_te = ~np.isnan(y_te)
        if keep_tr.sum() < 50 or keep_te.sum() < 5:
            log.warning("%s/%s/%s skip held=%s (tr=%d te=%d)",
                        target, model_name, held, keep_tr.sum(), keep_te.sum())
            continue
        Xtr = X[tr.values][keep_tr]
        Xte = X[te.values][keep_te]
        y_tr = y_tr[keep_tr]; y_te = y_te[keep_te]

        sc = StandardScaler().fit(Xtr)
        Xtr_s = sc.transform(Xtr); Xte_s = sc.transform(Xte)

        if model_name == "ridge":
            mdl = Ridge(alpha=alpha)
        elif model_name == "enet":
            mdl = ElasticNet(alpha=alpha, l1_ratio=l1, max_iter=5000)
        else:
            raise ValueError(model_name)
        mdl.fit(Xtr_s, y_tr)
        p = mdl.predict(Xte_s)
        # write back into preds
        te_idx = np.where(te.values)[0][keep_te]
        preds[te_idx] = p

        # metrics
        sp = safe_metric(lambda a, b: spearmanr(a, b).correlation, y_te, p)
        pe = safe_metric(lambda a, b: pearsonr(a, b)[0], y_te, p)
        r2 = safe_metric(r2_score, y_te, p)
        mae = safe_metric(mean_absolute_error, y_te, p)

        # DM1_high AUROC: train-defined threshold (top quartile of TRAIN)
        thr_tr = np.quantile(y_tr, 0.75)
        y_te_hi_train = (y_te >= thr_tr).astype(int)
        auc_tr_thr = (safe_metric(roc_auc_score, y_te_hi_train, p)
                      if 0 < y_te_hi_train.sum() < len(y_te_hi_train) else np.nan)
        # sensitivity: GLOBAL top quartile
        thr_g = np.quantile(meta[target].dropna().values, 0.75)
        y_te_hi_g = (y_te >= thr_g).astype(int)
        auc_g = (safe_metric(roc_auc_score, y_te_hi_g, p)
                 if 0 < y_te_hi_g.sum() < len(y_te_hi_g) else np.nan)

        rows.append({"target": target, "model": model_name, "fold_held": held,
                     "n_train": len(y_tr), "n_test": len(y_te),
                     "spearman_r": sp, "pearson_r": pe, "r2": r2, "mae": mae,
                     "auroc_DM1high_train_thr": auc_tr_thr,
                     "auroc_DM1high_global_thr": auc_g,
                     "thr_train_q75": thr_tr, "thr_global_q75": thr_g})
    df = pd.DataFrame(rows)
    return df, preds


def pooled_metrics(meta: pd.DataFrame, preds: np.ndarray, target: str) -> dict:
    y = meta[target].values.astype(float)
    keep = ~np.isnan(y) & ~np.isnan(preds)
    y = y[keep]; p = preds[keep]
    out = {
        "pooled_n": int(len(y)),
        "pooled_spearman_r": safe_metric(lambda a, b: spearmanr(a, b).correlation, y, p),
        "pooled_pearson_r": safe_metric(lambda a, b: pearsonr(a, b)[0], y, p),
        "pooled_r2": safe_metric(r2_score, y, p),
        "pooled_mae": safe_metric(mean_absolute_error, y, p),
    }
    thr_g = np.quantile(meta[target].dropna().values, 0.75)
    y_hi = (y >= thr_g).astype(int)
    out["pooled_auroc_DM1high_global_thr"] = (
        safe_metric(roc_auc_score, y_hi, p)
        if 0 < y_hi.sum() < len(y_hi) else np.nan)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--embeddings", required=True)
    ap.add_argument("--tile-size", type=int, required=True)
    ap.add_argument("--out-metrics",
                    default="project/results/03_pathology_poc/loso_metrics_resnet50.tsv")
    ap.add_argument("--out-preds",
                    default="project/results/03_pathology_poc/loso_predictions_resnet50.tsv.gz")
    ap.add_argument("--ridge-alpha", type=float, default=1.0)
    ap.add_argument("--append", action="store_true")
    args = ap.parse_args()

    meta_full = pd.read_csv(META, sep="\t")
    meta = meta_full[meta_full.tile_size == args.tile_size].reset_index(drop=True)
    log.info("meta: %d tiles size=%d", len(meta), args.tile_size)

    z = np.load(args.embeddings)
    emb_spots = z["spot_ids"]
    # align embeddings to meta order
    idx_map = {s: i for i, s in enumerate(emb_spots)}
    keep = meta.spot_id.isin(idx_map).values
    if keep.sum() != len(meta):
        log.warning("dropping %d tiles missing embeddings", (~keep).sum())
        meta = meta[keep].reset_index(drop=True)
    order = np.array([idx_map[s] for s in meta.spot_id.values])
    X = z["embeddings"][order]
    log.info("aligned embeddings: %s", X.shape)

    metrics_all = []
    preds_long = []
    for target in TARGETS:
        if target not in meta.columns:
            log.warning("target %s missing", target); continue
        for model_name in ["ridge", "enet"]:
            t0 = time.time()
            df, preds = loso_run(meta, X, target, model_name,
                                 alpha=(args.ridge_alpha if model_name == "ridge" else 0.001))
            df["tile_size"] = args.tile_size
            df["embedding"] = "resnet50_imagenet"
            pooled = pooled_metrics(meta, preds, target)
            df.attrs.update(pooled)
            log.info("%s/%s tile=%d  pooled_spearman=%.3f auroc_g=%.3f  (%.1fs)",
                     target, model_name, args.tile_size,
                     pooled["pooled_spearman_r"],
                     pooled["pooled_auroc_DM1high_global_thr"],
                     time.time() - t0)
            for k, v in pooled.items():
                df[k] = v
            metrics_all.append(df)
            preds_long.append(pd.DataFrame({
                "spot_id": meta.spot_id, "sample_id": meta.sample_id,
                "stage": meta.stage, "tile_size": args.tile_size,
                "target": target, "model": model_name,
                "y_obs": meta[target].values, "y_pred": preds,
            }))

    metrics = pd.concat(metrics_all, ignore_index=True)
    preds_df = pd.concat(preds_long, ignore_index=True)

    out_metrics = Path(args.out_metrics)
    out_preds = Path(args.out_preds)
    out_metrics.parent.mkdir(parents=True, exist_ok=True)
    if args.append and out_metrics.exists():
        prev = pd.read_csv(out_metrics, sep="\t")
        metrics = pd.concat([prev, metrics], ignore_index=True)
    metrics.to_csv(out_metrics, sep="\t", index=False)
    if args.append and out_preds.exists():
        prev = pd.read_csv(out_preds, sep="\t")
        preds_df = pd.concat([prev, preds_df], ignore_index=True)
    preds_df.to_csv(out_preds, sep="\t", index=False, compression="gzip")
    log.info("wrote %s and %s", out_metrics, out_preds)


if __name__ == "__main__":
    main()
