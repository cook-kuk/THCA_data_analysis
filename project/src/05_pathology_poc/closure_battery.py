#!/usr/bin/env python3
"""Closure battery — A, B, D, E, F.

Reads:
  embeddings_resnet50_224.npz                 (cached)
  tile_metadata_resid.tsv.gz                  (resid labels)
  all_spots_scored.tsv.gz                     (QC + raw scores)
  per-sample scored.h5ad                      (gene matrix for B residualization variants)

Writes:
  project/results/03_pathology_poc/closure_battery_metrics.tsv
  project/results/03_pathology_poc/closure_battery_predictions.tsv.gz
  project/results/03_pathology_poc/closure_battery_summary.png

Sections:
  A. Raw target sensitivity         (DM1/RAI/TDS raw, ResNet50 224)
  B. Residualization variants       (5 variants on DM1)
  D. Stage / morphology sanity      (stage ordinal, ATC vs non-ATC, epi/prolif raw)
  E. Slide-level aggregation        (re-aggregate existing predictions)
  F. Simple image features baseline (RGB mean/std/entropy)

C is in a separate script (tile_size_448.py) because it needs new tiles.
"""
from __future__ import annotations
import logging, os, time, glob
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
import statsmodels.api as sm
from PIL import Image
from scipy.stats import spearmanr, pearsonr, rankdata
from sklearn.linear_model import Ridge, ElasticNet, LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import (
    r2_score, mean_absolute_error, roc_auc_score,
    accuracy_score, balanced_accuracy_score,
)
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger("closure")

ROOT = Path(os.environ.get("THCA_ROOT", "/home/seungho/personal/THCA_data_analysis"))
EMB_FILE = ROOT/"project/results/03_pathology_poc/embeddings_resnet50_224.npz"
META = ROOT/"project/results/03_pathology_poc/tile_metadata_resid.tsv.gz"
SPOTS = ROOT/"project/results/01_spatial_score/all_spots_scored.tsv.gz"
PROCESSED = ROOT/"project/data/processed/GSE250521"
OUT_METRICS = ROOT/"project/results/03_pathology_poc/closure_battery_metrics.tsv"
OUT_PREDS   = ROOT/"project/results/03_pathology_poc/closure_battery_predictions.tsv.gz"
OUT_PNG     = ROOT/"project/results/03_pathology_poc/closure_battery_summary.png"


def safe(fn, *args, default=np.nan):
    try: return float(fn(*args))
    except Exception: return default


def loso_pooled_regression(meta: pd.DataFrame, X: np.ndarray, y: np.ndarray,
                           model_name: str = "ridge", alpha: float = 1.0):
    """Returns (preds_full_array, fold_metrics_df)."""
    samples = sorted(meta.sample_id.unique())
    preds = np.full(len(meta), np.nan)
    rows = []
    for held in samples:
        tr = (meta.sample_id != held).values
        te = ~tr
        y_tr = y[tr]; y_te = y[te]
        keep_tr = ~np.isnan(y_tr); keep_te = ~np.isnan(y_te)
        if keep_tr.sum() < 50 or keep_te.sum() < 5: continue
        Xtr = X[tr][keep_tr]; Xte = X[te][keep_te]
        sc = StandardScaler().fit(Xtr)
        Xtr_s = sc.transform(Xtr); Xte_s = sc.transform(Xte)
        if model_name == "ridge":
            mdl = Ridge(alpha=alpha)
        else:
            mdl = ElasticNet(alpha=alpha, l1_ratio=0.5, max_iter=3000)
        mdl.fit(Xtr_s, y_tr[keep_tr])
        p = mdl.predict(Xte_s)
        te_idx = np.where(te)[0][keep_te]
        preds[te_idx] = p
        sp = safe(lambda a, b: spearmanr(a, b).correlation, y_te[keep_te], p)
        pe = safe(lambda a, b: pearsonr(a, b)[0], y_te[keep_te], p)
        thr = np.quantile(y_tr[keep_tr], 0.75)
        yhi = (y_te[keep_te] >= thr).astype(int)
        auc = (safe(roc_auc_score, yhi, p)
               if 0 < yhi.sum() < len(yhi) else np.nan)
        rows.append({"fold_held": held, "n_train": int(keep_tr.sum()),
                     "n_test": int(keep_te.sum()),
                     "spearman_r": sp, "pearson_r": pe, "auroc_q75_train": auc})
    return preds, pd.DataFrame(rows)


def pooled_metrics(y, preds, label="target"):
    keep = ~np.isnan(y) & ~np.isnan(preds)
    if keep.sum() < 50: return {"pooled_n": 0}
    yk = y[keep]; pk = preds[keep]
    out = {
        "pooled_n": int(keep.sum()),
        "pooled_spearman_r": safe(lambda a, b: spearmanr(a, b).correlation, yk, pk),
        "pooled_pearson_r": safe(lambda a, b: pearsonr(a, b)[0], yk, pk),
        "pooled_r2": safe(r2_score, yk, pk),
        "pooled_mae": safe(mean_absolute_error, yk, pk),
    }
    thr = np.quantile(yk, 0.75)
    yhi = (yk >= thr).astype(int)
    out["pooled_auroc_q75"] = (safe(roc_auc_score, yhi, pk)
                               if 0 < yhi.sum() < len(yhi) else np.nan)
    thr50 = np.quantile(yk, 0.5)
    yhi50 = (yk >= thr50).astype(int)
    out["pooled_auroc_q50"] = (safe(roc_auc_score, yhi50, pk)
                                if 0 < yhi50.sum() < len(yhi50) else np.nan)
    return out


def load_all_spots():
    """Load merged spot table with QC + raw scores."""
    s = pd.read_csv(SPOTS, sep="\t")
    s["log_counts"] = np.log1p(s["total_counts"])
    s["log_ngenes"] = np.log1p(s["n_genes_by_counts"])
    return s


def residualize_per_sample(df, score_col, covars):
    out = np.full(len(df), np.nan)
    df = df.reset_index(drop=True)
    for sid, sub in df.groupby("sample_id"):
        y = sub[score_col].values.astype(float)
        if len(covars) == 0:
            out[sub.index] = y - np.nanmean(y)
            continue
        X = sub[covars].values.astype(float)
        X = sm.add_constant(X, has_constant="add")
        try:
            res = sm.OLS(y, X, missing="drop").fit()
            out[sub.index] = y - res.predict(X)
        except Exception:
            out[sub.index] = y - np.nanmean(y)
    return out


def rank_norm_per_sample(df, score_col):
    out = np.full(len(df), np.nan)
    df = df.reset_index(drop=True)
    for sid, sub in df.groupby("sample_id"):
        y = sub[score_col].values.astype(float)
        if np.isnan(y).all():
            out[sub.index] = np.nan; continue
        r = rankdata(y, nan_policy="omit") / max(np.sum(~np.isnan(y)), 1)
        out[sub.index] = r - r.mean()
    return out


def load_embeddings():
    z = np.load(EMB_FILE)
    return z["spot_ids"], z["sample_ids"], z["stages"], z["embeddings"]


def align_meta_emb(meta, emb_spots, emb_X):
    idx_map = {s: i for i, s in enumerate(emb_spots)}
    keep = meta.spot_id.isin(idx_map).values
    meta_a = meta[keep].reset_index(drop=True)
    order = np.array([idx_map[s] for s in meta_a.spot_id.values])
    return meta_a, emb_X[order]


# ============ Section A: raw targets ============
def section_A(meta, X, spots, all_metrics, all_preds):
    log.info("[A] raw target sensitivity")
    # merge raw scores from spots
    raw_cols = ["DM1_like_score", "RAI_8_score", "TDS_like_score",
                "Epithelial_score", "Proliferation_score"]
    extra = ["sample_id", "spot_id"] + raw_cols
    s = spots[extra].drop_duplicates(["spot_id", "sample_id"])
    m = meta.merge(s, on=["spot_id", "sample_id"], how="left", suffixes=("", "_y"))
    for tgt in ["DM1_like_score", "RAI_8_score", "TDS_like_score"]:
        # prefer the merged spot value (raw); fall back to meta
        col = tgt if tgt in m.columns else f"{tgt}_y"
        y = m[col].values.astype(float)
        for model_name in ["ridge", "enet"]:
            t0 = time.time()
            preds, fold = loso_pooled_regression(m, X, y, model_name=model_name)
            pm = pooled_metrics(y, preds, tgt)
            row = {"section": "A", "experiment": f"raw_{tgt}", "model": model_name,
                   "tile_size": 224, "label_kind": "raw", **pm}
            all_metrics.append(row)
            all_preds.append(pd.DataFrame({
                "section":"A", "experiment": f"raw_{tgt}", "model": model_name,
                "spot_id": m.spot_id, "sample_id": m.sample_id,
                "y_obs": y, "y_pred": preds}))
            log.info("  raw %s/%s sp=%.3f auc(q75)=%.3f r2=%.2f (%.1fs)",
                     tgt, model_name, pm.get("pooled_spearman_r", np.nan),
                     pm.get("pooled_auroc_q75", np.nan), pm.get("pooled_r2", np.nan),
                     time.time()-t0)


# ============ Section B: residualization variants ============
def section_B(meta, X, spots, all_metrics, all_preds):
    log.info("[B] residualization variants on DM1_like_score")
    base = "DM1_like_score"
    s = spots[["sample_id", "spot_id", base, "log_counts", "log_ngenes",
               "Epithelial_score"]].drop_duplicates(["spot_id", "sample_id"])
    m = meta.merge(s, on=["spot_id", "sample_id"], how="left", suffixes=("", "_y"))

    variants = {
        "B1_resid_log_counts_only": ("DM1_like_score", ["log_counts"]),
        "B2_resid_log_ngenes_only": ("DM1_like_score", ["log_ngenes"]),
        "B3_no_resid":              ("DM1_like_score", []),         # mean-centered per sample
        "B4_rank_norm":             ("DM1_like_score", "rank"),
        "B5_resid_full_eptop50":    ("DM1_like_score", ["log_counts","log_ngenes"]),
    }
    for name, (col, covars) in variants.items():
        if name == "B5_resid_full_eptop50":
            # filter to epithelial top 50% per sample
            keep_idx = []
            for sid, sub in m.groupby("sample_id"):
                if "Epithelial_score" not in sub.columns or sub["Epithelial_score"].isna().all():
                    continue
                thr = np.nanmedian(sub["Epithelial_score"].values)
                keep_idx.extend(sub[sub["Epithelial_score"] >= thr].index.tolist())
            m_sub = m.loc[sorted(set(keep_idx))].reset_index(drop=True)
            X_sub = X[m.index.isin(keep_idx)]
            y = residualize_per_sample(m_sub, col, ["log_counts", "log_ngenes"])
            X_use = X_sub
            m_use = m_sub
        elif covars == "rank":
            y = rank_norm_per_sample(m, col)
            X_use = X; m_use = m
        else:
            y = residualize_per_sample(m, col, covars)
            X_use = X; m_use = m
        for model_name in ["ridge"]:
            t0 = time.time()
            preds, _ = loso_pooled_regression(m_use, X_use, y, model_name=model_name)
            pm = pooled_metrics(y, preds, col)
            row = {"section":"B", "experiment": name, "model": model_name,
                   "tile_size": 224, "label_kind": "resid_variant", **pm,
                   "n_used": len(m_use)}
            all_metrics.append(row)
            log.info("  %s sp=%.3f auc=%.3f r2=%.2f n=%d (%.1fs)",
                     name, pm.get("pooled_spearman_r", np.nan),
                     pm.get("pooled_auroc_q75", np.nan),
                     pm.get("pooled_r2", np.nan), len(m_use), time.time()-t0)


# ============ Section D: stage / morphology sanity ============
def section_D(meta, X, spots, all_metrics, all_preds):
    log.info("[D] stage/morphology sanity")
    stage_to_int = {"PT": 0, "PTC": 1, "LPTC": 2, "ATC": 3}
    s = spots[["sample_id", "spot_id", "Epithelial_score",
               "Proliferation_score"]].drop_duplicates(["spot_id", "sample_id"])
    m = meta.merge(s, on=["spot_id", "sample_id"], how="left", suffixes=("", "_y"))
    m["stage_int"] = m["stage"].map(stage_to_int)
    m["is_atc"] = (m["stage"] == "ATC").astype(int)

    # D1: stage ordinal Spearman (regress as continuous 0..3)
    y = m["stage_int"].values.astype(float)
    preds, _ = loso_pooled_regression(m, X, y, model_name="ridge")
    pm = pooled_metrics(y, preds, "stage_int")
    pm["accuracy"] = float(accuracy_score(y, np.round(np.clip(preds, 0, 3))))
    all_metrics.append({"section":"D", "experiment":"D1_stage_ordinal",
                        "model":"ridge", "tile_size":224,
                        "label_kind":"stage_int", **pm})
    log.info("  D1 stage ordinal sp=%.3f acc=%.3f",
             pm.get("pooled_spearman_r", np.nan), pm.get("accuracy", np.nan))
    all_preds.append(pd.DataFrame({"section":"D", "experiment":"D1_stage_ordinal",
        "model":"ridge", "spot_id":m.spot_id, "sample_id":m.sample_id,
        "y_obs": y, "y_pred": preds}))

    # D2: ATC vs non-ATC binary, LOSO logistic
    samples = sorted(m.sample_id.unique())
    preds_atc = np.full(len(m), np.nan)
    fold_aucs = []
    for held in samples:
        tr = (m.sample_id != held).values; te = ~tr
        y_tr = m.loc[tr, "is_atc"].values; y_te = m.loc[te, "is_atc"].values
        if y_tr.sum() == 0 or y_tr.sum() == len(y_tr): continue
        if y_te.sum() == 0 or y_te.sum() == len(y_te):
            # held-out is single-class — still emit predictions but skip AUROC
            sc = StandardScaler().fit(X[tr])
            mdl = LogisticRegression(C=0.1, max_iter=2000, solver="lbfgs").fit(sc.transform(X[tr]), y_tr)
            preds_atc[np.where(te)[0]] = mdl.predict_proba(sc.transform(X[te]))[:,1]
            continue
        sc = StandardScaler().fit(X[tr])
        mdl = LogisticRegression(C=0.1, max_iter=2000, solver="lbfgs").fit(sc.transform(X[tr]), y_tr)
        p = mdl.predict_proba(sc.transform(X[te]))[:,1]
        preds_atc[np.where(te)[0]] = p
        fold_aucs.append(safe(roc_auc_score, y_te, p))
    y = m["is_atc"].values.astype(float)
    keep = ~np.isnan(preds_atc)
    pooled_atc_auc = safe(roc_auc_score, y[keep], preds_atc[keep]) if keep.sum() > 0 else np.nan
    all_metrics.append({"section":"D", "experiment":"D2_atc_binary",
                        "model":"logreg_C0.1", "tile_size":224,
                        "label_kind":"is_atc", "pooled_n": int(keep.sum()),
                        "pooled_auroc_atc": pooled_atc_auc,
                        "fold_auc_median": float(np.nanmedian(fold_aucs)) if fold_aucs else np.nan,
                        "fold_auc_min": float(np.nanmin(fold_aucs)) if fold_aucs else np.nan,
                        "fold_auc_max": float(np.nanmax(fold_aucs)) if fold_aucs else np.nan,
                        "n_folds_with_auc": len(fold_aucs)})
    log.info("  D2 ATC vs non-ATC pooled_auc=%.3f fold_med=%.3f (n_folds=%d)",
             pooled_atc_auc, np.nanmedian(fold_aucs) if fold_aucs else np.nan, len(fold_aucs))
    all_preds.append(pd.DataFrame({"section":"D", "experiment":"D2_atc_binary",
        "model":"logreg", "spot_id":m.spot_id, "sample_id":m.sample_id,
        "y_obs": y, "y_pred": preds_atc}))

    # D3: epithelial raw / proliferation raw
    for tgt in ["Epithelial_score", "Proliferation_score"]:
        y = m[tgt].values.astype(float)
        preds, _ = loso_pooled_regression(m, X, y, model_name="ridge")
        pm = pooled_metrics(y, preds, tgt)
        all_metrics.append({"section":"D", "experiment": f"D3_{tgt}",
                            "model":"ridge", "tile_size":224,
                            "label_kind":"raw_morphology", **pm})
        log.info("  D3 %s sp=%.3f auc=%.3f", tgt,
                 pm.get("pooled_spearman_r", np.nan), pm.get("pooled_auroc_q75", np.nan))
        all_preds.append(pd.DataFrame({"section":"D", "experiment": f"D3_{tgt}",
            "model":"ridge", "spot_id":m.spot_id, "sample_id":m.sample_id,
            "y_obs": y, "y_pred": preds}))


# ============ Section E: slide-level aggregation ============
def section_E(all_metrics):
    log.info("[E] slide-level aggregation of existing predictions")
    p = pd.read_csv(ROOT/"project/results/03_pathology_poc/loso_predictions_resnet50.tsv.gz", sep="\t")
    for tgt in ["DM1_like_score_resid", "RAI_8_score_resid", "TDS_like_score_resid"]:
        sub = p[(p.target == tgt) & (p.model == "ridge")].dropna(subset=["y_pred", "y_obs"])
        if sub.empty: continue
        for agg_name, agg_fn in [
            ("mean", lambda g: (g["y_obs"].mean(), g["y_pred"].mean())),
            ("top25_mean", lambda g: (g["y_obs"].mean(),
                                      g["y_pred"].nlargest(max(1, len(g)//4)).mean())),
            ("highrisk_frac", lambda g: (
                (g["y_obs"] >= np.quantile(p[p.target==tgt]["y_obs"].dropna(), 0.75)).mean(),
                (g["y_pred"] >= np.quantile(g["y_pred"], 0.75)).mean())),
        ]:
            slides = sub.groupby("sample_id").apply(lambda g: agg_fn(g)).to_dict()
            obs = np.array([v[0] for v in slides.values()])
            pred = np.array([v[1] for v in slides.values()])
            sp = safe(lambda a,b: spearmanr(a,b).correlation, obs, pred)
            pe = safe(lambda a,b: pearsonr(a,b)[0], obs, pred)
            all_metrics.append({"section":"E", "experiment": f"E_{tgt}_{agg_name}",
                                "model":"ridge", "tile_size":224,
                                "label_kind":"slide_level",
                                "pooled_spearman_r": sp, "pooled_pearson_r": pe,
                                "n_slides": len(obs)})
            log.info("  %s/%s sp=%.3f pe=%.3f n=%d", tgt, agg_name, sp, pe, len(obs))


# ============ Section F: simple image features ============
def section_F(meta, all_metrics, all_preds):
    log.info("[F] simple image features baseline (RGB mean/std/entropy)")
    feats = np.zeros((len(meta), 9), dtype=np.float32)  # R,G,B mean + std + entropy_R/G/B
    t0 = time.time()
    for i, p in enumerate(meta.tile_path.values):
        img = np.asarray(Image.open(p).convert("RGB"))
        for c in range(3):
            arr = img[:,:,c].astype(np.float32) / 255.0
            feats[i, c]   = arr.mean()
            feats[i, c+3] = arr.std()
            # crude entropy via histogram
            h, _ = np.histogram(img[:,:,c], bins=16, range=(0,255))
            h = h / h.sum() if h.sum() > 0 else h
            ent = -np.sum(h[h > 0] * np.log2(h[h > 0])) if h.sum() > 0 else 0.0
            feats[i, c+6] = ent
        if (i+1) % 800 == 0:
            log.info("  features %d/%d (%.1f tile/s)", i+1, len(meta),
                     (i+1)/max(time.time()-t0, 1e-6))
    np.save(ROOT/"project/results/03_pathology_poc/simple_features.npy", feats)

    # LOSO Ridge on resid + raw DM1
    spots = load_all_spots()
    s = spots[["sample_id", "spot_id", "DM1_like_score"]].drop_duplicates(["spot_id","sample_id"])
    m = meta.merge(s, on=["spot_id","sample_id"], how="left", suffixes=("","_y"))
    for tgt_name, y_arr in [
        ("DM1_like_score_resid_simplefeat", m["DM1_like_score_resid"].values.astype(float)),
        ("DM1_like_score_raw_simplefeat", m["DM1_like_score"].values.astype(float)),
    ]:
        preds, _ = loso_pooled_regression(m, feats, y_arr, model_name="ridge")
        pm = pooled_metrics(y_arr, preds, tgt_name)
        all_metrics.append({"section":"F", "experiment": tgt_name,
                            "model":"ridge_simplefeat_9d", "tile_size":224,
                            "label_kind":"simple_image_features", **pm})
        log.info("  %s sp=%.3f auc=%.3f", tgt_name,
                 pm.get("pooled_spearman_r", np.nan),
                 pm.get("pooled_auroc_q75", np.nan))


# ============ Plot summary ============
def plot_summary(metrics_df: pd.DataFrame):
    fig, ax = plt.subplots(figsize=(12, 9))
    plot = metrics_df[metrics_df.pooled_spearman_r.notna()][
        ["section","experiment","pooled_spearman_r","pooled_auroc_q75"]].copy()
    plot["label"] = plot["section"] + " " + plot["experiment"]
    plot = plot.sort_values("pooled_spearman_r")
    y = np.arange(len(plot))
    ax.barh(y, plot.pooled_spearman_r, color="steelblue", label="Spearman r")
    ax.scatter(plot.pooled_auroc_q75 - 0.5, y, color="firebrick", label="AUROC q75 - 0.5", s=22)
    ax.set_yticks(y); ax.set_yticklabels(plot.label.values, fontsize=7)
    ax.axvline(0.30, color="green", lw=0.6, ls="--", label="GO Spearman ≥ 0.30")
    ax.axvline(0.20, color="orange", lw=0.6, ls="--", label="BORDERLINE 0.20")
    ax.axvline(0, color="k", lw=0.4)
    ax.set_xlabel("metric value (Spearman r) ; AUROC plotted as (auroc - 0.5)")
    ax.set_title("Closure battery summary")
    ax.legend(loc="lower right", fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT_PNG, dpi=140, bbox_inches="tight")
    log.info("wrote %s", OUT_PNG)


def main():
    log.info("=== closure battery start ===")
    t_total = time.time()

    spot_ids, sample_ids, stages, X = load_embeddings()
    log.info("embeddings: %s", X.shape)

    meta = pd.read_csv(META, sep="\t")
    meta = meta[meta.tile_size == 224].reset_index(drop=True)
    meta, X = align_meta_emb(meta, spot_ids, X)
    log.info("aligned meta+emb: %d tiles", len(meta))

    spots = load_all_spots()

    all_metrics = []
    all_preds = []

    section_A(meta, X, spots, all_metrics, all_preds)
    section_B(meta, X, spots, all_metrics, all_preds)
    section_D(meta, X, spots, all_metrics, all_preds)
    section_E(all_metrics)
    section_F(meta, all_metrics, all_preds)

    metrics_df = pd.DataFrame(all_metrics)
    metrics_df.to_csv(OUT_METRICS, sep="\t", index=False)
    log.info("wrote %s (%d rows)", OUT_METRICS, len(metrics_df))
    if all_preds:
        preds_df = pd.concat(all_preds, ignore_index=True)
        preds_df.to_csv(OUT_PREDS, sep="\t", index=False, compression="gzip")
        log.info("wrote %s", OUT_PREDS)
    plot_summary(metrics_df)
    log.info("=== closure battery done in %.1fs ===", time.time() - t_total)


if __name__ == "__main__":
    main()
